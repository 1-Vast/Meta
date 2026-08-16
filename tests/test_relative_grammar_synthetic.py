"""Stage 1 gates for the signed relative-difference transport.

Each gate falsifies one property the Stage 4 admission failure showed the
previous transport did not have.
"""
import pytest
import torch

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from model.relative_grammar import RelativeDifferenceTransport, RelativeGrammarModel


PROTEIN_DIM = 16
RESIDUES = 6


def _model(seed: int = 9101, use_reliability: bool = True) -> RelativeGrammarModel:
    torch.manual_seed(seed)
    return RelativeGrammarModel(
        protein_dim=PROTEIN_DIM, hidden_dim=8, task_dim=4, ligand_layers=1,
        pair_dim=8, pair_latents=4, pair_heads=2,
        use_reliability=use_reliability)


def _episode(generator: torch.Generator, batch: int = 2, support: int = 3,
             query: int = 4, atoms: int = 5):
    def graph(count):
        raw = torch.rand(batch, count, atoms, ATOM_FEAT_DIM, generator=generator)
        bonds = torch.rand(batch, count, atoms, atoms, BOND_FEAT_DIM,
                           generator=generator)
        bonds = bonds * (bonds > 0.7)
        return raw, bonds, torch.ones(batch, count, atoms)
    support_graph = graph(support)
    query_graph = graph(query)
    return {
        "protein_pooled": torch.randn(batch, PROTEIN_DIM, generator=generator),
        "protein_tokens": torch.randn(batch, RESIDUES, PROTEIN_DIM,
                                      generator=generator),
        "protein_mask": torch.ones(batch, RESIDUES),
        "protein_chemistry": torch.rand(batch, RESIDUES, 4, generator=generator),
        "support_atoms": support_graph[0], "support_bonds": support_graph[1],
        "support_mask": support_graph[2],
        "support_y": torch.randn(batch, support, generator=generator),
        "query_atoms": query_graph[0], "query_bonds": query_graph[1],
        "query_mask": query_graph[2],
    }


def _forward(model, episode, **kwargs):
    return model(
        episode["protein_pooled"], episode["protein_tokens"],
        episode["protein_mask"], episode["support_atoms"],
        episode["support_bonds"], episode["support_mask"], episode["support_y"],
        episode["query_atoms"], episode["query_bonds"], episode["query_mask"],
        protein_chemistry=episode["protein_chemistry"], **kwargs)


# --------------------------------------------------------------------------
# Algebraic gates on the difference operator
# --------------------------------------------------------------------------

def test_difference_operator_is_exactly_antisymmetric_and_zero_on_the_diagonal():
    torch.manual_seed(9201)
    transport = RelativeDifferenceTransport(6, hidden_dim=16)
    left = torch.randn(3, 5, 6)
    right = torch.randn(3, 4, 6)
    forward = transport.delta(left, right)
    backward = transport.delta(right, left)
    assert torch.allclose(forward, -backward.transpose(1, 2), atol=1e-6, rtol=1e-5)
    self_delta = transport.delta(left, left)
    assert torch.allclose(self_delta.diagonal(dim1=1, dim2=2),
                          torch.zeros(3, 5), atol=1e-6)
    assert torch.allclose(self_delta, -self_delta.transpose(1, 2),
                          atol=1e-6, rtol=1e-5)


def test_transport_is_affine_in_the_support_labels():
    """Label locking: the residual enters linearly, the operator does not see it."""
    torch.manual_seed(9202)
    transport = RelativeDifferenceTransport(6, hidden_dim=16, use_reliability=False)
    support = torch.randn(3, 4, 6)
    query = torch.randn(3, 5, 6)
    residual = torch.randn(3, 4)
    with torch.no_grad():
        base, delta, weight, _ = transport(support, query, residual)
        scaled, _, weight_scaled, _ = transport(support, query, 2.0 * residual)
        zero, _, _, _ = transport(support, query, torch.zeros_like(residual))
    assert torch.allclose(weight, weight_scaled, atol=1e-6)
    assert torch.allclose(scaled - zero, 2.0 * (base - zero), atol=1e-5, rtol=1e-4)
    # With no residual evidence the transport is exactly the mean difference.
    assert torch.allclose(
        zero, torch.einsum("bqk,bqk->bq", weight, delta), atol=1e-6)


def test_zero_difference_operator_recovers_the_shrunken_support_mean():
    """Level-only abstention is the floor, not a special case."""
    torch.manual_seed(9203)
    transport = RelativeDifferenceTransport(6, hidden_dim=16, use_reliability=False)
    with torch.no_grad():
        for parameter in transport.pair.parameters():
            parameter.zero_()
        # Zeroing the key makes every similarity identical, hence flat weights.
        # (The temperature has a floor of 1.0 by construction and cannot do it.)
        transport.key.weight.zero_()
    support = torch.randn(3, 4, 6)
    query = torch.randn(3, 5, 6)
    residual = torch.randn(3, 4)
    with torch.no_grad():
        value, delta, weight, _ = transport(support, query, residual)
    assert torch.count_nonzero(delta) == 0
    assert torch.allclose(value, residual.mean(-1, keepdim=True).expand_as(value),
                          atol=1e-5, rtol=1e-4)


# --------------------------------------------------------------------------
# Whole-model structural gates
# --------------------------------------------------------------------------

def test_zero_support_is_exactly_the_zero_shot_endpoint():
    model = _model()
    episode = _episode(torch.Generator().manual_seed(7401))
    empty = dict(episode)
    for key, axis in (("support_atoms", 1), ("support_bonds", 1),
                      ("support_mask", 1), ("support_y", 1)):
        empty[key] = episode[key][:, :0]
    out = _forward(model, empty)
    assert torch.equal(out.prediction, out.zero_shot)
    assert torch.count_nonzero(out.adaptation) == 0
    frozen = _forward(model, episode, adapt=False)
    assert torch.equal(frozen.prediction, frozen.zero_shot)
    assert torch.equal(frozen.zero_shot, out.zero_shot)


def test_single_support_correction_is_query_specific():
    model = _model(9102)
    episode = _episode(torch.Generator().manual_seed(7402), support=1)
    out = _forward(model, episode)
    spread = float(out.adaptation.std(dim=-1).min())
    assert spread > 1e-4, "k=1 correction is a pure scalar level shift"
    shifted = dict(episode)
    shifted["support_y"] = episode["support_y"] + 1.0
    moved = _forward(model, shifted)
    delta = moved.prediction - out.prediction
    assert float(delta.abs().min()) > 0.0


def test_support_permutation_invariance():
    model = _model(9103)
    episode = _episode(torch.Generator().manual_seed(7403), support=4)
    base = _forward(model, episode)
    order = torch.tensor([2, 0, 3, 1])
    permuted = dict(episode)
    for key in ("support_atoms", "support_bonds", "support_mask", "support_y"):
        permuted[key] = episode[key].index_select(1, order)
    other = _forward(model, permuted)
    assert torch.allclose(base.prediction, other.prediction, atol=1e-5, rtol=1e-5)


def test_query_permutation_equivariance():
    model = _model(9104)
    episode = _episode(torch.Generator().manual_seed(7404), query=4)
    base = _forward(model, episode)
    order = torch.tensor([3, 1, 0, 2])
    permuted = dict(episode)
    for key in ("query_atoms", "query_bonds", "query_mask"):
        permuted[key] = episode[key].index_select(1, order)
    other = _forward(model, permuted)
    assert torch.allclose(base.prediction.index_select(1, order),
                          other.prediction, atol=1e-5, rtol=1e-5)


def test_label_permutation_changes_the_prediction():
    """The failure this mechanism exists to fix: support identity must matter."""
    model = _model(9105)
    episode = _episode(torch.Generator().manual_seed(7405), support=4)
    base = _forward(model, episode)
    permuted = dict(episode)
    permuted["support_y"] = episode["support_y"].roll(1, dims=-1)
    other = _forward(model, permuted)
    moved = (base.prediction - other.prediction).abs()
    assert float(moved.max()) > 1e-4, "support labels are interchangeable"


def test_padding_does_not_change_the_prediction():
    model = _model(9106)
    episode = _episode(torch.Generator().manual_seed(7406), atoms=5)
    base = _forward(model, episode)
    padded = dict(episode)
    for key in ("support_atoms", "query_atoms"):
        padded[key] = torch.nn.functional.pad(episode[key], (0, 0, 0, 3))
    for key in ("support_bonds", "query_bonds"):
        padded[key] = torch.nn.functional.pad(episode[key], (0, 0, 0, 3, 0, 3))
    for key in ("support_mask", "query_mask"):
        padded[key] = torch.nn.functional.pad(episode[key], (0, 3))
    other = _forward(model, padded)
    assert torch.allclose(base.prediction, other.prediction, atol=1e-5, rtol=1e-4)


def test_every_trainable_tensor_receives_gradient():
    model = _model(9107)
    episode = _episode(torch.Generator().manual_seed(7407), support=3)
    out = _forward(model, episode)
    out.prediction.square().mean().backward()
    dead = [name for name, parameter in model.named_parameters()
            if parameter.grad is None or float(parameter.grad.abs().sum()) == 0.0]
    assert dead == [], f"dead trainable branches: {dead}"


def test_query_labels_are_not_an_input():
    model = _model(9108)
    episode = _episode(torch.Generator().manual_seed(7408))
    with pytest.raises(TypeError):
        _forward(model, episode, query_y=torch.zeros(1))
    assert torch.isfinite(_forward(model, episode).prediction).all()


def test_coordinates_are_refused_because_no_deployment_pair_has_a_complex():
    model = _model(9109)
    episode = _episode(torch.Generator().manual_seed(7409))
    with pytest.raises(ValueError, match="common-frame"):
        _forward(model, episode,
                 geometry_coordinates=torch.zeros(2, 7, 11, 3),
                 geometry_edge_index=torch.zeros(2, 4, dtype=torch.long))


# --------------------------------------------------------------------------
# Synthetic mechanism gates
# --------------------------------------------------------------------------

def _signed_tasks(generator, kind, shift, tasks=192, supports=3, queries=6):
    types = shift.shape[0]
    support_type = torch.randint(types, (tasks, supports), generator=generator)
    query_type = torch.randint(types, (tasks, queries), generator=generator)
    support_embed = torch.nn.functional.one_hot(support_type, types).float()
    query_embed = torch.nn.functional.one_hot(query_type, types).float()
    level = torch.randn(tasks, 1, generator=generator)
    if kind == "signed":
        # Each ligand type has a signed offset from the task level.
        residual = level + shift[support_type]
        target = level + shift[query_type]
    elif kind == "level":
        residual = level.expand(tasks, supports).clone()
        target = level.expand(tasks, queries).clone()
    elif kind == "private":
        residual = level + shift[support_type]
        target = torch.randn(tasks, queries, generator=generator)
    else:
        raise ValueError(kind)
    return support_embed, query_embed, residual, target


def _fit(transport, generator, kind, shift, steps=500):
    optimizer = torch.optim.Adam(transport.parameters(), lr=3e-2)
    for _ in range(steps):
        support, query, residual, target = _signed_tasks(generator, kind, shift)
        value, _, _, _ = transport(support, query, residual)
        loss = (value - target).square().mean()
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()


def test_signed_relative_effects_are_recovered_beyond_the_level():
    torch.manual_seed(9301)
    generator = torch.Generator().manual_seed(9302)
    shift = torch.randn(5, generator=generator)
    transport = RelativeDifferenceTransport(5, hidden_dim=32, use_reliability=False)
    _fit(transport, torch.Generator().manual_seed(9303), "signed", shift)
    support, query, residual, target = _signed_tasks(
        torch.Generator().manual_seed(73201), "signed", shift, tasks=512)
    with torch.no_grad():
        value, _, _, _ = transport(support, query, residual)
    level = residual.mean(-1, keepdim=True).expand_as(target)
    assert (value - target).square().mean() < 0.5 * (level - target).square().mean()


def test_single_support_signed_effect_is_recovered():
    torch.manual_seed(9304)
    generator = torch.Generator().manual_seed(9305)
    shift = torch.randn(5, generator=generator)
    transport = RelativeDifferenceTransport(5, hidden_dim=32, use_reliability=False)
    _fit(transport, torch.Generator().manual_seed(9306), "signed", shift, steps=700)
    support, query, residual, target = _signed_tasks(
        torch.Generator().manual_seed(73202), "signed", shift, tasks=512,
        supports=1)
    with torch.no_grad():
        value, _, _, _ = transport(support, query, residual)
    level = residual.mean(-1, keepdim=True).expand_as(target)
    assert (value - target).square().mean() < 0.5 * (level - target).square().mean()


def test_level_only_task_is_matched_without_spurious_structure():
    torch.manual_seed(9307)
    generator = torch.Generator().manual_seed(9308)
    shift = torch.randn(5, generator=generator)
    transport = RelativeDifferenceTransport(5, hidden_dim=32, use_reliability=False)
    _fit(transport, torch.Generator().manual_seed(9309), "level", shift)
    support, query, residual, target = _signed_tasks(
        torch.Generator().manual_seed(73203), "level", shift, tasks=512)
    with torch.no_grad():
        value, _, _, _ = transport(support, query, residual)
    assert (value - target).square().mean() < 0.02 * target.square().mean()


def test_private_mechanism_is_rejected_on_heldout_tasks():
    torch.manual_seed(9310)
    generator = torch.Generator().manual_seed(9311)
    shift = torch.randn(5, generator=generator)
    transport = RelativeDifferenceTransport(5, hidden_dim=32, use_reliability=False)
    _fit(transport, torch.Generator().manual_seed(9312), "private", shift)
    support, query, residual, target = _signed_tasks(
        torch.Generator().manual_seed(73204), "private", shift, tasks=512)
    with torch.no_grad():
        value, _, _, _ = transport(support, query, residual)
    assert (value - target).square().mean() >= 0.95 * target.square().mean()


def test_reliability_downweights_a_corrupted_support_label():
    """AdaMBind-style task analysis without MAML: inconsistent labels lose credit."""
    torch.manual_seed(9313)
    generator = torch.Generator().manual_seed(9314)
    shift = torch.randn(5, generator=generator)
    transport = RelativeDifferenceTransport(5, hidden_dim=32, use_reliability=True)
    _fit(transport, torch.Generator().manual_seed(9315), "signed", shift, steps=500)
    support, query, residual, _ = _signed_tasks(
        torch.Generator().manual_seed(73205), "signed", shift, tasks=256,
        supports=4)
    corrupted = residual.clone()
    corrupted[:, 0] = corrupted[:, 0] + 6.0
    with torch.no_grad():
        _, _, _, clean_credit = transport(support, query, residual)
        _, _, _, dirty_credit = transport(support, query, corrupted)
    assert float(dirty_credit[:, 0].mean()) < float(clean_credit[:, 0].mean())
