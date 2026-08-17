"""Stage 1 held-out synthetic gates for the interaction-grammar candidate.

Every gate is a falsification test of one contract the Stage 0 audit showed the
retained baseline could not satisfy.
"""
import pytest
import torch

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from model.interaction_grammar import (
    InteractionGrammarModel, TransferabilityTransport,
)


PROTEIN_DIM = 16
RESIDUES = 6


def _model(seed: int = 5101) -> InteractionGrammarModel:
    torch.manual_seed(seed)
    return InteractionGrammarModel(
        protein_dim=PROTEIN_DIM, hidden_dim=8, task_dim=4, ligand_layers=1,
        pair_dim=8, pair_latents=4, pair_heads=2)


def _episode(generator: torch.Generator, batch: int = 2, support: int = 3,
             query: int = 4, atoms: int = 5):
    def graph(count):
        raw = torch.rand(batch, count, atoms, ATOM_FEAT_DIM, generator=generator)
        bonds = torch.rand(batch, count, atoms, atoms, BOND_FEAT_DIM,
                           generator=generator)
        bonds = bonds * (bonds > 0.7)
        mask = torch.ones(batch, count, atoms)
        return raw, bonds, mask
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
# Whole-model structural gates
# --------------------------------------------------------------------------

def test_zero_support_is_exactly_the_zero_shot_endpoint():
    model = _model()
    episode = _episode(torch.Generator().manual_seed(7301))
    empty = dict(episode)
    empty["support_atoms"] = episode["support_atoms"][:, :0]
    empty["support_bonds"] = episode["support_bonds"][:, :0]
    empty["support_mask"] = episode["support_mask"][:, :0]
    empty["support_y"] = episode["support_y"][:, :0]
    out = _forward(model, empty)
    assert torch.equal(out.prediction, out.zero_shot)
    assert torch.count_nonzero(out.adaptation) == 0
    frozen = _forward(model, episode, adapt=False)
    assert torch.equal(frozen.prediction, frozen.zero_shot)
    assert torch.equal(frozen.zero_shot, out.zero_shot)


def _k1_label_effect_spread(model, episode) -> float:
    """Std across queries of the response to a uniform k=1 label shift."""
    base = _forward(model, episode)
    shifted = dict(episode)
    shifted["support_y"] = episode["support_y"] + 1.0
    moved = _forward(model, shifted)
    delta = moved.prediction - base.prediction
    assert float(delta.abs().min()) > 0.0
    return float(delta.std(dim=-1).min())


def test_single_support_label_effect_is_query_specific():
    """k=1 must move different queries by different amounts, and be trainable."""
    model = _model(5102)
    episode = _episode(torch.Generator().manual_seed(7302), support=1)
    initial = _k1_label_effect_spread(model, episode)
    assert initial > 0.0, "k=1 correction is a pure scalar level shift"
    assert float(_forward(model, episode).sar_adaptation.abs().max()) > 1e-6
    # The query-specific k=1 channel must be reachable by gradient descent, not
    # merely nonzero at initialisation.
    # The gate is positive and bounded, so a feasible target rescales the
    # transported residual per query rather than flipping its sign.
    ramp = torch.linspace(0.4, 1.6, episode["query_atoms"].shape[1])
    target = (_forward(model, episode).adaptation.detach() * ramp)
    optimizer = torch.optim.Adam(model.transport.parameters(), lr=5e-2)
    for _ in range(120):
        out = _forward(model, episode)
        loss = (out.adaptation - target).square().mean()
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    trained = _k1_label_effect_spread(model, episode)
    assert trained > 20.0 * initial and trained > 1e-2


def test_support_permutation_invariance():
    model = _model(5103)
    episode = _episode(torch.Generator().manual_seed(7303), support=4)
    base = _forward(model, episode)
    order = torch.tensor([2, 0, 3, 1])
    permuted = dict(episode)
    for key in ("support_atoms", "support_bonds", "support_mask", "support_y"):
        permuted[key] = episode[key].index_select(1, order)
    other = _forward(model, permuted)
    assert torch.allclose(base.prediction, other.prediction, atol=1e-5, rtol=1e-5)


def test_query_permutation_equivariance():
    model = _model(5104)
    episode = _episode(torch.Generator().manual_seed(7304), query=4)
    base = _forward(model, episode)
    order = torch.tensor([3, 1, 0, 2])
    permuted = dict(episode)
    for key in ("query_atoms", "query_bonds", "query_mask"):
        permuted[key] = episode[key].index_select(1, order)
    other = _forward(model, permuted)
    assert torch.allclose(base.prediction.index_select(1, order),
                          other.prediction, atol=1e-5, rtol=1e-5)


def test_query_labels_never_enter_the_model():
    model = _model(5105)
    episode = _episode(torch.Generator().manual_seed(7305))
    signature = model.forward.__doc__ or ""
    assert "query_y" not in signature
    assert "query_y" not in _forward.__code__.co_consts.__str__()
    base = _forward(model, episode)
    # Support labels are the only affinity input; the model has no query-label
    # argument, so an intervention on the query labels cannot be expressed.
    with pytest.raises(TypeError):
        model(episode["protein_pooled"], episode["protein_tokens"],
              episode["protein_mask"], episode["support_atoms"],
              episode["support_bonds"], episode["support_mask"],
              episode["support_y"], episode["query_atoms"],
              episode["query_bonds"], episode["query_mask"],
              protein_chemistry=episode["protein_chemistry"],
              query_y=torch.zeros(1))
    assert torch.isfinite(base.prediction).all()


def test_support_labels_reach_every_transport_parameter():
    model = _model(5106)
    episode = _episode(torch.Generator().manual_seed(7306), support=3)
    out = _forward(model, episode)
    out.prediction.square().mean().backward()
    dead = [name for name, parameter in model.named_parameters()
            if parameter.grad is None or float(parameter.grad.abs().sum()) == 0.0]
    assert dead == [], f"dead trainable branches: {dead}"


def test_zero_residual_produces_exactly_zero_transport():
    """Level-only abstention: no residual evidence means no correction."""
    model = _model(5107)
    episode = _episode(torch.Generator().manual_seed(7307), support=3)
    frozen = _forward(model, episode, adapt=False)
    matched = dict(episode)
    matched["support_y"] = _forward(model, episode).support_evidence.new_zeros(
        episode["support_y"].shape)
    # Build labels that reproduce the endpoint exactly, so every residual is 0.
    support_zero = model.encode(
        episode["protein_pooled"], episode["protein_tokens"],
        episode["protein_mask"],
        torch.cat((episode["support_atoms"], episode["query_atoms"]), 1),
        torch.cat((episode["support_bonds"], episode["query_bonds"]), 1),
        torch.cat((episode["support_mask"], episode["query_mask"]), 1),
        episode["protein_chemistry"])[0][:, :episode["support_atoms"].shape[1]]
    matched["support_y"] = support_zero.detach()
    out = _forward(model, matched)
    assert torch.allclose(out.prediction, frozen.zero_shot, atol=1e-5, rtol=1e-5)
    assert float(out.adaptation.abs().max()) < 1e-5


# --------------------------------------------------------------------------
# Transport mechanism gates on held-out synthetic tasks
# --------------------------------------------------------------------------

def _sensitivity(types: int, generator: torch.Generator) -> torch.Tensor:
    return 0.6 + 0.8 * torch.rand(types, generator=generator)


def _tasks(generator, kind, sensitivity, tasks=192, supports=3, queries=6):
    types = sensitivity.shape[0]
    support_type = torch.randint(types, (tasks, supports), generator=generator)
    query_type = torch.randint(types, (tasks, queries), generator=generator)
    support_embed = torch.nn.functional.one_hot(support_type, types).float()
    query_embed = torch.nn.functional.one_hot(query_type, types).float()
    amplitude = torch.randn(tasks, 1, generator=generator)
    if kind == "shared":
        residual = amplitude * sensitivity[support_type]
        target = amplitude * sensitivity[query_type]
    elif kind == "level":
        residual = amplitude.expand(tasks, supports).clone()
        target = amplitude.expand(tasks, queries).clone()
    elif kind == "private":
        residual = amplitude * sensitivity[support_type]
        target = torch.randn(tasks, queries, generator=generator)
    else:
        raise ValueError(kind)
    return support_embed, query_embed, residual, target


def _fit(transport, generator, kind, sensitivity, steps=400):
    optimizer = torch.optim.Adam(transport.parameters(), lr=5e-2)
    for _ in range(steps):
        support_embed, query_embed, residual, target = _tasks(
            generator, kind, sensitivity)
        prediction, _, _ = transport(support_embed, query_embed, residual)
        loss = (prediction - target).square().mean()
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()


def test_shared_mechanism_recovery_beats_the_level_control():
    torch.manual_seed(6101)
    sensitivity = _sensitivity(5, torch.Generator().manual_seed(6102))
    transport = TransferabilityTransport(5, hidden_dim=32)
    _fit(transport, torch.Generator().manual_seed(6103), "shared", sensitivity)
    support_embed, query_embed, residual, target = _tasks(
        torch.Generator().manual_seed(73101), "shared", sensitivity, tasks=512)
    with torch.no_grad():
        prediction, _, _ = transport(support_embed, query_embed, residual)
    level = residual.mean(-1, keepdim=True).expand_as(target)
    assert (prediction - target).square().mean() < 0.6 * (
        level - target).square().mean()


def test_single_support_shared_mechanism_is_recovered():
    """The decisive k=1 gate: one support label, query-specific transfer."""
    torch.manual_seed(6104)
    sensitivity = _sensitivity(5, torch.Generator().manual_seed(6105))
    transport = TransferabilityTransport(5, hidden_dim=32)
    _fit(transport, torch.Generator().manual_seed(6106), "shared", sensitivity,
         steps=600)
    support_embed, query_embed, residual, target = _tasks(
        torch.Generator().manual_seed(73102), "shared", sensitivity, tasks=512,
        supports=1)
    with torch.no_grad():
        prediction, gate, _ = transport(support_embed, query_embed, residual)
    level = residual.mean(-1, keepdim=True).expand_as(target)
    assert (prediction - target).square().mean() < 0.6 * (
        level - target).square().mean()
    assert float(gate.std()) > 1e-3


def test_level_only_task_is_matched_without_spurious_structure():
    torch.manual_seed(6107)
    sensitivity = _sensitivity(5, torch.Generator().manual_seed(6108))
    transport = TransferabilityTransport(5, hidden_dim=32)
    _fit(transport, torch.Generator().manual_seed(6109), "level", sensitivity)
    support_embed, query_embed, residual, target = _tasks(
        torch.Generator().manual_seed(73103), "level", sensitivity, tasks=512)
    with torch.no_grad():
        prediction, _, _ = transport(support_embed, query_embed, residual)
    assert (prediction - target).square().mean() < 0.02 * target.square().mean()


def test_private_mechanism_is_rejected_on_heldout_tasks():
    torch.manual_seed(6110)
    sensitivity = _sensitivity(5, torch.Generator().manual_seed(6111))
    transport = TransferabilityTransport(5, hidden_dim=32)
    _fit(transport, torch.Generator().manual_seed(6112), "private", sensitivity)
    support_embed, query_embed, residual, target = _tasks(
        torch.Generator().manual_seed(73104), "private", sensitivity, tasks=512)
    with torch.no_grad():
        prediction, _, _ = transport(support_embed, query_embed, residual)
    assert (prediction - target).square().mean() >= 0.95 * target.square().mean()


# --------------------------------------------------------------------------
# Zero-shot trunk gate: can the endpoint express a protein-conditioned
# interaction at all?  The Stage 0 audit showed the retained trunk could not.
# --------------------------------------------------------------------------

def _bilinear_corpus(generator, tasks, types=4, atoms=6, residues=6,
                     queries=4, codebook=None, weights=None):
    atom_type = torch.randint(types, (tasks, queries, atoms), generator=generator)
    residue_type = torch.randint(types, (tasks, residues), generator=generator)
    raw_atoms = torch.zeros(tasks, queries, atoms, ATOM_FEAT_DIM)
    raw_atoms.scatter_(-1, atom_type.unsqueeze(-1), 1.0)
    tokens = codebook[residue_type]
    ligand_share = torch.nn.functional.one_hot(atom_type, types).float().mean(-2)
    pocket_share = torch.nn.functional.one_hot(residue_type, types).float().mean(-2)
    truth = (ligand_share * pocket_share[:, None, :] * weights).sum(-1)
    truth = (truth - truth.mean()) / truth.std()
    bonds = torch.zeros(tasks, queries, atoms, atoms, BOND_FEAT_DIM)
    bonds[..., 0] = 1.0
    return {
        "protein_pooled": tokens.mean(1),
        "protein_tokens": tokens,
        "protein_mask": torch.ones(tasks, residues),
        "protein_chemistry": torch.zeros(tasks, residues, 4),
        "support_atoms": raw_atoms[:, :0], "support_bonds": bonds[:, :0],
        "support_mask": torch.ones(tasks, 0, atoms),
        "support_y": torch.zeros(tasks, 0),
        "query_atoms": raw_atoms, "query_bonds": bonds,
        "query_mask": torch.ones(tasks, queries, atoms),
    }, truth


def test_zero_shot_trunk_learns_a_protein_conditioned_interaction():
    torch.manual_seed(8101)
    generator = torch.Generator().manual_seed(8102)
    types = 4
    codebook = torch.randn(types, PROTEIN_DIM, generator=generator)
    weights = torch.tensor([2.0, -1.5, 1.0, -0.5])
    model = InteractionGrammarModel(
        protein_dim=PROTEIN_DIM, hidden_dim=32, task_dim=8, ligand_layers=2,
        pair_dim=32, pair_latents=8, pair_heads=4)
    optimizer = torch.optim.Adam(model.parameters(), lr=3e-3)
    for _ in range(240):
        episode, truth = _bilinear_corpus(
            generator, 24, types, codebook=codebook, weights=weights)
        out = _forward(model, episode, adapt=False)
        loss = (out.prediction - truth).square().mean()
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()
    heldout, truth = _bilinear_corpus(
        torch.Generator().manual_seed(73106), 256, types,
        codebook=codebook, weights=weights)
    with torch.no_grad():
        out = _forward(model, heldout, adapt=False)
        shuffled = dict(heldout)
        order = torch.randperm(truth.shape[0],
                               generator=torch.Generator().manual_seed(73107))
        for key in ("protein_pooled", "protein_tokens", "protein_mask",
                    "protein_chemistry"):
            shuffled[key] = heldout[key].index_select(0, order)
        wrong = _forward(model, shuffled, adapt=False)
    variance = float(truth.var())
    fit = float((out.prediction - truth).square().mean())
    swapped = float((wrong.prediction - truth).square().mean())
    assert fit < 0.5 * variance, f"trunk cannot fit the interaction: {fit}"
    assert swapped > fit + 0.15 * variance, (
        "the trunk ignores the protein: "
        f"correct {fit:.4f} versus shuffled {swapped:.4f}")


def test_transport_is_linear_in_the_support_labels():
    """Label locking: doubling every residual doubles the correction."""
    torch.manual_seed(6113)
    transport = TransferabilityTransport(5, hidden_dim=16)
    generator = torch.Generator().manual_seed(73105)
    support_embed, query_embed, residual, _ = _tasks(
        generator, "shared", _sensitivity(5, generator), tasks=8)
    with torch.no_grad():
        single, _, _ = transport(support_embed, query_embed, residual)
        double, _, _ = transport(support_embed, query_embed, 2.0 * residual)
    assert torch.allclose(double, 2.0 * single, atol=1e-6, rtol=1e-5)
