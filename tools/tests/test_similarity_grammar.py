"""Stage 6 gates for chemistry-grounded support weighting."""
import pytest
import torch

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from model.similarity_grammar import (
    SimilarityGrammarModel, SimilarityTransport, tanimoto,
)


PROTEIN_DIM = 16
RESIDUES = 6
BITS = 32


def _model(seed: int = 4101, **kwargs) -> SimilarityGrammarModel:
    torch.manual_seed(seed)
    return SimilarityGrammarModel(
        protein_dim=PROTEIN_DIM, hidden_dim=8, task_dim=4, ligand_layers=1,
        pair_dim=8, pair_latents=4, pair_heads=2, **kwargs)


def _episode(generator: torch.Generator, batch: int = 2, support: int = 3,
             query: int = 4, atoms: int = 5):
    def graph(count):
        raw = torch.rand(batch, count, atoms, ATOM_FEAT_DIM, generator=generator)
        bonds = torch.rand(batch, count, atoms, atoms, BOND_FEAT_DIM,
                           generator=generator)
        return raw, bonds * (bonds > 0.7), torch.ones(batch, count, atoms)
    support_graph, query_graph = graph(support), graph(query)
    def fingerprint(count):
        return (torch.rand(batch, count, BITS, generator=generator) > 0.6).float()
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
        "support_fingerprint": fingerprint(support),
        "query_fingerprint": fingerprint(query),
    }


def _forward(model, episode, **kwargs):
    return model(
        episode["protein_pooled"], episode["protein_tokens"],
        episode["protein_mask"], episode["support_atoms"],
        episode["support_bonds"], episode["support_mask"], episode["support_y"],
        episode["query_atoms"], episode["query_bonds"], episode["query_mask"],
        protein_chemistry=episode["protein_chemistry"],
        support_fingerprint=episode["support_fingerprint"],
        query_fingerprint=episode["query_fingerprint"], **kwargs)


def test_tanimoto_matches_the_set_definition():
    a = torch.tensor([[[1.0, 1, 0, 0], [1, 1, 1, 1]]])
    b = torch.tensor([[[1.0, 0, 0, 0], [0, 0, 1, 1]]])
    value = tanimoto(a, b)
    assert torch.allclose(value[0, 0, 0], torch.tensor(0.5))
    assert torch.allclose(value[0, 0, 1], torch.tensor(0.0))
    assert torch.allclose(value[0, 1, 0], torch.tensor(0.25))
    assert torch.allclose(value[0, 1, 1], torch.tensor(0.5))
    assert torch.allclose(tanimoto(a, a).diagonal(dim1=1, dim2=2),
                          torch.ones(1, 2))


def test_zero_fingerprint_contributes_no_similarity():
    """Unparsable ligands must add no evidence rather than wrong evidence."""
    a = torch.zeros(1, 1, 4)
    b = torch.tensor([[[1.0, 1, 0, 0]]])
    assert float(tanimoto(a, b)) == 0.0


def test_zero_support_is_exactly_the_zero_shot_endpoint():
    model = _model()
    episode = _episode(torch.Generator().manual_seed(4201))
    empty = dict(episode)
    for key in ("support_atoms", "support_bonds", "support_mask", "support_y",
                "support_fingerprint"):
        empty[key] = episode[key][:, :0]
    out = _forward(model, empty)
    assert torch.equal(out.prediction, out.zero_shot)
    assert torch.count_nonzero(out.adaptation) == 0
    frozen = _forward(model, episode, adapt=False)
    assert torch.equal(frozen.prediction, frozen.zero_shot)


def test_flat_weighting_recovers_the_shrunken_support_mean():
    """Level-only abstention stays the floor."""
    torch.manual_seed(4202)
    transport = SimilarityTransport(6, hidden_dim=8, use_learned_key=False)
    with torch.no_grad():
        transport.similarity_scale.zero_()
    support = torch.randn(3, 4, 6)
    query = torch.randn(3, 5, 6)
    residual = torch.randn(3, 4)
    similarity = torch.rand(3, 5, 4)
    value, weight = transport(support, query, residual, similarity)
    assert torch.allclose(weight, torch.full_like(weight, 0.25), atol=1e-6)
    assert torch.allclose(
        value, residual.mean(-1, keepdim=True).expand_as(value), atol=1e-6)


def test_similarity_changes_the_weighting():
    torch.manual_seed(4203)
    transport = SimilarityTransport(6, hidden_dim=8, use_learned_key=False)
    support = torch.randn(2, 3, 6)
    query = torch.randn(2, 2, 6)
    residual = torch.tensor([[1.0, -1.0, 0.0], [2.0, 0.0, -2.0]])
    peaked = torch.zeros(2, 2, 3)
    peaked[:, :, 0] = 1.0
    value, weight = transport(support, query, residual, peaked)
    assert float(weight[..., 0].min()) > float(weight[..., 1].max())
    assert float(value[0, 0]) > float(residual[0].mean())


def test_transport_is_linear_in_the_support_labels():
    torch.manual_seed(4204)
    transport = SimilarityTransport(6, hidden_dim=8)
    support, query = torch.randn(2, 3, 6), torch.randn(2, 4, 6)
    residual, similarity = torch.randn(2, 3), torch.rand(2, 4, 3)
    with torch.no_grad():
        single, weight_a = transport(support, query, residual, similarity)
        double, weight_b = transport(support, query, 2 * residual, similarity)
    assert torch.allclose(weight_a, weight_b, atol=1e-7)
    assert torch.allclose(double, 2 * single, atol=1e-6)


def test_support_permutation_invariance():
    model = _model(4102)
    episode = _episode(torch.Generator().manual_seed(4205), support=4)
    base = _forward(model, episode)
    order = torch.tensor([2, 0, 3, 1])
    permuted = dict(episode)
    for key in ("support_atoms", "support_bonds", "support_mask", "support_y",
                "support_fingerprint"):
        permuted[key] = episode[key].index_select(1, order)
    assert torch.allclose(base.prediction, _forward(model, permuted).prediction,
                          atol=1e-5, rtol=1e-5)


def test_query_permutation_equivariance():
    model = _model(4103)
    episode = _episode(torch.Generator().manual_seed(4206), query=4)
    base = _forward(model, episode)
    order = torch.tensor([3, 1, 0, 2])
    permuted = dict(episode)
    for key in ("query_atoms", "query_bonds", "query_mask", "query_fingerprint"):
        permuted[key] = episode[key].index_select(1, order)
    assert torch.allclose(base.prediction.index_select(1, order),
                          _forward(model, permuted).prediction,
                          atol=1e-5, rtol=1e-5)


def test_label_permutation_changes_the_prediction():
    model = _model(4104)
    episode = _episode(torch.Generator().manual_seed(4207), support=4)
    base = _forward(model, episode)
    permuted = dict(episode)
    permuted["support_y"] = episode["support_y"].roll(1, dims=-1)
    moved = (base.prediction - _forward(model, permuted).prediction).abs()
    assert float(moved.max()) > 1e-4, "support labels are interchangeable"


def test_every_trainable_tensor_receives_gradient():
    model = _model(4105)
    episode = _episode(torch.Generator().manual_seed(4208), support=3)
    _forward(model, episode).prediction.square().mean().backward()
    dead = [name for name, parameter in model.named_parameters()
            if parameter.grad is None or float(parameter.grad.abs().sum()) == 0.0]
    assert dead == [], f"dead trainable branches: {dead}"


def test_similarity_only_has_no_dead_trainable_parameter():
    """`use_learned_key=False` must not leave gradient-free trainable tensors."""
    model = _model(4108, use_learned_key=False)
    episode = _episode(torch.Generator().manual_seed(4213), support=3)
    _forward(model, episode).prediction.square().mean().backward()
    trainable = [(n, p) for n, p in model.named_parameters() if p.requires_grad]
    dead = [n for n, p in trainable
            if p.grad is None or float(p.grad.abs().sum()) == 0.0]
    assert dead == [], f"dead trainable branches: {dead}"
    # The unused tensors are retained for state-dict compatibility but frozen.
    frozen = {n for n, p in model.named_parameters() if not p.requires_grad}
    assert frozen == {"transport.key.weight", "transport.log_temperature"}
    assert model.transport.similarity_scale.requires_grad
    assert model.transport.log_shrinkage.requires_grad


def test_similarity_only_checkpoints_remain_loadable():
    """Freezing must not change state-dict keys or initialisation."""
    torch.manual_seed(4109)
    learned = SimilarityGrammarModel(
        protein_dim=PROTEIN_DIM, hidden_dim=8, task_dim=4, ligand_layers=1,
        pair_dim=8, pair_latents=4, pair_heads=2, use_learned_key=True)
    torch.manual_seed(4109)
    frozen = SimilarityGrammarModel(
        protein_dim=PROTEIN_DIM, hidden_dim=8, task_dim=4, ligand_layers=1,
        pair_dim=8, pair_latents=4, pair_heads=2, use_learned_key=False)
    assert set(learned.state_dict()) == set(frozen.state_dict())
    for key, value in learned.state_dict().items():
        assert torch.equal(value, frozen.state_dict()[key]), key
    frozen.load_state_dict(learned.state_dict())


def test_missing_fingerprints_fail_loudly():
    model = _model(4106)
    episode = _episode(torch.Generator().manual_seed(4209))
    with pytest.raises(ValueError, match="fingerprint"):
        model(episode["protein_pooled"], episode["protein_tokens"],
              episode["protein_mask"], episode["support_atoms"],
              episode["support_bonds"], episode["support_mask"],
              episode["support_y"], episode["query_atoms"],
              episode["query_bonds"], episode["query_mask"],
              protein_chemistry=episode["protein_chemistry"])


def test_coordinates_are_refused():
    model = _model(4107)
    episode = _episode(torch.Generator().manual_seed(4210))
    with pytest.raises(ValueError, match="common-frame"):
        _forward(model, episode,
                 geometry_coordinates=torch.zeros(2, 7, 11, 3),
                 geometry_edge_index=torch.zeros(2, 4, dtype=torch.long))


def test_similarity_weighting_recovers_a_synthetic_local_sar():
    """The mechanism must beat the support mean when neighbours carry signal."""
    torch.manual_seed(4211)
    generator = torch.Generator().manual_seed(4212)
    transport = SimilarityTransport(4, hidden_dim=8, use_learned_key=False)
    tasks, supports, queries = 256, 4, 6
    cluster_support = torch.randint(3, (tasks, supports), generator=generator)
    cluster_query = torch.randint(3, (tasks, queries), generator=generator)
    offset = torch.randn(tasks, 3, generator=generator)
    residual = torch.gather(offset, 1, cluster_support)
    target = torch.gather(offset, 1, cluster_query)
    similarity = (cluster_query[:, :, None] == cluster_support[:, None, :]).float()
    with torch.no_grad():
        value, _ = transport(
            torch.zeros(tasks, supports, 4), torch.zeros(tasks, queries, 4),
            residual, similarity)
    level = residual.mean(-1, keepdim=True).expand_as(target)
    assert (value - target).square().mean() < 0.5 * (level - target).square().mean()
