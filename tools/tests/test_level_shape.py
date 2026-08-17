"""Structural gates for Core Innovation A, the level-shape factorized model.

These run before any training. Each one falsifies a specific way the
factorization could be cosmetic rather than real.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from model.level_shape import CHANNELS, LevelShapeModel

PROTEIN_DIM, SLOTS, ATOMS = 64, 16, 9


def build(seed: int = 0, **kwargs) -> LevelShapeModel:
    torch.manual_seed(seed)
    return LevelShapeModel(protein_dim=PROTEIN_DIM, hidden_dim=32, task_dim=16,
                           ligand_layers=2, pair_dim=24, pair_heads=4,
                           anchors=6, **kwargs).double()


def ligand(count: int, seed: int) -> tuple[torch.Tensor, ...]:
    generator = torch.Generator().manual_seed(seed)
    atoms = torch.rand(count, ATOMS, ATOM_FEAT_DIM,
                       generator=generator, dtype=torch.float64)
    bonds = torch.rand(count, ATOMS, ATOMS, BOND_FEAT_DIM,
                       generator=generator, dtype=torch.float64)
    bonds = bonds * (bonds > 0.7)
    bonds = 0.5 * (bonds + bonds.transpose(1, 2))
    mask = torch.ones(count, ATOMS, dtype=torch.float64)
    fingerprint = (torch.rand(count, 64, generator=generator,
                              dtype=torch.float64) > 0.7).double()
    fingerprint[:, 0] = 1.0                       # never an all-zero row
    return atoms, bonds, mask, fingerprint


def protein(seed: int) -> tuple[torch.Tensor, ...]:
    generator = torch.Generator().manual_seed(seed)
    pooled = torch.randn(PROTEIN_DIM, generator=generator, dtype=torch.float64)
    tokens = torch.randn(SLOTS, PROTEIN_DIM, generator=generator, dtype=torch.float64)
    mask = torch.ones(SLOTS, dtype=torch.float64)
    chemistry = torch.rand(SLOTS, 4, generator=generator, dtype=torch.float64)
    return pooled, tokens, mask, chemistry


def episode(support: int, query: int, *, protein_seed: int = 1,
            ligand_seed: int = 2) -> dict:
    pooled, tokens, mask, chemistry = protein(protein_seed)
    sa, sb, sm, sf = ligand(support, ligand_seed)
    qa, qb, qm, qf = ligand(query, ligand_seed + 100)
    return {
        "protein_pooled": pooled, "protein_tokens": tokens, "protein_mask": mask,
        "protein_chemistry": chemistry,
        "support_atoms": sa, "support_bonds": sb, "support_mask": sm,
        "support_y": torch.linspace(4.0, 9.0, support, dtype=torch.float64),
        "query_atoms": qa, "query_bonds": qb, "query_mask": qm,
        "support_fingerprint": sf, "query_fingerprint": qf,
    }


# ---------------------------------------------------------------- factorization


def test_centered_branch_has_no_constant_component():
    """mean over the anchors of centered_interaction is exactly 0, per protein.

    This is the gate that makes the factorization real: the centered branch
    cannot express a target-level constant, so the level must live in
    `target_level`.
    """
    model = build()
    pooled, tokens, mask, chemistry = protein(7)
    residues, _ = model.encode_protein(
        pooled[None], tokens[None], mask[None], chemistry[None])
    anchors = model.anchor.unsqueeze(0)
    centered, _ = model.centered_interaction(anchors, residues, mask[None])
    assert centered.mean(-1).abs().max() < 1e-10


def test_centered_branch_is_protein_specific():
    """Two different proteins give different anchor centerings and outputs."""
    model = build()
    left, right = episode(0, 5, protein_seed=1), episode(0, 5, protein_seed=9)
    a = model(**left, adapt=False)
    b = model(**{**right, "query_atoms": left["query_atoms"],
                 "query_bonds": left["query_bonds"],
                 "query_mask": left["query_mask"],
                 "query_fingerprint": left["query_fingerprint"]}, adapt=False)
    assert not torch.allclose(a.centered, b.centered, atol=1e-6)
    assert (a.anchor_centering - b.anchor_centering).abs().item() > 1e-8


def test_level_branch_cannot_change_within_target_ordering():
    """target_level is one number per protein, broadcast over every query."""
    output = build()(**episode(0, 6), adapt=False)
    assert output.target_level.std().item() == pytest.approx(0.0, abs=1e-12)


def test_ligand_prior_is_protein_blind():
    model = build()
    base = episode(0, 5, protein_seed=1)
    other = protein(9)
    swapped = model(**{**base, "protein_pooled": other[0],
                       "protein_tokens": other[1], "protein_mask": other[2],
                       "protein_chemistry": other[3]}, adapt=False)
    assert torch.allclose(model(**base, adapt=False).ligand_prior,
                          swapped.ligand_prior, atol=1e-12)


def test_endpoint_is_the_exact_sum_of_three_branches():
    """One scalar endpoint, no duplicate affinity path."""
    output = build()(**episode(0, 5), adapt=False)
    assert torch.allclose(
        output.endpoint,
        output.ligand_prior + output.target_level + output.centered, atol=1e-12)


# ---------------------------------------------------------------- leakage


def test_prediction_is_independent_of_the_other_queries():
    """Inductive: no query-panel mean or other transductive statistic."""
    model = build()
    full = episode(0, 6)
    joint = model(**full, adapt=False).endpoint
    for index in range(6):
        single = model(**{**full,
                          "query_atoms": full["query_atoms"][index:index + 1],
                          "query_bonds": full["query_bonds"][index:index + 1],
                          "query_mask": full["query_mask"][index:index + 1],
                          "query_fingerprint":
                              full["query_fingerprint"][index:index + 1]},
                       adapt=False).endpoint
        assert torch.allclose(joint[index], single[0], atol=1e-10)


def test_zero_support_returns_the_endpoint_exactly():
    output = build()(**episode(0, 5))
    assert torch.equal(output.prediction, output.endpoint)
    assert output.transport.abs().max().item() == 0.0


def test_query_labels_are_never_an_input():
    parameters = build().forward.__doc__ or ""
    del parameters
    import inspect
    signature = inspect.signature(LevelShapeModel.forward)
    assert "query_y" not in signature.parameters


def test_support_permutation_invariance():
    model = build()
    base = episode(4, 5)
    order = torch.tensor([2, 0, 3, 1])
    permuted = {**base,
                "support_atoms": base["support_atoms"][order],
                "support_bonds": base["support_bonds"][order],
                "support_mask": base["support_mask"][order],
                "support_y": base["support_y"][order],
                "support_fingerprint": base["support_fingerprint"][order]}
    assert torch.allclose(model(**base).prediction,
                          model(**permuted).prediction, atol=1e-10)


def test_query_permutation_equivariance():
    model = build()
    base = episode(3, 5)
    order = torch.tensor([4, 0, 2, 1, 3])
    permuted = {**base,
                "query_atoms": base["query_atoms"][order],
                "query_bonds": base["query_bonds"][order],
                "query_mask": base["query_mask"][order],
                "query_fingerprint": base["query_fingerprint"][order]}
    assert torch.allclose(model(**base).prediction[order],
                          model(**permuted).prediction, atol=1e-10)


def test_support_labels_enter_only_as_residuals():
    """Shifting every support label by c shifts the transport by shrink * c."""
    model = build()
    base = episode(3, 5)
    shifted = {**base, "support_y": base["support_y"] + 2.0}
    delta = (model(**shifted).transport - model(**base).transport)
    shrink = float(model.transport.shrinkage(3, torch.zeros(1, dtype=torch.float64)))
    assert torch.allclose(delta, torch.full_like(delta, 2.0 * shrink), atol=1e-9)


def test_geometry_input_is_refused():
    with pytest.raises(ValueError, match="common-frame"):
        build()(**episode(1, 2),
                geometry_available=torch.ones(1, dtype=torch.bool))


# ---------------------------------------------------------------- gradients


@pytest.mark.parametrize("support", [0, 1, 2, 5])
def test_every_trainable_parameter_receives_gradient(support):
    model = build()
    output = model(**episode(support, 6))
    output.prediction.sum().backward()
    missing = [name for name, parameter in model.named_parameters()
               if parameter.requires_grad
               and (parameter.grad is None or not parameter.grad.abs().sum())]
    if support == 0:
        # Inactive at k=0 by contract, and exercised by the k>=1 cases.
        missing = [n for n in missing if not n.startswith("transport.")]
    if support == 1:
        # A softmax over one support is identically 1, so the kernel scale has
        # no gradient at k=1. That degeneracy is a documented property of the
        # fixed Tanimoto kernel, not a dead branch: the shrinkage still trains
        # here, and the scale trains at k>=2.
        assert "transport.log_shrinkage" not in missing
        missing = [n for n in missing if n != "transport.similarity_scale"]
    assert missing == []


def test_level_and_centered_gradients_are_separable():
    """Detaching one branch removes gradient from its head and no other."""
    model = build()
    output = model(**episode(0, 5), adapt=False)
    (output.ligand_prior + output.centered).sum().backward()
    assert all(p.grad is None or not p.grad.abs().sum()
               for _, p in model.level_head.named_parameters())
    assert any(p.grad is not None and p.grad.abs().sum()
               for _, p in model.interaction.named_parameters())


def test_anchor_parameters_are_trained():
    model = build()
    model(**episode(0, 4), adapt=False).endpoint.sum().backward()
    assert model.anchor.grad is not None
    assert model.anchor.grad.abs().sum().item() > 0


def test_channel_count_matches_the_anchor_space():
    model = build()
    assert model.anchor.shape[1] == CHANNELS
    batch = episode(0, 3)
    channels = model.encode_ligand(*[batch[k].unsqueeze(0) for k in
                                     ("query_atoms", "query_bonds", "query_mask")])
    assert channels.shape == (1, 3, CHANNELS, model.hidden_dim)
