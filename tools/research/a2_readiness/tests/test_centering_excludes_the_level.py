"""Structural probes for the load-bearing property of Stage A3 (CPC).

The whole design rests on one algebraic claim: computing the protein
counterfactual on the **within-target centered** prediction makes the
level branch unable to satisfy it, because `protein_value(P)` is constant
across the queries of one target and cancels exactly under centering.

These probes verify that claim on the real `SimilarityGrammarModel`, with no
training and no dataset — small random tensors are enough, because the claim
is algebraic. They are research probes, not repository contracts, so they live
under `tools/research/`; they move to `tools/tests/` only if the family is
admitted.

Run: `conda run -n drug python -m pytest tools/research/a2_readiness/tests -q`
"""
from __future__ import annotations

from pathlib import Path
import sys

import pytest
import torch

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM  # noqa: E402
from model.similarity_grammar import SimilarityGrammarModel      # noqa: E402

PROTEIN_DIM, SLOTS, ATOMS, QUERIES = 32, 12, 7, 6


def build(seed: int = 0) -> SimilarityGrammarModel:
    torch.manual_seed(seed)
    return SimilarityGrammarModel(
        protein_dim=PROTEIN_DIM, hidden_dim=24, task_dim=12, ligand_layers=2,
        pair_dim=24, pair_latents=6, pair_heads=2,
        use_learned_key=False).double()


def protein(seed: int):
    generator = torch.Generator().manual_seed(seed)
    pooled = torch.randn(1, PROTEIN_DIM, generator=generator, dtype=torch.float64)
    tokens = torch.randn(1, SLOTS, PROTEIN_DIM, generator=generator,
                         dtype=torch.float64)
    mask = torch.ones(1, SLOTS, dtype=torch.float64)
    chemistry = torch.rand(1, SLOTS, 4, generator=generator, dtype=torch.float64)
    return pooled, tokens, mask, chemistry


def ligands(count: int, seed: int):
    generator = torch.Generator().manual_seed(seed)
    atoms = torch.rand(1, count, ATOMS, ATOM_FEAT_DIM, generator=generator,
                       dtype=torch.float64)
    bonds = torch.rand(1, count, ATOMS, ATOMS, BOND_FEAT_DIM,
                       generator=generator, dtype=torch.float64)
    bonds = bonds * (bonds > 0.7)
    bonds = 0.5 * (bonds + bonds.transpose(2, 3))
    mask = torch.ones(1, count, ATOMS, dtype=torch.float64)
    fingerprint = (torch.rand(1, count, 64, generator=generator,
                              dtype=torch.float64) > 0.7).double()
    return atoms, bonds, mask, fingerprint


def zero_shot(model, protein_parts, ligand_parts) -> torch.Tensor:
    pooled, tokens, mask, chemistry = protein_parts
    atoms, bonds, atom_mask, fingerprint = ligand_parts
    empty_atoms = atoms[:, :0]
    output = model(
        pooled, tokens, mask,
        empty_atoms, bonds[:, :0], atom_mask[:, :0],
        torch.zeros(1, 0, dtype=torch.float64),
        atoms, bonds, atom_mask, adapt=False,
        protein_chemistry=chemistry,
        support_fingerprint=fingerprint[:, :0],
        query_fingerprint=fingerprint)
    return output.zero_shot.squeeze(0)


def centered(values: torch.Tensor) -> torch.Tensor:
    return values - values.mean()


def cpc_loss(model, correct_protein, donor_protein, ligand_parts,
             labels: torch.Tensor, centered_form: bool) -> torch.Tensor:
    """The Stage A3 term; `centered_form=False` is A0's existing uncentered one."""
    p = zero_shot(model, correct_protein, ligand_parts)
    q = zero_shot(model, donor_protein, ligand_parts)
    if centered_form:
        p, q, target = centered(p), centered(q), centered(labels)
    else:
        target = labels
    correct = (p - target).square().mean()
    wrong = (q - target).square().mean()
    stacked = torch.stack((correct, wrong)) / 0.1
    return torch.log_softmax(-stacked, 0)[0].neg()


def test_centering_removes_the_protein_level_exactly():
    """`protein_value` is constant within a target, so it cancels under centering."""
    model = build()
    parts = ligands(QUERIES, seed=3)
    left = zero_shot(model, protein(1), parts)
    right = zero_shot(model, protein(2), parts)
    # Two proteins differ; the additive protein branch is a per-target constant,
    # so any difference it contributes is removed by centering.
    assert not torch.allclose(left, right)
    difference = centered(left) - centered(right)
    assert difference.abs().max() > 0, "the trunk has no protein-conditioned shape at all"


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_centered_term_gives_the_level_head_exactly_zero_gradient(seed):
    """The load-bearing property of Stage A3.

    `protein_head` is the dedicated level branch. Under the centered
    counterfactual it must receive *identically* zero gradient, so the level
    cannot satisfy the term and the gradient must reach the interaction path.
    """
    model = build(seed)
    parts = ligands(QUERIES, seed=10 + seed)
    labels = torch.randn(QUERIES, generator=torch.Generator().manual_seed(seed),
                         dtype=torch.float64)
    model.zero_grad(set_to_none=True)
    cpc_loss(model, protein(1), protein(2), parts, labels,
             centered_form=True).backward()
    for name, parameter in model.protein_head.named_parameters():
        assert parameter.grad is None or parameter.grad.abs().max() < 1e-12, (
            f"centered CPC leaked gradient into protein_head.{name}")


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_uncentered_term_does_reach_the_level_head(seed):
    """The contrast: A0's existing form is satisfiable by the level branch.

    This is why the incumbent's protein counterfactual never reached the
    ordering — see DATAFLOW_AUDIT.md F6.
    """
    model = build(seed)
    parts = ligands(QUERIES, seed=10 + seed)
    labels = torch.randn(QUERIES, generator=torch.Generator().manual_seed(seed),
                         dtype=torch.float64)
    model.zero_grad(set_to_none=True)
    cpc_loss(model, protein(1), protein(2), parts, labels,
             centered_form=False).backward()
    total = max(float(p.grad.abs().max()) for p in model.protein_head.parameters()
                if p.grad is not None)
    assert total > 1e-9, "uncentered CPC unexpectedly missed the level head"


def test_centered_term_reaches_the_interaction_path():
    """Zero gradient everywhere would make the term vacuous, not selective."""
    model = build(4)
    parts = ligands(QUERIES, seed=21)
    labels = torch.randn(QUERIES, generator=torch.Generator().manual_seed(4),
                         dtype=torch.float64)
    model.zero_grad(set_to_none=True)
    cpc_loss(model, protein(1), protein(2), parts, labels,
             centered_form=True).backward()
    reached = max(float(p.grad.abs().max())
                  for p in model.interaction_head.parameters()
                  if p.grad is not None)
    assert reached > 1e-9, "centered CPC delivers no gradient to the interaction head"


def test_centered_term_is_invariant_to_a_per_target_constant():
    """Adding any constant to every query must not change the term."""
    generator = torch.Generator().manual_seed(8)
    prediction = torch.randn(QUERIES, generator=generator, dtype=torch.float64)
    donor = torch.randn(QUERIES, generator=generator, dtype=torch.float64)
    labels = torch.randn(QUERIES, generator=generator, dtype=torch.float64)

    def term(offset: float) -> float:
        p, q = centered(prediction + offset), centered(donor + offset)
        target = centered(labels)
        correct = (p - target).square().mean()
        wrong = (q - target).square().mean()
        return float(correct - wrong)

    assert abs(term(0.0) - term(3.7)) < 1e-12


def test_centered_term_is_query_permutation_invariant():
    generator = torch.Generator().manual_seed(9)
    prediction = torch.randn(QUERIES, generator=generator, dtype=torch.float64)
    donor = torch.randn(QUERIES, generator=generator, dtype=torch.float64)
    labels = torch.randn(QUERIES, generator=generator, dtype=torch.float64)
    order = torch.randperm(QUERIES, generator=generator)

    def term(p, q, y) -> float:
        p, q, y = centered(p), centered(q), centered(y)
        return float((p - y).square().mean() - (q - y).square().mean())

    assert abs(term(prediction, donor, labels)
               - term(prediction[order], donor[order], labels[order])) < 1e-12


def test_single_query_panel_degenerates_to_zero_not_nan():
    """Centering a one-query panel gives 0; the term must be finite and inert."""
    one = torch.tensor([1.5], dtype=torch.float64)
    p, q, y = centered(one), centered(one * 2), centered(one * 3)
    value = (p - y).square().mean() - (q - y).square().mean()
    assert torch.isfinite(value) and float(value.abs()) < 1e-12
