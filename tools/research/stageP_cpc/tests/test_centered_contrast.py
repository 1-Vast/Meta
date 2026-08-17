"""Structural contract for Stage P's single training change.

`centered_protein_contrast` must have exactly one load-bearing property: the
additive protein level branch cannot satisfy it. `protein_value(P)` is constant
across the queries of one target, so subtracting the query mean removes it
identically and `d(loss)/d(protein_head)` is zero at every parameter value.

Without that, Stage P would repeat the incumbent's failure mode — the level
branch answers "which protein is this?" with a 0.215 pK shift and the gradient
never reaches the ligand-varying path (DATAFLOW_AUDIT F6).

These probes run on the loss function directly and on the real
`SimilarityGrammarModel`, in float64, with no dataset.
"""
from __future__ import annotations

from pathlib import Path
import sys

import pytest
import torch

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM     # noqa: E402
from model.similarity_grammar import SimilarityGrammarModel         # noqa: E402
from scripts.train_qpsmp import (                                   # noqa: E402
    TrainConfig, binding_contrastive_loss, centered_protein_contrast,
)

PROTEIN_DIM, SLOTS, ATOMS, QUERIES = 32, 12, 7, 6


def panel(seed: int):
    generator = torch.Generator().manual_seed(seed)
    return (torch.randn(QUERIES, generator=generator, dtype=torch.float64),
            torch.randn(QUERIES, generator=generator, dtype=torch.float64),
            torch.randn(QUERIES, generator=generator, dtype=torch.float64))


# --- algebra ---------------------------------------------------------------

@pytest.mark.parametrize("seed", [0, 1, 2])
def test_a_per_target_constant_on_either_side_changes_nothing(seed):
    """The defining property: level shifts are invisible to the term."""
    correct, wrong, truth = panel(seed)
    base = float(centered_protein_contrast(correct, wrong, truth, 0.1))
    for shift in (3.7, -2.1):
        assert float(centered_protein_contrast(
            correct + shift, wrong, truth, 0.1)) == pytest.approx(base, abs=1e-12)
        assert float(centered_protein_contrast(
            correct, wrong + shift, truth, 0.1)) == pytest.approx(base, abs=1e-12)
        assert float(centered_protein_contrast(
            correct, wrong, truth + shift, 0.1)) == pytest.approx(base, abs=1e-12)


def test_the_incumbent_form_is_not_invariant_to_a_level_shift():
    """The contrast that motivates the change."""
    correct, wrong, truth = panel(3)

    def uncentered(c, w):
        return float(binding_contrastive_loss(
            [(c - truth).square().mean(), (w - truth).square().mean()], 0.1))

    assert abs(uncentered(correct, wrong)
               - uncentered(correct + 3.7, wrong)) > 1e-3


def test_it_is_query_permutation_invariant():
    correct, wrong, truth = panel(4)
    order = torch.randperm(QUERIES, generator=torch.Generator().manual_seed(4))
    assert float(centered_protein_contrast(correct, wrong, truth, 0.1)) == \
        pytest.approx(float(centered_protein_contrast(
            correct[order], wrong[order], truth[order], 0.1)), abs=1e-12)


def test_a_one_query_panel_is_inert_rather_than_nan():
    one = torch.tensor([1.5], dtype=torch.float64)
    value = centered_protein_contrast(one, one * 2, one * 3, 0.1)
    assert torch.isfinite(value) and float(value) == 0.0


def test_it_is_minimised_by_a_better_correct_arm_not_only_a_worse_donor():
    """Both routes lower the loss; the gates, not the algebra, forbid one."""
    _, wrong, truth = panel(5)
    worse = float(centered_protein_contrast(truth * 0.0, wrong, truth, 0.1))
    better = float(centered_protein_contrast(truth, wrong, truth, 0.1))
    assert better < worse
    # Documented, so the analysis cannot forget it: degrading the donor also
    # lowers the loss, which is why Stage P reports both sides separately.
    ruined = float(centered_protein_contrast(truth, wrong * 50.0, truth, 0.1))
    assert ruined < worse


# --- on the real model -----------------------------------------------------

def build(seed: int = 0) -> SimilarityGrammarModel:
    torch.manual_seed(seed)
    return SimilarityGrammarModel(
        protein_dim=PROTEIN_DIM, hidden_dim=24, task_dim=12, ligand_layers=2,
        pair_dim=24, pair_latents=6, pair_heads=2,
        use_learned_key=False).double()


def protein(seed: int):
    generator = torch.Generator().manual_seed(seed)
    return (torch.randn(1, PROTEIN_DIM, generator=generator, dtype=torch.float64),
            torch.randn(1, SLOTS, PROTEIN_DIM, generator=generator, dtype=torch.float64),
            torch.ones(1, SLOTS, dtype=torch.float64),
            torch.rand(1, SLOTS, 4, generator=generator, dtype=torch.float64))


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


def zero_shot(model, protein_parts, ligand_parts):
    pooled, tokens, mask, chemistry = protein_parts
    atoms, bonds, atom_mask, fingerprint = ligand_parts
    return model(pooled, tokens, mask,
                 atoms[:, :0], bonds[:, :0], atom_mask[:, :0],
                 torch.zeros(1, 0, dtype=torch.float64),
                 atoms, bonds, atom_mask, adapt=False,
                 protein_chemistry=chemistry,
                 support_fingerprint=fingerprint[:, :0],
                 query_fingerprint=fingerprint).zero_shot.squeeze(0)


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_the_level_head_receives_exactly_zero_gradient(seed):
    """The property Stage P depends on, on the real model class."""
    model = build(seed)
    parts = ligands(QUERIES, seed=10 + seed)
    truth = torch.randn(QUERIES, generator=torch.Generator().manual_seed(seed),
                        dtype=torch.float64)
    model.zero_grad(set_to_none=True)
    centered_protein_contrast(
        zero_shot(model, protein(1), parts),
        zero_shot(model, protein(2), parts), truth, 0.1).backward()
    for name, parameter in model.protein_head.named_parameters():
        assert parameter.grad is None or parameter.grad.abs().max() < 1e-12, (
            f"the centered contrast leaked gradient into protein_head.{name}")


def test_it_does_reach_the_interaction_path():
    """Zero gradient everywhere would make the term vacuous, not selective."""
    model = build(4)
    parts = ligands(QUERIES, seed=21)
    truth = torch.randn(QUERIES, generator=torch.Generator().manual_seed(4),
                        dtype=torch.float64)
    model.zero_grad(set_to_none=True)
    centered_protein_contrast(
        zero_shot(model, protein(1), parts),
        zero_shot(model, protein(2), parts), truth, 0.1).backward()
    reached = max(float(p.grad.abs().max())
                  for p in model.interaction_head.parameters()
                  if p.grad is not None)
    assert reached > 1e-9


# --- the configuration switch ---------------------------------------------

def test_the_default_configuration_is_the_incumbent():
    """Every recorded arm must be unchanged by the existence of this flag."""
    assert TrainConfig.protein_contrast_form == "uncentered"


def test_the_trainer_exposes_exactly_the_two_forms():
    import subprocess
    finished = subprocess.run(
        [sys.executable, "-m", "scripts.train_qpsmp", "--help"],
        cwd=str(ROOT), capture_output=True, text=True, timeout=300)
    assert "--protein-contrast-form {uncentered,centered}" in finished.stdout


# --- the matched-control guard --------------------------------------------

def test_the_evaluator_refuses_an_unmatched_control(tmp_path):
    """A control that differs in anything but the intended change is not one.

    The guard runs before any metric is computed, so an accidentally mismatched
    budget, seed, architecture or split cannot produce an interpretable-looking
    contrast.
    """
    import json
    from tools.research.stageP_cpc.evaluate import (
        SEEDS, verify_arms_are_matched)

    def write(stage, seed, arm, **overrides):
        config = {"seed": seed, "steps": 1200, "arch": "similarity_only",
                  "split_directory": "x", "protein_contrast_form": "uncentered",
                  "protein_contrast_loss_weight": 0.5}
        config.update(overrides)
        path = stage / f"{arm}_seed{seed}"
        path.mkdir(parents=True, exist_ok=True)
        (path / "RESULT.json").write_text(json.dumps({"config": config}),
                                          encoding="utf-8")

    good = tmp_path / "good"
    for seed in SEEDS:
        write(good, seed, "A0repro")
        write(good, seed, "CPCoverdrive",
              protein_contrast_form="centered", protein_contrast_loss_weight=2.0)
    report = verify_arms_are_matched(good)
    assert set(report) == {str(s) for s in SEEDS}

    # An unintended difference must raise, not warn.
    bad = tmp_path / "bad"
    for seed in SEEDS:
        write(bad, seed, "A0repro")
        write(bad, seed, "CPCoverdrive", steps=600,
              protein_contrast_form="centered", protein_contrast_loss_weight=2.0)
    with pytest.raises(ValueError, match="not matched"):
        verify_arms_are_matched(bad)

    # A flag that silently failed to take effect must also raise.
    inert = tmp_path / "inert"
    for seed in SEEDS:
        write(inert, seed, "A0repro")
        write(inert, seed, "CPCoverdrive", protein_contrast_loss_weight=2.0)
    with pytest.raises(ValueError):
        verify_arms_are_matched(inert)


def test_the_guard_requires_a_complete_seed_pair(tmp_path):
    import json
    from tools.research.stageP_cpc.evaluate import verify_arms_are_matched
    stage = tmp_path / "partial"
    path = stage / "A0repro_seed20260815"
    path.mkdir(parents=True)
    (path / "RESULT.json").write_text(json.dumps({"config": {}}), encoding="utf-8")
    with pytest.raises(ValueError, match="no complete seed pair"):
        verify_arms_are_matched(stage)
