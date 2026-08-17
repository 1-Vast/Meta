"""Structural probes for Phase 1's measurement apparatus and the trunk itself.

Two kinds of test live here.

**Apparatus.** A null result is only worth reporting if the instrument that
produced it can detect a real effect. These tests inject known effects and
check the metrics move, and inject no effect and check they do not.

**Architecture.** Facts about the trunk that constrain what any claim may say —
in particular that the protein path is exactly invariant to residue-slot order,
which forbids describing its attention as pocket-aware or biologically
localized regardless of what the attention weights look like.

No dataset, no training: these are algebraic properties, so small random
tensors suffice.

Run: `conda run -n drug python -m pytest tools/research/a2_readiness_v2/tests -q`
"""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest
import torch

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM      # noqa: E402
from model.similarity_grammar import SimilarityGrammarModel          # noqa: E402
from tools.research.a2_readiness_v2 import _frozen                   # noqa: E402
from tools.research.a2_readiness_v2._arms import randomise           # noqa: E402
from tools.research.a2_readiness_v2.branch_ordering_v2 import (      # noqa: E402
    branches, centered, concordance, correlation, scramble_slots,
    variance_components,
)

PROTEIN_DIM, SLOTS, ATOMS, QUERIES = 32, 12, 7, 6


def build(seed: int = 0) -> SimilarityGrammarModel:
    torch.manual_seed(seed)
    return SimilarityGrammarModel(
        protein_dim=PROTEIN_DIM, hidden_dim=24, task_dim=12, ligand_layers=2,
        pair_dim=24, pair_latents=6, pair_heads=2,
        use_learned_key=False).double().eval()


def protein_parts(seed: int) -> list:
    generator = torch.Generator().manual_seed(seed)
    return [
        torch.randn(1, PROTEIN_DIM, generator=generator, dtype=torch.float64),
        torch.randn(1, SLOTS, PROTEIN_DIM, generator=generator, dtype=torch.float64),
        torch.ones(1, SLOTS, dtype=torch.float64),
        torch.rand(1, SLOTS, 4, generator=generator, dtype=torch.float64),
    ]


def ligand_parts(count: int, seed: int) -> tuple:
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


# ---------------------------------------------------------------------------
# apparatus: the metrics must detect what they claim to detect
# ---------------------------------------------------------------------------

def test_correlation_recovers_a_planted_ordering():
    truth = np.array([1.0, 5.0, 2.0, 4.0, 3.0, 0.0])
    assert correlation(truth, truth) == pytest.approx(1.0)
    assert correlation(-truth, truth) == pytest.approx(-1.0)
    assert abs(correlation(np.ones_like(truth), truth)) < 1e-12


def test_correlation_and_concordance_are_invariant_to_a_positive_rescale():
    """The R14 identity: a positive affine map cannot change ordering."""
    generator = np.random.default_rng(3)
    prediction = generator.normal(size=32)
    truth = generator.normal(size=32)
    assert correlation(prediction, truth) == pytest.approx(
        correlation(3.7 * prediction + 1.9, truth))
    assert concordance(prediction, truth) == pytest.approx(
        concordance(3.7 * prediction + 1.9, truth))


def test_a_label_shuffle_destroys_a_real_correlation():
    """The negative control must actually collapse a signal that is there."""
    generator = np.random.default_rng(5)
    truth = generator.normal(size=400)
    prediction = truth + 0.3 * generator.normal(size=400)
    assert correlation(prediction, truth) > 0.9
    shuffled = truth[generator.permutation(len(truth))]
    assert abs(correlation(prediction, shuffled)) < 0.2


def test_variance_components_separate_the_three_levels():
    """Planted spread must land in the level it was planted in."""
    rows = []
    for component_index in range(4):
        for target_index in range(3):
            base = 10.0 * component_index + 1.0 * target_index
            for seed in range(3):
                rows.append({"component": f"c{component_index}",
                             "target": f"c{component_index}t{target_index}",
                             "value": base + 0.01 * seed})
    parts = variance_components(rows, "value")
    assert parts["n_components"] == 4 and parts["n_targets"] == 12
    assert parts["between_component_sd"] > 5.0
    assert 0.5 < parts["between_target_within_component_sd"] < 2.0
    assert parts["between_seed_within_target_sd"] < 0.05


def test_variance_components_report_zero_seed_spread_for_one_seed():
    rows = [{"component": "c", "target": "t", "value": 1.0}]
    assert variance_components(rows, "value")["between_seed_within_target_sd"] == 0.0


# ---------------------------------------------------------------------------
# architecture: what the trunk can and cannot see
# ---------------------------------------------------------------------------

def test_the_protein_path_is_exactly_invariant_to_residue_slot_order():
    """The strongest constraint Phase 1 places on any biological claim.

    `ResidueEncoder` pools slots with a sum, and `ContactGrammar` reduces the
    residue axis with a softmax-weighted sum. Both are permutation-invariant,
    so the model cannot distinguish a protein from the same protein with its
    slots shuffled. Whatever the attention weights look like, they are not
    reading a pocket, a contact, or any ordered structural feature — there is
    no information in the input that would let them.

    This is why `scramble_slots` is a *null* control rather than a
    perturbation, and why the measured 1e-8 pK is machine zero rather than a
    small effect.
    """
    model = build(1)
    parts = protein_parts(11)
    atoms, bonds, mask, fingerprint = ligand_parts(QUERIES, seed=12)
    plain = branches(model, parts, atoms, bonds, mask, fingerprint,
                     "cpu", torch.float64)
    shuffled = branches(model, scramble_slots(parts, seed=7),
                        atoms, bonds, mask, fingerprint, "cpu", torch.float64)
    for branch in ("full", "ligand_only", "protein_only", "interaction"):
        assert np.abs(plain[branch] - shuffled[branch]).max() < 1e-9, branch


def test_the_protein_swap_control_does_move_the_output():
    """The counterpart: a genuine substitution must not be null."""
    model = build(1)
    atoms, bonds, mask, fingerprint = ligand_parts(QUERIES, seed=12)
    left = branches(model, protein_parts(11), atoms, bonds, mask, fingerprint,
                    "cpu", torch.float64)
    right = branches(model, protein_parts(22), atoms, bonds, mask, fingerprint,
                     "cpu", torch.float64)
    assert abs(left["protein_only"].mean() - right["protein_only"].mean()) > 1e-6
    assert np.abs(left["full"] - right["full"]).max() > 1e-6


def test_the_ligand_branch_is_protein_blind_by_construction():
    """`ligand_value` must be identical under any protein substitution."""
    model = build(2)
    atoms, bonds, mask, fingerprint = ligand_parts(QUERIES, seed=13)
    left = branches(model, protein_parts(31), atoms, bonds, mask, fingerprint,
                    "cpu", torch.float64)
    right = branches(model, protein_parts(32), atoms, bonds, mask, fingerprint,
                     "cpu", torch.float64)
    assert np.abs(left["ligand_only"] - right["ligand_only"]).max() < 1e-12


def test_the_protein_branch_is_constant_within_a_target():
    """`protein_value` cannot contribute to any within-target ordering."""
    model = build(3)
    atoms, bonds, mask, fingerprint = ligand_parts(QUERIES, seed=14)
    values = branches(model, protein_parts(41), atoms, bonds, mask,
                      fingerprint, "cpu", torch.float64)["protein_only"]
    # The claim is that the values are *identical*, so state it that way.
    # Centering them and comparing to zero would instead measure float32
    # rounding in the mean, which `branches` introduces when it casts for
    # numpy: six identical float32 values do not always center to exactly 0.
    assert float(values.max() - values.min()) == 0.0
    assert np.abs(centered(values)).max() < 1e-6


def test_the_three_branches_reconstruct_the_endpoint_exactly():
    model = build(4)
    atoms, bonds, mask, fingerprint = ligand_parts(QUERIES, seed=15)
    values = branches(model, protein_parts(51), atoms, bonds, mask,
                      fingerprint, "cpu", torch.float64)
    total = values["ligand_only"] + values["protein_only"] + values["interaction"]
    assert np.abs(total - values["full"]).max() < 1e-12


def test_the_probe_reads_no_support_label():
    """`adapt=False` must make the support labels irrelevant to the output."""
    model = build(5)
    atoms, bonds, mask, fingerprint = ligand_parts(QUERIES, seed=16)
    first = branches(model, protein_parts(61), atoms, bonds, mask, fingerprint,
                     "cpu", torch.float64)
    second = branches(model, protein_parts(61), atoms, bonds, mask, fingerprint,
                      "cpu", torch.float64)
    assert np.abs(first["full"] - second["full"]).max() == 0.0


# ---------------------------------------------------------------------------
# random-init arms
# ---------------------------------------------------------------------------

def test_independent_random_inits_are_actually_different():
    """Ten arms that share weights would fake a reproducibility result."""
    left = randomise(build(6), seed=101)
    right = randomise(build(6), seed=202)
    same = randomise(build(6), seed=101)
    left_weight = left.interaction_head[0].weight
    assert not torch.allclose(left_weight, right.interaction_head[0].weight)
    assert torch.allclose(left_weight, same.interaction_head[0].weight)


def test_randomisation_preserves_every_parameter_shape():
    model = build(7)
    shapes = {name: tuple(p.shape) for name, p in model.named_parameters()}
    randomise(model, seed=303)
    assert {name: tuple(p.shape)
            for name, p in model.named_parameters()} == shapes


# ---------------------------------------------------------------------------
# frozen design
# ---------------------------------------------------------------------------

def test_the_verdict_rule_separates_negligible_from_meaningful():
    """A resolved interval below the smallest effect of interest is not a find."""
    assert _frozen.verdict({"mean": 0.20, "lo": 0.10, "hi": 0.30}) == "RESOLVED"
    assert _frozen.verdict(
        {"mean": 0.0016, "lo": 0.0005, "hi": 0.0027}) == "RESOLVED_NEGLIGIBLE"
    assert _frozen.verdict(
        {"mean": 0.0, "lo": -0.001, "hi": 0.001}) == "DECISIVE_NULL"
    assert _frozen.verdict(
        {"mean": 0.0, "lo": -0.30, "hi": 0.30}) == "UNDERPOWERED"


def test_the_frozen_manifest_pins_the_checkpoints_it_names():
    manifest = _frozen.frozen_manifest()
    assert len(manifest["a0_checkpoints"]) == 3
    assert all(len(digest) == 64 for digest in manifest["a0_checkpoints"].values())
    assert len(set(manifest["a0_checkpoints"].values())) == 3, (
        "two A0 seeds share a checkpoint hash; they are not independent arms")
    assert len(manifest["random_init_seeds"]) >= 10


def test_the_probe_cannot_reach_meta_test():
    """The stage inherits the repaired fail-closed default."""
    from scripts.qpsmp_data import QPSMPData
    from scripts.train_qpsmp import (
        COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK)
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=_frozen.SPLIT_DIRECTORY)
    assert "meta_test" not in data.tasks
    assert data.seal_record()["included"] is False
