"""Regression tests for the S5D estimand and collapse diagnostics.

S5D trains nothing. These tests pin its registration, its two statistics and
the property that makes the conditional estimand worth measuring: pocket
membership cancels exactly inside the symmetric difference.
"""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "research" / "s7_l2b_r0r" / "s5d_diagnostics.py"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


s5d = _load("s5d_diagnostics_test", MODULE)


# ------------------------------------------------------------ registration
def test_registration_hash_matches_the_committed_document():
    document = (ROOT / "research" / "s7_l2b_r0r" /
                "PREREG_PHASE2B_S5D_ESTIMAND_AND_COLLAPSE_DIAGNOSTICS.md")
    assert s5d.sha_file(document) == s5d.PREREG_SHA


def test_frozen_constants_match_the_preregistration():
    assert s5d.MIN_PAIRS_FOR_RHO == 3
    assert s5d.RHO_COLLAPSE_MEDIAN == 0.80
    assert s5d.RHO_EXCESS_OVER_DATA == 0.10
    assert (s5d.E1_MARGIN, s5d.E2_MARGIN, s5d.E3_MARGIN) == (0.05, 0.03, 0.03)
    assert s5d.ARMS == ("candidate", "baseline41", "foreign", "permuted")
    assert s5d.SEED_BOOT == 20260903


def test_stage_introduces_no_trainable_parameter():
    source = MODULE.read_text(encoding="utf-8")
    for forbidden in ("optimizer", "backward", "requires_grad", "nn.Parameter",
                      "Adam"):
        assert forbidden not in source


# ------------------------------------------------------------ D1 statistic
def test_top_principal_energy_fraction_is_exactly_one_for_identical_fields():
    rng = np.random.default_rng(0)
    identical = np.repeat(rng.normal(size=(1, 40)), 12, axis=0)
    assert s5d.top_principal_energy_fraction(identical) == 1.0


def test_a_fully_collapsed_construct_clears_the_registered_collapse_threshold():
    """Rescaled copies of one direction normalize to within ~1e-16 of each
    other, so the centred spectrum is float noise rather than an exact zero and
    the statistic floors below 1.0. What the registered rule needs is only that
    total collapse still reads as collapse."""
    rng = np.random.default_rng(0)
    direction = rng.normal(size=(1, 40))
    scaled = rng.uniform(0.5, 4.0, size=(12, 1)) * direction
    assert s5d.top_principal_energy_fraction(scaled) > s5d.RHO_COLLAPSE_MEDIAN


def test_the_float_floor_only_understates_collapse_never_manufactures_it():
    """The noise floor biases rho downward in the degenerate limit, so it can
    hide collapse but cannot invent it. Diverse fields stay far below the
    threshold."""
    rng = np.random.default_rng(7)
    diverse = rng.normal(size=(12, 40))
    assert s5d.top_principal_energy_fraction(diverse) < s5d.RHO_COLLAPSE_MEDIAN


def test_top_principal_energy_fraction_is_small_for_isotropic_fields():
    rng = np.random.default_rng(1)
    assert s5d.top_principal_energy_fraction(rng.normal(size=(400, 40))) < 0.15


def test_top_principal_energy_fraction_ignores_sign_and_scale_only():
    """Unit normalization removes magnitude, so a sign flip is a genuine second
    direction rather than a rescaling."""
    rng = np.random.default_rng(2)
    direction = rng.normal(size=(1, 30))
    signs = np.array([[1.0]] * 6 + [[-1.0]] * 6)
    assert s5d.top_principal_energy_fraction(signs * direction) == pytest.approx(
        1.0, abs=1e-9)


def test_rho_is_undefined_below_the_registered_pair_floor():
    rng = np.random.default_rng(3)
    assert np.isnan(s5d.top_principal_energy_fraction(rng.normal(size=(2, 10))))


# ------------------------------------------------------------ D2 estimand
def test_conditional_estimand_ignores_unchanged_residues_entirely():
    """The whole point: residues outside the symmetric difference cannot move
    the statistic, so a pocket-membership field cancels exactly."""
    score = np.array([9.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    gain, loss = {0}, {1}
    reference = s5d.conditional_metrics(score, gain, loss, 6)
    shifted = score.copy()
    shifted[2:] = 100.0                      # unchanged residues dominate
    assert s5d.conditional_metrics(shifted, gain, loss, 6)["ap_cond"] == \
        pytest.approx(reference["ap_cond"])


def test_conditional_estimand_rewards_the_correct_direction():
    correct = np.array([5.0, -5.0, 0.0, 0.0])
    wrong = -correct
    gain, loss = {0}, {1}
    assert s5d.conditional_metrics(correct, gain, loss, 4)["ap_cond"] == 1.0
    assert s5d.conditional_metrics(wrong, gain, loss, 4)["ap_cond"] < \
        s5d.conditional_metrics(correct, gain, loss, 4)["ap_cond"]


def test_conditional_chance_is_the_constant_score_value():
    rng = np.random.default_rng(4)
    gain = set(range(0, 6))
    loss = set(range(6, 10))
    metrics = s5d.conditional_metrics(rng.normal(size=20), gain, loss, 20)
    constant = s5d.conditional_metrics(np.zeros(20), gain, loss, 20)
    assert metrics["chance_cond"] == pytest.approx(constant["ap_cond"])
    assert metrics["changed"] == 10
    assert metrics["gain_fraction"] == pytest.approx(0.6)


def test_conditional_estimand_is_ineligible_without_both_classes():
    assert s5d.conditional_metrics(np.zeros(5), {0, 1}, set(), 5) is None
    assert s5d.conditional_metrics(np.zeros(5), set(), {2, 3}, 5) is None
    assert s5d.conditional_metrics(np.zeros(5), {0}, set(), 5) is None


def test_conditional_estimand_is_invariant_to_a_constant_pocket_offset():
    rng = np.random.default_rng(5)
    score = rng.normal(size=30)
    gain = {1, 4, 7}
    loss = {2, 9, 11}
    base = s5d.conditional_metrics(score, gain, loss, 30)["ap_cond"]
    offset = score + 17.0
    assert s5d.conditional_metrics(offset, gain, loss, 30)["ap_cond"] == \
        pytest.approx(base)


# ------------------------------------------------------------ recorded result
def test_recorded_verdict_does_not_reopen_the_s4r_route():
    path = ROOT / "report" / "s7_l2b_r0r" / "PHASE2B_S5D_GATE.json"
    if not path.is_file():
        pytest.skip("S5D has not been executed in this checkout")
    gate = json.loads(path.read_text(encoding="utf-8"))
    assert gate["trainable_parameters_introduced"] == 0
    assert gate["affinity_value_reads"] == 0
    assert gate["heldoutB_status"] == "NOT_CREATED_AND_NOT_READ"
    assert gate["s4r_verdict_unchanged"] == "REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED"
    assert gate["TERMINAL_VERDICT"] in {
        "S5D_CONTRACT_FAIL_CLOSED",
        "LIGAND_DIRECTION_COLLAPSE_NOT_CONFIRMED",
        "POSE_FREE_LIGAND_CONDITIONED_DIRECTION_ABSENT_UNDER_CONDITIONAL_ESTIMAND",
        "CONDITIONAL_ESTIMAND_RECOVERS_LIGAND_SPECIFIC_DIRECTION_IN_DEVELOPMENT"}
    if gate["TERMINAL_VERDICT"] != \
            "CONDITIONAL_ESTIMAND_RECOVERS_LIGAND_SPECIFIC_DIRECTION_IN_DEVELOPMENT":
        assert gate["authorized_next_action"].startswith("none")
