"""Stage V structural and gate-verdict tests.

These tests do not train anything. They pin the two load-bearing facts of the
V0/V1 verdict: the core-inclusive exact key carries the shared core and the
frozen gates record exactly what the verified rerun produced — V0 fails on
target/component domination, the repeated-key internal surface is not
evaluable, and V1's interaction mean square sits below supervision noise.
"""
from __future__ import annotations

import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[4]
if str(REPO) not in sys.path:  # pragma: no cover
    sys.path.insert(0, str(REPO))

from tools.research.stageU_mmp_interaction.mmp import (  # noqa: E402
    fragment, strip_stereochemistry, transformation,
)
from tools.research.stageV_core_mmp.core_mmp import (  # noqa: E402
    DoubleDifference, Observation, TargetEffect, double_differences,
    target_effects,
)

STAGE = Path(__file__).resolve().parents[1]
SLOW = os.environ.get("RUN_SLOW") == "1"

PREREG_SHA = "c567f66066c301fefe293048a4643fe4f65158077c3540ce1bbb0beb5d5844d4"


# -- 1. frozen artifact contract --------------------------------------------

def test_preregistration_sha256_is_frozen():
    actual = hashlib.sha256(
        (STAGE / "PREREGISTRATION.md").read_bytes()).hexdigest()
    assert actual == PREREG_SHA


def test_v0v1_result_records_the_negative_verdict():
    report = json.loads((STAGE / "V0_V1_RESULT.json").read_text(
        encoding="utf-8"))
    assert report["v0_gate"]["all_pass"] is False
    assert report["v0_gate"]["checks"]["top1_target_share"]["measured"] > 0.25
    assert report["v0_gate"]["checks"]["top1_component_share"]["measured"] > 0.25
    assert report["v0b_evaluability"]["primary_surface_internal_repeated_evaluable"] is False
    assert report["double_difference_surfaces"]["internal_repeated"]["rows"] < 100
    assert report["v1_gate"]["fit"]["pass"] is False
    assert report["v1_gate"]["fit"]["theta_hi"] < 0.0


def test_no_neural_model_was_trained_in_stage_v():
    assert not list(STAGE.rglob("runs/*/RUN.json"))


def test_posthoc_sensitivity_is_descriptive_and_consistent():
    report = json.loads((STAGE / "POSTHOC_FORENSICS.json").read_text(
        encoding="utf-8"))
    assert "excluded from every gate" in report["disclosure"]
    top = report["sensitivity_degree_domination"]["after_removing_top_target"]
    assert top["top1_target_share"] < 0.25
    assert top["top1_component_share"] < 0.25
    cross = report["sensitivity_interaction_variance"][
        "cross_component_effects_only"]
    assert cross["MS_effect"] < report["noise_reference_sigma2_same"]
    assert report["evaluation_surfaces"]["internal_repeated_exact"]["rows"] == 32


# -- 2. the estimand is constructed correctly --------------------------------

def _one_transformation(left_smiles, right_smiles):
    for left in fragment(left_smiles):
        for right in fragment(right_smiles):
            built = transformation(left, right)
            if built is not None:
                return built
    raise AssertionError("no matched pair")


def test_exact_key_includes_core_and_inverse_differs():
    aromatic, _ = _one_transformation("CCc1ccc(Cl)cc1", "CCc1ccc(Br)cc1")
    aliphatic, _ = _one_transformation("CCCCCl", "CCCCBr")
    assert aromatic.core != aliphatic.core
    assert aromatic.exact_key != aliphatic.exact_key
    inverse = aromatic.inverse
    assert inverse.r_a == aromatic.r_b and inverse.r_b == aromatic.r_a
    assert inverse.exact_key != aromatic.exact_key


def test_attachment_stereo_charge_preserved():
    stereo, _ = _one_transformation("C[C@H](N)c1ccccc1",
                                    "C[C@@H](N)c1ccccc1")
    assert stereo.stereo_edit is True
    assert strip_stereochemistry(stereo.r_a) == strip_stereochemistry(stereo.r_b)
    charged, _ = _one_transformation("CCc1ccc(CC(=O)O)cc1",
                                     "CCc1ccc(CC(=O)[O-])cc1")
    assert charged.charge_change != 0


def test_double_differences_are_only_within_one_exact_key():
    left = TargetEffect("key_a", "coarse_a", "t1", "c1", 1.0, 1, 1, False, (0.0,))
    right = TargetEffect("key_b", "coarse_b", "t2", "c2", 2.0, 1, 1, False, (0.0,))
    assert double_differences([left, right]) == []


def test_double_difference_cancels_a_shared_constant():
    """A shared mu_tau and shared target level cannot enter D."""
    mu = 7.5
    effects = [
        TargetEffect("key", "coarse", "t1", "c1", mu + 0.3, 1, 1, False, (0.0,)),
        TargetEffect("key", "coarse", "t2", "c2", mu - 0.2, 1, 1, False, (0.0,)),
    ]
    rows = double_differences(effects)
    assert len(rows) == 1
    assert abs(rows[0].value - 0.5) < 1e-12


# -- 3. no forbidden mechanics in the stage source ---------------------------

def test_no_python_hash_in_stage_v_sources():
    for path in sorted(STAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                    and node.func.id == "hash"):
                raise AssertionError(f"{path.name}:{node.lineno} calls hash()")


def test_stage_v_never_names_the_development_validation_split():
    forbidden = "meta" + "_val"
    for path in sorted(STAGE.rglob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value == forbidden:
                raise AssertionError(f"{path.name}:{node.lineno}")


def test_v0_census_has_no_label_path_into_keys_or_splits():
    source = (STAGE / "core_mmp.py").read_text(encoding="utf-8")
    tree = ast.parse(source, filename="core_mmp.py")
    # `delta_y` is computed only in build_observations from already-governed
    # cells; key construction and split membership never touch it.
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name in {
                "target_effects", "double_differences"}:
            for child in ast.walk(node):
                if isinstance(child, ast.Attribute) and child.attr == "pK":
                    raise AssertionError("split/key function reads pK")
                if isinstance(child, ast.Constant) and child.value == "pK":
                    raise AssertionError("split/key function reads pK")


# -- 4. governed seal (slow) -------------------------------------------------

@pytest.mark.skipif(not SLOW, reason="set RUN_SLOW=1 to mount the corpus")
def test_physical_meta_test_seal():
    from tools.research.stageV_core_mmp.core_mmp import load_governed
    data, seal = load_governed()
    assert seal["included"] is False
    assert seal["evaluated"] is False
    assert seal["isolation"]["physically_isolated"] is True
    assert "meta_test" not in data.tasks


def test_pair_level_noise_audit_is_descriptive_and_consistent():
    report = json.loads((STAGE / "PAIR_LEVEL_NOISE_AUDIT.json").read_text(
        encoding="utf-8"))
    assert "excluded from every gate" in report["disclosure"]
    population = report["population"]
    assert population["repeated_shared_panel_pairs"] == 88
    assert population["zero_range_curation_duplicate_groups"] == 40
    assert population["disagreeing_groups"] == 48
    disagreeing = report["pair_level_variance"]["disagreeing_pairs_only"]
    assert disagreeing["point"] < 0.858
    cross = report["v1_against_alternative_noise_references"][
        "pair_level_disagreeing_only"]["cross_component_keys_only"]
    assert cross["pass"] is False
    assert cross["theta_lo"] < 0.0


def test_metaval_structure_census_is_consistent():
    report = json.loads((STAGE / "METAVAL_STRUCTURE_CENSUS.json").read_text(
        encoding="utf-8"))
    assert report["overlap"]["exact_keys_shared_train_to_development"] == 0
    assert report["splits"]["development_validation"]["targets"] == 41
    assert report["splits"]["development_validation"]["components"] == 19
    assert report["splits"]["development_validation"][
        "potential_D_rows"] > 1000
    assert "no pK value is accessed" in report["disclosure"]


def test_phase1_final_decision_is_bounded_and_consistent():
    report = json.loads((STAGE / "PHASE1_FINAL_DECISION.json").read_text(
        encoding="utf-8"))
    assert report["verdict"] == "BOUNDED_NEGATIVE_UNDER_CURRENT_BINDINGDB_KI_PROTOCOL"
    assert report["verification"]["meta_test_evaluated"] is False
    assert report["verification"]["neural_models_trained_for_phase1_verdict"] is False
    assert report["verification"]["stage_v_preregistration_sha256"] == PREREG_SHA
    text = (STAGE / "PHASE1_FINAL_DECISION.md").read_text(encoding="utf-8")
    assert "BOUNDED NEGATIVE" in text
    assert "not biological impossibility" in text.lower() or \
           "not closed, and not claimed" in text.lower()


def test_remaining_lanes_audit_completes_the_negative_branch():
    report = json.loads((STAGE / "REMAINING_LANES_AUDIT.json").read_text(
        encoding="utf-8"))
    assert report["routes"]["MSA_coevolution"]["status"] == \
        "BLOCKED_ON_EXTERNAL_ASSET"
    assert report["routes"]["MSA_coevolution"]["evidence"][
        "governed_uniref_snapshot_found"] is False
    assert report["routes"]["external_datasets_Davis_KIBA"]["status"] == \
        "PROMOTION_GATED"
    assert report["routes"]["sealed_confirmation_split"]["status"] == \
        "SEALED_ZERO_EVALUATIONS"
    assert "negative branch is therefore complete" in report["conclusion"]


def test_v1_synthetic_calibration_is_bounded_and_consistent():
    report = json.loads((STAGE / "V1_SYNTHETIC_CALIBRATION.json").read_text(
        encoding="utf-8"))
    assert report["observed"]["MS_effect"] == 0.45168040420690786
    assert report["observed"]["df"] == 7482
    implied = report["implied_interaction_variance_if_all_excess_is_signal"]
    assert 0.3 < implied["pair_level_disagreeing_only"]["delta_sd"] < 0.5
    assert 0.4 < implied["pair_level_all_repeated"]["delta_sd"] < 0.7
    assert "cannot reopen the frozen gate" in report["reading"]


def test_completion_evidence_manifest_hashes_match_current_files():
    import hashlib
    manifest_path = REPO / "report" / "COMPLETION_EVIDENCE_MANIFEST_CORE_TASK1.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["goal_status"].startswith("COMPLETE (bounded-negative branch")
    for item in manifest["artifacts"]:
        path = REPO / item["path"]
        assert path.exists(), item["path"]
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == item["sha256"], item["path"]
    goal = (REPO / "tools" / "research" / "GOAL_ACTIVE.md").read_text(
        encoding="utf-8")
    assert "COMPLETE（有界负结论分支）" in goal
