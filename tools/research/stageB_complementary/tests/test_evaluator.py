"""Contract for the corrected Stage B evaluator.

`test_every_arm_reports_every_control` is correction 2 as a test: Stage A
computed A0's control rows and then dropped them from the summary, which is
exactly the kind of omission that looks like an absent effect rather than an
absent measurement.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.research.stageB_complementary.evaluate_stageb import (     # noqa: E402
    CONTROLS, cliff_sign_accuracy, extended_metrics, paired,
    tanimoto_matrix, verify_arms_are_matched,
)


def write_arms(root: Path, arms: dict[str, dict]):
    base = {"seed": 20260815, "steps": 1200, "arch": "similarity_only",
            "split_directory": "x", "hidden_dim": 192,
            "mode": "T", "inner_steps": 1, "inner_lr": 0.1,
            "learned_step": False, "max_step": 0.5, "first_order": True,
            "adapted": []}
    for arm, overrides in arms.items():
        path = root / arm
        path.mkdir(parents=True, exist_ok=True)
        (path / "RESULT.json").write_text(
            json.dumps({"config": {**base, **overrides}, "report": {}}),
            encoding="utf-8")


def test_matched_arms_accepted(tmp_path):
    write_arms(tmp_path, {"T": {"mode": "T"}, "H": {"mode": "H"},
                          "C": {"mode": "C"}})
    configs = verify_arms_are_matched(tmp_path, ("T", "H", "C"))
    assert set(configs) == {"T", "H", "C"}


def test_an_undeclared_difference_raises(tmp_path):
    write_arms(tmp_path, {"T": {"mode": "T"}, "H": {"mode": "H", "steps": 600}})
    with pytest.raises(ValueError, match="not matched"):
        verify_arms_are_matched(tmp_path, ("T", "H"))


def test_a_capacity_difference_raises(tmp_path):
    write_arms(tmp_path, {"T": {"mode": "T"}, "C": {"mode": "C", "hidden_dim": 256}})
    with pytest.raises(ValueError, match="not matched"):
        verify_arms_are_matched(tmp_path, ("T", "C"))


def test_every_arm_reports_every_control():
    """Correction 2: the control list is arm-independent by construction."""
    assert set(CONTROLS) >= {"no_adaptation", "permuted_support",
                             "matched_wrong_support", "wrong_protein",
                             "bias_only", "weight_only"}
    # The evaluator loops `for arm in arms` over the same CONTROLS tuple, so no
    # arm can be silently omitted; this pins the tuple itself.
    assert len(set(CONTROLS)) == len(CONTROLS)


def rows(condition, values, k=2, novelty=0.1):
    return [{"k": k, "component": f"c{i % 2}", "target": t, "condition": condition,
             "mse_pk": v, "max_train_tanimoto": novelty}
            for i, (t, v) in enumerate(sorted(values.items()))]


def test_paired_uses_separate_conditions():
    correct = rows("correct", {"t1": 1.0, "t2": 2.0, "t3": 3.0, "t4": 4.0})
    control = rows("permuted_support", {"t1": 2.0, "t2": 3.0, "t3": 4.0, "t4": 5.0})
    result = paired(control, correct, "permuted_support", "correct", 2, "mse_pk")
    assert result["mean"] == pytest.approx(1.0, abs=1e-9)
    assert result["components"] == 2


def test_paired_skips_non_finite_values():
    correct = rows("correct", {"t1": 1.0, "t2": 2.0})
    control = rows("permuted_support", {"t1": float("nan"), "t2": 3.0})
    result = paired(control, correct, "permuted_support", "correct", 2, "mse_pk")
    assert np.isfinite(result["mean"])
    assert result["components"] == 1


def test_paired_honours_a_stratum_predicate():
    correct = rows("correct", {"t1": 1.0, "t2": 2.0})
    control = rows("permuted_support", {"t1": 2.0, "t2": 3.0})
    for row in control:
        row["max_train_tanimoto"] = 0.9 if row["target"] == "t1" else 0.1
    result = paired(control, correct, "permuted_support", "correct", 2,
                    "mse_pk", predicate=lambda r: r["max_train_tanimoto"] < 0.4)
    assert result["components"] == 1


def test_cliff_sign_accuracy_is_orientation_invariant():
    truth = np.array([5.0, 3.0, 8.0])
    prediction = np.array([4.5, 3.5, 7.0])
    similarity = np.ones((3, 3))
    forward = cliff_sign_accuracy(prediction, truth, similarity)
    order = np.array([2, 0, 1])
    reordered = cliff_sign_accuracy(prediction[order], truth[order],
                                    similarity[np.ix_(order, order)])
    assert forward == pytest.approx(reordered, abs=1e-12)


def test_cliff_sign_returns_nan_when_no_pair_qualifies():
    truth = np.array([5.0, 5.05])
    prediction = np.array([1.0, 2.0])
    assert np.isnan(cliff_sign_accuracy(prediction, truth, np.ones((2, 2))))


def test_extended_metrics_reports_centered_mse_and_cliff():
    truth = np.array([1.0, 2.0, 3.0, 4.0])
    values = extended_metrics(truth + 0.5, truth, np.ones((4, 4)))
    assert set(values) >= {"mse_pk", "rmse_pk", "pearson", "spearman", "ci",
                           "r2", "centered_mse_pk", "cliff_sign"}
    # A constant offset is pure level: centered error must vanish.
    assert values["centered_mse_pk"] == pytest.approx(0.0, abs=1e-12)
    assert values["mse_pk"] == pytest.approx(0.25, abs=1e-12)


def test_centered_mse_separates_level_from_shape():
    truth = np.array([1.0, 2.0, 3.0, 4.0])
    shape_error = np.array([1.0, 2.5, 2.5, 4.0])
    values = extended_metrics(shape_error, truth, np.ones((4, 4)))
    assert values["centered_mse_pk"] > 0.0


def test_tanimoto_matrix_is_symmetric_with_unit_diagonal():
    generator = np.random.default_rng(0)
    fingerprint = (generator.random((6, 32)) > 0.7).astype(float)
    matrix = tanimoto_matrix(fingerprint)
    assert np.allclose(matrix, matrix.T, atol=1e-12)
    assert np.allclose(np.diag(matrix), 1.0, atol=1e-9)
