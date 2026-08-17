"""Contract for the Stage A evaluator.

The matched-arm guard and the paired contrast are the two places where a silent
defect would produce an interpretable-looking but wrong result, so both are
pinned here. `test_paired_is_not_self_cancelling` exists because the first
version of `paired` took one condition for both sides, which made every
counterfactual an exact zero that looked like a clean null.
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

from tools.research.stageA_innerloop.evaluate import (                # noqa: E402
    ARMS, concordance, metrics_for, paired, pearson, r_squared, spearman,
    verify_arms_are_matched,
)


def write_arms(root: Path, **overrides):
    base = {"seed": 20260815, "steps": 1200, "arch": "similarity_only",
            "split_directory": "x", "inner_steps": 0, "inner_lr": 0.1,
            "task_selection": False, "hidden_dim": 192}
    settings = {
        "A0": {},
        "A1": {"inner_steps": 1},
        "A2": {"inner_steps": 1, "task_selection": True},
    }
    for arm in ARMS:
        config = {**base, "arm": arm, **settings[arm], **overrides.get(arm, {})}
        path = root / arm
        path.mkdir(parents=True, exist_ok=True)
        (path / "RESULT.json").write_text(
            json.dumps({"config": config, "report": {}}), encoding="utf-8")


def test_matched_arms_are_accepted(tmp_path):
    write_arms(tmp_path)
    configs = verify_arms_are_matched(tmp_path)
    assert set(configs) == set(ARMS)


def test_an_undeclared_difference_raises(tmp_path):
    write_arms(tmp_path, A1={"steps": 600})
    with pytest.raises(ValueError, match="not matched"):
        verify_arms_are_matched(tmp_path)


def test_a_capacity_difference_raises(tmp_path):
    write_arms(tmp_path, A2={"hidden_dim": 256})
    with pytest.raises(ValueError, match="not matched"):
        verify_arms_are_matched(tmp_path)


def test_an_arm_that_does_not_adapt_raises(tmp_path):
    write_arms(tmp_path, A1={"inner_steps": 0})
    with pytest.raises(ValueError, match="must actually adapt"):
        verify_arms_are_matched(tmp_path)


def test_a_baseline_that_adapts_raises(tmp_path):
    write_arms(tmp_path, A0={"inner_steps": 1})
    with pytest.raises(ValueError):
        verify_arms_are_matched(tmp_path)


def test_a_missing_arm_raises(tmp_path):
    write_arms(tmp_path)
    (tmp_path / "A2" / "RESULT.json").unlink()
    with pytest.raises(FileNotFoundError):
        verify_arms_are_matched(tmp_path)


def test_selection_flags_must_match_their_arm(tmp_path):
    write_arms(tmp_path, A1={"inner_steps": 1, "task_selection": True})
    with pytest.raises(ValueError, match="uniform"):
        verify_arms_are_matched(tmp_path)


# --- the paired contrast ---------------------------------------------------

def rows(condition: str, values: dict[str, float], k: int = 1):
    return [{"k": k, "component": f"c{i % 2}", "target": target,
             "condition": condition, "mse_pk": value}
            for i, (target, value) in enumerate(sorted(values.items()))]


def test_paired_is_not_self_cancelling():
    """Two different conditions must produce a non-zero contrast."""
    correct = rows("steps1", {"t1": 1.0, "t2": 2.0, "t3": 3.0, "t4": 4.0})
    control = rows("permuted_support", {"t1": 2.0, "t2": 3.0, "t3": 4.0, "t4": 5.0})
    result = paired(control, correct, "permuted_support", "steps1", 1, "mse_pk")
    assert result["mean"] == pytest.approx(1.0, abs=1e-9)
    assert result["components"] == 2


def test_paired_against_the_same_condition_is_exactly_zero():
    correct = rows("steps1", {"t1": 1.0, "t2": 2.0})
    result = paired(correct, correct, "steps1", "steps1", 1, "mse_pk")
    assert result["mean"] == pytest.approx(0.0, abs=1e-12)


def test_paired_ignores_targets_missing_from_one_side():
    left = rows("permuted_support", {"t1": 2.0, "t2": 3.0, "t9": 9.0})
    right = rows("steps1", {"t1": 1.0, "t2": 1.0})
    result = paired(left, right, "permuted_support", "steps1", 1, "mse_pk")
    assert result["components"] == 2
    assert np.isfinite(result["mean"])


def test_paired_matches_repeated_draws_seat_by_seat():
    """Two draws of one target must pair 1-to-1, not fan out to four pairs."""
    left = ([{"k": 1, "component": "c0", "target": "t1",
              "condition": "a", "mse_pk": value} for value in (3.0, 5.0)])
    right = ([{"k": 1, "component": "c0", "target": "t1",
               "condition": "b", "mse_pk": value} for value in (1.0, 2.0)])
    result = paired(left, right, "a", "b", 1, "mse_pk")
    assert result["mean"] == pytest.approx(2.5, abs=1e-9)   # (2 + 3) / 2


# --- metrics ---------------------------------------------------------------

def test_pearson_and_spearman_agree_on_a_monotone_linear_case():
    truth = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    prediction = 2.0 * truth + 1.0
    assert pearson(prediction, truth) == pytest.approx(1.0, abs=1e-12)
    assert spearman(prediction, truth) == pytest.approx(1.0, abs=1e-12)
    assert concordance(prediction, truth) == pytest.approx(1.0, abs=1e-12)


def test_spearman_survives_a_monotone_nonlinearity_that_pearson_loses():
    truth = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    prediction = np.exp(truth)
    assert spearman(prediction, truth) == pytest.approx(1.0, abs=1e-12)
    assert pearson(prediction, truth) < 0.95


def test_r2_is_negative_for_a_prediction_worse_than_the_mean():
    truth = np.array([1.0, 2.0, 3.0, 4.0])
    assert r_squared(np.full(4, 10.0), truth) < 0.0
    assert r_squared(truth, truth) == pytest.approx(1.0, abs=1e-12)
    assert r_squared(np.full(4, truth.mean()), truth) == pytest.approx(0.0, abs=1e-12)


def test_reversed_prediction_inverts_every_ranking_metric():
    truth = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    prediction = -truth
    assert pearson(prediction, truth) == pytest.approx(-1.0, abs=1e-12)
    assert spearman(prediction, truth) == pytest.approx(-1.0, abs=1e-12)
    assert concordance(prediction, truth) == pytest.approx(0.0, abs=1e-12)


def test_metrics_for_reports_every_preregistered_field():
    truth = np.array([1.0, 2.0, 3.0, 4.0])
    values = metrics_for(truth + 0.5, truth)
    assert set(values) == {"mse_pk", "rmse_pk", "pearson", "spearman", "ci", "r2"}
    assert values["mse_pk"] == pytest.approx(0.25, abs=1e-12)
    assert values["rmse_pk"] == pytest.approx(0.5, abs=1e-12)


def test_a_constant_truth_panel_does_not_crash_the_ranking_metrics():
    truth = np.full(4, 2.0)
    values = metrics_for(np.array([1.0, 2.0, 3.0, 4.0]), truth)
    assert np.isnan(values["ci"])
    assert values["pearson"] == 0.0
