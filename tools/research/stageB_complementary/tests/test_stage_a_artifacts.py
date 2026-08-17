"""Numerical regression tests against the recorded Stage A evidence.

These read `tools/research/stageA_innerloop/` and assert the facts the
correction audit rests on. They exist so a future reader cannot quietly restate
a withdrawn claim: if someone re-reports "A1 improves every metric at every k",
`test_k5_ci_regression_is_recorded` is the counter-evidence, taken from the leaf
artifact rather than from prose.

Stage A's artifacts are never modified by this suite.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[4]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

STAGE_A = ROOT / "tools/research/stageA_innerloop"
EVAL = STAGE_A / "STAGE_A_meta_val.json"


@pytest.fixture(scope="module")
def payload():
    if not EVAL.exists():
        pytest.skip("Stage A evaluation artifact is absent")
    return json.loads(EVAL.read_text(encoding="utf-8"))


def test_k5_ci_regression_is_recorded(payload):
    """Correction 7: A1's k=5 CI is BELOW A0's, so 'every metric' is false."""
    a0 = payload["arm_metrics"]["A0"]["5"]["steps0"]["ci"]
    a1 = payload["arm_metrics"]["A1"]["5"]["steps1"]["ci"]
    assert a1 < a0, "the recorded k=5 CI regression must remain visible"
    assert a0 == pytest.approx(0.6314, abs=5e-4)
    assert a1 == pytest.approx(0.6295, abs=5e-4)
    # Spearman disagrees with CI at k=5, which is itself the finding.
    assert (payload["arm_metrics"]["A1"]["5"]["steps1"]["spearman"] >
            payload["arm_metrics"]["A0"]["5"]["steps0"]["spearman"])


def test_a0_counterfactuals_were_never_summarized(payload):
    """Correction 2: the defect this stage repairs, pinned as a fact."""
    assert "A0" not in payload["counterfactuals"], (
        "Stage A omitted A0's controls; if this now passes, the artifact was "
        "edited and the audit needs revising")
    assert set(payload["counterfactuals"]) == {"A1", "A2"}


def test_the_rows_do_contain_a0_controls(payload):
    """The rows exist — only the summary omitted them, which is why a re-run
    is needed for correction 1 but not for the permutation control."""
    rows = STAGE_A / "STAGE_A_meta_val.rows.jsonl"
    if not rows.exists():
        pytest.skip("Stage A prediction rows are absent")
    conditions = set()
    with rows.open(encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if row["arm"] == "A0" and row["k"] > 1:
                conditions.add(row["condition"])
    assert "permuted_support" in conditions
    assert "matched_wrong_support" in conditions


def test_meta_test_status_is_the_audited_string(payload):
    """Correction 6: the defensible claim, not 'untouched'."""
    seal = payload["meta_test"]
    assert seal["included"] is False
    assert seal["evaluated"] is False
    assert seal["sealed_cells_withheld"] == 768
    assert seal["isolation"]["level"] == "logical_exclusion_after_parsing"
    assert seal["isolation"]["physically_isolated"] is False
    assert seal["isolation"]["labels_parsed_in_process"] is True


def test_stage_a_conditioning_used_the_uncorrected_formula():
    """Correction 3: pin the published value so the correction is auditable."""
    path = STAGE_A / "CONDITIONING.json"
    if not path.exists():
        pytest.skip("Stage A conditioning artifact is absent")
    recorded = json.loads(path.read_text(encoding="utf-8"))
    assert "2 * inner_lr * |h|^2" in recorded["definition"]
    assert "+ 1" not in recorded["definition"], (
        "the published Stage A definition omitted the adapted bias")
    assert recorded["arms"]["A0"]["alpha_mean"] == pytest.approx(1.514, abs=5e-3)
    lr = recorded["inner_lr"]
    # The corrected alpha adds exactly 2*lr for the adapted bias.
    corrected_a0 = recorded["arms"]["A0"]["alpha_mean"] + 2.0 * lr
    corrected_a1 = recorded["arms"]["A1"]["alpha_mean"] + 2.0 * lr
    assert corrected_a0 == pytest.approx(1.714, abs=5e-3)
    assert corrected_a1 == pytest.approx(0.441, abs=5e-3)
    # The qualitative conclusion survives: A0 overshoots, A1 does not.
    assert corrected_a0 > 1.0 and corrected_a1 < 1.0


def test_stage_a_artifacts_are_present_and_unmodified_by_this_stage():
    for name in ("PREREGISTRATION.md", "REPORT.md", "RESULT.json"):
        assert (STAGE_A / name).exists(), f"Stage A {name} must be preserved"
