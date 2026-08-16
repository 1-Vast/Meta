"""Keep the narrative research record consistent with its leaf artifacts.

These are cheap, no-training checks. They exist because the R13 gate-count
inconsistency (`15/16 passed` alongside two recorded xfails) and the mixed
k=0 "frontier" both survived several consolidation passes: a number that no
test recomputes will drift away from its evidence.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.audit_research_record import collect_arms, check_seals, pareto

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "report"
BOUNDARY = REPORT / "BOUNDARY_20260816.md"
R13 = REPORT / "meta_fewshot" / "stageR13_shape_direct_20260816"


def test_double_cold_meta_test_is_never_reported_as_evaluated():
    """The confirmation split opens once, and it has not opened yet."""
    seals = check_seals()
    assert seals["violations"] == [], seals["violations"]
    assert seals["sealed_explicit"], "no explicitly sealed double-cold artifact found"


def test_r13_gate_count_matches_the_suite():
    """The R13 record must state the count the suite actually collects."""
    collected = subprocess.run(
        [sys.executable, "-m", "pytest", "tests/test_shape_direct_synthetic.py",
         "--collect-only", "-q", "-p", "no:cacheprovider"],
        cwd=ROOT, capture_output=True, text=True, timeout=300)
    assert collected.returncode == 0, collected.stderr[-2000:]
    total = sum(1 for line in collected.stdout.splitlines()
                if line.startswith("tests/test_shape_direct_synthetic.py::"))

    record = json.loads((R13 / "RESULT.json").read_text(encoding="utf-8"))
    suite = record["gate_suite"]
    assert suite["collected"] == total, (
        f"RESULT.json says {suite['collected']} gates, suite collects {total}")
    assert suite["passed"] + suite["xfailed"] == total, (
        f"{suite['passed']} passed + {suite['xfailed']} xfailed != {total}")
    assert len(record["xfail_gates"]) == suite["xfailed"]

    narrative = (R13 / "REPORT.md").read_text(encoding="utf-8")
    assert f"**{suite['passed']} of {total} gates pass**" in narrative, (
        "REPORT.md does not state the pass/total the RESULT.json records")


def test_k0_frontier_in_the_boundary_document_is_the_real_pareto_set():
    """A frontier is a set of whole configurations, not a mix of metrics."""
    arms = collect_arms()
    front = pareto(arms)
    assert front == ["B3", "C2", "A0"], front

    text = BOUNDARY.read_text(encoding="utf-8")
    for label in front:
        mse = arms[label]["0"]["mse"]
        ci = arms[label]["0"]["ci"]
        assert f"{mse:.3f}" in text, f"{label} MSE {mse:.3f} absent from BOUNDARY"
        assert f"{ci:.3f}" in text, f"{label} CI {ci:.3f} absent from BOUNDARY"

    # Every arm on the frontier must be beaten on the other metric by another
    # frontier member; that is what makes the trade real rather than rhetorical.
    best_mse = min(arms[label]["0"]["mse"] for label in front)
    best_ci = max(arms[label]["0"]["ci"] for label in front)
    assert not any(arms[label]["0"]["mse"] == best_mse
                   and arms[label]["0"]["ci"] == best_ci for label in front), (
        "one arm dominates both metrics; the frontier language must be rewritten")


def test_best_cliff_record_is_scoped_to_meta_val_development():
    """0.782 is a development record on a Pareto-dominated arm; say so."""
    arms = collect_arms()
    best = max(arms.items(), key=lambda kv: kv[1].get("5", {}).get("cliff_sign") or 0)
    label, row = best
    assert label == "C1" and abs(row["5"]["cliff_sign"] - 0.782) < 5e-4, best

    assert label not in pareto(arms), (
        "C1 is no longer Pareto-dominated; the BOUNDARY caveat needs rewriting")

    text = BOUNDARY.read_text(encoding="utf-8")
    heading = "## The activity-cliff record, stated at its true scope"
    assert heading in text, "BOUNDARY has no scoped activity-cliff section"
    section = text[text.index(heading):]
    section = section[:section.index("\n## ", 1)]
    assert "0.782" in section, "the cliff record section omits the number"
    assert "meta_val" in section, "the cliff record is not scoped to meta_val"
    assert "development" in section, "the cliff record is not graded development"
    assert "dominated" in section, "the cliff record omits its dominated status"


@pytest.mark.parametrize("path,expected", [
    ("docs/PROJECT_FILE_ORGANIZATION.md", "R5-R13"),
    ("report/BOUNDARY_20260816.md", "sealed and never opened"),
])
def test_documents_carry_the_current_cycle_scope(path, expected):
    assert expected in (ROOT / path).read_text(encoding="utf-8")
