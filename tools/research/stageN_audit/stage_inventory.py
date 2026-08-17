"""One discovery rule for "which stages exist", shared by both audits.

Before 2026-08-18 the two audit scripts each carried their own hand-written
stage list: `final_audit.py` named seven, `completion_inventory.py` named
eight. Neither was wrong when written and both went stale the moment a stage
was added — Stage Q is exactly the failure that produced the seven-versus-eight
discrepancy the post-completion review recorded. Worse, both lists selected on
`*.rows.summary.json`, which Stages A, B and P_cpc never emitted, so three
trained stages were invisible to *both* counts.

The rule below reads the filesystem instead:

* a **stage directory** is any `tools/research/stage*/` directory;
* it is **preregistered** if it carries any `PREREGISTRATION*.md`;
* it is **trained/retained** if it also carries at least one evaluation row
  artifact (`*.rows.jsonl`) — the artifact every trained stage in this
  programme emits, with or without a summary sidecar;
* preregistration **precedes results** when the preregistration file's mtime is
  not later than the earliest row artifact's.

The mtime comparison is a weak check on a strong claim, and it is reported as
its own field rather than folded into a pass/fail: Git history is the real
ordering evidence. It exists to catch the one mistake it can catch — a
preregistration written after the fact.
"""
from __future__ import annotations

from pathlib import Path


def discover_stages(root: Path) -> list[dict]:
    """Every `tools/research/stage*` directory, classified by what it holds."""
    stages: list[dict] = []
    for directory in sorted((root / "tools/research").glob("stage*")):
        if not directory.is_dir():
            continue
        preregistrations = sorted(directory.glob("PREREGISTRATION*.md"))
        rows = sorted(directory.glob("*.rows.jsonl"))
        summaries = sorted(directory.glob("*.rows.summary.json"))
        results = sorted(directory.glob("RESULT.json"))
        reports = sorted(directory.glob("REPORT*.md"))
        json_artifacts = sorted(
            path for path in directory.glob("*.json")
            if path.name != "RESULT.json")

        precedes = None
        if preregistrations and rows:
            precedes = (min(p.stat().st_mtime for p in preregistrations)
                        <= min(r.stat().st_mtime for r in rows))

        stages.append({
            "stage": directory.name,
            "preregistered": bool(preregistrations),
            "preregistration_files": [p.name for p in preregistrations],
            "trained": bool(preregistrations and rows),
            "row_artifacts": len(rows),
            "summary_artifacts": len(summaries),
            "result_artifacts": len(results),
            "json_artifacts": len(json_artifacts),
            "reports": [p.name for p in reports],
            "preregistration_precedes_results": precedes,
        })
    return stages


def stage_table(stages: list[dict]) -> list[dict]:
    """Stable, JSON-serialisable ordering for the audit outputs."""
    return sorted(stages, key=lambda row: row["stage"])


def trained_stage_names(root: Path) -> list[str]:
    return [row["stage"] for row in discover_stages(root) if row["trained"]]
