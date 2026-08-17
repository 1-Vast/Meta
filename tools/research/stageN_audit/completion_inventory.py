"""Completion-evidence inventory: verify that every artifact referenced by the
final authorities exists and carries the expected schema/verdict fields.

Read-only. Output: tools/research/stageN_audit/COMPLETION_INVENTORY.json

The stage list is **discovered from the filesystem** by
`stage_inventory.discover_stages`, the same rule `final_audit.py` uses. Until
2026-08-18 this file hard-coded eight stages while the final audit hard-coded
seven, which is the discrepancy the post-completion review recorded; both lists
also missed Stages A, B and P_cpc, which emit row artifacts without a summary
sidecar.
"""
from __future__ import annotations

import datetime as _datetime
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from tools.research.stageN_audit.stage_inventory import discover_stages

CHECKS = [
    # (description, path, kind)
    ("post-completion review", "report/POST_COMPLETION_REVIEW_20260818.md", "file"),
    ("boundary doc", "report/BOUNDARY_20260817_NIGHT.md", "file"),
    ("final state doc", "report/FINAL_STATE_20260818.md", "file"),
    ("evidence ledger", "report/EVIDENCE_LEDGER.md", "file"),
    ("current model evidence", "report/CURRENT_MODEL_EVIDENCE.md", "file"),
    ("completion statement", "report/COMPLETION_STATEMENT_20260818.md", "file"),
    ("task contract", "task.md", "file"),
    ("history", "history.md", "file"),
    ("goal record", "tools/research/GOAL_ACTIVE.md", "file"),
    ("final audit", "tools/research/stageN_audit/FINAL_BOUNDARY_AUDIT.json", "json"),
    ("audit report", "tools/research/stageN_audit/AUDIT_REPORT.md", "file"),
    ("method-ladder closure map", "tools/research/method_ladder/CLOSURE_MAP.md", "file"),
    ("split isolation spec", "tools/research/a2_readiness_v2/SPLIT_ISOLATION_SPEC.md", "file"),
    ("governed split-view builder", "scripts/build_governed_split_views.py", "file"),
    ("physical seal contract", "tools/tests/test_physical_meta_test_seal.py", "file"),
    ("leak-free selection module", "scripts/internal_validation.py", "file"),
    ("D0 report", "tools/research/stageD_level_panel/D0_REPORT.md", "file"),
    ("D0b doc transfer", "tools/research/stageD_level_panel/D0b_DOC_TRANSFER.json", "json"),
    ("D0c journal", "tools/research/stageJ_assay/D0c_JOURNAL_IDENTIFIABILITY.json", "json"),
    ("T2 eval rows", "tools/research/stageD_level_panel/T2_meta_val.rows.jsonl", "file"),
    ("K2 pooled contrast", "tools/research/stageK_contrastive/K2_multiseed_contrast.json", "json"),
    ("L contrast", "tools/research/stageL_gated/L_vs_T2.contrast.json", "json"),
    ("Q contrast", "tools/research/stageQ_frozenhead/Q_vs_T2.contrast.json", "json"),
    ("Q0 probe", "tools/research/stageQ_frozenhead/Q0_JOINT_FROZEN_IDENTIFIABILITY.json", "json"),
    ("M0 probes", "tools/research/stageM_chemberta/M0_CHEMBERTA_PROBES.json", "json"),
    ("P0 probes", "tools/research/stageP_go/P0_GO_PROBES.json", "json"),
    ("Davis/KIBA frozen plan", "tools/research/stageR_daviskiba/PREREGISTRATION.md", "file"),
]


def main():
    report = {"schema": "MetaSieve.CompletionInventory.v2",
              "generated": _datetime.date.today().isoformat(),
              "checked": [], "stages": {}, "missing": []}
    for description, rel, kind in CHECKS:
        path = ROOT / rel
        ok = path.is_file()
        entry = {"path": rel, "exists": ok}
        if ok and kind == "json":
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
                entry["schema"] = payload.get("schema", payload.get("date", "n/a"))
            except Exception:  # noqa: BLE001
                entry["json_ok"] = False
        report["checked"].append(entry)
        if not ok:
            report["missing"].append(rel)
        print(("OK   " if ok else "MISS ") + description)

    discovered = discover_stages(ROOT)
    trained = [row for row in discovered if row["trained"]]
    for row in trained:
        report["stages"][row["stage"]] = {
            "preregistration": row["preregistered"],
            "preregistration_files": row["preregistration_files"],
            "preregistration_precedes_results": row[
                "preregistration_precedes_results"],
            "reports": row["reports"],
            "row_artifacts": row["row_artifacts"],
            "json_artifacts": row["json_artifacts"],
        }
        print(f"STAGE {row['stage']}: pre={row['preregistered']} "
              f"reports={len(row['reports'])} rows={row['row_artifacts']} "
              f"json={row['json_artifacts']}")

    report["stage_summary"] = {
        "trained_stages_retained": len(trained),
        "trained_and_preregistered": sum(1 for row in trained
                                         if row["preregistered"]),
        "preregistered_not_run": sorted(
            row["stage"] for row in discovered
            if row["preregistered"] and not row["trained"]),
        "measurement_only_directories": sorted(
            row["stage"] for row in discovered
            if not row["preregistered"] and not row["trained"]),
        "discovery_rule": "stage_inventory.discover_stages",
    }
    print("trained stages:", report["stage_summary"]["trained_stages_retained"])

    out = ROOT / "tools/research/stageN_audit/COMPLETION_INVENTORY.json"
    out.write_text(json.dumps(report, indent=1, sort_keys=True),
                   encoding="utf-8")
    print("missing:", report["missing"])
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
