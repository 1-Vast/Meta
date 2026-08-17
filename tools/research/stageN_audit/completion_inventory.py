"""Completion-evidence inventory: verify that every artifact referenced by the
final authorities exists and carries the expected schema/verdict fields.

Read-only. Output: tools/research/stageN_audit/COMPLETION_INVENTORY.json
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

CHECKS = [
    # (description, path, kind)
    ("boundary doc", "report/BOUNDARY_20260817_NIGHT.md", "file"),
    ("final state doc", "report/FINAL_STATE_20260818.md", "file"),
    ("evidence ledger", "report/EVIDENCE_LEDGER.md", "file"),
    ("current model evidence", "report/CURRENT_MODEL_EVIDENCE.md", "file"),
    ("task contract", "task.md", "file"),
    ("history", "history.md", "file"),
    ("goal record", "tools/research/GOAL_ACTIVE.md", "file"),
    ("final audit", "tools/research/stageN_audit/FINAL_BOUNDARY_AUDIT.json", "json"),
    ("audit report", "tools/research/stageN_audit/AUDIT_REPORT.md", "file"),
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

STAGES = ["stageD_level_panel", "stageF_pairwise", "stageG_esm650",
          "stageI_lm", "stageJ_assay", "stageK_contrastive",
          "stageL_gated", "stageQ_frozenhead"]


def main():
    report = {"checked": [], "stages": {}, "missing": []}
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
    for stage in STAGES:
        base = ROOT / "tools/research" / stage
        pre = (base / "PREREGISTRATION.md").is_file()
        report_files = sorted(base.glob("REPORT.md"))
        results = sorted(base.glob("*.json"))
        report["stages"][stage] = {"preregistration": pre,
                                   "reports": [p.name for p in report_files],
                                   "json_artifacts": len(results)}
        print(f"STAGE {stage}: pre={pre} reports={len(report_files)} "
              f"json={len(results)}")
    out = ROOT / "tools/research/stageN_audit/COMPLETION_INVENTORY.json"
    out.write_text(json.dumps(report, indent=1, sort_keys=True),
                   encoding="utf-8")
    print("missing:", report["missing"])
    print("wrote", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
