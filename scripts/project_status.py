"""Read-only status view over the governed project and archive manifests."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "archive/retired_research_20260811/ARCHIVE_MANIFEST.json"

STATE = {
    "schema_version": "5.0",
    "updated": "2026-08-11",
    "core_task": {
        "name": "unseen-target few-shot drug-target affinity prediction",
        "target_is_meta_task": True,
        "support_sizes": [1, 2, 3, 5],
        "query_labels_forbidden": True,
    },
    "current_stage": {
        "name": "R0C_EXACT_DISTANCE_CONFIRMATION_COMPLETE",
        "status": "MARGINAL_OR_SLOT_RECALIBRATION_ONLY",
        "training_authorized": False,
        "confirmation_labels_open": False,
    },
    "unresolved": [
        "FRESH_CLUSTER_CONFIRMATION_FOR_PROTEIN_SPECIFICITY",
        "INDEPENDENT_ASSAY_MATCHED_DENSE_CROSSED_SELECTIVITY_COHORT",
        "PREAGGREGATION_AND_LOCAL_INFORMATION_AUDIT_CONDITIONAL_ON_FRESH_A1_PASS",
        "END_TO_END_BIOLOGICAL_COORDINATE_TRAINING",
        "NEW_OBSERVABLE_TO_BREAK_ADDITIVE_ATOM_RESIDUE_SHORTCUT",
        "REAL_ASSAY_REPLICATE_UNCERTAINTY_FOR_ROBUST_META_OBJECTIVE",
        "PARTNER_SPECIFIC_AFFINITY_INCREMENT",
        "BIOLOGICAL_STATISTIC_ADMISSION_TO_Z",
        "END_TO_END_FEWSHOT_DTA",
        "CSMO_BIOLOGICAL_Z_BRIDGE",
    ],
    "archive_manifest": "archive/retired_research_20260811/ARCHIVE_MANIFEST.json",
}


def load_status(*, archive_only: bool = False) -> dict:
    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    if archive_only:
        return {
            "archive_manifest": ARCHIVE.relative_to(ROOT).as_posix(),
            "created": archive["created"],
            "purpose": archive["purpose"],
            "policy": archive["policy"],
            "production_impact": archive["production_impact"],
        }
    return {
        "schema_version": STATE["schema_version"],
        "updated": STATE["updated"],
        "core_task": STATE["core_task"],
        "current_stage": STATE["current_stage"],
        "unresolved": STATE["unresolved"],
        "archive_manifest": STATE["archive_manifest"],
    }


def _human_status(value: dict, *, archive_only: bool) -> str:
    if archive_only:
        return "\n".join((
            f"Archive: {value['archive_manifest']}",
            f"Created: {value['created']}",
            f"Policy: {value['policy']}",
            f"Production impact: {value['production_impact']}",
        ))
    stage = value["current_stage"]
    return "\n".join((
        f"Task: {value['core_task']['name']}",
        f"Stage: {stage['name']}",
        f"Status: {stage['status']}",
        f"Training authorized: {str(stage['training_authorized']).lower()}",
        f"Unresolved items: {len(value['unresolved'])}",
        f"Archive: {value['archive_manifest']}",
    ))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", action="store_true")
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    value = load_status(archive_only=args.archive)
    print(json.dumps(value, indent=2, sort_keys=True) if args.json
          else _human_status(value, archive_only=args.archive))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
