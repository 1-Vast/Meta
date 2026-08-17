"""Read-only status view over the governed cold-target DTA project."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STATE = {
    "schema_version": "7.0",
    "updated": "2026-08-16",
    "core_task": {
        "name": "cold-target few-shot drug-target affinity prediction",
        "target_is_meta_task": True,
        "support_sizes": [0, 1, 2, 3, 5],
        "query_labels_forbidden": True,
    },
    "current_stage": {
        "name": "A2_MOMENT_META_PREREGISTRATION",
        "status": "R14_BOUNDARY_ESTABLISHED_NEXT_FAMILY_NOT_RUN",
        "training_authorized": True,
        "confirmation_labels_open": False,
    },
    "unresolved": [
        "PROTEIN_CONDITIONED_SAR_COORDINATE_IDENTIFIABILITY",
        "K1_NON_SCALAR_ADAPTATION",
        "IMPROVEMENT_OVER_FIXED_TANIMOTO",
        "UNTOUCHED_DOUBLE_COLD_CONFIRMATION",
    ],
    "cleanup_record": "docs/REPOSITORY_CLEANUP_20260816.md",
    "recovery_index": "archive/README.md",
}


def load_status(*, archive_only: bool = False) -> dict:
    if archive_only:
        return {"cleanup_record": STATE["cleanup_record"],
                "recovery_index": STATE["recovery_index"],
                "policy": "recovery records are not active entry points"}
    return dict(STATE)


def _human_status(value: dict, *, archive_only: bool) -> str:
    if archive_only:
        return "\n".join((
            f"Cleanup: {value['cleanup_record']}",
            f"Recovery: {value['recovery_index']}",
        ))
    stage = value["current_stage"]
    return "\n".join((
        f"Task: {value['core_task']['name']}",
        f"Stage: {stage['name']}",
        f"Status: {stage['status']}",
        f"Development training authorized: {str(stage['training_authorized']).lower()}",
        f"Unresolved items: {len(value['unresolved'])}",
        f"Recovery: {value['recovery_index']}",
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
