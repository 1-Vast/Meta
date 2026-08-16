"""Read-only status view over the governed project and archive manifests."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "archive/retired_research_20260811/ARCHIVE_MANIFEST.json"

STATE = {
    "schema_version": "6.0",
    "updated": "2026-08-15",
    "core_task": {
        "name": "cold-target few-shot drug-target affinity prediction",
        "target_is_meta_task": True,
        "support_sizes": [0, 1, 2, 3, 5],
        "query_labels_forbidden": True,
    },
    "current_stage": {
        "name": "QPSMP_BPSF_DEVELOPMENT",
        "status": "TRAINABLE_INTERFACE_NOT_ADMITTED",
        "training_authorized": True,
        "confirmation_labels_open": False,
    },
    "unresolved": [
        "PREREGISTERED_COMPONENT_LEVEL_DEVELOPMENT_GAIN",
        "CORRECT_PROTEIN_AND_SAR_SPECIFICITY",
        "UNTOUCHED_CONFIRMATION_COMPONENTS",
        "CARTESIAN_PERFORMANCE_WITH_LEGAL_COMMON_FRAME_COORDINATES",
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
        f"Development training authorized: {str(stage['training_authorized']).lower()}",
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
