"""Independently verify EnergyPilot rows, task reconstruction, and hash chain."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.source_affinity.canonicalize import p_affinity, stable_hash, task_manifest
from scripts.source_affinity.common import canonical_json_bytes, sha256_file, write_canonical_json


def verified_rows(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            row = json.loads(line)
            prefix = f"row {line_number} activity {row.get('activity_id')}"
            if row.get("schema") != "MetaSieve.AffinityEnergyRow.v1":
                raise ValueError(f"{prefix}: schema mismatch")
            if row.get("chembl_release") != "37":
                raise ValueError(f"{prefix}: release mismatch")
            if row.get("endpoint_family") not in {"Ki", "Kd"} or row.get("standard_relation") != "=":
                raise ValueError(f"{prefix}: E0-Core endpoint/relation mismatch")
            sequence = row["protein_sequence"]
            if stable_hash(sequence) != row["protein_sequence_sha256"]:
                raise ValueError(f"{prefix}: protein hash mismatch")
            context_json = json.dumps(
                row["assay_context"], sort_keys=True, separators=(",", ":"), default=str,
            )
            if stable_hash(context_json) != row["assay_context_sha256"]:
                raise ValueError(f"{prefix}: context hash mismatch")
            expected_p = p_affinity(row["standard_value"], row["standard_units"])
            if not math.isclose(expected_p, row["p_affinity"], rel_tol=0, abs_tol=1e-12):
                raise ValueError(f"{prefix}: pAffinity mismatch")
            task_keys = [
                row["protein_sequence_sha256"], row["endpoint_family"],
                row["assay_chembl_id"], row["assay_context_sha256"],
            ]
            if task_keys != row["task_keys"]:
                raise ValueError(f"{prefix}: task keys mismatch")
            expected_task = hashlib.sha256(
                json.dumps(task_keys, separators=(",", ":")).encode()
            ).hexdigest()
            if expected_task != row["task_id"]:
                raise ValueError(f"{prefix}: task hash mismatch")
            if row["ligand_connectivity_key"] != row["standard_inchi_key"].split("-")[0]:
                raise ValueError(f"{prefix}: ligand connectivity mismatch")
            yield row


def verify(
    rows_path: Path, tasks_path: Path, corpus_manifest_path: Path,
    release_manifest_path: Path, sql_path: Path, schema_path: Path,
) -> dict:
    corpus = json.loads(corpus_manifest_path.read_text(encoding="utf-8"))
    release = json.loads(release_manifest_path.read_text(encoding="utf-8"))
    rows = verified_rows(rows_path)
    reconstructed_tasks = task_manifest(rows)
    reconstructed_task_hash = hashlib.sha256(canonical_json_bytes(reconstructed_tasks)).hexdigest()
    checks = {
        "canonical_rows_sha256": sha256_file(rows_path) == corpus["canonical_rows"]["sha256"],
        "canonical_row_count": sum(1 for _ in rows_path.open(encoding="utf-8"))
        == corpus["canonical_rows"]["rows"],
        "task_manifest_reconstruction_sha256": reconstructed_task_hash == sha256_file(tasks_path),
        "task_manifest_registered_sha256": sha256_file(tasks_path) == corpus["tasks"]["sha256"],
        "task_count": len(reconstructed_tasks) == corpus["tasks"]["count"],
        "sql_sha256": sha256_file(sql_path) == corpus["sql_sha256"],
        "normalizer_sha256": sha256_file(
            ROOT / "scripts/source_affinity/canonicalize.py"
        ) == corpus["normalizer_sha256"],
        "row_schema_sha256": sha256_file(schema_path) == corpus["row_schema_sha256"],
        "official_archive_checksum_verified": release["official_checksum_verified"] is True,
        "release_training_authorized_false": release["training_authorized"] is False,
        "recipient_labels_unread": (
            release["recipient_labels_read"] is False
            and corpus["recipient_labels_read"] is False
        ),
    }
    if not all(checks.values()):
        raise RuntimeError(f"corpus verification failed: {checks}")
    return {
        "schema": "MetaSieve.AffinityEnergyCorpusVerification.v1",
        "stage": "P1R2B-D0-C",
        "verdict": "PASS",
        "checks": checks,
        "canonical_rows": corpus["canonical_rows"]["rows"],
        "tasks": corpus["tasks"]["count"],
        "eligible_e0_core_tasks": corpus["tasks"]["eligible_e0_core"],
        "inputs": {
            "release_manifest_sha256": sha256_file(release_manifest_path),
            "corpus_manifest_sha256": sha256_file(corpus_manifest_path),
            "canonical_rows_sha256": sha256_file(rows_path),
            "task_manifest_sha256": sha256_file(tasks_path),
        },
        "recipient_labels_read": False,
        "training_authorized": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    pilot = ROOT / "dataset/processed/source_affinity/energy_pilot_v1"
    parser.add_argument("--rows", type=Path, default=pilot / "canonical_rows.jsonl")
    parser.add_argument("--tasks", type=Path, default=pilot / "task_manifest.json")
    parser.add_argument("--corpus-manifest", type=Path, default=pilot / "corpus_manifest.json")
    parser.add_argument(
        "--release-manifest", type=Path,
        default=ROOT / "dataset/raw/source_affinity/chembl37_sqlite_v1/release_manifest.json",
    )
    parser.add_argument(
        "--sql", type=Path,
        default=ROOT / "contracts/source_affinity/chembl37_e0_core.sql",
    )
    parser.add_argument(
        "--row-schema", type=Path,
        default=ROOT / "contracts/source_affinity/affinity_energy_row_v1.json",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(
        args.rows, args.tasks, args.corpus_manifest, args.release_manifest,
        args.sql, args.row_schema,
    )
    write_canonical_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
