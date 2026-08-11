"""Compile the frozen R0-C downloads with the P1B holo contract."""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path

from scripts.build_holo_complex_index import (
    _compile_candidate_worker,
    record_set_sha256,
)
from scripts.data_contract import read_jsonl, write_jsonl
from scripts.structure_sources.biolip import BioLiPEntry
from scripts.structure_sources.rcsb import sha256_file


def _entry(row: dict) -> BioLiPEntry:
    return BioLiPEntry(**{key: row[key] for key in (
        "pdb_id", "receptor_auth_asym_id", "resolution", "binding_site_id",
        "ligand_comp_id", "ligand_auth_asym_id", "ligand_serial",
        "ligand_auth_seq_id", "sequence",
    )})


def build(
    panel_path: str | Path,
    raw_root: str | Path,
    acquisition_manifest_path: str | Path,
    output_dir: str | Path,
    *,
    workers: int = 1,
) -> dict:
    panel_file = Path(panel_path)
    root = Path(raw_root)
    acquisition_path = Path(acquisition_manifest_path)
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"R0-C holo output already exists: {output}")
    if workers < 1:
        raise ValueError("workers must be positive")
    rows = read_jsonl(panel_file)
    row_by_entry = {str(row["source_entry_id"]): row for row in rows}
    if len(row_by_entry) != len(rows):
        raise ValueError("duplicate source_entry_id in R0-C panel")
    acquisition = json.loads(acquisition_path.read_text(encoding="utf-8"))
    if (
        acquisition.get("schema") != "MetaSieve.R0C.StructureAcquisition.v1"
        or acquisition.get("files_failed") != 0
        or acquisition.get("panel", {}).get("sha256") != sha256_file(panel_file)
    ):
        raise ValueError("R0-C acquisition manifest does not bind a complete panel")

    payloads = [(_entry(row), str(root)) for row in rows]
    if workers == 1:
        compiled = map(_compile_candidate_worker, payloads)
        pool = None
    else:
        pool = ProcessPoolExecutor(max_workers=workers)
        compiled = pool.map(_compile_candidate_worker, payloads, chunksize=8)
    accepted, exclusions = [], []
    try:
        for entry, record, reason in compiled:
            if record is not None and float(record["protein_mapping_coverage"]) < 0.999999:
                record, reason = None, "protein_mapping_not_exact"
            if record is None:
                exclusions.append({"source_entry_id": entry.source_entry_id, "reason": reason})
                continue
            source = row_by_entry[entry.source_entry_id]
            accepted.append(record | {
                "r0c_component_id": source["r0c_component_id"],
                "r0_split": "heldout_b",
            })
    finally:
        if pool is not None:
            pool.shutdown()
    accepted.sort(key=lambda row: row["source_entry_id"])
    exclusions.sort(key=lambda row: row["source_entry_id"])
    output.mkdir(parents=True)
    write_jsonl(output / "complexes.jsonl", accepted)
    write_jsonl(output / "exclusions.jsonl", exclusions)
    manifest = {
        "schema": "MetaSieve.R0C.HoloIndex.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_records": len(rows),
        "accepted_records": len(accepted),
        "excluded_records": len(exclusions),
        "exclusions": dict(Counter(row["reason"] for row in exclusions)),
        "records_sha256": record_set_sha256(accepted),
        "complexes_sha256": sha256_file(output / "complexes.jsonl"),
        "panel": {"path": str(panel_file.resolve()), "sha256": sha256_file(panel_file)},
        "acquisition_manifest": {
            "path": str(acquisition_path.resolve()),
            "sha256": sha256_file(acquisition_path),
        },
        "affinity_value_reads": 0,
        "distance_value_reads": 0,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("panel")
    parser.add_argument("raw_root")
    parser.add_argument("acquisition_manifest")
    parser.add_argument("output")
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    print(json.dumps(build(
        args.panel, args.raw_root, args.acquisition_manifest, args.output,
        workers=args.workers,
    ), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
