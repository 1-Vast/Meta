"""Select a fresh, coordinate-free BioLiP candidate pool for R0-C."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

from scripts.data_contract import write_jsonl
from scripts.structure_sources.biolip import (
    BioLiPEntry,
    iter_biolip,
    pilot_candidates,
    regular_ligand_ids,
)
from scripts.structure_sources.rcsb import sha256_file


P1B_CANDIDATE_LIMIT = 20_000
R0C_CANDIDATE_LIMIT = 3_000
SELECTION_NAMESPACE = "R0C-FRESH-CANDIDATE-v1"


def _coordinate_free_eligible(
    entry: BioLiPEntry, allowed_ligands: set[str]
) -> bool:
    return (
        0 < entry.resolution <= 3.0
        and 50 <= len(entry.sequence) <= 1022
        and entry.ligand_auth_seq_id.lstrip("-").isdigit()
        and entry.ligand_comp_id in allowed_ligands
    )


def _selection_key(entry: BioLiPEntry) -> str:
    return hashlib.sha256(
        f"{SELECTION_NAMESPACE}|{entry.source_entry_id}".encode("utf-8")
    ).hexdigest()


def select_fresh_candidates(
    entries: list[BioLiPEntry],
    old_candidates: list[BioLiPEntry],
    allowed_ligands: set[str],
    *,
    limit: int = R0C_CANDIDATE_LIMIT,
) -> tuple[list[BioLiPEntry], dict]:
    """Select metadata-novel entries without reading coordinates or labels."""
    if limit < 1:
        raise ValueError("R0-C candidate limit must be positive")
    old_entries = {entry.source_entry_id for entry in old_candidates}
    old_pdb = {entry.pdb_id for entry in old_candidates}
    old_sequences = {entry.sequence for entry in old_candidates}
    old_ligands = {entry.ligand_comp_id for entry in old_candidates}

    unique: dict[str, BioLiPEntry] = {}
    exclusions = {
        "not_coordinate_free_eligible": 0,
        "old_source_entry": 0,
        "old_pdb": 0,
        "old_exact_sequence": 0,
        "old_ligand_comp_id": 0,
    }
    for entry in entries:
        if not _coordinate_free_eligible(entry, allowed_ligands):
            exclusions["not_coordinate_free_eligible"] += 1
            continue
        if entry.source_entry_id in old_entries:
            exclusions["old_source_entry"] += 1
            continue
        if entry.pdb_id in old_pdb:
            exclusions["old_pdb"] += 1
            continue
        if entry.sequence in old_sequences:
            exclusions["old_exact_sequence"] += 1
            continue
        if entry.ligand_comp_id in old_ligands:
            exclusions["old_ligand_comp_id"] += 1
            continue
        unique.setdefault(entry.source_entry_id, entry)

    selected: list[BioLiPEntry] = []
    used_pdb: set[str] = set()
    used_sequences: set[str] = set()
    used_ligands: set[str] = set()
    for entry in sorted(unique.values(), key=_selection_key):
        if (
            entry.pdb_id in used_pdb
            or entry.sequence in used_sequences
            or entry.ligand_comp_id in used_ligands
        ):
            continue
        selected.append(entry)
        used_pdb.add(entry.pdb_id)
        used_sequences.add(entry.sequence)
        used_ligands.add(entry.ligand_comp_id)
        if len(selected) == limit:
            break
    if len(selected) != limit:
        raise RuntimeError(
            f"R0-C fresh candidate supply is {len(selected)}, below {limit}"
        )
    audit = {
        "schema": "MetaSieve.R0C.CoordinateFreeCandidates.v1",
        "selection_namespace": SELECTION_NAMESPACE,
        "requested_candidates": limit,
        "selected_candidates": len(selected),
        "selected_unique_pdb": len(used_pdb),
        "selected_unique_sequences": len(used_sequences),
        "selected_unique_ligand_comp_ids": len(used_ligands),
        "p1b_candidate_exposure": len(old_candidates),
        "eligible_after_old_exposure_filters": len(unique),
        "exclusions": exclusions,
        "affinity_value_reads": 0,
        "coordinate_value_reads": 0,
    }
    return selected, audit


def prepare(
    annotation_path: str | Path,
    ligand_summary_path: str | Path,
    output_dir: str | Path,
    *,
    limit: int = R0C_CANDIDATE_LIMIT,
) -> dict:
    annotation = Path(annotation_path)
    ligand_summary = Path(ligand_summary_path)
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"R0-C candidate output already exists: {output}")
    allowed = regular_ligand_ids(ligand_summary)
    old = pilot_candidates(
        annotation, limit=P1B_CANDIDATE_LIMIT, allowed_ligands=allowed
    )
    selected, audit = select_fresh_candidates(
        list(iter_biolip(annotation)), old, allowed, limit=limit
    )
    output.mkdir(parents=True)
    rows = [entry.to_dict() for entry in selected]
    write_jsonl(output / "candidates.jsonl", rows)
    audit = {
        **audit,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "annotation": {
            "path": str(annotation.resolve()),
            "sha256": sha256_file(annotation),
        },
        "ligand_summary": {
            "path": str(ligand_summary.resolve()),
            "sha256": sha256_file(ligand_summary),
        },
        "candidates_sha256": sha256_file(output / "candidates.jsonl"),
    }
    (output / "candidate_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("annotation")
    parser.add_argument("ligand_summary")
    parser.add_argument("output")
    parser.add_argument("--limit", type=int, default=R0C_CANDIDATE_LIMIT)
    args = parser.parse_args()
    print(
        json.dumps(
            prepare(
                args.annotation, args.ligand_summary, args.output, limit=args.limit
            ),
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
