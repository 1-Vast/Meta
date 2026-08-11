"""Prepare the label-free governed P1B panel for the R0-B distance audit."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

from scripts.data_contract import read_jsonl, write_jsonl
from scripts.structure_sources.rcsb import sha256_file


TRAIN_NAMESPACE = "R0B-TRAIN-REP-v1"
MIN_HELDOUT_COMPONENTS = 30
MAX_HELDOUT_COMPONENT_SHARE = 0.20


def _selection_key(record: dict) -> str:
    return hashlib.sha256(
        f"{TRAIN_NAMESPACE}|{record['source_entry_id']}".encode("utf-8")
    ).hexdigest()


def _chemical_key(record: dict) -> tuple[str, str]:
    return str(record["canonical_smiles"]), str(record["murcko_scaffold"])


def select_panel(records: list[dict]) -> tuple[list[dict], list[dict], dict]:
    """Select one train representative per family, then filter later splits.

    The operation reads no geometry or affinity value.  Heldout records are
    never used to choose among train representatives.
    """
    required = {
        "source_entry_id", "source_split", "homology_group_id", "pdb_id",
        "sequence", "sequence_sha256", "protein_mapping_coverage",
        "canonical_smiles", "murcko_scaffold", "ccd_sha256",
        "connectivity_sha256",
    }
    if not records:
        raise ValueError("R0-B panel selection received no records")
    if any(required - set(record) for record in records):
        raise ValueError("one or more governed records lack R0-B fields")

    eligible, exclusions = [], []
    seen_entries = set()
    for record in records:
        reason = None
        if record["source_entry_id"] in seen_entries:
            raise ValueError("duplicate source_entry_id in governed records")
        seen_entries.add(record["source_entry_id"])
        if record["source_split"] not in {"train", "val", "test"}:
            reason = "unrecognized_source_split"
        elif float(record["protein_mapping_coverage"]) < 0.999999:
            reason = "protein_mapping_not_exact"
        elif not record["canonical_smiles"]:
            reason = "missing_exact_ligand_graph"
        if reason:
            exclusions.append({"source_entry_id": record["source_entry_id"],
                               "source_split": record["source_split"],
                               "reason": reason})
        else:
            eligible.append(record)

    by_train_group: dict[str, list[dict]] = defaultdict(list)
    for record in eligible:
        if record["source_split"] == "train":
            by_train_group[str(record["homology_group_id"])].append(record)
    train = [min(group, key=_selection_key) for group in by_train_group.values()]
    selected_train = {record["source_entry_id"] for record in train}
    for record in eligible:
        if record["source_split"] == "train" and \
                record["source_entry_id"] not in selected_train:
            exclusions.append({"source_entry_id": record["source_entry_id"],
                               "source_split": "train",
                               "reason": "not_hash_selected_family_representative"})

    def filter_later(rows: list[dict], reference: list[dict], split: str) -> list[dict]:
        ccd = {record["ccd_sha256"] for record in reference}
        exact = {record["connectivity_sha256"] for record in reference}
        scaffold = {record["murcko_scaffold"] for record in reference}
        kept = []
        for record in rows:
            reasons = []
            if record["ccd_sha256"] in ccd:
                reasons.append("earlier_split_ccd")
            if record["connectivity_sha256"] in exact:
                reasons.append("earlier_split_exact_connectivity")
            if record["murcko_scaffold"] in scaffold:
                reasons.append("earlier_split_scaffold")
            if reasons:
                exclusions.append({"source_entry_id": record["source_entry_id"],
                                   "source_split": split,
                                   "reason": "+".join(reasons)})
            else:
                kept.append(record)
        return kept

    train_exposure = [record for record in records if record["source_split"] == "train"]
    val_candidates = [record for record in eligible if record["source_split"] == "val"]
    # The frozen prior saw every P1B train/validation record, not only the
    # representatives later selected for residual fitting.  Chemical novelty
    # is therefore filtered against the full earlier checkpoint exposure.
    val = filter_later(val_candidates, train_exposure, "val")
    test_candidates = [record for record in eligible if record["source_split"] == "test"]
    val_exposure = [record for record in records if record["source_split"] == "val"]
    heldout = filter_later(test_candidates, train_exposure + val_exposure, "test")

    selected = []
    for output_split, rows in (("train", train), ("val", val),
                               ("heldout_a", heldout)):
        selected.extend([{**record, "r0_split": output_split} for record in rows])

    split_groups = {
        split: {record["homology_group_id"] for record in selected
                if record["r0_split"] == split}
        for split in ("train", "val", "heldout_a")
    }
    if any(split_groups[left] & split_groups[right]
           for left, right in (("train", "val"), ("train", "heldout_a"),
                               ("val", "heldout_a"))):
        raise RuntimeError("homology group straddles the R0-B split")
    for earlier, later in (("train", "val"), ("train", "heldout_a"),
                           ("val", "heldout_a")):
        left = [record for record in selected if record["r0_split"] == earlier]
        right = [record for record in selected if record["r0_split"] == later]
        left_ccd = {record["ccd_sha256"] for record in left}
        left_exact = {record["connectivity_sha256"] for record in left}
        left_scaffold = {record["murcko_scaffold"] for record in left}
        if any(record["ccd_sha256"] in left_ccd or
               record["connectivity_sha256"] in left_exact or
               record["murcko_scaffold"] in left_scaffold for record in right):
            raise RuntimeError("ligand dependency straddles the R0-B split")

    heldout_sizes = Counter(record["homology_group_id"] for record in heldout)
    largest = max(heldout_sizes.values(), default=0)
    largest_share = largest / max(1, len(heldout))
    split_counts = Counter(record["r0_split"] for record in selected)
    audit = {
        "schema": "MetaSieve.R0B.LabelFreePanel.v1",
        "selection_namespace": TRAIN_NAMESPACE,
        "input_records": len(records),
        "exact_mapping_eligible": len(eligible),
        "selected_records": dict(split_counts),
        "selected_homology_components": {
            split: len(groups) for split, groups in split_groups.items()},
        "heldout_largest_component_records": largest,
        "heldout_largest_component_share": largest_share,
        "heldout_component_requirements": {
            "minimum_components": MIN_HELDOUT_COMPONENTS,
            "maximum_share": MAX_HELDOUT_COMPONENT_SHARE,
        },
        "split_gate_pass": (
            len(split_groups["heldout_a"]) >= MIN_HELDOUT_COMPONENTS
            and largest_share < MAX_HELDOUT_COMPONENT_SHARE
        ),
        "exclusions": dict(Counter(row["reason"] for row in exclusions)),
        "affinity_value_reads": 0,
        "geometry_value_reads": 0,
        "p1b_exposure": {
            "train": "used to fit the frozen P1B checkpoint",
            "val": "used for frozen P1B checkpoint selection",
            "heldout_a": "previously consumed by the P1B Gate; R0 development only",
        },
    }
    return sorted(selected, key=lambda row: (row["r0_split"], row["source_entry_id"])), \
        sorted(exclusions, key=lambda row: (row["source_split"], row["source_entry_id"])), \
        audit


def prepare(records_path: str | Path, output_dir: str | Path) -> dict:
    records_path = Path(records_path)
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"R0-B panel output already exists: {output}")
    selected, exclusions, audit = select_panel(read_jsonl(records_path))
    if not audit["split_gate_pass"]:
        raise RuntimeError("R0-B label-free split failed its component Gate")
    output.mkdir(parents=True)
    write_jsonl(output / "panel.jsonl", selected)
    write_jsonl(output / "exclusions.jsonl", exclusions)
    audit = {
        **audit,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_path": str(records_path.resolve()),
        "input_sha256": sha256_file(records_path),
    }
    (output / "panel_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("records")
    parser.add_argument("output")
    args = parser.parse_args()
    print(json.dumps(prepare(args.records, args.output), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
