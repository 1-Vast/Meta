"""Apply frozen chemistry quarantine to the exact-mapped R0-C records."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

from scripts.data_contract import read_jsonl, write_jsonl
from scripts.structure_sources.rcsb import sha256_file


MIN_COMPONENTS = 120
FINAL_NAMESPACE = "R0C-FINAL-REP-v1"


def _chemistry_keys(row: dict) -> set[tuple[str, str]]:
    keys = {
        ("ccd", str(row["ccd_sha256"])),
        ("connectivity", str(row["connectivity_sha256"])),
    }
    scaffold = str(row.get("murcko_scaffold", ""))
    if scaffold:
        keys.add(("scaffold", scaffold))
    return keys


def finalize(records: list[dict], exposure: list[dict]) -> tuple[list[dict], list[dict], dict]:
    exposure_keys = set().union(*(_chemistry_keys(row) for row in exposure)) if exposure else set()
    clean, exclusions = [], []
    for row in records:
        overlap = sorted(_chemistry_keys(row) & exposure_keys)
        if overlap:
            exclusions.append({
                "source_entry_id": row["source_entry_id"],
                "reason": "+".join(sorted({kind for kind, _ in overlap})),
            })
        else:
            clean.append(row)

    parent = {row["source_entry_id"]: row["source_entry_id"] for row in clean}
    def find(value: str) -> str:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value
    def union(left: str, right: str) -> None:
        left, right = find(left), find(right)
        if left != right:
            low, high = sorted((left, right))
            parent[high] = low
    owners: dict[tuple[str, str], str] = {}
    for row in clean:
        entry = row["source_entry_id"]
        for key in _chemistry_keys(row):
            if key in owners:
                union(entry, owners[key])
            else:
                owners[key] = entry
    components: dict[str, list[dict]] = {}
    for row in clean:
        components.setdefault(find(row["source_entry_id"]), []).append(row)
    selected = []
    for rows in components.values():
        chosen = min(rows, key=lambda row: hashlib.sha256(
            f"{FINAL_NAMESPACE}|{row['source_entry_id']}".encode("utf-8")
        ).hexdigest())
        selected.append(chosen | {"r0c_final_component_id": min(
            row["source_entry_id"] for row in rows
        )})
        for row in rows:
            if row["source_entry_id"] != chosen["source_entry_id"]:
                exclusions.append({
                    "source_entry_id": row["source_entry_id"],
                    "reason": "internal_chemistry_component_nonrepresentative",
                })
    selected.sort(key=lambda row: row["source_entry_id"])
    exclusions.sort(key=lambda row: row["source_entry_id"])
    audit = {
        "schema": "MetaSieve.R0C.FinalPanel.v1",
        "input_records": len(records),
        "exposure_records": len(exposure),
        "chemistry_clean_records": len(clean),
        "selected_records": len(selected),
        "selected_components": len(selected),
        "minimum_components": MIN_COMPONENTS,
        "largest_component_share": 1.0 / max(1, len(selected)),
        "component_gate_pass": len(selected) >= MIN_COMPONENTS,
        "exclusions": dict(Counter(row["reason"] for row in exclusions)),
        "affinity_value_reads": 0,
        "distance_value_reads": 0,
    }
    return selected, exclusions, audit


def prepare(records_path: str | Path, exposure_path: str | Path, output_dir: str | Path) -> dict:
    records_file, exposure_file = Path(records_path), Path(exposure_path)
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"R0-C final panel output already exists: {output}")
    selected, exclusions, audit = finalize(
        read_jsonl(records_file), read_jsonl(exposure_file)
    )
    output.mkdir(parents=True)
    write_jsonl(output / "panel.jsonl", selected)
    write_jsonl(output / "exclusions.jsonl", exclusions)
    audit = {
        **audit,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "records_input": {"path": str(records_file.resolve()), "sha256": sha256_file(records_file)},
        "exposure_input": {"path": str(exposure_file.resolve()), "sha256": sha256_file(exposure_file)},
        "panel_sha256": sha256_file(output / "panel.jsonl"),
    }
    (output / "panel_audit.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if not audit["component_gate_pass"]:
        raise RuntimeError("R0-C final panel failed the frozen component Gate")
    return audit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("records")
    parser.add_argument("exposure")
    parser.add_argument("output")
    args = parser.parse_args()
    print(json.dumps(prepare(args.records, args.exposure, args.output), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
