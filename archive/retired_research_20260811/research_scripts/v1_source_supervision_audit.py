"""Census measured source supervision for a partner-specific MetaSieve v1."""
from __future__ import annotations

import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
LABELS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/exact_labels.jsonl.gz"
OUT = ROOT / "report/meta_fewshot/v1_source_supervision_audit.json"
ALLOWED_SPLITS = ("meta_train", "meta_val")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_gzip_jsonl(path: Path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def source_key(row_id: str, target_id: str, ligand_id: str) -> tuple[str, str, str]:
    return str(row_id), target_id, ligand_id


def admitted_source_rows(cells: list[dict]) -> tuple[dict[tuple[str, str, str], str], set[tuple[str, str, str]]]:
    allowed, forbidden = {}, set()
    for cell in cells:
        destination = allowed if cell["split"] in ALLOWED_SPLITS else forbidden
        for row_id in cell["source_row_ids"]:
            key = source_key(row_id, cell["target_id"], cell["ligand_id"])
            if destination is allowed:
                previous = allowed.setdefault(key, cell["split"])
                if previous != cell["split"]:
                    raise ValueError("one source row is assigned to two development splits")
            else:
                forbidden.add(key)
    if set(allowed) & forbidden:
        raise ValueError("development and meta-test source rows overlap")
    return allowed, forbidden


def median_rows(rows: list[dict]) -> list[dict]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["split"], row["panel_id"], row["target_id"], row["ligand_id"])].append(
            float(row["pK"]))
    return [
        {"split": key[0], "panel_id": key[1], "target_id": key[2],
         "ligand_id": key[3], "pK": float(np.median(values)),
         "replicates": len(values)}
        for key, values in sorted(grouped.items())
    ]


def summarize(rows: list[dict], target_group: dict[str, str], split: str) -> dict:
    selected = [row for row in rows if row["split"] == split]
    within = defaultdict(list)
    partner = defaultdict(list)
    for row in selected:
        within[(row["panel_id"], row["target_id"])].append(row)
        partner[(row["panel_id"], row["ligand_id"])].append(row)

    within_eligible = [values for values in within.values()
                       if len({row["ligand_id"] for row in values}) >= 2]
    partner_eligible = []
    differences = []
    for values in partner.values():
        targets = {row["target_id"] for row in values}
        groups = {target_group[target] for target in targets}
        if len(targets) < 2 or len(groups) < 2:
            continue
        partner_eligible.append(values)
        for left in range(len(values) - 1):
            for right in range(left + 1, len(values)):
                if target_group[values[left]["target_id"]] != target_group[values[right]["target_id"]]:
                    differences.append(abs(values[left]["pK"] - values[right]["pK"]))

    within_targets = {row["target_id"] for values in within_eligible for row in values}
    partner_targets = {row["target_id"] for values in partner_eligible for row in values}
    partner_ligands = {row["ligand_id"] for values in partner_eligible for row in values}
    partner_groups = {target_group[target] for target in partner_targets}
    return {
        "panel_target_groups": len(within),
        "within_panel_groups_ge2_ligands": len(within_eligible),
        "within_panel_targets": len(within_targets),
        "within_panel_cdhit40_groups": len({target_group[target] for target in within_targets}),
        "within_panel_observations": sum(len(values) for values in within_eligible),
        "panel_ligand_groups": len(partner),
        "measured_partner_groups_cross_cdhit40": len(partner_eligible),
        "measured_partner_targets": len(partner_targets),
        "measured_partner_ligands": len(partner_ligands),
        "measured_partner_cdhit40_groups": len(partner_groups),
        "measured_partner_observations": sum(len(values) for values in partner_eligible),
        "measured_cross_family_pairs": len(differences),
        "absolute_delta_pK": {
            "median": float(np.median(differences)) if differences else None,
            "q25": float(np.quantile(differences, 0.25)) if differences else None,
            "q75": float(np.quantile(differences, 0.75)) if differences else None,
        },
    }


def audit(corpus: Path = CORPUS, labels: Path = LABELS) -> dict:
    cells = list(read_gzip_jsonl(corpus / "cells.jsonl.gz"))
    target_group = {row["target_id"]: row["protein_group_40"] for row in cells}
    allowed, forbidden = admitted_source_rows(cells)
    source_rows = []
    forbidden_seen = 0
    for row in read_gzip_jsonl(labels):
        if row["endpoint"] != "Ki":
            continue
        key = source_key(row["source_row_id"], row["target_id"], row["ligand_id"])
        if key in forbidden:
            forbidden_seen += 1
        split = allowed.get(key)
        if split is None:
            continue
        source_rows.append({**row, "split": split})
    medians = median_rows(source_rows)
    return {
        "schema": "MetaSieve.V1SourceSupervisionAudit.v1",
        "corpus_manifest_sha256": sha256(corpus / "manifest.json"),
        "exact_labels_sha256": sha256(labels),
        "allowed_splits": list(ALLOWED_SPLITS),
        "meta_test_values_used": 0,
        "meta_test_rows_encountered_but_rejected": forbidden_seen,
        "source_rows": len(source_rows),
        "panel_target_ligand_medians": len(medians),
        "replicated_panel_target_ligand_groups": sum(row["replicates"] > 1 for row in medians),
        "splits": {split: summarize(medians, target_group, split) for split in ALLOWED_SPLITS},
    }


def main() -> int:
    result = audit()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
