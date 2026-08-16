"""Audit the legal coordinate coverage of the active DTA task.

Answers one question with numbers rather than assumption: for how many
BindingDB Ki deployment pairs does a common-frame protein-ligand complex exist
in the governed holo corpus?  A Cartesian protein-ligand encoder is only
defensible where that number is non-trivial.
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


ROOT = Path(__file__).resolve().parents[1]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
HOLO = ROOT / "dataset/processed/open_structures/pilot20k_holo_governed_v2"
SUPERVISION = ROOT / "dataset/processed/open_structures/pilot20k_structure_supervision_v2"


def read_jsonl(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    holo = read_jsonl(HOLO / "complexes.jsonl")
    cells = read_jsonl(CORPUS / "cells.jsonl.gz")
    proteins = read_jsonl(CORPUS / "proteins.jsonl")
    ligands = read_jsonl(CORPUS / "ligands.jsonl")
    smiles_of = {row["drug_key"]: row.get("smiles") for row in ligands}

    holo_sequences = {row["sequence"] for row in holo}
    holo_sequence_sha = {row["sequence_sha256"] for row in holo}
    holo_smiles = {row["canonical_smiles"] for row in holo}
    holo_pairs = {(row["sequence_sha256"], row["canonical_smiles"]) for row in holo}

    dta_sequence_sha = {row["sequence_sha256"] for row in proteins}
    dta_smiles = {value for value in smiles_of.values() if value}

    holo_sequence_list = list(holo_sequences)
    containment = sum(
        any(row["sequence"] in other or other in row["sequence"]
            for other in holo_sequence_list)
        for row in proteins)

    joint_cells = 0
    joint_targets = set()
    for cell in cells:
        key = (cell["target_id"], smiles_of.get(cell["ligand_id"]))
        if key in holo_pairs:
            joint_cells += 1
            joint_targets.add(cell["target_id"])

    supervision = json.loads((SUPERVISION / "manifest.json").read_text(encoding="utf-8"))

    report = {
        "schema": "MetaSieve.DTAGeometryCoverageAudit.v1",
        "holo_corpus": {
            "path": str(HOLO),
            "complexes": len(holo),
            "unique_receptor_sequences": len(holo_sequence_sha),
            "unique_ligand_chemotypes": len(holo_smiles),
            "raw_coordinates_available": True,
            "raw_coordinate_format": "gzipped mmCIF under dataset/raw/open_structures/pilot20k/mmcif",
        },
        "processed_structure_supervision": {
            "path": str(SUPERVISION),
            "pairs": supervision.get("pairs"),
            "tensors": ["contact", "distance_bin", "atom_mask", "residue_mask"],
            "contains_cartesian_coordinates": False,
            "note": ("only rotation/translation-invariant contact and distance "
                     "bins are materialized; coordinates are not"),
        },
        "dta_corpus": {
            "path": str(CORPUS),
            "cells": len(cells),
            "targets": len(dta_sequence_sha),
            "ligands": len(dta_smiles),
        },
        "overlap": {
            "targets_with_exact_holo_sequence": len(dta_sequence_sha & holo_sequence_sha),
            "targets_with_containment_holo_sequence": containment,
            "ligands_with_exact_holo_smiles": len(dta_smiles & holo_smiles),
            "cells_with_common_frame_complex": joint_cells,
            "targets_with_common_frame_complex": len(joint_targets),
        },
        "verdict": (
            "No BindingDB Ki deployment pair has a solved common-frame "
            "protein-ligand complex. A joint Cartesian protein-ligand encoder "
            "has no legal input on this task; claiming Cartesian equivariance "
            "or atomic 3D recognition for the DTA path would require "
            "manufacturing cross-molecular geometry."),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(report["overlap"], indent=2, sort_keys=True))
    print(report["verdict"])


if __name__ == "__main__":
    main()
