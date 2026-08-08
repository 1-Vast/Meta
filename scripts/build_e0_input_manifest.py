"""Build a label-free, model-compatible E0 input manifest."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

from rdkit import Chem, rdBase

from contracts.ligand_graph import MAX_ATOMS
from scripts.audit_e0_input_feasibility import _read_jsonl, VALID_RESIDUES
from scripts.structure_sources.rcsb import sha256_file


def _write_jsonl(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def build_e0_input_manifest(rows_path: Path, splits_path: Path, output: Path) -> dict:
    assignments = {row["task_id"]: row for row in _read_jsonl(splits_path)}
    if not assignments:
        raise ValueError("no governed task assignments")

    candidate_rows = []
    task_compounds: dict[str, set[str]] = defaultdict(set)
    ligands: dict[str, dict] = {}
    proteins: dict[str, dict] = {}
    for row in _read_jsonl(rows_path):
        assignment = assignments.get(row.get("task_id"))
        if assignment is None:
            continue
        molecule = Chem.MolFromSmiles(row["canonical_smiles"])
        sequence = row["protein_sequence"]
        if (molecule is None or molecule.GetNumAtoms() > MAX_ATOMS or not sequence
                or set(sequence) - VALID_RESIDUES):
            continue
        state_key = row.get("standard_inchi_key") or row["ligand_connectivity_key"]
        ligand = {
            "atom_count": molecule.GetNumAtoms(),
            "canonical_smiles": row["canonical_smiles"],
            "ligand_connectivity_key": row["ligand_connectivity_key"],
            "ligand_state_key": state_key,
        }
        if state_key in ligands and ligands[state_key] != ligand:
            raise ValueError(f"inconsistent ligand state {state_key}")
        ligands[state_key] = ligand
        protein_key = row["protein_sequence_sha256"]
        protein = {"sequence": sequence, "sequence_sha256": protein_key}
        if protein_key in proteins and proteins[protein_key] != protein:
            raise ValueError(f"inconsistent protein state {protein_key}")
        proteins[protein_key] = protein
        task_compounds[row["task_id"]].add(row["ligand_connectivity_key"])
        candidate_rows.append({
            "activity_id": row["activity_id"],
            "closure_component_id": assignment["closure_component_id"],
            "endpoint_family": row["endpoint_family"],
            "ligand_connectivity_key": row["ligand_connectivity_key"],
            "ligand_state_key": state_key,
            "outer_oof_fold": int(assignment["outer_oof_fold"]),
            "protein_sequence_sha256": protein_key,
            "task_id": row["task_id"],
        })

    retained_tasks = {task for task, compounds in task_compounds.items()
                      if len(compounds) >= 20}
    retained_rows = [row for row in candidate_rows if row["task_id"] in retained_tasks]
    retained_ligand_keys = {row["ligand_state_key"] for row in retained_rows}
    retained_protein_keys = {row["protein_sequence_sha256"] for row in retained_rows}
    ligand_rows = [ligands[key] for key in sorted(retained_ligand_keys)]
    protein_rows = [proteins[key] for key in sorted(retained_protein_keys)]
    retained_rows.sort(key=lambda row: (row["outer_oof_fold"], row["task_id"],
                                        row["ligand_state_key"], row["activity_id"]))

    output.mkdir(parents=True, exist_ok=True)
    row_output = output / "rows.label_blind.jsonl"
    ligand_output = output / "ligands.jsonl"
    protein_output = output / "proteins.jsonl"
    _write_jsonl(row_output, retained_rows)
    _write_jsonl(ligand_output, ligand_rows)
    _write_jsonl(protein_output, protein_rows)
    manifest = {
        "schema": "MetaSieve.E0InputManifest.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P1R2B-E0",
        "label_blind": True,
        "affinity_values_present": False,
        "recipient_labels_read": False,
        "training_performed": False,
        "contract": {"max_atoms": MAX_ATOMS, "minimum_compounds_per_task": 20},
        "inputs": {"canonical_rows_sha256": sha256_file(rows_path),
                   "split_assignments_sha256": sha256_file(splits_path),
                   "builder_sha256": sha256_file(Path(__file__)),
                   "rdkit_version": rdBase.rdkitVersion},
        "counts": {"rows": len(retained_rows), "tasks": len(retained_tasks),
                   "ligand_states": len(ligand_rows), "proteins": len(protein_rows)},
        "outputs": {"rows": sha256_file(row_output), "ligands": sha256_file(ligand_output),
                    "proteins": sha256_file(protein_output)},
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    pilot = Path("dataset/processed/source_affinity/energy_pilot_v1")
    governance = Path("dataset/processed/source_affinity/energy_pilot_v1_governance")
    parser.add_argument("--rows", type=Path, default=pilot / "canonical_rows.jsonl")
    parser.add_argument("--splits", type=Path, default=governance / "split_assignments.jsonl")
    parser.add_argument("--output", type=Path,
                        default=Path("dataset/processed/source_affinity/e0_input_v1"))
    args = parser.parse_args()
    print(json.dumps(build_e0_input_manifest(args.rows, args.splits, args.output),
                     indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
