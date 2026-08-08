"""Label-blind feasibility audit for the frozen E0 source corpus."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
from statistics import median

from contracts.ligand_graph import MAX_ATOMS
from scripts.structure_sources.rcsb import sha256_file


AFFINITY_VALUE_FIELDS = frozenset({
    "p_affinity", "pchembl_value_reported", "published_value", "standard_value",
})
VALID_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _without_affinity_values(pairs):
    return {key: value for key, value in pairs if key not in AFFINITY_VALUE_FIELDS}


def _read_jsonl(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line, object_pairs_hook=_without_affinity_values)


def audit_e0_inputs(rows_path: Path, splits_path: Path) -> dict:
    from rdkit import Chem, rdBase

    assignments = list(_read_jsonl(splits_path))
    by_task = {row["task_id"]: row for row in assignments}
    if len(by_task) != len(assignments):
        raise ValueError("duplicate task assignment")

    row_counts = Counter()
    valid_row_counts = Counter()
    invalid_reasons = Counter()
    fold_rows = Counter()
    fold_valid_rows = Counter()
    fold_ligands: dict[int, set[str]] = defaultdict(set)
    fold_proteins: dict[int, set[str]] = defaultdict(set)
    task_valid_ligands: dict[str, set[str]] = defaultdict(set)
    ligands: dict[str, tuple[str, int]] = {}
    proteins: dict[str, str] = {}

    for row in _read_jsonl(rows_path):
        assignment = by_task.get(row.get("task_id"))
        if assignment is None:
            continue
        task_id = row["task_id"]
        fold = int(assignment["outer_oof_fold"])
        ligand_key = row["ligand_connectivity_key"]
        ligand_state_key = row.get("standard_inchi_key") or ligand_key
        protein_key = row["protein_sequence_sha256"]
        smiles = row["canonical_smiles"]
        sequence = row["protein_sequence"]
        row_counts[task_id] += 1
        fold_rows[fold] += 1

        molecule = Chem.MolFromSmiles(smiles)
        reason = None
        if molecule is None:
            reason = "invalid_smiles"
        elif molecule.GetNumAtoms() > MAX_ATOMS:
            reason = "atom_count_exceeds_128"
        elif not sequence or set(sequence) - VALID_RESIDUES:
            reason = "unsupported_protein_residue"
        if reason is not None:
            invalid_reasons[reason] += 1
            continue

        atom_count = molecule.GetNumAtoms()
        previous_ligand = ligands.setdefault(ligand_state_key, (smiles, atom_count))
        if previous_ligand != (smiles, atom_count):
            raise ValueError(f"inconsistent canonical ligand state {ligand_state_key}")
        previous_sequence = proteins.setdefault(protein_key, sequence)
        if previous_sequence != sequence:
            raise ValueError(f"inconsistent protein sequence {protein_key}")
        if assignment["protein_sequence_sha256"] != protein_key:
            raise ValueError(f"task/protein assignment mismatch {task_id}")

        valid_row_counts[task_id] += 1
        fold_valid_rows[fold] += 1
        fold_ligands[fold].add(ligand_state_key)
        fold_proteins[fold].add(protein_key)
        task_valid_ligands[task_id].add(ligand_key)

    missing_tasks = sorted(set(by_task) - set(row_counts))
    retained_tasks = {
        task_id for task_id, values in task_valid_ligands.items() if len(values) >= 20
    }
    retained_folds = Counter(int(by_task[task]["outer_oof_fold"])
                             for task in retained_tasks)
    retained_components = {by_task[task]["closure_component_id"] for task in retained_tasks}
    retained_proteins = {by_task[task]["protein_sequence_sha256"] for task in retained_tasks}
    atom_counts = sorted(value[1] for value in ligands.values())
    folds = {
        str(fold): {
            "tasks": sum(1 for row in assignments if int(row["outer_oof_fold"]) == fold),
            "rows": fold_rows[fold],
            "valid_rows": fold_valid_rows[fold],
            "unique_ligands": len(fold_ligands[fold]),
            "unique_proteins": len(fold_proteins[fold]),
            "retained_tasks_at_20_valid_compounds": retained_folds[fold],
        }
        for fold in range(5)
    }
    data_floor_pass = (
        len(retained_tasks) >= 2500 and len(retained_proteins) >= 500
        and len(retained_components) >= 200
        and min(retained_folds.get(fold, 0) for fold in range(5)) >= 500
    )
    return {
        "schema": "MetaSieve.E0InputFeasibility.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P1R2B-E0",
        "audit_mode": "LABEL_BLIND",
        "affinity_value_fields_materialized": False,
        "recipient_labels_read": False,
        "training_performed": False,
        "contract": {"max_atoms": MAX_ATOMS, "valid_residues": "ACDEFGHIKLMNPQRSTVWY"},
        "inputs": {
            "canonical_rows_sha256": sha256_file(rows_path),
            "split_assignments_sha256": sha256_file(splits_path),
            "rdkit_version": rdBase.rdkitVersion,
        },
        "counts": {
            "governed_tasks": len(by_task),
            "governed_rows": sum(row_counts.values()),
            "valid_rows": sum(valid_row_counts.values()),
            "unique_valid_ligands": len(ligands),
            "unique_valid_proteins": len(proteins),
            "missing_tasks": len(missing_tasks),
            "invalid_rows_by_reason": dict(sorted(invalid_reasons.items())),
        },
        "atom_counts": {
            "minimum": atom_counts[0] if atom_counts else None,
            "median": median(atom_counts) if atom_counts else None,
            "maximum": atom_counts[-1] if atom_counts else None,
        },
        "folds": folds,
        "post_contract_floor": {
            "tasks_with_at_least_20_valid_compounds": len(retained_tasks),
            "proteins": len(retained_proteins),
            "closure_components": len(retained_components),
            "minimum_tasks_per_fold": min(retained_folds.get(fold, 0) for fold in range(5)),
            "pass": data_floor_pass,
        },
        "missing_task_ids": missing_tasks,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    pilot = Path("dataset/processed/source_affinity/energy_pilot_v1")
    governance = Path("dataset/processed/source_affinity/energy_pilot_v1_governance")
    parser.add_argument("--rows", type=Path, default=pilot / "canonical_rows.jsonl")
    parser.add_argument("--splits", type=Path, default=governance / "split_assignments.jsonl")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit_e0_inputs(args.rows, args.splits)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["post_contract_floor"]["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
