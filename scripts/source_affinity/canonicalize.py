"""Deterministic ChEMBL E0-Core row canonicalization."""
from __future__ import annotations

from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import sqlite3
from statistics import median
from typing import Any, Iterable

from rdkit import Chem, rdBase

from scripts.source_affinity.common import sha256_file, write_canonical_json


UNIT_TO_MOLAR = {
    "M": 1.0,
    "mM": 1e-3,
    "uM": 1e-6,
    "um": 1e-6,
    "µM": 1e-6,
    "nM": 1e-9,
    "pM": 1e-12,
    "fM": 1e-15,
}


def stable_hash(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def canonical_smiles(value: str) -> str:
    molecule = Chem.MolFromSmiles(value)
    if molecule is None:
        raise ValueError("RDKit rejected canonical parent SMILES")
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def p_affinity(value: float, units: str) -> float:
    molar = float(value) * UNIT_TO_MOLAR[units]
    if not math.isfinite(molar) or molar <= 0:
        raise ValueError("affinity must be finite and positive")
    return -math.log10(molar)


def assay_parameters(connection: sqlite3.Connection) -> dict[int, list[dict[str, Any]]]:
    columns = {row[1] for row in connection.execute("PRAGMA table_info(assay_parameters)")}
    required = {
        "assay_id", "standard_type", "standard_relation", "standard_value",
        "standard_units", "standard_text_value",
    }
    if not required.issubset(columns):
        raise RuntimeError(f"assay_parameters schema missing: {sorted(required - columns)}")
    grouped: dict[int, list[dict[str, Any]]] = defaultdict(list)
    query = """
        SELECT assay_id, standard_type, standard_relation, standard_value,
               standard_units, standard_text_value
        FROM assay_parameters
        ORDER BY assay_id, standard_type, standard_value, standard_units,
                 standard_text_value
    """
    for row in connection.execute(query):
        grouped[int(row[0])].append({
            "type": row[1], "relation": row[2], "value": row[3],
            "units": row[4], "text_value": row[5],
        })
    return grouped


def context_payload(row: dict[str, Any], parameters: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "assay_organism": row["assay_organism"],
        "bao_format": row["bao_format"],
        "cell_id": row["cell_id"],
        "tissue_id": row["tissue_id"],
        "subcellular_fraction": row["subcellular_fraction"],
        "relationship_type": row["relationship_type"],
        "variant_id": row["variant_id"],
        "variant_mutation": row["target_variant_mutation"],
        "component_accession": row["target_accession"],
        "parameters": parameters,
    }


def canonicalize_row(row: dict[str, Any], parameters: list[dict[str, Any]]) -> dict[str, Any]:
    sequence = "".join(str(row["protein_sequence"]).split()).upper()
    smiles = canonical_smiles(str(row["canonical_smiles"]))
    inchi_key = str(row["standard_inchi_key"]).strip().upper()
    if len(inchi_key.split("-")[0]) != 14:
        raise ValueError("invalid standard InChIKey connectivity block")
    context = context_payload(row, parameters)
    context_json = json.dumps(context, sort_keys=True, separators=(",", ":"), default=str)
    endpoint = str(row["endpoint_family"])
    assay = str(row["assay_chembl_id"])
    sequence_sha = stable_hash(sequence)
    context_sha = stable_hash(context_json)
    task_key = [sequence_sha, endpoint, assay, context_sha]
    return {
        "schema": "MetaSieve.AffinityEnergyRow.v1",
        "chembl_release": "37",
        "activity_id": int(row["activity_id"]),
        "assay_chembl_id": assay,
        "document_chembl_id": row["document_chembl_id"] or None,
        "document_doi": row["document_doi"] or None,
        "document_pmid": row["document_pmid"],
        "document_year": row["document_year"],
        "molecule_chembl_id": str(row["molecule_chembl_id"]),
        "parent_molregno": int(row["parent_molregno"]),
        "canonical_smiles": smiles,
        "standard_inchi_key": inchi_key,
        "ligand_connectivity_key": inchi_key.split("-")[0],
        "target_chembl_id": str(row["target_chembl_id"]),
        "target_component_id": int(row["target_component_id"]),
        "target_accession": row["target_accession"],
        "protein_sequence": sequence,
        "protein_sequence_sha256": sequence_sha,
        "target_variant_id": row["target_variant_id"],
        "target_variant_mutation": row["target_variant_mutation"],
        "endpoint_family": endpoint,
        "standard_relation": "=",
        "standard_value": float(row["standard_value"]),
        "standard_units": str(row["standard_units"]),
        "p_affinity": p_affinity(float(row["standard_value"]), str(row["standard_units"])),
        "published_type": row["published_type"],
        "published_relation": row["published_relation"],
        "published_value": row["published_value"],
        "published_units": row["published_units"],
        "pchembl_value_reported": row["pchembl_value_reported"],
        "bao_endpoint": row["bao_endpoint"],
        "activity_comment": row["activity_comment"],
        "source_id": row["src_id"],
        "assay_type": row["assay_type"],
        "assay_confidence": row["assay_confidence"],
        "assay_description": row["assay_description"],
        "assay_context": context,
        "assay_context_sha256": context_sha,
        "task_id": stable_hash(json.dumps(task_key, separators=(",", ":"))),
        "task_keys": task_key,
    }


def _task_state(rows: Iterable[dict[str, Any]]) -> tuple[dict, dict]:
    tasks: dict[str, dict[str, Any]] = {}
    values: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        task_id = row["task_id"]
        if task_id not in tasks:
            tasks[task_id] = {
                "task_id": task_id,
                "task_keys": row["task_keys"],
                "protein_sequence_sha256": row["protein_sequence_sha256"],
                "endpoint_family": row["endpoint_family"],
                "assay_chembl_id": row["assay_chembl_id"],
                "assay_context_sha256": row["assay_context_sha256"],
                "document_ids": set(),
                "measurement_count": 0,
            }
        tasks[task_id]["measurement_count"] += 1
        if row["document_chembl_id"]:
            tasks[task_id]["document_ids"].add(row["document_chembl_id"])
        values[task_id][row["ligand_connectivity_key"]].append(row["p_affinity"])
    return tasks, values


def _manifest_from_state(tasks: dict, values: dict) -> list[dict[str, Any]]:
    output = []
    for task_id in sorted(tasks):
        task = tasks[task_id]
        compound_values = sorted(median(replicates) for replicates in values[task_id].values())
        tied = sum(
            compound_values[index] == compound_values[other]
            for index in range(len(compound_values))
            for other in range(index + 1, len(compound_values))
        )
        total = len(compound_values) * (len(compound_values) - 1) // 2
        output.append({
            **{key: value for key, value in task.items() if key != "document_ids"},
            "document_ids": sorted(task["document_ids"]),
            "exact_compound_count": len(compound_values),
            "pair_comparisons": total,
            "non_tied_pair_comparisons": total - tied,
            "eligible_e0_core": len(compound_values) >= 20 and total - tied >= 12,
        })
    return output


def task_manifest(rows: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    return _manifest_from_state(*_task_state(rows))


def build_corpus(database: Path, sql_path: Path, output: Path) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=True)
    sql = sql_path.read_text(encoding="utf-8")
    connection = sqlite3.connect(f"file:{database.as_posix()}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        parameters = assay_parameters(connection)
        rows_path = output / "canonical_rows.jsonl"
        rejected_path = output / "rejected_rows.jsonl"
        rows_partial = rows_path.with_suffix(".jsonl.partial")
        rejected_partial = rejected_path.with_suffix(".jsonl.partial")
        task_state: dict[str, dict[str, Any]] = {}
        value_state: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
        rejection_counts: Counter[str] = Counter()
        accepted_rows = rejected_rows = 0
        with rows_partial.open("wb") as accepted_handle, rejected_partial.open("wb") as rejected_handle:
            for processed_rows, raw in enumerate(connection.execute(sql), start=1):
                try:
                    row = canonicalize_row(
                        dict(raw), parameters.get(int(raw["assay_id"]), []),
                    )
                except ValueError as error:
                    reason = str(error)
                    rejected_rows += 1
                    rejection_counts[reason] += 1
                    rejected_handle.write((json.dumps({
                        "activity_id": int(raw["activity_id"]), "reason": reason,
                    }, sort_keys=True, separators=(",", ":")) + "\n").encode())
                    continue
                accepted_rows += 1
                accepted_handle.write(
                    (json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n").encode()
                )
                task_id = row["task_id"]
                if task_id not in task_state:
                    task_state[task_id] = {
                        "task_id": task_id,
                        "task_keys": row["task_keys"],
                        "protein_sequence_sha256": row["protein_sequence_sha256"],
                        "endpoint_family": row["endpoint_family"],
                        "assay_chembl_id": row["assay_chembl_id"],
                        "assay_context_sha256": row["assay_context_sha256"],
                        "document_ids": set(),
                        "measurement_count": 0,
                    }
                task_state[task_id]["measurement_count"] += 1
                if row["document_chembl_id"]:
                    task_state[task_id]["document_ids"].add(row["document_chembl_id"])
                value_state[task_id][row["ligand_connectivity_key"]].append(row["p_affinity"])
                if processed_rows % 50000 == 0:
                    print(
                        f"processed_rows={processed_rows} accepted_rows={accepted_rows} "
                        f"rejected_rows={rejected_rows}",
                        flush=True,
                    )
        rows_partial.replace(rows_path)
        rejected_partial.replace(rejected_path)
    finally:
        connection.close()
    tasks = _manifest_from_state(task_state, value_state)
    tasks_path = output / "task_manifest.json"
    write_canonical_json(tasks_path, tasks)
    manifest = {
        "schema": "MetaSieve.AffinityEnergyCorpus.v1",
        "name": "EnergyPilot.v1",
        "source": "ChEMBL37 static SQLite",
        "scope": "E0-Core Ki/Kd",
        "sql_sha256": sha256_file(sql_path),
        "normalizer_sha256": sha256_file(Path(__file__)),
        "row_schema_sha256": sha256_file(
            Path(__file__).resolve().parents[2]
            / "contracts/source_affinity/affinity_energy_row_v1.json"
        ),
        "rdkit_version": rdBase.rdkitVersion,
        "canonical_rows": {"rows": accepted_rows, "sha256": sha256_file(rows_path)},
        "rejected_rows": {
            "rows": rejected_rows,
            "reasons": dict(sorted(rejection_counts.items())),
            "sha256": sha256_file(rejected_path),
        },
        "tasks": {
            "count": len(tasks),
            "eligible_e0_core": sum(task["eligible_e0_core"] for task in tasks),
            "sha256": sha256_file(tasks_path),
        },
        "recipient_labels_read": False,
        "training_authorized": False,
    }
    write_canonical_json(output / "corpus_manifest.json", manifest)
    return manifest
