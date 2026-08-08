"""Build immutable source/metaval protocols from governed, label-separated rows."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from scripts.data_contract import read_jsonl, write_jsonl


PROTOCOLS = ("fewshot_core_v2", "novel_family_stress_v1")


def _digest(value: object) -> str:
    return hashlib.sha256(str(value).encode("utf-8")).hexdigest()


def _file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _recipient_targets(paths) -> set[str]:
    targets = set()
    for path in paths or ():
        path = Path(path)
        if path.suffix == ".tsv":
            with path.open(encoding="utf-8") as handle:
                header = handle.readline().rstrip("\n").split("\t")
                try:
                    index = header.index("sequence_sha256")
                except ValueError as error:
                    raise ValueError(f"recipient registry lacks sequence_sha256: {path}") from error
                targets.update(line.rstrip("\n").split("\t")[index] for line in handle if line.strip())
        else:
            targets.update(str(row["sequence_sha256"]) for row in read_jsonl(path))
    return targets


def _identity_key(row: dict) -> tuple:
    return (
        str(row["protein_sequence_sha256"]),
        str(row["ligand_connectivity_key"]),
        json.dumps(row.get("task_keys", {}), sort_keys=True),
    )


def _split(identifier: str) -> str:
    return "source" if int(_digest(identifier)[:16], 16) / 2**64 < 0.8 else "metaval"


def _protocol_assignments(assignments: list[dict], protocol: str) -> list[dict]:
    by_protein = {str(row["protein_sequence_sha256"]): dict(row) for row in assignments}
    result = []
    for protein in sorted(by_protein):
        row = by_protein[protein]
        component = str(row.get("homology_component_id", row.get("closure_component_id", protein)))
        key = protein if protocol == "fewshot_core_v2" else component
        row["split"] = _split(f"{protocol}:{key}")
        row["protocol"] = protocol
        result.append(row)
    return result


def _task_rows(identity_rows: list[dict], assignments: dict[str, dict]):
    grouped: dict[tuple[str, str], list[dict]] = {}
    for row in identity_rows:
        protein = str(row["protein_sequence_sha256"])
        task = json.dumps(row.get("task_keys", {}), sort_keys=True)
        grouped.setdefault((protein, task), []).append(row)
    tasks, emitted = [], []
    for (protein, task_key), rows in sorted(grouped.items()):
        by_ligand: dict[str, list[dict]] = {}
        for row in rows:
            by_ligand.setdefault(str(row["ligand_connectivity_key"]), []).append(row)
        ligands = sorted(by_ligand)
        support = set(ligands[:5])
        task_id = _digest(f"{protein}|{task_key}")
        for ligand in ligands:
            role = "support_pool" if ligand in support else "query_pool"
            for row in sorted(by_ligand[ligand], key=lambda value: str(value["row_id"])):
                value = dict(row)
                value["task_id"] = task_id
                value["split"] = assignments[protein]["split"]
                value["episode_role"] = role
                emitted.append(value)
        tasks.append({
            "task_id": task_id, "protein_sequence_sha256": protein,
            "split": assignments[protein]["split"], "support_compounds": len(support),
            "query_compounds": len(ligands) - len(support),
            "eligible_k5_q12": len(support) >= 5 and len(ligands) - len(support) >= 12,
        })
    return emitted, tasks


def build_protocols(identity_path: str | Path, labels_path: str | Path,
                    assignments_path: str | Path, output_dir: str | Path,
                    recipient_identity_paths=None, protein_annotations_path: str | Path | None = None) -> dict:
    """Emit the registered few-shot and homology-stress protocol artifacts."""
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"protocol output already exists: {output}")
    identities = read_jsonl(identity_path)
    labels = read_jsonl(labels_path)
    assignments = read_jsonl(assignments_path)
    by_id = {}
    for row in identities:
        identifier = str(row["row_id"])
        if identifier in by_id and by_id[identifier] != row:
            raise ValueError(f"conflicting identities for row_id {identifier}")
        by_id[identifier] = row
    label_by_id = {str(row["row_id"]): float(row["p_value"]) for row in labels}
    if set(by_id) != set(label_by_id):
        raise ValueError("identity and label row identifiers must match exactly")

    recipient = _recipient_targets(recipient_identity_paths)
    retained = [row for row in by_id.values()
                if str(row["protein_sequence_sha256"]) not in recipient]
    excluded_rows = len(by_id) - len(retained)
    groups: dict[tuple, list[dict]] = {}
    for row in retained:
        groups.setdefault(_identity_key(row), []).append(row)
    conflict_rows, deduplicated, conflicts = [], [], []
    for key, rows in sorted(groups.items(), key=lambda item: str(item[0])):
        representatives = {}
        for row in sorted(rows, key=lambda value: str(value["row_id"])):
            representatives.setdefault(label_by_id[str(row["row_id"])], row)
        values = sorted(representatives)
        if len(values) == 1:
            deduplicated.append(representatives[values[0]])
        else:
            conflict_rows.extend(representatives[value] for value in values)
            conflicts.append({"identity_key": list(key), "p_values": values,
                              "row_ids": sorted(str(row["row_id"]) for row in rows)})
    retained = deduplicated + conflict_rows
    retained.sort(key=lambda row: str(row["row_id"]))
    retained_proteins = {str(row["protein_sequence_sha256"]) for row in retained}
    filtered_assignments = [row for row in assignments
                            if str(row["protein_sequence_sha256"]) in retained_proteins]
    if len({str(row["protein_sequence_sha256"]) for row in filtered_assignments}) != len(retained_proteins):
        raise ValueError("retained identity rows lack homology assignments")
    annotations = read_jsonl(protein_annotations_path) if protein_annotations_path else []
    output.mkdir(parents=True)
    summaries = {}
    duplicate_rows = len(deduplicated) + len(conflict_rows)
    exact_duplicates_removed = len([row for rows in groups.values() for row in rows]) - duplicate_rows
    for protocol in PROTOCOLS:
        root = output / protocol
        protocol_assignments = _protocol_assignments(filtered_assignments, protocol)
        assignment_by_protein = {str(row["protein_sequence_sha256"]): row for row in protocol_assignments}
        protocol_identity, tasks = _task_rows(retained, assignment_by_protein)
        protocol_labels = [{"row_id": str(row["row_id"]),
                            "p_value": label_by_id[str(row["row_id"])]}
                           for row in protocol_identity]
        row_annotations = [{"row_id": str(row["row_id"]),
                            "murcko_scaffold": row.get("murcko_scaffold"),
                            "pocket_id": row.get("pocket_id")}
                           for row in protocol_identity]
        protein_annotations = []
        for row in annotations:
            protein = str(row["protein_sequence_sha256"])
            if protein in assignment_by_protein:
                value = dict(row)
                value["split"] = assignment_by_protein[protein]["split"]
                value["protocol"] = protocol
                protein_annotations.append(value)
        write_jsonl(root / "identity.jsonl", protocol_identity)
        write_jsonl(root / "labels.jsonl", protocol_labels)
        write_jsonl(root / "assignments.jsonl", protocol_assignments)
        write_jsonl(root / "tasks.jsonl", tasks)
        write_jsonl(root / "annotations.jsonl", row_annotations)
        write_jsonl(root / "protein_annotations.jsonl", protein_annotations)
        write_jsonl(root / "replicate_conflicts.jsonl", conflicts)
        artifacts = {name: _file_digest(root / name) for name in (
            "identity.jsonl", "labels.jsonl", "assignments.jsonl", "tasks.jsonl",
            "annotations.jsonl", "protein_annotations.jsonl", "replicate_conflicts.jsonl",
        )}
        summaries[protocol] = {
            "protocol": protocol, "rows": len(protocol_identity), "tasks": len(tasks),
            "recipient_rows_excluded": excluded_rows,
            "recipient_registry_exact_targets": len(recipient),
            "recipient_corpus_exact_targets_excluded": len({
                str(row["protein_sequence_sha256"]) for row in by_id.values()
                if str(row["protein_sequence_sha256"]) in recipient}),
            "exact_duplicate_rows_removed": exact_duplicates_removed,
            "conflicting_replicate_groups_retained": len(conflicts),
            "recipient_label_reads": 0, "homology_components_deleted": 0,
            "scaffolds_deleted": 0, "artifacts": artifacts,
        }
        (root / "manifest.json").write_text(json.dumps(summaries[protocol], indent=2,
                                                   sort_keys=True), encoding="utf-8")
    return summaries
