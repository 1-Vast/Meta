"""Collect accepted assay rows before protocol-specific split assignment."""
from __future__ import annotations

import json
from pathlib import Path

from scripts.data_contract import write_jsonl


def _read_rows(roots) -> list[dict]:
    rows = []
    for root in roots:
        for path in sorted(Path(root).rglob("*.jsonl")):
            with path.open(encoding="utf-8") as handle:
                rows.extend(json.loads(line) for line in handle if line.strip())
    return rows


def _components(sequences: dict[str, str], identity, threshold: float = 0.4) -> dict[str, str]:
    parent = {key: key for key in sequences}

    def find(key):
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    def union(left, right):
        left, right = find(left), find(right)
        if left != right:
            parent[max(left, right)] = min(left, right)

    keys = sorted(sequences)
    for index, left in enumerate(keys):
        for right in keys[index + 1:]:
            if float(identity(sequences[left], sequences[right])[0]) >= threshold:
                union(left, right)
    return {key: find(key) for key in keys}


def collect_governed_sources(roots, output_dir: str | Path, identity) -> dict:
    """Deduplicate eligible exact rows and derive homology components once."""
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"governed output already exists: {output}")
    seen, accepted, skipped, repeated = {}, [], 0, 0
    for row in _read_rows(roots):
        if not row.get("model_eligible", False):
            skipped += 1
            continue
        if row.get("standard_relation") != "=":
            skipped += 1
            continue
        identifier = str(row["source_row_id"])
        if identifier in seen:
            if seen[identifier] != row:
                raise ValueError(f"conflicting repeated source_row_id {identifier}")
            repeated += 1
            continue
        seen[identifier] = row
        accepted.append(row)
    sequences = {str(row["protein_sequence_sha256"]): str(row["protein_sequence"])
                 for row in accepted}
    components = _components(sequences, identity)
    identities, labels = [], []
    for row in sorted(accepted, key=lambda value: str(value["source_row_id"])):
        row_id = str(row["source_row_id"])
        identities.append({
            "row_id": row_id,
            "protein_sequence_sha256": str(row["protein_sequence_sha256"]),
            "ligand_connectivity_key": str(row["ligand_connectivity_key"]),
            "task_keys": row["task_keys"], "canonical_smiles": row["canonical_smiles"],
        })
        labels.append({"row_id": row_id, "p_value": float(row["p_value"])})
    assignments = [{"protein_sequence_sha256": key,
                    "homology_component_id": components[key]}
                   for key in sorted(components)]
    output.mkdir(parents=True)
    write_jsonl(output / "identity.jsonl", identities)
    write_jsonl(output / "labels.jsonl", labels)
    write_jsonl(output / "homology_assignments.jsonl", assignments)
    result = {
        "eligible_exact_rows": len(accepted), "proteins": len(sequences),
        "homology_components": len(set(components.values())),
        "repeated_identical_source_ids": repeated,
        "skipped_ineligible_rows": skipped,
        "rows_rejected_by_split_homology_scaffold": 0,
        "document_closure_applied": False,
    }
    (output / "manifest.json").write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result
