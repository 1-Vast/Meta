"""Materialize the governed BindingDB Ki quotient training corpus."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from research.crossed_interaction.bindingdb_cq_r0 import iter_projection, sha256_file
from research.crossed_interaction.bindingdb_cq_r1 import iter_labels, panel_residual


def stable_hash(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


class UnionFind:
    def __init__(self, values):
        self.parent = {value: value for value in values}

    def find(self, value):
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[max(left, right)] = min(left, right)


def sequence_identity(left: str, right: str) -> float:
    import parasail

    trace = parasail.nw_trace_striped_16(
        left, right, 10, 1, parasail.blosum62
    ).traceback
    matches = sum(
        a == b for a, b in zip(trace.query, trace.ref) if a != "-" and b != "-"
    )
    return matches / min(len(left), len(right))


def cluster_sequences(sequences: dict[str, str], workers: int = 8) -> dict[str, str]:
    keys = sorted(sequences)
    union = UnionFind(keys)

    def compare_left(index: int):
        left_key = keys[index]
        left = sequences[left_key]
        matches = []
        for right_key in keys[index + 1 :]:
            right = sequences[right_key]
            ratio = len(left) / len(right)
            if not 0.5 <= ratio <= 2.0:
                continue
            if sequence_identity(left, right) >= 0.40:
                matches.append(right_key)
        return left_key, matches

    with ThreadPoolExecutor(max_workers=workers) as executor:
        for left_key, matches in executor.map(compare_left, range(len(keys))):
            for right_key in matches:
                union.union(left_key, right_key)
    groups = defaultdict(list)
    for key in keys:
        groups[union.find(key)].append(key)
    result = {}
    for members in groups.values():
        group_id = stable_hash("protein40|" + "|".join(sorted(members)))
        for key in members:
            result[key] = group_id
    return result


def murcko_scaffold(smiles: str) -> str:
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return ""
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    return Chem.MolToSmiles(scaffold, canonical=True, isomericSmiles=True)


def canonical_smiles(smiles: str) -> str:
    from rdkit import Chem

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return ""
    return Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)


def _gzip_writer(path: Path):
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return raw, io.TextIOWrapper(compressed, encoding="utf-8", newline="\n")


def _write_jsonl_gz(path: Path, rows: list[dict]) -> None:
    raw, writer = _gzip_writer(path)
    try:
        for row in rows:
            writer.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    finally:
        writer.close()
        raw.close()


def build_corpus(projection: Path, labels: Path, output: Path, workers: int = 8) -> dict:
    metadata = {row["source_row_id"]: row for row in iter_projection(projection)}
    ki_rows = [row for row in iter_labels(labels) if row["endpoint"] == "Ki"]
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in ki_rows:
        grouped[row["panel_id"]].append(row)
    eligible_panels = {
        panel_id for panel_id, rows in grouped.items() if panel_residual(rows) is not None
    }

    cell_values: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in ki_rows:
        if row["panel_id"] in eligible_panels:
            cell_values[(row["panel_id"], row["target_id"], row["ligand_id"])].append(row)

    candidates = []
    smiles_by_ligand: dict[str, set[str]] = defaultdict(set)
    sequences_by_target: dict[str, set[str]] = defaultdict(set)
    conflicting_cells = 0
    for (panel_id, target_id, ligand_id), values in sorted(cell_values.items()):
        records = [metadata[value["source_row_id"]] for value in values]
        smiles_values = {canonical_smiles(record["ligand_smiles"]) for record in records}
        sequence_values = {record["target_sequence"] for record in records}
        if "" in smiles_values or len(smiles_values) != 1 or len(sequence_values) != 1:
            conflicting_cells += 1
            continue
        smiles = smiles_values.pop()
        sequence = sequence_values.pop()
        smiles_by_ligand[ligand_id].add(smiles)
        sequences_by_target[target_id].add(sequence)
        candidates.append(
            {
                "cell_id": stable_hash(f"{panel_id}|{target_id}|{ligand_id}"),
                "panel_id": panel_id,
                "document_id": values[0]["document_id"],
                "target_id": target_id,
                "ligand_id": ligand_id,
                "pK": float(np.mean([value["pK"] for value in values])),
                "replicates": len(values),
                "source_row_ids": sorted(value["source_row_id"] for value in values),
            }
        )

    conflicting_ligands = {
        key for key, values in smiles_by_ligand.items() if len(values) != 1
    }
    conflicting_targets = {
        key for key, values in sequences_by_target.items() if len(values) != 1
    }
    cells = [
        cell
        for cell in candidates
        if cell["ligand_id"] not in conflicting_ligands
        and cell["target_id"] not in conflicting_targets
    ]
    ligand_smiles = {
        key: next(iter(values))
        for key, values in smiles_by_ligand.items()
        if key not in conflicting_ligands
    }
    sequences = {
        key: next(iter(values))
        for key, values in sequences_by_target.items()
        if key not in conflicting_targets
    }

    scaffolds = {key: murcko_scaffold(smiles) for key, smiles in ligand_smiles.items()}
    invalid_ligands = {key for key, value in scaffolds.items() if not value}
    if invalid_ligands:
        cells = [cell for cell in cells if cell["ligand_id"] not in invalid_ligands]
        ligand_smiles = {
            key: value for key, value in ligand_smiles.items() if key not in invalid_ligands
        }
        scaffolds = {key: value for key, value in scaffolds.items() if value}

    protein_groups = cluster_sequences(sequences, workers=workers)
    panel_cells: dict[str, list[dict]] = defaultdict(list)
    for cell in cells:
        panel_cells[cell["panel_id"]].append(cell)
    # Scaffold filtering can break cycles, so recheck every panel.
    retained_panels = {}
    for panel_id, values in panel_cells.items():
        if panel_residual(values) is not None:
            retained_panels[panel_id] = values
    cells = [cell for panel_id in sorted(retained_panels) for cell in retained_panels[panel_id]]

    panel_ids = sorted(retained_panels)
    dependencies = UnionFind(panel_ids)
    first_seen = {}
    for panel_id, values in retained_panels.items():
        document = values[0]["document_id"]
        attributes = [("document", document)]
        attributes += [
            ("protein40", protein_groups[cell["target_id"]]) for cell in values
        ]
        attributes += [("scaffold", scaffolds[cell["ligand_id"]]) for cell in values]
        for attribute in set(attributes):
            if attribute in first_seen:
                dependencies.union(panel_id, first_seen[attribute])
            else:
                first_seen[attribute] = panel_id

    components: dict[str, list[str]] = defaultdict(list)
    for panel_id in panel_ids:
        components[dependencies.find(panel_id)].append(panel_id)
    component_items = []
    for members in components.values():
        component_id = stable_hash("component|" + "|".join(sorted(members)))
        edge_count = sum(len(retained_panels[panel]) for panel in members)
        component_items.append((component_id, sorted(members), edge_count))
    component_items.sort()
    split_by_panel = {}
    split_components = defaultdict(list)
    for component_id, members, edge_count in component_items:
        split = "development" if int(component_id[:8], 16) % 5 == 0 else "train"
        split_components[split].append(component_id)
        for panel_id in members:
            split_by_panel[panel_id] = split

    panels = []
    for panel_id, values in sorted(retained_panels.items()):
        result = panel_residual(values)
        assert result is not None
        component_id = next(
            cid for cid, members, _ in component_items if panel_id in members
        )
        panels.append(
            {
                "panel_id": panel_id,
                "document_id": values[0]["document_id"],
                "split": split_by_panel[panel_id],
                "dependency_component": component_id,
                "cell_ids": sorted(cell["cell_id"] for cell in values),
                "retained_rank": result["retained_rank"],
                "edges": result["edges"],
                "targets": result["targets"],
                "ligands": result["ligands"],
            }
        )
    for cell in cells:
        cell["split"] = split_by_panel[cell["panel_id"]]
        cell["protein_group_40"] = protein_groups[cell["target_id"]]
        cell["scaffold"] = scaffolds[cell["ligand_id"]]

    output.mkdir(parents=True, exist_ok=True)
    _write_jsonl_gz(output / "cells.jsonl.gz", sorted(cells, key=lambda row: row["cell_id"]))
    _write_jsonl_gz(output / "panels.jsonl.gz", panels)
    protein_rows = [
        {"sequence_sha256": key, "sequence": value}
        for key, value in sorted(sequences.items())
        if any(cell["target_id"] == key for cell in cells)
    ]
    ligand_rows = [
        {"drug_key": key, "smiles": value, "scaffold": scaffolds[key]}
        for key, value in sorted(ligand_smiles.items())
        if any(cell["ligand_id"] == key for cell in cells)
    ]
    (output / "proteins.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in protein_rows),
        encoding="utf-8",
    )
    (output / "ligands.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in ligand_rows),
        encoding="utf-8",
    )
    total_edges = sum(item[2] for item in component_items)
    split_summary = {}
    for split in ("train", "development"):
        split_panels = [panel for panel in panels if panel["split"] == split]
        split_summary[split] = {
            "components": len(split_components[split]),
            "panels": len(split_panels),
            "cells": sum(panel["edges"] for panel in split_panels),
            "retained_rank": sum(panel["retained_rank"] for panel in split_panels),
        }
    largest_share = max((item[2] for item in component_items), default=0) / max(total_edges, 1)
    ready = (
        split_summary["train"]["retained_rank"] > 0
        and split_summary["development"]["retained_rank"] > 0
        and split_summary["development"]["components"] >= 5
    )
    manifest = {
        "schema": "MetaSieve.BindingDB.CQTrainingCorpus.v1",
        "projection_sha256": sha256_file(projection),
        "labels_sha256": sha256_file(labels),
        "cells": len(cells),
        "panels": len(panels),
        "proteins": len(protein_rows),
        "ligands": len(ligand_rows),
        "protein_groups_40": len(set(protein_groups.values())),
        "dependency_components": len(component_items),
        "largest_component_share": largest_share,
        "identity_conflicting_cells_excluded": conflicting_cells,
        "identity_conflicting_ligands_excluded": len(conflicting_ligands),
        "identity_conflicting_targets_excluded": len(conflicting_targets),
        "invalid_or_scaffoldless_ligands_excluded": len(invalid_ligands),
        "splits": split_summary,
        "development_training_ready": ready,
        "biological_claim_ready": False,
        "files": {},
    }
    for filename in ("cells.jsonl.gz", "panels.jsonl.gz", "proteins.jsonl", "ligands.jsonl"):
        manifest["files"][filename] = sha256_file(output / filename)
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--projection", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    result = build_corpus(args.projection, args.labels, args.output, args.workers)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
