"""R2-E3: schema-explicit, label-free census of local crossed panels.

Affinity arrays are never loaded. Training interaction degrees of freedom and
dependency-closed confirmation units are reported separately for each panel;
they are never summed to manufacture a confirmation PASS.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from research.meta_fewshot.train_main_v0 import sha256

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT / "report/meta_fewshot/r2_crossed_census.json"
OVERLAP_REPORT = ROOT / "report/meta_fewshot/r2_crossed_overlap_audit.json"
SCHEMA = "MetaSieve.R2CrossedCensus.v2"
PATHS = {
    "BLK-METZ-XP2": ROOT / "dataset/processed/crossed_panels_xp2/blk_metz_xp2.npz",
    "BLK-BDB-PANELS": ROOT / "dataset/processed/multipanel/blk_bdb_panels.npz",
    "PDSP-CORE": ROOT / "dataset/processed/crossed_panels/pdsp_core.npz",
}


class DisjointSet:
    def __init__(self):
        self.parent = {}

    def find(self, value):
        self.parent.setdefault(value, value)
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[right] = left


def codes(values: np.ndarray) -> np.ndarray:
    mapping = {value: index for index, value in enumerate(sorted(set(values.tolist())))}
    return np.asarray([mapping[value] for value in values.tolist()], dtype=np.int64)


def connected_blocks(protein: np.ndarray, ligand: np.ndarray) -> list[np.ndarray]:
    dsu = DisjointSet()
    for p_value, l_value in zip(protein, ligand):
        dsu.union(("p", int(p_value)), ("l", int(l_value)))
    groups = defaultdict(list)
    for row, p_value in enumerate(protein):
        groups[dsu.find(("p", int(p_value)))].append(row)
    return [np.asarray(rows, dtype=np.int64) for rows in groups.values()]


def interaction_df(protein_raw: np.ndarray, ligand_raw: np.ndarray,
                   block_raw: np.ndarray) -> tuple[int, int]:
    total, connected = 0, 0
    for block in sorted(set(block_raw.tolist())):
        rows = np.flatnonzero(block_raw == block)
        protein, ligand = codes(protein_raw[rows]), codes(ligand_raw[rows])
        for component_rows in connected_blocks(protein, ligand):
            connected += 1
            n_protein = len(set(protein[component_rows].tolist()))
            n_ligand = len(set(ligand[component_rows].tolist()))
            total += max(0, len(component_rows) - (n_protein + n_ligand - 1))
    return int(total), connected


def dependency_components(block: np.ndarray, ligand: np.ndarray,
                          scaffold: np.ndarray | None,
                          cluster: np.ndarray | None) -> tuple[int, float]:
    dsu = DisjointSet()
    for row in range(len(block)):
        anchor = ("row", row)
        dsu.union(anchor, ("block", str(block[row])))
        dsu.union(anchor, ("ligand", str(ligand[row])))
        if scaffold is not None:
            dsu.union(anchor, ("scaffold", str(scaffold[row])))
        if cluster is not None:
            dsu.union(anchor, ("cluster", str(cluster[row])))
    sizes = defaultdict(int)
    for row in range(len(block)):
        sizes[dsu.find(("row", row))] += 1
    ordered = sorted(sizes.values(), reverse=True)
    return len(ordered), float(ordered[0] / len(block)) if ordered else 0.0


def dense_matrix_design(path: Path, name: str) -> dict:
    with np.load(path, allow_pickle=True) as stored:
        mask = stored["M"].astype(bool)
        ligand_index, protein_index = np.where(mask)
        if name == "BLK-METZ-XP2":
            ligand_values = stored["cmpd"][ligand_index]
            protein_values = stored["kinase"][protein_index]
            scaffold = stored["scaffold_component"][ligand_index]
            cluster = stored["group"][protein_index]
            freshness = "HISTORICALLY_CONSUMED_XP1_XP2_DEVELOPMENT"
        elif name == "PDSP-CORE":
            ligand_values = stored["ligands"][ligand_index]
            protein_values = stored["targets"][protein_index]
            scaffold = None
            cluster = stored["cluster"][protein_index]
            freshness = "HISTORICALLY_CONSUMED_XP1_DEVELOPMENT"
        else:
            raise ValueError(name)
    block = np.zeros(len(ligand_values), dtype=np.int8)
    return design_result(name, path, protein_values, ligand_values, block,
                         scaffold, cluster, freshness)


def bindingdb_design(path: Path) -> dict:
    with np.load(path, allow_pickle=True) as stored:
        # Deliberately do not access `pki`.
        protein = stored["uni"].copy()
        ligand = stored["smiles"].copy()
        block = stored["pmid"].copy()
        scaffold = stored["scaffold_component"].copy()
        cluster = stored["protein_cluster"].copy()
    return design_result(
        "BLK-BDB-PANELS", path, protein, ligand, block, scaffold, cluster,
        "KNOWN_EXACT_PROTEIN_LIGAND_DOCUMENT_OVERLAP_WITH_MAIN_V0")


def design_result(name: str, path: Path, protein: np.ndarray, ligand: np.ndarray,
                  block: np.ndarray, scaffold: np.ndarray | None,
                  cluster: np.ndarray | None, freshness: str) -> dict:
    value_df, connected = interaction_df(protein, ligand, block)
    components, largest = dependency_components(block, ligand, scaffold, cluster)
    closure_complete = scaffold is not None and cluster is not None
    training = value_df >= 5000
    confirmation_design = closure_complete and components >= 30 and largest <= 0.5
    return {
        "panel": name,
        "status": "OK",
        "outcome_values_read": 0,
        "rows": len(protein),
        "n_proteins": len(set(protein.tolist())),
        "n_ligands": len(set(ligand.tolist())),
        "assay_blocks": len(set(block.tolist())),
        "connected_blocks_within_assay": connected,
        "interaction_df": value_df,
        "dependency_closure_complete": closure_complete,
        "dependency_components": components,
        "largest_component_share": largest,
        "training_supply_eligible": training,
        "confirmation_design_eligible": confirmation_design,
        "freshness": freshness,
        "confirmation_supply_eligible": confirmation_design and freshness.startswith("FRESH_"),
        "artifact_sha256": sha256(path),
    }


def census(list_keys: bool = False) -> dict:
    if list_keys:
        panels = []
        for name, path in PATHS.items():
            with np.load(path, allow_pickle=True) as stored:
                panels.append({"panel": name, "keys": stored.files})
        return {"schema": SCHEMA, "mode": "KEYS_ONLY", "panels": panels}
    panels = [
        dense_matrix_design(PATHS["BLK-METZ-XP2"], "BLK-METZ-XP2"),
        bindingdb_design(PATHS["BLK-BDB-PANELS"]),
        dense_matrix_design(PATHS["PDSP-CORE"], "PDSP-CORE"),
    ]
    result = {
        "schema": SCHEMA,
        "declared_role": "LABEL_FREE_DESIGN_CENSUS_NOT_CONFIRMATORY",
        "outcome_values_read": 0,
        "panels": panels,
        "training_supplies": [row["panel"] for row in panels if row["training_supply_eligible"]],
        "confirmation_design_supplies": [
            row["panel"] for row in panels if row["confirmation_design_eligible"]],
        "fresh_confirmation_supplies": [
            row["panel"] for row in panels if row["confirmation_supply_eligible"]],
        "overlap_audit": {
            "path": str(OVERLAP_REPORT.relative_to(ROOT)),
            "sha256": sha256(OVERLAP_REPORT) if OVERLAP_REPORT.exists() else None,
        },
    }
    if result["fresh_confirmation_supplies"]:
        verdict = "LOCAL_FRESH_CONFIRMATION_SUPPLY_IDENTIFIED"
    elif result["training_supplies"]:
        verdict = "LOCAL_CROSSED_TRAINING_SUPPLY_EXISTS_NO_FRESH_CONFIRMATION"
    else:
        verdict = "LOCAL_CROSSED_TRAINING_SUPPLY_NOT_IDENTIFIED"
    result["verdict"] = verdict
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--list-keys", action="store_true")
    args = parser.parse_args()
    result = census(args.list_keys)
    if not args.list_keys:
        REPORT.parent.mkdir(parents=True, exist_ok=True)
        REPORT.write_text(json.dumps(result, allow_nan=False, indent=2, sort_keys=True) + "\n",
                          encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
