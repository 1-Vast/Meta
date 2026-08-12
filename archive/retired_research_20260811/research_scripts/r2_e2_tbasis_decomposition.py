"""R2-E2: label-free, design-aware decomposition of frozen T-BASIS features.

The result describes the observed sparse design; it does not prove architectural
capacity or absence of nonlinear biological information. In particular, the
additive residual is evaluated both on all rows and on the bipartite 2-core so
ligand singletons cannot manufacture an artificially tiny interaction residual.
"""
from __future__ import annotations

import argparse
import gzip
import json
from collections import defaultdict, deque
from pathlib import Path

import numpy as np

from research.meta_fewshot.train_main_v0 import sha256

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
STRUCTURAL_INDEX = CORPUS / "r2_structural_index.jsonl.gz"
FEATURES = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_tbasis_features.npz"
REPORT = ROOT / "report/meta_fewshot/r2_tbasis_decomposition.json"
SCHEMA = "MetaSieve.R2TbasisDecomposition.v2"


def read_cells() -> list[dict]:
    if not STRUCTURAL_INDEX.exists():
        raise FileNotFoundError(
            "run research.meta_fewshot.r2_build_structural_index first")
    with gzip.open(STRUCTURAL_INDEX, "rt", encoding="utf-8") as handle:
        rows = [json.loads(line) for line in handle if line.strip()]
    allowed = {
        "cell_id", "target_id", "ligand_id", "protein_group_40", "split", "panel_ids"
    }
    if any(set(row) != allowed for row in rows):
        raise ValueError("structural index is not exactly label-redacted")
    return rows


def reindex(values: np.ndarray) -> np.ndarray:
    mapping = {value: index for index, value in enumerate(sorted(set(values.tolist())))}
    return np.asarray([mapping[value] for value in values.tolist()], dtype=np.int64)


def group_mean(values: np.ndarray, index: np.ndarray, size: int) -> np.ndarray:
    result = np.zeros((size, values.shape[1]), dtype=np.float64)
    np.add.at(result, index, values)
    counts = np.bincount(index, minlength=size)[:, None]
    return np.divide(result, counts, out=np.zeros_like(result), where=counts > 0)


def additive_projection(phi: np.ndarray, protein: np.ndarray, ligand: np.ndarray,
                        max_iter: int = 100, tolerance: float = 1e-10) -> dict:
    """Least-squares additive projection by converged backfitting."""
    protein, ligand = reindex(protein), reindex(ligand)
    centered = phi - phi.mean(axis=0, keepdims=True)
    n_protein, n_ligand = protein.max() + 1, ligand.max() + 1
    alpha = np.zeros((n_protein, phi.shape[1]), dtype=np.float64)
    beta = np.zeros((n_ligand, phi.shape[1]), dtype=np.float64)
    fitted = np.zeros_like(centered)
    iterations = 0
    relative = float("inf")
    for iterations in range(1, max_iter + 1):
        alpha = group_mean(centered - beta[ligand], protein, n_protein)
        beta = group_mean(centered - alpha[protein], ligand, n_ligand)
        updated = alpha[protein] + beta[ligand]
        relative = np.linalg.norm(updated - fitted) / (np.linalg.norm(updated) + 1e-12)
        fitted = updated
        if relative < tolerance:
            break
    residual = centered - fitted
    total = float(np.square(centered).sum())
    protein_only = centered - group_mean(centered, protein, n_protein)[protein]
    ligand_only = centered - group_mean(centered, ligand, n_ligand)[ligand]
    residual_sse = float(np.square(residual).sum())
    return {
        "additive_explained_fraction": 1.0 - residual_sse / total,
        "interaction_residual_fraction": residual_sse / total,
        "protein_increment_beyond_ligand_fraction": (
            float(np.square(ligand_only).sum()) - residual_sse) / total,
        "ligand_increment_beyond_protein_fraction": (
            float(np.square(protein_only).sum()) - residual_sse) / total,
        "iterations": iterations,
        "converged": bool(relative < tolerance),
        "final_relative_change": float(relative),
        "tolerance": tolerance,
        "_alpha": alpha,
    }


def bipartite_two_core(protein: np.ndarray, ligand: np.ndarray) -> np.ndarray:
    """Rows remaining after recursively removing degree<2 protein/ligand nodes."""
    p_rows, l_rows = defaultdict(set), defaultdict(set)
    active = np.ones(len(protein), dtype=bool)
    for row, (p_value, l_value) in enumerate(zip(protein, ligand)):
        p_rows[int(p_value)].add(row)
        l_rows[int(l_value)].add(row)
    queue = deque(
        [("p", key) for key, rows in p_rows.items() if len(rows) < 2]
        + [("l", key) for key, rows in l_rows.items() if len(rows) < 2]
    )
    while queue:
        kind, key = queue.popleft()
        rows = p_rows[key] if kind == "p" else l_rows[key]
        for row in list(rows):
            if not active[row]:
                continue
            active[row] = False
            p_value, l_value = int(protein[row]), int(ligand[row])
            p_rows[p_value].discard(row)
            l_rows[l_value].discard(row)
            if len(p_rows[p_value]) == 1:
                queue.append(("p", p_value))
            if len(l_rows[l_value]) == 1:
                queue.append(("l", l_value))
    return active


def graph_component_membership(protein: np.ndarray, ligand: np.ndarray) -> dict[int, object]:
    parent = {}

    def find(value):
        parent.setdefault(value, value)
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left, right):
        left, right = find(left), find(right)
        if left != right:
            parent[right] = left

    for p_value, l_value in zip(protein, ligand):
        union(("p", int(p_value)), ("l", int(l_value)))
    return {int(value): find(("p", int(value))) for value in set(protein.tolist())}


def run() -> dict:
    cells = read_cells()
    with np.load(FEATURES, allow_pickle=False) as stored:
        if stored["cell_id"].tolist() != [row["cell_id"] for row in cells]:
            raise ValueError("feature rows do not match corpus cell order")
        correct = stored["correct"].astype(np.float64)
        wrong = stored["deranged_protein"].astype(np.float64)
    scale = correct.std(axis=0)
    scale[scale < 1e-9] = 1.0
    mean = correct.mean(axis=0)
    phi, phi_wrong = (correct - mean) / scale, (wrong - mean) / scale

    protein_keys = sorted({row["target_id"] for row in cells})
    ligand_keys = sorted({row["ligand_id"] for row in cells})
    proteins = {key: index for index, key in enumerate(protein_keys)}
    ligands = {key: index for index, key in enumerate(ligand_keys)}
    p_index = np.asarray([proteins[row["target_id"]] for row in cells])
    l_index = np.asarray([ligands[row["ligand_id"]] for row in cells])

    all_projection = additive_projection(phi, p_index, l_index)
    all_projection.pop("_alpha")
    core = bipartite_two_core(p_index, l_index)
    core_projection = additive_projection(phi[core], p_index[core], l_index[core])
    core_alpha = core_projection.pop("_alpha")
    core_p, core_l = p_index[core], l_index[core]
    component_of = graph_component_membership(core_p, core_l)
    core_components = len(set(component_of.values()))
    interaction_df = int(core.sum() - (
        len(set(core_p.tolist())) + len(set(core_l.tolist())) - core_components))

    by_ligand = defaultdict(lambda: defaultdict(list))
    for row, (ligand, protein) in enumerate(zip(l_index, p_index)):
        by_ligand[int(ligand)][int(protein)].append(row)
    within, shared = [], 0
    for targets in by_ligand.values():
        if len(targets) < 2:
            continue
        shared += 1
        block = np.stack([
            phi[rows].mean(axis=0) for rows in targets.values()
        ])
        within.append(float(np.square(block - block.mean(axis=0)).sum(axis=1).mean()))
    global_dispersion = float(np.square(phi - phi.mean(axis=0)).sum(axis=1).mean())
    partner_fraction = float(np.mean(within)) / global_dispersion if within else None

    cluster_of = {}
    for row in cells:
        cluster_of[proteins[row["target_id"]]] = row["protein_group_40"]
    core_protein_keys = sorted(set(core_p.tolist()))
    core_alpha_by_protein = {
        protein: core_alpha[position]
        for position, protein in enumerate(core_protein_keys)
    }
    within_pairs, between_pairs = [], []
    keys = core_protein_keys
    for position, left in enumerate(keys):
        for right in keys[position + 1:]:
            if component_of[left] != component_of[right]:
                continue
            distance = float(np.linalg.norm(
                core_alpha_by_protein[left] - core_alpha_by_protein[right]))
            (within_pairs if cluster_of[left] == cluster_of[right]
             else between_pairs).append(distance)
    homolog_ratio = (float(np.mean(within_pairs)) / float(np.mean(between_pairs))
                     if within_pairs and between_pairs else None)

    corruption = float(np.linalg.norm(phi - phi_wrong, axis=1).mean())
    natural = float(np.sqrt(np.mean(within))) if within else None
    ligand_counts = np.bincount(l_index, minlength=len(ligands))
    result = {
        "schema": SCHEMA,
        "declared_role": "LABEL_FREE_OBSERVED_DESIGN_DESCRIPTION_NOT_CAPACITY_PROOF",
        "affinity_values_read": 0,
        "structural_index_physically_label_redacted": True,
        "n_cells": len(cells),
        "n_proteins": len(proteins),
        "n_ligands": len(ligands),
        "ligand_singleton_fraction": float((ligand_counts == 1).mean()),
        "basis_shape": "8 atom channels x 6 residue classes x 6 radial shells = 288",
        "protein_path": (
            "ESM residue states condition bridge distance logits; the final radial "
            "aggregation also contracts against six residue-chemistry classes"),
        "all_rows_additive_projection": all_projection,
        "crossed_two_core": {
            "rows": int(core.sum()),
            "proteins": len(set(core_p.tolist())),
            "ligands": len(set(core_l.tolist())),
            "connected_components": core_components,
            "interaction_df": interaction_df,
            **core_projection,
        },
        "fixed_ligand_partner_dispersion": {
            "ligands_with_at_least_two_unique_proteins": shared,
            "partner_dispersion_fraction": partner_fraction,
        },
        "homolog_resolution": {
            "within_cluster_over_between_cluster_alpha_distance": homolog_ratio,
            "within_pairs": len(within_pairs),
            "between_pairs": len(between_pairs),
            "comparison_scope": "within_bipartite_component_only",
        },
        "deranged_partner_calibration": {
            "mean_feature_shift": corruption,
            "mean_natural_partner_shift": natural,
            "corruption_over_natural": corruption / natural if natural else None,
        },
        "inputs": {
            "features_sha256": sha256(FEATURES),
            "structural_index_sha256": sha256(STRUCTURAL_INDEX),
        },
    }
    verdict = []
    if core_projection["interaction_residual_fraction"] < 0.10:
        verdict.append("TBASIS_OBSERVED_CROSSED_DESIGN_NEAR_ADDITIVE")
    if partner_fraction is not None and partner_fraction < 0.10:
        verdict.append("TBASIS_FIXED_LIGAND_PARTNER_DISPERSION_LOW")
    if homolog_ratio is not None and homolog_ratio > 0.70:
        verdict.append("TBASIS_PROTEIN_MAIN_WEAKLY_SEPARATES_CDHIT40_CLUSTERS")
    result["verdict"] = verdict or ["NO_REGISTERED_REPRESENTATION_DEFECT_DETECTED"]
    return result


def main() -> int:
    argparse.ArgumentParser().parse_args()
    result = run()
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text(json.dumps(result, allow_nan=False, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
