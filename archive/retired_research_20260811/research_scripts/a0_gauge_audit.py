"""Decompose correct/wrong Meta-Section predictions into gauge-sensitive terms."""
from __future__ import annotations

import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
from scipy import stats

from research.meta_fewshot.train_main_v0 import MetaSectionRegressor, sha256
from research.meta_fewshot.train_v1_development import DEV, OUT as V1_OUT, SEEDS, load_sealed, read_jsonl_gz

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "report/meta_fewshot/a0_gauge_audit.json"


def ridge_kernel(support: np.ndarray, query: np.ndarray, ridge: float) -> np.ndarray:
    gram = support @ support.T
    cross = query @ support.T
    return np.linalg.solve(gram + ridge * np.eye(len(support)), cross.T).T


def procrustes_from_support(wrong: np.ndarray, correct: np.ndarray) -> np.ndarray:
    left, _, right = np.linalg.svd(wrong.T @ correct, full_matrices=False)
    return left @ right


def relative_error(actual: np.ndarray, reference: np.ndarray, eps: float = 1e-12) -> float:
    return float(np.linalg.norm(actual - reference) / (np.linalg.norm(reference) + eps))


def controls() -> dict:
    rng = np.random.default_rng(20260811)
    support = rng.normal(size=(5, 3))
    query = rng.normal(size=(7, 3))
    residual = rng.normal(size=5)
    rotation, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    base = ridge_kernel(support, query, 1.0)
    rotated = ridge_kernel(support @ rotation, query @ rotation, 1.0)
    scale_small = ridge_kernel(0.1 * support, 0.1 * query, 1.0)
    scale_large = ridge_kernel(10.0 * support, 10.0 * query, 1.0)
    shear = np.eye(3)
    shear[0, 1] = 1.0
    sheared = ridge_kernel(support @ shear, query @ shear, 1.0)
    rank_one = np.outer(rng.normal(size=5), rng.normal(size=3))
    rank_one_query = np.outer(rng.normal(size=7), rank_one[0])
    return {
        "orthogonal_H_max_abs": float(np.abs(base - rotated).max()),
        "orthogonal_correction_max_abs": float(np.abs(base @ residual - rotated @ residual).max()),
        "scale_0_1_H_relative": relative_error(scale_small, base),
        "scale_10_H_relative": relative_error(scale_large, base),
        "shear_H_relative": relative_error(sheared, base),
        "rank_one_support_rank": int(np.linalg.matrix_rank(rank_one)),
        "rank_one_kernel_finite": bool(np.isfinite(ridge_kernel(rank_one, rank_one_query, 1.0)).all()),
    }


def prediction_map(path: Path):
    allowed = {"v0_correct", "v0_sc_qw", "v0_sw_qc", "v0_sw_qw", "v0_permuted"}
    result = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            row = json.loads(line)
            if row["arm"] in allowed:
                key = (row["arm"], row["seed"], row["target_id"], row["draw"], row["cell_id"])
                if key in result:
                    raise ValueError(f"duplicate prediction key: {key}")
                result[key] = float(row["prediction_standardized"])
    return result


def summarize(values):
    array = np.asarray(values, dtype=np.float64)
    return {"median": float(np.median(array)), "q25": float(np.quantile(array, 0.25)),
            "q75": float(np.quantile(array, 0.75)), "mean": float(array.mean())}


def audit(device="cpu"):
    _, _, tensors, _, _, _, val_index, y_mean, y_scale = load_sealed(device)
    episodes = read_jsonl_gz(DEV / "metaval_episodes.jsonl.gz")
    predictions = prediction_map(V1_OUT / "predictions_before_query_labels.jsonl.gz")
    episode_rows, max_recompute_error = [], 0.0
    global_rows = []
    for seed in SEEDS:
        checkpoint_path = V1_OUT / "checkpoints" / f"v0_seed{seed}.pt"
        checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
        if (checkpoint["arm"] != "v0" or checkpoint["seed"] != seed
                or checkpoint["config"]["section_dim"] != 2
                or checkpoint["config"]["ridge"] != 1.0
                or checkpoint["development_seal"] != sha256(DEV / "manifest.json")):
            raise ValueError("v0 gauge checkpoint contract mismatch")
        model = MetaSectionRegressor(288, 2, 1.0).to(device)
        model.load_state_dict(checkpoint["model_state"])
        model.eval()
        basis = model.basis().detach().cpu().numpy()
        correct_all = tensors["val_pair"].cpu().numpy() @ basis
        wrong_all = tensors["val_wrong"].cpu().numpy() @ basis
        cluster_by_target = {row["target_id"]: row["protein_group_40"]
                             for row in read_jsonl_gz(DEV / "metaval_cells_without_labels.jsonl.gz")}
        cell_target = {row["cell_id"]: row["target_id"]
                       for row in read_jsonl_gz(DEV / "metaval_cells_without_labels.jsonl.gz")}
        # One map per held-out cluster tests whether a global rotation transfers.
        for heldout in sorted(set(cluster_by_target.values())):
            train_indices = [index for cell_id, index in val_index.items()
                             if cluster_by_target[cell_target[cell_id]] != heldout]
            test_indices = [index for cell_id, index in val_index.items()
                            if cluster_by_target[cell_target[cell_id]] == heldout]
            rotation = procrustes_from_support(wrong_all[train_indices], correct_all[train_indices])
            global_rows.append({"seed": seed, "cluster": heldout,
                                "heldout_query_residual": relative_error(
                                    wrong_all[test_indices] @ rotation, correct_all[test_indices])})
        for episode in episodes:
            if episode["seed"] != seed:
                continue
            support = np.asarray([val_index[cell_id] for cell_id in episode["support_cell_ids"]])
            query = np.asarray([val_index[cell_id] for cell_id in episode["query_cell_ids"]])
            ligand_s, ligand_q = tensors["val_ligand"][support], tensors["val_ligand"][query]
            pair_cs, pair_cq = tensors["val_pair"][support], tensors["val_pair"][query]
            pair_ws, pair_wq = tensors["val_wrong"][support], tensors["val_wrong"][query]
            with torch.no_grad():
                mu_cs, m_cs_t = model.components(ligand_s, pair_cs)
                mu_cq, m_cq_t = model.components(ligand_q, pair_cq)
                mu_ws, m_ws_t = model.components(ligand_s, pair_ws)
                mu_wq, m_wq_t = model.components(ligand_q, pair_wq)
            m_cs, m_cq = m_cs_t.cpu().numpy(), m_cq_t.cpu().numpy()
            m_ws, m_wq = m_ws_t.cpu().numpy(), m_wq_t.cpu().numpy()
            mu_cs, mu_cq = mu_cs.cpu().numpy(), mu_cq.cpu().numpy()
            mu_ws, mu_wq = mu_ws.cpu().numpy(), mu_wq.cpu().numpy()
            support_y = np.asarray([(value - y_mean) / y_scale for value in episode["support_pK"]])
            hcc = ridge_kernel(m_cs, m_cq, 1.0)
            hcw = ridge_kernel(m_cs, m_wq, 1.0)
            hwc = ridge_kernel(m_ws, m_cq, 1.0)
            hww = ridge_kernel(m_ws, m_wq, 1.0)
            pred = {
                "v0_correct": mu_cq + hcc @ (support_y - mu_cs),
                "v0_sc_qw": mu_wq + hcw @ (support_y - mu_cs),
                "v0_sw_qc": mu_cq + hwc @ (support_y - mu_ws),
                "v0_sw_qw": mu_wq + hww @ (support_y - mu_ws),
            }
            for arm, values in pred.items():
                stored = np.asarray([predictions[(arm, seed, episode["target_id"], episode["draw"], cell_id)]
                                     for cell_id in episode["query_cell_ids"]])
                max_recompute_error = max(max_recompute_error, float(np.abs(values - stored).max()))
            rotation = procrustes_from_support(m_ws, m_cs)
            gram_c, gram_w = m_cs @ m_cs.T, m_ws @ m_ws.T
            eigen_c = np.linalg.eigvalsh(gram_c)
            eigen_w = np.linalg.eigvalsh(gram_w)
            permuted = np.asarray([predictions[("v0_permuted", seed, episode["target_id"],
                                               episode["draw"], cell_id)]
                                   for cell_id in episode["query_cell_ids"]])
            episode_rows.append({
                "seed": seed, "target_id": episode["target_id"],
                "protein_group_40": episode["protein_group_40"], "draw": episode["draw"],
                "rank_correct": int(np.linalg.matrix_rank(m_cs)),
                "rank_wrong": int(np.linalg.matrix_rank(m_ws)),
                "support_procrustes": relative_error(m_ws @ rotation, m_cs),
                "query_transfer": relative_error(m_wq @ rotation, m_cq),
                "support_norm_ratio_wrong_correct": float(np.linalg.norm(m_ws) / (np.linalg.norm(m_cs) + 1e-12)),
                "gram_relative": relative_error(gram_w, gram_c),
                "gram_trace_ratio": float(np.trace(gram_w) / (np.trace(gram_c) + 1e-12)),
                "gram_eigen_relative": relative_error(eigen_w, eigen_c),
                "effective_ridge_correct": float(1.0 / (np.trace(gram_c) / len(m_cs) + 1e-12)),
                "effective_ridge_wrong": float(1.0 / (np.trace(gram_w) / len(m_ws) + 1e-12)),
                "H_CC_WW_relative": relative_error(hww, hcc),
                "H_CC_CW_relative": relative_error(hcw, hcc),
                "H_CC_WC_relative": relative_error(hwc, hcc),
                "mu_support_abs": float(np.mean(np.abs(mu_ws - mu_cs))),
                "mu_query_abs": float(np.mean(np.abs(mu_wq - mu_cq))),
                "section_CC_WW_abs": float(np.mean(np.abs(
                    hww @ (support_y - mu_ws) - hcc @ (support_y - mu_cs)))),
                "prediction_CC_WW_abs": float(np.mean(np.abs(pred["v0_sw_qw"] - pred["v0_correct"]))),
                "prediction_CC_permuted_abs": float(np.mean(np.abs(permuted - pred["v0_correct"]))),
            })
    numeric_keys = [key for key in episode_rows[0]
                    if key not in {"seed", "target_id", "protein_group_40", "draw"}]
    summary = {key: summarize([row[key] for row in episode_rows]) for key in numeric_keys}
    cluster_values = defaultdict(lambda: defaultdict(list))
    for row in episode_rows:
        for key in ("H_CC_WW_relative", "prediction_CC_WW_abs", "query_transfer", "gram_relative"):
            cluster_values[row["protein_group_40"]][key].append(row[key])
    cluster_rows = [{"cluster": cluster, **{key: float(np.mean(value)) for key, value in values.items()}}
                    for cluster, values in sorted(cluster_values.items())]
    relationship = {
        key: float(stats.spearmanr(
            [row[key] for row in cluster_rows],
            [row["prediction_CC_WW_abs"] for row in cluster_rows]).statistic)
        for key in ("H_CC_WW_relative", "query_transfer", "gram_relative")
    }
    if max_recompute_error > 2e-5:
        raise ValueError(f"prediction recomputation mismatch: {max_recompute_error}")
    global_targets = set(cluster_by_target)
    checkpoint_paths = [V1_OUT / "checkpoints" / f"v0_seed{seed}.pt" for seed in SEEDS]
    result = {
        "schema": "MetaSieve.A0GaugeAudit.v1",
        "TERMINAL_VERDICT": "NO_SINGLE_EXACT_ORTHOGONAL_GAUGE_IDENTIFIED",
        "query_affinity_values_used": 0,
        "episodes": len(episode_rows), "targets": len({row["target_id"] for row in episode_rows}),
        "clusters": len(cluster_rows), "max_prediction_recompute_error": max_recompute_error,
        "controls": controls(), "episode_summary": summary,
        "cluster_spearman_with_prediction_CC_WW_abs": relationship,
        "global_leave_one_cluster_out_query_residual": summarize(
            [row["heldout_query_residual"] for row in global_rows]),
        "global_audit_scope": {
            "cells": len(val_index),
            "targets": len(global_targets),
            "clusters": len(set(cluster_by_target.values())),
            "weighting": "cell-weighted Procrustes; descriptive only",
        },
        "interpretation_boundary": "descriptive local gauge-like audit; nine clusters; not confirmation",
        "inputs": {"prediction_sha256": sha256(V1_OUT / "predictions_before_query_labels.jsonl.gz"),
                   "episodes_sha256": sha256(DEV / "metaval_episodes.jsonl.gz"),
                   "correct_features_sha256": sha256(DEV / "metaval_features.npz"),
                   "wrong_features_sha256": sha256(DEV / "metaval_wrong_features.npz"),
                   "wrong_donor_map_sha256": sha256(DEV / "metaval_wrong_protein_map.json"),
                   "development_seal_sha256": sha256(DEV / "manifest.json"),
                   "checkpoint_sha256": {
                       path.name: sha256(path) for path in checkpoint_paths
                   }},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(audit(), indent=2, sort_keys=True))
