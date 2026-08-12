"""R2-E0/E1: post-hoc calibration and support-weight audit for main-v0.

This is a descriptive audit on already-consumed splits. It cannot open a Gate.
E0 compares the learned section with matched pair-prior and ligand-only
support-mean calibrators. E1 inspects the label-free support-weight operator.
The output is split-specific so a meta-val run cannot be overwritten by a
meta-test run.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch

from research.meta_fewshot.train_main_v0 import (
    OUT,
    SEEDS,
    MetaSectionRegressor,
    TrainConfig,
    draw_episode,
    load_data,
    sha256,
)
from research.meta_fewshot.r2_calibration_orthogonal_section import (
    calibration_orthogonal_prediction,
    ridge_support_weights,
)

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = "MetaSieve.R2RegimeAudit.v2"


def report_path(split: str) -> Path:
    return ROOT / f"report/meta_fewshot/r2_regime_audit_{split}.json"


def episode_rng(seed: int, target: str, draw: int) -> np.random.Generator:
    key = f"main-v0-episode|{seed}|{target}|{draw}"
    return np.random.default_rng(int(hashlib.sha256(key.encode()).hexdigest()[:16], 16))


def support_weights(coordinates_s: torch.Tensor, coordinates_q: torch.Tensor,
                    ridge: float) -> torch.Tensor:
    """Compatibility wrapper for the shared audited ridge primitive."""
    return ridge_support_weights(coordinates_s, coordinates_q, ridge)


def decompose(u: torch.Tensor) -> dict[str, torch.Tensor]:
    k = u.shape[-1]
    direction = torch.ones(k, dtype=u.dtype, device=u.device) / np.sqrt(k)
    parallel = (u @ direction)[:, None] * direction[None, :]
    orthogonal = u - parallel
    return {
        "calibration_norm": parallel.norm(dim=-1),
        "specific_norm": orthogonal.norm(dim=-1),
        "weight_sum": u.sum(dim=-1),
    }


def aggregate(per_target: dict, target_cluster: dict[str, str], fields: tuple[str, ...]) -> dict:
    target_mean = {
        target: {field: float(np.mean(values[field])) for field in fields}
        for target, values in per_target.items()
    }
    clusters = defaultdict(lambda: defaultdict(list))
    for target, values in target_mean.items():
        for field in fields:
            clusters[target_cluster[target]][field].append(values[field])
    cluster_mean = {
        cluster: {field: float(np.mean(values[field])) for field in fields}
        for cluster, values in clusters.items()
    }
    return {
        "target_macro": {
            field: float(np.mean([values[field] for values in target_mean.values()]))
            for field in fields
        },
        "cluster_macro": {
            field: float(np.mean([values[field] for values in cluster_mean.values()]))
            for field in fields
        },
        "per_target": target_mean,
        "per_cluster": cluster_mean,
        "n_clusters": len(cluster_mean),
    }


def load_checkpoint(model: MetaSectionRegressor, path: Path, seed: int, arm: str,
                    d: int, ridge: float) -> dict:
    state = torch.load(path, map_location=next(model.parameters()).device, weights_only=False)
    if (state["seed"] != seed or state["arm"] != arm or state["d"] != d
            or float(state["ridge"]) != ridge):
        raise ValueError(f"checkpoint contract mismatch: {path}")
    model.load_state_dict(state["model_state"])
    model.eval()
    return state


def run(split: str, device: str) -> dict:
    existing_path = OUT / "MAIN_V0_RESULT.json"
    existing = json.loads(existing_path.read_text(encoding="utf-8"))
    config = TrainConfig(**existing["config"])
    d, ridge = int(existing["selected"]["d"]), float(existing["selected"]["ridge"])
    cells, tensors, tasks, norm = load_data(device)
    if split not in tasks:
        raise ValueError(f"split {split} absent from corpus")
    targets = sorted(tasks[split])
    target_cluster = {}
    for row in cells:
        if row["target_id"] in tasks[split]:
            target_cluster[row["target_id"]] = row["protein_group_40"]
    if set(target_cluster) != set(targets):
        raise ValueError("target-to-cluster mapping is incomplete")

    loss_fields = (
        "full", "pair_intercept", "centered_section_posthoc", "pair_zero",
        "ligand_intercept", "ligand_zero",
    )
    structural_fields = (
        "calibration_share", "u_gap_wrong", "u_gap_uniform",
        "effective_ridge_ratio",
    )
    per_target = {
        target: {field: [] for field in loss_fields + structural_fields}
        for target in targets
    }
    checkpoint_hashes = {}

    for seed in SEEDS:
        full_model = MetaSectionRegressor(288, d, ridge).to(device)
        population_model = MetaSectionRegressor(288, 0, ridge).to(device)
        full_path = OUT / "checkpoints" / f"full_seed{seed}.pt"
        population_path = OUT / "checkpoints" / f"population_seed{seed}.pt"
        load_checkpoint(full_model, full_path, seed, "full", d, ridge)
        load_checkpoint(population_model, population_path, seed, "population", 0, ridge)
        checkpoint_hashes[full_path.name] = sha256(full_path)
        checkpoint_hashes[population_path.name] = sha256(population_path)

        with torch.no_grad():
            for target in targets:
                for draw in range(config.test_draws):
                    support, query = draw_episode(
                        tasks[split][target], episode_rng(seed, target, draw), config.k)
                    y_s, y_q = tensors["y"][support], tensors["y"][query]
                    ligand_s, ligand_q = tensors["ligand"][support], tensors["ligand"][query]
                    mu_s, coord_s = full_model.components(ligand_s, tensors["correct"][support])
                    mu_q, coord_q = full_model.components(ligand_q, tensors["correct"][query])
                    ligand_mu_s, _ = population_model.components(
                        ligand_s, tensors["correct"][support])
                    ligand_mu_q, _ = population_model.components(
                        ligand_q, tensors["correct"][query])
                    residual = y_s - mu_s
                    ligand_residual = y_s - ligand_mu_s

                    u_correct = support_weights(coord_s, coord_q, ridge)
                    full = mu_q + u_correct @ residual
                    pair_intercept = mu_q + residual.mean()
                    pair_zero = mu_q
                    ligand_intercept = ligand_mu_q + ligand_residual.mean()
                    ligand_zero = ligand_mu_q

                    centered_section, _, _ = calibration_orthogonal_prediction(
                        mu_s, coord_s, y_s, mu_q, coord_q, ridge)

                    predictions = {
                        "full": full,
                        "pair_intercept": pair_intercept,
                        "centered_section_posthoc": centered_section,
                        "pair_zero": pair_zero,
                        "ligand_intercept": ligand_intercept,
                        "ligand_zero": ligand_zero,
                    }
                    for name, prediction in predictions.items():
                        per_target[target][name].append(
                            float(((prediction - y_q) ** 2).mean()) * norm["y_scale"] ** 2)

                    _, wrong_s = full_model.components(ligand_s, tensors["wrong"][support])
                    _, wrong_q = full_model.components(ligand_q, tensors["wrong"][query])
                    u_wrong = support_weights(wrong_s, wrong_q, ridge)
                    parts = decompose(u_correct)
                    share = parts["calibration_norm"] / (
                        parts["calibration_norm"] + parts["specific_norm"] + 1e-12)
                    uniform = torch.full_like(u_correct, 1.0 / config.k)
                    gram = coord_s @ coord_s.T
                    per_target[target]["calibration_share"].append(float(share.mean()))
                    per_target[target]["u_gap_wrong"].append(
                        float((u_correct - u_wrong).norm(dim=-1).mean()))
                    per_target[target]["u_gap_uniform"].append(
                        float((u_correct - uniform).norm(dim=-1).mean()))
                    per_target[target]["effective_ridge_ratio"].append(
                        float(ridge / (torch.trace(gram).item() / config.k + 1e-12)))

    aggregated = aggregate(per_target, target_cluster, loss_fields + structural_fields)
    result = {
        "schema": SCHEMA,
        "split": split,
        "declared_role": "POSTHOC_DESCRIPTIVE_DIAGNOSIS_NOT_CONFIRMATORY",
        "query_labels_used_for_E0": True,
        "query_labels_used_for_E1_support_weight_metrics": False,
        "all_corpus_labels_loaded_by_legacy_loader": True,
        "d": d,
        "ridge": ridge,
        "k": config.k,
        "n_targets": len(targets),
        "n_clusters": aggregated["n_clusters"],
        "metrics": {key: aggregated[key] for key in ("target_macro", "cluster_macro")},
        "per_target": aggregated["per_target"],
        "per_cluster": aggregated["per_cluster"],
        "inputs": {
            "main_result_sha256": sha256(existing_path),
            "checkpoint_sha256": checkpoint_hashes,
        },
        "interpretation": {
            "pair_intercept_is_matched_to_full_pair_prior": True,
            "ligand_intercept_uses_independently_trained_d0_population": True,
            "centered_section_posthoc_was_not_used_for_model_selection": True,
        },
    }

    target = result["metrics"]["target_macro"]
    cluster = result["metrics"]["cluster_macro"]
    result["E0"] = {}
    for unit, values in (("target_macro", target), ("cluster_macro", cluster)):
        total_gain = values["pair_zero"] - values["full"]
        noncalibration_gain = values["pair_intercept"] - values["full"]
        result["E0"][unit] = {
            "pair_intercept_minus_full": noncalibration_gain,
            "pair_zero_minus_full": total_gain,
            "noncalibration_share_of_pair_support_gain": (
                noncalibration_gain / total_gain if total_gain > 0 else None),
            "ligand_intercept_minus_full": values["ligand_intercept"] - values["full"],
        }

    calibration_share = target["calibration_share"]
    gauge_ratio = target["u_gap_wrong"] / (target["u_gap_uniform"] + 1e-12)
    result["E1"] = {
        "target_macro_calibration_share": calibration_share,
        "target_macro_gauge_ratio": gauge_ratio,
        "registered_near_uniform_threshold": 0.9,
        "registered_near_gauge_threshold": 0.5,
        "registered_H0_falsifier_gauge_ratio_gt": 1.0,
    }

    shares = [result["E0"][unit]["noncalibration_share_of_pair_support_gain"]
              for unit in ("target_macro", "cluster_macro")]
    if all(value is not None and value <= 0 for value in shares):
        e0_verdict = "META_SECTION_EFFECT_IS_CALIBRATION"
    elif all(value is not None and value < 0.10 for value in shares):
        e0_verdict = "META_SECTION_EFFECT_PREDOMINANTLY_CALIBRATION"
    else:
        e0_verdict = "META_SECTION_EFFECT_NOT_EXPLAINED_BY_CALIBRATION"
    if gauge_ratio > 1.0:
        e1_verdict = "H0_REGIME_FALSIFIED_BY_REGISTERED_E1"
    elif gauge_ratio <= 0.5:
        e1_verdict = "PROTEIN_CHANNEL_NEAR_GAUGE_EQUIVALENT"
    else:
        e1_verdict = "PROTEIN_CHANNEL_NOT_NEAR_GAUGE_EQUIVALENT"
    result["verdict"] = [e0_verdict, e1_verdict, "RFMS_TRAINING_NOT_AUTHORIZED"]
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split", choices=("meta_val", "meta_test"), default="meta_val")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()
    result = run(args.split, args.device)
    path = report_path(args.split)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, allow_nan=False, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")
    print(json.dumps({key: value for key, value in result.items()
                      if key not in {"per_target", "per_cluster"}}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
