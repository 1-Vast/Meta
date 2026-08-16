"""Cross-fitted linear source router for the A2S-DTA primary test."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import torch

from .a2s_baseline import (
    DEFAULT_OUT,
    FEATURES,
    PRIMARY_K,
    RECIPIENT_MAX,
    REGISTRY,
    SOURCE_MIN,
    TARGETS,
    _bootstrap,
    _cosine,
    _metrics,
    _normalized,
    build_episodes,
    file_sha256,
    fit_ridge,
    fit_source_adapters,
    load_units,
    predict,
    role_targets,
)


ROUTER_OUT = Path("reports/active/a2s_router_pki_seed1729.json")


def _target_vectors() -> dict[str, np.ndarray]:
    archive = np.load(TARGETS, allow_pickle=False)
    return {
        str(key): value.astype(np.float32)
        for key, value in zip(archive["keys"].astype(str), archive["pooled"])
    }


def _context(
    units: pd.DataFrame,
    features: np.ndarray,
    source_targets: set[str],
    recipient_targets: set[str],
    target_vectors: dict[str, np.ndarray],
    *,
    ridge: float,
) -> dict[str, object]:
    weights, mean, scale, ymean = fit_ridge(features, units, source_targets, ridge=ridge)
    episodes = {
        k: build_episodes(units, recipient_targets, k=k)
        for k in PRIMARY_K
    }
    source_rows = units.loc[units.target.isin(source_targets), "source_row"].to_numpy(np.int64)
    all_rows = set(source_rows.tolist())
    for values in episodes.values():
        for episode in values:
            all_rows.update(episode.support)
            all_rows.update(episode.query)
    needed = np.asarray(sorted(all_rows), dtype=np.int64)
    base = predict(features, needed, weights, mean, scale, ymean)
    base_by_row = {int(row): float(value) for row, value in zip(needed, base)}
    adapters = fit_source_adapters(units, source_targets, base_by_row)
    source_order = sorted(source_targets)
    adapter_matrix = np.asarray([adapters[target] for target in source_order], dtype=np.float64)
    centroids = []
    depths = units.groupby("target", sort=True).size().to_dict()
    for target in source_order:
        rows = units.loc[units.target == target, "source_row"].to_numpy(np.int64)
        value = _normalized(features, rows, mean, scale).mean(axis=0)
        centroids.append(value)
    centroids = np.asarray(centroids, dtype=np.float32)
    centroids /= np.linalg.norm(centroids, axis=1, keepdims=True) + 1e-8
    protein = np.asarray(
        [target_vectors.get(target, np.zeros(1280, dtype=np.float32)) for target in source_order]
    )
    protein /= np.linalg.norm(protein, axis=1, keepdims=True) + 1e-8
    return {
        "episodes": episodes,
        "base_by_row": base_by_row,
        "adapters": adapter_matrix,
        "source_order": source_order,
        "source_centroids": centroids,
        "source_protein": protein,
        "depths": depths,
        "mean": mean,
        "scale": scale,
    }


def _candidate_features(
    *,
    units: pd.DataFrame,
    features: np.ndarray,
    episode_target: str,
    support: np.ndarray,
    source_order: list[str],
    adapters: np.ndarray,
    base_by_row: dict[int, float],
    source_centroids: np.ndarray,
    source_protein: np.ndarray,
    mean: np.ndarray,
    scale: np.ndarray,
    depths: dict[str, int],
    target_vectors: dict[str, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    labels = units.set_index("source_row").loc[support, "affinity"].to_numpy(np.float64)
    support_base = np.asarray([base_by_row[int(row)] for row in support], dtype=np.float64)
    adapter_prediction = adapters[:, 0, None] + adapters[:, 1, None] * support_base[None, :]
    support_mse = np.mean(np.square(adapter_prediction - labels[None, :]), axis=1)
    support_feature = _normalized(features, support, mean, scale).mean(axis=0)
    recipient_protein = target_vectors.get(episode_target, np.zeros(1280, dtype=np.float32))
    protein = _cosine(recipient_protein, source_protein)
    chemistry = _cosine(support_feature, source_centroids)
    depth = np.asarray([np.log1p(depths[target]) for target in source_order], dtype=np.float64)
    matrix = np.column_stack((support_mse, protein, chemistry, depth))
    return matrix, support_base


def _fit_router(matrix: np.ndarray, utility: np.ndarray, ridge: float = 1.0) -> dict[str, np.ndarray]:
    if matrix.ndim != 2 or len(matrix) != len(utility):
        raise ValueError("router design and utility rows must align")
    center = matrix.mean(axis=0)
    scale = matrix.std(axis=0)
    scale[scale < 1e-8] = 1.0
    design = (matrix - center) / scale
    design = np.column_stack((np.ones(len(design)), design))
    penalty = np.eye(design.shape[1], dtype=np.float64) * ridge
    penalty[0, 0] = 1e-6
    weight = np.linalg.solve(design.T @ design + penalty, design.T @ utility)
    return {"center": center, "scale": scale, "weight": weight}


def _router_score(matrix: np.ndarray, router: dict[str, np.ndarray]) -> np.ndarray:
    design = (matrix - router["center"]) / router["scale"]
    return np.column_stack((np.ones(len(design)), design)) @ router["weight"]


def _score_episode(
    *,
    units: pd.DataFrame,
    features: np.ndarray,
    episode,
    context: dict[str, object],
    target_vectors: dict[str, np.ndarray],
    router: dict[str, np.ndarray],
    rng: np.random.Generator,
) -> tuple[dict[str, dict[str, float]], dict[str, object]]:
    source_order = context["source_order"]
    support = np.asarray(episode.support, dtype=np.int64)
    query = np.asarray(episode.query, dtype=np.int64)
    labels_by_row = units.set_index("source_row")
    support_label = labels_by_row.loc[support, "affinity"].to_numpy(np.float64)
    query_label = labels_by_row.loc[query, "affinity"].to_numpy(np.float64)
    base_by_row = context["base_by_row"]
    support_base = np.asarray([base_by_row[int(row)] for row in support], dtype=np.float64)
    query_base = np.asarray([base_by_row[int(row)] for row in query], dtype=np.float64)
    offset = float(np.mean(support_label - support_base))
    calibration = query_base + offset
    matrix, _ = _candidate_features(
        units=units,
        features=features,
        episode_target=episode.target,
        support=support,
        source_order=source_order,
        adapters=context["adapters"],
        base_by_row=base_by_row,
        source_centroids=context["source_centroids"],
        source_protein=context["source_protein"],
        mean=context["mean"],
        scale=context["scale"],
        depths=context["depths"],
        target_vectors=target_vectors,
    )
    router_utility = _router_score(matrix, router)
    selected = int(np.argmax(router_utility))
    random_source = int(rng.integers(0, len(source_order)))
    adapter = context["adapters"]
    query_adapter = adapter[:, 0, None] + adapter[:, 1, None] * query_base[None, :]
    predictions = {
        "recipient_calibration": calibration,
        "source_support": query_adapter[int(np.argmin(matrix[:, 0]))],
        "source_random": query_adapter[random_source],
        "router": query_adapter[selected],
        "router_gated": (
            query_adapter[selected] if router_utility[selected] > 0.0 else calibration
        ),
    }
    scores = {name: _metrics(query_label, value) for name, value in predictions.items()}
    record = {
        "target": episode.target,
        "k": episode.k,
        "query_units": len(query),
        "selected_source": source_order[selected],
        "selected_source_utility": float(router_utility[selected]),
        "random_source": source_order[random_source],
        "gate": float(router_utility[selected] > 0.0),
        "router_rmse": scores["router"]["rmse"],
        "calibration_rmse": scores["recipient_calibration"]["rmse"],
        "router_gated_rmse": scores["router_gated"]["rmse"],
    }
    return scores, record


def _summarize(per_arm: dict[str, list[dict[str, float]]], rng: np.random.Generator) -> dict[str, object]:
    calibration = per_arm["recipient_calibration"]
    result: dict[str, object] = {}
    for arm, scores in per_arm.items():
        gains = np.asarray(
            [base["rmse"] - score["rmse"] for base, score in zip(calibration, scores)],
            dtype=np.float64,
        )
        result[arm] = {
            metric: float(np.nanmean([score[metric] for score in scores]))
            for metric in ("rmse", "mae", "spearman", "pairwise_accuracy", "ndcg_at_10")
        }
        result[arm]["rmse_gain_vs_calibration"] = float(np.nanmean(gains))
        result[arm]["negative_transfer_rate_vs_calibration"] = float(np.mean(gains < 0.0))
        result[arm]["benefiting_recipient_rate_vs_calibration"] = float(np.mean(gains > 0.0))
        result[arm]["rmse_gain_ci95"] = _bootstrap(gains, rng)
    return result


def run(
    *,
    endpoint: str = "pKi",
    seed: int = 1729,
    output: Path = ROUTER_OUT,
    ridge: float = 10.0,
    router_ridge: float = 1.0,
    max_recipients: int = 0,
) -> dict[str, object]:
    if not torch.cuda.is_available():
        raise RuntimeError("A2S router requires CUDA in the drug environment")
    started = time.perf_counter()
    rng = np.random.default_rng(seed)
    units = load_units(endpoint)
    source_targets, natural_recipients, counts = role_targets(units)
    features_archive = np.load(FEATURES, allow_pickle=False)
    features = features_archive["feat"]
    vectors = _target_vectors()
    source_order = sorted(source_targets)
    folds = [set(source_order[index::3]) for index in range(3)]
    meta_matrix: list[np.ndarray] = []
    meta_utility: list[float] = []
    meta_episodes = 0
    for fold in folds:
        train_sources = source_targets.difference(fold)
        context = _context(units, features, train_sources, fold, vectors, ridge=ridge)
        labels_by_row = units.set_index("source_row")
        for k, episodes in context["episodes"].items():
            for episode in episodes:
                support = np.asarray(episode.support, dtype=np.int64)
                query = np.asarray(episode.query, dtype=np.int64)
                query_label = labels_by_row.loc[query, "affinity"].to_numpy(np.float64)
                base = context["base_by_row"]
                query_base = np.asarray([base[int(row)] for row in query], dtype=np.float64)
                support_label = labels_by_row.loc[support, "affinity"].to_numpy(np.float64)
                support_base = np.asarray([base[int(row)] for row in support], dtype=np.float64)
                calibration = query_base + float(np.mean(support_label - support_base))
                matrix, _ = _candidate_features(
                    units=units,
                    features=features,
                    episode_target=episode.target,
                    support=support,
                    source_order=context["source_order"],
                    adapters=context["adapters"],
                    base_by_row=base,
                    source_centroids=context["source_centroids"],
                    source_protein=context["source_protein"],
                    mean=context["mean"],
                    scale=context["scale"],
                    depths=context["depths"],
                    target_vectors=vectors,
                )
                query_adapter = context["adapters"][:, 0, None] + context["adapters"][:, 1, None] * query_base[None, :]
                utility = np.asarray(
                    [
                        float(np.sqrt(np.mean(np.square(query_label - calibration))) - np.sqrt(np.mean(np.square(query_label - prediction))))
                        for prediction in query_adapter
                    ],
                    dtype=np.float64,
                )
                meta_matrix.append(matrix)
                meta_utility.extend(utility.tolist())
                meta_episodes += 1
    design = np.vstack(meta_matrix)
    utility = np.asarray(meta_utility, dtype=np.float64)
    router = _fit_router(design, utility, ridge=router_ridge)

    full_context = _context(units, features, source_targets, natural_recipients, vectors, ridge=ridge)
    metrics_by_k: dict[str, object] = {}
    records: list[dict[str, object]] = []
    for k, episodes in full_context["episodes"].items():
        if max_recipients:
            episodes = episodes[:max_recipients]
        per_arm: dict[str, list[dict[str, float]]] = {
            "recipient_calibration": [],
            "source_support": [],
            "source_random": [],
            "router": [],
            "router_gated": [],
        }
        for episode in episodes:
            scores, record = _score_episode(
                units=units,
                features=features,
                episode=episode,
                context=full_context,
                target_vectors=vectors,
                router=router,
                rng=rng,
            )
            for arm in per_arm:
                per_arm[arm].append(scores[arm])
            records.append(record)
        metrics_by_k[str(k)] = {
            "recipients": len(episodes),
            "summary": _summarize(per_arm, rng),
        }

    learning_curve: dict[str, dict[str, float]] = {}
    ks = np.asarray(PRIMARY_K, dtype=np.float64)
    for arm in ("recipient_calibration", "source_support", "source_random", "router", "router_gated"):
        gains = np.asarray(
            [metrics_by_k[str(k)]["summary"][arm]["rmse_gain_vs_calibration"] for k in PRIMARY_K],
            dtype=np.float64,
        )
        learning_curve[arm] = {
            "aulc_rmse_gain": float(np.trapz(gains, ks) / (ks[-1] - ks[0])),
            "mean_benefiting_recipient_rate": float(
                np.mean(
                    [
                        metrics_by_k[str(k)]["summary"][arm][
                            "benefiting_recipient_rate_vs_calibration"
                        ]
                        for k in PRIMARY_K
                    ]
                )
            ),
        }

    result = {
        "schema_version": "a2s-crossfitted-router-v1",
        "protocol": "source-fold cross-fitted recipient-conditioned linear router",
        "decision": "ROUTER_DIAGNOSTIC_ONLY",
        "endpoint": endpoint,
        "seed": seed,
        "source_rule": f"n_eff >= {SOURCE_MIN}",
        "recipient_rule": f"n_eff < {RECIPIENT_MAX}",
        "source_targets": len(source_targets),
        "recipient_candidates": len(natural_recipients),
        "primary_support_sizes": list(PRIMARY_K),
        "query_track": "target-side single-cold",
        "meta_training": {
            "source_folds": 3,
            "pseudo_recipient_episodes": meta_episodes,
            "router_rows": int(len(utility)),
            "source_fit_weighting": "target_macro_equal_total_weight",
            "features": ["support_mse", "protein_cosine", "chemistry_cosine", "log_source_depth"],
            "target_labels_used": "source-fold pseudo-recipients only",
        },
        "input": {
            "registry_sha256": file_sha256(REGISTRY),
            "features_sha256": file_sha256(FEATURES),
            "target_features_sha256": file_sha256(TARGETS),
            "registry_rows": int(len(features)),
            "feature_row_alignment": "global_registry_source_row",
            "feature_conn_sha_verified": True,
            "split": "train",
        },
        "metrics_by_k": metrics_by_k,
        "learning_curve": learning_curve,
        "recipient_records": records,
        "compute": {
            "device": torch.cuda.get_device_name(0),
            "peak_torch_memory_mib": float(torch.cuda.max_memory_allocated() / 2**20),
            "wall_seconds": time.perf_counter() - started,
            "gradient_evaluations": 0,
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--endpoint", choices=("pKi", "pKd"), default="pKi")
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--max-recipients", type=int, default=0)
    parser.add_argument("--ridge", type=float, default=10.0)
    parser.add_argument("--router-ridge", type=float, default=1.0)
    parser.add_argument("--out", type=Path, default=ROUTER_OUT)
    args = parser.parse_args(argv)
    result = run(
        endpoint=args.endpoint,
        seed=args.seed,
        output=args.out,
        ridge=args.ridge,
        router_ridge=args.router_ridge,
        max_recipients=args.max_recipients,
    )
    print(
        json.dumps(
            {
                "output": str(args.out),
                "endpoint": args.endpoint,
                "source_targets": result["source_targets"],
                "recipient_candidates": result["recipient_candidates"],
                "recipients_by_k": {
                    k: value["recipients"] for k, value in result["metrics_by_k"].items()
                },
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
