"""Low-cost A2S-DTA source-only baseline and routing smoke."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import torch


ROOT = Path("dataset/public/chembl_37/processed/dualcold")
REGISTRY = ROOT / "registry.parquet"
FEATURES = ROOT / "ligand_features.npz"
TARGETS = ROOT / "target_esm2.npz"
DEFAULT_OUT = Path("reports/active/a2s_pki_seed1729.json")
SOURCE_MIN = 100
RECIPIENT_MAX = 30
PRIMARY_K = (1, 3, 5)
MIN_QUERY = 5
COLUMNS = [
    "target",
    "conn",
    "endpoint",
    "affinity",
    "scaffold",
    "assays",
    "docs",
    "accession",
    "hcluster",
    "dual_cold_split",
]


@dataclass(frozen=True)
class Episode:
    target: str
    k: int
    support: tuple[int, ...]
    query: tuple[int, ...]
    homology: str


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _tokens(value: object) -> set[str]:
    return {part.strip() for part in str(value).split("|") if part.strip()}


def connection_hash(frame: pd.DataFrame) -> str:
    """Hash registry connectivity in its canonical global row order."""

    values = pd.util.hash_pandas_object(frame["conn"], index=False).values
    return sha256(values.tobytes()).hexdigest()


def _units(frame: pd.DataFrame) -> pd.DataFrame:
    """Collapse duplicate parent/provenance rows without mixing endpoints."""

    frame = frame.copy()
    if "source_row" not in frame.columns:
        frame["source_row"] = frame.index.astype(np.int64)
    else:
        frame["source_row"] = pd.to_numeric(frame["source_row"], errors="raise").astype(
            np.int64
        )
    frame["target"] = frame.target.astype(str)
    frame["endpoint"] = frame.endpoint.astype(str)
    frame["affinity"] = pd.to_numeric(frame.affinity, errors="coerce")
    frame = frame[np.isfinite(frame.affinity.to_numpy())].copy()
    if frame.empty:
        raise ValueError("no finite affinity rows remain")
    keys = ["target", "endpoint", "conn", "docs", "assays"]
    return (
        frame.sort_values(keys + ["source_row"], kind="stable")
        .drop_duplicates(keys, keep="first")
        .reset_index(drop=True)
    )


def load_units(endpoint: str) -> pd.DataFrame:
    # Feature rows are aligned to the complete registry, so filtering first
    # would silently replace global row ids with a local RangeIndex.
    frame = pd.read_parquet(REGISTRY, columns=COLUMNS).reset_index(names="source_row")
    archive = np.load(FEATURES, allow_pickle=False)
    if len(archive["feat"]) != len(frame):
        raise ValueError("ligand feature cache and registry have different row counts")
    cached_hash = str(archive["conn_sha"].item())
    actual_hash = connection_hash(frame)
    if cached_hash != actual_hash:
        raise ValueError("ligand feature cache is not aligned to registry connectivity rows")
    frame = frame[(frame.dual_cold_split == "train") & (frame.endpoint == endpoint)].copy()
    if frame.empty:
        raise ValueError(f"no TRAIN rows for endpoint {endpoint}")
    if set(frame.dual_cold_split.astype(str)) != {"train"}:
        raise ValueError("A2S source fitting received a non-TRAIN row")
    return _units(frame)


def role_targets(units: pd.DataFrame) -> tuple[set[str], set[str], dict[str, int]]:
    counts = units.groupby("target", sort=True).size().astype(int).to_dict()
    source = {target for target, count in counts.items() if count >= SOURCE_MIN}
    recipient = {target for target, count in counts.items() if count < RECIPIENT_MAX}
    if source.intersection(recipient):
        raise AssertionError("source and recipient targets overlap")
    return source, recipient, counts


def build_episodes(
    units: pd.DataFrame,
    recipients: set[str],
    *,
    k: int,
    min_query: int = MIN_QUERY,
    max_recipients: int = 0,
) -> list[Episode]:
    episodes: list[Episode] = []
    for target in sorted(recipients):
        group = units.loc[units.target == target].sort_values(
            ["scaffold", "conn", "docs", "assays", "source_row"], kind="stable"
        )
        support: list[int] = []
        scaffolds: set[str] = set()
        parents: set[str] = set()
        for row in group.itertuples(index=False):
            if str(row.scaffold) in scaffolds or str(row.conn) in parents:
                continue
            support.append(int(row.source_row))
            scaffolds.add(str(row.scaffold))
            parents.add(str(row.conn))
            if len(support) == k:
                break
        if len(support) < k:
            continue
        support_set = set(support)
        query = tuple(
            int(row.source_row)
            for row in group.itertuples(index=False)
            if int(row.source_row) not in support_set
        )
        if len(query) < min_query:
            continue
        homology = str(group.hcluster.iloc[0])
        episodes.append(
            Episode(target=target, k=k, support=tuple(support), query=query, homology=homology)
        )
        if max_recipients and len(episodes) >= max_recipients:
            break
    return episodes


def _feature_stats(
    features: np.ndarray, rows: np.ndarray, weights: np.ndarray | None = None
) -> tuple[np.ndarray, np.ndarray]:
    if weights is None:
        weights = np.ones(len(rows), dtype=np.float64)
    if len(weights) != len(rows) or np.any(weights < 0.0):
        raise ValueError("feature-stat weights must align with nonnegative rows")
    weight_total = float(weights.sum())
    if weight_total <= 0.0:
        raise ValueError("feature-stat weights must have positive total")
    total = np.zeros(features.shape[1], dtype=np.float64)
    squares = np.zeros(features.shape[1], dtype=np.float64)
    for start in range(0, len(rows), 4096):
        values = np.asarray(features[rows[start : start + 4096]], dtype=np.float64)
        chunk_weights = weights[start : start + len(values)]
        total += (values * chunk_weights[:, None]).sum(axis=0)
        squares += (np.square(values) * chunk_weights[:, None]).sum(axis=0)
    mean = total / weight_total
    variance = np.maximum(squares / weight_total - np.square(mean), 1e-8)
    return mean.astype(np.float32), np.sqrt(variance).astype(np.float32)


def _normalized(features: np.ndarray, rows: np.ndarray, mean: np.ndarray, scale: np.ndarray) -> np.ndarray:
    return (np.asarray(features[rows], dtype=np.float32) - mean) / scale


def fit_ridge(
    features: np.ndarray,
    frame: pd.DataFrame,
    source_targets: set[str],
    *,
    ridge: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    source_rows = frame.loc[frame.target.isin(source_targets), "source_row"].to_numpy(np.int64)
    labels = frame.loc[frame.target.isin(source_targets), "affinity"].to_numpy(np.float32)
    if len(source_rows) < 100:
        raise ValueError("source pool is too small for the pooled ridge")
    source_frame = frame.loc[frame.target.isin(source_targets), ["target", "source_row"]]
    target_counts = source_frame.groupby("target", sort=False).size().to_dict()
    target_count = len(target_counts)
    row_weights = np.asarray(
        [len(source_rows) / (target_count * target_counts[str(target)]) for target in source_frame.target],
        dtype=np.float64,
    )
    mean, scale = _feature_stats(features, source_rows, row_weights)
    device = torch.device("cuda")
    dimension = features.shape[1]
    xtx = torch.zeros((dimension, dimension), dtype=torch.float32, device=device)
    xty = torch.zeros(dimension, dtype=torch.float32, device=device)
    ymean = float(np.average(labels, weights=row_weights))
    eye = torch.eye(dimension, dtype=torch.float32, device=device)
    centered = labels - ymean
    for start in range(0, len(source_rows), 4096):
        rows = source_rows[start : start + 4096]
        x = torch.as_tensor(_normalized(features, rows, mean, scale), device=device)
        sqrt_weight = torch.as_tensor(
            np.sqrt(row_weights[start : start + len(rows)]), device=device, dtype=torch.float32
        )
        weighted_x = x * sqrt_weight[:, None]
        weighted_y = torch.as_tensor(centered[start : start + len(rows)], device=device) * sqrt_weight
        xtx.addmm_(weighted_x.T, weighted_x)
        xty.add_(weighted_x.T @ weighted_y)
    weights = torch.linalg.solve(xtx + ridge * eye, xty).cpu().numpy()
    return weights, mean, scale, ymean


def predict(
    features: np.ndarray,
    rows: np.ndarray,
    weights: np.ndarray,
    mean: np.ndarray,
    scale: np.ndarray,
    ymean: float,
) -> np.ndarray:
    result: list[np.ndarray] = []
    device = torch.device("cuda")
    weight = torch.as_tensor(weights, dtype=torch.float32, device=device)
    for start in range(0, len(rows), 4096):
        x = torch.as_tensor(
            _normalized(features, rows[start : start + 4096], mean, scale), device=device
        )
        result.append((x @ weight).cpu().numpy() + ymean)
    return np.concatenate(result) if result else np.empty(0, dtype=np.float32)


def fit_source_adapters(
    frame: pd.DataFrame,
    source_targets: set[str],
    base_by_row: dict[int, float],
) -> dict[str, tuple[float, float]]:
    adapters: dict[str, tuple[float, float]] = {}
    for target in sorted(source_targets):
        group = frame.loc[frame.target == target]
        base = np.asarray([base_by_row[int(row)] for row in group.source_row], dtype=np.float64)
        label = group.affinity.to_numpy(np.float64)
        matrix = np.column_stack((np.ones(len(base)), base))
        penalty = np.diag([1e-6, 1e-3])
        theta = np.linalg.solve(matrix.T @ matrix + penalty, matrix.T @ label)
        adapters[target] = (float(theta[0]), float(theta[1]))
    return adapters


def _cosine(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left_norm = np.linalg.norm(left) + 1e-8
    right_norm = np.linalg.norm(right, axis=1) + 1e-8
    return (right @ left) / (right_norm * left_norm)


def _metrics(label: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    error = label - prediction
    actual_rank = pd.Series(label).rank(method="average").to_numpy()
    predicted_rank = pd.Series(prediction).rank(method="average").to_numpy()
    spearman = (
        float(np.corrcoef(actual_rank, predicted_rank)[0, 1])
        if len(label) > 1 and np.std(actual_rank) > 0 and np.std(predicted_rank) > 0
        else float("nan")
    )
    left, right = np.triu_indices(len(label), 1)
    actual_diff = label[left] - label[right]
    predicted_diff = prediction[left] - prediction[right]
    usable = actual_diff != 0
    pairwise = (
        float(np.mean((actual_diff[usable] * predicted_diff[usable]) > 0))
        if np.any(usable)
        else float("nan")
    )
    cutoff = min(10, len(label))
    if cutoff == 0 or np.ptp(label) == 0.0:
        ndcg = float("nan")
    else:
        relevance = np.maximum(label - np.min(label), 0.0)
        discounts = 1.0 / np.log2(np.arange(2, cutoff + 2, dtype=np.float64))
        predicted_order = np.argsort(-prediction, kind="stable")[:cutoff]
        ideal_order = np.argsort(-label, kind="stable")[:cutoff]
        dcg = float(np.sum(relevance[predicted_order] * discounts))
        ideal = float(np.sum(relevance[ideal_order] * discounts))
        ndcg = dcg / ideal if ideal > 0.0 else float("nan")
    return {
        "rmse": float(np.sqrt(np.mean(np.square(error)))),
        "mae": float(np.mean(np.abs(error))),
        "spearman": spearman,
        "pairwise_accuracy": pairwise,
        "ndcg_at_10": ndcg,
    }


def _bootstrap(values: np.ndarray, rng: np.random.Generator) -> list[float]:
    values = values[np.isfinite(values)]
    if len(values) < 2:
        return [float("nan"), float("nan")]
    samples = values[rng.integers(0, len(values), size=(1000, len(values)))].mean(axis=1)
    return [float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))]


def run(
    *,
    endpoint: str,
    seed: int,
    output: Path,
    max_recipients: int,
    ridge: float,
) -> dict[str, object]:
    if not torch.cuda.is_available():
        raise RuntimeError("A2S numerical baseline requires CUDA in the drug environment")
    started = time.perf_counter()
    rng = np.random.default_rng(seed)
    units = load_units(endpoint)
    source_targets, recipients, counts = role_targets(units)
    features_archive = np.load(FEATURES, allow_pickle=False)
    features = features_archive["feat"]
    if len(features) <= int(units.source_row.max()):
        raise ValueError("ligand feature cache is not aligned to registry source rows")
    weights, mean, scale, ymean = fit_ridge(features, units, source_targets, ridge=ridge)

    source_rows = units.loc[units.target.isin(source_targets), "source_row"].to_numpy(np.int64)
    episodes_by_k = {
        k: build_episodes(units, recipients, k=k, max_recipients=max_recipients)
        for k in PRIMARY_K
    }
    all_rows = set(source_rows.tolist())
    for episodes in episodes_by_k.values():
        for episode in episodes:
            all_rows.update(episode.support)
            all_rows.update(episode.query)
    needed = np.asarray(sorted(all_rows), dtype=np.int64)
    base_values = predict(features, needed, weights, mean, scale, ymean)
    base_by_row = {int(row): float(value) for row, value in zip(needed, base_values)}
    adapters = fit_source_adapters(units, source_targets, base_by_row)

    source_order = sorted(source_targets)
    source_adapter = np.asarray([adapters[target] for target in source_order], dtype=np.float64)
    source_centroids = []
    for target in source_order:
        rows = units.loc[units.target == target, "source_row"].to_numpy(np.int64)
        source_centroids.append(_normalized(features, rows, mean, scale).mean(axis=0))
    source_centroids = np.asarray(source_centroids, dtype=np.float32)
    source_centroids /= np.linalg.norm(source_centroids, axis=1, keepdims=True) + 1e-8

    target_archive = np.load(TARGETS, allow_pickle=False)
    target_keys = target_archive["keys"].astype(str)
    target_vectors = target_archive["pooled"].astype(np.float32)
    target_vector = {key: vector for key, vector in zip(target_keys, target_vectors)}
    source_protein = np.asarray(
        [target_vector.get(target, np.zeros(target_vectors.shape[1], dtype=np.float32)) for target in source_order]
    )
    source_protein /= np.linalg.norm(source_protein, axis=1, keepdims=True) + 1e-8

    arms = ("b0", "recipient_calibration", "source_support", "source_random", "source_protein", "source_chemistry", "source_gated")
    metrics_by_k: dict[str, object] = {}
    recipient_records: list[dict[str, object]] = []
    for k, episodes in episodes_by_k.items():
        per_arm: dict[str, list[dict[str, float]]] = {arm: [] for arm in arms}
        for episode in episodes:
            support = np.asarray(episode.support, dtype=np.int64)
            query = np.asarray(episode.query, dtype=np.int64)
            support_label = units.set_index("source_row").loc[support, "affinity"].to_numpy(np.float64)
            query_label = units.set_index("source_row").loc[query, "affinity"].to_numpy(np.float64)
            support_base = np.asarray([base_by_row[int(row)] for row in support], dtype=np.float64)
            query_base = np.asarray([base_by_row[int(row)] for row in query], dtype=np.float64)
            calibration_offset = float(np.mean(support_label - support_base))
            calibration = query_base + calibration_offset
            support_calibration = support_base + calibration_offset
            adapter_support = source_adapter[:, 0, None] + source_adapter[:, 1, None] * support_base[None, :]
            support_mse = np.mean(np.square(adapter_support - support_label[None, :]), axis=1)
            best_support = int(np.argmin(support_mse))
            random_source = int(rng.integers(0, len(source_order)))
            recipient_vector = target_vector.get(episode.target, np.zeros(target_vectors.shape[1], dtype=np.float32))
            protein_source = int(np.argmax(_cosine(recipient_vector, source_protein)))
            support_feature = _normalized(features, support, mean, scale).mean(axis=0)
            chemistry_source = int(np.argmax(_cosine(support_feature, source_centroids)))
            query_adapter = source_adapter[:, 0, None] + source_adapter[:, 1, None] * query_base[None, :]
            source_prediction = query_adapter[best_support]
            random_prediction = query_adapter[random_source]
            protein_prediction = query_adapter[protein_source]
            chemistry_prediction = query_adapter[chemistry_source]
            gate = float(support_mse[best_support] <= np.mean(np.square(support_label - support_calibration)))
            gated_prediction = gate * source_prediction + (1.0 - gate) * calibration
            predictions = {
                "b0": query_base,
                "recipient_calibration": calibration,
                "source_support": source_prediction,
                "source_random": random_prediction,
                "source_protein": protein_prediction,
                "source_chemistry": chemistry_prediction,
                "source_gated": gated_prediction,
            }
            scores = {arm: _metrics(query_label, prediction) for arm, prediction in predictions.items()}
            for arm in arms:
                per_arm[arm].append(scores[arm])
            recipient_records.append(
                {
                    "target": episode.target,
                    "endpoint": endpoint,
                    "k": k,
                    "query_units": len(query),
                    "homology": episode.homology,
                    "selected_source_support": source_order[best_support],
                    "selected_source_random": source_order[random_source],
                    "selected_source_protein": source_order[protein_source],
                    "selected_source_chemistry": source_order[chemistry_source],
                    "gate": gate,
                    "source_support_rmse": scores["source_support"]["rmse"],
                    "recipient_calibration_rmse": scores["recipient_calibration"]["rmse"],
                    "source_gated_rmse": scores["source_gated"]["rmse"],
                }
            )
        summary: dict[str, object] = {}
        calibration_scores = per_arm["recipient_calibration"]
        for arm in arms:
            values = per_arm[arm]
            summary[arm] = {
                metric: float(np.nanmean([score[metric] for score in values]))
                for metric in ("rmse", "mae", "spearman", "pairwise_accuracy", "ndcg_at_10")
            }
            gains = np.asarray(
                [cal["rmse"] - score["rmse"] for cal, score in zip(calibration_scores, values)],
                dtype=np.float64,
            )
            summary[arm]["rmse_gain_vs_calibration"] = float(np.nanmean(gains))
            summary[arm]["negative_transfer_rate_vs_calibration"] = float(np.mean(gains < 0.0))
            summary[arm]["benefiting_recipient_rate_vs_calibration"] = float(np.mean(gains > 0.0))
            summary[arm]["rmse_gain_ci95"] = _bootstrap(
                gains, rng
            )
        metrics_by_k[str(k)] = {
            "recipients": len(episodes),
            "summary": summary,
        }

    learning_curve: dict[str, dict[str, float]] = {}
    ks = np.asarray(PRIMARY_K, dtype=np.float64)
    for arm in arms:
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
        "schema_version": "a2s-dta-baseline-v1",
        "protocol": "target-disjoint A2S target-side single-cold source-only ridge plus source adapters",
        "decision": "DESCRIPTIVE_BASELINE_ONLY",
        "endpoint": endpoint,
        "seed": seed,
        "source_rule": f"n_eff >= {SOURCE_MIN}",
        "recipient_rule": f"n_eff < {RECIPIENT_MAX}",
        "primary_support_sizes": list(PRIMARY_K),
        "query_track": "target-side single-cold; source chemical overlap permitted",
        "source_targets": len(source_targets),
        "recipient_candidates": len(recipients),
        "recipient_counts_by_n_eff": {
            "lt5": int(sum(value < 5 for value in counts.values())),
            "lt10": int(sum(value < 10 for value in counts.values())),
            "lt30": int(sum(value < 30 for value in counts.values())),
        },
        "input": {
            "registry_sha256": file_sha256(REGISTRY),
            "features_sha256": file_sha256(FEATURES),
            "target_features_sha256": file_sha256(TARGETS),
            "registry_rows": int(len(features)),
            "feature_row_alignment": "global_registry_source_row",
            "feature_conn_sha_verified": True,
            "split": "train",
            "source_labels_used_for_fit": True,
            "recipient_labels_used_only_support_query": True,
        },
        "ridge": {
            "lambda": ridge,
            "feature_width": int(features.shape[1]),
            "source_weighting": "target_macro_equal_total_weight",
        },
        "metrics_by_k": metrics_by_k,
        "learning_curve": learning_curve,
        "recipient_records": recipient_records,
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
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)
    result = run(
        endpoint=args.endpoint,
        seed=args.seed,
        output=args.out,
        max_recipients=args.max_recipients,
        ridge=args.ridge,
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
