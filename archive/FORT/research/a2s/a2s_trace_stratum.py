"""A2S-TRACE Q1: stratum-resolved support-information admission gate.

The v2 balanced ChEMBL gate returned ``NO_GO_INFORMATION_NOT_ADMITTED`` under a
single passive support policy that forces the support set and the query set to
be disjoint on scaffold, ligand connectivity, document and assay.  A separate
BindingDB branch reported a large fixed-Tanimoto-KRR gain under random
within-target support.  Those are different estimands, not contradictory
results.

This module measures the boundary directly.  It holds the corpus, the frozen
support-free base, the analytic smoother, the derangement control and the
component bootstrap fixed, and varies exactly two declared axes:

1. the **support policy** (how the k support compounds are drawn), and
2. the **support-query relation stratum** (nearest support Tanimoto of a query).

Nothing here is a proposed adaptation mechanism.  Every estimator is a fixed
analytic operator on frozen-base residuals.  Only the ``fit`` and ``probe``
source roles are opened; the ``locked`` role and the recipient roster are never
requested.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import rankdata

from research.a2s.a2s_information_gate import (
    DEVICE,
    MORGAN_BITS,
    canonical,
    component_oof,
    load_design,
    load_labeled_fit_probe,
    load_metadata,
    oof_folds,
    row_index,
    sha256_file,
    verify_lock,
)
from research.a2s.a2s_information_gate import FEATURES, PROTEINS, REGISTRY


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_LOCK = ROOT / "reports" / "active" / "a2s_source_information_gate_lock_v2_2026-08-01.json"
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_trace_q1_stratum_2026-08-01.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "a2s_trace_q1_records_2026-08-01.parquet"
DEFAULT_OOF = ROOT / "reports" / "active" / "a2s_trace_oof_2026-08-01.npz"

SEED = 20260801
EPISODE_SEEDS = (1729, 1730, 1731)
K_VALUES = (1, 3, 5)
MIN_QUERY = 8
MAX_QUERY = 32
MIN_STRATUM_QUERY = 6
RANDOM_DRAWS = 8
KRR_RIDGE = 0.1
BOOTSTRAP_DRAWS = 2000

POLICIES = ("random_within_target", "scaffold_disjoint", "provenance_disjoint")
DETERMINISTIC_POLICIES = frozenset({"provenance_disjoint"})
STRATUM_EDGES = (0.0, 0.20, 0.35, 0.55, 1.0001)
STRATUM_NAMES = ("t00_20", "t20_35", "t35_55", "t55_100")
ARMS = ("correct", "deranged", "null", "signflip", "wrong_target")


@dataclass(frozen=True)
class DrawnEpisode:
    policy: str
    seed: int
    draw: int
    k: int
    target: str
    component: str
    support_rows: tuple[int, ...]
    query_rows: tuple[int, ...]


def episode_rng(seed: int, policy: str, target: str, k: int, draw: int) -> np.random.Generator:
    digest = sha256(f"{seed}:{policy}:{target}:{k}:{draw}".encode("utf-8")).hexdigest()
    return np.random.default_rng(int(digest[:16], 16))


def _cap(rng: np.random.Generator, values: list[int], limit: int) -> list[int]:
    if len(values) <= limit:
        return sorted(values)
    chosen = rng.choice(np.asarray(values, dtype=np.int64), size=limit, replace=False)
    return sorted(int(value) for value in chosen)


def draw_random_within_target(
    group: pd.DataFrame, k: int, rng: np.random.Generator
) -> tuple[tuple[int, ...], tuple[int, ...]] | None:
    """BindingDB-style policy: support drawn uniformly at random inside a target."""

    rows = group.source_row.to_numpy(dtype=np.int64)
    if len(rows) < k + MIN_QUERY:
        return None
    support = rng.choice(rows, size=k, replace=False)
    support_set = set(int(value) for value in support)
    remainder = [int(value) for value in rows if int(value) not in support_set]
    if len(remainder) < MIN_QUERY:
        return None
    return tuple(sorted(support_set)), tuple(_cap(rng, remainder, MAX_QUERY))


def draw_scaffold_disjoint(
    group: pd.DataFrame, k: int, rng: np.random.Generator
) -> tuple[tuple[int, ...], tuple[int, ...]] | None:
    """Random support, but every query is scaffold- and connectivity-cold to it."""

    rows = group.source_row.to_numpy(dtype=np.int64)
    if len(rows) < k + MIN_QUERY:
        return None
    scaffolds = group.set_index("source_row").scaffold.astype(str).to_dict()
    conns = group.set_index("source_row").conn.astype(str).to_dict()
    order = rng.permutation(rows)
    support: list[int] = []
    used_scaffold: set[str] = set()
    for value in order:
        row = int(value)
        if scaffolds[row] in used_scaffold:
            continue
        support.append(row)
        used_scaffold.add(scaffolds[row])
        if len(support) == k:
            break
    if len(support) < k:
        return None
    support_set = set(support)
    used_conn = {conns[row] for row in support}
    remainder = [
        int(value)
        for value in rows
        if int(value) not in support_set
        and scaffolds[int(value)] not in used_scaffold
        and conns[int(value)] not in used_conn
    ]
    if len(remainder) < MIN_QUERY:
        return None
    return tuple(sorted(support_set)), tuple(_cap(rng, remainder, MAX_QUERY))


def draw_provenance_disjoint(
    group: pd.DataFrame, k: int, rng: np.random.Generator
) -> tuple[tuple[int, ...], tuple[int, ...]] | None:
    """The frozen v2 policy: support and query disjoint on four provenance axes."""

    ordered = group.sort_values(["scaffold", "conn", "docs", "assays", "source_row"], kind="stable")
    fields = ("conn", "scaffold", "docs", "assays")
    used = {field: set() for field in fields}
    support: list[int] = []
    for row in ordered.itertuples(index=False):
        if any(str(getattr(row, field)) in used[field] for field in fields):
            continue
        support.append(int(row.source_row))
        for field in fields:
            used[field].add(str(getattr(row, field)))
        if len(support) == k:
            break
    if len(support) < k:
        return None
    support_set = set(support)
    query: list[int] = []
    for row in ordered.itertuples(index=False):
        source_row = int(row.source_row)
        if source_row in support_set:
            continue
        if any(str(getattr(row, field)) in used[field] for field in fields):
            continue
        query.append(source_row)
        if len(query) == MAX_QUERY:
            break
    if len(query) < MIN_QUERY:
        return None
    return tuple(sorted(support_set)), tuple(sorted(query))


POLICY_FUNCTIONS = {
    "random_within_target": draw_random_within_target,
    "scaffold_disjoint": draw_scaffold_disjoint,
    "provenance_disjoint": draw_provenance_disjoint,
}


def build_episodes(labeled: pd.DataFrame, role: str) -> list[DrawnEpisode]:
    episodes: list[DrawnEpisode] = []
    active = labeled.loc[labeled.role == role]
    for target, group in active.groupby("target", sort=True):
        component = str(group.component.iloc[0])
        for policy in POLICIES:
            draws = 1 if policy in DETERMINISTIC_POLICIES else RANDOM_DRAWS
            seeds = (EPISODE_SEEDS[0],) if policy in DETERMINISTIC_POLICIES else EPISODE_SEEDS
            for seed in seeds:
                for draw in range(draws):
                    for k in K_VALUES:
                        rng = episode_rng(seed, policy, str(target), k, draw)
                        selected = POLICY_FUNCTIONS[policy](group, k, rng)
                        if selected is None:
                            continue
                        support, query = selected
                        episodes.append(
                            DrawnEpisode(
                                policy=policy,
                                seed=int(seed),
                                draw=int(draw),
                                k=int(k),
                                target=str(target),
                                component=component,
                                support_rows=support,
                                query_rows=query,
                            )
                        )
    return episodes


def tanimoto_matrix(query_bits: np.ndarray, support_bits: np.ndarray) -> np.ndarray:
    intersection = query_bits @ support_bits.T
    query_count = query_bits.sum(axis=1, keepdims=True)
    support_count = support_bits.sum(axis=1, keepdims=True).T
    union = query_count + support_count - intersection
    return (intersection / np.maximum(union, 1.0)).astype(np.float64)


def stratum_of(values: np.ndarray) -> np.ndarray:
    index = np.digitize(values, np.asarray(STRATUM_EDGES[1:-1], dtype=np.float64), right=False)
    return np.asarray([STRATUM_NAMES[int(position)] for position in index], dtype=object)


def metric_loss(label: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    """Identical semantics to ``a2s_information_gate.metric_loss``, vectorised.

    The gate version calls ``pandas.Series.rank`` per metric, which dominates the
    runtime of a stratified sweep.  ``scipy.stats.rankdata`` gives the same
    average-tie ranks at a fraction of the cost, and the pairwise CI, NDCG@10 and
    RMSE definitions are copied unchanged.
    """

    label = np.asarray(label, dtype=np.float64)
    prediction = np.asarray(prediction, dtype=np.float64)
    if len(label) < 2:
        return {"rmse": float("nan"), "ci": float("nan"), "spearman": float("nan"), "ndcg10": float("nan"), "rank_loss": float("nan")}
    error = label - prediction
    left, right = np.triu_indices(len(label), k=1)
    truth = np.sign(label[left] - label[right])
    active = truth != 0
    pred = np.sign(prediction[left] - prediction[right])
    if active.any():
        ci = float((pred[active] == truth[active]).mean() + 0.5 * (pred[active] == 0).mean())
    else:
        ci = float("nan")
    label_rank = rankdata(label)
    prediction_rank = rankdata(prediction)
    if np.std(label_rank) and np.std(prediction_rank):
        spearman = float(np.corrcoef(label_rank, prediction_rank)[0, 1])
    else:
        spearman = float("nan")
    cutoff = min(10, len(label))
    order = np.argsort(-prediction, kind="stable")[:cutoff]
    ideal = np.argsort(-label, kind="stable")[:cutoff]
    gain = label - label.min() + 1e-6
    discounts = 1.0 / np.log2(np.arange(2, cutoff + 2))
    dcg = float(np.sum(gain[order] * discounts))
    ideal_dcg = float(np.sum(gain[ideal] * discounts))
    ndcg = dcg / ideal_dcg if ideal_dcg > 0 else float("nan")
    return {"rmse": float(np.sqrt(np.mean(error**2))), "ci": ci, "spearman": spearman, "ndcg10": ndcg, "rank_loss": 1.0 - ndcg}


def derangement(size: int, rng: np.random.Generator) -> np.ndarray:
    if size < 2:
        raise ValueError("a derangement needs at least two support compounds")
    while True:
        candidate = rng.permutation(size)
        if not np.any(candidate == np.arange(size)):
            return candidate


def krr_weights(support_bits: np.ndarray, query_bits: np.ndarray, ridge: float) -> np.ndarray:
    gram = tanimoto_matrix(support_bits, support_bits)
    cross = tanimoto_matrix(query_bits, support_bits)
    eye = np.eye(gram.shape[0], dtype=np.float64)
    return cross @ np.linalg.inv(gram + ridge * eye)


def nw_weights(cross: np.ndarray) -> np.ndarray:
    return cross / np.maximum(cross.sum(axis=1, keepdims=True), 1e-9)


def evaluate_episodes(
    episodes: list[DrawnEpisode],
    labeled: pd.DataFrame,
    raw_features: np.ndarray,
    base: np.ndarray,
) -> pd.DataFrame:
    position = row_index(labeled)
    affinity = labeled.affinity.to_numpy(dtype=np.float64)
    residual_all = affinity - base.astype(np.float64)
    bits = raw_features[:, :MORGAN_BITS].astype(np.float64)

    # Wrong-target support residuals are drawn from a different probe target so
    # that the residual multiset is real but the compound-evidence assignment is
    # biologically meaningless.
    by_target: dict[str, np.ndarray] = {}
    for target, group in labeled.groupby("target", sort=True):
        by_target[str(target)] = group.source_row.to_numpy(dtype=np.int64)
    target_names = sorted(by_target)

    records: list[dict[str, object]] = []
    for episode in episodes:
        support = np.asarray(episode.support_rows, dtype=np.int64)
        query = np.asarray(episode.query_rows, dtype=np.int64)
        support_index = np.asarray([position[int(row)] for row in support], dtype=np.int64)
        query_index = np.asarray([position[int(row)] for row in query], dtype=np.int64)
        support_bits = bits[support]
        query_bits = bits[query]
        cross = tanimoto_matrix(query_bits, support_bits)
        nearest = cross.max(axis=1)
        strata = stratum_of(nearest)

        residual = residual_all[support_index]
        base_query = base.astype(np.float64)[query_index]
        label_query = affinity[query_index]

        rng = episode_rng(episode.seed + 7, "arms:" + episode.policy, episode.target, episode.k, episode.draw)
        arm_residual: dict[str, np.ndarray] = {
            "correct": residual,
            "null": np.zeros_like(residual),
            "signflip": -residual,
        }
        if episode.k >= 3:
            arm_residual["deranged"] = residual[derangement(episode.k, rng)]
        other = [name for name in target_names if name != episode.target]
        donor = str(other[int(rng.integers(0, len(other)))])
        donor_rows = by_target[donor]
        donor_pick = rng.choice(donor_rows, size=min(episode.k, len(donor_rows)), replace=False)
        donor_index = np.asarray([position[int(row)] for row in donor_pick], dtype=np.int64)
        donor_residual = residual_all[donor_index]
        if len(donor_residual) < episode.k:
            donor_residual = np.resize(donor_residual, episode.k)
        # Norm-match the wrong-target residual so magnitude cannot explain a gap.
        donor_norm = float(np.linalg.norm(donor_residual))
        if donor_norm > 1e-9:
            donor_residual = donor_residual * (float(np.linalg.norm(residual)) / donor_norm)
        arm_residual["wrong_target"] = donor_residual

        weights = {
            "nw": nw_weights(cross),
            "krr": krr_weights(support_bits, query_bits, KRR_RIDGE),
            "level": np.full((len(query), episode.k), 1.0 / episode.k, dtype=np.float64),
        }

        predictions: dict[str, np.ndarray] = {"base__correct": base_query}
        for estimator, weight in weights.items():
            for arm, values in arm_residual.items():
                predictions[f"{estimator}__{arm}"] = base_query + weight @ values

        for stratum in (*STRATUM_NAMES, "all"):
            mask = np.ones(len(query), dtype=bool) if stratum == "all" else (strata == stratum)
            if int(mask.sum()) < MIN_STRATUM_QUERY:
                continue
            truth = label_query[mask]
            if float(np.std(truth)) < 1e-9:
                continue
            row: dict[str, object] = {
                "policy": episode.policy,
                "seed": episode.seed,
                "draw": episode.draw,
                "k": episode.k,
                "target": episode.target,
                "component": episode.component,
                "stratum": stratum,
                "n_query": int(mask.sum()),
                "nearest_tanimoto_mean": float(nearest[mask].mean()),
                "support_residual_sd": float(np.std(residual)),
            }
            for name, values in predictions.items():
                metrics = metric_loss(truth, values[mask])
                for metric, value in metrics.items():
                    row[f"{name}__{metric}"] = float(value)
            records.append(row)
    return pd.DataFrame.from_records(records)


def paired_bootstrap(frame: pd.DataFrame, column: str, *, seed: int = SEED, draws: int = BOOTSTRAP_DRAWS) -> dict[str, float]:
    """Aggregate draws -> seed/target mean -> component mean -> paired bootstrap."""

    usable = frame[["component", "target", column]].dropna()
    if usable.empty:
        return {"components": 0, "mean": float("nan"), "lower95": float("nan"), "upper95": float("nan")}
    per_target = usable.groupby(["component", "target"], sort=True)[column].mean().reset_index()
    per_component = per_target.groupby("component", sort=True)[column].mean().to_numpy(dtype=np.float64)
    per_component = per_component[np.isfinite(per_component)]
    if len(per_component) == 0:
        return {"components": 0, "mean": float("nan"), "lower95": float("nan"), "upper95": float("nan")}
    rng = np.random.default_rng(seed)
    sample = rng.integers(0, len(per_component), size=(draws, len(per_component)))
    means = per_component[sample].mean(axis=1)
    return {
        "components": int(len(per_component)),
        "targets": int(len(per_target)),
        "mean": float(per_component.mean()),
        "lower95": float(np.quantile(means, 0.025)),
        "upper95": float(np.quantile(means, 0.975)),
    }


CONTRASTS = {
    "krr_minus_base": ("krr__correct", "base__correct"),
    "krr_minus_level": ("krr__correct", "level__correct"),
    "krr_correct_minus_deranged": ("krr__correct", "krr__deranged"),
    "krr_correct_minus_wrong_target": ("krr__correct", "krr__wrong_target"),
    "krr_correct_minus_signflip": ("krr__correct", "krr__signflip"),
    "nw_minus_base": ("nw__correct", "base__correct"),
    "nw_correct_minus_deranged": ("nw__correct", "nw__deranged"),
    "level_minus_base": ("level__correct", "base__correct"),
}
RANK_METRICS = ("ci", "ndcg10", "spearman")
ALL_METRICS = ("ci", "ndcg10", "spearman", "rmse")


def summarise(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    for policy in sorted(records.policy.unique()):
        policy_frame = records.loc[records.policy == policy]
        policy_summary: dict[str, object] = {}
        for k in sorted(policy_frame.k.unique()):
            k_frame = policy_frame.loc[policy_frame.k == k]
            k_summary: dict[str, object] = {}
            for stratum in sorted(k_frame.stratum.unique()):
                stratum_frame = k_frame.loc[k_frame.stratum == stratum].copy()
                entry: dict[str, object] = {
                    "episodes": int(len(stratum_frame)),
                    "targets": int(stratum_frame.target.nunique()),
                    "components": int(stratum_frame.component.nunique()),
                    "mean_query_per_episode": float(stratum_frame.n_query.mean()),
                    "nearest_tanimoto_mean": float(stratum_frame.nearest_tanimoto_mean.mean()),
                    "absolute": {},
                    "contrasts": {},
                }
                for estimator in ("base__correct", "level__correct", "nw__correct", "krr__correct"):
                    entry["absolute"][estimator] = {
                        metric: float(np.nanmean(stratum_frame[f"{estimator}__{metric}"].to_numpy()))
                        for metric in ALL_METRICS
                        if f"{estimator}__{metric}" in stratum_frame
                    }
                for name, (left, right) in CONTRASTS.items():
                    for metric in ALL_METRICS:
                        left_column = f"{left}__{metric}"
                        right_column = f"{right}__{metric}"
                        if left_column not in stratum_frame or right_column not in stratum_frame:
                            continue
                        sign = -1.0 if metric == "rmse" else 1.0
                        stratum_frame["_delta"] = sign * (
                            stratum_frame[left_column] - stratum_frame[right_column]
                        )
                        entry["contrasts"].setdefault(name, {})[metric] = paired_bootstrap(
                            stratum_frame.rename(columns={"_delta": "value"}), "value"
                        )
                k_summary[stratum] = entry
            policy_summary[f"k{int(k)}"] = k_summary
        summary[policy] = policy_summary
    return summary


ADMISSION_MDE = 0.005


def admission(summary: dict[str, object]) -> dict[str, object]:
    """Preregistered stratum admission rule (see the TRACE exploration prompt)."""

    admitted: list[dict[str, object]] = []
    for policy, policy_summary in summary.items():
        for k_label, k_summary in policy_summary.items():
            k = int(k_label[1:])
            for stratum, entry in k_summary.items():
                contrasts = entry.get("contrasts", {})
                gain = contrasts.get("krr_minus_base", {}).get("ci")
                assign = contrasts.get("krr_correct_minus_deranged", {}).get("ci")
                if gain is None:
                    continue
                gain_pass = bool(gain["lower95"] > ADMISSION_MDE)
                assign_pass = bool(assign is not None and assign["lower95"] > 0.0) if k >= 3 else None
                record = {
                    "policy": policy,
                    "k": k,
                    "stratum": stratum,
                    "components": gain["components"],
                    "ci_gain_lower95": gain["lower95"],
                    "ci_gain_mean": gain["mean"],
                    "assignment_lower95": None if assign is None else assign["lower95"],
                    "gain_pass": gain_pass,
                    "assignment_pass": assign_pass,
                    "admitted": bool(gain_pass and (assign_pass is not False)),
                }
                admitted.append(record)
    ranked = sorted(admitted, key=lambda item: (-float(item["ci_gain_lower95"]), item["policy"], item["k"]))
    return {
        "rule": {
            "ci_gain_lower95_above": ADMISSION_MDE,
            "assignment_lower95_above": 0.0,
            "assignment_applies_for_k": [3, 5],
        },
        "records": ranked,
        "any_admitted": any(bool(item["admitted"]) for item in ranked),
    }


def resolve_base(labeled: pd.DataFrame, values: np.ndarray, metadata: pd.DataFrame, cache: Path) -> tuple[np.ndarray, dict[str, object]]:
    if cache.exists():
        stored = np.load(cache, allow_pickle=False)
        if np.array_equal(stored["source_row"], labeled.source_row.to_numpy(np.int64)):
            base = np.asarray(stored["base"], dtype=np.float32)
            if np.isfinite(base).all():
                return base, {"reused": True, "path": str(cache)}
    base, stats = component_oof(metadata, labeled, values)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache,
        source_row=labeled.source_row.to_numpy(np.int64),
        base=base,
        residual=labeled.affinity.to_numpy(np.float32) - base,
    )
    stats["reused"] = False
    stats["path"] = str(cache)
    return base, stats


def run(lock_path: Path, output: Path, records_path: Path, oof_cache: Path, *, role: str = "probe") -> dict[str, object]:
    if DEVICE.type != "cuda":
        raise RuntimeError("run this gate with D:\\anaconda\\envs\\drug\\python.exe")
    lock = verify_lock(lock_path)
    metadata = load_metadata(lock)
    labeled = load_labeled_fit_probe(metadata)
    values, raw_features, _ = load_design(metadata, labeled)
    base, oof_stats = resolve_base(labeled, values, metadata, oof_cache)

    episodes = build_episodes(labeled, role)
    if not episodes:
        raise RuntimeError("no episodes could be drawn for the requested role")
    records = evaluate_episodes(episodes, labeled, raw_features, base)
    if records.empty:
        raise RuntimeError("no evaluable episode/stratum cells")
    records_path.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)

    summary = summarise(records)
    verdict = admission(summary)
    result: dict[str, object] = {
        "schema": "a2s-trace-q1-stratum-v1",
        "status": "SOURCE_ONLY_DIAGNOSTIC",
        "question": "Q1: under which support policy and support-query relation stratum does correctly assigned support carry transferable ranking information?",
        "role_measured": role,
        "labels": {
            "opened_roles": ["fit", "probe"],
            "locked_labels_requested": False,
            "recipient_labels_requested": False,
        },
        "device": torch.cuda.get_device_name(0),
        "lock": {"path": str(lock_path), "content_sha256": lock["content_sha256"]},
        "inputs": {
            "registry_sha256": sha256_file(REGISTRY),
            "features_sha256": sha256_file(FEATURES),
            "proteins_sha256": sha256_file(PROTEINS),
        },
        "protocol": {
            "policies": list(POLICIES),
            "episode_seeds": list(EPISODE_SEEDS),
            "random_draws_per_target": RANDOM_DRAWS,
            "k_values": list(K_VALUES),
            "query_cap": MAX_QUERY,
            "min_query": MIN_QUERY,
            "min_stratum_query": MIN_STRATUM_QUERY,
            "stratum_edges": list(STRATUM_EDGES),
            "krr_ridge": KRR_RIDGE,
            "bootstrap_draws": BOOTSTRAP_DRAWS,
            "aggregation": "episode draws -> seed/target mean -> component mean -> paired component bootstrap",
            "arms": list(ARMS),
        },
        "oof": oof_stats,
        "episodes": {
            "total": len(episodes),
            "by_policy": {
                policy: int(sum(1 for episode in episodes if episode.policy == policy))
                for policy in POLICIES
            },
            "evaluable_cells": int(len(records)),
        },
        "summary": summary,
        "admission": verdict,
        "interpretation": {
            "fact": "Every estimator here is a fixed analytic operator on frozen-base residuals; no adaptation parameter was fitted.",
            "inference": "A stratum where the KRR gain lower bound clears the MDE while the derangement contrast is positive is one where compound-to-residual assignment carries transferable ranking information.",
            "hypothesis": "A learned transport rule can beat the fixed analytic smoother only inside an admitted stratum.",
        },
    }
    result["content_sha256"] = sha256(canonical(result).encode("utf-8")).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=True, default=float) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    parser.add_argument("--oof-cache", type=Path, default=DEFAULT_OOF)
    parser.add_argument("--role", type=str, default="probe", choices=("fit", "probe"))
    args = parser.parse_args()
    result = run(
        args.lock.resolve(),
        args.out.resolve(),
        args.records.resolve(),
        args.oof_cache.resolve(),
        role=args.role,
    )
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.out.resolve()),
                "episodes": result["episodes"],
                "any_admitted": result["admission"]["any_admitted"],
                "top_strata": result["admission"]["records"][:8],
            },
            indent=2,
            sort_keys=True,
            default=float,
        )
    )


if __name__ == "__main__":
    main()
