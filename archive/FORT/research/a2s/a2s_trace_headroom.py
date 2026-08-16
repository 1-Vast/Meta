"""A2S-TRACE headroom oracles: an upper bound on learned transport adaptation.

A null result for a learned mechanism is only interpretable next to the ceiling
it failed to reach.  This module measures, on the same probe episodes and the
same frozen base, how far the *oracle* version of each adaptation class can go.
Every oracle here reads the query labels it is scored on, so none of them is a
method; they are upper bounds.

Three classes, in increasing order of what they are allowed to condition on:

``episode_scale``
    one transport scale per episode, chosen with hindsight from a grid.
    Upper bound for any episode-level magnitude router (the TAMSK family).

``episode_subset``
    one support subset per episode, chosen with hindsight.
    Upper bound for any episode-level support-selection rule.

``query_subset``
    a support subset per *query*, chosen with hindsight by coordinate ascent.
    Upper bound for the per-pair reliability class that A2S-TRACE claims.

If ``query_subset`` headroom is large and the learned model captures none of it,
the limit is transfer.  If it is small, the limit is structural and no per-pair
mechanism can pay in this construction.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from research.a2s.a2s_information_gate import canonical
from research.a2s.a2s_trace import (
    ADMITTED_K,
    ADMITTED_STRATA,
    BATCH_EPISODES,
    EVALUATION_EPISODE_SEED,
    TRAIN_POLICIES,
    Substrate,
    analytic_delta,
    build_episodes,
    episode_ci,
    gather_batch,
    group_by_k,
    load_substrate,
    make_batch,
    select_krr_ridge,
    split_fit_episodes,
    tanimoto,
)
from research.a2s.a2s_trace_stratum import (
    DEFAULT_LOCK,
    DEFAULT_OOF,
    MIN_STRATUM_QUERY,
    DrawnEpisode,
    paired_bootstrap,
)


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_trace_headroom_2026-08-01.json"
SCALE_GRID = (0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0)
ASCENT_SWEEPS = 3
ADMITTED_MIN_TANIMOTO = 0.55


def subsets(k: int) -> list[tuple[int, ...]]:
    """All non-empty support subsets, smallest first."""

    return [
        combination
        for size in range(1, k + 1)
        for combination in combinations(range(k), size)
    ]


def episode_scale_oracle(cross: np.ndarray, solved: np.ndarray, base: np.ndarray, label: np.ndarray, scale: float) -> np.ndarray:
    return base + scale * (cross @ solved)


def query_subset_ascent(
    cross: np.ndarray,
    solved: np.ndarray,
    base: np.ndarray,
    label: np.ndarray,
    scale: float,
    options: list[tuple[int, ...]],
) -> tuple[float, np.ndarray]:
    """Coordinate ascent over a per-query support subset, scored by episode CI."""

    n_query, k = cross.shape
    candidate = np.empty((n_query, len(options)), dtype=np.float64)
    for index, option in enumerate(options):
        selector = np.zeros(k, dtype=np.float64)
        selector[list(option)] = 1.0
        candidate[:, index] = scale * (cross * selector) @ solved
    candidate += base[:, None]

    truth = np.sign(label[:, None] - label[None, :])
    active = truth != 0
    total_pairs = float(active.sum()) / 2.0
    if total_pairs <= 0:
        return float("nan"), np.zeros(n_query, dtype=np.int64)

    choice = np.full(n_query, len(options) - 1, dtype=np.int64)  # start from the full support
    prediction = candidate[np.arange(n_query), choice]
    for _ in range(ASCENT_SWEEPS):
        improved = False
        for query in range(n_query):
            # Changing one query's prediction only moves the pairs it belongs to,
            # so the coordinate step is an exact argmax over its own contribution.
            others = np.ones(n_query, dtype=bool)
            others[query] = False
            gap = candidate[query][:, None] - prediction[None, others]
            row_truth = truth[query][others]
            row_active = active[query][others]
            score = ((np.sign(gap) == row_truth) + 0.5 * (np.sign(gap) == 0)) * row_active
            best_index = int(np.argmax(score.sum(axis=1)))
            if best_index != choice[query]:
                improved = True
                choice[query] = best_index
                prediction[query] = candidate[query, best_index]
        if not improved:
            break
    return episode_ci(prediction, label), choice


def measure(substrate: Substrate, episodes: list[DrawnEpisode], ridge: float) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    by_k = group_by_k(episodes, substrate)
    for k, (block, support, query) in sorted(by_k.items()):
        if k not in ADMITTED_K:
            continue
        options = subsets(k)
        for start in range(0, len(block), BATCH_EPISODES):
            rows = np.arange(start, min(start + BATCH_EPISODES, len(block)))
            batch = make_batch(support, query, rows)
            data = gather_batch(substrate, batch)
            cross_t = tanimoto(data["query_bits"], data["support_bits"])
            gram = tanimoto(data["support_bits"], data["support_bits"])
            eye = torch.eye(k, device=gram.device, dtype=gram.dtype)
            solved_t = torch.linalg.solve(
                gram + ridge * eye, data["residual"].unsqueeze(-1)
            ).squeeze(-1)
            cross_all = cross_t.cpu().numpy().astype(np.float64)
            solved_all = solved_t.cpu().numpy().astype(np.float64)
            base_all = data["base"].cpu().numpy().astype(np.float64)
            label_all = data["label"].cpu().numpy().astype(np.float64)
            mask_all = data["mask"].cpu().numpy()
            for offset, row in enumerate(rows):
                episode = block[int(row)]
                nearest = cross_all[offset].max(axis=1)
                active = mask_all[offset] & (nearest >= ADMITTED_MIN_TANIMOTO)
                if int(active.sum()) < MIN_STRATUM_QUERY:
                    continue
                label = label_all[offset][active]
                if float(np.std(label)) < 1e-9:
                    continue
                cross = cross_all[offset][active]
                solved = solved_all[offset]
                base = base_all[offset][active]
                fixed = episode_ci(base + cross @ solved, label)
                if not np.isfinite(fixed):
                    continue
                scale_scores = {
                    scale: episode_ci(episode_scale_oracle(cross, solved, base, label, scale), label)
                    for scale in SCALE_GRID
                }
                best_scale = max(SCALE_GRID, key=lambda value: scale_scores[value])
                subset_best = -np.inf
                for option in options:
                    selector = np.zeros(cross.shape[1], dtype=np.float64)
                    selector[list(option)] = 1.0
                    score = episode_ci(base + (cross * selector) @ solved, label)
                    subset_best = max(subset_best, score)
                query_best, _ = query_subset_ascent(cross, solved, base, label, 1.0, options)
                records.append(
                    {
                        "policy": episode.policy,
                        "k": episode.k,
                        "target": episode.target,
                        "component": episode.component,
                        "n_query": int(active.sum()),
                        "base_ci": episode_ci(base, label),
                        "fixed_krr_ci": fixed,
                        "oracle_episode_scale_ci": float(scale_scores[best_scale]),
                        "oracle_episode_scale_value": float(best_scale),
                        "oracle_episode_subset_ci": float(subset_best),
                        "oracle_query_subset_ci": float(query_best),
                    }
                )
    return pd.DataFrame.from_records(records)


def summarise(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    for policy in sorted(records.policy.unique()):
        for k in sorted(records.k.unique()):
            frame = records.loc[(records.policy == policy) & (records.k == k)].copy()
            if frame.empty:
                continue
            cell: dict[str, object] = {
                "episodes": int(len(frame)),
                "targets": int(frame.target.nunique()),
                "components": int(frame.component.nunique()),
                "absolute_ci": {
                    name: float(frame[f"{name}_ci"].mean())
                    for name in (
                        "base", "fixed_krr", "oracle_episode_scale",
                        "oracle_episode_subset", "oracle_query_subset",
                    )
                },
                "mean_oracle_scale": float(frame.oracle_episode_scale_value.mean()),
                "headroom_over_fixed_krr": {},
            }
            for name in ("oracle_episode_scale", "oracle_episode_subset", "oracle_query_subset"):
                frame["value"] = frame[f"{name}_ci"] - frame["fixed_krr_ci"]
                cell["headroom_over_fixed_krr"][name] = paired_bootstrap(frame, "value")
            summary.setdefault(policy, {})[f"k{int(k)}"] = cell
    return summary


def run(lock_path: Path, output: Path, oof_cache: Path) -> dict[str, object]:
    substrate, context = load_substrate(lock_path, oof_cache)
    fit_all = [e for e in build_episodes(substrate.labeled, "fit") if e.policy in TRAIN_POLICIES]
    _, fit_val = split_fit_episodes(fit_all, substrate)
    ridge, _ = select_krr_ridge(substrate, fit_val)
    probe = [
        episode
        for episode in build_episodes(substrate.labeled, "probe")
        if episode.policy in TRAIN_POLICIES and episode.seed == EVALUATION_EPISODE_SEED
    ]
    records = measure(substrate, probe, ridge)
    summary = summarise(records)
    result = {
        "schema": "a2s-trace-headroom-v1",
        "status": "SOURCE_ONLY_DIAGNOSTIC",
        "note": "Every oracle reads the query labels it is scored on. These are upper bounds, not methods.",
        "lock": {"path": str(lock_path), "content_sha256": context["lock"]["content_sha256"]},
        "protocol": {
            "stratum": list(ADMITTED_STRATA),
            "min_nearest_tanimoto": ADMITTED_MIN_TANIMOTO,
            "k_values": list(ADMITTED_K),
            "krr_ridge_selected_on_fit": ridge,
            "scale_grid": list(SCALE_GRID),
            "ascent_sweeps": ASCENT_SWEEPS,
            "episode_seed": EVALUATION_EPISODE_SEED,
        },
        "summary": summary,
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
    parser.add_argument("--oof-cache", type=Path, default=DEFAULT_OOF)
    args = parser.parse_args()
    result = run(args.lock.resolve(), args.out.resolve(), args.oof_cache.resolve())
    print(json.dumps(result["summary"], indent=2, sort_keys=True, default=float))


if __name__ == "__main__":
    main()
