"""Offline compiled-dataset fixtures used to verify input isolation."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch

from scripts.data_contract import read_jsonl
from scripts.data import EpisodeTarget, ModelEpisode


def build_episodes(rows_path: str | Path, output_path: str | Path, *, split: str = "train",
                   k: int = 5, draws_per_task: int = 8, seed: int = 0) -> dict:
    """Compile disjoint support/query indices without serializing query labels."""
    if k < 1 or draws_per_task < 1:
        raise ValueError("k and draws_per_task must be positive")
    rows = read_jsonl(rows_path)
    by_task: dict[str, list[int]] = {}
    for index, row in enumerate(rows):
        if row["split"] == split:
            by_task.setdefault(row["task_key"], []).append(index)
    rng = np.random.default_rng(seed)
    query, support, task = [], [], []
    skipped = {"too_few_rows": 0}
    task_keys = sorted(by_task)
    for task_index, task_key in enumerate(task_keys):
        indices = np.asarray(by_task[task_key], dtype=np.int64)
        if len(indices) < k + 1:
            skipped["too_few_rows"] += 1
            continue
        for _ in range(draws_per_task):
            selected = rng.choice(indices, size=k + 1, replace=False)
            query.append(selected[0])
            support.append(selected[1:])
            task.append(task_index)
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(target, q_row_idx=np.asarray(query, dtype=np.int64),
                        s_row_idx=np.asarray(support, dtype=np.int64),
                        task_idx=np.asarray(task, dtype=np.int64),
                        task_keys=np.asarray(task_keys))
    return {"split": split, "k": k, "episodes": len(query), "tasks": len(task_keys),
            "skipped": skipped, "contains_query_label": False}


def build_index_maps(rows_path: str | Path) -> tuple[list[dict], dict[str, int], dict[str, int]]:
    rows = read_jsonl(rows_path)
    drugs = {key: index for index, key in enumerate(sorted({row["drug_key"] for row in rows}))}
    targets = {key: index for index, key in enumerate(sorted({row["target_key"] for row in rows}))}
    return rows, drugs, targets


def batch_from_episode_indices(rows_path: str | Path, episode_path: str | Path, *,
                               device=None, dtype=torch.float64) -> tuple[ModelEpisode, torch.Tensor]:
    """Return model inputs and query targets separately for leakage checks."""
    rows, drugs, targets = build_index_maps(rows_path)
    episode = np.load(episode_path)
    query_rows = episode["q_row_idx"]
    support_rows = episode["s_row_idx"]
    query = [rows[int(index)] for index in query_rows]
    support = [[rows[int(index)] for index in indices] for indices in support_rows]
    return ModelEpisode(
        target_idx=torch.as_tensor([targets[row["target_key"]] for row in query], dtype=torch.long, device=device),
        query_pair_idx=torch.as_tensor([drugs[row["drug_key"]] for row in query], dtype=torch.long, device=device),
        support_pair_idx=torch.as_tensor([[drugs[row["drug_key"]] for row in group] for group in support], dtype=torch.long, device=device),
        support_y=torch.as_tensor([[row["y"] for row in group] for group in support], dtype=dtype, device=device),
        support_mask=torch.ones((len(query), support_rows.shape[1]), dtype=dtype, device=device),
        context_id=torch.as_tensor([row["context_id"] for row in query], dtype=torch.long, device=device),
        context_cont=torch.as_tensor([row["context_cont"] for row in query], dtype=dtype, device=device),
        context_mask=torch.as_tensor([row["context_mask"] for row in query], dtype=dtype, device=device),
    ), torch.as_tensor([row["y"] for row in query], dtype=dtype, device=device)


def _band(labels: list[float], grid: np.ndarray, minimum: int) -> tuple[np.ndarray, np.ndarray, str]:
    if len(labels) < minimum:
        return np.zeros_like(grid), np.ones_like(grid), "vacuous"
    values = np.asarray(labels, dtype=float)
    cdf = np.asarray([(values <= point).mean() for point in grid])
    radius = float(np.sqrt(np.log(20.0) / (2.0 * len(values))))
    return np.clip(cdf - radius, 0.0, 1.0), np.clip(cdf + radius, 0.0, 1.0), "empirical"


def freeze_population_bands(rows_path: str | Path, output_path: str | Path, *, grid_size: int = 33,
                            minimum_count: int = 30) -> dict:
    """Freeze source-only empirical bands for offline regression checks."""
    rows = [row for row in read_jsonl(rows_path) if row["split"] == "train"]
    grid = np.linspace(0.0, 1.0, grid_size)
    by_context: dict[str, list[float]] = {}
    for row in rows:
        by_context.setdefault(row["endpoint_key"], []).append(float(row["y"]))
    lower, upper, status = {}, {}, {}
    for context_key, labels in sorted(by_context.items()):
        lo, hi, state = _band(labels, grid, minimum_count)
        lower[context_key], upper[context_key], status[context_key] = lo.tolist(), hi.tolist(), state
    result = {"schema": "MetaSieve.PopulationBands.v1", "grid": grid.tolist(),
              "minimum_count": minimum_count, "context_bands":
              {key: {"lower": lower[key], "upper": upper[key], "status": status[key],
                     "source_rows": len(by_context[key])} for key in lower},
              "fallback": "vacuous [0,1] within the requested endpoint only"}
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result
