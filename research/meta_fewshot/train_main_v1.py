"""Train and evaluate experimental MetaSieve V1 on cold-target few-shot DTA.

V1 keeps the admitted main-v0 representation and uncentered positive ridge.
It tests two training-only interventions inspired by AdaMBind: a source-task
scheduler and support-label perturbation.  Neither is loaded for inference.
"""
from __future__ import annotations

import argparse
import copy
from collections import defaultdict
from dataclasses import asdict, dataclass, replace
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
import torch
import torch.nn.functional as F
from scipy import stats

from model.metasieve_v1 import MetaSieveV1, TaskScheduler, uniform_label_noise
from model.runtime import require_cuda
from research.meta_fewshot.train_main_v0 import (
    CORPUS,
    FEATURES,
    concordance,
    json_safe,
    load_data,
    sha256,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "report/meta_fewshot/main_v1"
SEEDS = (20260831, 20260832, 20260833)
VARIANTS = {
    "uniform_clean": ("uniform", False),
    "uniform_support_noise": ("uniform", True),
    "ats_clean": ("ats", False),
    "ats_support_noise": ("ats", True),
    "ats_null_clean": ("ats_null", False),
    "ats_null_support_noise": ("ats_null", True),
}
@dataclass(frozen=True)
class V1TrainConfig:
    support_sizes: tuple[int, ...] = (1, 2, 3, 5)
    min_query: int = 3
    section_dim: int = 2
    ridge: float = 1.0
    support_only_section: bool = False
    population_hidden_dim: int = 0
    pair_hidden_dim: int = 0
    population_pretrain_steps: int = 0
    population_pretrain_batch_size: int = 1024
    population_pretrain_learning_rate: float = 0.001
    candidate_tasks: int = 15
    selected_tasks: int = 8
    query_per_task: int = 32
    steps: int = 1000
    learning_rate: float = 0.01
    weight_decay: float = 1e-4
    scheduler_hidden_dim: int = 8
    scheduler_learning_rate: float = 0.003
    scheduler_target_temperature: float = 0.25
    support_noise_half_width_pk: float = 0.2
    val_interval: int = 100
    val_draws: int = 2
    test_draws: int = 5
    test_max_query: int = 128
    bootstrap_draws: int = 9999


def build_tasks(cells: list[dict], split: str, k: int,
                min_query: int) -> dict[str, np.ndarray]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(cells):
        if row["split"] == split:
            grouped[row["target_id"]].append(index)
    return {target: np.asarray(indices, dtype=np.int64)
            for target, indices in sorted(grouped.items())
            if len(indices) >= k + min_query}


def cluster_tasks(cells: list[dict], tasks: dict[str, np.ndarray]) -> dict[str, list[str]]:
    target_cluster = {row["target_id"]: row["protein_group_40"] for row in cells}
    grouped: dict[str, list[str]] = defaultdict(list)
    for target in tasks:
        grouped[target_cluster[target]].append(target)
    return {cluster: sorted(targets) for cluster, targets in sorted(grouped.items())}


def draw_episode(indices: np.ndarray, rng: np.random.Generator, k: int,
                 max_query: int | None) -> tuple[np.ndarray, np.ndarray]:
    order = rng.permutation(indices)
    support, query = order[:k], order[k:]
    if max_query is not None and len(query) > max_query:
        query = rng.choice(query, size=max_query, replace=False)
    return np.asarray(support), np.asarray(query)


def draw_nested_episode(indices: np.ndarray, rng: np.random.Generator, k: int,
                        max_k: int, max_query: int | None) -> tuple[np.ndarray, np.ndarray]:
    """Draw nested support prefixes and a k-invariant query from one permutation."""
    if not 1 <= k <= max_k or len(indices) <= max_k:
        raise ValueError("invalid nested few-shot episode")
    order = rng.permutation(indices)
    support, query = order[:k], order[max_k:]
    if max_query is not None and len(query) > max_query:
        query = rng.choice(query, size=max_query, replace=False)
    return np.asarray(support), np.asarray(query)


def sample_cluster_balanced(clusters: dict[str, list[str]], rng: np.random.Generator,
                            count: int) -> list[str]:
    keys = sorted(clusters)
    if count > len(keys):
        raise ValueError("candidate task count exceeds eligible protein clusters")
    selected_clusters = rng.choice(len(keys), size=count, replace=False)
    return [clusters[keys[int(index)]][
        int(rng.integers(0, len(clusters[keys[int(index)]])))]
        for index in selected_clusters]


def _flat_gradient(loss: torch.Tensor, parameters: tuple[torch.nn.Parameter, ...], *,
                   retain_graph: bool) -> torch.Tensor:
    gradients = torch.autograd.grad(
        loss, parameters, retain_graph=retain_graph, allow_unused=True)
    return torch.cat([
        (torch.zeros_like(parameter) if gradient is None else gradient).reshape(-1)
        for parameter, gradient in zip(parameters, gradients)
    ])


def _flat_batched_gradients(
        losses: torch.Tensor, parameters: tuple[torch.nn.Parameter, ...], *,
        retain_graph: bool) -> torch.Tensor:
    """Return one flattened parameter gradient per scalar loss."""
    identity = torch.eye(
        len(losses), device=losses.device, dtype=losses.dtype)
    gradients = torch.autograd.grad(
        losses, parameters, grad_outputs=identity, retain_graph=retain_graph,
        allow_unused=True, is_grads_batched=True)
    return torch.cat([
        ((torch.zeros(
            (len(losses),) + tuple(parameter.shape),
            device=parameter.device, dtype=parameter.dtype)
          if gradient is None else gradient).flatten(start_dim=1))
        for parameter, gradient in zip(parameters, gradients)
    ], dim=1)


def gradient_cosine(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
    denominator = torch.linalg.vector_norm(left) * torch.linalg.vector_norm(right)
    safe_denominator = torch.where(
        denominator > 0, denominator, torch.ones_like(denominator))
    return torch.where(
        denominator > 0, torch.dot(left, right) / safe_denominator,
        torch.zeros_like(denominator))


def gradient_cosine_rows(left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
    if left.ndim != 2 or left.shape != right.shape:
        raise ValueError("batched gradients need the same [tasks,parameters] shape")
    denominator = (torch.linalg.vector_norm(left, dim=1)
                   * torch.linalg.vector_norm(right, dim=1))
    numerator = (left * right).sum(dim=1)
    safe_denominator = torch.where(
        denominator > 0, denominator, torch.ones_like(denominator))
    return torch.where(
        denominator > 0, numerator / safe_denominator,
        torch.zeros_like(numerator))


def pack_episode_indices(
        episodes: list[tuple[np.ndarray, np.ndarray]],
        device: torch.device | str) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Pack fixed-k, variable-query episode indices for one GPU batch."""
    if not episodes or any(len(query) < 1 for _, query in episodes):
        raise ValueError("cannot pack empty episodes or queries")
    support_size = len(episodes[0][0])
    if any(len(support) != support_size for support, _ in episodes):
        raise ValueError("batched episodes require one shared support size")
    support = np.stack([item[0] for item in episodes]).astype(np.int64, copy=False)
    max_query = max(len(item[1]) for item in episodes)
    query = np.empty((len(episodes), max_query), dtype=np.int64)
    mask = np.zeros((len(episodes), max_query), dtype=np.bool_)
    for index, (_, values) in enumerate(episodes):
        query[index] = values[0]
        query[index, :len(values)] = values
        mask[index, :len(values)] = True
    return (
        torch.as_tensor(support, device=device),
        torch.as_tensor(query, device=device),
        torch.as_tensor(mask, device=device),
    )


def add_episodewise_label_noise(
        support_y: torch.Tensor, standard_deviation: float, *,
        generator: torch.Generator) -> torch.Tensor:
    """Preserve the scalar-episode RNG stream while batching model compute."""
    if support_y.ndim != 2:
        raise ValueError("episodewise support labels need shape [tasks,k]")
    return torch.stack([
        uniform_label_noise(row, standard_deviation, generator=generator)
        for row in support_y
    ])


def _episode_loss(model: MetaSieveV1, tensors: dict, family: str,
                  support: np.ndarray, query: np.ndarray,
                  support_y: torch.Tensor) -> torch.Tensor:
    prediction = model.episode(
        tensors["ligand"][support], tensors[family][support], support_y,
        tensors["ligand"][query], tensors[family][query])
    return (prediction - tensors["y"][query]).square().mean()


def validation_metrics(model: MetaSieveV1, cells: list[dict], tensors: dict,
                       family: str, split: str, config: V1TrainConfig,
                       seed: int, y_scale: float, *,
                       noisy_support: bool = False) -> dict[int, float]:
    model.eval()
    result = {}
    max_k = max(config.support_sizes)
    tasks = build_tasks(cells, split, max_k, config.min_query)
    with torch.no_grad():
        for k in config.support_sizes:
            episodes = []
            identities = []
            for target, indices in tasks.items():
                for draw in range(config.val_draws):
                    episode_seed = _episode_seed(seed, target, draw)
                    support, query = draw_nested_episode(
                        indices, np.random.default_rng(episode_seed), k, max_k,
                        config.query_per_task)
                    episodes.append((support, query))
                    identities.append((target, draw))
            support_index, query_index, query_mask = pack_episode_indices(
                episodes, tensors["y"].device)
            support_y = tensors["y"][support_index]
            if noisy_support:
                noise_sd = (
                    config.support_noise_half_width_pk / math.sqrt(3.0) / y_scale)
                support_y = torch.stack([
                    uniform_label_noise(
                        row, noise_sd,
                        generator=torch.Generator(
                            device=tensors["y"].device).manual_seed(
                                _episode_seed(seed + 2_017_000, target, draw)))
                    for row, (target, draw) in zip(support_y, identities)
                ])
            prediction, _ = model.batched_episode(
                tensors["ligand"][support_index], tensors[family][support_index],
                support_y, tensors["ligand"][query_index],
                tensors[family][query_index])
            squared = ((prediction - tensors["y"][query_index]) * y_scale).square()
            episode_loss = (
                (squared * query_mask).sum(dim=1)
                / query_mask.sum(dim=1).to(squared.dtype))
            target_loss = episode_loss.reshape(len(tasks), config.val_draws).mean(dim=1)
            result[k] = float(target_loss.mean())
    return result


def _mean_validation_score(clean: dict[int, float], noisy: dict[int, float]) -> float:
    return float(np.mean([clean[k] + noisy[k] for k in sorted(clean)]))


def train_one(model: MetaSieveV1, cells: list[dict], tensors: dict, *,
              family: str, schedule: str, add_support_noise: bool,
              seed: int, config: V1TrainConfig, y_scale: float) -> tuple[dict, TaskScheduler | None]:
    if schedule not in {"uniform", "ats", "ats_null"}:
        raise ValueError("unknown V1 task schedule")
    optimizer = torch.optim.Adam(
        model.parameters(), lr=config.learning_rate,
        weight_decay=config.weight_decay)
    scheduler = (TaskScheduler(config.scheduler_hidden_dim).to(tensors["y"].device)
                 if schedule != "uniform" else None)
    scheduler_optimizer = (torch.optim.Adam(
        scheduler.parameters(), lr=config.scheduler_learning_rate)
        if scheduler is not None else None)
    parameters = tuple(model.parameters())
    source_rng = np.random.default_rng(seed)
    validation_rng = np.random.default_rng(seed + 17_171)
    noise_rng = torch.Generator(device=tensors["y"].device).manual_seed(seed + 31_337)
    sample_rng = torch.Generator(device=tensors["y"].device).manual_seed(seed + 91_919)
    null_rng = torch.Generator(device=tensors["y"].device).manual_seed(seed + 131_071)
    tasks_by_k = {k: build_tasks(cells, "meta_train", k, config.min_query)
                  for k in config.support_sizes}
    val_tasks_by_k = {k: build_tasks(cells, "meta_val", k, config.min_query)
                      for k in config.support_sizes}
    clusters_by_k = {k: cluster_tasks(cells, tasks) for k, tasks in tasks_by_k.items()}
    val_clusters_by_k = {k: cluster_tasks(cells, tasks)
                         for k, tasks in val_tasks_by_k.items()}
    best_state, best_scheduler, best_score, best_step = None, None, float("inf"), 0
    losses, scheduler_losses = [], []

    for step in range(1, config.steps + 1):
        model.train()
        k = config.support_sizes[(step - 1) % len(config.support_sizes)]
        candidate_targets = sample_cluster_balanced(
            clusters_by_k[k], source_rng, config.candidate_tasks)
        episodes = [draw_episode(
            tasks_by_k[k][target], source_rng, k, config.query_per_task)
            for target in candidate_targets]
        support_index, query_index, query_mask = pack_episode_indices(
            episodes, tensors["y"].device)
        noise_sd = config.support_noise_half_width_pk / math.sqrt(3.0) / y_scale
        support_y = tensors["y"][support_index]
        if add_support_noise:
            support_y = add_episodewise_label_noise(
                support_y, noise_sd, generator=noise_rng)
        prediction, support_population = model.batched_episode(
            tensors["ligand"][support_index], tensors[family][support_index],
            support_y, tensors["ligand"][query_index],
            tensors[family][query_index])
        squared_error = (prediction - tensors["y"][query_index]).square()
        episode_losses = (
            (squared_error * query_mask).sum(dim=1)
            / query_mask.sum(dim=1).to(squared_error.dtype))

        if scheduler is not None:
            support_losses = (support_population - support_y).square().mean(dim=1)
            support_gradients = _flat_batched_gradients(
                support_losses, parameters, retain_graph=True).detach()
            query_gradients = _flat_batched_gradients(
                episode_losses, parameters, retain_graph=True).detach()
            progress = episode_losses.new_full(
                (len(episode_losses),),
                (step - 1) / max(1, config.steps - 1))
            statistics = torch.stack((
                torch.log1p(episode_losses.detach()),
                gradient_cosine_rows(support_gradients, query_gradients),
                progress,
            ), dim=1)

        if scheduler is None:
            selected = torch.arange(config.selected_tasks, device=tensors["y"].device)
        else:
            val_target = sample_cluster_balanced(
                val_clusters_by_k[k], validation_rng, 1)[0]
            val_support, val_query = draw_episode(
                val_tasks_by_k[k][val_target], validation_rng, k,
                config.query_per_task)
            val_loss = _episode_loss(
                model, tensors, family, val_support, val_query,
                tensors["y"][val_support])
            val_gradient = _flat_gradient(
                val_loss, parameters, retain_graph=False).detach()
            utility = gradient_cosine_rows(
                query_gradients, val_gradient.expand_as(query_gradients))
            scheduler_input = statistics
            if schedule == "ats_null":
                scheduler_input = scheduler_input[
                    torch.randperm(len(scheduler_input), generator=null_rng,
                                   device=scheduler_input.device)]
            scheduler_optimizer.zero_grad()
            logits = scheduler(scheduler_input)
            target_probability = torch.softmax(
                utility / config.scheduler_target_temperature, dim=0)
            scheduler_loss = -(target_probability * torch.log_softmax(logits, dim=0)).sum()
            scheduler_loss.backward()
            scheduler_optimizer.step()
            scheduler_losses.append(scheduler_loss.detach())
            selected = scheduler.sample(
                scheduler_input.detach(), config.selected_tasks, sample_rng)

        optimizer.zero_grad()
        loss = episode_losses[selected].mean()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
        optimizer.step()
        losses.append(loss.detach())

        if step % config.val_interval == 0 or step == config.steps:
            clean = validation_metrics(
                model, cells, tensors, family, "meta_val", config,
                seed + step, y_scale, noisy_support=False)
            noisy = validation_metrics(
                model, cells, tensors, family, "meta_val", config,
                seed + step, y_scale, noisy_support=True)
            score = _mean_validation_score(clean, noisy)
            if score < best_score:
                best_score, best_step = score, step
                best_state = copy.deepcopy(model.state_dict())
                best_scheduler = (copy.deepcopy(scheduler.state_dict())
                                  if scheduler is not None else None)
    if best_state is None:
        raise RuntimeError("V1 training produced no validation checkpoint")
    model.load_state_dict(best_state)
    if scheduler is not None:
        scheduler.load_state_dict(best_scheduler)
    final_clean = validation_metrics(
        model, cells, tensors, family, "meta_val", config,
        seed + 700_000, y_scale, noisy_support=False)
    final_noisy = validation_metrics(
        model, cells, tensors, family, "meta_val", config,
        seed + 700_000, y_scale, noisy_support=True)
    diagnostics = {
        "best_combined_val_score": best_score,
        "best_step": best_step,
        "clean_val_mse_by_k": final_clean,
        "noisy_support_val_mse_by_k": final_noisy,
        "initial_loss": float(torch.stack(losses[:min(50, len(losses))]).mean()),
        "final_loss": float(torch.stack(losses[-min(50, len(losses)):]).mean()),
        "mean_scheduler_loss": (float(torch.stack(scheduler_losses).mean())
                                if scheduler_losses else None),
    }
    return diagnostics, scheduler


def pretrain_population(model: MetaSieveV1, cells: list[dict], tensors: dict, *,
                        seed: int, config: V1TrainConfig) -> dict:
    """Initialize only the target-independent ligand population on source rows."""
    if config.population_pretrain_steps < 0:
        raise ValueError("population pretrain steps cannot be negative")
    if config.population_pretrain_steps == 0:
        return {"steps": 0, "initial_loss": None, "final_loss": None}
    source = torch.tensor(
        [index for index, row in enumerate(cells) if row["split"] == "meta_train"],
        dtype=torch.long, device=tensors["y"].device)
    if len(source) == 0:
        raise ValueError("population pretraining has no source rows")
    optimizer = torch.optim.Adam(
        model.population.parameters(),
        lr=config.population_pretrain_learning_rate,
        weight_decay=config.weight_decay)
    generator = torch.Generator(device=tensors["y"].device).manual_seed(
        seed + 2_000_003)
    losses = []
    model.train()
    for _ in range(config.population_pretrain_steps):
        selected = source[torch.randint(
            len(source), (min(config.population_pretrain_batch_size, len(source)),),
            generator=generator, device=source.device)]
        optimizer.zero_grad()
        loss = F.mse_loss(
            model.population(tensors["ligand"][selected]).squeeze(-1),
            tensors["y"][selected])
        loss.backward()
        optimizer.step()
        losses.append(loss.detach())
    edge = min(50, len(losses))
    return {
        "steps": config.population_pretrain_steps,
        "initial_loss": float(torch.stack(losses[:edge]).mean()),
        "final_loss": float(torch.stack(losses[-edge:]).mean()),
    }


def _episode_seed(seed: int, target: str, draw: int) -> int:
    key = f"main-v1-cold-target|{seed}|{target}|{draw}"
    return int(hashlib.sha256(key.encode("utf-8")).hexdigest()[:16], 16)


def different_cluster_donors(cells: list[dict], targets: list[str]) -> dict[str, str]:
    target_cluster = {row["target_id"]: row["protein_group_40"] for row in cells}
    donors = {}
    for position, target in enumerate(targets):
        for offset in range(1, len(targets)):
            candidate = targets[(position + offset) % len(targets)]
            if target_cluster[candidate] != target_cluster[target]:
                donors[target] = candidate
                break
        if target not in donors:
            raise ValueError("foreign-support control needs another protein cluster")
    return donors


def predict(model: MetaSieveV1, cells: list[dict], tensors: dict,
            tasks: dict[str, np.ndarray], *, family: str, seed: int,
            k: int, draws: int, max_query: int, arm: str,
            y_scale: float, max_support_k: int, support_mode: str = "correct",
            query_family: str | None = None,
            noisy_support: bool = False,
            noise_half_width_pk: float = 0.2) -> list[dict]:
    targets = sorted(tasks)
    donor = different_cluster_donors(cells, targets)
    rows = []
    model.eval()
    with torch.no_grad():
        episodes = []
        support_values = []
        metadata = []
        for target in targets:
            for draw in range(draws):
                rng = np.random.default_rng(_episode_seed(seed, target, draw))
                support, query = draw_nested_episode(
                    tasks[target], rng, k, max_support_k, max_query)
                source_support = support
                source_family = family
                support_y = tensors["y"][support]
                if support_mode == "foreign":
                    donor_rng = np.random.default_rng(
                        _episode_seed(seed, donor[target], draw))
                    source_support, _ = draw_nested_episode(
                        tasks[donor[target]], donor_rng, k, max_support_k,
                        max_query)
                    support_y = tensors["y"][source_support]
                elif support_mode == "label_control":
                    if k > 1:
                        support_y = support_y.roll(1)
                    else:
                        donor_rng = np.random.default_rng(
                            _episode_seed(seed, donor[target], draw))
                        donor_support, _ = draw_nested_episode(
                            tasks[donor[target]], donor_rng, k, max_support_k,
                            max_query)
                        support_y = tensors["y"][donor_support]
                if noisy_support:
                    noise_seed = _episode_seed(seed + 9_000_000, target, draw)
                    torch_rng = torch.Generator(
                        device=tensors["y"].device).manual_seed(noise_seed)
                    support_y = uniform_label_noise(
                        support_y, noise_half_width_pk / math.sqrt(3.0) / y_scale,
                        generator=torch_rng)
                episodes.append((source_support, query))
                support_values.append(support_y)
                metadata.append((target, draw, source_support, query))
        support_index, query_index, _ = pack_episode_indices(
            episodes, tensors["y"].device)
        support_y = torch.stack(support_values)
        query_feature = query_family or family
        if support_mode == "zero":
            prediction, _ = model.batched_components(
                tensors["ligand"][query_index], tensors[query_feature][query_index])
        else:
            prediction, _ = model.batched_episode(
                tensors["ligand"][support_index],
                tensors[source_family][support_index], support_y,
                tensors["ligand"][query_index], tensors[query_feature][query_index])
        prediction = prediction.cpu().numpy()
        for episode_index, (target, draw, source_support, query) in enumerate(metadata):
            support_hash = _cell_hash(cells, source_support)
            query_hash = _cell_hash(cells, query)
            for index, value in zip(query, prediction[episode_index, :len(query)]):
                rows.append({
                    "arm": arm, "seed": seed, "k": k,
                    "target_id": target, "draw": draw,
                    "cell_id": cells[int(index)]["cell_id"],
                    "support_hash": support_hash, "query_hash": query_hash,
                    "prediction_standardized": float(value),
                })
    return rows


def _cell_hash(cells: list[dict], indices: np.ndarray) -> str:
    content = "|".join(sorted(cells[int(index)]["cell_id"] for index in indices))
    return hashlib.sha256(content.encode("utf-8")).hexdigest()


def write_predictions(path: Path, rows: list[dict]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    writer = io.TextIOWrapper(compressed, encoding="utf-8", newline="\n")
    try:
        for row in rows:
            writer.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
    finally:
        writer.close()
        raw.close()
    return sha256(path)


def _difference_mse(y: np.ndarray, prediction: np.ndarray) -> float:
    if len(y) < 2:
        return float("nan")
    left, right = np.triu_indices(len(y), k=1)
    return float(np.square(
        (prediction[left] - prediction[right]) - (y[left] - y[right])).mean())


def score_predictions(rows: list[dict], cells: list[dict],
                      y_mean: float, y_scale: float):
    truth = {row["cell_id"]: float(row["pK"]) for row in cells}
    target_cluster = {row["target_id"]: row["protein_group_40"] for row in cells}
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(row["k"], row["arm"], row["seed"],
                 row["target_id"], row["draw"])].append(row)
    draw_metrics = []
    for (k, arm, seed, target, draw), values in grouped.items():
        y = np.asarray([truth[row["cell_id"]] for row in values])
        prediction = np.asarray([
            row["prediction_standardized"] * y_scale + y_mean for row in values])
        squared = np.square(y - prediction)
        variance = np.square(y - y.mean()).sum()
        draw_metrics.append({
            "k": k, "arm": arm, "seed": seed, "target_id": target,
            "protein_group_40": target_cluster[target], "draw": draw,
            "mse": float(squared.mean()),
            "rmse": float(np.sqrt(squared.mean())),
            "r2": float(1.0 - squared.sum() / variance) if variance > 0 else np.nan,
            "ci": concordance(y, prediction),
            "pearson": float(stats.pearsonr(y, prediction).statistic)
            if y.std() > 0 and prediction.std() > 0 else np.nan,
            "spearman": float(stats.spearmanr(y, prediction).statistic)
            if y.std() > 0 and prediction.std() > 0 else np.nan,
            "difference_mse": _difference_mse(y, prediction),
        })
    target_metrics: dict[tuple, dict[str, float]] = {}
    target_rows: dict[tuple, list[dict]] = defaultdict(list)
    for row in draw_metrics:
        target_rows[(row["k"], row["arm"], row["target_id"])].append(row)
    metric_names = ("mse", "rmse", "r2", "ci", "pearson", "spearman",
                    "difference_mse")
    for key, values in target_rows.items():
        target_metrics[key] = {
            name: float(np.nanmean([row[name] for row in values]))
            for name in metric_names}
    summary: dict[str, dict[str, dict[str, float]]] = defaultdict(dict)
    for k in sorted({row["k"] for row in draw_metrics}):
        for arm in sorted({row["arm"] for row in draw_metrics if row["k"] == k}):
            values = [metrics for (value_k, value_arm, _), metrics in target_metrics.items()
                      if value_k == k and value_arm == arm]
            summary[str(k)][arm] = {
                name: float(np.nanmean([value[name] for value in values]))
                for name in metric_names}
    return dict(summary), target_metrics, draw_metrics, target_cluster


def cluster_contrast(target_metrics: dict, target_cluster: dict, *,
                     k: int, correct: str, control: str, metric: str,
                     draws: int, seed: int) -> dict:
    targets = sorted(
        {target for value_k, arm, target in target_metrics
         if value_k == k and arm == correct}
        & {target for value_k, arm, target in target_metrics
           if value_k == k and arm == control})
    grouped: dict[str, list[float]] = defaultdict(list)
    for target in targets:
        grouped[target_cluster[target]].append(
            target_metrics[(k, control, target)][metric]
            - target_metrics[(k, correct, target)][metric])
    component = np.asarray([np.mean(values) for _, values in sorted(grouped.items())])
    if len(component) < 2:
        raise ValueError("cold-target contrast needs at least two protein components")
    rng = np.random.default_rng(seed)
    samples = component[
        rng.integers(0, len(component), size=(draws, len(component)))].mean(axis=1)
    return {
        "k": k, "correct": correct, "control": control, "metric": metric,
        "components": len(component),
        "cluster_macro_reduction": float(component.mean()),
        "one_sided_95_lcb": float(np.quantile(samples, 0.05)),
        "pass": bool(np.quantile(samples, 0.05) > 0),
    }


def run(config: V1TrainConfig = V1TrainConfig(), *, device: str = "cuda",
        seeds: tuple[int, ...] = SEEDS, output: Path = OUT) -> dict:
    run_started = time.perf_counter()
    if output.exists():
        raise FileExistsError(f"V1 output already exists: {output}")
    device = str(require_cuda(device))
    output.mkdir(parents=True, exist_ok=True)
    cells, tensors, _, normalization = load_data(device)
    y_scale = normalization["y_scale"]
    runner_sha256 = sha256(Path(__file__))
    model_sha256 = sha256(ROOT / "model/metasieve_v1.py")
    models: dict[tuple[str, int], MetaSieveV1] = {}
    diagnostics, checkpoints = [], []

    specifications = list(VARIANTS.items()) + [
        ("population_d0", ("uniform", False)),
        ("ligand_only", ("uniform", False)),
    ]
    for seed in seeds:
        for name, (schedule, noise) in specifications:
            family = "ligand" if name == "ligand_only" else "correct"
            dimension = 0 if name == "population_d0" else config.section_dim
            torch.manual_seed(seed)
            model = MetaSieveV1(
                288, dimension, config.ridge,
                support_only_section=config.support_only_section,
                population_hidden_dim=config.population_hidden_dim,
                pair_hidden_dim=config.pair_hidden_dim).to(device)
            checkpoint_path = output / "checkpoints" / f"{name}_seed{seed}.pt"
            checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
            expected = {
                "seed": seed, "arm": name, "config": asdict(config),
                "corpus_sha256": sha256(CORPUS / "manifest.json"),
                "features_sha256": sha256(FEATURES),
                "runner_sha256": runner_sha256,
                "model_sha256": model_sha256,
            }
            if checkpoint_path.exists():
                checkpoint = torch.load(
                    checkpoint_path, map_location=device, weights_only=False)
                if any(checkpoint.get(key) != value for key, value in expected.items()):
                    raise ValueError(f"incompatible V1 resume checkpoint: {checkpoint_path}")
                model.load_state_dict(checkpoint["model_state"])
                detail = checkpoint["train_diagnostics"] | {"resumed_checkpoint": True}
            else:
                training_started = time.perf_counter()
                population_pretrain = pretrain_population(
                    model, cells, tensors, seed=seed, config=config)
                detail, scheduler = train_one(
                    model, cells, tensors, family=family, schedule=schedule,
                    add_support_noise=noise, seed=seed, config=config,
                    y_scale=y_scale)
                detail["population_pretrain"] = population_pretrain
                detail["training_elapsed_seconds"] = (
                    time.perf_counter() - training_started)
                torch.save({
                    **expected,
                    "model_state": model.state_dict(),
                    "scheduler_state": (scheduler.state_dict()
                                        if scheduler is not None else None),
                    "train_diagnostics": detail,
                    "inference_requires_scheduler": False,
                }, checkpoint_path)
            models[name, seed] = model
            diagnostics.append({"seed": seed, "arm": name, **detail})
            checkpoints.append({
                "path": str(checkpoint_path.relative_to(ROOT)),
                "sha256": sha256(checkpoint_path),
            })

    variant_scores = {
        name: float(np.mean([
            next(row["best_combined_val_score"] for row in diagnostics
                 if row["arm"] == name and row["seed"] == seed)
            for seed in seeds]))
        for name in VARIANTS}
    ats_meta_validation_gate = {
        "clean": variant_scores["ats_clean"] < variant_scores["ats_null_clean"],
        "support_noise": (
            variant_scores["ats_support_noise"]
            < variant_scores["ats_null_support_noise"]),
    }
    eligible_variants = ["uniform_clean", "uniform_support_noise"]
    if ats_meta_validation_gate["clean"]:
        eligible_variants.append("ats_clean")
    if ats_meta_validation_gate["support_noise"]:
        eligible_variants.append("ats_support_noise")
    selected_variant = min(
        eligible_variants, key=lambda name: (variant_scores[name], name))

    all_rows = []
    max_support_k = max(config.support_sizes)
    evaluation_tasks = build_tasks(
        cells, "meta_test", max_support_k, config.min_query)
    for seed in seeds:
        for k in config.support_sizes:
            for variant in VARIANTS:
                model = models[variant, seed]
                all_rows += predict(
                    model, cells, tensors, evaluation_tasks, family="correct", seed=seed,
                    k=k, draws=config.test_draws, max_query=config.test_max_query,
                    arm=f"{variant}_correct", y_scale=y_scale,
                    max_support_k=max_support_k)
                all_rows += predict(
                    model, cells, tensors, evaluation_tasks, family="correct", seed=seed,
                    k=k, draws=config.test_draws, max_query=config.test_max_query,
                    arm=f"{variant}_support_noisy", y_scale=y_scale,
                    max_support_k=max_support_k,
                    noisy_support=True,
                    noise_half_width_pk=config.support_noise_half_width_pk)
            selected = models[selected_variant, seed]
            for arm, kwargs in (
                    ("selected_zero", {"support_mode": "zero"}),
                    ("selected_label_control", {"support_mode": "label_control"}),
                    ("selected_foreign", {"support_mode": "foreign"}),
                    ("selected_support_wrong", {
                        "family": "wrong", "query_family": "correct"}),
                    ("selected_query_wrong", {"query_family": "wrong"}),
                    ("selected_both_wrong", {"family": "wrong", "query_family": "wrong"})):
                arguments = {"family": "correct", **kwargs}
                all_rows += predict(
                    selected, cells, tensors, evaluation_tasks, seed=seed, k=k,
                    draws=config.test_draws, max_query=config.test_max_query,
                    arm=arm, y_scale=y_scale, max_support_k=max_support_k,
                    **arguments)
            all_rows += predict(
                models["population_d0", seed], cells, tensors, evaluation_tasks,
                family="correct", seed=seed, k=k, draws=config.test_draws,
                max_query=config.test_max_query, arm="population_d0", y_scale=y_scale,
                max_support_k=max_support_k)
            all_rows += predict(
                models["ligand_only", seed], cells, tensors, evaluation_tasks,
                family="ligand", seed=seed, k=k, draws=config.test_draws,
                max_query=config.test_max_query, arm="ligand_only", y_scale=y_scale,
                max_support_k=max_support_k)

    prediction_path = output / "predictions_without_query_labels.jsonl.gz"
    prediction_hash = write_predictions(prediction_path, all_rows)
    summary, target_metrics, draw_metrics, target_cluster = score_predictions(
        all_rows, cells, normalization["y_mean"], y_scale)
    selected_arm = f"{selected_variant}_correct"
    comparisons = [
        ("population_d0", "mse"), ("ligand_only", "mse"),
        ("selected_zero", "mse"), ("selected_label_control", "mse"),
        ("selected_foreign", "mse"),
        ("selected_support_wrong", "difference_mse"),
        ("selected_query_wrong", "difference_mse"),
        ("selected_both_wrong", "difference_mse"),
    ]
    contrasts = []
    for k in config.support_sizes:
        for index, (control, metric) in enumerate(comparisons):
            contrasts.append(cluster_contrast(
                target_metrics, target_cluster, k=k, correct=selected_arm,
                control=control, metric=metric, draws=config.bootstrap_draws,
                seed=seeds[0] + 1000 * k + index))
    contrast_lookup = {(row["k"], row["control"], row["metric"]): row
                       for row in contrasts}
    mechanism_specs = (
        ("ats_clean_correct", "ats_null_clean_correct", "ats_clean_vs_null"),
        ("ats_support_noise_correct", "ats_null_support_noise_correct",
         "ats_support_noise_vs_null"),
        ("uniform_support_noise_correct", "uniform_clean_correct",
         "uniform_noise_on_clean_support"),
        ("uniform_support_noise_support_noisy", "uniform_clean_support_noisy",
         "uniform_noise_on_noisy_support"),
        ("ats_support_noise_correct", "ats_clean_correct",
         "ats_noise_on_clean_support"),
        ("ats_support_noise_support_noisy", "ats_clean_support_noisy",
         "ats_noise_on_noisy_support"),
    )
    mechanism_contrasts = []
    for k in config.support_sizes:
        for index, (correct, control, mechanism) in enumerate(mechanism_specs):
            row = cluster_contrast(
                target_metrics, target_cluster, k=k, correct=correct,
                control=control, metric="mse", draws=config.bootstrap_draws,
                seed=seeds[0] + 20_000 + 1000 * k + index)
            row["mechanism"] = mechanism
            mechanism_contrasts.append(row)
    mechanism_lookup = {
        (row["k"], row["mechanism"]): row for row in mechanism_contrasts}
    scheduler_test_gate = {
        condition: ats_meta_validation_gate[condition] and all(
            mechanism_lookup[(k, f"ats_{condition}_vs_null")]["pass"]
            for k in config.support_sizes)
        for condition in ("clean", "support_noise")
    }
    noise_test_gate = {
        schedule: all(
            mechanism_lookup[(k, f"{schedule}_noise_on_clean_support")]["pass"]
            and mechanism_lookup[(k, f"{schedule}_noise_on_noisy_support")]["pass"]
            for k in config.support_sizes)
        for schedule in ("uniform", "ats")
    }
    performance_thresholds = {
        "rmse_pk_at_most": 1.0,
        "concordance_at_least": 0.60,
        "primary_k": [1, 2, 3],
        "k5_must_not_exceed_rmse_pk": 1.1,
    }
    performance_by_k = {}
    for k in config.support_sizes:
        point = summary[str(k)][selected_arm]
        performance_by_k[str(k)] = {
            "rmse_pass": point["rmse"] <= (1.1 if k == 5 else 1.0),
            "ci_pass": point["ci"] >= 0.60,
            "beats_ligand_only_component_lcb": contrast_lookup[
                (k, "ligand_only", "mse")]["pass"],
            "beats_population_component_lcb": contrast_lookup[
                (k, "population_d0", "mse")]["pass"],
        }
    primary_good = all(
        all(performance_by_k[str(k)].values()) for k in (1, 2, 3)) \
        and performance_by_k["5"]["rmse_pass"]
    biology_gate = all(
        contrast_lookup[(k, control, "difference_mse")]["pass"]
        for k in config.support_sizes
        for control in ("selected_support_wrong", "selected_query_wrong",
                        "selected_both_wrong"))
    support_gate = all(
        contrast_lookup[(k, control, "mse")]["pass"]
        for k in config.support_sizes
        for control in ("selected_zero", "selected_label_control", "selected_foreign"))
    selected_mechanism_gate = True
    if selected_variant.startswith("ats_"):
        selected_mechanism_gate = scheduler_test_gate[
            "support_noise" if selected_variant.endswith("support_noise") else "clean"]
    if selected_variant.endswith("support_noise"):
        selected_mechanism_gate = selected_mechanism_gate and noise_test_gate[
            "ats" if selected_variant.startswith("ats_") else "uniform"]

    result = json_safe({
        "schema": "MetaSieve.MainV1ColdTargetResult.v1",
        "config": asdict(config),
        "seeds": list(seeds),
        "environment": {
            "python": subprocess.check_output(["python", "--version"], text=True).strip(),
            "torch": torch.__version__, "device": device,
            "git_head": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
            "runner_sha256": runner_sha256,
            "model_sha256": model_sha256,
            "elapsed_seconds_this_run": time.perf_counter() - run_started,
            "rng_contract": (
                "paired model initialization, source episodes, scheduler sampling, "
                "scheduler-null permutation, support noise, and validation use "
                "separate deterministic streams"),
        },
        "data": {
            "corpus_sha256": sha256(CORPUS / "manifest.json"),
            "features_sha256": sha256(FEATURES),
            "cold_target_tasks_by_k": {
                str(k): len(build_tasks(cells, "meta_test", k, config.min_query))
                for k in config.support_sizes},
            "evaluation_common_cohort_tasks": len(evaluation_tasks),
            "evaluation_episode_contract": (
                "same max-k-eligible targets, shared permutation, nested support "
                "prefixes, and k-invariant query for k=1/2/3/5"),
            "split": "complete CD-HIT-40 target clusters; consumed development test",
            "split_limitations": (
                "homology-cold target only; ligand scaffold and document closures "
                "are not enforced across splits"),
        },
        "adambind_translation": {
            "core_reference": "Wan et al., Nature Communications 2026",
            "task_scheduler": "source-only first-order meta-validation gradient alignment",
            "candidate_tasks": config.candidate_tasks,
            "selected_without_replacement": config.selected_tasks,
            "support_noise_distribution_pk": (
                f"Uniform(-{config.support_noise_half_width_pk},"
                f"+{config.support_noise_half_width_pk})"),
            "query_label_noise": False,
            "scheduler_loaded_at_test": False,
        },
        "variant_meta_validation_score": variant_scores,
        "ats_meta_validation_beats_matched_null": ats_meta_validation_gate,
        "eligible_variants_before_test_scoring": eligible_variants,
        "selected_variant_before_test_scoring": selected_variant,
        "prediction_artifact": {
            "path": str(prediction_path.relative_to(ROOT)),
            "sha256_before_scoring": prediction_hash,
            "query_labels_in_artifact": False,
            "firewall_scope": (
                "development code-path firewall; query labels are loaded by the "
                "runner but are not indexed by predict() or serialized"),
            "rows": len(all_rows),
        },
        "point_metrics_target_macro": summary,
        "cluster_bootstrap_contrasts": contrasts,
        "training_intervention_cluster_contrasts": mechanism_contrasts,
        "performance_thresholds": performance_thresholds,
        "performance_by_k": performance_by_k,
        "gates": {
            "cold_target_performance_good": primary_good,
            "support_specificity_all_k": support_gate,
            "partner_identity_difference_all_k": biology_gate,
            "scheduler_beats_matched_null_all_k": scheduler_test_gate,
            "support_noise_beats_clean_training_all_k": noise_test_gate,
            "selected_training_mechanism_identified": selected_mechanism_gate,
        },
        "train_diagnostics": diagnostics,
        "checkpoints": checkpoints,
        "development_only": True,
        "production_migration_authorized": False,
        "TERMINAL_VERDICT": (
            "COLD_TARGET_FEWSHOT_V1_GOOD_DEVELOPMENT"
            if primary_good and support_gate and biology_gate and selected_mechanism_gate
            else "COLD_TARGET_FEWSHOT_V1_NOT_YET_GOOD"
        ),
    })
    (output / "MAIN_V1_COLD_TARGET_RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    (output / "draw_metrics.json").write_text(
        json.dumps(json_safe(draw_metrics), indent=2, sort_keys=True,
                   allow_nan=False) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--support-only-section", action="store_true")
    parser.add_argument("--population-hidden-dim", type=int, default=0)
    parser.add_argument("--pair-hidden-dim", type=int, default=0)
    parser.add_argument("--section-dim", type=int)
    parser.add_argument("--ridge", type=float)
    parser.add_argument("--population-pretrain-steps", type=int, default=0)
    args = parser.parse_args()
    config = V1TrainConfig()
    if (args.support_only_section or args.population_hidden_dim
            or args.pair_hidden_dim or args.section_dim is not None
            or args.ridge is not None or args.population_pretrain_steps):
        config = replace(
            config, support_only_section=args.support_only_section,
            population_hidden_dim=args.population_hidden_dim,
            pair_hidden_dim=args.pair_hidden_dim,
            section_dim=(args.section_dim if args.section_dim is not None
                         else config.section_dim),
            ridge=(args.ridge if args.ridge is not None else config.ridge))
        config = replace(
            config, population_pretrain_steps=args.population_pretrain_steps)
    seeds = SEEDS
    output = args.output
    if not output.is_absolute():
        output = ROOT / output
    output = output.resolve()
    allowed_output = (ROOT / "report/meta_fewshot").resolve()
    if output == allowed_output or not output.is_relative_to(allowed_output):
        raise ValueError(
            f"V1 output must be a new child of {allowed_output}: {output}")
    if args.smoke:
        config = replace(
            config, steps=4, val_interval=4, val_draws=1,
            test_draws=1, bootstrap_draws=199,
            candidate_tasks=8, selected_tasks=4)
        seeds = (SEEDS[0],)
    result = run(config, device=args.device, seeds=seeds, output=output)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
