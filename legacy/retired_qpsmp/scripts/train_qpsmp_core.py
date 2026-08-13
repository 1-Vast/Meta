"""Train the QPSMP core on the main-v0 cold-target few-shot corpus.

This is the first endpoint-scalar training smoke for the CURRENT_THEORY QPSMP
architecture. It is intentionally small and compares against previous V1
metrics rather than claiming admission to G2/G3.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass
import gzip
import hashlib
import json
import time
from pathlib import Path
import sys

import numpy as np
import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from model.qpsmp import QPSMPCore
from research.meta_fewshot.train_main_v0 import (
    CORPUS,
    FEATURES,
    json_safe,
    load_data,
    sha256,
)
from research.meta_fewshot.train_main_v1 import (
    build_tasks,
    cluster_contrast,
    draw_episode,
    draw_nested_episode,
    pack_episode_indices,
    score_predictions,
)


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "report/meta_fewshot/qpsmp_core_smoke"
SEEDS = (20260841,)


@dataclass(frozen=True)
class QPSMPTrainConfig:
    support_sizes: tuple[int, ...] = (1, 2, 3, 5)
    min_query: int = 3
    task_dim: int = 2
    ridge: float = 1.0
    section_radius_bound: float = 1.0
    tasks_per_step: int = 8
    query_per_task: int = 32
    steps: int = 300
    learning_rate: float = 0.005
    weight_decay: float = 1e-4
    val_interval: int = 50
    val_draws: int = 1
    test_draws: int = 3
    test_max_query: int = 128
    bootstrap_draws: int = 999


def feature_tensors(tensors: dict, family: str) -> torch.Tensor:
    return tensors[family]


def baseline_tensors(tensors: dict) -> torch.Tensor:
    return tensors["ligand"]


def validation_mse(
        model: QPSMPCore, cells: list[dict], tensors: dict,
        split: str, config: QPSMPTrainConfig, *,
        seed: int, y_scale: float) -> float:
    model.eval()
    tasks = build_tasks(cells, split, max(config.support_sizes), config.min_query)
    losses = []
    with torch.no_grad():
        for target, indices in tasks.items():
            target_losses = []
            for draw in range(config.val_draws):
                episode_seed = stable_seed("qpsmp-val", seed, target, draw)
                support, query = draw_nested_episode(
                    indices, np.random.default_rng(episode_seed),
                    max(config.support_sizes), max(config.support_sizes),
                    config.query_per_task)
                output = model(
                    feature_tensors(tensors, "correct")[support],
                    tensors["y"][support],
                    feature_tensors(tensors, "correct")[query],
                    support_baseline_features=baseline_tensors(tensors)[support],
                    query_baseline_features=baseline_tensors(tensors)[query])
                target_losses.append(float(
                    ((output.ridge_prediction - tensors["y"][query]) * y_scale)
                    .square().mean()))
            losses.append(float(np.mean(target_losses)))
    return float(np.mean(losses))


def stable_seed(*parts: object) -> int:
    digest = hashlib.sha256("|".join(map(str, parts)).encode("utf-8")).hexdigest()
    return int(digest[:16], 16)


def train_one(
        cells: list[dict], tensors: dict, config: QPSMPTrainConfig,
        *, seed: int, y_scale: float) -> tuple[QPSMPCore, dict]:
    torch.manual_seed(seed)
    model = QPSMPCore(
        feature_dim=tensors["correct"].shape[1],
        baseline_dim=tensors["ligand"].shape[1],
        task_dim=config.task_dim,
        ridge=config.ridge,
        section_radius_bound=config.section_radius_bound).to(tensors["y"].device)
    optimizer = torch.optim.Adam(
        model.parameters(), lr=config.learning_rate, weight_decay=config.weight_decay)
    rng = np.random.default_rng(seed)
    tasks_by_k = {
        k: build_tasks(cells, "meta_train", k, config.min_query)
        for k in config.support_sizes
    }
    best_state, best_val, best_step = None, float("inf"), 0
    losses = []
    for step in range(1, config.steps + 1):
        k = config.support_sizes[(step - 1) % len(config.support_sizes)]
        targets = sorted(tasks_by_k[k])
        selected = rng.choice(len(targets), size=config.tasks_per_step, replace=True)
        optimizer.zero_grad()
        loss = torch.zeros((), device=tensors["y"].device)
        for target_index in selected:
            indices = tasks_by_k[k][targets[int(target_index)]]
            support, query = draw_episode(
                indices, rng, k, max_query=config.query_per_task)
            output = model(
                feature_tensors(tensors, "correct")[support],
                tensors["y"][support],
                feature_tensors(tensors, "correct")[query],
                support_baseline_features=baseline_tensors(tensors)[support],
                query_baseline_features=baseline_tensors(tensors)[query])
            loss = loss + (
                output.ridge_prediction - tensors["y"][query]).square().mean()
        loss = loss / config.tasks_per_step
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
        optimizer.step()
        losses.append(float(loss.detach()))
        if step % config.val_interval == 0 or step == config.steps:
            value = validation_mse(
                model, cells, tensors, "meta_val", config,
                seed=seed, y_scale=y_scale)
            if value < best_val:
                best_state = copy.deepcopy(model.state_dict())
                best_val, best_step = value, step
    assert best_state is not None
    model.load_state_dict(best_state)
    return model, {
        "best_step": best_step,
        "best_val_mse_pk": best_val,
        "initial_loss": float(np.mean(losses[:min(25, len(losses))])),
        "final_loss": float(np.mean(losses[-min(25, len(losses)):]))}


def predict_rows(
        model: QPSMPCore, cells: list[dict], tensors: dict,
        tasks: dict[str, np.ndarray], *, seed: int, k: int,
        draws: int, max_query: int, arm: str, family: str = "correct",
        query_family: str | None = None,
        support_mode: str = "real", max_support_k: int = 5) -> list[dict]:
    rows = []
    query_family = family if query_family is None else query_family
    targets = sorted(tasks)
    donor = {target: targets[(index + 1) % len(targets)] for index, target in enumerate(targets)}
    model.eval()
    with torch.no_grad():
        for target in targets:
            for draw in range(draws):
                episode_seed = stable_seed("qpsmp-test", seed, target, draw)
                support, query = draw_nested_episode(
                    tasks[target], np.random.default_rng(episode_seed),
                    k, max_support_k, max_query)
                source_support = support
                source_y = tensors["y"][support]
                source_family = family
                adapt = support_mode != "zero"
                task_state_override = None
                if support_mode == "foreign":
                    donor_seed = stable_seed("qpsmp-donor", seed, target, draw)
                    donor_support, _ = draw_episode(
                        tasks[donor[target]], np.random.default_rng(donor_seed),
                        k, max_query=max_query)
                    task_state_override = model.centered_state(
                        feature_tensors(tensors, family)[donor_support],
                        tensors["y"][donor_support],
                        baseline_tensors(tensors)[donor_support],
                    )[0]
                prediction = model(
                    feature_tensors(tensors, source_family)[source_support],
                    source_y,
                    feature_tensors(tensors, query_family)[query],
                    support_baseline_features=baseline_tensors(tensors)[source_support],
                    query_baseline_features=baseline_tensors(tensors)[query],
                    adapt=adapt,
                    task_state_override=task_state_override).ridge_prediction
                for cell_index, value in zip(query, prediction):
                    cell = cells[int(cell_index)]
                    rows.append({
                        "k": k,
                        "arm": arm,
                        "seed": seed,
                        "target_id": target,
                        "draw": draw,
                        "cell_id": cell["cell_id"],
                        "prediction_standardized": float(value),
                    })
    return rows


def write_predictions(path: Path, rows: list[dict]) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", compresslevel=6) as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    return sha256(path)


def run(
        config: QPSMPTrainConfig = QPSMPTrainConfig(), *,
        device: str = "cuda", output: Path = OUT,
        seeds: tuple[int, ...] = SEEDS) -> dict:
    if output.exists():
        raise FileExistsError(f"QPSMP output already exists: {output}")
    started = time.perf_counter()
    resolved = torch.device(device)
    if resolved.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    device = str(resolved)
    output.mkdir(parents=True, exist_ok=False)
    cells, tensors, _, normalization = load_data(device)
    y_scale = normalization["y_scale"]
    models = {}
    diagnostics = []
    for seed in seeds:
        model, detail = train_one(cells, tensors, config, seed=seed, y_scale=y_scale)
        models[seed] = model
        diagnostics.append({"seed": seed, **detail})
    evaluation_tasks = build_tasks(
        cells, "meta_test", max(config.support_sizes), config.min_query)
    rows = []
    for seed, model in models.items():
        for k in config.support_sizes:
            rows += predict_rows(
                model, cells, tensors, evaluation_tasks, seed=seed, k=k,
                draws=config.test_draws, max_query=config.test_max_query,
                arm="qpsmp_correct", max_support_k=max(config.support_sizes))
            rows += predict_rows(
                model, cells, tensors, evaluation_tasks, seed=seed, k=k,
                draws=config.test_draws, max_query=config.test_max_query,
                arm="qpsmp_zero_support", support_mode="zero",
                max_support_k=max(config.support_sizes))
            rows += predict_rows(
                model, cells, tensors, evaluation_tasks, seed=seed, k=k,
                draws=config.test_draws, max_query=config.test_max_query,
                arm="qpsmp_foreign_support", support_mode="foreign",
                max_support_k=max(config.support_sizes))
            rows += predict_rows(
                model, cells, tensors, evaluation_tasks, seed=seed, k=k,
                draws=config.test_draws, max_query=config.test_max_query,
                arm="qpsmp_wrong_query", family="correct", query_family="wrong",
                max_support_k=max(config.support_sizes))
            rows += level_rows(
                cells, tensors, evaluation_tasks, seed=seed, k=k,
                draws=config.test_draws, max_query=config.test_max_query,
                max_support_k=max(config.support_sizes))
    prediction_path = output / "predictions_without_query_labels.jsonl.gz"
    prediction_hash = write_predictions(prediction_path, rows)
    summary, target_metrics, _, target_cluster = score_predictions(
        rows, cells, normalization["y_mean"], y_scale)
    contrasts = []
    for k in config.support_sizes:
        for index, control in enumerate(
                ("level", "qpsmp_zero_support", "qpsmp_foreign_support",
                 "qpsmp_wrong_query")):
            contrasts.append(cluster_contrast(
                target_metrics, target_cluster, k=k, correct="qpsmp_correct",
                control=control, metric="mse", draws=config.bootstrap_draws,
                seed=seeds[0] + 2000 * k + index))
    previous = previous_metrics()
    gates = {
        "qpsmp_beats_level_all_k": all(
            row["pass"] for row in contrasts if row["control"] == "level"),
        "qpsmp_beats_zero_support_all_k": all(
            row["pass"] for row in contrasts if row["control"] == "qpsmp_zero_support"),
        "qpsmp_beats_foreign_support_all_k": all(
            row["pass"] for row in contrasts if row["control"] == "qpsmp_foreign_support"),
        "qpsmp_beats_wrong_query_all_k": all(
            row["pass"] for row in contrasts if row["control"] == "qpsmp_wrong_query"),
    }
    verdict = (
        "QPSMP_CORE_SMOKE_PASS"
        if all(gates.values())
        else "QPSMP_CORE_SMOKE_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.QPSMPCoreSmoke.v1",
        "hypothesis": (
            "The CURRENT_THEORY QPSMP retained scalar path improves cold-target "
            "few-shot endpoint prediction against matched level/support/protein controls."),
        "config": asdict(config),
        "seeds": list(seeds),
        "data": {
            "corpus_sha256": sha256(CORPUS / "manifest.json"),
            "features_sha256": sha256(FEATURES),
            "test_tasks": len(evaluation_tasks),
            "y_scale": y_scale,
        },
        "training_diagnostics": diagnostics,
        "prediction_artifact": {
            "path": str(prediction_path.resolve().relative_to(ROOT)),
            "sha256": prediction_hash,
            "rows": len(rows),
            "query_labels_in_artifact": False,
        },
        "summary": summary,
        "contrasts": contrasts,
        "previous_metrics": previous,
        "gates": gates,
        "elapsed_seconds": time.perf_counter() - started,
        "g2_authorized": False,
        "g3_authorized": False,
        "v1_integration_authorized": False,
        "TERMINAL_VERDICT": verdict,
    }
    (output / "RESULT.json").write_text(
        json.dumps(json_safe(result), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    return result


def level_rows(
        cells: list[dict], tensors: dict, tasks: dict[str, np.ndarray], *,
        seed: int, k: int, draws: int, max_query: int,
        max_support_k: int) -> list[dict]:
    rows = []
    for target, indices in sorted(tasks.items()):
        for draw in range(draws):
            episode_seed = stable_seed("qpsmp-test", seed, target, draw)
            support, query = draw_nested_episode(
                indices, np.random.default_rng(episode_seed),
                k, max_support_k, max_query)
            prediction = torch.full_like(
                tensors["y"][query], float(tensors["y"][support].mean()))
            for cell_index, value in zip(query, prediction):
                rows.append({
                    "k": k,
                    "arm": "level",
                    "seed": seed,
                    "target_id": target,
                    "draw": draw,
                    "cell_id": cells[int(cell_index)]["cell_id"],
                    "prediction_standardized": float(value),
                })
    return rows


def previous_metrics() -> dict:
    path = ROOT / "report/meta_fewshot/v1_targeted_repairs/RESULT.json"
    if not path.exists():
        return {}
    result = json.loads(path.read_text(encoding="utf-8"))
    return {
        "v1_targeted_final_rmse_by_k": result["final_candidate"]["rmse_by_k"],
        "v1_vectorized_baseline_rmse_by_k": result["baseline"]["rmse_by_k"],
        "v1_terminal_verdict": result["terminal_verdict"],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--steps", type=int, default=300)
    parser.add_argument("--task-dim", type=int, default=2)
    parser.add_argument("--ridge", type=float, default=1.0)
    parser.add_argument("--test-draws", type=int, default=3)
    parser.add_argument("--bootstrap-draws", type=int, default=999)
    args = parser.parse_args()
    result = run(
        QPSMPTrainConfig(
            steps=args.steps,
            task_dim=args.task_dim,
            ridge=args.ridge,
            test_draws=args.test_draws,
            bootstrap_draws=args.bootstrap_draws),
        device=args.device,
        output=args.output)
    print(json.dumps(json_safe(result), indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
