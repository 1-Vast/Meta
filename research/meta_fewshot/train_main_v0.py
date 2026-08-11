"""Real BindingDB MetaSieve-main v0 episodic training and controlled evaluation."""
from __future__ import annotations

import copy
import gzip
import hashlib
import io
import json
import math
import random
import subprocess
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from scipy import stats
from torch import nn

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
FEATURES = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_tbasis_features.npz"
OUT = ROOT / "report/meta_fewshot/main_v0"
SEEDS = (20260811, 20260812, 20260813, 20260814, 20260815)


@dataclass(frozen=True)
class TrainConfig:
    k: int = 5
    min_query: int = 3
    tasks_per_step: int = 8
    query_per_task: int = 32
    selection_steps: int = 500
    final_steps: int = 1000
    learning_rate: float = 0.01
    weight_decay: float = 1e-4
    val_interval: int = 100
    val_draws: int = 2
    test_draws: int = 5
    bootstrap_draws: int = 9999


class MetaSectionRegressor(nn.Module):
    def __init__(self, input_dim: int, section_dim: int, ridge: float):
        super().__init__()
        if not 0 <= section_dim <= 5 or ridge <= 0:
            raise ValueError("section_dim must be in [0,5] and ridge positive")
        self.section_dim = section_dim
        self.ridge = float(ridge)
        self.population = nn.Linear(input_dim, 1)
        if section_dim:
            self.raw_basis = nn.Parameter(torch.randn(input_dim, section_dim) / math.sqrt(input_dim))
            self.population_coordinate = nn.Parameter(torch.zeros(section_dim))

    def basis(self):
        return torch.linalg.qr(self.raw_basis, mode="reduced").Q

    def components(self, ligand, family):
        population = self.population(ligand).squeeze(-1)
        if not self.section_dim:
            return population, None
        coordinates = family @ self.basis()
        return population + coordinates @ self.population_coordinate, coordinates

    def episode(self, support_ligand, support_family, support_y,
                query_ligand, query_family):
        support_population, support_coordinates = self.components(support_ligand, support_family)
        query_population, query_coordinates = self.components(query_ligand, query_family)
        if not self.section_dim:
            return query_population
        residual = support_y - support_population
        identity = torch.eye(len(support_y), device=support_y.device, dtype=support_y.dtype)
        dual = torch.linalg.solve(
            support_coordinates @ support_coordinates.T + self.ridge * identity,
            residual,
        )
        coefficient = support_coordinates.T @ dual
        return query_population + query_coordinates @ coefficient


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def json_safe(value):
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    return value


def read_cells(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def ligand_fingerprints(cells: list[dict]) -> np.ndarray:
    ligands = {
        row["drug_key"]: row["smiles"]
        for row in map(json.loads, (CORPUS / "ligands.jsonl").read_text().splitlines())
    }
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=288)
    cache = {}
    for key, smiles in ligands.items():
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"invalid admitted ligand {key}")
        cache[key] = generator.GetFingerprintAsNumPy(molecule).astype(np.float32)
    return np.stack([cache[row["ligand_id"]] for row in cells])


def standardize(values: np.ndarray, indices: np.ndarray) -> tuple[np.ndarray, dict]:
    mean = values[indices].mean(axis=0)
    scale = values[indices].std(axis=0)
    scale[scale < 1e-6] = 1.0
    return ((values - mean) / scale).astype(np.float32), {"mean": mean, "scale": scale}


def load_data(device: str):
    cells = read_cells(CORPUS / "cells.jsonl.gz")
    with np.load(FEATURES, allow_pickle=False) as stored:
        if stored["cell_id"].tolist() != [row["cell_id"] for row in cells]:
            raise ValueError("feature rows do not match corpus cell order")
        correct = stored["correct"].astype(np.float32)
        wrong = stored["deranged_protein"].astype(np.float32)
    if not np.isfinite(correct).all() or not np.isfinite(wrong).all():
        raise ValueError("features contain non-finite values")
    ligand = ligand_fingerprints(cells)
    train_indices = np.asarray([i for i, row in enumerate(cells) if row["split"] == "meta_train"])
    correct, full_norm = standardize(correct, train_indices)
    wrong = ((wrong - full_norm["mean"]) / full_norm["scale"]).astype(np.float32)
    ligand, ligand_norm = standardize(ligand, train_indices)
    y_raw = np.asarray([row["pK"] for row in cells], dtype=np.float32)
    y_mean, y_scale = float(y_raw[train_indices].mean()), float(y_raw[train_indices].std())
    y = (y_raw - y_mean) / y_scale
    tensors = {
        "correct": torch.from_numpy(correct).to(device),
        "wrong": torch.from_numpy(wrong).to(device),
        "ligand": torch.from_numpy(ligand).to(device),
        "y": torch.from_numpy(y).to(device),
    }
    tasks = defaultdict(dict)
    by_split_target = defaultdict(lambda: defaultdict(list))
    for index, row in enumerate(cells):
        by_split_target[row["split"]][row["target_id"]].append(index)
    for split, split_tasks in by_split_target.items():
        tasks[split] = {
            target: np.asarray(indices, dtype=np.int64)
            for target, indices in sorted(split_tasks.items()) if len(indices) >= 8
        }
    metadata = {
        "y_mean": y_mean, "y_scale": y_scale,
        "full_mean": full_norm["mean"], "full_scale": full_norm["scale"],
        "ligand_mean": ligand_norm["mean"], "ligand_scale": ligand_norm["scale"],
    }
    return cells, tensors, tasks, metadata


def draw_episode(indices: np.ndarray, rng: np.random.Generator, k: int,
                 max_query: int | None = None):
    order = rng.permutation(indices)
    support, query = order[:k], order[k:]
    if max_query is not None and len(query) > max_query:
        query = rng.choice(query, size=max_query, replace=False)
    return support, np.asarray(query)


def validation_mse(model, family_name, tensors, tasks, seed, draws, y_scale):
    rng = np.random.default_rng(seed)
    losses = []
    model.eval()
    with torch.no_grad():
        for target, indices in tasks.items():
            target_losses = []
            for _ in range(draws):
                support, query = draw_episode(indices, rng, 5)
                prediction = model.episode(
                    tensors["ligand"][support], tensors[family_name][support], tensors["y"][support],
                    tensors["ligand"][query], tensors[family_name][query],
                )
                target_losses.append(float(((prediction - tensors["y"][query]) * y_scale).square().mean()))
            losses.append(float(np.mean(target_losses)))
    return float(np.mean(losses))


def train_model(tensors, train_tasks, val_tasks, family_name, section_dim, ridge,
                seed, steps, config: TrainConfig):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = MetaSectionRegressor(288, section_dim, ridge).to(tensors["y"].device)
    optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate,
                                 weight_decay=config.weight_decay)
    targets = sorted(train_tasks)
    best_state, best_val, best_step = None, float("inf"), 0
    losses = []
    for step in range(1, steps + 1):
        model.train()
        optimizer.zero_grad()
        loss = torch.zeros((), device=tensors["y"].device)
        selected = rng.choice(len(targets), size=config.tasks_per_step, replace=True)
        for target_index in selected:
            indices = train_tasks[targets[int(target_index)]]
            support, query = draw_episode(indices, rng, config.k, config.query_per_task)
            prediction = model.episode(
                tensors["ligand"][support], tensors[family_name][support], tensors["y"][support],
                tensors["ligand"][query], tensors[family_name][query],
            )
            loss = loss + (prediction - tensors["y"][query]).square().mean()
        loss = loss / config.tasks_per_step
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
        optimizer.step()
        losses.append(float(loss.detach()))
        if step % config.val_interval == 0 or step == steps:
            value = validation_mse(model, family_name, tensors, val_tasks,
                                   seed + step, config.val_draws, 1.0)
            if value < best_val:
                best_val, best_step = value, step
                best_state = copy.deepcopy(model.state_dict())
    model.load_state_dict(best_state)
    return model, {"best_standardized_val_mse": best_val, "best_step": best_step,
                   "initial_loss": float(np.mean(losses[:50])),
                   "final_loss": float(np.mean(losses[-50:]))}


def save_checkpoint(model, path: Path, payload: dict) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    torch.save({"model_state": model.state_dict(), **payload}, path)
    return {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}


def predict_arm(model, family_name, tensors, tasks, cells, seed, draws, arm,
                support_mode="correct", query_family=None):
    targets = sorted(tasks)
    donor = {target: targets[(index + 1) % len(targets)] for index, target in enumerate(targets)}
    rows = []
    model.eval()
    with torch.no_grad():
        for target in targets:
            for draw in range(draws):
                episode_key = f"main-v0-episode|{seed}|{target}|{draw}"
                episode_seed = int(hashlib.sha256(episode_key.encode()).hexdigest()[:16], 16)
                support, query = draw_episode(tasks[target], np.random.default_rng(episode_seed), 5)
                source_support = support
                source_family = family_name
                source_y = tensors["y"][support]
                if support_mode == "foreign":
                    donor_key = f"main-v0-donor|{seed}|{target}|{draw}|{donor[target]}"
                    donor_seed = int(hashlib.sha256(donor_key.encode()).hexdigest()[:16], 16)
                    source_support, _ = draw_episode(
                        tasks[donor[target]], np.random.default_rng(donor_seed), 5)
                    source_y = tensors["y"][source_support]
                elif support_mode == "permuted":
                    source_y = source_y.roll(1)
                query_feature_name = query_family or family_name
                if support_mode == "zero":
                    prediction, _ = model.components(
                        tensors["ligand"][query], tensors[query_feature_name][query])
                else:
                    prediction = model.episode(
                        tensors["ligand"][source_support], tensors[source_family][source_support], source_y,
                        tensors["ligand"][query], tensors[query_feature_name][query],
                    )
                prediction = prediction.detach().cpu().numpy()
                support_hash = stable_cell_hash(cells, source_support)
                query_hash = stable_cell_hash(cells, query)
                for index, value in zip(query, prediction):
                    rows.append({"arm": arm, "seed": seed, "target_id": target,
                                 "draw": draw, "cell_id": cells[int(index)]["cell_id"],
                                 "support_hash": support_hash, "query_hash": query_hash,
                                 "prediction_standardized": float(value)})
    return rows


def stable_cell_hash(cells, indices):
    value = "|".join(sorted(cells[int(index)]["cell_id"] for index in indices))
    return hashlib.sha256(value.encode()).hexdigest()


def write_prediction_rows(path: Path, rows: list[dict]) -> str:
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


def concordance(y, prediction):
    concordant, comparable = 0.0, 0
    for index in range(len(y) - 1):
        difference = y[index + 1:] - y[index]
        mask = difference != 0
        if not mask.any():
            continue
        pred_difference = prediction[index + 1:] - prediction[index]
        product = difference[mask] * pred_difference[mask]
        concordant += float((product > 0).sum()) + 0.5 * float((product == 0).sum())
        comparable += int(mask.sum())
    return concordant / comparable if comparable else float("nan")


def score_predictions(rows, cells, y_mean, y_scale):
    truth = {row["cell_id"]: float(row["pK"]) for row in cells}
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["arm"], row["seed"], row["target_id"], row["draw"])].append(row)
    draw_metrics, losses = [], defaultdict(list)
    for (arm, seed, target, draw), values in grouped.items():
        y = np.asarray([truth[row["cell_id"]] for row in values])
        prediction = np.asarray([row["prediction_standardized"] * y_scale + y_mean for row in values])
        squared = (y - prediction) ** 2
        variance = ((y - y.mean()) ** 2).sum()
        pearson = stats.pearsonr(y, prediction).statistic if len(y) > 1 and y.std() > 0 and prediction.std() > 0 else np.nan
        spearman = stats.spearmanr(y, prediction).statistic if len(y) > 1 and y.std() > 0 and prediction.std() > 0 else np.nan
        metric = {
            "arm": arm, "seed": seed, "target_id": target, "draw": draw,
            "mse": float(squared.mean()), "rmse": float(np.sqrt(squared.mean())),
            "r2": float(1.0 - squared.sum() / variance) if variance > 0 else np.nan,
            "ci": concordance(y, prediction), "pearson": float(pearson),
            "spearman": float(spearman),
        }
        draw_metrics.append(metric)
        losses[(arm, target)].append(metric["mse"])
    target_loss = {key: float(np.mean(value)) for key, value in losses.items()}
    summary = {}
    for arm in sorted({row["arm"] for row in draw_metrics}):
        target_metrics = defaultdict(list)
        for metric in draw_metrics:
            if metric["arm"] == arm:
                for name in ("mse", "rmse", "r2", "ci", "pearson", "spearman"):
                    target_metrics[(metric["target_id"], name)].append(metric[name])
        summary[arm] = {
            name: float(np.nanmean([np.nanmean(values) for (target, metric_name), values
                                    in target_metrics.items() if metric_name == name]))
            for name in ("mse", "rmse", "r2", "ci", "pearson", "spearman")
        }
    return summary, target_loss, draw_metrics


def bootstrap_contrast(target_loss, correct, control, draws, seed):
    targets = sorted({target for arm, target in target_loss if arm == correct}
                     & {target for arm, target in target_loss if arm == control})
    delta = np.asarray([target_loss[(control, target)] - target_loss[(correct, target)]
                        for target in targets])
    rng = np.random.default_rng(seed)
    boot = delta[rng.integers(0, len(delta), size=(draws, len(delta)))].mean(axis=1)
    return {"correct": correct, "control": control, "targets": len(targets),
            "mean_mse_reduction": float(delta.mean()),
            "one_sided_95_lcb": float(np.quantile(boot, 0.05)),
            "pass": bool(np.quantile(boot, 0.05) > 0)}


def cluster_bootstrap_contrast(target_loss, target_cluster, correct, control, draws, seed):
    targets = sorted({target for arm, target in target_loss if arm == correct}
                     & {target for arm, target in target_loss if arm == control})
    by_cluster = defaultdict(list)
    for target in targets:
        by_cluster[target_cluster[target]].append(
            target_loss[(control, target)] - target_loss[(correct, target)])
    cluster_delta = np.asarray([np.mean(values) for _, values in sorted(by_cluster.items())])
    rng = np.random.default_rng(seed)
    boot = cluster_delta[
        rng.integers(0, len(cluster_delta), size=(draws, len(cluster_delta)))
    ].mean(axis=1)
    return {"correct": correct, "control": control, "clusters": len(cluster_delta),
            "cluster_macro_mse_reduction": float(cluster_delta.mean()),
            "one_sided_95_lcb": float(np.quantile(boot, 0.05)),
            "sensitivity_pass": bool(np.quantile(boot, 0.05) > 0)}


def run(config: TrainConfig = TrainConfig(), device: str = "cuda") -> dict:
    if device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    OUT.mkdir(parents=True, exist_ok=True)
    cells, tensors, tasks, normalization = load_data(device)
    grid = []
    for section_dim in range(1, 6):
        for ridge in (0.01, 0.1, 1.0):
            model, diagnostics = train_model(
                tensors, tasks["meta_train"], tasks["meta_val"], "correct",
                section_dim, ridge, SEEDS[0], config.selection_steps, config,
            )
            value = validation_mse(model, "correct", tensors, tasks["meta_val"],
                                   SEEDS[0] + 991, 5, normalization["y_scale"])
            grid.append({"d": section_dim, "ridge": ridge,
                         "validation_target_macro_mse": value, **diagnostics})
    selected = min(grid, key=lambda row: (row["validation_target_macro_mse"],
                                          row["d"], -row["ridge"]))

    checkpoints, all_predictions, train_diagnostics = [], [], []
    for seed in SEEDS:
        models = {}
        for name, family, d in (("population", "correct", 0),
                                ("full", "correct", selected["d"]),
                                ("ligand_only", "ligand", selected["d"])):
            model, diagnostics = train_model(
                tensors, tasks["meta_train"], tasks["meta_val"], family,
                d, selected["ridge"], seed, config.final_steps, config,
            )
            models[name] = model
            checkpoints.append(save_checkpoint(
                model, OUT / "checkpoints" / f"{name}_seed{seed}.pt",
                {"seed": seed, "arm": name, "d": d, "ridge": selected["ridge"],
                 "corpus_sha256": sha256(CORPUS / "manifest.json"),
                 "features_sha256": sha256(FEATURES)},
            ))
            train_diagnostics.append({"seed": seed, "arm": name, **diagnostics})

        all_predictions += predict_arm(models["population"], "correct", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws, "population_d0")
        all_predictions += predict_arm(models["ligand_only"], "ligand", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws, "ligand_only")
        all_predictions += predict_arm(models["full"], "correct", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws, "full_correct")
        all_predictions += predict_arm(models["full"], "correct", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws,
                                       "full_zero", support_mode="zero")
        all_predictions += predict_arm(models["full"], "correct", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws,
                                       "full_foreign", support_mode="foreign")
        all_predictions += predict_arm(models["full"], "correct", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws,
                                       "full_permuted", support_mode="permuted")
        all_predictions += predict_arm(models["full"], "wrong", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws,
                                       "full_wrong_protein")

    prediction_path = OUT / "predictions_before_query_labels.jsonl.gz"
    prediction_hash = write_prediction_rows(prediction_path, all_predictions)
    summary, target_loss, draw_metrics = score_predictions(
        all_predictions, cells, normalization["y_mean"], normalization["y_scale"])
    controls = ("population_d0", "full_zero", "full_foreign", "full_permuted",
                "ligand_only", "full_wrong_protein")
    contrasts = [bootstrap_contrast(target_loss, "full_correct", control,
                                    config.bootstrap_draws, SEEDS[0] + index)
                 for index, control in enumerate(controls)]
    target_cluster = {row["target_id"]: row["protein_group_40"] for row in cells}
    cluster_sensitivity = [cluster_bootstrap_contrast(
        target_loss, target_cluster, "full_correct", control,
        config.bootstrap_draws, SEEDS[0] + 100 + index)
        for index, control in enumerate(controls)]
    gates = {
        "M1_meta_effect": next(row["pass"] for row in contrasts if row["control"] == "population_d0"),
        "M2_support_specificity": all(row["pass"] for row in contrasts
                                      if row["control"] in {"full_zero", "full_foreign", "full_permuted"}),
        "M3_biological_specificity": all(row["pass"] for row in contrasts
                                         if row["control"] in {"ligand_only", "full_wrong_protein"}),
    }
    if all(gates.values()):
        verdict = "REAL_BIOLOGICAL_META_SECTION_V0_PASS"
    elif not gates["M1_meta_effect"]:
        verdict = "META_EFFECT_NOT_IDENTIFIED"
    elif not gates["M2_support_specificity"]:
        verdict = "SUPPORT_SPECIFICITY_NOT_IDENTIFIED"
    else:
        verdict = "BIOLOGICAL_SPECIFICITY_NOT_IDENTIFIED"
    result = {
        "schema": "MetaSieve.MainV0Result.v1", "config": asdict(config),
        "environment": {"python": subprocess.check_output(["python", "--version"], text=True).strip(),
                        "torch": torch.__version__, "device": device,
                        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()},
        "data": {"corpus_sha256": sha256(CORPUS / "manifest.json"),
                 "features_sha256": sha256(FEATURES),
                 "train_tasks_k5": len(tasks["meta_train"]),
                 "val_tasks_k5": len(tasks["meta_val"]),
                 "test_tasks_k5": len(tasks["meta_test"])},
        "selection_grid": grid, "selected": {"d": selected["d"], "ridge": selected["ridge"]},
        "train_diagnostics": train_diagnostics, "checkpoints": checkpoints,
        "prediction_artifact": {"path": str(prediction_path.relative_to(ROOT)),
                                "sha256_before_scoring": prediction_hash,
                                "query_labels_in_artifact": False,
                                "rows": len(all_predictions)},
        "point_metrics_target_macro": summary, "paired_mse_reductions": contrasts,
        "cdhit_cluster_bootstrap_sensitivity": cluster_sensitivity,
        "gates": gates, "law_metrics": "NA_NOT_ADMITTED",
        "frozen_strict_confirmation_opened": False,
        "TERMINAL_VERDICT": verdict,
    }
    result = json_safe(result)
    (OUT / "MAIN_V0_RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "draw_metrics.json").write_text(
        json.dumps(json_safe(draw_metrics), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    return result


def evaluate_existing(device: str = "cuda") -> dict:
    result_path = OUT / "MAIN_V0_RESULT.json"
    existing = json.loads(result_path.read_text(encoding="utf-8"))
    config = TrainConfig(**existing["config"])
    selected = existing["selected"]
    cells, tensors, tasks, normalization = load_data(device)
    all_predictions = []
    for seed in SEEDS:
        models = {}
        for name, family, d in (("population", "correct", 0),
                                ("full", "correct", selected["d"]),
                                ("ligand_only", "ligand", selected["d"])):
            model = MetaSectionRegressor(288, d, selected["ridge"]).to(device)
            checkpoint = torch.load(
                OUT / "checkpoints" / f"{name}_seed{seed}.pt",
                map_location=device, weights_only=False)
            model.load_state_dict(checkpoint["model_state"])
            models[name] = model
        all_predictions += predict_arm(models["population"], "correct", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws, "population_d0")
        all_predictions += predict_arm(models["ligand_only"], "ligand", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws, "ligand_only")
        all_predictions += predict_arm(models["full"], "correct", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws, "full_correct")
        all_predictions += predict_arm(models["full"], "correct", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws,
                                       "full_zero", support_mode="zero")
        all_predictions += predict_arm(models["full"], "correct", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws,
                                       "full_foreign", support_mode="foreign")
        all_predictions += predict_arm(models["full"], "correct", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws,
                                       "full_permuted", support_mode="permuted")
        all_predictions += predict_arm(models["full"], "wrong", tensors,
                                       tasks["meta_test"], cells, seed, config.test_draws,
                                       "full_wrong_protein")

    prediction_path = OUT / "predictions_before_query_labels.jsonl.gz"
    prediction_hash = write_prediction_rows(prediction_path, all_predictions)
    summary, target_loss, draw_metrics = score_predictions(
        all_predictions, cells, normalization["y_mean"], normalization["y_scale"])
    controls = ("population_d0", "full_zero", "full_foreign", "full_permuted",
                "ligand_only", "full_wrong_protein")
    contrasts = [bootstrap_contrast(target_loss, "full_correct", control,
                                    config.bootstrap_draws, SEEDS[0] + index)
                 for index, control in enumerate(controls)]
    target_cluster = {row["target_id"]: row["protein_group_40"] for row in cells}
    cluster_sensitivity = [cluster_bootstrap_contrast(
        target_loss, target_cluster, "full_correct", control,
        config.bootstrap_draws, SEEDS[0] + 100 + index)
        for index, control in enumerate(controls)]
    gates = {
        "M1_meta_effect": next(row["pass"] for row in contrasts if row["control"] == "population_d0"),
        "M2_support_specificity": all(row["pass"] for row in contrasts
                                      if row["control"] in {"full_zero", "full_foreign", "full_permuted"}),
        "M3_biological_specificity": all(row["pass"] for row in contrasts
                                         if row["control"] in {"ligand_only", "full_wrong_protein"}),
    }
    if all(gates.values()):
        verdict = "REAL_BIOLOGICAL_META_SECTION_V0_PASS"
    elif not gates["M1_meta_effect"]:
        verdict = "META_EFFECT_NOT_IDENTIFIED"
    elif not gates["M2_support_specificity"]:
        verdict = "SUPPORT_SPECIFICITY_NOT_IDENTIFIED"
    else:
        verdict = "BIOLOGICAL_SPECIFICITY_NOT_IDENTIFIED"
    existing.update({
        "evaluation_protocol": "paired_episode_v2_deterministic_per_target_draw",
        "prediction_artifact": {"path": str(prediction_path.relative_to(ROOT)),
                                "sha256_before_scoring": prediction_hash,
                                "query_labels_in_artifact": False,
                                "rows": len(all_predictions)},
        "point_metrics_target_macro": summary,
        "paired_mse_reductions": contrasts,
        "cdhit_cluster_bootstrap_sensitivity": cluster_sensitivity,
        "gates": gates,
        "registered_target_level_verdict": verdict,
        "scientific_admission_verdict": (
            "REAL_BIOLOGICAL_META_SECTION_V0_PASS"
            if all(row["sensitivity_pass"] for row in cluster_sensitivity)
            else "BIOLOGICAL_SPECIFICITY_NOT_IDENTIFIED_CLUSTER_SENSITIVITY"
        ),
        "production_migration_authorized": all(
            row["sensitivity_pass"] for row in cluster_sensitivity),
        "TERMINAL_VERDICT": (
            "REAL_BIOLOGICAL_META_SECTION_V0_PASS"
            if all(row["sensitivity_pass"] for row in cluster_sensitivity)
            else "BIOLOGICAL_SPECIFICITY_NOT_IDENTIFIED_CLUSTER_SENSITIVITY"
        ),
    })
    existing = json_safe(existing)
    result_path.write_text(json.dumps(existing, indent=2, sort_keys=True,
                                      allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "draw_metrics.json").write_text(
        json.dumps(json_safe(draw_metrics), indent=2, sort_keys=True,
                   allow_nan=False) + "\n", encoding="utf-8")
    return existing


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--evaluate-existing", action="store_true")
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    result = evaluate_existing(args.device) if args.evaluate_existing else run(device=args.device)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
