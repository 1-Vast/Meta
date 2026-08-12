"""Source-only MetaSieve v1 development experiment with sealed meta-validation."""
from __future__ import annotations

import copy
import gzip
import hashlib
import io
import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import torch
from scipy import stats
from torch import nn

from research.meta_fewshot.train_main_v0 import (
    MetaSectionRegressor,
    cluster_bootstrap_contrast,
    concordance,
    json_safe,
    sha256,
)

ROOT = Path(__file__).resolve().parents[2]
DEV = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_v1_development"
OUT = ROOT / "report/meta_fewshot/v1_development"
SEEDS = (20260821, 20260822, 20260823, 20260824, 20260825)


@dataclass(frozen=True)
class V1Config:
    k: int = 5
    section_dim: int = 2
    ridge: float = 1.0
    bottleneck: int = 32
    tasks_per_step: int = 8
    query_per_task: int = 32
    steps: int = 1000
    learning_rate: float = 0.01
    weight_decay: float = 1e-4
    auxiliary_weight: float = 0.1
    auxiliary_groups_per_step: int = 8
    bootstrap_draws: int = 9999


class PairPriorMetaSection(nn.Module):
    def __init__(self, input_dim: int, section_dim: int, ridge: float, bottleneck: int = 32):
        super().__init__()
        if not 0 <= section_dim <= 5 or ridge <= 0 or bottleneck <= 0:
            raise ValueError("invalid pair-prior Meta-Section dimensions")
        self.section_dim = section_dim
        self.ridge = float(ridge)
        self.bottleneck = bottleneck
        self.ligand_population = nn.Linear(input_dim, 1)
        self.adapter_down = nn.Linear(input_dim, bottleneck, bias=False)
        self.adapter_up = nn.Linear(bottleneck, input_dim, bias=False)
        nn.init.zeros_(self.adapter_up.weight)
        self.pair_population = nn.Linear(input_dim, 1, bias=False)
        if section_dim:
            self.raw_basis = nn.Parameter(torch.randn(input_dim, section_dim) / math.sqrt(input_dim))

    def encoded_pair(self, pair):
        return pair + self.adapter_up(torch.nn.functional.silu(self.adapter_down(pair))) / math.sqrt(self.bottleneck)

    def basis(self):
        return torch.linalg.qr(self.raw_basis, mode="reduced").Q

    def components(self, ligand, pair):
        encoded = self.encoded_pair(pair)
        population = self.ligand_population(ligand).squeeze(-1) + self.pair_population(encoded).squeeze(-1)
        coordinates = encoded @ self.basis() if self.section_dim else None
        return population, coordinates

    def episode(self, support_ligand, support_pair, support_y, query_ligand, query_pair):
        support_population, support_coordinates = self.components(support_ligand, support_pair)
        query_population, query_coordinates = self.components(query_ligand, query_pair)
        if not self.section_dim:
            return query_population
        residual = support_y - support_population
        identity = torch.eye(len(support_y), device=support_y.device, dtype=support_y.dtype)
        dual = torch.linalg.solve(
            support_coordinates @ support_coordinates.T + self.ridge * identity, residual)
        coefficient = support_coordinates.T @ dual
        return query_population + query_coordinates @ coefficient


def read_jsonl_gz(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def standardize(values: np.ndarray, mean=None, scale=None):
    if mean is None:
        mean = values.mean(axis=0)
    if scale is None:
        scale = values.std(axis=0)
    scale = np.asarray(scale).copy()
    scale[scale < 1e-6] = 1.0
    return ((values - mean) / scale).astype(np.float32), mean, scale


def load_sealed(device: str):
    manifest = json.loads((DEV / "manifest.json").read_text())
    if manifest["main_v0_test_values_used"] != 0:
        raise ValueError("development seal used main-v0 test values")
    source = read_jsonl_gz(DEV / "source_cells.jsonl.gz")
    validation = read_jsonl_gz(DEV / "metaval_cells_without_labels.jsonl.gz")
    if any("pK" in row for row in validation):
        raise ValueError("meta-validation public rows contain query labels")
    blacklist = json.loads((DEV / "main_v0_test_blacklist.json").read_text())
    if ({row["target_id"] for row in source + validation} & set(blacklist["targets"]) or
            {row["cell_id"] for row in source + validation} & set(blacklist["cells"])):
        raise ValueError("development loader crossed the main-v0 test blacklist")
    with np.load(DEV / "source_features.npz", allow_pickle=False) as stored:
        if stored["cell_id"].tolist() != [row["cell_id"] for row in source]:
            raise ValueError("source feature order mismatch")
        source_pair = stored["correct"].astype(np.float32)
        source_ligand = stored["ligand"].astype(np.float32)
    with np.load(DEV / "metaval_features.npz", allow_pickle=False) as stored:
        if stored["cell_id"].tolist() != [row["cell_id"] for row in validation]:
            raise ValueError("meta-validation feature order mismatch")
        val_pair = stored["correct"].astype(np.float32)
        val_ligand = stored["ligand"].astype(np.float32)
    with np.load(DEV / "metaval_wrong_features.npz", allow_pickle=False) as stored:
        if stored["cell_id"].tolist() != [row["cell_id"] for row in validation]:
            raise ValueError("wrong feature order mismatch")
        val_wrong = stored["wrong"].astype(np.float32)
    source_pair, pair_mean, pair_scale = standardize(source_pair)
    val_pair, _, _ = standardize(val_pair, pair_mean, pair_scale)
    val_wrong, _, _ = standardize(val_wrong, pair_mean, pair_scale)
    source_ligand, ligand_mean, ligand_scale = standardize(source_ligand)
    val_ligand, _, _ = standardize(val_ligand, ligand_mean, ligand_scale)
    source_y_raw = np.asarray([row["pK"] for row in source], dtype=np.float32)
    y_mean, y_scale = float(source_y_raw.mean()), float(source_y_raw.std())
    source_y = ((source_y_raw - y_mean) / y_scale).astype(np.float32)
    tensors = {
        "source_pair": torch.from_numpy(source_pair).to(device),
        "source_ligand": torch.from_numpy(source_ligand).to(device),
        "source_y": torch.from_numpy(source_y).to(device),
        "val_pair": torch.from_numpy(val_pair).to(device),
        "val_wrong": torch.from_numpy(val_wrong).to(device),
        "val_ligand": torch.from_numpy(val_ligand).to(device),
    }
    source_index = {row["cell_id"]: index for index, row in enumerate(source)}
    val_index = {row["cell_id"]: index for index, row in enumerate(validation)}
    tasks = defaultdict(list)
    target_cluster = {}
    for index, row in enumerate(source):
        tasks[row["target_id"]].append(index)
        target_cluster[row["target_id"]] = row["protein_group_40"]
    tasks = {target: np.asarray(indices, dtype=np.int64)
             for target, indices in tasks.items() if len(indices) >= 8}
    clusters = defaultdict(list)
    for target in tasks:
        clusters[target_cluster[target]].append(target)
    contrast_groups = read_jsonl_gz(DEV / "source_contrast_groups.jsonl.gz")
    contrasts = defaultdict(list)
    for group in contrast_groups:
        members = [{"index": source_index[row["cell_id"]],
                    "y": (float(row["pK"]) - y_mean) / y_scale,
                    "protein_group_40": row["protein_group_40"]} for row in group["members"]]
        contrasts[group["kind"]].append(members)
    return source, validation, tensors, tasks, dict(clusters), contrasts, val_index, y_mean, y_scale


def draw_episode(indices, rng, k, max_query):
    order = rng.permutation(indices)
    support, query = order[:k], order[k:]
    if len(query) > max_query:
        query = rng.choice(query, size=max_query, replace=False)
    return support, np.asarray(query)


def sample_cluster_targets(clusters, rng, count):
    cluster_keys = sorted(clusters)
    chosen = rng.choice(len(cluster_keys), size=count, replace=True)
    return [clusters[cluster_keys[int(index)]][int(rng.integers(0, len(clusters[cluster_keys[int(index)]])))]
            for index in chosen]


def auxiliary_loss(model, tensors, groups, rng, count):
    total = torch.zeros((), device=tensors["source_y"].device)
    for _ in range(count):
        members = groups[int(rng.integers(0, len(groups)))]
        left_index, right_index = rng.choice(len(members), size=2, replace=False)
        left, right = members[int(left_index)], members[int(right_index)]
        if "protein_group_40" in left and len({row["protein_group_40"] for row in members}) > 1:
            eligible = [row for row in members
                        if row["protein_group_40"] != left["protein_group_40"]]
            right = eligible[int(rng.integers(0, len(eligible)))]
        indices = torch.as_tensor([left["index"], right["index"]], device=tensors["source_y"].device)
        population, _ = model.components(tensors["source_ligand"][indices], tensors["source_pair"][indices])
        observed = torch.tensor(left["y"] - right["y"], device=population.device, dtype=population.dtype)
        total = total + ((population[0] - population[1]) - observed).square()
    return total / count


def train(model, tensors, tasks, clusters, contrasts, seed, config, auxiliaries=False):
    optimizer = torch.optim.Adam(model.parameters(), lr=config.learning_rate,
                                 weight_decay=config.weight_decay)
    episode_rng = np.random.default_rng(seed)
    auxiliary_rng = np.random.default_rng(seed + 100000)
    losses = []
    for _ in range(config.steps):
        optimizer.zero_grad()
        loss = torch.zeros((), device=tensors["source_y"].device)
        for target in sample_cluster_targets(clusters, episode_rng, config.tasks_per_step):
            support, query = draw_episode(tasks[target], episode_rng, config.k, config.query_per_task)
            prediction = model.episode(
                tensors["source_ligand"][support], tensors["source_pair"][support],
                tensors["source_y"][support], tensors["source_ligand"][query],
                tensors["source_pair"][query])
            loss = loss + (prediction - tensors["source_y"][query]).square().mean()
        loss = loss / config.tasks_per_step
        if auxiliaries:
            within = auxiliary_loss(model, tensors, contrasts["within_panel_ligand"],
                                    auxiliary_rng, config.auxiliary_groups_per_step)
            partner = auxiliary_loss(model, tensors, contrasts["measured_partner"],
                                     auxiliary_rng, config.auxiliary_groups_per_step)
            loss = loss + config.auxiliary_weight * (within + partner)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
        optimizer.step()
        losses.append(float(loss.detach()))
    return {"initial_loss": float(np.mean(losses[:50])),
            "final_loss": float(np.mean(losses[-50:]))}


def predict(model, episodes, tensors, val_index, seed, arm, support_pair="correct",
            query_pair="correct", support_mode="correct", y_mean=0.0, y_scale=1.0):
    pair_tensor = {"correct": tensors["val_pair"], "wrong": tensors["val_wrong"]}
    rows = []
    model.eval()
    with torch.no_grad():
        for episode in episodes:
            if episode["seed"] != seed:
                continue
            support = np.asarray([val_index[cell_id] for cell_id in episode["support_cell_ids"]])
            query = np.asarray([val_index[cell_id] for cell_id in episode["query_cell_ids"]])
            support_y = torch.tensor(
                [(value - y_mean) / y_scale for value in episode["support_pK"]],
                device=tensors["val_pair"].device)
            if support_mode == "permuted":
                support_y = support_y.roll(1)
            if support_mode == "zero":
                prediction, _ = model.components(tensors["val_ligand"][query], pair_tensor[query_pair][query])
            else:
                prediction = model.episode(
                    tensors["val_ligand"][support], pair_tensor[support_pair][support], support_y,
                    tensors["val_ligand"][query], pair_tensor[query_pair][query])
            for cell_id, value in zip(episode["query_cell_ids"], prediction.cpu().numpy()):
                rows.append({"arm": arm, "seed": seed, "target_id": episode["target_id"],
                             "protein_group_40": episode["protein_group_40"],
                             "draw": episode["draw"], "cell_id": cell_id,
                             "support_hash": episode["support_hash"],
                             "query_hash": episode["query_hash"],
                             "prediction_standardized": float(value)})
    return rows


def predict_foreign(model, episodes, tensors, val_index, seed, arm, y_mean, y_scale):
    active = [row for row in episodes if row["seed"] == seed]
    targets = sorted({row["target_id"] for row in active})
    donor = {target: targets[(index + 1) % len(targets)] for index, target in enumerate(targets)}
    lookup = {(row["target_id"], row["draw"]): row for row in active}
    rows = []
    model.eval()
    with torch.no_grad():
        for episode in active:
            source = lookup[(donor[episode["target_id"]], episode["draw"])]
            support = np.asarray([val_index[cell_id] for cell_id in source["support_cell_ids"]])
            query = np.asarray([val_index[cell_id] for cell_id in episode["query_cell_ids"]])
            support_y = torch.tensor(
                [(value - y_mean) / y_scale for value in source["support_pK"]],
                device=tensors["val_pair"].device)
            prediction = model.episode(
                tensors["val_ligand"][support], tensors["val_pair"][support], support_y,
                tensors["val_ligand"][query], tensors["val_pair"][query])
            for cell_id, value in zip(episode["query_cell_ids"], prediction.cpu().numpy()):
                rows.append({"arm": arm, "seed": seed, "target_id": episode["target_id"],
                             "protein_group_40": episode["protein_group_40"],
                             "draw": episode["draw"], "cell_id": cell_id,
                             "support_hash": source["support_hash"],
                             "query_hash": episode["query_hash"],
                             "prediction_standardized": float(value)})
    return rows


def write_predictions(path, rows):
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


def score(rows, truth_path, y_mean, y_scale):
    truth = {row["cell_id"]: float(row["pK"]) for row in read_jsonl_gz(truth_path)}
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["arm"], row["seed"], row["target_id"], row["draw"])].append(row)
    metrics, losses = [], defaultdict(list)
    for (arm, seed, target, draw), values in grouped.items():
        y = np.asarray([truth[row["cell_id"]] for row in values])
        prediction = np.asarray([row["prediction_standardized"] * y_scale + y_mean for row in values])
        squared = (y - prediction) ** 2
        variance = ((y - y.mean()) ** 2).sum()
        metrics.append({
            "arm": arm, "seed": seed, "target_id": target, "draw": draw,
            "protein_group_40": values[0]["protein_group_40"],
            "mse": float(squared.mean()), "rmse": float(np.sqrt(squared.mean())),
            "r2": float(1 - squared.sum() / variance) if variance > 0 else np.nan,
            "ci": concordance(y, prediction),
            "pearson": float(stats.pearsonr(y, prediction).statistic)
            if y.std() > 0 and prediction.std() > 0 else np.nan,
            "spearman": float(stats.spearmanr(y, prediction).statistic)
            if y.std() > 0 and prediction.std() > 0 else np.nan,
        })
        losses[(arm, target)].append(float(squared.mean()))
    target_loss = {key: float(np.mean(value)) for key, value in losses.items()}
    summary = {}
    for arm in sorted({row["arm"] for row in metrics}):
        summary[arm] = {}
        for name in ("mse", "rmse", "r2", "ci", "pearson", "spearman"):
            by_target = defaultdict(list)
            for row in metrics:
                if row["arm"] == arm:
                    by_target[row["target_id"]].append(row[name])
            summary[arm][name] = float(np.nanmean([np.nanmean(value) for value in by_target.values()]))
    target_cluster = {row["target_id"]: row["protein_group_40"] for row in metrics}
    return summary, target_loss, target_cluster, metrics


def run(config=V1Config(), device="cuda"):
    if device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    OUT.mkdir(parents=True, exist_ok=True)
    _, _, tensors, tasks, clusters, contrasts, val_index, y_mean, y_scale = load_sealed(device)
    episodes = read_jsonl_gz(DEV / "metaval_episodes.jsonl.gz")
    all_rows, diagnostics, checkpoints = [], [], []
    for seed in SEEDS:
        models = {}
        specs = (("ligand_d0", "legacy", 0, False),
                 ("v0", "legacy", config.section_dim, False),
                 ("pair_d0", "pair", 0, False),
                 ("v1a", "pair", config.section_dim, False),
                 ("v1b", "pair", config.section_dim, True))
        for name, kind, dimension, auxiliaries in specs:
            torch.manual_seed(seed)
            model = (MetaSectionRegressor(288, dimension, config.ridge) if kind == "legacy"
                     else PairPriorMetaSection(288, dimension, config.ridge, config.bottleneck)).to(device)
            path = OUT / "checkpoints" / f"{name}_seed{seed}.pt"
            path.parent.mkdir(parents=True, exist_ok=True)
            seal_hash = sha256(DEV / "manifest.json")
            if path.exists():
                checkpoint = torch.load(path, map_location=device, weights_only=False)
                if (checkpoint.get("seed") != seed or checkpoint.get("arm") != name
                        or checkpoint.get("config") != asdict(config)
                        or checkpoint.get("development_seal") != seal_hash):
                    raise ValueError(f"incompatible resume checkpoint {path}")
                model.load_state_dict(checkpoint["model_state"])
                detail = {"resumed_checkpoint": True}
            else:
                detail = train(model, tensors, tasks, clusters, contrasts, seed, config, auxiliaries)
                torch.save({"model_state": model.state_dict(), "seed": seed, "arm": name,
                            "config": asdict(config), "development_seal": seal_hash,
                            "train_diagnostics": detail}, path)
            models[name] = model
            checkpoints.append({"path": str(path.relative_to(ROOT)), "sha256": sha256(path)})
            diagnostics.append({"seed": seed, "arm": name, **detail})

        all_rows += predict(models["ligand_d0"], episodes, tensors, val_index, seed,
                            "ligand_d0", y_mean=y_mean, y_scale=y_scale)
        all_rows += predict(models["pair_d0"], episodes, tensors, val_index, seed,
                            "pair_d0", y_mean=y_mean, y_scale=y_scale)
        for name in ("v0", "v1a", "v1b"):
            model = models[name]
            common = {"y_mean": y_mean, "y_scale": y_scale}
            all_rows += predict(model, episodes, tensors, val_index, seed, f"{name}_correct", **common)
            all_rows += predict(model, episodes, tensors, val_index, seed, f"{name}_zero",
                                support_mode="zero", **common)
            all_rows += predict(model, episodes, tensors, val_index, seed, f"{name}_permuted",
                                support_mode="permuted", **common)
            all_rows += predict(model, episodes, tensors, val_index, seed, f"{name}_sc_qw",
                                query_pair="wrong", **common)
            all_rows += predict(model, episodes, tensors, val_index, seed, f"{name}_sw_qc",
                                support_pair="wrong", **common)
            all_rows += predict(model, episodes, tensors, val_index, seed, f"{name}_sw_qw",
                                support_pair="wrong", query_pair="wrong", **common)
            all_rows += predict_foreign(model, episodes, tensors, val_index, seed,
                                        f"{name}_foreign", y_mean, y_scale)

    prediction_path = OUT / "predictions_before_query_labels.jsonl.gz"
    prediction_hash = write_predictions(prediction_path, all_rows)
    summary, target_loss, target_cluster, draw_metrics = score(
        all_rows, DEV / "metaval_query_truth.jsonl.gz", y_mean, y_scale)
    comparisons = [
        ("pair_d0", "ligand_d0"),
        ("v0_correct", "v0_zero"), ("v0_correct", "v0_permuted"),
        ("v0_correct", "v0_sc_qw"), ("v0_correct", "v0_sw_qc"),
        ("v0_correct", "v0_sw_qw"), ("v0_correct", "v0_foreign"),
        ("v1a_correct", "v1a_zero"), ("v1a_correct", "v1a_permuted"),
        ("v1a_correct", "v1a_sc_qw"), ("v1a_correct", "v1a_sw_qc"),
        ("v1a_correct", "v1a_sw_qw"), ("v1a_correct", "v1a_foreign"),
        ("v1b_correct", "v1b_zero"), ("v1b_correct", "v1b_permuted"),
        ("v1b_correct", "v1b_sc_qw"), ("v1b_correct", "v1b_sw_qc"),
        ("v1b_correct", "v1b_sw_qw"), ("v1b_correct", "v1b_foreign"),
        ("v1b_correct", "ligand_d0"),
    ]
    cluster_contrasts = [cluster_bootstrap_contrast(
        target_loss, target_cluster, correct, control, config.bootstrap_draws,
        SEEDS[0] + index) for index, (correct, control) in enumerate(comparisons)]
    lookup = {(row["correct"], row["control"]): row for row in cluster_contrasts}
    candidate = (
        summary["v1b_correct"]["mse"] <= summary["v1a_correct"]["mse"]
        and lookup[("v1b_correct", "v1b_permuted")]["cluster_macro_mse_reduction"]
        > lookup[("v1a_correct", "v1a_permuted")]["cluster_macro_mse_reduction"]
        and lookup[("v1b_correct", "v1b_sc_qw")]["cluster_macro_mse_reduction"]
        > lookup[("v1a_correct", "v1a_sc_qw")]["cluster_macro_mse_reduction"]
        and lookup[("v1b_correct", "v1b_zero")]["one_sided_95_lcb"] > 0
        and lookup[("v1b_correct", "ligand_d0")]["one_sided_95_lcb"] > 0
    )
    result = json_safe({
        "schema": "MetaSieve.V1DevelopmentResult.v1",
        "TERMINAL_VERDICT": ("V1_DEVELOPMENT_CANDIDATE_SELECTED" if candidate
                             else "V1_DEVELOPMENT_REPAIR_NOT_SELECTED"),
        "config": asdict(config), "seeds": list(SEEDS),
        "development_only": True, "main_v0_test_values_used": 0,
        "scientific_confirmation": False, "production_migration_authorized": False,
        "data": {"development_seal_sha256": sha256(DEV / "manifest.json"),
                 "wrong_features_manifest_sha256": sha256(DEV / "metaval_wrong_features.manifest.json"),
                 "source_k5_tasks": len(tasks), "source_cdhit40_clusters": len(clusters),
                 "metaval_k5_tasks": 37, "metaval_eligible_clusters": 9},
        "prediction_artifact": {"path": str(prediction_path.relative_to(ROOT)),
                                "sha256_before_scoring": prediction_hash,
                                "query_labels_in_artifact": False, "rows": len(all_rows)},
        "point_metrics_target_macro": summary,
        "cluster_bootstrap_contrasts": cluster_contrasts,
        "candidate_criteria_pass": candidate,
        "train_diagnostics": diagnostics, "checkpoints": checkpoints,
        "law_metrics": "NA_NOT_ADMITTED",
    })
    (OUT / "V1_DEVELOPMENT_RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    (OUT / "draw_metrics.json").write_text(
        json.dumps(json_safe(draw_metrics), indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    print(json.dumps(run(device=args.device), indent=2, sort_keys=True))
