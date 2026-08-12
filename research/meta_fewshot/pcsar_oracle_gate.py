"""Source-only PC-SAR oracle Gate.

This fail-fast audit tests whether source targets contain target-specific SAR
geometry before training a protein-conditioned V2 representation. It compares
three arms on identical target-internal support/query episodes:

LEVEL   support mean only
GLOBAL  one source-trained ligand kernel shared by all targets
ORACLE  one free target-specific ligand kernel fitted from that target support

ORACLE is not a production model. It is an upper-bound diagnostic: if free
target-specific SAR does not beat GLOBAL/LEVEL, PC-SAR should stop.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from research.meta_fewshot.train_main_v0 import CORPUS, FEATURES, sha256


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "report/meta_fewshot/pcsar_oracle_gate"


def read_cells(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def stable_seed(*parts: object) -> int:
    text = "|".join(str(part) for part in parts)
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:16], 16)


def draw_episode(
        indices: np.ndarray, *, seed: int, k: int, max_query: int | None) -> tuple[np.ndarray, np.ndarray]:
    if len(indices) <= k:
        raise ValueError("episode needs more rows than support size")
    rng = np.random.default_rng(seed)
    order = rng.permutation(indices)
    support = order[:k]
    query = order[k:]
    if max_query is not None and len(query) > max_query:
        query = rng.choice(query, size=max_query, replace=False)
    return np.asarray(support, dtype=np.int64), np.asarray(query, dtype=np.int64)


def fit_ridge(x: np.ndarray, y: np.ndarray, ridge: float, *,
              center: bool = True) -> dict:
    if ridge <= 0:
        raise ValueError("ridge must be strictly positive")
    if center:
        mean = x.mean(axis=0)
        y_mean = float(y.mean())
    else:
        mean = np.zeros(x.shape[1], dtype=np.float64)
        y_mean = 0.0
    scale = x.std(axis=0)
    scale[scale < 1e-6] = 1.0
    xs = (x - mean) / scale
    yc = y - y_mean
    identity = np.eye(xs.shape[1], dtype=np.float64)
    weights = np.linalg.solve(xs.T @ xs + ridge * identity, xs.T @ yc)
    return {
        "mean": mean,
        "scale": scale,
        "y_mean": y_mean,
        "weights": weights,
        "ridge": ridge,
        "feature_dim": int(x.shape[1]),
    }


def predict(model: dict, x: np.ndarray) -> np.ndarray:
    return ((x - model["mean"]) / model["scale"]) @ model["weights"] + model["y_mean"]


def fit_pca_basis(x: np.ndarray, dim: int) -> dict:
    if dim <= 0 or dim > x.shape[1]:
        raise ValueError("PCA dim must be in [1, feature_dim]")
    mean = x.mean(axis=0)
    centered = x - mean
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    basis = vt[:dim].T
    return {"mean": mean, "basis": basis, "dim": int(dim)}


def apply_basis(model: dict, x: np.ndarray) -> np.ndarray:
    return (x - model["mean"]) @ model["basis"]


def load_source_data(corpus: Path, features: Path) -> tuple[list[dict], np.ndarray, np.ndarray, dict]:
    cells = read_cells(corpus / "cells.jsonl.gz")
    with np.load(features, allow_pickle=False) as stored:
        if stored["cell_id"].tolist() != [row["cell_id"] for row in cells]:
            raise ValueError("feature rows do not match corpus cell order")
        ligand = stored["correct"].astype(np.float64)
    train_indices = np.asarray([
        index for index, row in enumerate(cells) if row["split"] == "meta_train"
    ], dtype=np.int64)
    mean = ligand[train_indices].mean(axis=0)
    scale = ligand[train_indices].std(axis=0)
    scale[scale < 1e-6] = 1.0
    ligand = (ligand - mean) / scale
    y = np.asarray([row["pK"] for row in cells], dtype=np.float64)
    y_mean = float(y[train_indices].mean())
    y_scale = float(y[train_indices].std())
    y = (y - y_mean) / y_scale
    metadata = {
        "feature_source": "main_v0_tbasis_correct_as_ligand_kernel_proxy",
        "feature_dim": int(ligand.shape[1]),
        "y_mean": y_mean,
        "y_scale": y_scale,
        "corpus_manifest_sha256": sha256(corpus / "manifest.json"),
        "features_sha256": sha256(features),
    }
    return cells, ligand, y, metadata


def source_tasks(cells: list[dict], *, split: str, min_rows: int) -> dict[str, np.ndarray]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(cells):
        if row["split"] == split:
            grouped[row["target_id"]].append(index)
    return {
        target: np.asarray(indices, dtype=np.int64)
        for target, indices in sorted(grouped.items())
        if len(indices) >= min_rows
    }


def component_contrast(rows: list[dict], correct: str, control: str, *,
                       draws: int, seed: int) -> dict:
    by_target: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        by_target[row["target_id"]][row["arm"]].append(row["squared_error"])
    targets = sorted([
        target for target, arms in by_target.items()
        if correct in arms and control in arms
    ])
    if len(targets) < 2:
        raise ValueError("target bootstrap needs at least two targets")
    delta = np.asarray([
        float(np.mean(by_target[target][control]) - np.mean(by_target[target][correct]))
        for target in targets
    ], dtype=np.float64)
    rng = np.random.default_rng(seed)
    samples = delta[rng.integers(0, len(delta), size=(draws, len(delta)))].mean(axis=1)
    lcb = float(np.quantile(samples, 0.05))
    return {
        "correct": correct,
        "control": control,
        "targets": len(targets),
        "target_macro_reduction": float(delta.mean()),
        "one_sided_95_lcb": lcb,
        "pass": bool(lcb > 0.0),
    }


def summarize(rows: list[dict], arm: str) -> dict:
    selected = [row for row in rows if row["arm"] == arm]
    return {
        "predictions": len(selected),
        "mse": float(np.mean([row["squared_error"] for row in selected])),
        "target_macro_mse": float(np.mean([
            np.mean([row["squared_error"] for row in selected if row["target_id"] == target])
            for target in sorted({row["target_id"] for row in selected})
        ])),
    }


def run(
        corpus: Path = CORPUS, features: Path = FEATURES, output: Path = OUT,
        k: int = 5, min_query: int = 3, draws: int = 5,
        max_query: int = 128, ridge: float = 1.0, oracle_dim: int = 5,
        bootstrap_draws: int = 9999, seed: int = 20260812) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    cells, x, y, metadata = load_source_data(corpus, features)
    tasks = source_tasks(cells, split="meta_train", min_rows=k + min_query)
    if len(tasks) < 2:
        raise ValueError("oracle Gate needs at least two source tasks")
    train_indices = np.concatenate(list(tasks.values()))
    basis_model = fit_pca_basis(x[train_indices], oracle_dim)
    z = apply_basis(basis_model, x)
    global_model = fit_ridge(z[train_indices], y[train_indices], ridge)
    rows = []
    for target, indices in tasks.items():
        for draw in range(draws):
            support, query = draw_episode(
                indices, seed=stable_seed(seed, target, draw), k=k,
                max_query=max_query)
            oracle_model = fit_ridge(z[support], y[support], ridge)
            level_prediction = np.full(len(query), float(y[support].mean()))
            predictions = {
                "LEVEL": level_prediction,
                "GLOBAL": predict(global_model, z[query]),
                "ORACLE": predict(oracle_model, z[query]),
            }
            for arm, values in predictions.items():
                for cell_index, estimate in zip(query, values):
                    rows.append({
                        "target_id": target,
                        "draw": draw,
                        "arm": arm,
                        "cell_id": cells[int(cell_index)]["cell_id"],
                        "squared_error": float(((y[int(cell_index)] - estimate) * metadata["y_scale"]) ** 2),
                    })
    contrasts = [
        component_contrast(rows, "ORACLE", "GLOBAL", draws=bootstrap_draws, seed=seed),
        component_contrast(rows, "ORACLE", "LEVEL", draws=bootstrap_draws, seed=seed + 1),
        component_contrast(rows, "GLOBAL", "LEVEL", draws=bootstrap_draws, seed=seed + 2),
    ]
    gates = {
        "source_targets_ge_10": len(tasks) >= 10,
        "oracle_beats_global": contrasts[0]["pass"],
        "oracle_beats_level": contrasts[1]["pass"],
    }
    verdict = (
        "PCSAR_ORACLE_GATE_PASS_HEADROOM"
        if all(gates.values())
        else "PCSAR_ORACLE_GATE_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.PCSAROracleGate.v1",
        "hypothesis": (
            "Source targets contain target-specific SAR kernel headroom: a free "
            "target-ID ligand kernel fitted only from support beats protein-blind "
            "and level-only controls on held-out query ligands."),
        "corpus": {
            "tasks": len(tasks),
            "cells": len(cells),
            "split": "meta_train",
            **metadata,
        },
        "config": {
            "k": k,
            "min_query": min_query,
            "draws": draws,
            "max_query": max_query,
            "ridge": ridge,
            "oracle_dim": oracle_dim,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
        },
        "train_summary": {
            "global_train_mse_standardized": float(np.square(
                y[train_indices] - predict(global_model, z[train_indices])).mean()),
            "global_feature_dim": global_model["feature_dim"],
        },
        "development_summary": {
            arm: summarize(rows, arm) for arm in ("LEVEL", "GLOBAL", "ORACLE")
        },
        "development_contrasts": contrasts,
        "gates": gates,
        "development_training_authorized": verdict.endswith("PASS_HEADROOM"),
        "v1_integration_authorized": False,
        "biological_claim_authorized": False,
        "TERMINAL_VERDICT": verdict,
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    (output / "development_rows.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
        encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--features", type=Path, default=FEATURES)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--min-query", type=int, default=3)
    parser.add_argument("--draws", type=int, default=5)
    parser.add_argument("--max-query", type=int, default=128)
    parser.add_argument("--ridge", type=float, default=1.0)
    parser.add_argument("--oracle-dim", type=int, default=5)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, features=args.features, output=args.output,
        k=args.k, min_query=args.min_query, draws=args.draws,
        max_query=args.max_query, ridge=args.ridge, oracle_dim=args.oracle_dim,
        bootstrap_draws=args.bootstrap_draws, seed=args.seed)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
