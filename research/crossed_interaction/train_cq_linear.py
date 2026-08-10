"""Fit the preregistered linear witness on complete-panel quotient features."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from research.crossed_interaction.quotient_operator import nuisance_basis, panel_design
from scripts.structure_sources.rcsb import sha256_file


RIDGE_GRID = (1e-6, 1e-4, 1e-2, 1.0)


def read_jsonl_gz(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def panel_quotient(rows: list[dict], features: np.ndarray) -> tuple[np.ndarray, np.ndarray, int]:
    design = panel_design(
        [row["target_id"] for row in rows], [row["ligand_id"] for row in rows]
    )
    basis, retained_rank = nuisance_basis(design)
    response = np.asarray([row["pK"] for row in rows], dtype=np.float64)
    response_q = response - basis @ (basis.T @ response)
    feature_q = features - basis @ (basis.T @ features)
    return response_q, feature_q, retained_rank


def fit_ridge(panels: list[dict], ridge: float) -> np.ndarray:
    if not panels:
        raise ValueError("at least one panel is required")
    dimension = panels[0]["X"].shape[1]
    gram = np.zeros((dimension, dimension), dtype=np.float64)
    target = np.zeros(dimension, dtype=np.float64)
    for panel in panels:
        weight = 1.0 / (len(panels) * panel["rank"])
        gram += weight * panel["X"].T @ panel["X"]
        target += weight * panel["X"].T @ panel["y"]
    return np.linalg.solve(gram + ridge * np.eye(dimension), target)


def panel_losses(panels: list[dict], weights: np.ndarray, arm: str) -> dict[str, float]:
    result = {}
    for panel in panels:
        residual = panel["y"] - panel[arm] @ weights
        result[panel["panel_id"]] = float(residual @ residual / panel["rank"])
    return result


def macro_metrics(panels: list[dict], losses: dict[str, float]) -> dict:
    by_component: dict[str, list[float]] = defaultdict(list)
    for panel in panels:
        by_component[panel["component"]].append(losses[panel["panel_id"]])
    values = np.asarray([np.mean(value) for value in by_component.values()], dtype=np.float64)
    return {"loss": float(values.mean()), "rmse": float(np.sqrt(values.mean())),
            "components": len(values)}


def paired_bootstrap(
    panels: list[dict], left: dict[str, float], right: dict[str, float],
    draws: int = 10000, seed: int = 20260810,
) -> dict:
    by_component: dict[str, list[float]] = defaultdict(list)
    for panel in panels:
        panel_id = panel["panel_id"]
        by_component[panel["component"]].append(right[panel_id] - left[panel_id])
    values = np.asarray([np.mean(value) for value in by_component.values()], dtype=np.float64)
    rng = np.random.default_rng(seed)
    sampled = values[rng.integers(0, len(values), size=(draws, len(values)))].mean(axis=1)
    return {"point_loss_reduction": float(values.mean()),
            "ci95": [float(np.quantile(sampled, 0.025)), float(np.quantile(sampled, 0.975))]}


def run(corpus: Path, features_path: Path, output: Path) -> dict:
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    panels_meta = read_jsonl_gz(corpus / "panels.jsonl.gz")
    by_cell = {row["cell_id"]: row for row in cells}
    with np.load(features_path, allow_pickle=False) as stored:
        ids = [str(value) for value in stored["cell_id"]]
        if ids != [row["cell_id"] for row in cells]:
            raise ValueError("feature rows do not exactly match corpus cells")
        arms = {name: stored[name].astype(np.float64) for name in
                ("correct", "foreign_ligand", "deranged_protein")}
    feature_index = {cell_id: index for index, cell_id in enumerate(ids)}
    prepared = []
    for panel in panels_meta:
        rows = [by_cell[cell_id] for cell_id in panel["cell_ids"]]
        indices = [feature_index[cell_id] for cell_id in panel["cell_ids"]]
        yq, correct_q, rank = panel_quotient(rows, arms["correct"][indices])
        if rank != panel["retained_rank"]:
            raise ValueError("stored and recomputed quotient ranks differ")
        prepared.append({
            "panel_id": panel["panel_id"], "component": panel["dependency_component"],
            "split": panel["split"], "rank": rank, "y": yq, "X": correct_q,
            "correct": correct_q,
            "foreign_ligand": panel_quotient(rows, arms["foreign_ligand"][indices])[1],
            "deranged_protein": panel_quotient(rows, arms["deranged_protein"][indices])[1],
        })
    train = [panel for panel in prepared if panel["split"] == "train"]
    development = [panel for panel in prepared if panel["split"] == "development"]
    components = sorted({panel["component"] for panel in train})
    folds = {component: int(hashlib.sha256(component.encode()).hexdigest()[:8], 16) % 5
             for component in components}
    cv = {}
    for ridge in RIDGE_GRID:
        fold_losses = []
        for fold in range(5):
            fit = [panel for panel in train if folds[panel["component"]] != fold]
            held = [panel for panel in train if folds[panel["component"]] == fold]
            if not fit or not held:
                continue
            weights = fit_ridge(fit, ridge)
            fold_losses.append(macro_metrics(held, panel_losses(held, weights, "correct"))["loss"])
        cv[str(ridge)] = float(np.mean(fold_losses))
    selected = min(RIDGE_GRID, key=lambda value: (cv[str(value)], value))
    weights = fit_ridge(train, selected)
    correct = panel_losses(development, weights, "correct")
    foreign = panel_losses(development, weights, "foreign_ligand")
    deranged = panel_losses(development, weights, "deranged_protein")
    zero = {panel["panel_id"]: float(panel["y"] @ panel["y"] / panel["rank"])
            for panel in development}
    correct_metrics = macro_metrics(development, correct)
    zero_metrics = macro_metrics(development, zero)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez(output, weights=weights, ridge=np.asarray(selected))
    result = {
        "schema": "MetaSieve.BindingDB.CQLinearWitness.v1",
        "verdict": "CQ_LINEAR_DEVELOPMENT_WITNESS_EVALUATED",
        "selected_ridge": selected,
        "train_component_cv": cv,
        "development": {
            "correct": correct_metrics,
            "zero_interaction": zero_metrics,
            "foreign_ligand": macro_metrics(development, foreign),
            "deranged_protein": macro_metrics(development, deranged),
            "explained_fraction_vs_zero": 1.0 - correct_metrics["loss"] / zero_metrics["loss"],
            "correct_vs_zero": paired_bootstrap(development, correct, zero),
            "correct_vs_foreign_ligand": paired_bootstrap(development, correct, foreign),
            "correct_vs_deranged_protein": paired_bootstrap(development, correct, deranged),
        },
        "weights_norm": float(np.linalg.norm(weights)),
        "nonzero_coefficients": int(np.count_nonzero(np.abs(weights) > 1e-12)),
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
        "features_sha256": sha256_file(features_path),
        "weights_sha256": sha256_file(output),
        "claim_boundary": "development witness only; not affinity or few-shot admission",
    }
    output.with_suffix(".report.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args.corpus, args.features, args.output), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
