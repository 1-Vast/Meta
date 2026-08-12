"""Train a source-only quotient observable for crossed Ki panels.

This is an admission Gate for a biological pair coordinate, not a V1 predictor.
Every panel is projected away from target-only and ligand-only main effects
before any feature model is trained or scored.
"""
from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from research.crossed_interaction.bindingdb_cq_r0 import sha256_file


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
FEATURES = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_tbasis_features.npz"
OUT = ROOT / "report/crossed_interaction/cq_observable_gate1"
ARMS = ("correct", "deranged_protein", "foreign_ligand")


@dataclass(frozen=True)
class QuotientBlock:
    panel_id: str
    split: str
    dependency_component: str
    retained_rank: int
    y: np.ndarray
    features: dict[str, np.ndarray]


def read_jsonl_gz(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def additive_residual(
        target_ids: list[str], ligand_ids: list[str],
        values: np.ndarray) -> tuple[np.ndarray, int, float]:
    if values.ndim == 1:
        values = values[:, None]
    targets = {value: index for index, value in enumerate(sorted(set(target_ids)))}
    ligands = {value: index for index, value in enumerate(sorted(set(ligand_ids)))}
    design = np.zeros((len(target_ids), 1 + len(targets) + len(ligands)), dtype=np.float64)
    for row, (target, ligand) in enumerate(zip(target_ids, ligand_ids)):
        design[row, 0] = 1.0
        design[row, 1 + targets[target]] = 1.0
        design[row, 1 + len(targets) + ligands[ligand]] = 1.0
    fitted = design @ np.linalg.lstsq(design, values, rcond=None)[0]
    residual = values - fitted
    rank = int(np.linalg.matrix_rank(design))
    retained_rank = len(target_ids) - rank
    orthogonality = float(np.max(np.abs(design.T @ residual))) if residual.size else 0.0
    return residual.squeeze(), retained_rank, orthogonality


def load_blocks(corpus: Path, features_path: Path) -> tuple[list[QuotientBlock], dict]:
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
    with np.load(features_path, allow_pickle=False) as stored:
        cell_ids = stored["cell_id"].tolist()
        feature_by_arm = {arm: stored[arm].astype(np.float64) for arm in ARMS}
    cell_index = {cell_id: index for index, cell_id in enumerate(cell_ids)}
    cell_by_id = {cell["cell_id"]: cell for cell in cells}
    if set(cell_by_id) != set(cell_index):
        raise ValueError("corpus cells and feature rows disagree")

    blocks: list[QuotientBlock] = []
    max_orthogonality = 0.0
    for panel in panels:
        ordered_cells = [cell_by_id[cell_id] for cell_id in panel["cell_ids"]]
        target_ids = [cell["target_id"] for cell in ordered_cells]
        ligand_ids = [cell["ligand_id"] for cell in ordered_cells]
        y_raw = np.asarray([cell["pK"] for cell in ordered_cells], dtype=np.float64)
        y_residual, retained_rank, y_orthogonality = additive_residual(
            target_ids, ligand_ids, y_raw)
        if retained_rank <= 0:
            continue
        arm_features = {}
        for arm in ARMS:
            raw = feature_by_arm[arm][[cell_index[cell["cell_id"]] for cell in ordered_cells]]
            residual, feature_rank, feature_orthogonality = additive_residual(
                target_ids, ligand_ids, raw)
            if feature_rank != retained_rank:
                raise ValueError(f"feature rank mismatch in panel {panel['panel_id']}")
            max_orthogonality = max(max_orthogonality, feature_orthogonality)
            arm_features[arm] = residual
        max_orthogonality = max(max_orthogonality, y_orthogonality)
        blocks.append(QuotientBlock(
            panel_id=panel["panel_id"],
            split=panel["split"],
            dependency_component=panel["dependency_component"],
            retained_rank=retained_rank,
            y=np.asarray(y_residual, dtype=np.float64),
            features=arm_features,
        ))
    metadata = {
        "cells": len(cells),
        "panels": len(panels),
        "blocks": len(blocks),
        "max_projection_orthogonality": max_orthogonality,
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
        "features_sha256": sha256_file(features_path),
    }
    return blocks, metadata


def _stack(blocks: list[QuotientBlock], arm: str) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.concatenate([block.features[arm] for block in blocks], axis=0),
        np.concatenate([block.y for block in blocks], axis=0),
    )


def fit_ridge(blocks: list[QuotientBlock], arm: str, ridge: float) -> dict:
    if ridge <= 0:
        raise ValueError("observable ridge must be strictly positive")
    x, y = _stack(blocks, arm)
    scale = x.std(axis=0)
    scale[scale < 1e-6] = 1.0
    x_scaled = x / scale
    identity = np.eye(x_scaled.shape[1], dtype=np.float64)
    weights = np.linalg.solve(x_scaled.T @ x_scaled + ridge * identity, x_scaled.T @ y)
    prediction = x_scaled @ weights
    return {
        "arm": arm,
        "ridge": ridge,
        "scale": scale,
        "weights": weights,
        "train_rank_weighted_mse": float(np.square(y - prediction).mean()),
    }


def predict(model: dict, block: QuotientBlock) -> np.ndarray:
    x = block.features[model["arm"]] / model["scale"]
    return x @ model["weights"]


def score_blocks(
        blocks: list[QuotientBlock], models: dict[str, dict]) -> tuple[list[dict], dict]:
    rows = []
    for block in blocks:
        predictions = {"zero_additive": np.zeros_like(block.y)}
        predictions.update({arm: predict(model, block) for arm, model in models.items()})
        for arm, prediction in predictions.items():
            squared = np.square(block.y - prediction)
            rows.append({
                "panel_id": block.panel_id,
                "split": block.split,
                "dependency_component": block.dependency_component,
                "arm": arm,
                "retained_rank": block.retained_rank,
                "rank_normalized_mse": float(squared.sum() / block.retained_rank),
                "row_mse": float(squared.mean()),
            })
    summary = {}
    for arm in sorted({row["arm"] for row in rows}):
        selected = [row for row in rows if row["arm"] == arm]
        rank_total = sum(row["retained_rank"] for row in selected)
        summary[arm] = {
            "panels": len(selected),
            "rank_weighted_mse": (
                sum(row["rank_normalized_mse"] * row["retained_rank"] for row in selected)
                / rank_total
            ),
            "panel_macro_rank_normalized_mse": float(
                np.mean([row["rank_normalized_mse"] for row in selected])),
        }
    return rows, summary


def component_metric(rows: list[dict], arm: str) -> dict[str, float]:
    grouped: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        if row["arm"] == arm:
            grouped[row["dependency_component"]].append(row["rank_normalized_mse"])
    return {component: float(np.mean(values)) for component, values in sorted(grouped.items())}


def bootstrap_contrast(
        rows: list[dict], correct: str, control: str, *,
        draws: int, seed: int) -> dict:
    correct_by_component = component_metric(rows, correct)
    control_by_component = component_metric(rows, control)
    components = sorted(set(correct_by_component) & set(control_by_component))
    if len(components) < 2:
        raise ValueError("component bootstrap needs at least two shared components")
    delta = np.asarray([
        control_by_component[component] - correct_by_component[component]
        for component in components
    ], dtype=np.float64)
    rng = np.random.default_rng(seed)
    samples = delta[rng.integers(0, len(delta), size=(draws, len(delta)))].mean(axis=1)
    lcb = float(np.quantile(samples, 0.05))
    return {
        "correct": correct,
        "control": control,
        "components": len(components),
        "component_macro_reduction": float(delta.mean()),
        "one_sided_95_lcb": lcb,
        "pass": bool(lcb > 0.0),
    }


def run(
        corpus: Path = CORPUS, features: Path = FEATURES, output: Path = OUT,
        ridge: float = 10.0, bootstrap_draws: int = 9999,
        seed: int = 20260812) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    blocks, metadata = load_blocks(corpus, features)
    train_blocks = [block for block in blocks if block.split == "train"]
    development_blocks = [block for block in blocks if block.split == "development"]
    if not train_blocks or not development_blocks:
        raise ValueError("train and development quotient blocks are required")
    models = {arm: fit_ridge(train_blocks, arm, ridge) for arm in ARMS}
    train_rows, train_summary = score_blocks(train_blocks, models)
    development_rows, development_summary = score_blocks(development_blocks, models)
    contrasts = [
        bootstrap_contrast(
            development_rows, "correct", control,
            draws=bootstrap_draws, seed=seed + index)
        for index, control in enumerate(("zero_additive", "deranged_protein", "foreign_ligand"))
    ]
    gates = {
        "projection_orthogonality": metadata["max_projection_orthogonality"] <= 1e-7,
        "development_components_ge_5": len({
            block.dependency_component for block in development_blocks}) >= 5,
        "correct_beats_zero_additive": contrasts[0]["pass"],
        "correct_beats_deranged_protein": contrasts[1]["pass"],
        "correct_beats_foreign_ligand": contrasts[2]["pass"],
    }
    verdict = (
        "CQ_OBSERVABLE_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "CQ_OBSERVABLE_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.CQObservableGate1.v1",
        "hypothesis": (
            "A source-trained quotient-space linear observable on the current "
            "pair features carries affinity-directed non-additive partner signal."),
        "literature_mechanism": {
            "hodge_cycle_space": (
                "cycle residuals remove node/main-effect potentials before testing "
                "interaction signal"),
            "bindingdb_panel_provenance": (
                "document/protocol panels are treated as source units, not IID rows"),
            "adambind_boundary": (
                "target-as-task meta-learning remains downstream; this Gate does "
                "not import MAML, scheduler, or support-noise mechanisms"),
        },
        "corpus": metadata,
        "config": {
            "ridge": ridge,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
            "arms": list(ARMS),
            "train_split": "train",
            "evaluation_split": "development",
        },
        "train_summary": train_summary,
        "development_summary": development_summary,
        "development_contrasts": contrasts,
        "gates": gates,
        "development_training_authorized": verdict.endswith("PASS_DEVELOPMENT"),
        "v1_integration_authorized": False,
        "biological_claim_authorized": False,
        "TERMINAL_VERDICT": verdict,
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    (output / "development_panel_metrics.json").write_text(
        json.dumps(development_rows, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--features", type=Path, default=FEATURES)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ridge", type=float, default=10.0)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, features=args.features, output=args.output,
        ridge=args.ridge, bootstrap_draws=args.bootstrap_draws, seed=args.seed)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
