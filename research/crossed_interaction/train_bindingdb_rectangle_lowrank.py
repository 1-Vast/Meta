"""Source-only low-rank oracle for BindingDB rectangle interactions.

X2 follows the label-only X1 Gate. It tests whether 2x2 rectangle
double-differences have a transferable low-rank target x transformation
structure. No protein encoder is trained here.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from research.crossed_interaction.audit_bindingdb_rectangle_interaction import (
    CORPUS,
    OUT as RECTANGLE_OUT,
    build_rectangles,
)


OUT = RECTANGLE_OUT.parent / "bindingdb_rectangle_lowrank_x2"


def transformation_key(row: dict) -> str:
    left, right = sorted([row["ligand_a"], row["ligand_b"]])
    return f"{left}>{right}"


def materialize_matrix(rows: list[dict], split: str, *,
                       transforms: list[str] | None = None,
                       min_targets: int = 2) -> tuple[np.ndarray, dict]:
    selected = [row for row in rows if row["split"] == split]
    by_target_transform: dict[tuple[str, str], list[float]] = defaultdict(list)
    component_by_target: dict[str, str] = {}
    for row in selected:
        target_pair = "|".join(sorted([row["target_a"], row["target_b"]]))
        transform = transformation_key(row)
        by_target_transform[(target_pair, transform)].append(float(row["rectangle"]))
        component_by_target[target_pair] = row["dependency_component"]
    target_pairs = sorted({key[0] for key in by_target_transform})
    transform_filter = None if transforms is None else set(transforms)
    if transforms is None:
        transforms = sorted({key[1] for key in by_target_transform})
    t_index = {value: index for index, value in enumerate(target_pairs)}
    m_index = {value: index for index, value in enumerate(transforms)}
    matrix = np.full((len(target_pairs), len(transforms)), np.nan, dtype=np.float64)
    for (target_pair, transform), values in by_target_transform.items():
        if transform_filter is not None and transform not in transform_filter:
            continue
        matrix[t_index[target_pair], m_index[transform]] = float(np.mean(values))
    if min_targets > 1:
        valid_transform = np.sum(np.isfinite(matrix), axis=0) >= min_targets
        matrix = matrix[:, valid_transform]
        transforms = [value for value, keep in zip(transforms, valid_transform) if keep]
    metadata = {
        "split": split,
        "target_pairs": len(target_pairs),
        "transforms": len(transforms),
        "observed_cells": int(np.isfinite(matrix).sum()),
        "density": float(np.isfinite(matrix).mean()) if matrix.size else 0.0,
        "components": len(set(component_by_target.values())),
        "transform_keys": transforms,
        "target_pair_components": component_by_target,
    }
    return matrix, metadata


def double_center_observed(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    observed = np.isfinite(matrix)
    filled = np.where(observed, matrix, 0.0)
    row_count = observed.sum(axis=1, keepdims=True)
    col_count = observed.sum(axis=0, keepdims=True)
    row_mean = np.divide(
        filled.sum(axis=1, keepdims=True), row_count,
        out=np.zeros_like(row_count, dtype=np.float64), where=row_count > 0)
    col_mean = np.divide(
        filled.sum(axis=0, keepdims=True), col_count,
        out=np.zeros_like(col_count, dtype=np.float64), where=col_count > 0)
    grand = float(filled.sum() / max(observed.sum(), 1))
    residual = matrix - row_mean - col_mean + grand
    return residual, grand


def lowrank_reconstruct(train: np.ndarray, dev: np.ndarray, rank: int) -> tuple[np.ndarray, dict]:
    train_resid, grand = double_center_observed(train)
    observed_train = np.isfinite(train_resid)
    filled = np.where(observed_train, train_resid, 0.0)
    if not observed_train.any():
        raise ValueError("low-rank oracle needs observed train cells")
    u, s, vt = np.linalg.svd(filled, full_matrices=False)
    basis = vt[:rank].T
    dev_resid, _ = double_center_observed(dev)
    prediction = np.zeros_like(dev_resid)
    observed_dev = np.isfinite(dev_resid)
    for row in range(dev_resid.shape[0]):
        mask = observed_dev[row]
        if not mask.any():
            continue
        b = basis[mask]
        y = dev_resid[row, mask]
        coef = np.linalg.lstsq(b, y, rcond=None)[0]
        prediction[row, mask] = b @ coef
    metadata = {
        "rank": rank,
        "train_singular_values": [float(value) for value in s[:min(10, len(s))]],
        "train_grand_mean": grand,
    }
    return prediction, metadata


def score(dev: np.ndarray, prediction: np.ndarray) -> dict:
    dev_resid, _ = double_center_observed(dev)
    mask = np.isfinite(dev_resid)
    if not mask.any():
        raise ValueError("score needs observed development cells")
    zero_error = np.square(dev_resid[mask])
    pred_error = np.square(dev_resid[mask] - prediction[mask])
    return {
        "observed_cells": int(mask.sum()),
        "zero_mse": float(zero_error.mean()),
        "prediction_mse": float(pred_error.mean()),
        "reduction": float(zero_error.mean() - pred_error.mean()),
    }


def run(
        corpus: Path = CORPUS, output: Path = OUT,
        ranks: tuple[int, ...] = (1, 2, 3), seed: int = 20260812) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    rows, rectangle_metadata = build_rectangles(corpus)
    train_all, train_all_meta = materialize_matrix(rows, "train")
    dev_all, dev_all_meta = materialize_matrix(rows, "development")
    shared = sorted(
        set(train_all_meta["transform_keys"]) & set(dev_all_meta["transform_keys"]))
    train, train_meta = materialize_matrix(rows, "train", transforms=shared, min_targets=1)
    dev, dev_meta = materialize_matrix(rows, "development", transforms=shared, min_targets=1)
    train_meta["pre_shared_transforms"] = train_all_meta["transforms"]
    dev_meta["pre_shared_transforms"] = dev_all_meta["transforms"]
    train_meta["shared_transforms"] = len(shared)
    dev_meta["shared_transforms"] = len(shared)
    if not shared or train_meta["observed_cells"] == 0 or dev_meta["observed_cells"] == 0:
        result = {
            "schema": "MetaSieve.BindingDBRectangleLowRankX2.v1",
            "hypothesis": (
                "Label-side rectangle interactions contain a transferable low-rank "
                "target-pair x ligand-transformation structure."),
            "rectangle_corpus": rectangle_metadata,
            "matrix": {
                "train": train_meta,
                "development": dev_meta,
            },
            "config": {
                "ranks": list(ranks),
                "seed": seed,
                "labels_used_for_training": True,
                "protein_encoder_used": False,
            },
            "rank_results": [],
            "best_rank": None,
            "gates": {
                "x2_train_observed_ge_1000": train_meta["observed_cells"] >= 1000,
                "x2_development_observed_ge_100": dev_meta["observed_cells"] >= 100,
                "x2_shared_transforms_ge_10": len(shared) >= 10,
                "x2_best_lowrank_beats_zero": False,
            },
            "failure_reason": "no_shared_train_development_ligand_transformations",
            "development_training_authorized": False,
            "v1_integration_authorized": False,
            "biological_claim_authorized": False,
            "TERMINAL_VERDICT": "BINDINGDB_RECTANGLE_LOWRANK_X2_FAIL_CLOSED",
        }
        output.mkdir(parents=True, exist_ok=False)
        (output / "RESULT.json").write_text(
            json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8")
        return result
    rank_results = []
    for rank in ranks:
        prediction, metadata = lowrank_reconstruct(train, dev, rank)
        rank_results.append({
            **metadata,
            "development": score(dev, prediction),
        })
    best = min(rank_results, key=lambda row: row["development"]["prediction_mse"])
    gates = {
        "x2_train_observed_ge_1000": train_meta["observed_cells"] >= 1000,
        "x2_development_observed_ge_100": dev_meta["observed_cells"] >= 100,
        "x2_shared_transforms_ge_10": len(shared) >= 10,
        "x2_best_lowrank_beats_zero": best["development"]["reduction"] > 0.0,
    }
    verdict = (
        "BINDINGDB_RECTANGLE_LOWRANK_X2_PASS"
        if all(gates.values())
        else "BINDINGDB_RECTANGLE_LOWRANK_X2_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.BindingDBRectangleLowRankX2.v1",
        "hypothesis": (
            "Label-side rectangle interactions contain a transferable low-rank "
            "target-pair x ligand-transformation structure."),
        "rectangle_corpus": rectangle_metadata,
        "matrix": {
            "train": train_meta,
            "development": dev_meta,
        },
        "config": {
            "ranks": list(ranks),
            "seed": seed,
            "labels_used_for_training": True,
            "protein_encoder_used": False,
        },
        "rank_results": rank_results,
        "best_rank": best["rank"],
        "gates": gates,
        "development_training_authorized": verdict.endswith("PASS"),
        "v1_integration_authorized": False,
        "biological_claim_authorized": False,
        "TERMINAL_VERDICT": verdict,
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ranks", type=int, nargs="+", default=[1, 2, 3])
    parser.add_argument("--seed", type=int, default=20260812)
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, output=args.output,
        ranks=tuple(args.ranks), seed=args.seed)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
