"""Audit whether BindingDB SAR-delta pairs remain shortcut-prone after reversal.

Gate 0 for the UniPert-inspired route: every forward matched pair is augmented
with its reverse pair. A target-main predictor should then lose to zero because
the empirical distribution enforces delta antisymmetry for each target.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from research.crossed_interaction.bindingdb_cq_r0 import sha256_file
from research.crossed_interaction.train_bindingdb_sardelta_attribution import (
    CORPUS,
    OUT as ATTRIBUTION_OUT,
    build_pairs,
    component_contrast,
    fit_positive_ridge_no_intercept,
    predict,
    stack,
    summarize,
)
from research.source_affinity.train_chembl_assay_sardelta import PAIR_SIMILARITY


OUT = ATTRIBUTION_OUT.parent / "bindingdb_sardelta_symmetry_gate0"


def augment_reverse_pairs(pairs: list[dict]) -> list[dict]:
    augmented = []
    for pair in pairs:
        augmented.append(pair)
        reversed_pair = dict(pair)
        reversed_pair["left_cell_id"] = pair["right_cell_id"]
        reversed_pair["right_cell_id"] = pair["left_cell_id"]
        reversed_pair["delta_pK"] = -float(pair["delta_pK"])
        reversed_pair["ligand_delta"] = -pair["ligand_delta"]
        reversed_pair["concat"] = np.concatenate([pair["protein"], reversed_pair["ligand_delta"]])
        reversed_pair["interaction"] = -pair["interaction"]
        reversed_pair["reversed_from"] = f"{pair['left_cell_id']}->{pair['right_cell_id']}"
        augmented.append(reversed_pair)
    return augmented


def antisymmetry_error(rows: list[dict], arm: str) -> dict:
    by_pair = {}
    errors = []
    for row in rows:
        key = tuple(sorted((row["left_cell_id"], row["right_cell_id"])))
        if key in by_pair:
            errors.append(abs(by_pair.pop(key) + row[f"{arm}_prediction"]))
        else:
            by_pair[key] = row[f"{arm}_prediction"]
    if not errors:
        return {"pairs": 0, "mean_abs_sum": 0.0, "max_abs_sum": 0.0}
    return {
        "pairs": len(errors),
        "mean_abs_sum": float(np.mean(errors)),
        "max_abs_sum": float(np.max(errors)),
    }


def run(
        corpus: Path = CORPUS, output: Path = OUT, ridge: float = 100.0,
        max_pairs_per_group: int | None = 100,
        bootstrap_draws: int = 9999, seed: int = 20260812) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    train_forward, train_meta = build_pairs(
        corpus, split="train", max_pairs_per_group=max_pairs_per_group)
    dev_forward, dev_meta = build_pairs(
        corpus, split="development", max_pairs_per_group=max_pairs_per_group)
    train_pairs = augment_reverse_pairs(train_forward)
    dev_pairs = augment_reverse_pairs(dev_forward)
    y_train = np.asarray([pair["delta_pK"] for pair in train_pairs], dtype=np.float64)
    y_dev = np.asarray([pair["delta_pK"] for pair in dev_pairs], dtype=np.float64)
    models = {
        "L": fit_positive_ridge_no_intercept(stack(train_pairs, "ligand_delta"), y_train, ridge),
        "P": fit_positive_ridge_no_intercept(stack(train_pairs, "protein"), y_train, ridge),
        "A": fit_positive_ridge_no_intercept(stack(train_pairs, "concat"), y_train, ridge),
    }
    predictions = {
        "Z": np.zeros_like(y_dev),
        "L": predict(models["L"], stack(dev_pairs, "ligand_delta")),
        "P": predict(models["P"], stack(dev_pairs, "protein")),
        "A": predict(models["A"], stack(dev_pairs, "concat")),
    }
    rows = []
    for index, (pair, true) in enumerate(zip(dev_pairs, y_dev)):
        row = {
            "dependency_component": pair["dependency_component"],
            "target_id": pair["target_id"],
            "scaffold": pair["scaffold"],
            "left_cell_id": pair["left_cell_id"],
            "right_cell_id": pair["right_cell_id"],
            "delta_pK": float(true),
        }
        for arm, values in predictions.items():
            estimate = float(values[index])
            row[f"{arm}_prediction"] = estimate
            row[f"{arm}_squared_error"] = float((true - estimate) ** 2)
        rows.append(row)
    contrasts = [
        component_contrast(rows, "P", "Z", draws=bootstrap_draws, seed=seed),
        component_contrast(rows, "A", "Z", draws=bootstrap_draws, seed=seed + 1),
        component_contrast(rows, "A", "L", draws=bootstrap_draws, seed=seed + 2),
        component_contrast(rows, "L", "Z", draws=bootstrap_draws, seed=seed + 3),
    ]
    antisymmetry = {arm: antisymmetry_error(rows, arm) for arm in ("L", "P", "A")}
    gates = {
        "development_components_ge_5": len({row["dependency_component"] for row in rows}) >= 5,
        "target_main_not_better_than_zero": not contrasts[0]["pass"],
        "ligand_delta_antisymmetry": antisymmetry["L"]["max_abs_sum"] <= 1e-9,
        "additive_concat_antisymmetry": antisymmetry["A"]["max_abs_sum"] <= 1e-9,
    }
    verdict = (
        "BINDINGDB_SARDELTA_SYMMETRY_GATE0_PASS"
        if all(gates.values())
        else "BINDINGDB_SARDELTA_SYMMETRY_GATE0_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.BindingDBSARDeltaSymmetryGate0.v1",
        "hypothesis": (
            "Forward+reverse SAR-delta augmentation removes target-main pair-order "
            "shortcuts before testing any richer UniPert-inspired bridge."),
        "corpus": {
            "manifest_sha256": sha256_file(corpus / "manifest.json"),
            "train_forward_pairs": train_meta,
            "development_forward_pairs": dev_meta,
            "train_augmented_pairs": len(train_pairs),
            "development_augmented_pairs": len(dev_pairs),
        },
        "config": {
            "ridge": ridge,
            "max_pairs_per_group": max_pairs_per_group,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
            "pair_similarity_threshold": PAIR_SIMILARITY,
        },
        "train_summary": {
            arm: {
                "feature_dim": models[arm]["feature_dim"],
                "train_mse": models[arm]["train_mse"],
            }
            for arm in ("L", "P", "A")
        },
        "development_summary": {
            arm: summarize(rows, arm)
            for arm in ("Z", "L", "P", "A")
        },
        "development_contrasts": contrasts,
        "antisymmetry_audit": antisymmetry,
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
    (output / "development_rows.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
        encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ridge", type=float, default=100.0)
    parser.add_argument("--max-pairs-per-group", type=int, default=100)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, output=args.output, ridge=args.ridge,
        max_pairs_per_group=args.max_pairs_per_group,
        bootstrap_draws=args.bootstrap_draws, seed=args.seed)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
