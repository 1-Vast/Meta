"""Attribute BindingDB SAR-delta transfer to ligand and target-interaction terms.

This is a diagnostic Gate for the F-153 interpretation. It reuses the same
same-target, same-scaffold BindingDB pair protocol, but changes only the model
family:

Z   zero delta
L   ligand-delta ridge
P   target-main ridge
A   additive concat ridge [target; ligand_delta]
I   bilinear interaction ridge target x ligand_delta
IW  I with a wrong protein descriptor at development time
IS  I with a shuffled protein descriptor at development time

The L and I arms use no intercept and no feature centering so pair
antisymmetry can be audited directly.
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from research.crossed_interaction.bindingdb_cq_r0 import sha256_file
from research.crossed_interaction.train_bindingdb_sardelta_cq_bridge import CORPUS
from research.crossed_interaction.train_cq_observable import OUT as CQ_OUT
from research.crossed_interaction.train_seqchem_cq_observable import (
    donor_maps,
    ligand_descriptor,
    protein_descriptor,
    read_jsonl,
    read_jsonl_gz,
)
from research.source_affinity.train_chembl_assay_sardelta import PAIR_SIMILARITY


OUT = CQ_OUT.parent / "bindingdb_sardelta_attribution_gate1"


def interaction_feature(protein: np.ndarray, delta: np.ndarray) -> np.ndarray:
    return np.multiply.outer(protein, delta).ravel()


def fit_positive_ridge_no_intercept(x: np.ndarray, y: np.ndarray, ridge: float) -> dict:
    if ridge <= 0:
        raise ValueError("ridge must be strictly positive")
    scale = x.std(axis=0)
    scale[scale < 1e-6] = 1.0
    x_scaled = x / scale
    identity = np.eye(x_scaled.shape[1], dtype=np.float64)
    weights = np.linalg.solve(x_scaled.T @ x_scaled + ridge * identity, x_scaled.T @ y)
    prediction = x_scaled @ weights
    return {
        "ridge": ridge,
        "scale": scale,
        "weights": weights,
        "train_mse": float(np.square(y - prediction).mean()),
        "feature_dim": int(x.shape[1]),
    }


def predict(model: dict, x: np.ndarray) -> np.ndarray:
    return (x / model["scale"]) @ model["weights"]


def shuffled_targets(targets: list[str], *, seed: int) -> dict[str, str]:
    ordered = sorted(targets)
    rng = np.random.default_rng(seed)
    shuffled = ordered.copy()
    for _ in range(100):
        rng.shuffle(shuffled)
        if all(left != right for left, right in zip(ordered, shuffled)):
            return dict(zip(ordered, shuffled))
    rotated = ordered[1:] + ordered[:1]
    return dict(zip(ordered, rotated))


def build_pairs(
        corpus: Path, *, split: str,
        max_pairs_per_group: int | None) -> tuple[list[dict], dict]:
    cells = [row for row in read_jsonl_gz(corpus / "cells.jsonl.gz") if row["split"] == split]
    cell_component = {}
    for panel in read_jsonl_gz(corpus / "panels.jsonl.gz"):
        if panel["split"] != split:
            continue
        for cell_id in panel["cell_ids"]:
            cell_component[cell_id] = panel["dependency_component"]
    proteins = {
        row["sequence_sha256"]: protein_descriptor(row["sequence"])
        for row in read_jsonl(corpus / "proteins.jsonl")
    }
    ligands = {
        row["drug_key"]: {
            "descriptor": ligand_descriptor(row["smiles"]),
            "smiles": row["smiles"],
            "scaffold": row["scaffold"],
        }
        for row in read_jsonl(corpus / "ligands.jsonl")
    }
    by_target_scaffold: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for cell in cells:
        by_target_scaffold[(cell["target_id"], cell["scaffold"])].append(cell)

    pairs = []
    skipped_similarity = 0
    for (target, scaffold), group in sorted(by_target_scaffold.items()):
        if len(group) < 2:
            continue
        ordered = sorted(group, key=lambda row: (row["pK"], row["ligand_id"]))
        group_pairs = []
        for i in range(len(ordered)):
            for j in range(i + 1, len(ordered)):
                left = ordered[i]
                right = ordered[j]
                left_ligand = ligands[left["ligand_id"]]
                right_ligand = ligands[right["ligand_id"]]
                left_fp = AllChem.GetMorganFingerprintAsBitVect(
                    Chem.MolFromSmiles(left_ligand["smiles"]), 2, 1024)
                right_fp = AllChem.GetMorganFingerprintAsBitVect(
                    Chem.MolFromSmiles(right_ligand["smiles"]), 2, 1024)
                similarity = DataStructs.TanimotoSimilarity(left_fp, right_fp)
                if similarity < PAIR_SIMILARITY:
                    skipped_similarity += 1
                    continue
                delta = left_ligand["descriptor"] - right_ligand["descriptor"]
                protein = proteins[target]
                group_pairs.append({
                    "target_id": target,
                    "scaffold": scaffold,
                    "dependency_component": cell_component[left["cell_id"]],
                    "left_cell_id": left["cell_id"],
                    "right_cell_id": right["cell_id"],
                    "delta_pK": float(left["pK"] - right["pK"]),
                    "protein": protein,
                    "ligand_delta": delta,
                    "concat": np.concatenate([protein, delta]),
                    "interaction": interaction_feature(protein, delta),
                })
        if max_pairs_per_group is not None and len(group_pairs) > max_pairs_per_group:
            group_pairs = group_pairs[:max_pairs_per_group]
        pairs.extend(group_pairs)
    metadata = {
        "split": split,
        "cells": len(cells),
        "target_scaffold_groups": len(by_target_scaffold),
        "pairs": len(pairs),
        "skipped_similarity": skipped_similarity,
        "pair_similarity_threshold": PAIR_SIMILARITY,
    }
    return pairs, metadata


def protein_controls(corpus: Path, *, seed: int) -> tuple[dict[str, str], dict[str, str], dict[str, np.ndarray]]:
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    wrong, _ = donor_maps(cells)
    shuffle = shuffled_targets(sorted({cell["target_id"] for cell in cells}), seed=seed)
    proteins = {
        row["sequence_sha256"]: protein_descriptor(row["sequence"])
        for row in read_jsonl(corpus / "proteins.jsonl")
    }
    return wrong, shuffle, proteins


def stack(pairs: list[dict], field: str) -> np.ndarray:
    return np.stack([pair[field] for pair in pairs]).astype(np.float64)


def controlled_interactions(
        pairs: list[dict], target_map: dict[str, str],
        proteins: dict[str, np.ndarray]) -> np.ndarray:
    return np.stack([
        interaction_feature(proteins[target_map[pair["target_id"]]], pair["ligand_delta"])
        for pair in pairs
    ]).astype(np.float64)


def summarize(rows: list[dict], arm: str) -> dict:
    return {
        "pairs": len(rows),
        "mse": float(np.mean([row[f"{arm}_squared_error"] for row in rows])),
    }


def component_mse(rows: list[dict], arm: str) -> dict[str, float]:
    grouped: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        grouped[row["dependency_component"]].append(row[f"{arm}_squared_error"])
    return {key: float(np.mean(values)) for key, values in sorted(grouped.items())}


def component_contrast(rows: list[dict], correct: str, control: str, *, draws: int, seed: int) -> dict:
    correct_by_component = component_mse(rows, correct)
    control_by_component = component_mse(rows, control)
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


def antisymmetry_audit(model: dict, x: np.ndarray, reversed_x: np.ndarray) -> dict:
    error = np.abs(predict(model, x) + predict(model, reversed_x))
    return {
        "mean_abs_sum": float(error.mean()),
        "max_abs_sum": float(error.max()),
    }


def run(
        corpus: Path = CORPUS, output: Path = OUT, ridge: float = 100.0,
        max_pairs_per_group: int | None = 100,
        bootstrap_draws: int = 9999, seed: int = 20260812) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    train_pairs, train_meta = build_pairs(
        corpus, split="train", max_pairs_per_group=max_pairs_per_group)
    development_pairs, development_meta = build_pairs(
        corpus, split="development", max_pairs_per_group=max_pairs_per_group)
    if len(train_pairs) < 2 or len(development_pairs) < 2:
        raise ValueError("train and development pair sets are required")

    y_train = np.asarray([pair["delta_pK"] for pair in train_pairs], dtype=np.float64)
    y_dev = np.asarray([pair["delta_pK"] for pair in development_pairs], dtype=np.float64)
    wrong_targets, shuffled, proteins = protein_controls(corpus, seed=seed)

    models = {
        "L": fit_positive_ridge_no_intercept(stack(train_pairs, "ligand_delta"), y_train, ridge),
        "P": fit_positive_ridge_no_intercept(stack(train_pairs, "protein"), y_train, ridge),
        "A": fit_positive_ridge_no_intercept(stack(train_pairs, "concat"), y_train, ridge),
        "I": fit_positive_ridge_no_intercept(stack(train_pairs, "interaction"), y_train, ridge),
    }
    dev_features = {
        "L": stack(development_pairs, "ligand_delta"),
        "P": stack(development_pairs, "protein"),
        "A": stack(development_pairs, "concat"),
        "I": stack(development_pairs, "interaction"),
        "IW": controlled_interactions(development_pairs, wrong_targets, proteins),
        "IS": controlled_interactions(development_pairs, shuffled, proteins),
    }
    predictions = {
        "Z": np.zeros_like(y_dev),
        "L": predict(models["L"], dev_features["L"]),
        "P": predict(models["P"], dev_features["P"]),
        "A": predict(models["A"], dev_features["A"]),
        "I": predict(models["I"], dev_features["I"]),
        "IW": predict(models["I"], dev_features["IW"]),
        "IS": predict(models["I"], dev_features["IS"]),
    }
    rows = []
    for index, (pair, true) in enumerate(zip(development_pairs, y_dev)):
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
        component_contrast(rows, "I", "L", draws=bootstrap_draws, seed=seed),
        component_contrast(rows, "I", "IW", draws=bootstrap_draws, seed=seed + 1),
        component_contrast(rows, "I", "IS", draws=bootstrap_draws, seed=seed + 2),
        component_contrast(rows, "A", "L", draws=bootstrap_draws, seed=seed + 3),
        component_contrast(rows, "A", "Z", draws=bootstrap_draws, seed=seed + 4),
        component_contrast(rows, "L", "Z", draws=bootstrap_draws, seed=seed + 5),
    ]
    reversed_features = {
        "L": -dev_features["L"],
        "P": dev_features["P"],
        "A": np.concatenate([dev_features["P"], -dev_features["L"]], axis=1),
        "I": -dev_features["I"],
    }
    antisymmetry = {
        arm: antisymmetry_audit(models[arm], dev_features[arm], reversed_features[arm])
        for arm in ("L", "P", "A", "I")
    }
    gates = {
        "development_components_ge_5": development_meta["pairs"] > 0 and len({
            row["dependency_component"] for row in rows}) >= 5,
        "interaction_beats_ligand_delta": contrasts[0]["pass"],
        "interaction_beats_wrong_target": contrasts[1]["pass"],
        "interaction_beats_shuffled_target": contrasts[2]["pass"],
        "interaction_antisymmetry": antisymmetry["I"]["max_abs_sum"] <= 1e-9,
        "ligand_delta_antisymmetry": antisymmetry["L"]["max_abs_sum"] <= 1e-9,
    }
    verdict = (
        "BINDINGDB_SARDELTA_ATTRIBUTION_GATE1_PASS_TARGET_CONDITIONING"
        if all(gates.values())
        else "BINDINGDB_SARDELTA_ATTRIBUTION_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.BindingDBSARDeltaAttributionGate1.v1",
        "hypothesis": (
            "A true bilinear target x ligand-delta model explains BindingDB "
            "SAR-delta transfer beyond ligand-delta and wrong-target controls."),
        "corpus": {
            "manifest_sha256": sha256_file(corpus / "manifest.json"),
            "train_pairs": train_meta,
            "development_pairs": development_meta,
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
            for arm in ("L", "P", "A", "I")
        },
        "development_summary": {
            arm: summarize(rows, arm)
            for arm in ("Z", "L", "P", "A", "I", "IW", "IS")
        },
        "development_contrasts": contrasts,
        "antisymmetry_audit": antisymmetry,
        "gates": gates,
        "development_training_authorized": verdict.endswith("PASS_TARGET_CONDITIONING"),
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
