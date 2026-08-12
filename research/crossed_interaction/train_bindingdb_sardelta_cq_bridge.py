"""Bridge target-conditioned SAR-delta supervision to BindingDB.

F-153 tests whether the F-152 source-only mechanism survives on the governed
BindingDB corpus: train on train-split same-target same-scaffold ligand-pair
deltas and score held-out development components against zero-delta.
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
from research.crossed_interaction.train_cq_observable import OUT as CQ_OUT
from research.crossed_interaction.train_seqchem_cq_observable import (
    ligand_descriptor,
    protein_descriptor,
    read_jsonl,
    read_jsonl_gz,
)
from research.source_affinity.train_chembl_assay_sardelta import (
    PAIR_SIMILARITY,
    fit_ridge,
    pair_feature,
    predict,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
OUT = CQ_OUT.parent / "bindingdb_sardelta_cq_bridge_gate1"


def build_cell_pairs(
        corpus: Path, *, split: str, feature_mode: str,
        max_pairs_per_group: int | None = None) -> tuple[list[dict], dict]:
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
                group_pairs.append({
                    "target_id": target,
                    "scaffold": scaffold,
                    "dependency_component": cell_component[left["cell_id"]],
                    "left_cell_id": left["cell_id"],
                    "right_cell_id": right["cell_id"],
                    "delta_pK": float(left["pK"] - right["pK"]),
                    "feature": pair_feature(
                        left_ligand["descriptor"],
                        right_ligand["descriptor"],
                        target=proteins[target],
                        mode=feature_mode),
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


def _summary(rows: list[dict], field: str) -> dict:
    return {
        "pairs": len(rows),
        "mse": float(np.mean([row[field] for row in rows])),
    }


def component_contrast(rows: list[dict], *, draws: int, seed: int) -> dict:
    correct: dict[str, list[float]] = defaultdict(list)
    zero: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        component = row["dependency_component"]
        correct[component].append(row["squared_error"])
        zero[component].append(row["zero_squared_error"])
    components = sorted(set(correct) & set(zero))
    if len(components) < 2:
        raise ValueError("component bootstrap needs at least two components")
    delta = np.asarray([
        float(np.mean(zero[component]) - np.mean(correct[component]))
        for component in components
    ], dtype=np.float64)
    rng = np.random.default_rng(seed)
    samples = delta[rng.integers(0, len(delta), size=(draws, len(delta)))].mean(axis=1)
    return {
        "components": len(components),
        "control": "zero_delta",
        "component_macro_reduction": float(delta.mean()),
        "one_sided_95_lcb": float(np.quantile(samples, 0.05)),
        "pass": bool(np.quantile(samples, 0.05) > 0.0),
    }


def run(
        corpus: Path = CORPUS, output: Path = OUT, ridge: float = 100.0,
        feature_mode: str = "delta_target", max_pairs_per_group: int | None = 200,
        bootstrap_draws: int = 9999, seed: int = 20260812) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    train_pairs, train_meta = build_cell_pairs(
        corpus, split="train", feature_mode=feature_mode,
        max_pairs_per_group=max_pairs_per_group)
    development_pairs, development_meta = build_cell_pairs(
        corpus, split="development", feature_mode=feature_mode,
        max_pairs_per_group=max_pairs_per_group)
    if len(train_pairs) < 2 or len(development_pairs) < 2:
        raise ValueError("train and development pair sets are required")
    x_train = np.stack([row["feature"] for row in train_pairs]).astype(np.float64)
    y_train = np.asarray([row["delta_pK"] for row in train_pairs], dtype=np.float64)
    model = fit_ridge(x_train, y_train, ridge)
    x_dev = np.stack([row["feature"] for row in development_pairs]).astype(np.float64)
    y_dev = np.asarray([row["delta_pK"] for row in development_pairs], dtype=np.float64)
    pred = predict(model, x_dev)
    rows = []
    for pair, true, estimate in zip(development_pairs, y_dev, pred):
        rows.append({
            "dependency_component": pair["dependency_component"],
            "target_id": pair["target_id"],
            "scaffold": pair["scaffold"],
            "left_cell_id": pair["left_cell_id"],
            "right_cell_id": pair["right_cell_id"],
            "delta_pK": float(true),
            "prediction": float(estimate),
            "squared_error": float((true - estimate) ** 2),
            "zero_squared_error": float(true ** 2),
        })
    contrast = component_contrast(rows, draws=bootstrap_draws, seed=seed)
    gates = {
        "development_components_ge_5": contrast["components"] >= 5,
        "correct_beats_zero_delta": contrast["pass"],
    }
    verdict = (
        "BINDINGDB_SARDELTA_CQ_BRIDGE_GATE1_PASS"
        if all(gates.values())
        else "BINDINGDB_SARDELTA_CQ_BRIDGE_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.BindingDBSARDeltaCQBridgeGate1.v1",
        "hypothesis": (
            "The target-conditioned SAR-delta mechanism admitted on ChEMBL "
            "transfers to governed BindingDB same-target same-scaffold deltas."),
        "corpus": {
            "manifest_sha256": sha256_file(corpus / "manifest.json"),
            "train_pairs": train_meta,
            "development_pairs": development_meta,
        },
        "config": {
            "ridge": ridge,
            "feature_mode": feature_mode,
            "max_pairs_per_group": max_pairs_per_group,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
        },
        "train_summary": {
            "pairs": len(train_pairs),
            "ridge_train_mse": model["train_mse"],
            "feature_dim": model["feature_dim"],
        },
        "development_summary": {
            "correct": _summary(rows, "squared_error"),
            "zero": _summary(rows, "zero_squared_error"),
        },
        "development_contrast": contrast,
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
    parser.add_argument("--feature-mode", choices=("delta", "delta_target", "concat"), default="delta_target")
    parser.add_argument("--max-pairs-per-group", type=int, default=200)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, output=args.output, ridge=args.ridge,
        feature_mode=args.feature_mode, max_pairs_per_group=args.max_pairs_per_group,
        bootstrap_draws=args.bootstrap_draws, seed=args.seed)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
