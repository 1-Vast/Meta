"""Test a UniPert-inspired pair-score difference SAR-delta arm.

After symmetry Gate 0, this U1 diagnostic learns a low-capacity interaction
score s(P,L) from forward+reverse pairs and predicts
s(P,L_i) - s(P,L_j). The architecture is deliberately simple and closed-form:
the score feature is [ligand_descriptor(L); protein_descriptor(P) outer
ligand_descriptor(L)].
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from research.crossed_interaction.audit_bindingdb_sardelta_symmetry import (
    augment_reverse_pairs,
)
from research.crossed_interaction.bindingdb_cq_r0 import sha256_file
from research.crossed_interaction.train_bindingdb_sardelta_attribution import (
    CORPUS,
    OUT as ATTRIBUTION_OUT,
    build_pairs,
    component_contrast,
    controlled_interactions,
    fit_positive_ridge_no_intercept,
    interaction_feature,
    predict,
    protein_controls,
    stack,
    summarize,
)
from research.crossed_interaction.train_seqchem_cq_observable import product_feature
from research.source_affinity.train_chembl_assay_sardelta import PAIR_SIMILARITY


OUT = ATTRIBUTION_OUT.parent / "bindingdb_pair_score_difference_gate_u1"


def score_difference_feature(pair: dict, protein: np.ndarray | None = None) -> np.ndarray:
    descriptor = pair["protein"] if protein is None else protein
    left = np.concatenate([pair["left_ligand"], product_feature(descriptor, pair["left_ligand"])])
    right = np.concatenate([pair["right_ligand"], product_feature(descriptor, pair["right_ligand"])])
    return left - right


def build_pairs_with_endpoints(
        corpus: Path, *, split: str,
        max_pairs_per_group: int | None) -> tuple[list[dict], dict]:
    pairs, metadata = build_pairs(
        corpus, split=split, max_pairs_per_group=max_pairs_per_group)
    # build_pairs stores only the delta. Recover endpoints from delta-compatible
    # descriptors by extending it in-place through a local corpus lookup.
    from research.crossed_interaction.train_seqchem_cq_observable import (
        ligand_descriptor,
        read_jsonl,
    )

    ligands = {
        row["drug_key"]: ligand_descriptor(row["smiles"])
        for row in read_jsonl(corpus / "ligands.jsonl")
    }
    # Cell IDs do not encode ligand IDs in all corpora, so fill endpoint
    # descriptors by joining against cells.
    from research.crossed_interaction.train_seqchem_cq_observable import read_jsonl_gz

    cell_ligand = {
        row["cell_id"]: row["ligand_id"]
        for row in read_jsonl_gz(corpus / "cells.jsonl.gz")
    }
    for pair in pairs:
        pair["left_ligand"] = ligands[cell_ligand[pair["left_cell_id"]]]
        pair["right_ligand"] = ligands[cell_ligand[pair["right_cell_id"]]]
    return pairs, metadata


def score_difference_matrix(pairs: list[dict], target_map: dict[str, str] | None,
                            proteins: dict[str, np.ndarray] | None) -> np.ndarray:
    features = []
    for pair in pairs:
        protein = None if target_map is None else proteins[target_map[pair["target_id"]]]
        features.append(score_difference_feature(pair, protein=protein))
    return np.stack(features).astype(np.float64)


def run(
        corpus: Path = CORPUS, output: Path = OUT, ridge: float = 100.0,
        max_pairs_per_group: int | None = 100,
        bootstrap_draws: int = 9999, seed: int = 20260812) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    train_forward, train_meta = build_pairs_with_endpoints(
        corpus, split="train", max_pairs_per_group=max_pairs_per_group)
    dev_forward, dev_meta = build_pairs_with_endpoints(
        corpus, split="development", max_pairs_per_group=max_pairs_per_group)
    train_pairs = augment_reverse_pairs(train_forward)
    dev_pairs = augment_reverse_pairs(dev_forward)
    # Reverse augmentation flips delta and interaction; endpoint order also
    # needs to be reversed for score-difference features.
    for pair in train_pairs + dev_pairs:
        if "reversed_from" in pair:
            pair["left_ligand"], pair["right_ligand"] = pair["right_ligand"], pair["left_ligand"]
    y_train = np.asarray([pair["delta_pK"] for pair in train_pairs], dtype=np.float64)
    y_dev = np.asarray([pair["delta_pK"] for pair in dev_pairs], dtype=np.float64)
    wrong_targets, shuffled, proteins = protein_controls(corpus, seed=seed)
    models = {
        "L": fit_positive_ridge_no_intercept(stack(train_pairs, "ligand_delta"), y_train, ridge),
        "B": fit_positive_ridge_no_intercept(stack(train_pairs, "interaction"), y_train, ridge),
        "R": fit_positive_ridge_no_intercept(score_difference_matrix(train_pairs, None, None), y_train, ridge),
    }
    dev_features = {
        "L": stack(dev_pairs, "ligand_delta"),
        "B": stack(dev_pairs, "interaction"),
        "R": score_difference_matrix(dev_pairs, None, None),
        "RW": score_difference_matrix(dev_pairs, wrong_targets, proteins),
        "RS": score_difference_matrix(dev_pairs, shuffled, proteins),
        "BW": controlled_interactions(dev_pairs, wrong_targets, proteins),
        "BS": controlled_interactions(dev_pairs, shuffled, proteins),
    }
    predictions = {
        "Z": np.zeros_like(y_dev),
        "L": predict(models["L"], dev_features["L"]),
        "B": predict(models["B"], dev_features["B"]),
        "BW": predict(models["B"], dev_features["BW"]),
        "BS": predict(models["B"], dev_features["BS"]),
        "R": predict(models["R"], dev_features["R"]),
        "RW": predict(models["R"], dev_features["RW"]),
        "RS": predict(models["R"], dev_features["RS"]),
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
        component_contrast(rows, "R", "L", draws=bootstrap_draws, seed=seed),
        component_contrast(rows, "R", "RW", draws=bootstrap_draws, seed=seed + 1),
        component_contrast(rows, "R", "RS", draws=bootstrap_draws, seed=seed + 2),
        component_contrast(rows, "B", "L", draws=bootstrap_draws, seed=seed + 3),
        component_contrast(rows, "B", "BW", draws=bootstrap_draws, seed=seed + 4),
        component_contrast(rows, "B", "BS", draws=bootstrap_draws, seed=seed + 5),
    ]
    gates = {
        "development_components_ge_5": len({row["dependency_component"] for row in rows}) >= 5,
        "score_difference_beats_ligand_delta": contrasts[0]["pass"],
        "score_difference_beats_wrong_target": contrasts[1]["pass"],
        "score_difference_beats_shuffled_target": contrasts[2]["pass"],
    }
    verdict = (
        "BINDINGDB_PAIR_SCORE_DIFFERENCE_GATE_U1_PASS"
        if all(gates.values())
        else "BINDINGDB_PAIR_SCORE_DIFFERENCE_GATE_U1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.BindingDBPairScoreDifferenceGateU1.v1",
        "hypothesis": (
            "A UniPert-inspired score difference s(P,L_i)-s(P,L_j) identifies "
            "target-specific SAR-delta signal beyond ligand and wrong-target controls."),
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
            for arm in ("L", "B", "R")
        },
        "development_summary": {
            arm: summarize(rows, arm)
            for arm in ("Z", "L", "B", "BW", "BS", "R", "RW", "RS")
        },
        "development_contrasts": contrasts,
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
