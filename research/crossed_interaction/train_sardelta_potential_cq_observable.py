"""Train a SAR-delta potential observable for crossed Ki quotients.

F-154 lifts the successful target-conditioned pair-delta bridge into the
original panel CQ Gate. A train-split SAR-delta ridge predicts local ligand
differences; each cell is represented by the mean predicted delta against
same-target same-scaffold train neighbors. The final admission test is the
standard additive-residual CQ Gate with wrong-partner controls.
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

from research.crossed_interaction.bindingdb_cq_r0 import sha256_file
from research.crossed_interaction.train_bindingdb_sardelta_cq_bridge import (
    build_cell_pairs,
)
from research.crossed_interaction.train_cq_observable import (
    ARMS,
    OUT as CQ_OUT,
    QuotientBlock,
    additive_residual,
    bootstrap_contrast,
    fit_ridge,
    score_blocks,
)
from research.crossed_interaction.train_seqchem_cq_observable import (
    donor_maps,
    ligand_descriptor,
    protein_descriptor,
    read_jsonl,
    read_jsonl_gz,
)
from research.source_affinity.train_chembl_assay_sardelta import (
    fit_ridge as fit_delta_ridge,
    pair_feature,
    predict as predict_delta,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
OUT = CQ_OUT.parent / "sardelta_potential_cq_observable_gate1"


def train_delta_model(
        corpus: Path, *, ridge: float, feature_mode: str,
        max_pairs_per_group: int | None) -> tuple[dict, dict]:
    pairs, metadata = build_cell_pairs(
        corpus, split="train", feature_mode=feature_mode,
        max_pairs_per_group=max_pairs_per_group)
    x = np.stack([pair["feature"] for pair in pairs]).astype(np.float64)
    y = np.asarray([pair["delta_pK"] for pair in pairs], dtype=np.float64)
    model = fit_delta_ridge(x, y, ridge)
    metadata.update({
        "delta_model_ridge": ridge,
        "delta_model_train_mse": model["train_mse"],
        "delta_model_feature_dim": model["feature_dim"],
    })
    return model, metadata


def sardelta_potential(
        target_descriptor: np.ndarray, ligand: np.ndarray,
        neighbor_ligands: list[np.ndarray], model: dict, *,
        feature_mode: str) -> np.ndarray:
    if not neighbor_ligands:
        return np.zeros(1, dtype=np.float64)
    features = np.stack([
        pair_feature(ligand, neighbor, target=target_descriptor, mode=feature_mode)
        for neighbor in neighbor_ligands
    ]).astype(np.float64)
    return np.asarray([float(np.mean(predict_delta(model, features)))], dtype=np.float64)


def materialize_features(
        corpus: Path, model: dict, *, feature_mode: str) -> tuple[
            dict[str, dict[str, np.ndarray]], dict]:
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    proteins_json = read_jsonl(corpus / "proteins.jsonl")
    ligands_json = read_jsonl(corpus / "ligands.jsonl")
    proteins = {
        row["sequence_sha256"]: protein_descriptor(row["sequence"])
        for row in proteins_json
    }
    ligands = {
        row["drug_key"]: ligand_descriptor(row["smiles"])
        for row in ligands_json
    }
    train_neighbors: dict[tuple[str, str], list[np.ndarray]] = defaultdict(list)
    for cell in cells:
        if cell["split"] == "train":
            train_neighbors[(cell["target_id"], cell["scaffold"])].append(
                ligands[cell["ligand_id"]])
    protein_donor, ligand_donor = donor_maps(cells)
    features = {}
    for cell in cells:
        target = cell["target_id"]
        ligand = cell["ligand_id"]
        protein_control = protein_donor[target]
        ligand_control = ligand_donor[ligand]
        scaffold = cell["scaffold"]
        features[cell["cell_id"]] = {
            "correct": sardelta_potential(
                proteins[target], ligands[ligand],
                train_neighbors.get((target, scaffold), []),
                model, feature_mode=feature_mode),
            "deranged_protein": sardelta_potential(
                proteins[protein_control], ligands[ligand],
                train_neighbors.get((protein_control, scaffold), []),
                model, feature_mode=feature_mode),
            "foreign_ligand": sardelta_potential(
                proteins[target], ligands[ligand_control],
                train_neighbors.get((target, scaffold), []),
                model, feature_mode=feature_mode),
        }
    metadata = {
        "feature_source": "bindingdb_train_sardelta_neighbor_potential",
        "feature_mode": feature_mode,
        "feature_dim": 1,
        "cells": len(cells),
        "proteins": len(proteins),
        "ligands": len(ligands),
        "train_target_scaffold_neighbor_groups": len(train_neighbors),
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
    }
    return features, metadata


def load_blocks(
        corpus: Path, model: dict, *, feature_mode: str) -> tuple[list[QuotientBlock], dict]:
    cells = {row["cell_id"]: row for row in read_jsonl_gz(corpus / "cells.jsonl.gz")}
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
    features, metadata = materialize_features(corpus, model, feature_mode=feature_mode)
    blocks = []
    max_orthogonality = 0.0
    for panel in panels:
        ordered = [cells[cell_id] for cell_id in panel["cell_ids"]]
        target_ids = [cell["target_id"] for cell in ordered]
        ligand_ids = [cell["ligand_id"] for cell in ordered]
        y_raw = np.asarray([cell["pK"] for cell in ordered], dtype=np.float64)
        y, retained_rank, y_orthogonality = additive_residual(target_ids, ligand_ids, y_raw)
        if retained_rank <= 0:
            continue
        arm_features = {}
        for arm in ARMS:
            raw = np.stack([features[cell["cell_id"]][arm] for cell in ordered])
            residual, feature_rank, feature_orthogonality = additive_residual(
                target_ids, ligand_ids, raw)
            if residual.ndim == 1:
                residual = residual[:, None]
            if feature_rank != retained_rank:
                raise ValueError(f"feature rank mismatch in panel {panel['panel_id']}")
            arm_features[arm] = residual
            max_orthogonality = max(max_orthogonality, feature_orthogonality)
        max_orthogonality = max(max_orthogonality, y_orthogonality)
        blocks.append(QuotientBlock(
            panel_id=panel["panel_id"],
            split=panel["split"],
            dependency_component=panel["dependency_component"],
            retained_rank=retained_rank,
            y=np.asarray(y, dtype=np.float64),
            features=arm_features,
        ))
    metadata["blocks"] = len(blocks)
    metadata["max_projection_orthogonality"] = max_orthogonality
    return blocks, metadata


def run(
        corpus: Path = CORPUS, output: Path = OUT, ridge: float = 10000.0,
        delta_ridge: float = 100.0, feature_mode: str = "delta_target",
        max_pairs_per_group: int | None = 100,
        bootstrap_draws: int = 9999, seed: int = 20260812) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    delta_model, delta_metadata = train_delta_model(
        corpus, ridge=delta_ridge, feature_mode=feature_mode,
        max_pairs_per_group=max_pairs_per_group)
    blocks, metadata = load_blocks(corpus, delta_model, feature_mode=feature_mode)
    train_blocks = [block for block in blocks if block.split == "train"]
    development_blocks = [block for block in blocks if block.split == "development"]
    models = {arm: fit_ridge(train_blocks, arm, ridge) for arm in ARMS}
    train_rows, train_summary = score_blocks(train_blocks, models)
    development_rows, development_summary = score_blocks(development_blocks, models)
    controls = ("zero_additive", "deranged_protein", "foreign_ligand")
    contrasts = [
        bootstrap_contrast(
            development_rows, "correct", control,
            draws=bootstrap_draws, seed=seed + index)
        for index, control in enumerate(controls)
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
        "SARDELTA_POTENTIAL_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "SARDELTA_POTENTIAL_CQ_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.SARDeltaPotentialCQObservableGate1.v1",
        "hypothesis": (
            "A train-split target-conditioned SAR-delta field can be lifted to "
            "a pair-level potential that survives the original BindingDB CQ Gate."),
        "source_delta_model": delta_metadata,
        "corpus": metadata,
        "config": {
            "ridge": ridge,
            "delta_ridge": delta_ridge,
            "feature_mode": feature_mode,
            "max_pairs_per_group": max_pairs_per_group,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
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
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ridge", type=float, default=10000.0)
    parser.add_argument("--delta-ridge", type=float, default=100.0)
    parser.add_argument("--feature-mode", choices=("delta", "delta_target", "concat"), default="delta_target")
    parser.add_argument("--max-pairs-per-group", type=int, default=100)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, output=args.output, ridge=args.ridge,
        delta_ridge=args.delta_ridge, feature_mode=args.feature_mode,
        max_pairs_per_group=args.max_pairs_per_group,
        bootstrap_draws=args.bootstrap_draws, seed=args.seed)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
