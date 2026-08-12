"""Train a family-context ESM2 observable on crossed Ki quotients.

This source-only Gate changes the biological information source from whole
protein summaries to a train-only sequence-family context. Family centroids are
estimated only from train-panel targets, then frozen for development scoring.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from research.crossed_interaction.bindingdb_cq_r0 import sha256_file
from research.crossed_interaction.train_cq_observable import (
    ARMS,
    OUT as CQ_OUT,
    QuotientBlock,
    additive_residual,
    bootstrap_contrast,
    fit_ridge,
    score_blocks,
)
from research.crossed_interaction.train_physchem_cq_observable import (
    ligand_pharmacophore_descriptor,
)
from research.crossed_interaction.train_plm_slot_cq_observable import (
    PROTEIN_BANK,
    load_protein_bank,
)
from research.crossed_interaction.train_seqchem_cq_observable import (
    donor_maps,
    product_feature,
    read_jsonl,
    read_jsonl_gz,
)
from research.crossed_interaction.train_slot_localizer_cq_observable import (
    protein_slot_blocks,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
OUT = CQ_OUT.parent / "family_context_cq_observable_gate1"


@dataclass(frozen=True)
class FamilyContext:
    global_centroid: np.ndarray
    family_centroids: dict[str, np.ndarray]
    target_group: dict[str, str]
    train_targets: set[str]
    min_train_family_size: int
    mode: str


def pooled_slot_descriptor(blocks: np.ndarray, mask: np.ndarray) -> np.ndarray:
    if blocks.ndim != 2:
        raise ValueError("blocks must be 2D")
    if mask.ndim != 1 or mask.shape[0] != blocks.shape[0]:
        raise ValueError("mask must match slot axis")
    if np.any(mask):
        return blocks[mask].mean(axis=0)
    return np.zeros(blocks.shape[1], dtype=np.float64)


def train_targets_from_panels(cells: dict[str, dict], panels: list[dict]) -> set[str]:
    targets = set()
    for panel in panels:
        if panel["split"] != "train":
            continue
        for cell_id in panel["cell_ids"]:
            targets.add(cells[cell_id]["target_id"])
    return targets


def build_family_context(
        cells: dict[str, dict], panels: list[dict],
        protein_descriptors: dict[str, np.ndarray], *,
        min_train_family_size: int = 2,
        mode: str = "family_contrast") -> tuple[FamilyContext, dict]:
    if min_train_family_size <= 0:
        raise ValueError("min_train_family_size must be positive")
    if mode not in {"family_contrast", "global_contrast", "raw_pooled"}:
        raise ValueError(f"unknown context mode: {mode}")
    target_group = {}
    for cell in cells.values():
        target_group[cell["target_id"]] = cell["protein_group_40"]
    train_targets = train_targets_from_panels(cells, panels)
    train_vectors = [protein_descriptors[target] for target in sorted(train_targets)]
    if not train_vectors:
        raise ValueError("no train targets available for family context")
    global_centroid = np.mean(train_vectors, axis=0)
    by_group: dict[str, list[np.ndarray]] = {}
    for target in sorted(train_targets):
        by_group.setdefault(target_group[target], []).append(protein_descriptors[target])
    family_centroids = {
        group: np.mean(vectors, axis=0)
        for group, vectors in by_group.items()
        if len(vectors) >= min_train_family_size
    }
    metadata = {
        "train_targets": len(train_targets),
        "target_groups": len(set(target_group.values())),
        "groups_with_train_centroid": len(family_centroids),
        "targets_with_family_centroid": sum(
            1 for target, group in target_group.items() if group in family_centroids),
        "targets_using_global_centroid": sum(
            1 for target, group in target_group.items() if group not in family_centroids),
    }
    return FamilyContext(
        global_centroid=np.asarray(global_centroid, dtype=np.float64),
        family_centroids=family_centroids,
        target_group=target_group,
        train_targets=train_targets,
        min_train_family_size=min_train_family_size,
        mode=mode,
    ), metadata


def family_context_descriptor(
        target: str, protein_descriptors: dict[str, np.ndarray],
        context: FamilyContext) -> np.ndarray:
    pooled = protein_descriptors[target]
    if context.mode == "raw_pooled":
        return pooled
    centroid = context.global_centroid
    if context.mode == "family_contrast":
        centroid = context.family_centroids.get(
            context.target_group[target], context.global_centroid)
    contrast = pooled - centroid
    return np.concatenate([contrast, centroid])


def materialize_inputs(
        corpus: Path, protein_bank: Path, *,
        hidden_blocks: int) -> tuple[
            dict[str, dict], list[dict], dict[str, np.ndarray], dict[str, np.ndarray], dict]:
    cells_list = read_jsonl_gz(corpus / "cells.jsonl.gz")
    cells = {cell["cell_id"]: cell for cell in cells_list}
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
    protein_keys = {row["sequence_sha256"] for row in read_jsonl(corpus / "proteins.jsonl")}
    protein_rows, protein_manifest = load_protein_bank(protein_bank, protein_keys)
    protein_descriptors = {
        key: pooled_slot_descriptor(
            protein_slot_blocks(row["residues"], row["mask"], hidden_blocks=hidden_blocks),
            row["mask"].astype(bool))
        for key, row in protein_rows.items()
    }
    ligands = {
        row["drug_key"]: ligand_pharmacophore_descriptor(row["smiles"])
        for row in read_jsonl(corpus / "ligands.jsonl")
    }
    metadata = {
        "cells": len(cells_list),
        "panels": len(panels),
        "proteins": len(protein_descriptors),
        "ligands": len(ligands),
        "protein_bank_model_id": protein_manifest.get("model_id", ""),
        "protein_bank_model_revision": protein_manifest.get("model_revision", ""),
        "protein_bank_slot_policy": protein_manifest.get("slot_policy", ""),
    }
    return cells, panels, protein_descriptors, ligands, metadata


def materialize_features(
        cells: dict[str, dict], protein_descriptors: dict[str, np.ndarray],
        ligands: dict[str, np.ndarray], context: FamilyContext) -> dict[str, dict[str, np.ndarray]]:
    protein_donor, ligand_donor = donor_maps(list(cells.values()))
    target_features = {
        target: family_context_descriptor(target, protein_descriptors, context)
        for target in protein_descriptors
    }
    features = {}
    for cell_id, cell in cells.items():
        target = cell["target_id"]
        ligand = cell["ligand_id"]
        features[cell_id] = {
            "correct": product_feature(target_features[target], ligands[ligand]),
            "deranged_protein": product_feature(
                target_features[protein_donor[target]], ligands[ligand]),
            "foreign_ligand": product_feature(
                target_features[target], ligands[ligand_donor[ligand]]),
        }
    return features


def load_blocks(
        cells: dict[str, dict], panels: list[dict],
        features: dict[str, dict[str, np.ndarray]]) -> tuple[list[QuotientBlock], float]:
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
    return blocks, max_orthogonality


def run(
        corpus: Path = CORPUS, protein_bank: Path = PROTEIN_BANK,
        output: Path = OUT, ridge: float = 10000.0,
        bootstrap_draws: int = 9999, seed: int = 20260812,
        hidden_blocks: int = 8, min_train_family_size: int = 2,
        context_mode: str = "family_contrast") -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    cells, panels, protein_descriptors, ligands, input_metadata = materialize_inputs(
        corpus, protein_bank, hidden_blocks=hidden_blocks)
    context, context_metadata = build_family_context(
        cells, panels, protein_descriptors,
        min_train_family_size=min_train_family_size, mode=context_mode)
    features = materialize_features(cells, protein_descriptors, ligands, context)
    blocks, max_projection_orthogonality = load_blocks(cells, panels, features)
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
        "projection_orthogonality": max_projection_orthogonality <= 1e-7,
        "development_components_ge_5": len({
            block.dependency_component for block in development_blocks}) >= 5,
        "correct_beats_zero_additive": contrasts[0]["pass"],
        "correct_beats_deranged_protein": contrasts[1]["pass"],
        "correct_beats_foreign_ligand": contrasts[2]["pass"],
    }
    verdict = (
        "FAMILY_CONTEXT_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "FAMILY_CONTEXT_CQ_GATE1_FAIL_CLOSED"
    )
    ligand_dim = int(next(iter(ligands.values())).shape[0])
    protein_dim = int(next(iter(features.values()))["correct"].shape[0] / ligand_dim)
    result = {
        "schema": "MetaSieve.FamilyContextCQObservableGate1.v1",
        "hypothesis": (
            "Train-only sequence-family ESM2 centroids provide a transferable "
            "family/domain context for quotient interaction residuals."),
        "literature_mechanism": {
            "proteochemometrics": (
                "target-family context and compound descriptors can improve "
                "cross-target affinity modelling"),
            "homology_domain_surrogate": (
                "sequence-family grouping is a governed surrogate for domain "
                "context when pocket structures are unavailable"),
            "hodge_cycle_space": (
                "final scoring removes target and ligand main effects with the "
                "same additive quotient Gate"),
        },
        "corpus": {
            **input_metadata,
            "blocks": len(blocks),
            "feature_source": "train_only_family_context_esm2_x_ligand_estate",
            "feature_dim": int(protein_dim * ligand_dim),
            "protein_descriptor_dim": protein_dim,
            "ligand_descriptor_dim": ligand_dim,
            "hidden_blocks": hidden_blocks,
            "max_projection_orthogonality": max_projection_orthogonality,
            "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
            "protein_bank_manifest_sha256": sha256_file(protein_bank / "manifest.json"),
        },
        "family_context": {
            **context_metadata,
            "mode": context_mode,
            "centroid_train_split": "train",
            "min_train_family_size": min_train_family_size,
        },
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
    parser.add_argument("--protein-bank", type=Path, default=PROTEIN_BANK)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ridge", type=float, default=10000.0)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--hidden-blocks", type=int, default=8)
    parser.add_argument("--min-train-family-size", type=int, default=2)
    parser.add_argument(
        "--context-mode",
        choices=("family_contrast", "global_contrast", "raw_pooled"),
        default="family_contrast")
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, protein_bank=args.protein_bank, output=args.output,
        ridge=args.ridge, bootstrap_draws=args.bootstrap_draws, seed=args.seed,
        hidden_blocks=args.hidden_blocks,
        min_train_family_size=args.min_train_family_size,
        context_mode=args.context_mode)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
