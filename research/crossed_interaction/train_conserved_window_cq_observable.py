"""Train a family-conserved ESM2 window observable on crossed Ki quotients.

Family conserved windows are selected only from train-panel targets using
within-family ESM2 slot stability. Development targets are transformed with the
frozen train-derived windows and centroids before the shared quotient Gate.
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
from research.crossed_interaction.train_family_context_cq_observable import (
    train_targets_from_panels,
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
OUT = CQ_OUT.parent / "conserved_window_cq_observable_gate1"


@dataclass(frozen=True)
class ConservedWindowContext:
    target_group: dict[str, str]
    train_targets: set[str]
    global_slots: np.ndarray
    global_centroid: np.ndarray
    family_slots: dict[str, np.ndarray]
    family_centroids: dict[str, np.ndarray]
    top_windows: int
    min_train_family_size: int
    mode: str


def selected_window_descriptor(
        slot_blocks: np.ndarray, selected_slots: np.ndarray, centroid: np.ndarray,
        *, mode: str) -> np.ndarray:
    if mode == "raw_windows":
        return slot_blocks[selected_slots].reshape(-1)
    contrast = slot_blocks[selected_slots] - centroid
    return np.concatenate([contrast.reshape(-1), centroid.reshape(-1)])


def _top_stable_slots(vectors: np.ndarray, top_windows: int) -> tuple[np.ndarray, np.ndarray]:
    if vectors.ndim != 3:
        raise ValueError("family vectors must be target x slot x channel")
    centroid = vectors.mean(axis=0)
    stability = vectors.var(axis=0).mean(axis=1)
    selected = np.argsort(stability)[:top_windows]
    selected.sort()
    return selected.astype(int), centroid[selected].astype(np.float64)


def build_conserved_window_context(
        cells: dict[str, dict], panels: list[dict],
        protein_slots: dict[str, np.ndarray], *, top_windows: int = 4,
        min_train_family_size: int = 2,
        mode: str = "family_conserved") -> tuple[ConservedWindowContext, dict]:
    if top_windows <= 0:
        raise ValueError("top_windows must be positive")
    slot_count = next(iter(protein_slots.values())).shape[0]
    if top_windows > slot_count:
        raise ValueError("top_windows cannot exceed available slots")
    if min_train_family_size <= 0:
        raise ValueError("min_train_family_size must be positive")
    if mode not in {"family_conserved", "global_conserved", "raw_windows"}:
        raise ValueError(f"unknown conserved-window mode: {mode}")
    target_group = {cell["target_id"]: cell["protein_group_40"] for cell in cells.values()}
    train_targets = train_targets_from_panels(cells, panels)
    if not train_targets:
        raise ValueError("no train targets available for conserved windows")
    train_matrix = np.stack([protein_slots[target] for target in sorted(train_targets)])
    global_slots, global_centroid = _top_stable_slots(train_matrix, top_windows)

    by_group: dict[str, list[np.ndarray]] = {}
    for target in sorted(train_targets):
        by_group.setdefault(target_group[target], []).append(protein_slots[target])
    family_slots = {}
    family_centroids = {}
    for group, values in by_group.items():
        if len(values) < min_train_family_size:
            continue
        selected, centroid = _top_stable_slots(np.stack(values), top_windows)
        family_slots[group] = selected
        family_centroids[group] = centroid
    metadata = {
        "train_targets": len(train_targets),
        "target_groups": len(set(target_group.values())),
        "groups_with_conserved_windows": len(family_slots),
        "targets_with_family_windows": sum(
            1 for target, group in target_group.items() if group in family_slots),
        "targets_using_global_windows": sum(
            1 for target, group in target_group.items() if group not in family_slots),
        "global_slots": global_slots.astype(int).tolist(),
    }
    return ConservedWindowContext(
        target_group=target_group,
        train_targets=train_targets,
        global_slots=global_slots,
        global_centroid=global_centroid,
        family_slots=family_slots,
        family_centroids=family_centroids,
        top_windows=top_windows,
        min_train_family_size=min_train_family_size,
        mode=mode,
    ), metadata


def conserved_window_descriptor(
        target: str, protein_slots: dict[str, np.ndarray],
        context: ConservedWindowContext) -> np.ndarray:
    selected = context.global_slots
    centroid = context.global_centroid
    if context.mode == "family_conserved":
        group = context.target_group[target]
        selected = context.family_slots.get(group, context.global_slots)
        centroid = context.family_centroids.get(group, context.global_centroid)
    return selected_window_descriptor(
        protein_slots[target], selected, centroid, mode=context.mode)


def materialize_inputs(
        corpus: Path, protein_bank: Path, *,
        hidden_blocks: int) -> tuple[
            dict[str, dict], list[dict], dict[str, np.ndarray], dict[str, np.ndarray], dict]:
    cells_list = read_jsonl_gz(corpus / "cells.jsonl.gz")
    cells = {cell["cell_id"]: cell for cell in cells_list}
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
    protein_keys = {row["sequence_sha256"] for row in read_jsonl(corpus / "proteins.jsonl")}
    protein_rows, protein_manifest = load_protein_bank(protein_bank, protein_keys)
    protein_slots = {
        key: protein_slot_blocks(row["residues"], row["mask"], hidden_blocks=hidden_blocks)
        for key, row in protein_rows.items()
    }
    ligands = {
        row["drug_key"]: ligand_pharmacophore_descriptor(row["smiles"])
        for row in read_jsonl(corpus / "ligands.jsonl")
    }
    metadata = {
        "cells": len(cells_list),
        "panels": len(panels),
        "proteins": len(protein_slots),
        "ligands": len(ligands),
        "protein_bank_model_id": protein_manifest.get("model_id", ""),
        "protein_bank_model_revision": protein_manifest.get("model_revision", ""),
        "protein_bank_slot_policy": protein_manifest.get("slot_policy", ""),
    }
    return cells, panels, protein_slots, ligands, metadata


def materialize_features(
        cells: dict[str, dict], protein_slots: dict[str, np.ndarray],
        ligands: dict[str, np.ndarray],
        context: ConservedWindowContext) -> dict[str, dict[str, np.ndarray]]:
    protein_donor, ligand_donor = donor_maps(list(cells.values()))
    target_features = {
        target: conserved_window_descriptor(target, protein_slots, context)
        for target in protein_slots
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
        hidden_blocks: int = 8, top_windows: int = 4,
        min_train_family_size: int = 2,
        context_mode: str = "family_conserved") -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    cells, panels, protein_slots, ligands, input_metadata = materialize_inputs(
        corpus, protein_bank, hidden_blocks=hidden_blocks)
    context, context_metadata = build_conserved_window_context(
        cells, panels, protein_slots, top_windows=top_windows,
        min_train_family_size=min_train_family_size, mode=context_mode)
    features = materialize_features(cells, protein_slots, ligands, context)
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
        "CONSERVED_WINDOW_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "CONSERVED_WINDOW_CQ_GATE1_FAIL_CLOSED"
    )
    ligand_dim = int(next(iter(ligands.values())).shape[0])
    protein_dim = int(next(iter(features.values()))["correct"].shape[0] / ligand_dim)
    result = {
        "schema": "MetaSieve.ConservedWindowCQObservableGate1.v1",
        "hypothesis": (
            "Train-only family-conserved ESM2 windows provide a transferable "
            "local family/domain context for quotient interaction residuals."),
        "literature_mechanism": {
            "proteochemometrics": (
                "family/domain target descriptors and ligand-target cross terms "
                "are established PCM mechanisms"),
            "domain_binding_region_dti": (
                "domain or binding-region protein descriptors are more local "
                "than whole-sequence summaries"),
            "hodge_cycle_space": (
                "final scoring uses the same additive quotient and wrong-partner Gate"),
        },
        "corpus": {
            **input_metadata,
            "blocks": len(blocks),
            "feature_source": "train_only_family_conserved_esm2_windows_x_ligand_estate",
            "feature_dim": int(protein_dim * ligand_dim),
            "protein_descriptor_dim": protein_dim,
            "ligand_descriptor_dim": ligand_dim,
            "hidden_blocks": hidden_blocks,
            "top_windows": top_windows,
            "max_projection_orthogonality": max_projection_orthogonality,
            "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
            "protein_bank_manifest_sha256": sha256_file(protein_bank / "manifest.json"),
        },
        "conserved_window_context": {
            **context_metadata,
            "mode": context_mode,
            "window_train_split": "train",
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
    parser.add_argument("--top-windows", type=int, default=4)
    parser.add_argument("--min-train-family-size", type=int, default=2)
    parser.add_argument(
        "--context-mode",
        choices=("family_conserved", "global_conserved", "raw_windows"),
        default="family_conserved")
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, protein_bank=args.protein_bank, output=args.output,
        ridge=args.ridge, bootstrap_draws=args.bootstrap_draws, seed=args.seed,
        hidden_blocks=args.hidden_blocks, top_windows=args.top_windows,
        min_train_family_size=args.min_train_family_size,
        context_mode=args.context_mode)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
