"""Train a source-only supervised slot localizer for crossed Ki quotients.

The localizer is learned only from train quotient blocks: each ESM2 residue
slot is scored by its standardized residual correlation with the correct-pair
quotient label after crossing with ligand E-state/pharmacophore chemistry.
The selected slots are then frozen before development scoring.
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
    read_jsonl,
    read_jsonl_gz,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
OUT = CQ_OUT.parent / "slot_localizer_cq_observable_gate1"


def protein_slot_blocks(
        residues: np.ndarray, mask: np.ndarray, *, hidden_blocks: int = 8) -> np.ndarray:
    if residues.ndim != 2:
        raise ValueError("residues must be a 2D array")
    if mask.ndim != 1 or mask.shape[0] != residues.shape[0]:
        raise ValueError("mask must match residue slots")
    if hidden_blocks <= 0:
        raise ValueError("hidden_blocks must be positive")
    if residues.shape[1] % hidden_blocks != 0:
        raise ValueError("hidden dimension must divide hidden_blocks")
    hidden_edges = np.linspace(0, residues.shape[1], hidden_blocks + 1, dtype=int)
    blocks = np.zeros((residues.shape[0], hidden_blocks), dtype=np.float64)
    for block_index, (start, end) in enumerate(zip(hidden_edges[:-1], hidden_edges[1:])):
        blocks[:, block_index] = residues[:, start:end].mean(axis=1)
    blocks[~mask] = 0.0
    return blocks


def _slot_products(
        ordered: list[dict], proteins: dict[str, np.ndarray],
        ligands: dict[str, np.ndarray], *, protein_key,
        ligand_key, slots: np.ndarray | None = None) -> np.ndarray:
    protein_values = np.stack([proteins[protein_key(cell)] for cell in ordered])
    if slots is not None:
        protein_values = protein_values[:, slots, :]
    ligand_values = np.stack([ligands[ligand_key(cell)] for cell in ordered])
    product = np.einsum("nsh,nl->nshl", protein_values, ligand_values, optimize=True)
    return product.reshape(product.shape[0], product.shape[1], -1)


def learn_slot_localizer(
        cells: dict[str, dict], panels: list[dict], proteins: dict[str, np.ndarray],
        ligands: dict[str, np.ndarray], *, top_slots: int,
        seed: int = 20260812, mode: str = "supervised") -> dict:
    if top_slots <= 0:
        raise ValueError("top_slots must be positive")
    slot_count = next(iter(proteins.values())).shape[0]
    if top_slots > slot_count:
        raise ValueError("top_slots cannot exceed available protein slots")
    if mode not in {"supervised", "uniform", "ligand_shuffled", "protein_only"}:
        raise ValueError(f"unknown localizer mode: {mode}")
    if mode == "uniform":
        selected = np.linspace(0, slot_count - 1, top_slots, dtype=int)
        return {
            "mode": mode,
            "selected_slots": selected.astype(int).tolist(),
            "slot_scores": [0.0] * slot_count,
            "train_panels_used": 0,
            "max_localizer_projection_orthogonality": 0.0,
        }

    score = np.zeros(slot_count, dtype=np.float64)
    panel_count = 0
    max_orthogonality = 0.0
    ligand_shuffle = {}
    if mode == "ligand_shuffled":
        ligand_keys = sorted(ligands)
        offset = int(np.random.default_rng(seed).integers(1, len(ligand_keys)))
        ligand_shuffle = {
            ligand: ligand_keys[(index + offset) % len(ligand_keys)]
            for index, ligand in enumerate(ligand_keys)
        }
    for panel in panels:
        if panel["split"] != "train":
            continue
        ordered = [cells[cell_id] for cell_id in panel["cell_ids"]]
        target_ids = [cell["target_id"] for cell in ordered]
        ligand_ids = [cell["ligand_id"] for cell in ordered]
        y_raw = np.asarray([cell["pK"] for cell in ordered], dtype=np.float64)
        y, retained_rank, y_orthogonality = additive_residual(target_ids, ligand_ids, y_raw)
        if retained_rank <= 0:
            continue
        if mode == "protein_only":
            raw = np.stack([proteins[cell["target_id"]] for cell in ordered])
        else:
            raw = _slot_products(
                ordered, proteins, ligands,
                protein_key=lambda cell: cell["target_id"],
                ligand_key=lambda cell: ligand_shuffle.get(
                    cell["ligand_id"], cell["ligand_id"]))
        flat = raw.reshape(raw.shape[0], -1)
        residual, feature_rank, feature_orthogonality = additive_residual(
            target_ids, ligand_ids, flat)
        if feature_rank != retained_rank:
            raise ValueError(f"feature rank mismatch in panel {panel['panel_id']}")
        residual = residual.reshape(raw.shape)
        y_energy = float(np.dot(y, y))
        if y_energy > 1e-12:
            x_energy = np.square(residual).sum(axis=0)
            xty = np.einsum("nsd,n->sd", residual, y, optimize=True)
            score += np.sum(np.square(xty) / (x_energy * y_energy + 1e-12), axis=1)
        max_orthogonality = max(max_orthogonality, y_orthogonality, feature_orthogonality)
        panel_count += 1
    if panel_count == 0:
        raise ValueError("no train quotient panels available for localizer")
    tie_breaker = np.random.default_rng(seed).uniform(0.0, 1e-12, size=slot_count)
    selected = np.argsort(-(score + tie_breaker))[:top_slots]
    selected.sort()
    return {
        "mode": mode,
        "selected_slots": selected.astype(int).tolist(),
        "slot_scores": score.tolist(),
        "train_panels_used": panel_count,
        "max_localizer_projection_orthogonality": max_orthogonality,
    }


def materialize_inputs(
        corpus: Path, protein_bank: Path, *,
        hidden_blocks: int) -> tuple[
            dict[str, dict], list[dict], dict[str, np.ndarray], dict[str, np.ndarray], dict]:
    cells_list = read_jsonl_gz(corpus / "cells.jsonl.gz")
    cells = {cell["cell_id"]: cell for cell in cells_list}
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
    proteins_json = read_jsonl(corpus / "proteins.jsonl")
    protein_keys = {row["sequence_sha256"] for row in proteins_json}
    protein_rows, protein_manifest = load_protein_bank(protein_bank, protein_keys)
    proteins = {
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
        "proteins": len(proteins),
        "ligands": len(ligands),
        "protein_bank_model_id": protein_manifest.get("model_id", ""),
        "protein_bank_model_revision": protein_manifest.get("model_revision", ""),
        "protein_bank_slot_policy": protein_manifest.get("slot_policy", ""),
    }
    return cells, panels, proteins, ligands, metadata


def materialize_features(
        cells: dict[str, dict], proteins: dict[str, np.ndarray],
        ligands: dict[str, np.ndarray], selected_slots: list[int]) -> dict[str, dict[str, np.ndarray]]:
    protein_donor, ligand_donor = donor_maps(list(cells.values()))
    slots = np.asarray(selected_slots, dtype=int)
    features = {}
    for cell_id, cell in cells.items():
        singleton = [cell]
        target = cell["target_id"]
        ligand = cell["ligand_id"]
        if target not in proteins or ligand not in ligands:
            raise ValueError("cell references missing protein or ligand")
        features[cell_id] = {
            "correct": _slot_products(
                singleton, proteins, ligands,
                protein_key=lambda row: row["target_id"],
                ligand_key=lambda row: row["ligand_id"],
                slots=slots).reshape(-1),
            "deranged_protein": _slot_products(
                singleton, proteins, ligands,
                protein_key=lambda row: protein_donor[row["target_id"]],
                ligand_key=lambda row: row["ligand_id"],
                slots=slots).reshape(-1),
            "foreign_ligand": _slot_products(
                singleton, proteins, ligands,
                protein_key=lambda row: row["target_id"],
                ligand_key=lambda row: ligand_donor[row["ligand_id"]],
                slots=slots).reshape(-1),
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
        hidden_blocks: int = 8, top_slots: int = 8,
        localizer_mode: str = "supervised") -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    cells, panels, proteins, ligands, input_metadata = materialize_inputs(
        corpus, protein_bank, hidden_blocks=hidden_blocks)
    localizer = learn_slot_localizer(
        cells, panels, proteins, ligands, top_slots=top_slots, seed=seed,
        mode=localizer_mode)
    features = materialize_features(cells, proteins, ligands, localizer["selected_slots"])
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
        "localizer_projection_orthogonality": (
            localizer["max_localizer_projection_orthogonality"] <= 1e-7),
        "development_components_ge_5": len({
            block.dependency_component for block in development_blocks}) >= 5,
        "correct_beats_zero_additive": contrasts[0]["pass"],
        "correct_beats_deranged_protein": contrasts[1]["pass"],
        "correct_beats_foreign_ligand": contrasts[2]["pass"],
    }
    verdict = (
        "SLOT_LOCALIZER_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "SLOT_LOCALIZER_CQ_GATE1_FAIL_CLOSED"
    )
    ligand_dim = int(next(iter(ligands.values())).shape[0])
    result = {
        "schema": "MetaSieve.SlotLocalizerCQObservableGate1.v1",
        "hypothesis": (
            "A train-only supervised localizer over frozen ESM2 residue slots, "
            "crossed with ligand E-state pharmacophore chemistry, carries "
            "dependency-transferable quotient interaction signal."),
        "literature_mechanism": {
            "plm_dta": "PLM residue states provide trainable protein representations",
            "attention_localization": (
                "ligand-relevant target regions should be selected rather than "
                "whole-protein averaged"),
            "hodge_cycle_space": (
                "all slot scoring and ridge fitting are performed after additive "
                "target+ligand quotient projection"),
        },
        "corpus": {
            **input_metadata,
            "blocks": len(blocks),
            "feature_source": "train_only_esm2_slot_localizer_x_ligand_estate",
            "feature_dim": int(top_slots * hidden_blocks * ligand_dim),
            "protein_descriptor_dim": int(top_slots * hidden_blocks),
            "ligand_descriptor_dim": ligand_dim,
            "hidden_blocks": hidden_blocks,
            "top_slots": top_slots,
            "max_projection_orthogonality": max_projection_orthogonality,
            "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
            "protein_bank_manifest_sha256": sha256_file(protein_bank / "manifest.json"),
        },
        "localizer": localizer,
        "config": {
            "ridge": ridge,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
            "arms": list(ARMS),
            "localizer_split": "train",
            "localizer_mode": localizer_mode,
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
    parser.add_argument("--top-slots", type=int, default=8)
    parser.add_argument(
        "--localizer-mode",
        choices=("supervised", "uniform", "ligand_shuffled", "protein_only"),
        default="supervised")
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, protein_bank=args.protein_bank, output=args.output,
        ridge=args.ridge, bootstrap_draws=args.bootstrap_draws, seed=args.seed,
        hidden_blocks=args.hidden_blocks, top_slots=args.top_slots,
        localizer_mode=args.localizer_mode)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
