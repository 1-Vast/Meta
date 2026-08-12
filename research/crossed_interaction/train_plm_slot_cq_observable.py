"""Train an ESM slot-region observable on crossed Ki quotients.

This source-only Gate tests whether existing frozen protein-language-model
residue slots, summarized into deterministic sequence regions, carry a
partner-specific quotient interaction signal when crossed with ligand
E-state/pharmacophore chemistry. It does not use the failed T-BASIS bridge.
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
from research.crossed_interaction.train_seqchem_cq_observable import (
    donor_maps,
    product_feature,
    read_jsonl,
    read_jsonl_gz,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
PROTEIN_BANK = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_protein_bank"
OUT = CQ_OUT.parent / "plm_slot_cq_observable_gate1"


def load_protein_bank(path: Path, required: set[str]) -> tuple[dict[str, dict[str, np.ndarray]], dict]:
    manifest = json.loads((path / "manifest.json").read_text(encoding="utf-8"))
    rows = {}
    for shard in manifest["shards"]:
        with np.load(path / shard["path"], allow_pickle=False) as stored:
            keys = stored["keys"].tolist()
            residues = stored["residues"]
            masks = stored["mask"]
            for index, key in enumerate(keys):
                if key in required:
                    rows[key] = {
                        "residues": residues[index].astype(np.float64),
                        "mask": masks[index].astype(bool),
                    }
    missing = sorted(required - set(rows))
    if missing:
        raise ValueError(f"protein bank is missing {len(missing)} required proteins")
    return rows, manifest


def protein_slot_descriptor(
        residues: np.ndarray, mask: np.ndarray, *,
        slot_segments: int = 4, hidden_blocks: int = 8) -> np.ndarray:
    if residues.ndim != 2:
        raise ValueError("residues must be a 2D array")
    if mask.ndim != 1 or mask.shape[0] != residues.shape[0]:
        raise ValueError("mask must match residue slots")
    if slot_segments <= 0 or hidden_blocks <= 0:
        raise ValueError("slot_segments and hidden_blocks must be positive")
    if residues.shape[1] % hidden_blocks != 0:
        raise ValueError("hidden dimension must divide hidden_blocks")

    slot_edges = np.linspace(0, residues.shape[0], slot_segments + 1, dtype=int)
    hidden_edges = np.linspace(0, residues.shape[1], hidden_blocks + 1, dtype=int)
    features = []
    for slot_start, slot_end in zip(slot_edges[:-1], slot_edges[1:]):
        segment_mask = mask[slot_start:slot_end]
        segment = residues[slot_start:slot_end]
        if not np.any(segment_mask):
            features.extend([0.0] * hidden_blocks)
            continue
        active = segment[segment_mask]
        for hidden_start, hidden_end in zip(hidden_edges[:-1], hidden_edges[1:]):
            features.append(float(active[:, hidden_start:hidden_end].mean()))
    return np.asarray(features, dtype=np.float64)


def materialize_features(
        corpus: Path, protein_bank: Path, *,
        slot_segments: int, hidden_blocks: int) -> tuple[dict[str, dict[str, np.ndarray]], dict]:
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    proteins_json = read_jsonl(corpus / "proteins.jsonl")
    protein_keys = {row["sequence_sha256"] for row in proteins_json}
    protein_rows, protein_manifest = load_protein_bank(protein_bank, protein_keys)
    proteins = {
        key: protein_slot_descriptor(
            row["residues"], row["mask"],
            slot_segments=slot_segments, hidden_blocks=hidden_blocks)
        for key, row in protein_rows.items()
    }
    ligands = {
        row["drug_key"]: ligand_pharmacophore_descriptor(row["smiles"])
        for row in read_jsonl(corpus / "ligands.jsonl")
    }
    protein_donor, ligand_donor = donor_maps(cells)
    features = {}
    for cell in cells:
        target = cell["target_id"]
        ligand = cell["ligand_id"]
        features[cell["cell_id"]] = {
            "correct": product_feature(proteins[target], ligands[ligand]),
            "deranged_protein": product_feature(
                proteins[protein_donor[target]], ligands[ligand]),
            "foreign_ligand": product_feature(
                proteins[target], ligands[ligand_donor[ligand]]),
        }

    first = next(iter(features.values()))["correct"]
    metadata = {
        "feature_source": "esm2_slot_region_x_ligand_estate_pharmacophore",
        "feature_dim": int(first.shape[0]),
        "protein_descriptor_dim": int(next(iter(proteins.values())).shape[0]),
        "ligand_descriptor_dim": int(next(iter(ligands.values())).shape[0]),
        "slot_segments": slot_segments,
        "hidden_blocks": hidden_blocks,
        "cells": len(cells),
        "proteins": len(proteins),
        "ligands": len(ligands),
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
        "protein_bank_manifest_sha256": sha256_file(protein_bank / "manifest.json"),
        "protein_bank_model_id": protein_manifest.get("model_id", ""),
        "protein_bank_model_revision": protein_manifest.get("model_revision", ""),
        "protein_bank_slot_policy": protein_manifest.get("slot_policy", ""),
    }
    return features, metadata


def load_blocks(
        corpus: Path, protein_bank: Path, *,
        slot_segments: int, hidden_blocks: int) -> tuple[list[QuotientBlock], dict]:
    cells = {row["cell_id"]: row for row in read_jsonl_gz(corpus / "cells.jsonl.gz")}
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
    features, metadata = materialize_features(
        corpus, protein_bank,
        slot_segments=slot_segments, hidden_blocks=hidden_blocks)
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
    metadata["blocks"] = len(blocks)
    metadata["max_projection_orthogonality"] = max_orthogonality
    return blocks, metadata


def run(
        corpus: Path = CORPUS, protein_bank: Path = PROTEIN_BANK,
        output: Path = OUT, ridge: float = 10000.0,
        bootstrap_draws: int = 9999, seed: int = 20260812,
        slot_segments: int = 4, hidden_blocks: int = 8) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    blocks, metadata = load_blocks(
        corpus, protein_bank,
        slot_segments=slot_segments, hidden_blocks=hidden_blocks)
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
        "PLM_SLOT_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "PLM_SLOT_CQ_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.PLMSlotCQObservableGate1.v1",
        "hypothesis": (
            "Frozen ESM2 slot-region descriptors crossed with ligand E-state "
            "pharmacophore chemistry carry dependency-transferable quotient "
            "interaction signal."),
        "literature_mechanism": {
            "plm_dta": (
                "protein language model embeddings can serve as structure-free "
                "protein representations for DTA"),
            "plm_binding_site": (
                "residue-level PLM states can localize ligand-relevant protein "
                "regions without holo structures"),
            "hodge_cycle_space": (
                "quotient scoring removes target and ligand main effects"),
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
    parser.add_argument("--protein-bank", type=Path, default=PROTEIN_BANK)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ridge", type=float, default=10000.0)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--slot-segments", type=int, default=4)
    parser.add_argument("--hidden-blocks", type=int, default=8)
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, protein_bank=args.protein_bank, output=args.output,
        ridge=args.ridge, bootstrap_draws=args.bootstrap_draws, seed=args.seed,
        slot_segments=args.slot_segments, hidden_blocks=args.hidden_blocks)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
