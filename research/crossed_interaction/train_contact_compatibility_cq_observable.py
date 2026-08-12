"""Train an external contact-compatibility observable on crossed Ki quotients.

F-148 keeps the external BioLiP2/RCSB contact teacher from F-146 but changes
the BindingDB feature: the ligand-conditioned contact score distribution is
used directly as the pair observable instead of as a pooling weight.
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
from research.crossed_interaction.train_slot_localizer_cq_observable import (
    protein_slot_blocks,
)
from research.crossed_interaction.train_structure_pocket_prior_cq_observable import (
    STRUCTURE_PROTEIN_BANK,
    STRUCTURE_RECORDS,
    STRUCTURE_SUPERVISION,
    _sigmoid,
    fit_structure_pocket_prior,
    predict_contact_scores,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
OUT = CQ_OUT.parent / "contact_compatibility_cq_observable_gate1"


def contact_compatibility_descriptor(
        blocks: np.ndarray, mask: np.ndarray, ligand: np.ndarray, prior: dict) -> np.ndarray:
    if blocks.ndim != 2:
        raise ValueError("blocks must be a 2D slot matrix")
    if mask.ndim != 1 or mask.shape[0] != blocks.shape[0]:
        raise ValueError("mask must match slot blocks")
    active = mask.astype(bool)
    if not np.any(active):
        return np.zeros(17, dtype=np.float64)
    probabilities = _sigmoid(predict_contact_scores(prior, blocks, ligand))[active]
    positions = np.arange(len(active), dtype=np.float64)[active] / max(len(active) - 1, 1)
    ordered = np.sort(probabilities)[::-1]
    top_means = [
        float(ordered[:min(k, len(ordered))].mean())
        for k in (1, 3, 8, 16)
    ]
    total = float(probabilities.sum())
    if total <= 1e-12 or len(probabilities) == 1:
        entropy = 0.0
    else:
        weights = probabilities / total
        entropy = float(-np.sum(weights * np.log(weights + 1e-12)) / np.log(len(weights)))
    if total <= 1e-12:
        position_mean = 0.0
        position_std = 0.0
    else:
        weights = probabilities / total
        position_mean = float(np.sum(weights * positions))
        position_std = float(np.sqrt(np.sum(weights * np.square(positions - position_mean))))
    top_position = float(positions[int(np.argmax(probabilities))])
    mean = float(probabilities.mean())
    std = float(probabilities.std())
    quantiles = np.quantile(probabilities, [0.1, 0.25, 0.5, 0.75, 0.9])
    return np.asarray([
        mean,
        std,
        float(probabilities.min()),
        float(probabilities.max()),
        *top_means,
        *[float(value) for value in quantiles],
        entropy,
        position_mean,
        position_std,
        top_position,
    ], dtype=np.float64)


def materialize_features(
        corpus: Path, protein_bank: Path, prior: dict, *,
        hidden_blocks: int) -> tuple[dict[str, dict[str, np.ndarray]], dict]:
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    proteins_json = read_jsonl(corpus / "proteins.jsonl")
    protein_keys = {row["sequence_sha256"] for row in proteins_json}
    protein_rows, protein_manifest = load_protein_bank(protein_bank, protein_keys)
    proteins = {}
    protein_masks = {}
    for key, row in protein_rows.items():
        proteins[key] = protein_slot_blocks(
            row["residues"], row["mask"].astype(bool), hidden_blocks=hidden_blocks)
        protein_masks[key] = row["mask"].astype(bool)
    ligands = {
        row["drug_key"]: ligand_pharmacophore_descriptor(row["smiles"])
        for row in read_jsonl(corpus / "ligands.jsonl")
    }
    protein_donor, ligand_donor = donor_maps(cells)
    features = {}

    def descriptor(protein_key: str, ligand_key: str) -> np.ndarray:
        return contact_compatibility_descriptor(
            proteins[protein_key],
            protein_masks[protein_key],
            ligands[ligand_key],
            prior,
        )

    for cell in cells:
        target = cell["target_id"]
        ligand = cell["ligand_id"]
        features[cell["cell_id"]] = {
            "correct": descriptor(target, ligand),
            "deranged_protein": descriptor(protein_donor[target], ligand),
            "foreign_ligand": descriptor(target, ligand_donor[ligand]),
        }
    first = next(iter(features.values()))["correct"]
    metadata = {
        "feature_dim": int(first.shape[0]),
        "hidden_blocks": hidden_blocks,
        "cells": len(cells),
        "proteins": len(proteins),
        "ligands": len(ligands),
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
        "bindingdb_protein_bank_manifest_sha256": sha256_file(protein_bank / "manifest.json"),
        "bindingdb_protein_bank_model_id": protein_manifest.get("model_id", ""),
        "bindingdb_protein_bank_model_revision": protein_manifest.get("model_revision", ""),
        "bindingdb_protein_bank_slot_policy": protein_manifest.get("slot_policy", ""),
    }
    return features, metadata


def load_blocks(
        corpus: Path, features: dict[str, dict[str, np.ndarray]]) -> tuple[list[QuotientBlock], dict]:
    cells = {row["cell_id"]: row for row in read_jsonl_gz(corpus / "cells.jsonl.gz")}
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
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
    return blocks, {
        "blocks": len(blocks),
        "max_projection_orthogonality": max_orthogonality,
    }


def run(
        corpus: Path = CORPUS, protein_bank: Path = PROTEIN_BANK,
        structure_supervision: Path = STRUCTURE_SUPERVISION,
        structure_records: Path = STRUCTURE_RECORDS,
        structure_protein_bank: Path = STRUCTURE_PROTEIN_BANK,
        output: Path = OUT, ridge: float = 10000.0,
        pocket_ridge: float = 10.0, max_structure_records: int | None = 4096,
        bootstrap_draws: int = 9999, seed: int = 20260812,
        hidden_blocks: int = 8) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    prior, prior_metadata = fit_structure_pocket_prior(
        structure_supervision, structure_records, structure_protein_bank,
        hidden_blocks=hidden_blocks,
        pocket_ridge=pocket_ridge,
        max_records=max_structure_records,
        prior_feature_mode="ligand_conditioned",
    )
    features, feature_metadata = materialize_features(
        corpus, protein_bank, prior, hidden_blocks=hidden_blocks)
    blocks, block_metadata = load_blocks(corpus, features)
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
        "projection_orthogonality": block_metadata["max_projection_orthogonality"] <= 1e-7,
        "development_components_ge_5": len({
            block.dependency_component for block in development_blocks}) >= 5,
        "correct_beats_zero_additive": contrasts[0]["pass"],
        "correct_beats_deranged_protein": contrasts[1]["pass"],
        "correct_beats_foreign_ligand": contrasts[2]["pass"],
    }
    verdict = (
        "CONTACT_COMPATIBILITY_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "CONTACT_COMPATIBILITY_CQ_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.ContactCompatibilityCQObservableGate1.v1",
        "hypothesis": (
            "A frozen ligand-conditioned contact teacher trained on external "
            "holo structures yields transferable pair-specific compatibility "
            "statistics in BindingDB quotient space."),
        "literature_mechanism": {
            "ligand_aware_binding_site": (
                "LaMPSite/LABind-style mechanisms condition binding-residue "
                "prediction on both protein residue states and ligand chemistry"),
            "external_contact_teacher": (
                "BioLiP2/RCSB holo contacts provide label-free structure "
                "supervision independent of BindingDB affinity labels"),
            "quotient_gate": (
                "target and ligand main effects are removed before positive-ridge "
                "interaction scoring and wrong-partner controls"),
        },
        "corpus": prior_metadata | feature_metadata | block_metadata,
        "config": {
            "ridge": ridge,
            "pocket_ridge": pocket_ridge,
            "max_structure_records": max_structure_records,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
            "hidden_blocks": hidden_blocks,
            "arms": list(ARMS),
            "train_split": "train",
            "evaluation_split": "development",
        },
        "pocket_prior": {
            "ridge": prior["ridge"],
            "slot_rows": prior["slot_rows"],
            "positive_contact_rate": prior["positive_contact_rate"],
            "train_mse": prior["train_mse"],
            "prior_feature_mode": prior["prior_feature_mode"],
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--protein-bank", type=Path, default=PROTEIN_BANK)
    parser.add_argument("--structure-supervision", type=Path, default=STRUCTURE_SUPERVISION)
    parser.add_argument("--structure-records", type=Path, default=STRUCTURE_RECORDS)
    parser.add_argument("--structure-protein-bank", type=Path, default=STRUCTURE_PROTEIN_BANK)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--ridge", type=float, default=10000.0)
    parser.add_argument("--pocket-ridge", type=float, default=10.0)
    parser.add_argument("--max-structure-records", type=int, default=4096)
    parser.add_argument("--all-structure-records", action="store_true")
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--hidden-blocks", type=int, default=8)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    max_records = None if args.all_structure_records else args.max_structure_records
    result = run(
        corpus=args.corpus,
        protein_bank=args.protein_bank,
        structure_supervision=args.structure_supervision,
        structure_records=args.structure_records,
        structure_protein_bank=args.structure_protein_bank,
        output=args.output,
        ridge=args.ridge,
        pocket_ridge=args.pocket_ridge,
        max_structure_records=max_records,
        bootstrap_draws=args.bootstrap_draws,
        seed=args.seed,
        hidden_blocks=args.hidden_blocks,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
