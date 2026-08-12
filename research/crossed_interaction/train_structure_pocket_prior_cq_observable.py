"""Train an external-structure pocket prior for crossed Ki quotients.

F-146 is a source-only admission Gate. A pocket propensity prior is learned
from external BioLiP2/RCSB holo contact supervision and frozen before any
BindingDB quotient scoring. BindingDB affinity labels are used only by the
shared positive-ridge observable and the development Gate.
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
    product_feature,
    read_jsonl,
    read_jsonl_gz,
)
from research.crossed_interaction.train_slot_localizer_cq_observable import (
    protein_slot_blocks,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
STRUCTURE_SUPERVISION = (
    ROOT / "dataset/processed/open_structures/pilot20k_structure_supervision_v2"
)
STRUCTURE_PROTEIN_BANK = (
    ROOT / "dataset/processed/open_structures/pilot20k_esm2_t30_slots128_v1"
)
STRUCTURE_RECORDS = ROOT / "dataset/processed/open_structures/pilot20k_holo_governed_v2/complexes.jsonl"
OUT = CQ_OUT.parent / "structure_pocket_prior_cq_observable_gate1"


def _sigmoid(values: np.ndarray) -> np.ndarray:
    clipped = np.clip(values, -30.0, 30.0)
    return 1.0 / (1.0 + np.exp(-clipped))


def fit_contact_prior_from_arrays(blocks: np.ndarray, labels: np.ndarray, ridge: float) -> dict:
    if blocks.ndim != 2:
        raise ValueError("blocks must be a 2D matrix")
    if labels.ndim != 1 or labels.shape[0] != blocks.shape[0]:
        raise ValueError("labels must be a vector aligned to blocks")
    if ridge <= 0:
        raise ValueError("pocket prior ridge must be strictly positive")
    if blocks.shape[0] < 2:
        raise ValueError("pocket prior needs at least two slot rows")
    mean = blocks.mean(axis=0)
    scale = blocks.std(axis=0)
    scale[scale < 1e-6] = 1.0
    standardized = (blocks - mean) / scale
    design = np.column_stack([np.ones(blocks.shape[0], dtype=np.float64), standardized])
    penalty = np.eye(design.shape[1], dtype=np.float64)
    penalty[0, 0] = 0.0
    coefficients = np.linalg.solve(
        design.T @ design + ridge * penalty,
        design.T @ labels.astype(np.float64),
    )
    prediction = design @ coefficients
    return {
        "mean": mean,
        "scale": scale,
        "coefficients": coefficients,
        "ridge": float(ridge),
        "slot_rows": int(blocks.shape[0]),
        "positive_contact_rate": float(labels.mean()),
        "train_mse": float(np.square(labels - prediction).mean()),
    }


def contact_prior_features(
        blocks: np.ndarray, ligand: np.ndarray | None = None, *,
        prior_feature_mode: str = "protein_only") -> np.ndarray:
    if prior_feature_mode == "protein_only":
        return blocks
    if prior_feature_mode != "ligand_conditioned":
        raise ValueError(f"unknown prior feature mode: {prior_feature_mode}")
    if ligand is None:
        raise ValueError("ligand descriptor is required for ligand-conditioned prior")
    return np.einsum("sh,l->shl", blocks, ligand, optimize=True).reshape(blocks.shape[0], -1)


def predict_contact_scores(
        prior: dict, blocks: np.ndarray, ligand: np.ndarray | None = None) -> np.ndarray:
    features = contact_prior_features(
        blocks, ligand, prior_feature_mode=prior.get("prior_feature_mode", "protein_only"))
    standardized = (features - prior["mean"]) / prior["scale"]
    design = np.column_stack([np.ones(features.shape[0], dtype=np.float64), standardized])
    return design @ prior["coefficients"]


def pocket_weighted_descriptor(
        blocks: np.ndarray, mask: np.ndarray, prior: dict, *,
        ligand: np.ndarray | None = None,
        mode: str = "structure_prior", seed: int = 20260812) -> np.ndarray:
    if blocks.ndim != 2:
        raise ValueError("blocks must be a 2D slot matrix")
    if mask.ndim != 1 or mask.shape[0] != blocks.shape[0]:
        raise ValueError("mask must match slot blocks")
    if mode not in {"structure_prior", "uniform", "shuffled_prior"}:
        raise ValueError(f"unknown pocket mode: {mode}")
    active = mask.astype(bool)
    if not np.any(active):
        return np.zeros(blocks.shape[1], dtype=np.float64)
    if mode == "uniform":
        weights = active.astype(np.float64)
    else:
        propensity = _sigmoid(predict_contact_scores(prior, blocks, ligand))
        if mode == "shuffled_prior":
            offset = int(np.random.default_rng(seed).integers(1, len(propensity)))
            propensity = np.roll(propensity, offset)
        weights = propensity * active
    total = float(weights.sum())
    if total <= 1e-12:
        weights = active.astype(np.float64)
        total = float(weights.sum())
    return (weights[:, None] * blocks).sum(axis=0) / total


def _record_sequence_index(records_path: Path) -> dict[str, str]:
    records = read_jsonl(records_path)
    index = {}
    for record in records:
        entry = str(record["source_entry_id"])
        sequence = str(record["sequence_sha256"])
        if entry in index and index[entry] != sequence:
            raise ValueError(f"conflicting sequence for structure entry {entry}")
        index[entry] = sequence
    return index


def _record_sequence_ligand_index(records_path: Path) -> dict[str, tuple[str, np.ndarray]]:
    records = read_jsonl(records_path)
    index = {}
    for record in records:
        entry = str(record["source_entry_id"])
        try:
            ligand = ligand_pharmacophore_descriptor(str(record["canonical_smiles"]))
        except Exception:
            continue
        sequence = str(record["sequence_sha256"])
        index[entry] = (sequence, ligand)
    return index


def load_structure_contact_labels(
        supervision: Path, records_path: Path, *,
        max_records: int | None = None) -> tuple[dict[str, np.ndarray], dict]:
    sequence_by_entry = _record_sequence_index(records_path)
    pairs = read_jsonl(supervision / "pairs.jsonl")
    if max_records is not None:
        if max_records <= 0:
            raise ValueError("max_records must be positive when provided")
        pairs = pairs[:max_records]
    by_shard: dict[str, list[dict]] = {}
    for pair in pairs:
        by_shard.setdefault(str(pair["shard"]), []).append(pair)
    labels: dict[str, np.ndarray] = {}
    records_used = 0
    for shard, shard_pairs in sorted(by_shard.items()):
        with np.load(supervision / shard, allow_pickle=False) as stored:
            contact = stored["contact"]
            residue_mask = stored["residue_mask"].astype(bool)
            for pair in shard_pairs:
                entry = str(pair["source_entry_id"])
                sequence = sequence_by_entry.get(entry)
                if sequence is None:
                    continue
                index = int(pair["shard_index"])
                slot_label = np.any(contact[index].astype(bool), axis=0) & residue_mask[index]
                previous = labels.get(sequence)
                labels[sequence] = slot_label if previous is None else (previous | slot_label)
                records_used += 1
    metadata = {
        "structure_pairs_scanned": len(pairs),
        "structure_pairs_used": records_used,
        "structure_sequences_with_contact_labels": len(labels),
        "structure_supervision_manifest_sha256": sha256_file(supervision / "manifest.json"),
        "structure_records_sha256": sha256_file(records_path),
    }
    return labels, metadata


def load_structure_contact_examples(
        supervision: Path, records_path: Path, *,
        max_records: int | None = None) -> tuple[list[dict], dict]:
    sequence_ligand_by_entry = _record_sequence_ligand_index(records_path)
    pairs = read_jsonl(supervision / "pairs.jsonl")
    if max_records is not None:
        if max_records <= 0:
            raise ValueError("max_records must be positive when provided")
        pairs = pairs[:max_records]
    by_shard: dict[str, list[dict]] = {}
    for pair in pairs:
        by_shard.setdefault(str(pair["shard"]), []).append(pair)
    examples = []
    parsed_pairs = 0
    for shard, shard_pairs in sorted(by_shard.items()):
        with np.load(supervision / shard, allow_pickle=False) as stored:
            contact = stored["contact"]
            residue_mask = stored["residue_mask"].astype(bool)
            for pair in shard_pairs:
                entry = str(pair["source_entry_id"])
                sequence_ligand = sequence_ligand_by_entry.get(entry)
                if sequence_ligand is None:
                    continue
                sequence, ligand = sequence_ligand
                index = int(pair["shard_index"])
                slot_label = np.any(contact[index].astype(bool), axis=0) & residue_mask[index]
                examples.append({
                    "sequence": sequence,
                    "ligand": ligand,
                    "label": slot_label,
                })
                parsed_pairs += 1
    metadata = {
        "structure_pairs_scanned": len(pairs),
        "structure_pairs_used": parsed_pairs,
        "structure_pairs_skipped_ligand_parse": len(pairs) - parsed_pairs,
        "structure_sequences_with_contact_labels": len({
            example["sequence"] for example in examples}),
        "structure_supervision_manifest_sha256": sha256_file(supervision / "manifest.json"),
        "structure_records_sha256": sha256_file(records_path),
    }
    return examples, metadata


def fit_structure_pocket_prior(
        supervision: Path, records_path: Path, structure_protein_bank: Path, *,
        hidden_blocks: int, pocket_ridge: float, max_records: int | None,
        prior_feature_mode: str = "protein_only") -> tuple[dict, dict]:
    if prior_feature_mode == "protein_only":
        labels_by_sequence, metadata = load_structure_contact_labels(
            supervision, records_path, max_records=max_records)
        examples_by_sequence = {
            sequence: [{"ligand": None, "label": label}]
            for sequence, label in labels_by_sequence.items()
        }
    elif prior_feature_mode == "ligand_conditioned":
        examples, metadata = load_structure_contact_examples(
            supervision, records_path, max_records=max_records)
        examples_by_sequence: dict[str, list[dict]] = {}
        for example in examples:
            examples_by_sequence.setdefault(example["sequence"], []).append(example)
    else:
        raise ValueError(f"unknown prior feature mode: {prior_feature_mode}")
    required = set(examples_by_sequence)
    if not required:
        raise ValueError("no external structure contact labels are available")
    manifest = json.loads((structure_protein_bank / "manifest.json").read_text(encoding="utf-8"))
    block_rows = []
    label_rows = []
    sequences_used = 0
    for shard in manifest["shards"]:
        with np.load(structure_protein_bank / shard["path"], allow_pickle=False) as stored:
            keys = stored["keys"].tolist()
            residues = stored["residues"]
            masks = stored["mask"]
            for index, key in enumerate(keys):
                if key not in required:
                    continue
                blocks = protein_slot_blocks(
                    residues[index].astype(np.float64),
                    masks[index].astype(bool),
                    hidden_blocks=hidden_blocks,
                )
                active = masks[index].astype(bool)
                for example in examples_by_sequence[key]:
                    labels = example["label"].astype(np.float64)
                    features = contact_prior_features(
                        blocks, example["ligand"],
                        prior_feature_mode=prior_feature_mode)
                    block_rows.append(features[active])
                    label_rows.append(labels[active])
                sequences_used += 1
    if not block_rows:
        raise ValueError("structure protein bank had no rows for selected contact labels")
    blocks = np.concatenate(block_rows, axis=0)
    labels = np.concatenate(label_rows, axis=0)
    prior = fit_contact_prior_from_arrays(blocks, labels, pocket_ridge)
    prior["prior_feature_mode"] = prior_feature_mode
    metadata.update({
        "feature_source": "external_structure_supervised_esm2_pocket_prior",
        "hidden_blocks": hidden_blocks,
        "max_structure_records": max_records,
        "prior_feature_mode": prior_feature_mode,
        "structure_sequences_used_for_prior": sequences_used,
        "structure_slot_rows_used_for_prior": int(blocks.shape[0]),
        "structure_protein_bank_manifest_sha256": sha256_file(
            structure_protein_bank / "manifest.json"),
        "structure_protein_bank_model_id": manifest.get("model_id", ""),
        "structure_protein_bank_model_revision": manifest.get("model_revision", ""),
        "structure_protein_bank_slot_policy": manifest.get("slot_policy", ""),
        "pocket_prior_positive_contact_rate": prior["positive_contact_rate"],
        "pocket_prior_train_mse": prior["train_mse"],
    })
    return prior, metadata


def materialize_features(
        corpus: Path, protein_bank: Path, prior: dict, *,
        hidden_blocks: int, pocket_mode: str, seed: int) -> tuple[dict[str, dict[str, np.ndarray]], dict]:
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
        return pocket_weighted_descriptor(
            proteins[protein_key],
            protein_masks[protein_key],
            prior,
            ligand=ligands[ligand_key],
            mode=pocket_mode,
            seed=seed,
        )

    for cell in cells:
        target = cell["target_id"]
        ligand = cell["ligand_id"]
        protein_control = protein_donor[target]
        ligand_control = ligand_donor[ligand]
        correct_protein = descriptor(target, ligand)
        deranged_protein = descriptor(protein_control, ligand)
        foreign_ligand = descriptor(target, ligand_control)
        features[cell["cell_id"]] = {
            "correct": product_feature(correct_protein, ligands[ligand]),
            "deranged_protein": product_feature(
                deranged_protein, ligands[ligand]),
            "foreign_ligand": product_feature(
                foreign_ligand, ligands[ligand_control]),
        }
    first = next(iter(features.values()))["correct"]
    metadata = {
        "feature_dim": int(first.shape[0]),
        "protein_descriptor_dim": int(hidden_blocks),
        "ligand_descriptor_dim": int(next(iter(ligands.values())).shape[0]),
        "pocket_mode": pocket_mode,
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
        hidden_blocks: int = 8, pocket_mode: str = "structure_prior",
        prior_feature_mode: str = "protein_only") -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    prior, prior_metadata = fit_structure_pocket_prior(
        structure_supervision, structure_records, structure_protein_bank,
        hidden_blocks=hidden_blocks,
        pocket_ridge=pocket_ridge,
        max_records=max_structure_records,
        prior_feature_mode=prior_feature_mode,
    )
    features, feature_metadata = materialize_features(
        corpus, protein_bank, prior,
        hidden_blocks=hidden_blocks, pocket_mode=pocket_mode, seed=seed)
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
        "STRUCTURE_POCKET_PRIOR_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "STRUCTURE_POCKET_PRIOR_CQ_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.StructurePocketPriorCQObservableGate1.v1",
        "hypothesis": (
            "A pocket propensity prior trained only on external holo contact "
            "supervision can focus frozen ESM2 target states on ligand-relevant "
            "regions and yield transferable quotient interaction signal."),
        "literature_mechanism": {
            "binding_site_prior": (
                "structure-supervised binding-site predictors learn residue-level "
                "pocket propensity from holo ligand contacts"),
            "protein_language_model_slots": (
                "frozen PLM residue states can supply protein descriptors without "
                "BindingDB label leakage"),
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
            "pocket_mode": pocket_mode,
            "prior_feature_mode": prior_feature_mode,
            "arms": list(ARMS),
            "train_split": "train",
            "evaluation_split": "development",
        },
        "pocket_prior": {
            "ridge": prior["ridge"],
            "slot_rows": prior["slot_rows"],
            "positive_contact_rate": prior["positive_contact_rate"],
            "train_mse": prior["train_mse"],
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
    parser.add_argument(
        "--prior-feature-mode",
        choices=("protein_only", "ligand_conditioned"),
        default="protein_only",
    )
    parser.add_argument(
        "--pocket-mode",
        choices=("structure_prior", "uniform", "shuffled_prior"),
        default="structure_prior",
    )
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
        pocket_mode=args.pocket_mode,
        prior_feature_mode=args.prior_feature_mode,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
