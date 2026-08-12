"""Train an external pocket-prototype observable on crossed Ki quotients.

F-147 tests whether BioLiP2/RCSB holo pockets can define a frozen pocket-family
dictionary. BindingDB pairs are represented only by similarity to that external
dictionary before the shared positive-ridge quotient Gate is applied.
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
    _record_sequence_ligand_index,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
OUT = CQ_OUT.parent / "pocket_prototype_cq_observable_gate1"


def _unit_rows(values: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    norms[norms < 1e-12] = 1.0
    return values / norms


def _unit_vector(value: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(value))
    if norm < 1e-12:
        return np.zeros_like(value, dtype=np.float64)
    return value.astype(np.float64) / norm


def farthest_first_indices(features: np.ndarray, count: int, seed: int) -> np.ndarray:
    if features.ndim != 2:
        raise ValueError("features must be a 2D matrix")
    if count <= 0:
        raise ValueError("prototype count must be positive")
    if count >= features.shape[0]:
        return np.arange(features.shape[0], dtype=int)
    rng = np.random.default_rng(seed)
    selected = [int(rng.integers(0, features.shape[0]))]
    closest = 2.0 - 2.0 * (features @ features[selected[0]])
    for _ in range(1, count):
        index = int(np.argmax(closest))
        selected.append(index)
        distance = 2.0 - 2.0 * (features @ features[index])
        closest = np.minimum(closest, distance)
    return np.asarray(selected, dtype=int)


def load_structure_pocket_examples(
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
                labels = np.any(contact[index].astype(bool), axis=0) & residue_mask[index]
                if not np.any(labels):
                    continue
                examples.append({
                    "source_entry_id": entry,
                    "sequence": sequence,
                    "ligand": ligand,
                    "label": labels,
                })
                parsed_pairs += 1
    metadata = {
        "structure_pairs_scanned": len(pairs),
        "structure_pairs_used": parsed_pairs,
        "structure_pairs_skipped": len(pairs) - parsed_pairs,
        "structure_sequences_with_pockets": len({example["sequence"] for example in examples}),
        "structure_supervision_manifest_sha256": sha256_file(supervision / "manifest.json"),
        "structure_records_sha256": sha256_file(records_path),
    }
    return examples, metadata


def build_pocket_prototypes(
        supervision: Path, records_path: Path, structure_protein_bank: Path, *,
        hidden_blocks: int, prototype_count: int, max_records: int | None,
        ligand_weight: float, seed: int) -> tuple[dict, dict]:
    if hidden_blocks <= 0:
        raise ValueError("hidden_blocks must be positive")
    if ligand_weight < 0:
        raise ValueError("ligand_weight cannot be negative")
    examples, metadata = load_structure_pocket_examples(
        supervision, records_path, max_records=max_records)
    examples_by_sequence: dict[str, list[dict]] = {}
    for example in examples:
        examples_by_sequence.setdefault(example["sequence"], []).append(example)
    required = set(examples_by_sequence)
    if not required:
        raise ValueError("no external pocket examples available")
    manifest = json.loads((structure_protein_bank / "manifest.json").read_text(encoding="utf-8"))
    pocket_vectors = []
    ligand_vectors = []
    source_entries = []
    contact_slot_counts = []
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
                    contact = example["label"].astype(bool) & active
                    if not np.any(contact):
                        continue
                    pocket_vectors.append(blocks[contact].mean(axis=0))
                    ligand_vectors.append(example["ligand"])
                    source_entries.append(example["source_entry_id"])
                    contact_slot_counts.append(int(contact.sum()))
    if not pocket_vectors:
        raise ValueError("structure protein bank yielded no pocket vectors")
    pockets = np.stack(pocket_vectors).astype(np.float64)
    ligands = np.stack(ligand_vectors).astype(np.float64)
    selection_features = np.concatenate([
        _unit_rows(pockets),
        ligand_weight * _unit_rows(ligands),
    ], axis=1)
    selection_features = _unit_rows(selection_features)
    selected = farthest_first_indices(selection_features, prototype_count, seed)
    prototypes = {
        "pockets": pockets[selected],
        "ligands": ligands[selected],
        "source_entry_ids": [source_entries[int(index)] for index in selected],
        "contact_slot_counts": [contact_slot_counts[int(index)] for index in selected],
    }
    metadata.update({
        "feature_source": "external_holo_pocket_prototype_dictionary",
        "hidden_blocks": hidden_blocks,
        "prototype_count_requested": prototype_count,
        "prototype_count": int(len(selected)),
        "candidate_pockets": int(len(pockets)),
        "ligand_weight": ligand_weight,
        "max_structure_records": max_records,
        "mean_contact_slots_per_candidate": float(np.mean(contact_slot_counts)),
        "structure_protein_bank_manifest_sha256": sha256_file(
            structure_protein_bank / "manifest.json"),
        "structure_protein_bank_model_id": manifest.get("model_id", ""),
        "structure_protein_bank_model_revision": manifest.get("model_revision", ""),
        "structure_protein_bank_slot_policy": manifest.get("slot_policy", ""),
    })
    return prototypes, metadata


def protein_prototype_similarity(
        blocks: np.ndarray, mask: np.ndarray, prototype_pockets: np.ndarray) -> np.ndarray:
    active = mask.astype(bool)
    if not np.any(active):
        return np.zeros(prototype_pockets.shape[0], dtype=np.float64)
    slots = _unit_rows(blocks[active].astype(np.float64))
    prototypes = _unit_rows(prototype_pockets.astype(np.float64))
    return np.max(slots @ prototypes.T, axis=0)


def ligand_prototype_similarity(ligand: np.ndarray, prototype_ligands: np.ndarray) -> np.ndarray:
    return _unit_vector(ligand) @ _unit_rows(prototype_ligands).T


def materialize_features(
        corpus: Path, protein_bank: Path, prototypes: dict, *,
        hidden_blocks: int) -> tuple[dict[str, dict[str, np.ndarray]], dict]:
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    proteins_json = read_jsonl(corpus / "proteins.jsonl")
    protein_keys = {row["sequence_sha256"] for row in proteins_json}
    protein_rows, protein_manifest = load_protein_bank(protein_bank, protein_keys)
    protein_similarity = {}
    for key, row in protein_rows.items():
        blocks = protein_slot_blocks(
            row["residues"], row["mask"].astype(bool), hidden_blocks=hidden_blocks)
        protein_similarity[key] = protein_prototype_similarity(
            blocks, row["mask"].astype(bool), prototypes["pockets"])
    ligand_similarity_by_key = {
        row["drug_key"]: ligand_prototype_similarity(
            ligand_pharmacophore_descriptor(row["smiles"]),
            prototypes["ligands"])
        for row in read_jsonl(corpus / "ligands.jsonl")
    }
    protein_donor, ligand_donor = donor_maps(cells)
    features = {}
    for cell in cells:
        target = cell["target_id"]
        ligand = cell["ligand_id"]
        protein_control = protein_donor[target]
        ligand_control = ligand_donor[ligand]
        features[cell["cell_id"]] = {
            "correct": protein_similarity[target] * ligand_similarity_by_key[ligand],
            "deranged_protein": (
                protein_similarity[protein_control] * ligand_similarity_by_key[ligand]
            ),
            "foreign_ligand": (
                protein_similarity[target] * ligand_similarity_by_key[ligand_control]
            ),
        }
    first = next(iter(features.values()))["correct"]
    metadata = {
        "feature_dim": int(first.shape[0]),
        "protein_descriptor_dim": int(first.shape[0]),
        "ligand_descriptor_dim": int(first.shape[0]),
        "cells": len(cells),
        "proteins": len(protein_similarity),
        "ligands": len(ligand_similarity_by_key),
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
        max_structure_records: int | None = 4096,
        prototype_count: int = 64, ligand_weight: float = 1.0,
        bootstrap_draws: int = 9999, seed: int = 20260812,
        hidden_blocks: int = 8) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    prototypes, prototype_metadata = build_pocket_prototypes(
        structure_supervision, structure_records, structure_protein_bank,
        hidden_blocks=hidden_blocks,
        prototype_count=prototype_count,
        max_records=max_structure_records,
        ligand_weight=ligand_weight,
        seed=seed,
    )
    features, feature_metadata = materialize_features(
        corpus, protein_bank, prototypes, hidden_blocks=hidden_blocks)
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
        "POCKET_PROTOTYPE_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "POCKET_PROTOTYPE_CQ_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.PocketPrototypeCQObservableGate1.v1",
        "hypothesis": (
            "A frozen external holo pocket-family dictionary can convert ESM2 "
            "target slots and ligand descriptors into transferable quotient "
            "interaction coordinates."),
        "literature_mechanism": {
            "pocket_similarity": (
                "protein-ligand binding pockets can be compared as reusable "
                "local biological contexts rather than whole-protein summaries"),
            "prototype_learning": (
                "few-shot prototype mechanisms motivate support-free dictionary "
                "coordinates that are later admitted by support/query-isolated Gates"),
            "quotient_gate": (
                "target and ligand main effects are removed before positive-ridge "
                "interaction scoring and wrong-partner controls"),
        },
        "corpus": prototype_metadata | feature_metadata | block_metadata,
        "config": {
            "ridge": ridge,
            "max_structure_records": max_structure_records,
            "prototype_count": prototype_count,
            "ligand_weight": ligand_weight,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
            "hidden_blocks": hidden_blocks,
            "arms": list(ARMS),
            "train_split": "train",
            "evaluation_split": "development",
        },
        "prototype_examples": prototypes["source_entry_ids"][:16],
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
    parser.add_argument("--max-structure-records", type=int, default=4096)
    parser.add_argument("--all-structure-records", action="store_true")
    parser.add_argument("--prototype-count", type=int, default=64)
    parser.add_argument("--ligand-weight", type=float, default=1.0)
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
        max_structure_records=max_records,
        prototype_count=args.prototype_count,
        ligand_weight=args.ligand_weight,
        bootstrap_draws=args.bootstrap_draws,
        seed=args.seed,
        hidden_blocks=args.hidden_blocks,
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
