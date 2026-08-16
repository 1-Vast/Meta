"""Descriptor-based G2 Gate for BindingDB rectangle quotients.

This replaces the failed exact transformation-ID X2 arm with a deployment-
computable quotient model. For each complete 2x2 rectangle it predicts

    R = y(P_a,L_b) - y(P_a,L_a) - y(P_b,L_b) + y(P_b,L_a)

from the crossed descriptor

    (protein(P_a) - protein(P_b)) x (ligand(L_b) - ligand(L_a)).

The arm is an admission Gate only. It does not train an end-to-end DTA model and
does not claim latent non-additivity without replicate/noise correction.
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

from research.crossed_interaction.audit_bindingdb_rectangle_interaction import (
    CORPUS,
    OUT as RECTANGLE_OUT,
    build_rectangles,
)
from research.crossed_interaction.bindingdb_cq_r0 import sha256_file
from research.crossed_interaction.train_cq_observable import bootstrap_contrast
from research.crossed_interaction.train_plm_slot_cq_observable import (
    PROTEIN_BANK,
    load_protein_bank,
    protein_slot_descriptor,
)
from research.crossed_interaction.train_seqchem_cq_observable import (
    ligand_descriptor,
    protein_descriptor,
    read_jsonl,
)


OUT = RECTANGLE_OUT.parent / "bindingdb_rectangle_descriptor_g2"


def protein_plm_slot_descriptors(
        corpus: Path, protein_bank: Path, *,
        slot_segments: int, hidden_blocks: int) -> tuple[dict[str, np.ndarray], dict]:
    proteins_json = read_jsonl(corpus / "proteins.jsonl")
    required = {row["sequence_sha256"] for row in proteins_json}
    protein_rows, manifest = load_protein_bank(protein_bank, required)
    descriptors = {
        key: protein_slot_descriptor(
            row["residues"], row["mask"],
            slot_segments=slot_segments,
            hidden_blocks=hidden_blocks).astype(np.float64)
        for key, row in protein_rows.items()
    }
    metadata = {
        "protein_descriptor": "frozen_esm2_slot_region_means",
        "protein_descriptor_dim": int(next(iter(descriptors.values())).shape[0]),
        "protein_bank_manifest_sha256": sha256_file(protein_bank / "manifest.json"),
        "protein_bank_model_id": manifest.get("model_id", ""),
        "protein_bank_model_revision": manifest.get("model_revision", ""),
        "protein_bank_slot_policy": manifest.get("slot_policy", ""),
        "slot_segments": slot_segments,
        "hidden_blocks": hidden_blocks,
    }
    return descriptors, metadata


def crossed_feature(protein_delta: np.ndarray, ligand_delta: np.ndarray) -> np.ndarray:
    return np.multiply.outer(protein_delta, ligand_delta).ravel()


def load_descriptors(
        corpus: Path, *, protein_descriptor_mode: str,
        protein_bank: Path, slot_segments: int,
        hidden_blocks: int) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict]:
    if protein_descriptor_mode == "composition":
        proteins = {
            row["sequence_sha256"]: protein_descriptor(row["sequence"]).astype(np.float64)
            for row in read_jsonl(corpus / "proteins.jsonl")
        }
        protein_metadata = {
            "protein_descriptor": "amino_acid_composition_grouped_length",
            "protein_descriptor_dim": int(next(iter(proteins.values())).shape[0]),
        }
    elif protein_descriptor_mode == "plm_slots":
        proteins, protein_metadata = protein_plm_slot_descriptors(
            corpus, protein_bank,
            slot_segments=slot_segments, hidden_blocks=hidden_blocks)
    else:
        raise ValueError(f"unknown protein descriptor mode: {protein_descriptor_mode}")
    ligands = {
        row["drug_key"]: ligand_descriptor(row["smiles"]).astype(np.float64)
        for row in read_jsonl(corpus / "ligands.jsonl")
    }
    first_protein = next(iter(proteins.values()))
    first_ligand = next(iter(ligands.values()))
    metadata = {
        "ligand_descriptor": "morgan_radius2_96_bits_plus_physchem",
        "ligand_descriptor_dim": int(first_ligand.shape[0]),
        "crossed_descriptor_dim": int(first_protein.shape[0] * first_ligand.shape[0]),
        "proteins": len(proteins),
        "ligands": len(ligands),
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
    }
    metadata.update(protein_metadata)
    return proteins, ligands, metadata


def materialize_examples(
        corpus: Path, *, protein_descriptor_mode: str,
        protein_bank: Path, slot_segments: int,
        hidden_blocks: int) -> tuple[list[dict], dict]:
    rows, rectangle_metadata = build_rectangles(corpus)
    proteins, ligands, descriptor_metadata = load_descriptors(
        corpus,
        protein_descriptor_mode=protein_descriptor_mode,
        protein_bank=protein_bank,
        slot_segments=slot_segments,
        hidden_blocks=hidden_blocks)
    examples = []
    for row in rows:
        protein_delta = proteins[row["target_a"]] - proteins[row["target_b"]]
        ligand_delta = ligands[row["ligand_b"]] - ligands[row["ligand_a"]]
        examples.append({
            **row,
            "target_pair": "|".join([row["target_a"], row["target_b"]]),
            "feature_correct": crossed_feature(protein_delta, ligand_delta),
            "feature_ligand_only": ligand_delta,
            "protein_delta": protein_delta,
            "ligand_delta": ligand_delta,
            "y": float(row["rectangle"]),
        })
    metadata = {
        "rectangle_corpus": rectangle_metadata,
        "descriptors": descriptor_metadata,
        "examples": len(examples),
    }
    return examples, metadata


def fit_ridge(x: np.ndarray, y: np.ndarray, ridge: float) -> dict:
    if ridge <= 0:
        raise ValueError("ridge must be strictly positive")
    scale = x.std(axis=0)
    scale[scale < 1e-6] = 1.0
    x_scaled = x / scale
    identity = np.eye(x_scaled.shape[1], dtype=np.float64)
    weights = np.linalg.solve(x_scaled.T @ x_scaled + ridge * identity, x_scaled.T @ y)
    train_prediction = x_scaled @ weights
    return {
        "ridge": ridge,
        "scale": scale,
        "weights": weights,
        "train_mse": float(np.square(y - train_prediction).mean()),
    }


def predict(model: dict, x: np.ndarray) -> np.ndarray:
    return (x / model["scale"]) @ model["weights"]


def stack(examples: list[dict], key: str) -> tuple[np.ndarray, np.ndarray]:
    return (
        np.stack([row[key] for row in examples]).astype(np.float64),
        np.asarray([row["y"] for row in examples], dtype=np.float64),
    )


def shuffled_protein_features(examples: list[dict], *, seed: int) -> tuple[np.ndarray, dict]:
    rng = np.random.default_rng(seed)
    protein_deltas = [row["protein_delta"] for row in examples]
    ligand_deltas = [row["ligand_delta"] for row in examples]
    order = np.arange(len(examples))
    if len(order) > 1:
        for _ in range(100):
            rng.shuffle(order)
            if np.all(order != np.arange(len(examples))):
                break
        if np.any(order == np.arange(len(examples))):
            order = np.roll(np.arange(len(examples)), 1)
    features = np.stack([
        crossed_feature(protein_deltas[int(source)], ligand_delta)
        for source, ligand_delta in zip(order, ligand_deltas)
    ]).astype(np.float64)
    return features, {
        "seed": seed,
        "fixed_points": int(np.sum(order == np.arange(len(examples)))),
    }


def wrong_target_features(examples: list[dict]) -> tuple[np.ndarray, dict]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(examples):
        grouped[row["dependency_component"]].append(index)
    donor = {}
    components = sorted(grouped)
    for component_index, component in enumerate(components):
        candidates = [
            index
            for offset in range(1, len(components) + 1)
            for index in grouped[components[(component_index + offset) % len(components)]]
            if components[(component_index + offset) % len(components)] != component
        ]
        if not candidates:
            candidates = [index for index in range(len(examples)) if index not in grouped[component]]
        for local_index, index in enumerate(grouped[component]):
            donor[index] = candidates[local_index % len(candidates)]
    features = np.stack([
        crossed_feature(examples[donor[index]]["protein_delta"], row["ligand_delta"])
        for index, row in enumerate(examples)
    ]).astype(np.float64)
    return features, {
        "dependency_components": len(components),
        "self_donors": int(sum(index == source for index, source in donor.items())),
    }


def rows_for_scoring(examples: list[dict], predictions: dict[str, np.ndarray]) -> list[dict]:
    rows = []
    y = np.asarray([row["y"] for row in examples], dtype=np.float64)
    for arm, prediction in predictions.items():
        for example, squared in zip(examples, np.square(y - prediction)):
            rows.append({
                "panel_id": example["panel_id"],
                "split": example["split"],
                "dependency_component": example["dependency_component"],
                "arm": arm,
                "retained_rank": 1,
                "rank_normalized_mse": float(squared),
                "row_mse": float(squared),
            })
    return rows


def summarize(rows: list[dict]) -> dict:
    result = {}
    for arm in sorted({row["arm"] for row in rows}):
        selected = [row for row in rows if row["arm"] == arm]
        result[arm] = {
            "rows": len(selected),
            "mse": float(np.mean([row["row_mse"] for row in selected])),
            "component_macro_mse": float(np.mean([
                np.mean([
                    row["row_mse"] for row in selected
                    if row["dependency_component"] == component
                ])
                for component in sorted({row["dependency_component"] for row in selected})
            ])),
        }
    return result


def split_summary(examples: list[dict]) -> dict:
    return {
        "rectangles": len(examples),
        "panels": len({row["panel_id"] for row in examples}),
        "target_pairs": len({row["target_pair"] for row in examples}),
        "ligand_pairs": len({
            tuple(sorted([row["ligand_a"], row["ligand_b"]]))
            for row in examples
        }),
        "dependency_components": len({row["dependency_component"] for row in examples}),
    }


def run(
        corpus: Path = CORPUS, protein_bank: Path = PROTEIN_BANK,
        output: Path = OUT, protein_descriptor_mode: str = "composition",
        slot_segments: int = 4, hidden_blocks: int = 8,
        ridge: float = 100.0, bootstrap_draws: int = 9999,
        seed: int = 20260812) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    examples, metadata = materialize_examples(
        corpus,
        protein_descriptor_mode=protein_descriptor_mode,
        protein_bank=protein_bank,
        slot_segments=slot_segments,
        hidden_blocks=hidden_blocks)
    train = [row for row in examples if row["split"] == "train"]
    development = [row for row in examples if row["split"] == "development"]
    if not train or not development:
        raise ValueError("train and development examples are required")

    x_train, y_train = stack(train, "feature_correct")
    model = fit_ridge(x_train, y_train, ridge)
    x_ligand_train, _ = stack(train, "feature_ligand_only")
    ligand_model = fit_ridge(x_ligand_train, y_train, ridge)

    x_dev, y_dev = stack(development, "feature_correct")
    x_ligand_dev, _ = stack(development, "feature_ligand_only")
    wrong_x, wrong_meta = wrong_target_features(development)
    shuffled_x, shuffled_meta = shuffled_protein_features(development, seed=seed + 17)
    predictions = {
        "zero": np.zeros_like(y_dev),
        "ligand_only": predict(ligand_model, x_ligand_dev),
        "correct": predict(model, x_dev),
        "wrong_protein": predict(model, wrong_x),
        "shuffled_protein": predict(model, shuffled_x),
    }
    development_rows = rows_for_scoring(development, predictions)
    development_summary = summarize(development_rows)
    contrasts = [
        bootstrap_contrast(development_rows, "correct", control, draws=bootstrap_draws, seed=seed + i)
        for i, control in enumerate(("zero", "ligand_only", "wrong_protein", "shuffled_protein"))
    ]
    contrast_by_control = {row["control"]: row for row in contrasts}
    gates = {
        "train_rectangles_ge_10000": len(train) >= 10000,
        "development_rectangles_ge_1000": len(development) >= 1000,
        "development_components_ge_5": split_summary(development)["dependency_components"] >= 5,
        "correct_beats_zero": contrast_by_control["zero"]["pass"],
        "correct_beats_ligand_only": contrast_by_control["ligand_only"]["pass"],
        "correct_beats_wrong_protein": contrast_by_control["wrong_protein"]["pass"],
        "correct_beats_shuffled_protein": contrast_by_control["shuffled_protein"]["pass"],
    }
    verdict = (
        "BINDINGDB_RECTANGLE_DESCRIPTOR_G2_PASS"
        if all(gates.values())
        else "BINDINGDB_RECTANGLE_DESCRIPTOR_G2_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.BindingDBRectangleDescriptorG2.v1",
        "hypothesis": (
            "A deployment-computable crossed descriptor, "
            "(protein_a-protein_b) x (ligand_b-ligand_a), transfers a cold "
            "protein-conditioned quotient signal across BindingDB components."),
        "literature_mechanism": {
            "matched_molecular_pair_analysis": (
                "chemical transformations provide interpretable SAR perturbation units"),
            "activity_cliffs": (
                "near-neighbor potency changes are a stringent SAR generalization test"),
            "double_mutant_cycle": (
                "2x2 cycles algebraically remove main effects before interaction testing"),
            "unipert_g2cp_boundary": (
                "retains cross-modal relation supervision only; no full UniPert model is imported"),
        },
        "metadata": metadata,
        "config": {
            "ridge": ridge,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
            "protein_descriptor_mode": protein_descriptor_mode,
            "slot_segments": slot_segments if protein_descriptor_mode == "plm_slots" else None,
            "hidden_blocks": hidden_blocks if protein_descriptor_mode == "plm_slots" else None,
            "labels_used_for_training": True,
            "support_query_isolation": "train_split_fit_development_split_score",
            "development_training_authorized": False,
            "latent_nonadditivity_claim_authorized": False,
        },
        "split_summary": {
            "train": split_summary(train),
            "development": split_summary(development),
        },
        "models": {
            "correct_train_mse": model["train_mse"],
            "ligand_only_train_mse": ligand_model["train_mse"],
            "wrong_protein": wrong_meta,
            "shuffled_protein": shuffled_meta,
        },
        "development_summary": development_summary,
        "contrasts": contrasts,
        "gates": gates,
        "g3a_authorized": verdict.endswith("PASS"),
        "v1_integration_authorized": False,
        "biological_claim_authorized": False,
        "TERMINAL_VERDICT": verdict,
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--protein-bank", type=Path, default=PROTEIN_BANK)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--protein-descriptor-mode",
        choices=["composition", "plm_slots"],
        default="composition")
    parser.add_argument("--slot-segments", type=int, default=4)
    parser.add_argument("--hidden-blocks", type=int, default=8)
    parser.add_argument("--ridge", type=float, default=100.0)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    args = parser.parse_args()
    result = run(
        corpus=args.corpus,
        protein_bank=args.protein_bank,
        output=args.output,
        protein_descriptor_mode=args.protein_descriptor_mode,
        slot_segments=args.slot_segments,
        hidden_blocks=args.hidden_blocks,
        ridge=args.ridge,
        bootstrap_draws=args.bootstrap_draws,
        seed=args.seed)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
