"""Train a local k-mer protein motif observable on crossed Ki quotients.

This is a source-only admission Gate for a biological pair coordinate. It
keeps the quotient target unchanged and only replaces the failed global
sequence-composition axis with a local protein k-mer motif hash crossed with
ligand chemistry.
"""
from __future__ import annotations

import argparse
import hashlib
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
from research.crossed_interaction.train_seqchem_cq_observable import (
    AA,
    donor_maps,
    ligand_descriptor,
    product_feature,
    read_jsonl,
    read_jsonl_gz,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
OUT = CQ_OUT.parent / "kmer_cq_observable_gate1"


def stable_int(value: str) -> int:
    return int(hashlib.sha256(value.encode("utf-8")).hexdigest()[:16], 16)


def protein_kmer_descriptor(sequence: str, *, k: int = 3, bins: int = 32) -> np.ndarray:
    if k <= 0:
        raise ValueError("k-mer size must be positive")
    if bins <= 0:
        raise ValueError("k-mer bins must be positive")
    cleaned = "".join(aa for aa in sequence.upper() if aa in AA)
    if not cleaned:
        raise ValueError("empty protein sequence")

    descriptor = np.zeros(bins, dtype=np.float64)
    if len(cleaned) < k:
        descriptor[stable_int(cleaned) % bins] = 1.0
        return descriptor

    total = len(cleaned) - k + 1
    for index in range(total):
        descriptor[stable_int(cleaned[index:index + k]) % bins] += 1.0
    return descriptor / total


def materialize_features(
        corpus: Path, *, kmer_size: int, kmer_bins: int,
        ligand_fp_bits: int) -> tuple[dict[str, dict[str, np.ndarray]], dict]:
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    proteins = {
        row["sequence_sha256"]: protein_kmer_descriptor(
            row["sequence"], k=kmer_size, bins=kmer_bins)
        for row in read_jsonl(corpus / "proteins.jsonl")
    }
    ligands = {
        row["drug_key"]: ligand_descriptor(row["smiles"], fp_bits=ligand_fp_bits)
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
        "feature_source": "protein_kmer_motif_hash_x_ligand_morgan_descriptors",
        "feature_dim": int(first.shape[0]),
        "protein_descriptor_dim": int(next(iter(proteins.values())).shape[0]),
        "ligand_descriptor_dim": int(next(iter(ligands.values())).shape[0]),
        "kmer_size": kmer_size,
        "kmer_bins": kmer_bins,
        "ligand_fp_bits": ligand_fp_bits,
        "cells": len(cells),
        "proteins": len(proteins),
        "ligands": len(ligands),
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
    }
    return features, metadata


def load_blocks(
        corpus: Path, *, kmer_size: int, kmer_bins: int,
        ligand_fp_bits: int) -> tuple[list[QuotientBlock], dict]:
    cells = {row["cell_id"]: row for row in read_jsonl_gz(corpus / "cells.jsonl.gz")}
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
    features, metadata = materialize_features(
        corpus, kmer_size=kmer_size, kmer_bins=kmer_bins,
        ligand_fp_bits=ligand_fp_bits)
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
        corpus: Path = CORPUS, output: Path = OUT, ridge: float = 10000.0,
        bootstrap_draws: int = 9999, seed: int = 20260812,
        kmer_size: int = 3, kmer_bins: int = 64,
        ligand_fp_bits: int = 64) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    blocks, metadata = load_blocks(
        corpus, kmer_size=kmer_size, kmer_bins=kmer_bins,
        ligand_fp_bits=ligand_fp_bits)
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
        "KMER_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "KMER_CQ_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.KmerCQObservableGate1.v1",
        "hypothesis": (
            "Local protein k-mer motif frequencies crossed with ligand chemistry "
            "carry dependency-transferable quotient interaction signal."),
        "literature_mechanism": {
            "deepdta": (
                "protein sequence and ligand sequence encoders can support DTA"),
            "widedta": (
                "protein domains or motifs and ligand substructure words are "
                "useful interaction coordinates"),
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
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ridge", type=float, default=10000.0)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--kmer-size", type=int, default=3)
    parser.add_argument("--kmer-bins", type=int, default=64)
    parser.add_argument("--ligand-fp-bits", type=int, default=64)
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, output=args.output, ridge=args.ridge,
        bootstrap_draws=args.bootstrap_draws, seed=args.seed,
        kmer_size=args.kmer_size, kmer_bins=args.kmer_bins,
        ligand_fp_bits=args.ligand_fp_bits)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
