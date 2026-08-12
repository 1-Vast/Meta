"""Train a sequence-chemistry product observable on crossed Ki quotients.

This replaces the failed T-BASIS input with a simple, falsifiable biological
coordinate: amino-acid composition statistics crossed with ligand chemistry.
It remains an admission Gate only; V1 integration is forbidden unless the Gate
passes against wrong-partner controls.
"""
from __future__ import annotations

import argparse
import gzip
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


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
OUT = CQ_OUT.parent / "seqchem_cq_observable_gate1"
AA = "ACDEFGHIKLMNPQRSTVWY"


def read_jsonl(path: Path) -> list[dict]:
    with path.open("rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def read_jsonl_gz(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def stable_int(value: str) -> int:
    return int(hashlib.sha256(value.encode("utf-8")).hexdigest()[:16], 16)


def protein_descriptor(sequence: str) -> np.ndarray:
    cleaned = "".join(aa for aa in sequence.upper() if aa in AA)
    if not cleaned:
        raise ValueError("empty protein sequence")
    length = len(cleaned)
    counts = np.asarray([cleaned.count(aa) for aa in AA], dtype=np.float64) / length
    groups = {
        "hydrophobic": "AILMFWYV",
        "polar": "STNQCY",
        "positive": "KRH",
        "negative": "DE",
        "gly_pro": "GP",
        "aromatic": "FWY",
    }
    grouped = np.asarray([
        sum(cleaned.count(aa) for aa in members) / length
        for members in groups.values()
    ], dtype=np.float64)
    log_length = np.asarray([np.log1p(length) / 10.0], dtype=np.float64)
    return np.concatenate([counts, grouped, log_length])


def ligand_descriptor(smiles: str, *, fp_bits: int = 96) -> np.ndarray:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid ligand SMILES: {smiles}")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=fp_bits)
    fingerprint = generator.GetFingerprintAsNumPy(molecule).astype(np.float64)
    descriptor = np.asarray([
        Descriptors.MolWt(molecule) / 500.0,
        Crippen.MolLogP(molecule) / 5.0,
        Descriptors.TPSA(molecule) / 200.0,
        Descriptors.NumHDonors(molecule) / 10.0,
        Descriptors.NumHAcceptors(molecule) / 15.0,
        Descriptors.NumRotatableBonds(molecule) / 20.0,
        Descriptors.RingCount(molecule) / 10.0,
        molecule.GetNumHeavyAtoms() / 80.0,
    ], dtype=np.float64)
    return np.concatenate([fingerprint, descriptor])


def product_feature(protein: np.ndarray, ligand: np.ndarray) -> np.ndarray:
    return np.multiply.outer(protein, ligand).ravel()


def donor_maps(cells: list[dict]) -> tuple[dict[str, str], dict[str, str]]:
    target_group = {}
    ligand_scaffold = {}
    for cell in cells:
        target_group[cell["target_id"]] = cell["protein_group_40"]
        ligand_scaffold[cell["ligand_id"]] = cell["scaffold"]
    targets = sorted(target_group)
    ligands = sorted(ligand_scaffold)
    protein_donor = {}
    for index, target in enumerate(targets):
        for offset in range(1, len(targets)):
            candidate = targets[(index + offset) % len(targets)]
            if target_group[candidate] != target_group[target]:
                protein_donor[target] = candidate
                break
        if target not in protein_donor:
            raise ValueError("no different-protein-group donor available")
    ligand_donor = {}
    for index, ligand in enumerate(ligands):
        for offset in range(1, len(ligands)):
            candidate = ligands[(index + offset) % len(ligands)]
            if ligand_scaffold[candidate] != ligand_scaffold[ligand]:
                ligand_donor[ligand] = candidate
                break
        if ligand not in ligand_donor:
            raise ValueError("no different-scaffold donor available")
    return protein_donor, ligand_donor


def materialize_features(corpus: Path) -> tuple[dict[str, dict[str, np.ndarray]], dict]:
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    proteins = {
        row["sequence_sha256"]: protein_descriptor(row["sequence"])
        for row in read_jsonl(corpus / "proteins.jsonl")
    }
    ligands = {
        row["drug_key"]: ligand_descriptor(row["smiles"])
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
        "feature_source": "sequence_composition_x_ligand_morgan_descriptors",
        "feature_dim": int(first.shape[0]),
        "protein_descriptor_dim": int(next(iter(proteins.values())).shape[0]),
        "ligand_descriptor_dim": int(next(iter(ligands.values())).shape[0]),
        "cells": len(cells),
        "proteins": len(proteins),
        "ligands": len(ligands),
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
    }
    return features, metadata


def load_blocks(corpus: Path) -> tuple[list[QuotientBlock], dict]:
    cells = {row["cell_id"]: row for row in read_jsonl_gz(corpus / "cells.jsonl.gz")}
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
    features, metadata = materialize_features(corpus)
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
        corpus: Path = CORPUS, output: Path = OUT, ridge: float = 100.0,
        bootstrap_draws: int = 9999, seed: int = 20260812) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    blocks, metadata = load_blocks(corpus)
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
        "SEQUENCE_CHEM_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "SEQUENCE_CHEM_CQ_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.SequenceChemCQObservableGate1.v1",
        "hypothesis": (
            "A source-trained product of protein sequence composition and ligand "
            "chemistry carries dependency-transferable quotient interaction signal."),
        "literature_mechanism": {
            "deepdta_graphdta": (
                "sequence and ligand graph/SMILES features are valid DTA inputs"),
            "hodge_cycle_space": (
                "quotient scoring removes target and ligand main effects"),
            "adambind_boundary": (
                "few-shot target-as-task remains downstream and unchanged"),
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
    parser.add_argument("--ridge", type=float, default=100.0)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, output=args.output, ridge=args.ridge,
        bootstrap_draws=args.bootstrap_draws, seed=args.seed)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
