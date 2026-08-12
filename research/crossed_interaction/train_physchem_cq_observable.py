"""Train a physicochemical autocorrelation observable on crossed Ki quotients.

This source-only admission Gate tests a more typed biological coordinate than
sequence text products: short-lag protein residue-property autocorrelations
crossed with ligand pharmacophore and E-state-like descriptors.
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
from research.crossed_interaction.train_seqchem_cq_observable import (
    AA,
    donor_maps,
    product_feature,
    read_jsonl,
    read_jsonl_gz,
)


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
OUT = CQ_OUT.parent / "physchem_cq_observable_gate1"

AA_PROPERTIES = {
    "hydropathy": {
        "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
        "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
        "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
        "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
    },
    "volume": {
        "A": 88.6, "C": 108.5, "D": 111.1, "E": 138.4, "F": 189.9,
        "G": 60.1, "H": 153.2, "I": 166.7, "K": 168.6, "L": 166.7,
        "M": 162.9, "N": 114.1, "P": 112.7, "Q": 143.8, "R": 173.4,
        "S": 89.0, "T": 116.1, "V": 140.0, "W": 227.8, "Y": 193.6,
    },
    "polarity": {
        "A": 8.1, "C": 5.5, "D": 13.0, "E": 12.3, "F": 5.2,
        "G": 9.0, "H": 10.4, "I": 5.2, "K": 11.3, "L": 4.9,
        "M": 5.7, "N": 11.6, "P": 8.0, "Q": 10.5, "R": 10.5,
        "S": 9.2, "T": 8.6, "V": 5.9, "W": 5.4, "Y": 6.2,
    },
    "charge": {
        "A": 0.0, "C": 0.0, "D": -1.0, "E": -1.0, "F": 0.0,
        "G": 0.0, "H": 0.1, "I": 0.0, "K": 1.0, "L": 0.0,
        "M": 0.0, "N": 0.0, "P": 0.0, "Q": 0.0, "R": 1.0,
        "S": 0.0, "T": 0.0, "V": 0.0, "W": 0.0, "Y": 0.0,
    },
    "aromatic": {
        "A": 0.0, "C": 0.0, "D": 0.0, "E": 0.0, "F": 1.0,
        "G": 0.0, "H": 1.0, "I": 0.0, "K": 0.0, "L": 0.0,
        "M": 0.0, "N": 0.0, "P": 0.0, "Q": 0.0, "R": 0.0,
        "S": 0.0, "T": 0.0, "V": 0.0, "W": 1.0, "Y": 1.0,
    },
}


def _normalized_property_values(property_map: dict[str, float]) -> dict[str, float]:
    values = np.asarray([property_map[aa] for aa in AA], dtype=np.float64)
    mean = values.mean()
    std = values.std()
    if std < 1e-12:
        std = 1.0
    return {aa: float((property_map[aa] - mean) / std) for aa in AA}


NORMALIZED_PROPERTIES = {
    name: _normalized_property_values(values)
    for name, values in AA_PROPERTIES.items()
}


def protein_physchem_descriptor(sequence: str, *, max_lag: int = 8) -> np.ndarray:
    if max_lag <= 0:
        raise ValueError("max_lag must be positive")
    cleaned = "".join(aa for aa in sequence.upper() if aa in AA)
    if not cleaned:
        raise ValueError("empty protein sequence")

    features = []
    for property_values in NORMALIZED_PROPERTIES.values():
        values = np.asarray([property_values[aa] for aa in cleaned], dtype=np.float64)
        centered = values - values.mean()
        variance = float(np.mean(np.square(centered)))
        for lag in range(1, max_lag + 1):
            if lag >= len(values) or variance < 1e-12:
                features.append(0.0)
            else:
                numerator = float(np.mean(centered[:-lag] * centered[lag:]))
                features.append(numerator / variance)
    return np.asarray(features, dtype=np.float64)


def ligand_pharmacophore_descriptor(smiles: str) -> np.ndarray:
    from rdkit import Chem
    from rdkit.Chem import Crippen, Descriptors
    from rdkit.Chem.EState import Fingerprinter

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid ligand SMILES: {smiles}")
    formal_charge = sum(atom.GetFormalCharge() for atom in molecule.GetAtoms())
    hetero_atoms = sum(1 for atom in molecule.GetAtoms() if atom.GetAtomicNum() not in (1, 6))
    estate_counts, estate_sums = Fingerprinter.FingerprintMol(molecule)
    estate_bins = np.zeros(32, dtype=np.float64)
    for index, value in enumerate(estate_counts):
        estate_bins[index % 16] += float(value)
    for index, value in enumerate(estate_sums):
        estate_bins[16 + index % 16] += float(value)
    heavy_atoms = max(float(molecule.GetNumHeavyAtoms()), 1.0)
    estate_bins[:16] /= heavy_atoms
    estate_bins[16:] /= max(heavy_atoms * 10.0, 1.0)
    scalars = np.asarray([
        Descriptors.MolWt(molecule) / 500.0,
        Crippen.MolLogP(molecule) / 5.0,
        Descriptors.TPSA(molecule) / 200.0,
        Descriptors.NumHDonors(molecule) / 10.0,
        Descriptors.NumHAcceptors(molecule) / 15.0,
        Descriptors.NumRotatableBonds(molecule) / 20.0,
        Descriptors.RingCount(molecule) / 10.0,
        hetero_atoms / 30.0 + formal_charge / 10.0,
    ], dtype=np.float64)
    return np.concatenate([estate_bins, scalars])


def materialize_features(
        corpus: Path, *, max_lag: int) -> tuple[dict[str, dict[str, np.ndarray]], dict]:
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    proteins = {
        row["sequence_sha256"]: protein_physchem_descriptor(row["sequence"], max_lag=max_lag)
        for row in read_jsonl(corpus / "proteins.jsonl")
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
        "feature_source": "protein_physchem_autocorrelation_x_ligand_pharmacophore",
        "feature_dim": int(first.shape[0]),
        "protein_descriptor_dim": int(next(iter(proteins.values())).shape[0]),
        "ligand_descriptor_dim": int(next(iter(ligands.values())).shape[0]),
        "protein_property_channels": list(NORMALIZED_PROPERTIES),
        "max_lag": max_lag,
        "cells": len(cells),
        "proteins": len(proteins),
        "ligands": len(ligands),
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
    }
    return features, metadata


def load_blocks(corpus: Path, *, max_lag: int) -> tuple[list[QuotientBlock], dict]:
    cells = {row["cell_id"]: row for row in read_jsonl_gz(corpus / "cells.jsonl.gz")}
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
    features, metadata = materialize_features(corpus, max_lag=max_lag)
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
        max_lag: int = 8) -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    blocks, metadata = load_blocks(corpus, max_lag=max_lag)
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
        "PHYSCHEM_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "PHYSCHEM_CQ_GATE1_FAIL_CLOSED"
    )
    result = {
        "schema": "MetaSieve.PhyschemCQObservableGate1.v1",
        "hypothesis": (
            "Short-lag protein physicochemical autocorrelations crossed with "
            "ligand pharmacophore descriptors carry dependency-transferable "
            "quotient interaction signal."),
        "literature_mechanism": {
            "classical_qsar_dti": (
                "sequence autocorrelation and ligand pharmacophore/E-state "
                "descriptors encode typed biochemical neighborhoods"),
            "hodge_cycle_space": (
                "quotient scoring removes target and ligand main effects"),
            "adambind_boundary": (
                "few-shot target-as-task V1 remains downstream and unchanged"),
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
    parser.add_argument("--max-lag", type=int, default=8)
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, output=args.output, ridge=args.ridge,
        bootstrap_draws=args.bootstrap_draws, seed=args.seed,
        max_lag=args.max_lag)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
