"""Generate label-blind, source-donor wrong-protein controls for v1 meta-validation."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import torch

from research.crossed_interaction.generate_tbasis_features import (
    _compute_arm,
    _ligand_states,
    _protein_states,
    ligand_channels_smiles,
)
from research.e0_identifiability.run_tdir_pilot import _load_frozen_model, _load_protein_rows
from research.meta_fewshot.seal_v1_development import read_gzip_jsonl, sha256, stable_hash
from scripts.build_ligand_bank import load_ligand_bank

ROOT = Path(__file__).resolve().parents[2]
MAIN = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
DEV = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_v1_development"
PROTEIN_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank"
LIGAND_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank"
CHECKPOINT = ROOT / "report/mechanism_refactor/p1b_pilot20k_seed17_v1/best.pt"
CALIBRATION = ROOT / "research/e0_identifiability/artifacts/tbasis_r0_v1/basis_values.npz"
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def composition(sequence: str) -> np.ndarray:
    count = np.asarray([sequence.count(code) for code in AMINO_ACIDS], dtype=np.float64)
    return count / max(len(sequence), 1)


def matched_donor_map(source: list[dict], validation: list[dict], sequences: dict[str, str]) -> dict[str, str]:
    source_target = {row["target_id"]: row["protein_group_40"] for row in source}
    val_target = {row["target_id"]: row["protein_group_40"] for row in validation}
    source_profiles = {
        target: (len(sequences[target]), composition(sequences[target])) for target in source_target
    }
    result = {}
    for target in sorted(val_target):
        length, amino = len(sequences[target]), composition(sequences[target])
        candidates = []
        for donor, donor_group in source_target.items():
            if donor_group == val_target[target]:
                continue
            donor_length, donor_amino = source_profiles[donor]
            score = abs(math.log(donor_length / length)) + 0.5 * float(np.abs(donor_amino - amino).sum())
            candidates.append((score, abs(donor_length - length), donor))
        if not candidates:
            raise ValueError(f"no source donor for {target}")
        result[target] = min(candidates)[2]
    return result


def generate(device: str = "cuda", pair_batch_size: int = 64) -> dict:
    if not device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("v1 wrong-feature generation requires CUDA")
    source = list(read_gzip_jsonl(DEV / "source_cells.jsonl.gz"))
    validation = list(read_gzip_jsonl(DEV / "metaval_cells_without_labels.jsonl.gz"))
    if any("pK" in row for row in validation):
        raise ValueError("meta-validation feature generation received query labels")
    proteins_json = [json.loads(line) for line in (MAIN / "proteins.jsonl").read_text().splitlines()]
    ligands_json = [json.loads(line) for line in (MAIN / "ligands.jsonl").read_text().splitlines()]
    sequences_all = {row["sequence_sha256"]: row["sequence"] for row in proteins_json}
    selected_targets = {row["target_id"] for row in source + validation}
    sequences = {key: sequences_all[key] for key in selected_targets}
    donor_map = matched_donor_map(source, validation, sequences)
    donor_targets = set(donor_map.values())

    protein_rows = _load_protein_rows(PROTEIN_BANK, selected_targets | donor_targets)
    protein_dim = int(next(iter(protein_rows.values()))["residues"].shape[-1])
    model, _ = _load_frozen_model(CHECKPOINT, protein_dim, device)
    protein_states = _protein_states(model, protein_rows, device)
    graphs = load_ligand_bank(LIGAND_BANK)
    required_ligands = {row["ligand_id"] for row in validation}
    graphs = {key: graphs[key] for key in required_ligands}
    ligand_states = _ligand_states(model, graphs, device)
    smiles = {row["drug_key"]: row["smiles"] for row in ligands_json if row["drug_key"] in required_ligands}
    channels = {key: ligand_channels_smiles(value) for key, value in smiles.items()}
    with np.load(CALIBRATION, allow_pickle=False) as stored:
        calibration = {
            "coef": stored["calibration_coef"], "intercept": stored["calibration_intercept"],
            "mean": stored["train_mean"], "scale": stored["train_scale"],
            "active": stored["active"], "bin": stored["bin_rbf_expectation"],
        }
    wrong = _compute_arm(
        validation, lambda row: donor_map[row["target_id"]], lambda row: row["ligand_id"],
        model, protein_states, ligand_states, channels, sequences, calibration,
        device, pair_batch_size,
    )
    output = DEV / "metaval_wrong_features.npz"
    np.savez_compressed(output,
                        cell_id=np.asarray([row["cell_id"] for row in validation]), wrong=wrong)
    map_path = DEV / "metaval_wrong_protein_map.json"
    map_path.write_text(json.dumps(donor_map, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    source_targets = {row["target_id"] for row in source}
    validation_targets = {row["target_id"] for row in validation}
    lengths = [abs(len(sequences[target]) - len(sequences[donor]))
               for target, donor in donor_map.items()]
    manifest = {
        "schema": "MetaSieve.V1WrongProteinFeatures.v1",
        "cells": len(validation), "targets": len(donor_map),
        "query_affinity_values_used": 0,
        "donor_pool": "meta_train_only",
        "all_donors_in_source": set(donor_map.values()) <= source_targets,
        "donors_in_metaval": len(set(donor_map.values()) & validation_targets),
        "different_cdhit40_group": True,
        "matching": "minimum abs(log-length-ratio) + 0.5*amino-composition-L1",
        "median_absolute_length_difference": float(np.median(lengths)),
        "max_absolute_length_difference": int(max(lengths)),
        "donor_map_hash": stable_hash(f"{key}:{value}" for key, value in donor_map.items()),
        "files": {"metaval_wrong_features.npz": sha256(output),
                  "metaval_wrong_protein_map.json": sha256(map_path)},
        "checkpoint_sha256": sha256(CHECKPOINT), "calibration_sha256": sha256(CALIBRATION),
    }
    manifest_path = DEV / "metaval_wrong_features.manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    print(json.dumps(generate(), indent=2, sort_keys=True))
