"""Generate frozen 288D T-BASIS features for a governed quotient corpus."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import torch

from research.e0_identifiability.run_tbasis_radial import (
    ATOM_CHANNELS,
    aggregate_basis,
    slot_composition,
)
from research.e0_identifiability.run_tdir_pilot import (
    _load_frozen_model,
    _load_protein_rows,
)
from scripts.build_ligand_bank import load_ligand_bank
from scripts.structure_sources.rcsb import sha256_file


def read_jsonl_gz(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def ligand_channels_smiles(smiles: str) -> np.ndarray:
    from rdkit import Chem, RDConfig
    from rdkit.Chem import ChemicalFeatures

    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None or any(atom.GetAtomicNum() == 1 for atom in molecule.GetAtoms()):
        raise ValueError("ligand must be a valid implicit-hydrogen SMILES")
    result = np.zeros((molecule.GetNumAtoms(), len(ATOM_CHANNELS)), dtype=np.float64)
    factory = ChemicalFeatures.BuildFeatureFactory(
        str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef")
    )
    family_to_index = {
        "Hydrophobe": 0,
        "LumpedHydrophobe": 0,
        "Donor": 2,
        "Acceptor": 3,
        "PosIonizable": 4,
        "NegIonizable": 5,
    }
    for feature in factory.GetFeaturesForMol(molecule):
        output_index = family_to_index.get(feature.GetFamily())
        if output_index is not None:
            result[list(feature.GetAtomIds()), output_index] = 1.0
    for atom in molecule.GetAtoms():
        index = atom.GetIdx()
        result[index, 1] = float(atom.GetIsAromatic())
        result[index, 6] = float(atom.GetAtomicNum() in {9, 17, 35, 53})
    result[:, 7] = (result[:, :7].sum(axis=1) == 0).astype(np.float64)
    if np.any(result.sum(axis=1) == 0):
        raise ValueError("ligand channel mapping left an untyped atom")
    return result


def apply_frozen_calibration(raw: np.ndarray, calibration: dict[str, np.ndarray]) -> np.ndarray:
    shape = raw.shape
    calibrated = raw.reshape(-1, shape[-1]) @ calibration["coef"].T
    calibrated += calibration["intercept"]
    normalized = (
        calibrated.reshape(shape).reshape(-1) - calibration["mean"]
    ) / calibration["scale"]
    return normalized[calibration["active"]]


def deterministic_control_map(keys: list[str], incompatible) -> dict[str, str]:
    ordered = sorted(keys)
    result = {}
    for key in ordered:
        start = int(hashlib.sha256(f"CQ-control|{key}".encode()).hexdigest()[:8], 16)
        for offset in range(1, len(ordered) + 1):
            candidate = ordered[(start + offset) % len(ordered)]
            if candidate != key and not incompatible(key, candidate):
                result[key] = candidate
                break
        if key not in result:
            raise ValueError(f"no eligible deterministic control for {key}")
    return result


def _protein_states(model, proteins: dict, device: str, batch_size: int = 64):
    result = {}
    keys = sorted(proteins)
    with torch.inference_mode():
        for start in range(0, len(keys), batch_size):
            batch = keys[start : start + batch_size]
            pooled = torch.from_numpy(
                np.stack([proteins[key]["pooled"] for key in batch])
            ).float().to(device)
            residues = torch.from_numpy(
                np.stack([proteins[key]["residues"] for key in batch])
            ).float().to(device)
            _, states = model.protein(pooled, residues)
            for index, key in enumerate(batch):
                result[key] = {
                    "states": states[index].cpu(),
                    "mask": torch.from_numpy(proteins[key]["mask"].astype(bool)),
                }
    return result


def _ligand_states(model, graphs: dict, device: str, batch_size: int = 32):
    result = {}
    keys = sorted(graphs)
    with torch.inference_mode():
        for start in range(0, len(keys), batch_size):
            batch = keys[start : start + batch_size]
            X = torch.from_numpy(np.stack([graphs[key]["X"] for key in batch])).to(device)
            A = torch.from_numpy(np.stack([graphs[key]["A"] for key in batch])).to(device)
            mask = torch.from_numpy(np.stack([graphs[key]["mask"] for key in batch])).to(device)
            _, states = model.ligand(X.float(), A.float(), mask.float())
            for index, key in enumerate(batch):
                result[key] = {
                    "states": states[index].cpu(),
                    "mask": mask[index].bool().cpu(),
                }
    return result


def _compute_arm(
    cells: list[dict], protein_key, ligand_key, model, protein_states, ligand_states,
    channels, sequences, calibration, device: str, pair_batch_size: int,
) -> np.ndarray:
    features = np.empty((len(cells), int(calibration["active"].sum())), dtype=np.float32)
    with torch.inference_mode():
        for start in range(0, len(cells), pair_batch_size):
            batch = cells[start : start + pair_batch_size]
            p_keys = [protein_key(row) for row in batch]
            l_keys = [ligand_key(row) for row in batch]
            atom_states = torch.stack([ligand_states[key]["states"] for key in l_keys]).to(device)
            atom_mask = torch.stack([ligand_states[key]["mask"] for key in l_keys]).to(device)
            residue_states = torch.stack([protein_states[key]["states"] for key in p_keys]).to(device)
            residue_mask = torch.stack([protein_states[key]["mask"] for key in p_keys]).to(device)
            prediction = model.bridge(atom_states, atom_mask, residue_states, residue_mask)
            distance = torch.softmax(prediction.distance_logits, dim=-1).cpu().numpy()
            for local, (p_key, l_key) in enumerate(zip(p_keys, l_keys)):
                count = int(atom_mask[local].sum().item())
                atom_channels = channels[l_key]
                if len(atom_channels) != count:
                    raise ValueError("ligand chemistry and graph atom order disagree")
                radial = np.einsum(
                    "isb,bk->isk", distance[local, :count], calibration["bin"], optimize=True
                )
                raw = aggregate_basis(
                    atom_channels,
                    slot_composition(sequences[p_key]),
                    radial,
                    residue_mask[local].cpu().numpy(),
                )
                features[start + local] = apply_frozen_calibration(raw, calibration)
    return features


def generate(
    corpus: Path,
    protein_bank: Path,
    ligand_bank: Path,
    checkpoint: Path,
    calibration_path: Path,
    output: Path,
    device: str = "cuda",
    pair_batch_size: int = 64,
) -> dict:
    if not device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("T-BASIS feature generation requires CUDA")
    cells = read_jsonl_gz(corpus / "cells.jsonl.gz")
    proteins_json = [json.loads(line) for line in (corpus / "proteins.jsonl").read_text().splitlines()]
    ligands_json = [json.loads(line) for line in (corpus / "ligands.jsonl").read_text().splitlines()]
    sequences = {row["sequence_sha256"]: row["sequence"] for row in proteins_json}
    smiles = {row["drug_key"]: row["smiles"] for row in ligands_json}
    protein_rows = _load_protein_rows(protein_bank, set(sequences))
    protein_dim = int(next(iter(protein_rows.values()))["residues"].shape[-1])
    model, checkpoint_value = _load_frozen_model(checkpoint, protein_dim, device)
    if checkpoint_value.get("protein_dim", protein_dim) != protein_dim:
        raise ValueError("checkpoint and protein bank dimensions differ")
    started = time.perf_counter()
    protein_states = _protein_states(model, protein_rows, device)
    graphs = load_ligand_bank(ligand_bank)
    ligand_states = _ligand_states(model, graphs, device)
    del graphs
    channels = {key: ligand_channels_smiles(value) for key, value in smiles.items()}
    with np.load(calibration_path, allow_pickle=False) as stored:
        calibration = {
            "coef": stored["calibration_coef"],
            "intercept": stored["calibration_intercept"],
            "mean": stored["train_mean"],
            "scale": stored["train_scale"],
            "active": stored["active"],
            "bin": stored["bin_rbf_expectation"],
        }
    scaffolds = {row["drug_key"]: row["scaffold"] for row in ligands_json}
    groups = {}
    for row in cells:
        groups.setdefault(row["target_id"], row["protein_group_40"])
    foreign_ligand = deterministic_control_map(
        list(smiles), lambda left, right: scaffolds[left] == scaffolds[right]
    )
    wrong_protein = deterministic_control_map(
        list(sequences),
        lambda left, right: groups[left] == groups[right]
        or not 0.5 <= len(sequences[left]) / len(sequences[right]) <= 2.0,
    )
    correct = _compute_arm(
        cells, lambda row: row["target_id"], lambda row: row["ligand_id"], model,
        protein_states, ligand_states, channels, sequences, calibration, device, pair_batch_size,
    )
    foreign = _compute_arm(
        cells, lambda row: row["target_id"], lambda row: foreign_ligand[row["ligand_id"]], model,
        protein_states, ligand_states, channels, sequences, calibration, device, pair_batch_size,
    )
    deranged = _compute_arm(
        cells, lambda row: wrong_protein[row["target_id"]], lambda row: row["ligand_id"], model,
        protein_states, ligand_states, channels, sequences, calibration, device, pair_batch_size,
    )
    elapsed = time.perf_counter() - started
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        cell_id=np.asarray([row["cell_id"] for row in cells]),
        correct=correct,
        foreign_ligand=foreign,
        deranged_protein=deranged,
    )
    manifest = {
        "schema": "MetaSieve.BindingDB.TBasisFeatures.v1",
        "cells": len(cells),
        "dimensions": int(correct.shape[1]),
        "arms": ["correct", "foreign_ligand", "deranged_protein"],
        "seconds": elapsed,
        "cells_per_second": len(cells) / elapsed,
        "device": device,
        "checkpoint_sha256": sha256_file(checkpoint),
        "calibration_sha256": sha256_file(calibration_path),
        "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
        "protein_manifest_sha256": sha256_file(protein_bank / "manifest.json"),
        "ligand_manifest_sha256": sha256_file(ligand_bank / "manifest.json"),
        "features_sha256": sha256_file(output),
        "affinity_values_used_by_feature_generator": 0,
    }
    output.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--corpus", type=Path, required=True)
    parser.add_argument("--protein-bank", type=Path, required=True)
    parser.add_argument("--ligand-bank", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--pair-batch-size", type=int, default=64)
    args = parser.parse_args()
    print(json.dumps(generate(
        corpus=args.corpus,
        protein_bank=args.protein_bank,
        ligand_bank=args.ligand_bank,
        checkpoint=args.checkpoint,
        calibration_path=args.calibration,
        output=args.output,
        device=args.device,
        pair_batch_size=args.pair_batch_size,
    ), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
