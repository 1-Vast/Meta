"""Stage R2: frozen ligand/protein representation discriminator, double-cold.

Criteria are fixed in
`report/meta_fewshot/stageR2_representation_discriminator_20260815/PREREGISTRATION.md`.
No training and no model. Neighbour counts are fixed a priori so that no arm can
win by having its bandwidth tuned.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.evaluate_qpsmp import concordance_index, spearman
from scripts.qpsmp_data import QPSMPData
from scripts.stageR0_retrieval_falsification import (
    component_bootstrap, component_target_mean, tanimoto_rows,
)
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
)

LIGAND_NEIGHBORS = 10
PROTEIN_NEIGHBORS = 16
ROBUSTNESS_NEIGHBORS = (1, 25)
LOW_TIER = 0.4


def unit(matrix: np.ndarray) -> np.ndarray:
    return matrix / np.maximum(np.linalg.norm(matrix, axis=-1, keepdims=True), 1e-9)


def rdkit_descriptors(smiles: dict[str, str | None], keys: list[str]) -> np.ndarray:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    rows = []
    for key in keys:
        molecule = Chem.MolFromSmiles(smiles.get(key) or "")
        if molecule is None:
            rows.append([0.0] * 8)
            continue
        rows.append([
            Descriptors.MolWt(molecule), Crippen.MolLogP(molecule),
            rdMolDescriptors.CalcTPSA(molecule),
            rdMolDescriptors.CalcNumHBD(molecule),
            rdMolDescriptors.CalcNumHBA(molecule),
            rdMolDescriptors.CalcNumRotatableBonds(molecule),
            rdMolDescriptors.CalcNumRings(molecule),
            rdMolDescriptors.CalcFractionCSP3(molecule),
        ])
    return np.asarray(rows, dtype=np.float32)


def chemberta_embeddings(smiles: dict[str, str | None], keys: list[str],
                         device: str) -> np.ndarray:
    import torch
    from transformers import AutoModel, AutoTokenizer
    name = "DeepChem/ChemBERTa-77M-MLM"
    tokenizer = AutoTokenizer.from_pretrained(name)
    model = AutoModel.from_pretrained(name).to(device).eval()
    out = []
    with torch.no_grad():
        for start in range(0, len(keys), 128):
            batch = [smiles.get(key) or "C" for key in keys[start:start + 128]]
            encoded = tokenizer(batch, padding=True, truncation=True,
                                max_length=256, return_tensors="pt").to(device)
            hidden = model(**encoded).last_hidden_state
            mask = encoded["attention_mask"].unsqueeze(-1).to(hidden.dtype)
            out.append(((hidden * mask).sum(1)
                        / mask.sum(1).clamp_min(1.0)).float().cpu().numpy())
    del model
    return np.concatenate(out, 0)


def kmer_matrix(sequences: list[str], k: int = 3) -> np.ndarray:
    vocabulary: dict[str, int] = {}
    counts = []
    for sequence in sequences:
        local: dict[int, float] = defaultdict(float)
        for index in range(len(sequence) - k + 1):
            token = sequence[index:index + k]
            if token not in vocabulary:
                vocabulary[token] = len(vocabulary)
            local[vocabulary[token]] += 1.0
        counts.append(local)
    matrix = np.zeros((len(sequences), len(vocabulary)), dtype=np.float32)
    for row, local in enumerate(counts):
        for column, value in local.items():
            matrix[row, column] = value
    return matrix


def neighbor_mean(similarity: np.ndarray, values: np.ndarray, count: int
                  ) -> np.ndarray:
    """Mean label of the `count` most similar bank entries, per query row."""
    count = min(count, similarity.shape[-1])
    index = np.argpartition(-similarity, count - 1, axis=-1)[:, :count]
    return values[index].mean(-1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--split-directory", type=Path, required=True)
    parser.add_argument("--split", default="meta_val")
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=args.split_directory)
    fingerprints = data.fingerprints

    train_cells = [c for c in data.cells if c["split"] == "meta_train"]
    eval_cells = [c for c in data.cells if c["split"] == args.split]
    ligand_values: dict[str, list[float]] = defaultdict(list)
    target_values: dict[str, list[float]] = defaultdict(list)
    for cell in train_cells:
        ligand_values[cell["ligand_id"]].append(float(cell["pK"]))
        target_values[cell["target_id"]].append(float(cell["pK"]))
    train_ligands = sorted(ligand_values)
    train_mean = np.asarray([float(np.mean(ligand_values[k])) for k in train_ligands])
    train_targets = sorted(target_values)
    train_target_mean = np.asarray(
        [float(np.mean(target_values[t])) for t in train_targets])

    eval_ligands = sorted({c["ligand_id"] for c in eval_cells})
    all_ligands = train_ligands + eval_ligands
    component_of = {c["target_id"]: c["protein_group_40"] for c in data.cells}

    # ---------------- ligand representations -------------------------------
    morgan = np.stack([fingerprints[k].numpy() for k in all_ligands])
    ligand_similarity: dict[str, np.ndarray] = {}
    ligand_similarity["morgan"] = None                     # Tanimoto, handled below
    descriptors = rdkit_descriptors(data._ligand_smiles, all_ligands)
    train_slice = descriptors[:len(train_ligands)]
    mean, deviation = train_slice.mean(0), train_slice.std(0) + 1e-6
    ligand_vectors = {"rdkit_desc": unit((descriptors - mean) / deviation)}
    try:
        ligand_vectors["chemberta"] = unit(chemberta_embeddings(
            data._ligand_smiles, all_ligands, args.device))
        chemberta_available = True
    except Exception as error:                              # noqa: BLE001
        chemberta_available = False
        chemberta_error = f"{type(error).__name__}: {error}"

    def ligand_pair_similarity(arm: str, left: np.ndarray,
                               right: np.ndarray) -> np.ndarray:
        if arm == "morgan":
            return tanimoto_rows(morgan[left], morgan[right])
        return ligand_vectors[arm][left] @ ligand_vectors[arm][right].T

    ligand_index = {key: i for i, key in enumerate(all_ligands)}
    train_index = np.arange(len(train_ligands))
    ligand_arms = ["morgan", "rdkit_desc"] + (
        ["chemberta"] if chemberta_available else [])

    # ---------------- protein representations ------------------------------
    eval_targets = sorted({c["target_id"] for c in eval_cells})
    pooled = {t: np.asarray(data.protein_for_target(t)[0], dtype=np.float32)
              for t in (*train_targets, *eval_targets)}
    train_pooled = np.stack([pooled[t] for t in train_targets])
    eval_pooled = np.stack([pooled[t] for t in eval_targets])
    center = train_pooled.mean(0, keepdims=True)
    deviation_matrix = train_pooled - center
    covariance = deviation_matrix.T @ deviation_matrix / max(len(train_pooled) - 1, 1)
    values, vectors = np.linalg.eigh(covariance.astype(np.float64))
    whiten = ((vectors / np.sqrt(np.maximum(values, 1e-3))) @ vectors.T).astype(np.float32)
    train_kmer = kmer_matrix([data._protein_sequences[t] for t in train_targets])
    eval_kmer_full = kmer_matrix(
        [data._protein_sequences[t] for t in (*train_targets, *eval_targets)])
    protein_bank = {
        "esm_pooled": (unit(train_pooled), unit(eval_pooled)),
        "esm_centered": (unit(train_pooled - center), unit(eval_pooled - center)),
        "esm_whitened": (unit((train_pooled - center) @ whiten.T),
                         unit((eval_pooled - center) @ whiten.T)),
        "kmer3": (unit(eval_kmer_full[:len(train_targets)]),
                  unit(eval_kmer_full[len(train_targets):])),
    }
    del train_kmer

    # ---------------- evaluation -------------------------------------------
    by_target: dict[str, list[dict]] = defaultdict(list)
    for cell in eval_cells:
        by_target[cell["target_id"]].append(cell)

    similarity_to_train = {}
    for target, cells in by_target.items():
        rows = np.stack([fingerprints[c["ligand_id"]].numpy() for c in cells])
        similarity_to_train[target] = tanimoto_rows(
            rows, morgan[:len(train_ligands)]).max(-1)

    ligand_report: dict[str, dict] = {}
    for arm in ligand_arms:
        continuity, low_tier, overall, ci_rows, rho_rows = [], [], [], [], []
        robustness = {n: [] for n in ROBUSTNESS_NEIGHBORS}
        for target, cells in by_target.items():
            component = component_of[target]
            index = np.asarray([ligand_index[c["ligand_id"]] for c in cells])
            truth = np.asarray([float(c["pK"]) for c in cells])
            if len(cells) >= 3:
                pair = ligand_pair_similarity(arm, index, index)
                rows, columns = np.triu_indices(len(cells), 1)
                gap = -np.abs(truth[rows] - truth[columns])
                value = spearman(pair[rows, columns], gap)
                if value is not None and np.isfinite(value):
                    continuity.append((component, target, float(value)))
            bank = ligand_pair_similarity(arm, index, train_index)
            prediction = neighbor_mean(bank, train_mean, LIGAND_NEIGHBORS)
            squared = (prediction - truth) ** 2
            low = similarity_to_train[target] < LOW_TIER
            overall.extend((component, target, float(v)) for v in squared)
            low_tier.extend((component, target, float(v)) for v in squared[low])
            for count in ROBUSTNESS_NEIGHBORS:
                other = neighbor_mean(bank, train_mean, count)
                robustness[count].extend(
                    (component, target, float(v)) for v in (other - truth) ** 2)
            if len(cells) >= 2:
                value, comparable = concordance_index(prediction, truth)
                if comparable:
                    ci_rows.append((component, target, value))
                rho = spearman(prediction, truth)
                if rho is not None and np.isfinite(rho):
                    rho_rows.append((component, target, float(rho)))
        ligand_report[arm] = {
            "L1_continuity": component_target_mean(continuity),
            "L1_bootstrap": component_bootstrap(
                continuity, args.bootstrap_draws, 20260815),
            "L2_low_tier_mse": component_target_mean(low_tier),
            "L2_all_mse": component_target_mean(overall),
            "L3_ci": component_target_mean(ci_rows),
            "L3_spearman": component_target_mean(rho_rows),
            "robustness_mse": {str(n): component_target_mean(robustness[n])
                               for n in ROBUSTNESS_NEIGHBORS},
            "low_tier_cells": len(low_tier), "cells": len(overall),
        }

    protein_report: dict[str, dict] = {}
    rng = np.random.default_rng(20260815)
    for arm, (bank, query) in protein_bank.items():
        similarity = query @ bank.T
        shuffled = similarity[:, rng.permutation(similarity.shape[1])]
        level = neighbor_mean(similarity, train_target_mean, PROTEIN_NEIGHBORS)
        level_shuffled = neighbor_mean(shuffled, train_target_mean, PROTEIN_NEIGHBORS)
        correct, control = [], []
        for position, target in enumerate(eval_targets):
            component = component_of[target]
            truth = np.asarray([float(c["pK"]) for c in by_target[target]])
            correct.extend((component, target, float((level[position] - v) ** 2))
                           for v in truth)
            control.extend(
                (component, target, float((level_shuffled[position] - v) ** 2))
                for v in truth)
        p1 = component_target_mean(correct)
        protein_report[arm] = {
            "P1_level_mse": p1,
            "P2_specificity": component_target_mean(control) - p1,
            "P2_bootstrap": component_bootstrap(
                [(c, t, control[i][2] - correct[i][2])
                 for i, (c, t, _) in enumerate(correct)],
                args.bootstrap_draws, 20260815),
            "top16_similarity_spread": float(np.mean(
                np.sort(similarity, -1)[:, -PROTEIN_NEIGHBORS:].ptp(-1))),
        }

    # ---------------- decision ---------------------------------------------
    incumbent_ligand, incumbent_protein = "morgan", "esm_pooled"
    ligand_choice = max(ligand_arms, key=lambda a: ligand_report[a]["L1_continuity"])
    if (ligand_report[ligand_choice]["L2_low_tier_mse"]
            > ligand_report[incumbent_ligand]["L2_low_tier_mse"]):
        ligand_choice = incumbent_ligand
    protein_choice = min(protein_bank, key=lambda a: protein_report[a]["P1_level_mse"])
    if protein_report[protein_choice]["P2_specificity"] <= 0:
        protein_choice = incumbent_protein

    payload = {
        "schema": "MetaSieve.StageR2RepresentationDiscriminator.v1",
        "split_directory": str(args.split_directory), "split": args.split,
        "population": {"targets": len(by_target), "cells": len(eval_cells),
                       "components": len({component_of[t] for t in by_target})},
        "fixed_neighbors": {"ligand": LIGAND_NEIGHBORS,
                            "protein": PROTEIN_NEIGHBORS},
        "ligand": ligand_report, "protein": protein_report,
        "chemberta_available": chemberta_available,
        **({} if chemberta_available else {"chemberta_error": chemberta_error}),
        "decision": {
            "ligand": ligand_choice, "protein": protein_choice,
            "ligand_incumbent": incumbent_ligand,
            "protein_incumbent": incumbent_protein,
            "ligand_changed": ligand_choice != incumbent_ligand,
            "protein_changed": protein_choice != incumbent_protein,
        },
        "external_data_disclosure": {
            "chemberta": "DeepChem/ChemBERTa-77M-MLM, masked-LM pretraining on "
                         "public SMILES; no affinity labels",
            "esm": "ESM2-t30-150M, already part of the retained pipeline",
        },
        "not_run": ["facebook/esm2_t33_650M_UR50D (incomplete local cache)",
                    "MoLFormer, GraphMVP, PMMR, TM-Vec (not available offline)",
                    "trained GINE (every checkpoint saw double-cold meta_val "
                    "ligands under the old split)"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")

    print("ligand arms (%d targets, %d cells, %d low-tier cells)" % (
        len(by_target), len(eval_cells), ligand_report["morgan"]["low_tier_cells"]))
    print("%-14s %10s %10s %10s %8s %8s" % (
        "arm", "L1 contin", "L2 lt40", "L2 all", "L3 CI", "L3 rho"))
    for arm in ligand_arms:
        entry = ligand_report[arm]
        print("%-14s %10.4f %10.4f %10.4f %8.4f %8.4f" % (
            arm, entry["L1_continuity"], entry["L2_low_tier_mse"],
            entry["L2_all_mse"], entry["L3_ci"], entry["L3_spearman"]))
    print("\n%-14s %10s %12s %14s" % ("arm", "P1 level", "P2 specific", "top16 spread"))
    for arm, entry in protein_report.items():
        print("%-14s %10.4f %12.4f %14.4f" % (
            arm, entry["P1_level_mse"], entry["P2_specificity"],
            entry["top16_similarity_spread"]))
    print(f"\ndecision: ligand={ligand_choice} protein={protein_choice}")


if __name__ == "__main__":
    main()
