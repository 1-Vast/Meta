"""Run the preregistered label-free E0 synthetic trainability gate."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import random

import numpy as np
import torch

from scripts.data_contract import read_jsonl
from research.e0_identifiability.mechanistic_affinity import (
    LocalMechanisticAffinityPotential, e0_loss,
)
from research.e0_identifiability.metrics import concordance as _concordance
from scripts.pretrain_mechanistic_bridge import MechanismPretrainer, TrainConfig
from scripts.structure_sources.rcsb import sha256_file


SEED = 17
TASKS_PER_FOLD = 8
LIGANDS_PER_TASK = 20
EPOCHS = 60
BATCH_SIZE = 4
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 1e-4
DISTANCE_WEIGHTS = np.asarray([1.0, 0.7, 0.2, -0.2, -0.6], dtype=np.float32)


def _load_states(root: Path, ligand_keys: set[str], protein_keys: set[str]):
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    proteins, ligands = {}, {}
    for item in manifest["protein_shards"]:
        path = root / item["path"]
        if sha256_file(path) != item["sha256"]:
            raise ValueError(f"protein state shard hash mismatch: {path}")
        with np.load(path, allow_pickle=False) as shard:
            for index, raw_key in enumerate(shard["keys"]):
                key = str(raw_key)
                if key in protein_keys:
                    proteins[key] = {
                        "residues": shard["residues"][index].astype(np.float32),
                        "chemistry": shard["chemistry"][index].astype(np.float32),
                        "mask": shard["mask"][index].astype(np.float32),
                    }
    for item in manifest["ligand_shards"]:
        path = root / item["path"]
        if sha256_file(path) != item["sha256"]:
            raise ValueError(f"ligand state shard hash mismatch: {path}")
        with np.load(path, allow_pickle=False) as shard:
            offsets = shard["offsets"]
            for index, raw_key in enumerate(shard["keys"]):
                key = str(raw_key)
                if key in ligand_keys:
                    left, right = int(offsets[index]), int(offsets[index + 1])
                    ligands[key] = {
                        "atoms": shard["atoms"][left:right].astype(np.float32),
                        "chemistry": shard["chemistry"][left:right].astype(np.float32),
                        "pooled": shard["pooled"][index].astype(np.float32),
                    }
    if set(proteins) != protein_keys or set(ligands) != ligand_keys:
        raise ValueError("local-state cache does not cover synthetic selection")
    return proteins, ligands


def _select_rows(path: Path) -> list[dict]:
    grouped: dict[tuple[int, str], dict[str, dict]] = defaultdict(dict)
    task_meta = {}
    for row in read_jsonl(path):
        key = (int(row["outer_oof_fold"]), row["task_id"])
        grouped[key].setdefault(row["ligand_state_key"], row)
        task_meta[key] = row
    selected = []
    for fold in range(5):
        tasks = [key for key in sorted(grouped) if key[0] == fold
                 and len(grouped[key]) >= LIGANDS_PER_TASK][:TASKS_PER_FOLD]
        if len(tasks) != TASKS_PER_FOLD:
            raise ValueError(f"fold {fold} lacks synthetic tasks")
        for key in tasks:
            selected.extend(grouped[key][ligand] for ligand in sorted(grouped[key])[:20])
    return selected


def _load_bridge(checkpoint_path: Path, device: str):
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    model = MechanismPretrainer(int(checkpoint["protein_dim"]),
                                TrainConfig(**checkpoint["config"])).to(device)
    model.load_state_dict(checkpoint["model_state"])
    model.eval()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    return model.bridge


def _batch(rows, proteins, ligands, device, geometry=None):
    width = max(len(ligands[row["ligand_state_key"]]["atoms"]) for row in rows)
    atom_state = np.zeros((len(rows), width, 128), dtype=np.float32)
    atom_chemistry = np.zeros((len(rows), width, 40), dtype=np.float32)
    atom_mask = np.zeros((len(rows), width), dtype=np.float32)
    residue_state, residue_chemistry, residue_mask = [], [], []
    for index, row in enumerate(rows):
        ligand = ligands[row["ligand_state_key"]]
        count = len(ligand["atoms"])
        atom_state[index, :count] = ligand["atoms"]
        atom_chemistry[index, :count] = ligand["chemistry"]
        atom_mask[index, :count] = 1
        protein = proteins[row["active_protein_key"]]
        residue_state.append(protein["residues"])
        residue_chemistry.append(protein["chemistry"])
        residue_mask.append(protein["mask"])
    values = {
        "atom_state": torch.from_numpy(atom_state).to(device),
        "atom_chemistry": torch.from_numpy(atom_chemistry).to(device),
        "atom_mask": torch.from_numpy(atom_mask).to(device),
        "residue_state": torch.from_numpy(np.stack(residue_state)).to(device),
        "residue_chemistry": torch.from_numpy(np.stack(residue_chemistry)).to(device),
        "residue_mask": torch.from_numpy(np.stack(residue_mask)).to(device),
    }
    if geometry is not None:
        contact = np.zeros((len(rows), width, 128), dtype=np.float32)
        distance = np.zeros((len(rows), width, 128, 5), dtype=np.float32)
        for index, row in enumerate(rows):
            cached = geometry[row["example_id"]]
            count = cached["contact"].shape[0]
            contact[index, :count] = cached["contact"]
            distance[index, :count] = cached["distance"]
        values["contact_prob"] = torch.from_numpy(contact).to(device)
        values["distance_prob"] = torch.from_numpy(distance).to(device)
    return values


def _geometry(rows, bridge, proteins, ligands, device):
    output = {}
    with torch.inference_mode():
        for start in range(0, len(rows), BATCH_SIZE):
            batch_rows = rows[start:start + BATCH_SIZE]
            values = _batch(batch_rows, proteins, ligands, device)
            prediction = bridge(values["atom_state"], values["atom_mask"],
                                values["residue_state"], values["residue_mask"])
            contact = prediction.contact_prob.cpu().numpy().astype(np.float16)
            distance = prediction.distance_prob.cpu().numpy().astype(np.float16)
            for index, row in enumerate(batch_rows):
                count = int(values["atom_mask"][index].sum().item())
                output[row["example_id"]] = {
                    "contact": contact[index, :count], "distance": distance[index, :count]}
    return output


def _teacher(rows, proteins, ligands, geometry):
    rng = np.random.default_rng(SEED)
    weights = rng.normal(0, 1, size=(8, 6)).astype(np.float32)
    raw, baselines = [], []
    for row in rows:
        ligand = ligands[row["ligand_state_key"]]
        protein = proteins[row["protein_sequence_sha256"]]
        atom = ligand["chemistry"][:, 32:40]
        residue = protein["chemistry"]
        contact = geometry[row["example_id"]]["contact"].astype(np.float32)
        distance = geometry[row["example_id"]]["distance"].astype(np.float32)
        compatibility = np.einsum("na,ar,lr->nl", atom, weights, residue)
        distance_score = np.einsum("nld,d->nl", distance, DISTANCE_WEIGHTS)
        raw.append(float(np.sum(contact * compatibility * distance_score) /
                         max(float(contact.sum()), 1e-6)))
        mean_pharmacophore = atom.mean(0)
        baselines.append(float(0.4 * mean_pharmacophore[1]
                               + 0.2 * mean_pharmacophore[2]
                               - 0.15 * mean_pharmacophore[7]))
    raw = np.asarray(raw)
    train = np.asarray([row["outer_oof_fold"] < 4 for row in rows])
    normalized = (raw - raw[train].mean()) / max(raw[train].std(), 1e-6)
    return normalized.astype(np.float32), np.asarray(baselines, dtype=np.float32), weights


def _macro_ci(rows, labels, predictions):
    grouped = defaultdict(list)
    for index, row in enumerate(rows):
        grouped[row["task_id"]].append(index)
    return float(np.mean([_concordance(labels[indices], predictions[indices])
                          for indices in grouped.values()]))


def run_synthetic(input_root: Path, cache_root: Path, checkpoint: Path,
                  output: Path, device: str) -> dict:
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    rows = _select_rows(input_root / "rows.label_blind.jsonl")
    for index, row in enumerate(rows):
        row["example_id"] = index
        row["active_protein_key"] = row["protein_sequence_sha256"]
    proteins_by_closure = {}
    for row in rows:
        proteins_by_closure.setdefault(row["closure_component_id"],
                                       row["protein_sequence_sha256"])
    closures = sorted(proteins_by_closure)
    deranged = {}
    for index, closure in enumerate(closures):
        for offset in range(1, len(closures)):
            candidate = closures[(index + offset) % len(closures)]
            if candidate != closure:
                deranged[closure] = proteins_by_closure[candidate]
                break
    ligand_keys = {row["ligand_state_key"] for row in rows}
    protein_keys = {row["protein_sequence_sha256"] for row in rows} | set(deranged.values())
    proteins, ligands = _load_states(cache_root, ligand_keys, protein_keys)
    bridge = _load_bridge(checkpoint, device)
    correct_geometry = _geometry(rows, bridge, proteins, ligands, device)
    residuals, baselines, teacher_weights = _teacher(
        rows, proteins, ligands, correct_geometry)
    labels = baselines + residuals

    train_rows = [row for row in rows if row["outer_oof_fold"] < 4]
    model = LocalMechanisticAffinityPotential().to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=LEARNING_RATE,
                                  weight_decay=WEIGHT_DECAY)
    by_task = defaultdict(list)
    for row in train_rows:
        by_task[row["task_id"]].append(row)
    for epoch in range(EPOCHS):
        epoch_losses = []
        task_order = sorted(by_task)
        random.shuffle(task_order)
        for task_id in task_order:
            task_rows = list(by_task[task_id])
            random.shuffle(task_rows)
            for start in range(0, len(task_rows), BATCH_SIZE):
                batch_rows = task_rows[start:start + BATCH_SIZE]
                values = _batch(batch_rows, proteins, ligands, device, correct_geometry)
                indices = torch.tensor([row["example_id"] for row in batch_rows], device=device)
                output_values = model(**values).potential
                target = torch.from_numpy(residuals).to(device)[indices]
                full_labels = torch.from_numpy(labels).to(device)[indices]
                task_index = torch.zeros(len(batch_rows), dtype=torch.long, device=device)
                loss = e0_loss(output_values, target, full_labels, task_index)
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                optimizer.step()
                epoch_losses.append(float(loss.detach().cpu()))
        if (epoch + 1) % 10 == 0:
            print(f"synthetic_epoch={epoch + 1} mean_loss={np.mean(epoch_losses):.6f}",
                  flush=True)

    evaluation = [dict(row) for row in rows if row["outer_oof_fold"] == 4]
    deranged_rows = [dict(row, active_protein_key=deranged[row["closure_component_id"]])
                     for row in evaluation]
    deranged_geometry = _geometry(deranged_rows, bridge, proteins, ligands, device)
    model.eval()
    correct_scores, deranged_scores = [], []
    with torch.inference_mode():
        for batch_rows, geometry, destination in (
                (evaluation, correct_geometry, correct_scores),
                (deranged_rows, deranged_geometry, deranged_scores)):
            for start in range(0, len(batch_rows), BATCH_SIZE):
                subset = batch_rows[start:start + BATCH_SIZE]
                destination.extend(model(**_batch(
                    subset, proteins, ligands, device, geometry)).potential.cpu().tolist())
    eval_indices = np.asarray([row["example_id"] for row in evaluation])
    eval_labels = labels[eval_indices]
    eval_baseline = baselines[eval_indices]
    correct_prediction = eval_baseline + np.asarray(correct_scores)
    deranged_prediction = eval_baseline + np.asarray(deranged_scores)
    ligand_ci = _macro_ci(evaluation, eval_labels, eval_baseline)
    correct_ci = _macro_ci(evaluation, eval_labels, correct_prediction)
    deranged_ci = _macro_ci(evaluation, eval_labels, deranged_prediction)

    example = evaluation[0]
    values = _batch([example], proteins, ligands, device, correct_geometry)
    original = model(**values).potential
    permutation = torch.arange(values["atom_state"].shape[1] - 1, -1, -1, device=device)
    for key in ("atom_state", "atom_chemistry", "atom_mask", "contact_prob",
                "distance_prob"):
        values[key] = values[key][:, permutation]
    permutation_error = float(torch.abs(original - model(**values).potential).item())
    gate = {
        "correct_ci_at_least_0_80": correct_ci >= 0.80,
        "correct_minus_ligand_at_least_0_10": correct_ci - ligand_ci >= 0.10,
        "correct_minus_deranged_at_least_0_10": correct_ci - deranged_ci >= 0.10,
        "permutation_error_at_most_1e_6": permutation_error <= 1e-6,
    }
    result = {
        "schema": "MetaSieve.E0SyntheticPregate.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P1R2B-E0",
        "affinity_labels_read": False,
        "recipient_labels_read": False,
        "selection": {"seed": SEED, "tasks_per_fold": TASKS_PER_FOLD,
                      "ligands_per_task": LIGANDS_PER_TASK,
                      "train_folds": [0, 1, 2, 3], "holdout_fold": 4},
        "optimization": {"epochs": EPOCHS, "batch_size": BATCH_SIZE,
                         "learning_rate": LEARNING_RATE, "weight_decay": WEIGHT_DECAY},
        "teacher_weights": teacher_weights.tolist(),
        "metrics": {"ligand_ci": ligand_ci, "correct_ci": correct_ci,
                    "deranged_ci": deranged_ci,
                    "correct_minus_ligand": correct_ci - ligand_ci,
                    "correct_minus_deranged": correct_ci - deranged_ci,
                    "permutation_error": permutation_error},
        "gate": {**gate, "pass": all(gate.values())},
        "inputs": {"input_manifest_sha256": sha256_file(input_root / "manifest.json"),
                   "cache_manifest_sha256": sha256_file(cache_root / "manifest.json"),
                   "checkpoint_sha256": sha256_file(checkpoint)},
    }
    output.mkdir(parents=True, exist_ok=True)
    torch.save({"model_state": model.state_dict(), "result": result}, output / "model.pt")
    result["model_sha256"] = sha256_file(output / "model.pt")
    (output / "synthetic_gate.json").write_text(
        json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path,
                        default=Path("dataset/processed/source_affinity/e0_input_v1"))
    parser.add_argument("--cache", type=Path,
                        default=Path("dataset/processed/source_affinity/e0_local_states_v1"))
    parser.add_argument("--checkpoint", type=Path,
                        default=Path("report/mechanism_refactor/p1b_pilot20k_seed17_v1/best.pt"))
    parser.add_argument("--output", type=Path,
                        default=Path("research/e0_identifiability/artifacts/e0_local_map_v1"))
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    result = run_synthetic(args.input, args.cache, args.checkpoint, args.output, args.device)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["gate"]["pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
