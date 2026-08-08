"""Train the isolated P1B contact/distance bridge; no affinity or CSMO code."""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import random

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM, MAX_ATOMS
from contracts.mechanism import (CONTACT_THRESHOLD_ANGSTROM,
    DISTANCE_BINS_ANGSTROM, MECHANISM_RESIDUE_SLOTS, MECHANISM_SCHEMA)
from model.encoders import LigandEncoder, ProteinEncoder
from model.mechanism import MechanisticInteractionBridge
from scripts.data_contract import read_jsonl, write_jsonl
from scripts.structure_sources.rcsb import sha256_file


@dataclass(frozen=True)
class TrainConfig:
    seed: int = 17
    hidden_dim: int = 128
    bridge_rank: int = 32
    gine_layers: int = 3
    batch_size: int = 8
    epochs: int = 5
    learning_rate: float = 3e-4
    weight_decay: float = 1e-4
    distance_loss_weight: float = 1.0


class MechanismPretrainer(nn.Module):
    def __init__(self, protein_dim: int, config: TrainConfig):
        super().__init__()
        self.protein = ProteinEncoder(protein_dim, config.hidden_dim, dtype=torch.float32)
        self.ligand = LigandEncoder(config.hidden_dim, n_layers=config.gine_layers,
                                    dtype=torch.float32)
        self.bridge = MechanisticInteractionBridge(
            config.hidden_dim, config.hidden_dim, rank=config.bridge_rank,
            dtype=torch.float32)

    def forward(self, X, A, atom_mask, protein_pooled, protein_residues, residue_mask):
        _, atoms = self.ligand(X, A, atom_mask)
        _, residues = self.protein(protein_pooled, protein_residues)
        return self.bridge(atoms, atom_mask, residues, residue_mask)


class MechanismCorpus:
    """Memory-resident immutable banks for deterministic single-GPU training."""

    def __init__(self, records_path: str | Path, supervision_dir: str | Path,
                 protein_bank_dir: str | Path, ligand_bank_path: str | Path):
        self.records = read_jsonl(records_path)
        self.index_by_entry = {record["source_entry_id"]: index
                               for index, record in enumerate(self.records)}
        if len(self.index_by_entry) != len(self.records):
            raise ValueError("duplicate source_entry_id in structure split")

        protein_root = Path(protein_bank_dir)
        protein_manifest = json.loads((protein_root / "manifest.json").read_text(encoding="utf-8"))
        self.protein_dim = int(protein_manifest["hidden_dim"])
        self.protein_shards: list[dict[str, np.ndarray]] = []
        self.protein_lookup: dict[str, tuple[int, int]] = {}
        for item in protein_manifest["shards"]:
            path = protein_root / item["path"]
            if sha256_file(path) != item["sha256"]:
                raise ValueError(f"protein shard hash mismatch: {path}")
            with np.load(path, allow_pickle=False) as shard:
                arrays = {key: shard[key].copy() for key in
                          ("keys", "pooled", "residues", "mask")}
            shard_index = len(self.protein_shards)
            self.protein_shards.append(arrays)
            for index, key in enumerate(arrays["keys"]):
                self.protein_lookup[str(key)] = (shard_index, index)

        supervision_root = Path(supervision_dir)
        pairs = read_jsonl(supervision_root / "pairs.jsonl")
        self.pairs = {row["source_entry_id"]: row for row in pairs}
        if set(self.pairs) != set(self.index_by_entry):
            raise ValueError("split records and geometry pairs differ")
        manifest = json.loads((supervision_root / "manifest.json").read_text(encoding="utf-8"))
        self.geometry_shards = {}
        for item in manifest["shards"]:
            path = supervision_root / item["path"]
            if sha256_file(path) != item["sha256"]:
                raise ValueError(f"geometry shard hash mismatch: {path}")
            with np.load(path, allow_pickle=False) as shard:
                self.geometry_shards[item["path"]] = {
                    key: shard[key].copy() for key in
                    ("contact", "distance_bin", "atom_mask", "residue_mask")}
        self.ligands = torch.load(ligand_bank_path, map_location="cpu", weights_only=False)
        self.split_indices = {name: [index for index, record in enumerate(self.records)
                                     if record["source_split"] == name]
                              for name in ("train", "val", "test")}
        if any(not values for values in self.split_indices.values()):
            raise ValueError("structure corpus must contain train, val, and test records")

    def geometry(self, index: int) -> dict[str, np.ndarray]:
        pair = self.pairs[self.records[index]["source_entry_id"]]
        shard = self.geometry_shards[pair["shard"]]
        return {key: value[pair["shard_index"]] for key, value in shard.items()}

    def batch(self, indices: list[int], device: str, *,
              protein_indices: list[int] | None = None,
              ligand_indices: list[int] | None = None) -> dict[str, torch.Tensor]:
        protein_indices = indices if protein_indices is None else protein_indices
        ligand_indices = indices if ligand_indices is None else ligand_indices
        if len(protein_indices) != len(indices) or len(ligand_indices) != len(indices):
            raise ValueError("control index lists must match the label batch")
        batch_size = len(indices)
        X = torch.zeros(batch_size, MAX_ATOMS, ATOM_FEAT_DIM)
        A = torch.zeros(batch_size, MAX_ATOMS, MAX_ATOMS, BOND_FEAT_DIM)
        ligand_mask = torch.zeros(batch_size, MAX_ATOMS)
        pooled, residues, geometry = [], [], []
        for batch_index, (index, protein_index, ligand_index) in enumerate(zip(
                indices, protein_indices, ligand_indices)):
            record = self.records[index]
            graph = self.ligands[self.records[ligand_index]["ccd_sha256"]]
            X[batch_index] = graph["X"]
            ligand_mask[batch_index] = graph["mask"]
            edge = graph["edge_index"].long()
            A[batch_index, edge[0], edge[1]] = graph["edge_attr"]
            protein_shard, protein_row = self.protein_lookup[
                self.records[protein_index]["sequence_sha256"]]
            protein = self.protein_shards[protein_shard]
            pooled.append(protein["pooled"][protein_row])
            residues.append(protein["residues"][protein_row])
            geometry.append(self.geometry(index))
        result = {
            "X": X, "A": A, "atom_mask": ligand_mask,
            "protein_pooled": torch.from_numpy(np.stack(pooled)).float(),
            "protein_residues": torch.from_numpy(np.stack(residues)).float(),
            "residue_mask": torch.from_numpy(np.stack(
                [value["residue_mask"] for value in geometry])).float(),
            "contact": torch.from_numpy(np.stack(
                [value["contact"] for value in geometry])).float(),
            "distance": torch.from_numpy(np.stack(
                [value["distance_bin"] for value in geometry])).long(),
        }
        if not torch.equal(result["atom_mask"].to(torch.uint8), torch.from_numpy(np.stack(
                [value["atom_mask"] for value in geometry]))):
            raise ValueError("CCD graph atom mask differs from geometry atom mask")
        return {key: value.to(device, non_blocking=True) for key, value in result.items()}

    def class_weights(self, indices: list[int]) -> tuple[float, torch.Tensor, dict]:
        positive = negative = 0
        distance_counts = torch.zeros(len(DISTANCE_BINS_ANGSTROM) - 1, dtype=torch.long)
        for index in indices:
            value = self.geometry(index)
            valid = value["atom_mask"][:, None] * value["residue_mask"][None, :]
            contacts = value["contact"][valid.astype(bool)]
            positive += int(contacts.sum())
            negative += int(len(contacts) - contacts.sum())
            bins = torch.from_numpy(value["distance_bin"][valid.astype(bool)].astype(np.int64))
            distance_counts += torch.bincount(bins, minlength=len(distance_counts))
        if positive == 0 or (distance_counts == 0).any():
            raise ValueError("training labels do not cover contact and every distance class")
        pos_weight = min(negative / positive, 20.0)
        distance_weights = distance_counts.sum() / (len(distance_counts) * distance_counts.float())
        audit = {"contact_positive": positive, "contact_negative": negative,
                 "contact_pos_weight": pos_weight,
                 "distance_counts": distance_counts.tolist(),
                 "distance_weights": distance_weights.tolist()}
        return pos_weight, distance_weights, audit


def _loss(output, batch, pos_weight: torch.Tensor, distance_weights: torch.Tensor,
          distance_scale: float) -> tuple[torch.Tensor, dict[str, float]]:
    valid = output.pair_mask.bool()
    contact = F.binary_cross_entropy_with_logits(
        output.contact_logits[valid], batch["contact"][valid], pos_weight=pos_weight)
    distance = F.cross_entropy(output.distance_logits[valid], batch["distance"][valid],
                               weight=distance_weights)
    total = contact + distance_scale * distance
    return total, {"contact_loss": float(contact.detach()),
                   "distance_loss": float(distance.detach()), "loss": float(total.detach())}


def _evaluate_loss(model, corpus, indices, batch_size, device, pos_weight, distance_weights,
                   distance_scale) -> dict:
    totals = {"loss": 0.0, "contact_loss": 0.0, "distance_loss": 0.0}
    count = 0
    model.eval()
    with torch.inference_mode():
        for start in range(0, len(indices), batch_size):
            batch_indices = indices[start:start + batch_size]
            batch = corpus.batch(batch_indices, device)
            with torch.autocast(device_type="cuda", dtype=torch.float16,
                                enabled=device.startswith("cuda")):
                output = model(batch["X"], batch["A"], batch["atom_mask"],
                               batch["protein_pooled"], batch["protein_residues"],
                               batch["residue_mask"])
                _, metrics = _loss(output, batch, pos_weight, distance_weights, distance_scale)
            for key in totals:
                totals[key] += metrics[key] * len(batch_indices)
            count += len(batch_indices)
    return {key: value / count for key, value in totals.items()}


def train_mechanistic_bridge(records_path: str | Path, supervision_dir: str | Path,
                             protein_bank_dir: str | Path, ligand_bank_path: str | Path,
                             output_dir: str | Path, config: TrainConfig, *,
                             device: str = "cuda", max_train_records: int | None = None,
                             max_eval_records: int | None = None) -> dict:
    if not device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("P1B is registered for CUDA and fails closed without it")
    random.seed(config.seed)
    np.random.seed(config.seed)
    torch.manual_seed(config.seed)
    torch.cuda.manual_seed_all(config.seed)
    output = Path(output_dir)
    if output.exists():
        raise FileExistsError(f"P1B output already exists: {output}")
    output.mkdir(parents=True)
    corpus = MechanismCorpus(records_path, supervision_dir, protein_bank_dir, ligand_bank_path)
    train_indices = list(corpus.split_indices["train"])
    if max_train_records is not None:
        train_indices = train_indices[:max_train_records]
    val_indices = list(corpus.split_indices["val"])
    if max_eval_records is not None:
        val_indices = val_indices[:max_eval_records]
    pos_weight_value, distance_weights_cpu, class_audit = corpus.class_weights(train_indices)
    pos_weight = torch.tensor(pos_weight_value, device=device)
    distance_weights = distance_weights_cpu.to(device)
    model = MechanismPretrainer(corpus.protein_dim, config).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=config.learning_rate,
                                  weight_decay=config.weight_decay)
    scaler = torch.amp.GradScaler("cuda")
    history, best_loss = [], math.inf
    generator = torch.Generator().manual_seed(config.seed)
    for epoch in range(config.epochs):
        order = torch.randperm(len(train_indices), generator=generator).tolist()
        model.train()
        totals = {"loss": 0.0, "contact_loss": 0.0, "distance_loss": 0.0}
        seen = 0
        for start in range(0, len(order), config.batch_size):
            indices = [train_indices[value] for value in order[start:start + config.batch_size]]
            batch = corpus.batch(indices, device)
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type="cuda", dtype=torch.float16):
                prediction = model(batch["X"], batch["A"], batch["atom_mask"],
                                   batch["protein_pooled"], batch["protein_residues"],
                                   batch["residue_mask"])
                loss, metrics = _loss(prediction, batch, pos_weight, distance_weights,
                                      config.distance_loss_weight)
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(optimizer)
            scaler.update()
            for key in totals:
                totals[key] += metrics[key] * len(indices)
            seen += len(indices)
        train_metrics = {key: value / seen for key, value in totals.items()}
        val_metrics = _evaluate_loss(
            model, corpus, val_indices, config.batch_size, device,
            pos_weight, distance_weights, config.distance_loss_weight)
        row = {"epoch": epoch + 1, "train": train_metrics, "val": val_metrics}
        history.append(row)
        if val_metrics["loss"] < best_loss:
            best_loss = val_metrics["loss"]
            torch.save({"schema": "MetaSieve.MechanismCheckpoint.v1",
                        "model_state": model.state_dict(), "protein_dim": corpus.protein_dim,
                        "config": asdict(config), "epoch": epoch + 1,
                        "val_loss": best_loss}, output / "best.pt")
    write_jsonl(output / "history.jsonl", history)
    mechanism_schema = {"schema": MECHANISM_SCHEMA,
                        "contact_threshold_angstrom": CONTACT_THRESHOLD_ANGSTROM,
                        "distance_bins_angstrom": list(DISTANCE_BINS_ANGSTROM),
                        "residue_slots": MECHANISM_RESIDUE_SLOTS}
    (output / "mechanism_schema.json").write_text(
        json.dumps(mechanism_schema, indent=2, sort_keys=True), encoding="utf-8")
    manifest = {
        "schema": "MetaSieve.MechanismPretrainingRun.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config": asdict(config), "class_balance": class_audit,
        "train_records": len(train_indices),
        "val_records": len(val_indices),
        "test_records": len(corpus.split_indices["test"]),
        "best_val_loss": best_loss,
        "checkpoint_sha256": sha256_file(output / "best.pt"),
        "inputs": {
            "records": sha256_file(records_path),
            "supervision_manifest": sha256_file(Path(supervision_dir) / "manifest.json"),
            "protein_manifest": sha256_file(Path(protein_bank_dir) / "manifest.json"),
            "ligand_bank": sha256_file(ligand_bank_path),
        },
        "affinity_labels_used": False, "csmo_used": False,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("records")
    parser.add_argument("supervision_dir")
    parser.add_argument("protein_bank_dir")
    parser.add_argument("ligand_bank")
    parser.add_argument("output")
    parser.add_argument("--epochs", type=int, default=5)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--max-train-records", type=int)
    parser.add_argument("--max-eval-records", type=int)
    args = parser.parse_args()
    config = TrainConfig(seed=args.seed, epochs=args.epochs, batch_size=args.batch_size)
    result = train_mechanistic_bridge(
        args.records, args.supervision_dir, args.protein_bank_dir,
        args.ligand_bank, args.output, config, device=args.device,
        max_train_records=args.max_train_records,
        max_eval_records=args.max_eval_records)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
