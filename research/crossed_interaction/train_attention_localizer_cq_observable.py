"""Train a ligand-conditioned ESM2 slot attention observable on Ki quotients.

This is a source-only admission Gate. A small differentiable localizer is
trained only on train quotient blocks, then frozen. The neural readout used for
localizer training is discarded; final development admission still uses the
shared positive-ridge quotient observable and wrong-partner controls.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch

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


ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/crossed_interaction/bindingdb_202608/cq_ki_corpus"
OUT = CQ_OUT.parent / "attention_localizer_cq_observable_gate1"


@dataclass(frozen=True)
class TrainBlock:
    panel_id: str
    protein_keys: list[str]
    ligand_keys: list[str]
    y: np.ndarray
    residual_matrix: np.ndarray


class SlotAttentionObservable(torch.nn.Module):
    def __init__(self, slot_dim: int, ligand_dim: int, attention_dim: int):
        super().__init__()
        self.slot_projection = torch.nn.Linear(slot_dim, attention_dim, bias=False)
        self.ligand_projection = torch.nn.Linear(ligand_dim, attention_dim, bias=False)
        self.score = torch.nn.Linear(attention_dim, 1, bias=False)
        self.readout = torch.nn.Linear(slot_dim * ligand_dim, 1, bias=False)

    def attention_weights(
            self, slots: torch.Tensor, mask: torch.Tensor,
            attention_ligand: torch.Tensor) -> torch.Tensor:
        hidden = self.slot_projection(slots) + self.ligand_projection(attention_ligand)[:, None, :]
        logits = self.score(torch.tanh(hidden)).squeeze(-1)
        logits = logits.masked_fill(~mask, -1.0e9)
        return torch.softmax(logits, dim=1)

    def feature(
            self, slots: torch.Tensor, mask: torch.Tensor,
            attention_ligand: torch.Tensor, feature_ligand: torch.Tensor) -> torch.Tensor:
        weights = self.attention_weights(slots, mask, attention_ligand)
        pooled = torch.sum(weights[:, :, None] * slots, dim=1)
        return torch.einsum("nh,nl->nhl", pooled, feature_ligand).flatten(start_dim=1)

    def forward(
            self, slots: torch.Tensor, mask: torch.Tensor,
            attention_ligand: torch.Tensor, feature_ligand: torch.Tensor) -> torch.Tensor:
        return self.readout(self.feature(slots, mask, attention_ligand, feature_ligand)).squeeze(-1)


def projection_residual_matrix(target_ids: list[str], ligand_ids: list[str]) -> np.ndarray:
    targets = {value: index for index, value in enumerate(sorted(set(target_ids)))}
    ligands = {value: index for index, value in enumerate(sorted(set(ligand_ids)))}
    design = np.zeros((len(target_ids), 1 + len(targets) + len(ligands)), dtype=np.float64)
    for row, (target, ligand) in enumerate(zip(target_ids, ligand_ids)):
        design[row, 0] = 1.0
        design[row, 1 + targets[target]] = 1.0
        design[row, 1 + len(targets) + ligands[ligand]] = 1.0
    return np.eye(len(target_ids), dtype=np.float64) - design @ np.linalg.pinv(design)


def materialize_inputs(
        corpus: Path, protein_bank: Path, *,
        hidden_blocks: int) -> tuple[
            dict[str, dict], list[dict], dict[str, np.ndarray], dict[str, np.ndarray],
            dict[str, np.ndarray], dict]:
    cells_list = read_jsonl_gz(corpus / "cells.jsonl.gz")
    cells = {cell["cell_id"]: cell for cell in cells_list}
    panels = read_jsonl_gz(corpus / "panels.jsonl.gz")
    proteins_json = read_jsonl(corpus / "proteins.jsonl")
    protein_keys = {row["sequence_sha256"] for row in proteins_json}
    protein_rows, protein_manifest = load_protein_bank(protein_bank, protein_keys)
    proteins = {}
    protein_masks = {}
    for key, row in protein_rows.items():
        proteins[key] = protein_slot_blocks(
            row["residues"], row["mask"], hidden_blocks=hidden_blocks)
        protein_masks[key] = row["mask"].astype(bool)
    ligands = {
        row["drug_key"]: ligand_pharmacophore_descriptor(row["smiles"])
        for row in read_jsonl(corpus / "ligands.jsonl")
    }
    metadata = {
        "cells": len(cells_list),
        "panels": len(panels),
        "proteins": len(proteins),
        "ligands": len(ligands),
        "protein_bank_model_id": protein_manifest.get("model_id", ""),
        "protein_bank_model_revision": protein_manifest.get("model_revision", ""),
        "protein_bank_slot_policy": protein_manifest.get("slot_policy", ""),
    }
    return cells, panels, proteins, protein_masks, ligands, metadata


def build_train_blocks(cells: dict[str, dict], panels: list[dict]) -> tuple[list[TrainBlock], float]:
    blocks = []
    max_orthogonality = 0.0
    for panel in panels:
        if panel["split"] != "train":
            continue
        ordered = [cells[cell_id] for cell_id in panel["cell_ids"]]
        target_ids = [cell["target_id"] for cell in ordered]
        ligand_ids = [cell["ligand_id"] for cell in ordered]
        y_raw = np.asarray([cell["pK"] for cell in ordered], dtype=np.float64)
        y, retained_rank, y_orthogonality = additive_residual(target_ids, ligand_ids, y_raw)
        if retained_rank <= 0:
            continue
        blocks.append(TrainBlock(
            panel_id=panel["panel_id"],
            protein_keys=target_ids,
            ligand_keys=ligand_ids,
            y=np.asarray(y, dtype=np.float64),
            residual_matrix=projection_residual_matrix(target_ids, ligand_ids),
        ))
        max_orthogonality = max(max_orthogonality, y_orthogonality)
    if not blocks:
        raise ValueError("no train quotient panels available")
    return blocks, max_orthogonality


def shuffled_ligand_map(ligands: dict[str, np.ndarray], seed: int) -> dict[str, str]:
    keys = sorted(ligands)
    if len(keys) < 2:
        raise ValueError("ligand-shuffled control needs at least two ligands")
    offset = int(np.random.default_rng(seed).integers(1, len(keys)))
    return {key: keys[(index + offset) % len(keys)] for index, key in enumerate(keys)}


def _attention_ligand_keys(
        ligand_keys: list[str], mode: str, ligand_shuffle: dict[str, str]) -> list[str]:
    if mode == "protein_only":
        return ["__zero__"] * len(ligand_keys)
    if mode == "ligand_shuffled":
        return [ligand_shuffle[key] for key in ligand_keys]
    return ligand_keys


def _torch_batch(
        protein_keys: list[str], ligand_keys: list[str], proteins: dict[str, np.ndarray],
        protein_masks: dict[str, np.ndarray], ligands: dict[str, np.ndarray],
        mode: str, ligand_shuffle: dict[str, str], device: torch.device) -> tuple[
            torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    attention_keys = _attention_ligand_keys(ligand_keys, mode, ligand_shuffle)
    slot_values = torch.tensor(
        np.stack([proteins[key] for key in protein_keys]),
        dtype=torch.float32, device=device)
    masks = torch.tensor(
        np.stack([protein_masks[key] for key in protein_keys]),
        dtype=torch.bool, device=device)
    feature_ligands = torch.tensor(
        np.stack([ligands[key] for key in ligand_keys]),
        dtype=torch.float32, device=device)
    zero_ligand = np.zeros_like(next(iter(ligands.values())))
    attention_ligands = torch.tensor(
        np.stack([
            zero_ligand if key == "__zero__" else ligands[key]
            for key in attention_keys
        ]),
        dtype=torch.float32, device=device)
    return slot_values, masks, attention_ligands, feature_ligands


def train_attention_localizer(
        blocks: list[TrainBlock], proteins: dict[str, np.ndarray],
        protein_masks: dict[str, np.ndarray], ligands: dict[str, np.ndarray], *,
        mode: str, attention_dim: int, epochs: int, learning_rate: float,
        weight_decay: float, seed: int, device: str) -> tuple[SlotAttentionObservable, dict]:
    if mode not in {"attention", "ligand_shuffled", "protein_only"}:
        raise ValueError(f"trainable localizer mode expected, got {mode}")
    torch.manual_seed(seed)
    np.random.seed(seed)
    device_value = torch.device(device)
    slot_dim = int(next(iter(proteins.values())).shape[1])
    ligand_dim = int(next(iter(ligands.values())).shape[0])
    model = SlotAttentionObservable(slot_dim, ligand_dim, attention_dim).to(device_value)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    ligand_shuffle = shuffled_ligand_map(ligands, seed) if mode == "ligand_shuffled" else {}
    history = []
    for epoch in range(epochs):
        total_loss = 0.0
        total_rows = 0
        for block in blocks:
            optimizer.zero_grad(set_to_none=True)
            slots, masks, attention_ligands, feature_ligands = _torch_batch(
                block.protein_keys, block.ligand_keys, proteins, protein_masks,
                ligands, mode, ligand_shuffle, device_value)
            raw_prediction = model(slots, masks, attention_ligands, feature_ligands)
            residual_matrix = torch.tensor(
                block.residual_matrix, dtype=torch.float32, device=device_value)
            y = torch.tensor(block.y, dtype=torch.float32, device=device_value)
            prediction = residual_matrix @ raw_prediction
            loss = torch.mean(torch.square(prediction - y))
            loss.backward()
            optimizer.step()
            total_loss += float(loss.detach().cpu()) * len(block.y)
            total_rows += len(block.y)
        mean_loss = total_loss / max(total_rows, 1)
        if epoch in {0, epochs - 1}:
            history.append({"epoch": epoch + 1, "train_row_mse": mean_loss})
    return model.eval(), {
        "mode": mode,
        "attention_dim": attention_dim,
        "epochs": epochs,
        "learning_rate": learning_rate,
        "weight_decay": weight_decay,
        "train_blocks": len(blocks),
        "loss_history": history,
    }


def uniform_feature(
        protein: np.ndarray, mask: np.ndarray, ligand: np.ndarray) -> np.ndarray:
    if np.any(mask):
        pooled = protein[mask].mean(axis=0)
    else:
        pooled = np.zeros(protein.shape[1], dtype=np.float64)
    return np.multiply.outer(pooled, ligand).ravel()


def model_feature(
        model: SlotAttentionObservable, protein: np.ndarray, mask: np.ndarray,
        attention_ligand: np.ndarray, feature_ligand: np.ndarray,
        device: torch.device) -> np.ndarray:
    with torch.inference_mode():
        feature = model.feature(
            torch.tensor(protein[None], dtype=torch.float32, device=device),
            torch.tensor(mask[None], dtype=torch.bool, device=device),
            torch.tensor(attention_ligand[None], dtype=torch.float32, device=device),
            torch.tensor(feature_ligand[None], dtype=torch.float32, device=device),
        )
    return feature.squeeze(0).detach().cpu().numpy().astype(np.float64)


def attention_entropy(
        model: SlotAttentionObservable, cells: dict[str, dict], panels: list[dict],
        proteins: dict[str, np.ndarray], protein_masks: dict[str, np.ndarray],
        ligands: dict[str, np.ndarray], *, mode: str,
        ligand_shuffle: dict[str, str], device: str) -> dict:
    if mode == "uniform" or model is None:
        slot_count = next(iter(proteins.values())).shape[0]
        return {
            "train_mean": float(np.log(slot_count)),
            "development_mean": float(np.log(slot_count)),
            "normalized_by_log_slots": 1.0,
        }
    device_value = torch.device(device)
    zero_ligand = np.zeros_like(next(iter(ligands.values())))
    rows = {"train": [], "development": []}
    with torch.inference_mode():
        for panel in panels:
            if panel["split"] not in rows:
                continue
            ordered = [cells[cell_id] for cell_id in panel["cell_ids"]]
            protein_keys = [cell["target_id"] for cell in ordered]
            ligand_keys = [cell["ligand_id"] for cell in ordered]
            if mode == "protein_only":
                attention_values = np.stack([zero_ligand for _ in ligand_keys])
            elif mode == "ligand_shuffled":
                attention_values = np.stack([ligands[ligand_shuffle[key]] for key in ligand_keys])
            else:
                attention_values = np.stack([ligands[key] for key in ligand_keys])
            weights = model.attention_weights(
                torch.tensor(
                    np.stack([proteins[key] for key in protein_keys]),
                    dtype=torch.float32, device=device_value),
                torch.tensor(
                    np.stack([protein_masks[key] for key in protein_keys]),
                    dtype=torch.bool, device=device_value),
                torch.tensor(attention_values, dtype=torch.float32, device=device_value),
            ).detach().cpu().numpy()
            entropy = -np.sum(weights * np.log(np.clip(weights, 1e-12, 1.0)), axis=1)
            rows[panel["split"]].extend(entropy.tolist())
    train_mean = float(np.mean(rows["train"]))
    development_mean = float(np.mean(rows["development"]))
    slot_count = next(iter(proteins.values())).shape[0]
    return {
        "train_mean": train_mean,
        "development_mean": development_mean,
        "normalized_by_log_slots": float(development_mean / np.log(slot_count)),
    }


def materialize_features(
        cells: dict[str, dict], proteins: dict[str, np.ndarray],
        protein_masks: dict[str, np.ndarray], ligands: dict[str, np.ndarray], *,
        mode: str, model: SlotAttentionObservable | None,
        ligand_shuffle: dict[str, str], device: str) -> dict[str, dict[str, np.ndarray]]:
    protein_donor, ligand_donor = donor_maps(list(cells.values()))
    device_value = torch.device(device)
    zero_ligand = np.zeros_like(next(iter(ligands.values())))

    def feature_for(protein_key: str, ligand_key: str) -> np.ndarray:
        if mode == "uniform":
            return uniform_feature(proteins[protein_key], protein_masks[protein_key], ligands[ligand_key])
        if model is None:
            raise ValueError("trainable localizer mode requires a model")
        if mode == "protein_only":
            attention_ligand = zero_ligand
        elif mode == "ligand_shuffled":
            attention_ligand = ligands[ligand_shuffle[ligand_key]]
        else:
            attention_ligand = ligands[ligand_key]
        return model_feature(
            model, proteins[protein_key], protein_masks[protein_key],
            attention_ligand, ligands[ligand_key], device_value)

    features = {}
    for cell_id, cell in cells.items():
        target = cell["target_id"]
        ligand = cell["ligand_id"]
        features[cell_id] = {
            "correct": feature_for(target, ligand),
            "deranged_protein": feature_for(protein_donor[target], ligand),
            "foreign_ligand": feature_for(target, ligand_donor[ligand]),
        }
    return features


def load_blocks(
        cells: dict[str, dict], panels: list[dict],
        features: dict[str, dict[str, np.ndarray]]) -> tuple[list[QuotientBlock], float]:
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
    return blocks, max_orthogonality


def run(
        corpus: Path = CORPUS, protein_bank: Path = PROTEIN_BANK,
        output: Path = OUT, ridge: float = 10000.0,
        bootstrap_draws: int = 9999, seed: int = 20260812,
        hidden_blocks: int = 8, attention_dim: int = 16,
        epochs: int = 80, learning_rate: float = 0.01,
        weight_decay: float = 0.001, localizer_mode: str = "attention",
        device: str = "cuda") -> dict:
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    if device == "cuda" and not torch.cuda.is_available():
        device = "cpu"
    if localizer_mode not in {"attention", "uniform", "ligand_shuffled", "protein_only"}:
        raise ValueError(f"unknown localizer mode: {localizer_mode}")
    cells, panels, proteins, protein_masks, ligands, input_metadata = materialize_inputs(
        corpus, protein_bank, hidden_blocks=hidden_blocks)
    train_blocks_for_localizer, localizer_projection_orthogonality = build_train_blocks(
        cells, panels)
    localizer: dict
    model = None
    ligand_shuffle = {}
    if localizer_mode == "uniform":
        localizer = {
            "mode": localizer_mode,
            "attention_dim": 0,
            "epochs": 0,
            "learning_rate": 0.0,
            "weight_decay": 0.0,
            "train_blocks": 0,
            "loss_history": [],
        }
    else:
        model, localizer = train_attention_localizer(
            train_blocks_for_localizer, proteins, protein_masks, ligands,
            mode=localizer_mode, attention_dim=attention_dim, epochs=epochs,
            learning_rate=learning_rate, weight_decay=weight_decay,
            seed=seed, device=device)
        if localizer_mode == "ligand_shuffled":
            ligand_shuffle = shuffled_ligand_map(ligands, seed)
    features = materialize_features(
        cells, proteins, protein_masks, ligands, mode=localizer_mode,
        model=model, ligand_shuffle=ligand_shuffle, device=device)
    entropy = attention_entropy(
        model, cells, panels, proteins, protein_masks, ligands,
        mode=localizer_mode, ligand_shuffle=ligand_shuffle, device=device)
    blocks, max_projection_orthogonality = load_blocks(cells, panels, features)
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
        "projection_orthogonality": max_projection_orthogonality <= 1e-7,
        "localizer_projection_orthogonality": localizer_projection_orthogonality <= 1e-7,
        "development_components_ge_5": len({
            block.dependency_component for block in development_blocks}) >= 5,
        "correct_beats_zero_additive": contrasts[0]["pass"],
        "correct_beats_deranged_protein": contrasts[1]["pass"],
        "correct_beats_foreign_ligand": contrasts[2]["pass"],
    }
    verdict = (
        "ATTENTION_LOCALIZER_CQ_GATE1_PASS_DEVELOPMENT"
        if all(gates.values())
        else "ATTENTION_LOCALIZER_CQ_GATE1_FAIL_CLOSED"
    )
    ligand_dim = int(next(iter(ligands.values())).shape[0])
    result = {
        "schema": "MetaSieve.AttentionLocalizerCQObservableGate1.v1",
        "hypothesis": (
            "A ligand-conditioned differentiable localizer over frozen ESM2 "
            "residue slots carries dependency-transferable quotient interaction "
            "signal when frozen and scored by positive ridge."),
        "literature_mechanism": {
            "plm_dta": "protein language model residue states support DTA representations",
            "attention_cpi": (
                "compound-conditioned target-region weighting can be an "
                "interaction mechanism"),
            "hodge_cycle_space": (
                "localizer training and final ridge scoring use additive "
                "target+ligand quotient residuals"),
        },
        "corpus": {
            **input_metadata,
            "blocks": len(blocks),
            "feature_source": "train_only_ligand_conditioned_esm2_slot_attention_x_ligand_estate",
            "feature_dim": int(hidden_blocks * ligand_dim),
            "protein_descriptor_dim": hidden_blocks,
            "ligand_descriptor_dim": ligand_dim,
            "hidden_blocks": hidden_blocks,
            "max_projection_orthogonality": max_projection_orthogonality,
            "corpus_manifest_sha256": sha256_file(corpus / "manifest.json"),
            "protein_bank_manifest_sha256": sha256_file(protein_bank / "manifest.json"),
        },
        "localizer": {
            **localizer,
            "train_split": "train",
            "slot_weight_inputs": (
                "protein_slots_only" if localizer_mode == "protein_only"
                else "shuffled_ligand_estate" if localizer_mode == "ligand_shuffled"
                else "none_uniform" if localizer_mode == "uniform"
                else "protein_slots_plus_ligand_estate"
            ),
            "max_localizer_projection_orthogonality": localizer_projection_orthogonality,
            "attention_entropy": entropy,
            "neural_readout_used_for_admission": False,
        },
        "config": {
            "ridge": ridge,
            "bootstrap_draws": bootstrap_draws,
            "seed": seed,
            "arms": list(ARMS),
            "localizer_mode": localizer_mode,
            "localizer_split": "train",
            "train_split": "train",
            "evaluation_split": "development",
            "device": device,
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
    parser.add_argument("--protein-bank", type=Path, default=PROTEIN_BANK)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--ridge", type=float, default=10000.0)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260812)
    parser.add_argument("--hidden-blocks", type=int, default=8)
    parser.add_argument("--attention-dim", type=int, default=16)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--learning-rate", type=float, default=0.01)
    parser.add_argument("--weight-decay", type=float, default=0.001)
    parser.add_argument(
        "--localizer-mode",
        choices=("attention", "uniform", "ligand_shuffled", "protein_only"),
        default="attention")
    parser.add_argument("--device", choices=("cuda", "cpu"), default="cuda")
    args = parser.parse_args()
    result = run(
        corpus=args.corpus, protein_bank=args.protein_bank, output=args.output,
        ridge=args.ridge, bootstrap_draws=args.bootstrap_draws, seed=args.seed,
        hidden_blocks=args.hidden_blocks, attention_dim=args.attention_dim,
        epochs=args.epochs, learning_rate=args.learning_rate,
        weight_decay=args.weight_decay, localizer_mode=args.localizer_mode,
        device=args.device)
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
