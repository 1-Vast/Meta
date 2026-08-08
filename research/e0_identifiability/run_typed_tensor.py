"""Run the synthetic-only E0R0 typed-tensor identifiability experiment."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import random

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.stats import pearsonr, spearmanr

from research.e0_identifiability.typed_energy_tensor_contract import (
    CP_RANK, CP_TENSOR_PARAMETERS, DISTANCE_BIN_DIM, FULL_TENSOR_PARAMETERS,
    LIGAND_TYPE_DIM, RESIDUE_TYPE_DIM,
)
from research.e0_identifiability.audit_evidence import (
    _derangement, _macro_ci, _teacher_values,
)
from research.e0_identifiability.mechanistic_affinity import pairwise_rank_loss
from research.e0_identifiability.run_synthetic_pregate import (
    BATCH_SIZE, DISTANCE_WEIGHTS, EPOCHS, LEARNING_RATE, SEED, WEIGHT_DECAY,
    _geometry, _load_bridge, _load_states, _select_rows,
)
from scripts.structure_sources.rcsb import sha256_file


class FullTypedTensor(nn.Module):
    """A centered 8 x 6 x 5 energy tensor with exactly 240 parameters."""

    def __init__(self):
        super().__init__()
        self.energy = nn.Parameter(torch.zeros(
            LIGAND_TYPE_DIM, RESIDUE_TYPE_DIM, DISTANCE_BIN_DIM))

    def forward(self, tensor: torch.Tensor) -> torch.Tensor:
        return torch.einsum("nard,ard->n", tensor, self.energy)


class CPTypedTensor(nn.Module):
    """A rank-6 CP energy tensor with exactly 114 parameters."""

    def __init__(self, seed: int = SEED):
        super().__init__()
        generator = torch.Generator().manual_seed(seed)
        self.ligand = nn.Parameter(torch.randn(
            LIGAND_TYPE_DIM, CP_RANK, generator=generator) * 0.1)
        self.residue = nn.Parameter(torch.randn(
            RESIDUE_TYPE_DIM, CP_RANK, generator=generator) * 0.1)
        self.distance = nn.Parameter(torch.randn(
            DISTANCE_BIN_DIM, CP_RANK, generator=generator) * 0.1)

    def energy_tensor(self) -> torch.Tensor:
        return torch.einsum("ak,rk,dk->ard", self.ligand, self.residue, self.distance)

    def forward(self, tensor: torch.Tensor) -> torch.Tensor:
        return torch.einsum("nard,ard->n", tensor, self.energy_tensor())


def _typed_statistics(rows: list[dict], proteins: dict, ligands: dict,
                      geometry: dict) -> np.ndarray:
    values = []
    for row in rows:
        atom = ligands[row["ligand_state_key"]]["chemistry"][:, 32:40].astype(np.float64)
        residue = proteins[row["active_protein_key"]]["chemistry"].astype(np.float64)
        cached = geometry[row["example_id"]]
        contact = cached["contact"].astype(np.float64)
        distance = cached["distance"].astype(np.float64)
        denominator = max(float(contact.sum()), 1e-6)
        values.append(np.einsum(
            "na,lr,nl,nld->ard", atom, residue, contact, distance) / denominator)
    return np.asarray(values, dtype=np.float64)


def _analytic_cp_tensor(weights: np.ndarray, distance_weights: np.ndarray) -> tuple[np.ndarray, float]:
    left, singular, right = np.linalg.svd(weights.astype(np.float64), full_matrices=False)
    ligand = left * singular[None, :]
    residue = right.T
    distance = np.repeat(distance_weights[:, None], len(singular), axis=1)
    reconstructed = np.einsum("ak,rk,dk->ard", ligand, residue, distance)
    target = weights[:, :, None] * distance_weights[None, None, :]
    return reconstructed, float(np.max(np.abs(reconstructed - target)))


def _gate(ligand_ci: float, correct_ci: float, deranged_ci: float,
          permutation_error: float) -> dict:
    checks = {
        "correct_ci_at_least_0_80": correct_ci >= 0.80,
        "correct_minus_ligand_at_least_0_10": correct_ci - ligand_ci >= 0.10,
        "correct_minus_deranged_at_least_0_10": correct_ci - deranged_ci >= 0.10,
        "permutation_error_at_most_1e_6": permutation_error <= 1e-6,
    }
    return {**checks, "pass": all(checks.values())}


def _metrics(model: nn.Module, features: torch.Tensor, deranged_features: torch.Tensor,
             rows: list[dict], labels: np.ndarray, baselines: np.ndarray,
             residuals: np.ndarray, train_mask: np.ndarray,
             holdout_mask: np.ndarray) -> dict:
    with torch.inference_mode():
        prediction = model(features).detach().cpu().numpy()
        deranged = model(deranged_features).detach().cpu().numpy()
    train_rows = [row for row, keep in zip(rows, train_mask) if keep]
    holdout_rows = [row for row, keep in zip(rows, holdout_mask) if keep]
    train_ci, _ = _macro_ci(
        train_rows, labels[train_mask], baselines[train_mask] + prediction[train_mask])
    correct_ci, _ = _macro_ci(
        holdout_rows, labels[holdout_mask], baselines[holdout_mask] + prediction[holdout_mask])
    deranged_ci, _ = _macro_ci(
        holdout_rows, labels[holdout_mask], baselines[holdout_mask] + deranged[holdout_mask])
    holdout_prediction = prediction[holdout_mask]
    holdout_target = residuals[holdout_mask]
    return {
        "train_correct_ci": train_ci,
        "correct_ci": correct_ci,
        "deranged_ci": deranged_ci,
        "correct_minus_deranged": correct_ci - deranged_ci,
        "residual_pearson": float(pearsonr(holdout_prediction, holdout_target).statistic),
        "residual_spearman": float(spearmanr(holdout_prediction, holdout_target).statistic),
    }


def _train(name: str, model: nn.Module, features: torch.Tensor,
           deranged_features: torch.Tensor, rows: list[dict], labels: np.ndarray,
           baselines: np.ndarray, residuals: np.ndarray, train_mask: np.ndarray,
           holdout_mask: np.ndarray, device: str) -> tuple[dict, list[dict]]:
    random.seed(SEED)
    np.random.seed(SEED)
    torch.manual_seed(SEED)
    model = model.to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=LEARNING_RATE,
                                  weight_decay=WEIGHT_DECAY)
    by_task: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        if train_mask[index]:
            by_task[row["task_id"]].append(index)
    target_tensor = torch.from_numpy(residuals).float().to(device)
    label_tensor = torch.from_numpy(labels).float().to(device)
    trace = []
    for epoch in range(EPOCHS):
        total_losses, residual_losses, rank_losses = [], [], []
        gradient_values: dict[str, list[float]] = defaultdict(list)
        task_order = sorted(by_task)
        random.shuffle(task_order)
        for task_id in task_order:
            indices = list(by_task[task_id])
            random.shuffle(indices)
            for start in range(0, len(indices), BATCH_SIZE):
                selected = indices[start:start + BATCH_SIZE]
                batch_index = torch.tensor(selected, dtype=torch.long, device=device)
                prediction = model(features[batch_index])
                residual_loss = F.huber_loss(prediction, target_tensor[batch_index])
                task_index = torch.zeros(len(selected), dtype=torch.long, device=device)
                rank_loss = pairwise_rank_loss(
                    prediction, label_tensor[batch_index], task_index)
                loss = residual_loss + rank_loss
                optimizer.zero_grad(set_to_none=True)
                loss.backward()
                for parameter_name, parameter in model.named_parameters():
                    if parameter.grad is not None:
                        gradient_values[parameter_name].append(
                            float(parameter.grad.detach().norm().cpu()))
                optimizer.step()
                total_losses.append(float(loss.detach().cpu()))
                residual_losses.append(float(residual_loss.detach().cpu()))
                rank_losses.append(float(rank_loss.detach().cpu()))
        epoch_metrics = _metrics(
            model, features, deranged_features, rows, labels, baselines, residuals,
            train_mask, holdout_mask)
        trace.append({
            "epoch": epoch + 1,
            "loss": float(np.mean(total_losses)),
            "residual_loss": float(np.mean(residual_losses)),
            "rank_loss": float(np.mean(rank_losses)),
            "gradient_l2_mean": {key: float(np.mean(value))
                                 for key, value in gradient_values.items()},
            "gradient_l2_max": {key: float(np.max(value))
                                for key, value in gradient_values.items()},
            "parameter_l2": {key: float(value.detach().norm().cpu())
                             for key, value in model.named_parameters()},
            **epoch_metrics,
        })
    final = _metrics(model, features, deranged_features, rows, labels, baselines,
                     residuals, train_mask, holdout_mask)
    final["parameter_count"] = sum(parameter.numel() for parameter in model.parameters())
    final["model"] = name
    return {"model": model, "metrics": final}, trace


def run(input_root: Path, cache_root: Path, checkpoint: Path,
        synthetic_root: Path, output: Path, device: str) -> dict:
    gate_report = json.loads((synthetic_root / "synthetic_gate.json").read_text())
    rows = _select_rows(input_root / "rows.label_blind.jsonl")
    for index, row in enumerate(rows):
        row["example_id"] = index
        row["active_protein_key"] = row["protein_sequence_sha256"]
    deranged_rows, _ = _derangement(rows)
    ligand_keys = {row["ligand_state_key"] for row in rows}
    protein_keys = ({row["protein_sequence_sha256"] for row in rows}
                    | {row["active_protein_key"] for row in deranged_rows})
    proteins, ligands = _load_states(cache_root, ligand_keys, protein_keys)
    bridge = _load_bridge(checkpoint, device)
    correct_geometry = _geometry(rows, bridge, proteins, ligands, device)
    deranged_geometry = _geometry(deranged_rows, bridge, proteins, ligands, device)
    weights = np.asarray(gate_report["teacher_weights"], dtype=np.float64)
    raw, _, baselines = _teacher_values(rows, proteins, ligands, correct_geometry, weights)
    train_mask = np.asarray([row["outer_oof_fold"] < 4 for row in rows])
    holdout_mask = ~train_mask
    center, scale = float(raw[train_mask].mean()), max(float(raw[train_mask].std()), 1e-6)
    residuals = (raw - center) / scale
    labels = baselines + residuals

    correct = _typed_statistics(rows, proteins, ligands, correct_geometry)
    deranged = _typed_statistics(deranged_rows, proteins, ligands, deranged_geometry)
    raw_from_tensor = np.einsum(
        "nard,ar,d->n", correct, weights, DISTANCE_WEIGHTS.astype(np.float64))
    raw_reconstruction_error = float(np.max(np.abs(raw_from_tensor - raw)))
    feature_center = correct[train_mask].mean(axis=0)
    centered = correct - feature_center
    centered_deranged = deranged - feature_center
    analytic_tensor = weights[:, :, None] * DISTANCE_WEIGHTS[None, None, :] / scale
    analytic_prediction = np.einsum("nard,ard->n", centered, analytic_tensor)
    analytic_deranged = np.einsum("nard,ard->n", centered_deranged, analytic_tensor)
    analytic_error = float(np.max(np.abs(analytic_prediction - residuals)))
    cp_tensor, cp_error = _analytic_cp_tensor(weights / scale, DISTANCE_WEIGHTS)
    cp_prediction = np.einsum("nard,ard->n", centered, cp_tensor)
    cp_numeric_error = float(np.max(np.abs(cp_prediction - residuals)))

    holdout_rows = [row for row, keep in zip(rows, holdout_mask) if keep]
    ligand_ci, _ = _macro_ci(holdout_rows, labels[holdout_mask], baselines[holdout_mask])
    analytic_correct_ci, _ = _macro_ci(
        holdout_rows, labels[holdout_mask],
        baselines[holdout_mask] + analytic_prediction[holdout_mask])
    analytic_deranged_ci, _ = _macro_ci(
        holdout_rows, labels[holdout_mask],
        baselines[holdout_mask] + analytic_deranged[holdout_mask])

    features = torch.from_numpy(centered).float().to(device)
    deranged_features = torch.from_numpy(centered_deranged).float().to(device)
    full, full_trace = _train(
        "full_240", FullTypedTensor(), features, deranged_features, rows, labels,
        baselines, residuals, train_mask, holdout_mask, device)
    cp, cp_trace = _train(
        "cp_rank_6", CPTypedTensor(), features, deranged_features, rows, labels,
        baselines, residuals, train_mask, holdout_mask, device)

    example = rows[np.flatnonzero(holdout_mask)[0]]
    atom = ligands[example["ligand_state_key"]]["chemistry"][:, 32:40].astype(np.float64)
    residue = proteins[example["active_protein_key"]]["chemistry"].astype(np.float64)
    cached = correct_geometry[example["example_id"]]
    contact = cached["contact"].astype(np.float64)
    distance = cached["distance"].astype(np.float64)
    original = np.einsum("na,lr,nl,nld->ard", atom, residue, contact, distance)
    permuted = np.einsum("na,lr,nl,nld->ard", atom[::-1], residue,
                         contact[::-1], distance[::-1])
    permutation_error = float(np.max(np.abs(original - permuted)))

    for arm in (full, cp):
        metrics = arm["metrics"]
        metrics["ligand_ci"] = ligand_ci
        metrics["correct_minus_ligand"] = metrics["correct_ci"] - ligand_ci
        metrics["permutation_error"] = permutation_error
        metrics["gate"] = _gate(
            ligand_ci, metrics["correct_ci"], metrics["deranged_ci"], permutation_error)

    analytic_gate = _gate(
        ligand_ci, analytic_correct_ci, analytic_deranged_ci, permutation_error)
    if raw_reconstruction_error > 1e-6 or cp_error > 1e-10 or permutation_error > 1e-6:
        verdict = "ANALYTIC_REALIZATION_IMPLEMENTATION_DEFECT"
    elif full["metrics"]["gate"]["pass"] and cp["metrics"]["gate"]["pass"]:
        verdict = "TYPED_TENSOR_REALIZATION_IDENTIFIED"
    elif full["metrics"]["gate"]["pass"]:
        verdict = "FULL_TENSOR_IDENTIFIED_CP_RANK6_TRAINING_FAIL"
    elif cp["metrics"]["gate"]["pass"]:
        verdict = "CP_RANK6_IDENTIFIED_FULL_TENSOR_TRAINING_FAIL"
    else:
        verdict = "TYPED_TENSOR_LEARNED_HEADS_FAIL"

    output.mkdir(parents=True, exist_ok=True)
    for name, trace in (("full_240_trace.jsonl", full_trace),
                        ("cp_rank_6_trace.jsonl", cp_trace)):
        (output / name).write_text("".join(json.dumps(row, sort_keys=True) + "\n"
                                           for row in trace), encoding="utf-8")
    for name, arm in (("full_240_model.pt", full), ("cp_rank_6_model.pt", cp)):
        torch.save({"model_state": arm["model"].state_dict(),
                    "metrics": arm["metrics"]}, output / name)
    result = {
        "schema": "MetaSieve.E0R0TypedTensor.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P1R2B-E0R0_TYPED_TENSOR_IDENTIFIABILITY",
        "verdict": verdict,
        "selection": gate_report["selection"],
        "optimization": gate_report["optimization"],
        "analytic_witness": {
            "maximum_residual_error": analytic_error,
            "maximum_raw_teacher_reconstruction_error": raw_reconstruction_error,
            "cp_rank_6_tensor_error": cp_error,
            "cp_rank_6_residual_error": cp_numeric_error,
            "ligand_ci": ligand_ci,
            "correct_ci": analytic_correct_ci,
            "deranged_ci": analytic_deranged_ci,
            "correct_minus_ligand": analytic_correct_ci - ligand_ci,
            "correct_minus_deranged": analytic_correct_ci - analytic_deranged_ci,
            "gate": analytic_gate,
        },
        "learned_full_240": full["metrics"],
        "learned_cp_rank_6": cp["metrics"],
        "frozen_map": gate_report["metrics"],
        "parameter_contract": {
            "full_expected": FULL_TENSOR_PARAMETERS,
            "cp_rank_6_expected": CP_TENSOR_PARAMETERS,
        },
        "inputs": {
            "input_manifest_sha256": sha256_file(input_root / "manifest.json"),
            "cache_manifest_sha256": sha256_file(cache_root / "manifest.json"),
            "checkpoint_sha256": sha256_file(checkpoint),
            "synthetic_gate_sha256": sha256_file(synthetic_root / "synthetic_gate.json"),
        },
        "affinity_labels_read": False,
        "recipient_labels_read": False,
        "davis_accessed": False,
        "typed_interaction_training_performed": False,
        "production_model_modified": False,
        "downstream_authorized": False,
    }
    result["model_hashes"] = {
        "full_240": sha256_file(output / "full_240_model.pt"),
        "cp_rank_6": sha256_file(output / "cp_rank_6_model.pt"),
    }
    (output / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result


def _write_report(output: Path, result: dict) -> None:
    analytic = result["analytic_witness"]
    full = result["learned_full_240"]
    cp = result["learned_cp_rank_6"]
    frozen = result["frozen_map"]
    text = f"""# E0R0 Typed Tensor Identifiability

Decision: `{result['verdict']}`.

| Arm | Parameters | Correct CI | Deranged CI | Correct-Ligand | Correct-Deranged | Gate |
|---|---:|---:|---:|---:|---:|---|
| Analytic tensor | 240 | {analytic['correct_ci']:.5f} | {analytic['deranged_ci']:.5f} | {analytic['correct_minus_ligand']:+.5f} | {analytic['correct_minus_deranged']:+.5f} | {analytic['gate']['pass']} |
| Learned full tensor | {full['parameter_count']} | {full['correct_ci']:.5f} | {full['deranged_ci']:.5f} | {full['correct_minus_ligand']:+.5f} | {full['correct_minus_deranged']:+.5f} | {full['gate']['pass']} |
| Learned CP rank 6 | {cp['parameter_count']} | {cp['correct_ci']:.5f} | {cp['deranged_ci']:.5f} | {cp['correct_minus_ligand']:+.5f} | {cp['correct_minus_deranged']:+.5f} | {cp['gate']['pass']} |
| Frozen generic MAP | comparison | {frozen['correct_ci']:.5f} | {frozen['deranged_ci']:.5f} | {frozen['correct_minus_ligand']:+.5f} | {frozen['correct_minus_deranged']:+.5f} | False |

The analytic 240D residual error is `{analytic['maximum_residual_error']:.3g}`.
The numerical CP-rank-6 tensor/residual errors are
`{analytic['cp_rank_6_tensor_error']:.3g}` and
`{analytic['cp_rank_6_residual_error']:.3g}`.

No affinity labels or DAVIS data were read. PLIP/typed-interaction training,
production integration, CSMO/Band changes and downstream authorization remain
outside this stage.
"""
    (output / "STAGE_REPORT.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path,
                        default=Path("dataset/processed/source_affinity/e0_input_v1"))
    parser.add_argument("--cache", type=Path,
                        default=Path("dataset/processed/source_affinity/e0_local_states_v1"))
    parser.add_argument("--checkpoint", type=Path,
                        default=Path("report/mechanism_refactor/p1b_pilot20k_seed17_v1/best.pt"))
    parser.add_argument("--synthetic", type=Path,
                        default=Path("research/e0_identifiability/artifacts/e0_local_map_v1"))
    parser.add_argument("--output", type=Path,
                        default=Path("research/e0_identifiability/artifacts/e0r0_typed_tensor_v1"))
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    result = run(args.input, args.cache, args.checkpoint, args.synthetic,
                 args.output, args.device)
    _write_report(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
