"""Source-only gate for a protein-conditioned two-channel adaptation state.

The gate opens only source ``fit`` and ``probe`` roles. It trains interaction
coordinates on component-held-out fit tasks and evaluates them on scaffold-
disjoint probe tasks. ``locked`` and recipient labels are never loaded.

This is an admission experiment, not promoted model code. The empirical-Bayes
state solve is an instrument for testing representation recoverability; a final
meta-adaptation claim additionally requires a learned support operator to beat
this solve in a later registered gate.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
import time
from typing import Iterable

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch.nn import functional as F

from research.a2s.a2s_information_gate import sha256_file
from research.a2s.a2s_mode_generalization import TargetSplit, target_splits
from research.a2s.a2s_trace import DEVICE, Substrate, load_substrate, tanimoto
from research.a2s.a2s_trace_stratum import DEFAULT_LOCK, DEFAULT_OOF, metric_loss, paired_bootstrap


ROOT = Path(__file__).resolve().parents[2]
SEGMENT_ARCHIVE = (
    ROOT / "dataset" / "public" / "chembl_37" / "processed" / "dualcold"
    / "target_esm2_segments32.npz"
)
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_interaction_state_gate_2026-08-02.json"
DEFAULT_RECORDS = (
    ROOT / "reports" / "active" / "a2s_interaction_state_gate_records_2026-08-02.parquet"
)
DEFAULT_WEIGHTS = (
    ROOT / "reports" / "active" / "a2s_interaction_state_gate_weights_2026-08-02.pt"
)

SEEDS = (1729, 1731, 1733)
MODES = ("segment", "pooled", "ligand", "random")
SUPPORT_K = (1, 3, 5)
TRAIN_K = (3, 5)
SEGMENT_DIM = 32
HIDDEN = 32
STATE_DIM = 2
STATE_NORM = math.sqrt(2.0)
RIDGE = 1.0
ORACLE_RIDGE = 0.03
DELTA_BOUND = 1.5
EPOCHS = 30
PATIENCE = 7
LEARNING_RATE = 3e-3
WEIGHT_DECAY = 1e-4
MICRO_BATCH = 12
MAX_TRAIN_QUERY = 48
DRAWS = 8
MIN_QUERY = 5
LOW_SIMILARITY = 0.35
BOOTSTRAP_DRAWS = 2000
MDE = 0.005
PCA_SEED = 20260802


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def stable_seed(*parts: object) -> int:
    return int(sha256(canonical(parts).encode()).hexdigest()[:8], 16)


def assert_source_roles(frame: pd.DataFrame) -> None:
    roles = set(frame.role.astype(str).unique())
    if roles - {"fit", "probe"}:
        raise AssertionError(f"unauthorized source role entered gate: {sorted(roles)}")
    if roles != {"fit", "probe"}:
        raise AssertionError(f"fit/probe source roles are both required, got {sorted(roles)}")


@dataclass(frozen=True)
class SegmentFeatures:
    by_target: dict[str, torch.Tensor]
    projection_shape: tuple[int, int]
    fit_segments: int


def build_segment_features(substrate: Substrate) -> SegmentFeatures:
    """Fit a label-free ESM segment projection on fit targets only."""

    archive = np.load(SEGMENT_ARCHIVE, allow_pickle=False)
    keys = [str(value) for value in archive["keys"]]
    values = np.asarray(archive["segments"], dtype=np.float32)
    key_to_row = {key: row for row, key in enumerate(keys)}
    targets = sorted(set(substrate.labeled.target.astype(str)))
    missing = sorted(set(targets) - set(key_to_row))
    if missing:
        raise RuntimeError(f"segment archive misses {len(missing)} source targets")

    fit_targets = sorted(
        set(substrate.labeled.loc[substrate.labeled.role == "fit", "target"].astype(str))
    )
    fit_rows = np.asarray([key_to_row[target] for target in fit_targets], dtype=np.int64)
    fit = torch.as_tensor(values[fit_rows].reshape(-1, values.shape[-1]), device=DEVICE)
    centre = fit.mean(dim=0, keepdim=True)
    torch.manual_seed(PCA_SEED)
    _, _, right = torch.pca_lowrank(fit - centre, q=SEGMENT_DIM, center=False, niter=4)

    selected = torch.as_tensor(
        values[[key_to_row[target] for target in targets]], device=DEVICE
    )
    scores = (selected - centre.reshape(1, 1, -1)) @ right[:, :SEGMENT_DIM]
    fit_scores = (fit - centre) @ right[:, :SEGMENT_DIM]
    mean = fit_scores.mean(dim=0, keepdim=True)
    scale = fit_scores.std(dim=0, keepdim=True).clamp(min=1e-6)
    scores = (scores - mean.reshape(1, 1, -1)) / scale.reshape(1, 1, -1)
    return SegmentFeatures(
        by_target={target: scores[row] for row, target in enumerate(targets)},
        projection_shape=(int(values.shape[-1]), SEGMENT_DIM),
        fit_segments=int(len(fit)),
    )


class InteractionChannels(nn.Module):
    """Two bounded coordinates from ligand motifs and protein sequence segments."""

    def __init__(self, ligand_dim: int, mode: str, hidden: int = HIDDEN) -> None:
        super().__init__()
        if mode not in MODES and mode != "random":
            raise ValueError(f"unknown interaction mode: {mode}")
        self.mode = mode
        self.ligand = nn.Sequential(
            nn.Linear(ligand_dim, hidden),
            nn.SiLU(),
            nn.LayerNorm(hidden),
        )
        self.output = nn.Linear(hidden, STATE_DIM, bias=False)
        if mode in {"segment", "pooled", "random"}:
            self.segment = nn.Linear(SEGMENT_DIM, hidden, bias=False)
            self.query = nn.Linear(hidden, hidden, bias=False)
            self.key = nn.Linear(hidden, hidden, bias=False)
            self.value = nn.Linear(hidden, hidden, bias=False)

    def forward(self, ligand: torch.Tensor, segments: torch.Tensor) -> torch.Tensor:
        ligand_hidden = self.ligand(ligand)
        if self.mode == "ligand":
            interaction = ligand_hidden
        else:
            segment_hidden = torch.tanh(self.segment(segments))
            if self.mode == "pooled":
                context = self.value(segment_hidden.mean(dim=0, keepdim=True)).expand_as(
                    ligand_hidden
                )
            else:
                logits = self.query(ligand_hidden) @ self.key(segment_hidden).T
                logits = logits / math.sqrt(ligand_hidden.shape[-1])
                context = torch.softmax(logits, dim=-1) @ self.value(segment_hidden)
            interaction = ligand_hidden * context
        return torch.tanh(self.output(interaction))


def active_dimensions(k: int) -> int:
    if k == 1:
        return 0
    if k == 3:
        return 1
    if k == 5:
        return 2
    raise ValueError(f"unsupported budget: {k}")


def solve_state(
    support_phi: torch.Tensor,
    support_residual: torch.Tensor,
    k: int,
    ridge: float = RIDGE,
) -> tuple[torch.Tensor, torch.Tensor, float]:
    """Centered ridge state with an orthogonally invariant norm bound."""

    dimensions = active_dimensions(k)
    state = torch.zeros(STATE_DIM, dtype=support_phi.dtype, device=support_phi.device)
    centre = support_phi.mean(dim=0)
    if dimensions == 0:
        return state, centre, 0.0
    design = support_phi[:, :dimensions] - centre[:dimensions]
    target = support_residual - support_residual.mean()
    gram = design.T @ design + ridge * torch.eye(
        dimensions, dtype=design.dtype, device=design.device
    )
    raw = torch.linalg.solve(gram, design.T @ target)
    norm = torch.linalg.vector_norm(raw).clamp(min=1e-12)
    raw = raw * torch.clamp(torch.as_tensor(STATE_NORM, device=raw.device) / norm, max=1.0)
    state[:dimensions] = raw
    smoother = design @ torch.linalg.solve(gram, design.T)
    effective_dof = float(torch.trace(smoother).detach())
    return state, centre, effective_dof


def rank_delta(query_phi: torch.Tensor, state: torch.Tensor, centre: torch.Tensor) -> torch.Tensor:
    raw = (query_phi - centre) @ state
    return DELTA_BOUND * torch.tanh(raw / DELTA_BOUND)


def adapted_prediction(
    support_phi: torch.Tensor,
    query_phi: torch.Tensor,
    support_residual: torch.Tensor,
    query_base: torch.Tensor,
    k: int,
    ridge: float = RIDGE,
) -> tuple[torch.Tensor, torch.Tensor, float]:
    state, centre, dof = solve_state(support_phi, support_residual, k, ridge)
    return query_base + rank_delta(query_phi, state, centre), state, dof


def pairwise_proper_loss(label: torch.Tensor, prediction: torch.Tensor) -> torch.Tensor:
    left, right = torch.triu_indices(len(label), len(label), offset=1, device=label.device)
    difference = label[left] - label[right]
    active = difference != 0
    if not bool(active.any()):
        return prediction.sum() * 0.0
    logits = prediction[left[active]] - prediction[right[active]]
    target = (difference[active] > 0).to(logits.dtype)
    return F.binary_cross_entropy_with_logits(logits, target)


def pairwise_proper_value(label: np.ndarray, prediction: np.ndarray) -> float:
    label_t = torch.as_tensor(label, dtype=torch.float64)
    prediction_t = torch.as_tensor(prediction, dtype=torch.float64)
    return float(pairwise_proper_loss(label_t, prediction_t))


def norm_matched_transplant(correct: torch.Tensor, wrong: torch.Tensor) -> torch.Tensor:
    correct_c = correct - correct.mean()
    wrong_c = wrong - wrong.mean()
    wrong_norm = torch.linalg.vector_norm(wrong_c)
    if float(wrong_norm) < 1e-12:
        return torch.zeros_like(correct_c)
    return wrong_c * (torch.linalg.vector_norm(correct_c) / wrong_norm)


def ligand_matrix(substrate: Substrate) -> torch.Tensor:
    return torch.cat((substrate.bits, substrate.desc), dim=1)


def scaffold_splits(substrate: Substrate, role: str) -> list[TargetSplit]:
    return [split for split in target_splits(substrate, role) if split.split == "scaffold_disjoint"]


def split_fit_components(splits: list[TargetSplit], seed: int) -> tuple[list[TargetSplit], list[TargetSplit]]:
    components = sorted(set(split.component for split in splits))
    rng = np.random.default_rng(stable_seed("fit-component-split", seed))
    validation = set(rng.choice(components, size=max(1, len(components) // 5), replace=False))
    return (
        [split for split in splits if split.component not in validation],
        [split for split in splits if split.component in validation],
    )


def select_rows(rows: np.ndarray, count: int, rng: np.random.Generator) -> np.ndarray:
    if len(rows) <= count:
        return rows
    return rng.choice(rows, size=count, replace=False)


def episode_loss(
    model: InteractionChannels,
    split: TargetSplit,
    k: int,
    rng: np.random.Generator,
    ligand: torch.Tensor,
    segments: SegmentFeatures,
    substrate: Substrate,
) -> torch.Tensor:
    support = rng.choice(split.train_rows, size=k, replace=False)
    query = select_rows(split.eval_rows, MAX_TRAIN_QUERY, rng)
    rows = np.concatenate((support, query))
    phi = model(ligand[torch.as_tensor(rows, device=DEVICE)], segments.by_target[split.target])
    support_phi, query_phi = phi[:k], phi[k:]
    residual = substrate.affinity[torch.as_tensor(support, device=DEVICE)] - substrate.base[
        torch.as_tensor(support, device=DEVICE)
    ]
    prediction, _, _ = adapted_prediction(
        support_phi,
        query_phi,
        residual,
        substrate.base[torch.as_tensor(query, device=DEVICE)],
        k,
    )
    label = substrate.affinity[torch.as_tensor(query, device=DEVICE)]
    loss = pairwise_proper_loss(label, prediction)
    # Prevent a zero or duplicated channel solution without dictating chemistry.
    channel_sd = phi.std(dim=0, unbiased=False)
    scale_penalty = ((channel_sd - 0.5) ** 2).mean()
    centred = phi - phi.mean(dim=0, keepdim=True)
    covariance = (centred[:, 0] * centred[:, 1]).mean()
    return loss + 0.01 * scale_penalty + 0.005 * covariance.square()


@torch.no_grad()
def validation_loss(
    model: InteractionChannels,
    splits: list[TargetSplit],
    ligand: torch.Tensor,
    segments: SegmentFeatures,
    substrate: Substrate,
    seed: int,
) -> float:
    model.eval()
    losses: list[float] = []
    for split in splits:
        for k in TRAIN_K:
            rng = np.random.default_rng(stable_seed("validation", seed, split.target, k))
            losses.append(float(episode_loss(model, split, k, rng, ligand, segments, substrate)))
    return float(np.mean(losses)) if losses else float("inf")


def train_representation(
    mode: str,
    seed: int,
    train_splits: list[TargetSplit],
    validation_splits: list[TargetSplit],
    ligand: torch.Tensor,
    segments: SegmentFeatures,
    substrate: Substrate,
    epochs: int = EPOCHS,
) -> tuple[InteractionChannels, dict[str, object]]:
    torch.manual_seed(seed)
    model = InteractionChannels(ligand.shape[1], mode).to(DEVICE)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY
    )
    best_state = deepcopy(model.state_dict())
    best_loss = float("inf")
    best_epoch = -1
    stale = 0
    history: list[dict[str, float]] = []

    for epoch in range(epochs):
        model.train()
        rng = np.random.default_rng(stable_seed("train", seed, mode, epoch))
        order = rng.permutation(len(train_splits))
        optimizer.zero_grad(set_to_none=True)
        pending: list[torch.Tensor] = []
        train_values: list[float] = []
        for position in order:
            split = train_splits[int(position)]
            for k in TRAIN_K:
                loss = episode_loss(model, split, k, rng, ligand, segments, substrate)
                pending.append(loss)
                train_values.append(float(loss.detach()))
                if len(pending) == MICRO_BATCH:
                    torch.stack(pending).mean().backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
                    optimizer.step()
                    optimizer.zero_grad(set_to_none=True)
                    pending.clear()
        if pending:
            torch.stack(pending).mean().backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            optimizer.zero_grad(set_to_none=True)

        value = validation_loss(
            model, validation_splits, ligand, segments, substrate, seed
        )
        history.append(
            {"epoch": float(epoch), "train": float(np.mean(train_values)), "validation": value}
        )
        if value < best_loss - 1e-4:
            best_loss = value
            best_epoch = epoch
            best_state = deepcopy(model.state_dict())
            stale = 0
        else:
            stale += 1
        if stale >= PATIENCE:
            break

    model.load_state_dict(best_state)
    model.eval()
    return model, {
        "mode": mode,
        "seed": seed,
        "best_epoch": int(best_epoch),
        "best_validation_loss": float(best_loss),
        "epochs_run": int(len(history)),
        "history": history,
    }


def fixed_rotation(seed: int) -> torch.Tensor:
    rng = np.random.default_rng(stable_seed("rotation", seed))
    matrix, _ = np.linalg.qr(rng.normal(size=(STATE_DIM, STATE_DIM)))
    return torch.as_tensor(matrix.astype(np.float32), device=DEVICE)


def prediction_metrics(label: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    result = metric_loss(label, prediction)
    result["proper"] = pairwise_proper_value(label, prediction)
    return result


def add_metrics(row: dict[str, object], name: str, label: np.ndarray, prediction: np.ndarray) -> None:
    for metric, value in prediction_metrics(label, prediction).items():
        row[f"{name}__{metric}"] = float(value)


@torch.no_grad()
def evaluate_model(
    model: InteractionChannels,
    mode: str,
    seed: int,
    splits: list[TargetSplit],
    ligand: torch.Tensor,
    segments: SegmentFeatures,
    substrate: Substrate,
) -> pd.DataFrame:
    model.eval()
    residual_all = substrate.affinity - substrate.base
    records: list[dict[str, object]] = []
    rotation = fixed_rotation(seed)

    for split_index, split in enumerate(splits):
        other = splits[(split_index + 7) % len(splits)]
        query_rows = torch.as_tensor(split.eval_rows, device=DEVICE)
        query_phi = model(ligand[query_rows], segments.by_target[split.target])
        zero_segments = torch.zeros_like(segments.by_target[split.target])
        query_phi_zero_protein = model(ligand[query_rows], zero_segments)
        query_phi_wrong_protein = model(ligand[query_rows], segments.by_target[other.target])
        query_base = substrate.base[query_rows]
        query_label = substrate.affinity[query_rows]

        train_rows = torch.as_tensor(split.train_rows, device=DEVICE)
        train_phi = model(ligand[train_rows], segments.by_target[split.target])
        oracle_prediction, _, _ = adapted_prediction(
            train_phi,
            query_phi,
            residual_all[train_rows],
            query_base,
            5,
            ORACLE_RIDGE,
        )

        for draw in range(DRAWS):
            current_rng = np.random.default_rng(stable_seed("eval", seed, split.target, draw))
            other_rng = np.random.default_rng(stable_seed("eval", seed, other.target, draw))
            current_order = current_rng.permutation(split.train_rows)
            other_order = other_rng.permutation(other.train_rows)
            for k in SUPPORT_K:
                support = current_order[:k]
                other_support = other_order[:k]
                support_rows = torch.as_tensor(support, device=DEVICE)
                other_rows = torch.as_tensor(other_support, device=DEVICE)
                support_phi = model(ligand[support_rows], segments.by_target[split.target])
                support_phi_zero_protein = model(ligand[support_rows], zero_segments)
                support_phi_wrong_protein = model(
                    ligand[support_rows], segments.by_target[other.target]
                )
                residual = residual_all[support_rows]
                wrong_residual = norm_matched_transplant(residual, residual_all[other_rows])

                correct, state, dof = adapted_prediction(
                    support_phi, query_phi, residual, query_base, k
                )
                wrong, _, _ = adapted_prediction(
                    support_phi, query_phi, wrong_residual, query_base, k
                )
                if k == 1:
                    deranged_residual = residual
                else:
                    deranged_residual = torch.roll(residual, shifts=1)
                deranged, _, _ = adapted_prediction(
                    support_phi, query_phi, deranged_residual, query_base, k
                )
                protein_zero, _, _ = adapted_prediction(
                    support_phi_zero_protein,
                    query_phi_zero_protein,
                    residual,
                    query_base,
                    k,
                )
                protein_transplant, _, _ = adapted_prediction(
                    support_phi_wrong_protein,
                    query_phi_wrong_protein,
                    residual,
                    query_base,
                    k,
                )
                rotated, _, _ = adapted_prediction(
                    support_phi @ rotation,
                    query_phi @ rotation,
                    residual,
                    query_base,
                    k,
                )

                nearest = tanimoto(
                    substrate.bits[query_rows].unsqueeze(0),
                    substrate.bits[support_rows].unsqueeze(0),
                ).amax(dim=-1).squeeze(0)
                masks = {
                    "all": torch.ones_like(nearest, dtype=torch.bool),
                    "low": nearest < LOW_SIMILARITY,
                    "local": nearest >= LOW_SIMILARITY,
                }
                predictions = {
                    "base": query_base,
                    "oracle": oracle_prediction,
                    "correct": correct,
                    "deranged": deranged,
                    "wrong_residual": wrong,
                    "protein_zero": protein_zero,
                    "protein_transplant": protein_transplant,
                    "rotated": rotated,
                }
                for stratum, mask in masks.items():
                    if int(mask.sum()) < MIN_QUERY:
                        continue
                    label_np = query_label[mask].cpu().numpy()
                    if float(np.std(label_np)) < 1e-9:
                        continue
                    row: dict[str, object] = {
                        "mode": mode,
                        "model_seed": seed,
                        "target": split.target,
                        "component": split.component,
                        "draw": draw,
                        "k": k,
                        "stratum": stratum,
                        "n_query": int(mask.sum()),
                        "nearest_tanimoto_mean": float(nearest[mask].mean()),
                        "state_0": float(state[0]),
                        "state_1": float(state[1]),
                        "state_norm": float(torch.linalg.vector_norm(state)),
                        "effective_dof": dof,
                    }
                    for name, prediction in predictions.items():
                        add_metrics(row, name, label_np, prediction[mask].cpu().numpy())
                    records.append(row)
    return pd.DataFrame.from_records(records)


METRICS = ("ci", "ndcg10", "spearman", "proper", "rmse")
ARMS = (
    "base",
    "oracle",
    "correct",
    "deranged",
    "wrong_residual",
    "protein_zero",
    "protein_transplant",
    "rotated",
)


def contrast_values(frame: pd.DataFrame, left: str, right: str, metric: str) -> pd.Series:
    if metric in {"proper", "rmse"}:
        return frame[f"{right}__{metric}"] - frame[f"{left}__{metric}"]
    return frame[f"{left}__{metric}"] - frame[f"{right}__{metric}"]


def summarise_records(records: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    contrasts = {
        "oracle_minus_base": ("oracle", "base"),
        "correct_minus_base": ("correct", "base"),
        "correct_minus_deranged": ("correct", "deranged"),
        "correct_minus_wrong_residual": ("correct", "wrong_residual"),
        "correct_minus_protein_zero": ("correct", "protein_zero"),
        "correct_minus_protein_transplant": ("correct", "protein_transplant"),
        "rotated_minus_correct": ("rotated", "correct"),
    }
    for mode in sorted(records["mode"].unique()):
        mode_frame = records.loc[records["mode"] == mode]
        mode_summary: dict[str, object] = {}
        for k in SUPPORT_K:
            k_frame = mode_frame.loc[mode_frame.k == k]
            k_summary: dict[str, object] = {}
            for stratum in sorted(k_frame.stratum.unique()):
                cell = k_frame.loc[k_frame.stratum == stratum].copy()
                entry: dict[str, object] = {
                    "components": int(cell.component.nunique()),
                    "targets": int(cell.target.nunique()),
                    "records": int(len(cell)),
                    "mean_query": float(cell.n_query.mean()),
                    "mean_dof": float(cell.effective_dof.mean()),
                    "mean_state_norm": float(cell.state_norm.mean()),
                    "absolute": {},
                    "contrasts": {},
                }
                for arm in ARMS:
                    entry["absolute"][arm] = {
                        metric: float(cell[f"{arm}__{metric}"].mean()) for metric in METRICS
                    }
                for name, (left, right) in contrasts.items():
                    entry["contrasts"][name] = {}
                    for metric in METRICS:
                        cell["_delta"] = contrast_values(cell, left, right, metric)
                        entry["contrasts"][name][metric] = paired_bootstrap(
                            cell, "_delta", seed=stable_seed(mode, k, stratum, name, metric),
                            draws=BOOTSTRAP_DRAWS,
                        )
                k_summary[stratum] = entry
            mode_summary[f"k{k}"] = k_summary
        summary[mode] = mode_summary
    return summary


def compare_modes(records: pd.DataFrame, left: str, right: str) -> dict[str, object]:
    keys = ["model_seed", "target", "component", "draw", "k", "stratum"]
    columns = keys + [f"correct__{metric}" for metric in METRICS]
    merged = records.loc[records["mode"] == left, columns].merge(
        records.loc[records["mode"] == right, columns],
        on=keys,
        suffixes=("__left", "__right"),
        validate="one_to_one",
    )
    result: dict[str, object] = {}
    for k in SUPPORT_K:
        result[f"k{k}"] = {}
        for stratum in sorted(merged.stratum.unique()):
            cell = merged.loc[(merged.k == k) & (merged.stratum == stratum)].copy()
            if cell.empty:
                continue
            result[f"k{k}"][stratum] = {}
            for metric in METRICS:
                left_column = f"correct__{metric}__left"
                right_column = f"correct__{metric}__right"
                if metric in {"proper", "rmse"}:
                    cell["_delta"] = cell[right_column] - cell[left_column]
                else:
                    cell["_delta"] = cell[left_column] - cell[right_column]
                result[f"k{k}"][stratum][metric] = paired_bootstrap(
                    cell,
                    "_delta",
                    seed=stable_seed("mode", left, right, k, stratum, metric),
                    draws=BOOTSTRAP_DRAWS,
                )
    return result


def synthetic_control(seed: int = 99173) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    records: list[dict[str, object]] = []
    targets = 80
    observations = 36
    noise = 0.8
    for target in range(targets):
        phi = rng.normal(size=(observations, STATE_DIM)).astype(np.float32)
        state = rng.normal(scale=1.6, size=STATE_DIM).astype(np.float32)
        label = phi @ state + rng.normal(scale=noise, size=observations)
        for draw in range(8):
            order = rng.permutation(observations)
            query = order[5:]
            for k in SUPPORT_K:
                support = order[:k]
                other_state = rng.normal(scale=1.6, size=STATE_DIM).astype(np.float32)
                wrong_label = phi[support] @ other_state + rng.normal(scale=noise, size=k)
                phi_s = torch.as_tensor(phi[support])
                phi_q = torch.as_tensor(phi[query])
                base = torch.zeros(len(query))
                correct, _, _ = adapted_prediction(
                    phi_s, phi_q, torch.as_tensor(label[support], dtype=torch.float32), base, k
                )
                wrong, _, _ = adapted_prediction(
                    phi_s, phi_q, torch.as_tensor(wrong_label, dtype=torch.float32), base, k
                )
                truth = label[query]
                records.append(
                    {
                        "component": f"c{target:03d}",
                        "target": f"t{target:03d}",
                        "draw": draw,
                        "k": k,
                        "correct": metric_loss(truth, correct.numpy())["ci"],
                        "wrong": metric_loss(truth, wrong.numpy())["ci"],
                        "base": metric_loss(truth, np.zeros_like(truth))["ci"],
                        "max_abs_k1": float(torch.max(torch.abs(correct))) if k == 1 else 0.0,
                    }
                )
    frame = pd.DataFrame.from_records(records)
    result: dict[str, object] = {}
    for k in SUPPORT_K:
        cell = frame.loc[frame.k == k].copy()
        cell["gain"] = cell.correct - cell.base
        cell["specificity"] = cell.correct - cell.wrong
        result[f"k{k}"] = {
            "gain": paired_bootstrap(cell, "gain", seed=seed + k, draws=BOOTSTRAP_DRAWS),
            "correct_minus_wrong": paired_bootstrap(
                cell, "specificity", seed=seed + 100 + k, draws=BOOTSTRAP_DRAWS
            ),
            "max_abs_k1": float(cell.max_abs_k1.max()),
        }
    passed = (
        result["k1"]["max_abs_k1"] == 0.0
        and result["k5"]["correct_minus_wrong"]["lower95"] > MDE
        and result["k5"]["gain"]["lower95"] > MDE
        and result["k5"]["gain"]["mean"] >= result["k3"]["gain"]["mean"]
    )
    return {"noise_sd": noise, "targets": targets, "result": result, "pass": bool(passed)}


def value_at(summary: dict[str, object], path: Iterable[str], default: float = float("nan")) -> float:
    current: object = summary
    for key in path:
        if not isinstance(current, dict) or key not in current:
            return default
        current = current[key]
    return float(current) if isinstance(current, (int, float)) else default


def decide(
    synthetic: dict[str, object],
    summary: dict[str, object],
    comparisons: dict[str, object],
) -> dict[str, object]:
    primary = summary.get("segment", {})
    k3 = primary.get("k3", {}).get("all", {})
    k5 = primary.get("k5", {}).get("all", {})
    low5 = primary.get("k5", {}).get("low", {})
    checks = {
        "synthetic_positive_control": bool(synthetic.get("pass", False)),
        "pooled_components_at_least_47": int(k5.get("components", 0)) >= 47,
        "oracle_pooled_lcb_above_mde": value_at(
            k5, ("contrasts", "oracle_minus_base", "ci", "lower95")
        ) > MDE,
        "oracle_low_similarity_lcb_above_mde": value_at(
            low5, ("contrasts", "oracle_minus_base", "ci", "lower95")
        ) > MDE,
        "k5_gain_lcb_above_mde": value_at(
            k5, ("contrasts", "correct_minus_base", "ci", "lower95")
        ) > MDE,
        "k3_to_k5_monotone": value_at(
            k5, ("contrasts", "correct_minus_base", "ci", "mean")
        ) >= value_at(k3, ("contrasts", "correct_minus_base", "ci", "mean")),
        "correct_beats_deranged": value_at(
            k5, ("contrasts", "correct_minus_deranged", "ci", "lower95")
        ) > 0.0,
        "correct_beats_wrong_residual": value_at(
            k5, ("contrasts", "correct_minus_wrong_residual", "ci", "lower95")
        ) > 0.0,
        "correct_beats_protein_zero": value_at(
            k5, ("contrasts", "correct_minus_protein_zero", "ci", "lower95")
        ) > 0.0,
        "correct_beats_protein_transplant": value_at(
            k5, ("contrasts", "correct_minus_protein_transplant", "ci", "lower95")
        ) > 0.0,
        "rotation_invariant": abs(
            value_at(k5, ("contrasts", "rotated_minus_correct", "ci", "mean"))
        ) < 1e-8,
        "segment_beats_ligand": value_at(
            comparisons.get("segment_minus_ligand", {}), ("k5", "all", "ci", "lower95")
        ) > 0.0,
        "segment_beats_pooled": value_at(
            comparisons.get("segment_minus_pooled", {}), ("k5", "all", "ci", "lower95")
        ) > 0.0,
        "segment_beats_frozen_random": value_at(
            comparisons.get("segment_minus_random", {}), ("k5", "all", "ci", "lower95")
        ) > 0.0,
    }
    passed = all(checks.values())
    return {
        "checks": checks,
        "pass": bool(passed),
        "verdict": (
            "INTERACTION_STATE_REPRESENTATION_ADMITTED_PROCEED_R1"
            if passed
            else "INTERACTION_STATE_REPRESENTATION_NOT_ADMITTED"
        ),
        "authorizes": (
            "R1 learned support-to-state operator only" if passed else "no learned operator or promotion"
        ),
    }


def run(
    output: Path,
    records_path: Path,
    weights_path: Path,
    *,
    seeds: tuple[int, ...] = SEEDS,
    modes: tuple[str, ...] = MODES,
    epochs: int = EPOCHS,
) -> dict[str, object]:
    started = time.time()
    substrate, substrate_metadata = load_substrate(DEFAULT_LOCK, DEFAULT_OOF)
    assert_source_roles(substrate.labeled)
    segment_features = build_segment_features(substrate)
    ligand = ligand_matrix(substrate)
    fit_splits = scaffold_splits(substrate, "fit")
    probe_splits = scaffold_splits(substrate, "probe")
    if len(set(split.component for split in probe_splits)) < 47:
        raise RuntimeError("probe gate does not meet the preregistered component floor")

    synthetic = synthetic_control()
    all_records: list[pd.DataFrame] = []
    training: list[dict[str, object]] = []
    weights: dict[str, dict[str, torch.Tensor]] = {}
    for seed in seeds:
        train_splits, validation_splits = split_fit_components(fit_splits, seed)
        for mode in modes:
            if mode == "random":
                torch.manual_seed(seed)
                model = InteractionChannels(ligand.shape[1], mode).to(DEVICE).eval()
                train_record = {
                    "mode": mode,
                    "seed": seed,
                    "best_epoch": None,
                    "best_validation_loss": None,
                    "epochs_run": 0,
                    "history": [],
                    "frozen_random_control": True,
                }
            else:
                model, train_record = train_representation(
                    mode,
                    seed,
                    train_splits,
                    validation_splits,
                    ligand,
                    segment_features,
                    substrate,
                    epochs,
                )
            training.append(
                {
                    **train_record,
                    "train_components": len(set(split.component for split in train_splits)),
                    "validation_components": len(
                        set(split.component for split in validation_splits)
                    ),
                }
            )
            weights[f"{mode}_{seed}"] = {
                name: value.detach().cpu() for name, value in model.state_dict().items()
            }
            all_records.append(
                evaluate_model(
                    model,
                    mode,
                    seed,
                    probe_splits,
                    ligand,
                    segment_features,
                    substrate,
                )
            )

    records = pd.concat(all_records, ignore_index=True)
    summary = summarise_records(records)
    comparisons: dict[str, object] = {}
    if {"segment", "ligand"}.issubset(set(modes)):
        comparisons["segment_minus_ligand"] = compare_modes(records, "segment", "ligand")
    if {"segment", "pooled"}.issubset(set(modes)):
        comparisons["segment_minus_pooled"] = compare_modes(records, "segment", "pooled")
    if {"segment", "random"}.issubset(set(modes)):
        comparisons["segment_minus_random"] = compare_modes(records, "segment", "random")
    decision = decide(synthetic, summary, comparisons)

    output.parent.mkdir(parents=True, exist_ok=True)
    records_path.parent.mkdir(parents=True, exist_ok=True)
    weights_path.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    torch.save(weights, weights_path)
    payload: dict[str, object] = {
        "schema": "a2s-interaction-state-gate-v1",
        "status": "SOURCE_ONLY",
        "protocol": {
            "roles_opened": ["fit", "probe"],
            "locked_labels_requested": False,
            "recipient_labels_requested": False,
            "support_k": list(SUPPORT_K),
            "state_dimensions": {"k1": 0, "k3": 1, "k5": 2},
            "seeds": list(seeds),
            "modes": list(modes),
            "epochs_max": epochs,
            "ridge": RIDGE,
            "oracle_ridge": ORACLE_RIDGE,
            "delta_bound": DELTA_BOUND,
            "low_similarity_threshold": LOW_SIMILARITY,
            "probe_role_used_once_by_this_gate": True,
        },
        "data": {
            "rows": int(len(substrate.labeled)),
            "fit_rows": int((substrate.labeled.role == "fit").sum()),
            "probe_rows": int((substrate.labeled.role == "probe").sum()),
            "fit_splits": len(fit_splits),
            "fit_components": len(set(split.component for split in fit_splits)),
            "probe_splits": len(probe_splits),
            "probe_components": len(set(split.component for split in probe_splits)),
            "segment_projection_shape": list(segment_features.projection_shape),
            "fit_segments": segment_features.fit_segments,
        },
        "synthetic": synthetic,
        "training": training,
        "summary": summary,
        "mode_comparisons": comparisons,
        "decision": decision,
        "substrate": substrate_metadata,
        "runtime": {
            "device": str(DEVICE),
            "cuda_name": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
            "wall_clock_seconds": round(time.time() - started, 3),
        },
        "artifacts": {
            "lock_sha256": sha256_file(DEFAULT_LOCK),
            "oof_sha256": sha256_file(DEFAULT_OOF),
            "segments_sha256": sha256_file(SEGMENT_ARCHIVE),
            "records": str(records_path.relative_to(ROOT)),
            "records_sha256": sha256_file(records_path),
            "weights": str(weights_path.relative_to(ROOT)),
            "weights_sha256": sha256_file(weights_path),
        },
    }
    payload["content_sha256"] = sha256(canonical(payload).encode()).hexdigest()
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def parse_values(value: str, cast: type) -> tuple:
    return tuple(cast(item.strip()) for item in value.split(",") if item.strip())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    parser.add_argument("--weights", type=Path, default=DEFAULT_WEIGHTS)
    parser.add_argument("--seeds", default=",".join(str(seed) for seed in SEEDS))
    parser.add_argument("--modes", default=",".join(MODES))
    parser.add_argument("--epochs", type=int, default=EPOCHS)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    if args.quick:
        parser.error(
            "--quick is disabled because the shared substrate loader opens probe labels; "
            "use the unit tests before the single full registered gate"
        )
    seeds = parse_values(args.seeds, int)
    modes = parse_values(args.modes, str)
    payload = run(
        args.output,
        args.records,
        args.weights,
        seeds=seeds,
        modes=modes,
        epochs=args.epochs,
    )
    print(payload["decision"]["verdict"])
    print(canonical(payload["decision"]["checks"]))


if __name__ == "__main__":
    main()
