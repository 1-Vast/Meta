"""F1B section-conditioned bilinear mechanism operator.

This module is an isolated research pilot.  It trains only on PKIS1, uses
PKIS2 and Anastassiadis2011 as already-consumed development panels, and never
imports or modifies ``model/``.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import torch
from torch import nn

from ..pkis_mechanism_pilot.mechanism import (
    double_center,
    protein_model_features,
    protein_pair_properties,
    stable_fold,
)
from .ceiling_probe import BOOTSTRAP_DRAWS, SEEDS, _load


SEED = 20260808
K_PRIMARY = 5
LIGAND_PCA_DIM = 24
PROTEIN_PCA_DIM = 12
BILINEAR_WIDTH = 16
PAIR_WIDTH = 32
SET_WIDTH = 32
EPISODE_BATCH = 24
QUERY_SIZE = 24
MAX_STEPS = 1500
VALIDATE_EVERY = 50
PATIENCE = 8
LOCATION_RIDGE = 10.0
ADDITIVE_RIDGE = 1000.0

VIEW_NAMES = ("hinge_polar", "dfg_back", "front_solvent")
LIGAND_CHANNELS = ((0, 2), (2, 3, 4), (0, 1, 4))
PROTEIN_PROPERTIES = ((0, 1, 4), (4, 5, 6), (0, 1, 2, 3, 6))
SUBPOCKET_NAMES = ("hinge_region", "dfg_region", "front_pocket")


def _write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True, default=str) + "\n",
        encoding="utf-8",
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _subpocket_weights(path: str | Path, scale: float = 6.0) -> dict[str, np.ndarray]:
    """Build fixed soft residue masks from KiSSim's public distance ranges."""
    frame = pd.read_csv(path).set_index(["subpocket", "min_max"])
    out = {}
    for name in SUBPOCKET_NAMES:
        low = frame.loc[(name, "min")].to_numpy(dtype=np.float64)
        high = frame.loc[(name, "max")].to_numpy(dtype=np.float64)
        midpoint = 0.5 * (low + high)
        weight = np.exp(-midpoint / float(scale))
        out[name] = (weight / max(float(weight.max()), 1e-12)).astype(np.float32)
    return out


def _ligand_raw_views(records) -> list[np.ndarray]:
    descriptors = np.stack([record.nuisance_features[-12:] for record in records])
    fingerprints = np.stack([record.channel_fingerprints for record in records])
    views = []
    for channels in LIGAND_CHANNELS:
        selected = fingerprints[:, channels].reshape(len(records), -1)
        views.append(np.concatenate([selected, descriptors], axis=1).astype(np.float32))
    return views


def _protein_raw_views(target_items, weights: dict[str, np.ndarray]) -> list[np.ndarray]:
    rows = []
    for view, properties in zip(SUBPOCKET_NAMES, PROTEIN_PROPERTIES):
        encoded = []
        weight = weights[view]
        for item in target_items:
            values = protein_pair_properties(item["record"].pocket)
            encoded.append((values[list(properties)] * weight[None, :]).reshape(-1))
        rows.append(np.asarray(encoded, dtype=np.float32))
    return rows


@dataclass
class ViewTransforms:
    ligand: list[object]
    protein: list[object]

    def ligand_apply(self, raw: list[np.ndarray]) -> list[np.ndarray]:
        return [pipeline.transform(value).astype(np.float32)
                for pipeline, value in zip(self.ligand, raw)]

    def protein_apply(self, raw: list[np.ndarray]) -> list[np.ndarray]:
        return [pipeline.transform(value).astype(np.float32)
                for pipeline, value in zip(self.protein, raw)]


def _pca_pipeline(value: np.ndarray, dimension: int, seed: int):
    n_components = min(int(dimension), value.shape[0] - 1, value.shape[1])
    pipeline = make_pipeline(
        PCA(n_components=n_components, svd_solver="randomized", random_state=seed),
        StandardScaler(),
    )
    return pipeline.fit(value)


def _fit_view_transforms(ligand_raw: list[np.ndarray], protein_raw: list[np.ndarray],
                         ligand_indices=None, protein_indices=None) -> ViewTransforms:
    ligand_indices = (np.arange(len(ligand_raw[0])) if ligand_indices is None
                      else np.asarray(ligand_indices, dtype=np.int64))
    protein_indices = (np.arange(len(protein_raw[0])) if protein_indices is None
                       else np.asarray(protein_indices, dtype=np.int64))
    ligand = [_pca_pipeline(value[ligand_indices], LIGAND_PCA_DIM, SEED + axis)
              for axis, value in enumerate(ligand_raw)]
    protein = [_pca_pipeline(value[protein_indices], PROTEIN_PCA_DIM, SEED + 17 + axis)
               for axis, value in enumerate(protein_raw)]
    return ViewTransforms(ligand=ligand, protein=protein)


@dataclass
class AdditiveModel:
    mu: float
    ligand: object
    protein: object

    def predict(self, ligand_features: np.ndarray,
                protein_features: np.ndarray) -> np.ndarray:
        ligand = np.asarray(self.ligand.predict(ligand_features), dtype=np.float64)
        protein = np.asarray(self.protein.predict(protein_features), dtype=np.float64)
        return self.mu + ligand[:, None] + protein[None, :]


def _fit_additive(y: np.ndarray, ligand_features: np.ndarray,
                  protein_features: np.ndarray) -> AdditiveModel:
    mu = float(np.nanmean(y))
    ligand_effect = np.nanmean(y, axis=1) - mu
    protein_effect = np.nanmean(y, axis=0) - mu
    ligand = make_pipeline(
        StandardScaler(), Ridge(alpha=ADDITIVE_RIDGE, solver="lsqr", tol=1e-8)
    ).fit(ligand_features, ligand_effect)
    protein = make_pipeline(
        StandardScaler(), Ridge(alpha=ADDITIVE_RIDGE, solver="lsqr", tol=1e-8)
    ).fit(protein_features, protein_effect)
    return AdditiveModel(mu=mu, ligand=ligand, protein=protein)


class SectionConditionedBilinear(nn.Module):
    """Three joint pair surfaces and a four-coordinate DeepSets section."""

    def __init__(self, ligand_dims, protein_dims):
        super().__init__()
        self.ligand_projection = nn.ModuleList([
            nn.Linear(int(width), BILINEAR_WIDTH, bias=False) for width in ligand_dims
        ])
        self.protein_projection = nn.ModuleList([
            nn.Linear(int(width), BILINEAR_WIDTH, bias=False) for width in protein_dims
        ])
        self.pair_hidden = nn.Linear(
            len(VIEW_NAMES) * BILINEAR_WIDTH, PAIR_WIDTH, bias=False)
        self.surface_head = nn.Linear(PAIR_WIDTH, 3, bias=False)
        self.element = nn.Sequential(
            nn.Linear(4, SET_WIDTH), nn.SiLU(),
            nn.Linear(SET_WIDTH, SET_WIDTH), nn.SiLU(),
        )
        self.section = nn.Sequential(
            nn.Linear(2 * SET_WIDTH, SET_WIDTH), nn.SiLU(),
            nn.Linear(SET_WIDTH, 4),
        )

    def surfaces(self, ligand_views, protein_views):
        blocks = []
        for ligand, protein, ligand_layer, protein_layer in zip(
                ligand_views, protein_views,
                self.ligand_projection, self.protein_projection):
            blocks.append(torch.tanh(ligand_layer(ligand))
                          * torch.tanh(protein_layer(protein)))
        hidden = torch.tanh(self.pair_hidden(torch.cat(blocks, dim=-1)))
        return 0.35 * torch.tanh(self.surface_head(hidden))

    def adapt(self, support_surfaces, support_residual):
        if support_surfaces.ndim != 3 or support_residual.ndim != 2:
            raise ValueError("support tensors must be batched as (episode, support, feature)")
        value = torch.cat([
            support_surfaces,
            support_residual.unsqueeze(-1) - support_surfaces[..., :1],
        ], dim=-1)
        encoded = self.element(value)
        aggregate = torch.cat([encoded.mean(dim=1), encoded.square().mean(dim=1)], dim=-1)
        raw = self.section(aggregate)
        tau = 0.35 * torch.tanh(raw[:, 0])
        tangent = torch.tanh(raw[:, 1:3])
        confidence = torch.sigmoid(raw[:, 3])
        return torch.cat([tau[:, None], tangent, confidence[:, None]], dim=-1)

    @staticmethod
    def apply_section(query_surfaces, coordinates):
        tau = coordinates[:, None, 0]
        tangent = coordinates[:, None, 1:3]
        confidence = coordinates[:, None, 3]
        correction = tau + (query_surfaces[..., 1:3] * tangent).sum(dim=-1)
        return query_surfaces[..., 0] + confidence * correction


def _tensor_views(value: list[np.ndarray]) -> list[torch.Tensor]:
    return [torch.as_tensor(item, dtype=torch.float32) for item in value]


class EpisodeSampler:
    def __init__(self, y: np.ndarray, scaffolds, k: int, query_size: int):
        self.y = np.asarray(y, dtype=np.float32)
        self.scaffolds = np.asarray(scaffolds, dtype=object)
        self.k = int(k)
        self.query_size = int(query_size)
        self.by_target = []
        for target in range(self.y.shape[1]):
            by_scaffold = {}
            for ligand in np.flatnonzero(np.isfinite(self.y[:, target])):
                by_scaffold.setdefault(str(self.scaffolds[ligand]), []).append(int(ligand))
            self.by_target.append(by_scaffold)
        self.valid_targets = np.asarray([
            target for target, groups in enumerate(self.by_target) if len(groups) >= self.k + 1
        ], dtype=np.int64)
        if not len(self.valid_targets):
            raise ValueError("no target has enough distinct scaffolds for an episode")

    def batch(self, rng: np.random.Generator, batch_size: int):
        targets, supports, queries = [], [], []
        attempts = 0
        while len(targets) < batch_size:
            attempts += 1
            if attempts > batch_size * 100:
                raise RuntimeError("could not sample enough valid episodes")
            target = int(rng.choice(self.valid_targets))
            groups = self.by_target[target]
            chosen = rng.choice(sorted(groups), size=self.k, replace=False)
            support = np.asarray([rng.choice(groups[str(key)]) for key in chosen], dtype=np.int64)
            forbidden = {str(value) for value in chosen}
            eligible = np.asarray([
                ligand for ligand in np.flatnonzero(np.isfinite(self.y[:, target]))
                if str(self.scaffolds[ligand]) not in forbidden
            ], dtype=np.int64)
            if len(eligible) < 3:
                continue
            size = min(self.query_size, len(eligible))
            query = rng.choice(eligible, size=size, replace=False)
            if size < self.query_size:
                query = np.resize(query, self.query_size)
            targets.append(target)
            supports.append(support)
            queries.append(query)
        return (np.asarray(targets, dtype=np.int64),
                np.asarray(supports, dtype=np.int64),
                np.asarray(queries, dtype=np.int64))


def _pair(model, ligand_views, protein_views, ligand_indices, target_indices):
    ligand_indices = torch.as_tensor(ligand_indices, dtype=torch.long)
    target_indices = torch.as_tensor(target_indices, dtype=torch.long)
    return model.surfaces(
        [value[ligand_indices] for value in ligand_views],
        [value[target_indices] for value in protein_views],
    )


def _episode_loss(model, ligand_views, protein_views, residual, sampler, rng,
                  batch_size=EPISODE_BATCH):
    targets, support, query = sampler.batch(rng, batch_size)
    repeated_support_target = np.repeat(targets, sampler.k)
    repeated_query_target = np.repeat(targets, sampler.query_size)
    support_surface = _pair(
        model, ligand_views, protein_views, support.reshape(-1), repeated_support_target
    ).reshape(batch_size, sampler.k, 3)
    query_surface = _pair(
        model, ligand_views, protein_views, query.reshape(-1), repeated_query_target
    ).reshape(batch_size, sampler.query_size, 3)
    support_y = torch.as_tensor(
        residual[support, targets[:, None]], dtype=torch.float32)
    query_y = torch.as_tensor(
        residual[query, targets[:, None]], dtype=torch.float32)
    coordinates = model.adapt(support_surface, support_y)
    prediction = model.apply_section(query_surface, coordinates)
    return torch.mean((prediction - query_y) ** 2), coordinates


def _auxiliary_losses(model, ligand_views, protein_views, residual,
                      rng, pair_count=576):
    n_ligand, n_target = residual.shape
    ligand = rng.integers(0, n_ligand, size=pair_count)
    target = rng.integers(0, n_target, size=pair_count)
    surface = _pair(model, ligand_views, protein_views, ligand, target)
    outcome = torch.as_tensor(residual[ligand, target], dtype=torch.float32)
    zero_loss = torch.mean((surface[:, 0] - outcome) ** 2)

    rectangles = max(32, pair_count // 4)
    i = rng.integers(0, n_ligand, size=rectangles)
    j = rng.integers(0, n_ligand, size=rectangles)
    t = rng.integers(0, n_target, size=rectangles)
    u = rng.integers(0, n_target, size=rectangles)
    predicted = (
        _pair(model, ligand_views, protein_views, i, t)[:, 0]
        - _pair(model, ligand_views, protein_views, i, u)[:, 0]
        - _pair(model, ligand_views, protein_views, j, t)[:, 0]
        + _pair(model, ligand_views, protein_views, j, u)[:, 0]
    )
    observed = torch.as_tensor(
        residual[i, t] - residual[i, u] - residual[j, t] + residual[j, u],
        dtype=torch.float32,
    )
    dd_loss = torch.mean((predicted - observed) ** 2)
    return zero_loss, dd_loss


@torch.no_grad()
def _validation_loss(model, ligand_views, protein_views, residual, sampler,
                     seed=SEED + 991):
    model.eval()
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(8):
        loss, _ = _episode_loss(
            model, ligand_views, protein_views, residual, sampler, rng, batch_size=16)
        values.append(float(loss))
    model.train()
    return float(np.mean(values))


def _train(ligand_views, protein_views, residual, scaffolds, steps,
           validation=None, seed=SEED):
    torch.manual_seed(seed)
    np.random.seed(seed)
    model = SectionConditionedBilinear(
        [value.shape[1] for value in ligand_views],
        [value.shape[1] for value in protein_views],
    )
    ligand_tensor = _tensor_views(ligand_views)
    protein_tensor = _tensor_views(protein_views)
    sampler = EpisodeSampler(residual, scaffolds, K_PRIMARY, QUERY_SIZE)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=1e-3, weight_decay=1e-4)
    rng = np.random.default_rng(seed)
    history = []
    best = {"step": 0, "loss": np.inf, "state": None}
    stale = 0

    if validation is not None:
        valid_ligand, valid_protein, valid_residual, valid_scaffold = validation
        valid_ligand_tensor = _tensor_views(valid_ligand)
        valid_protein_tensor = _tensor_views(valid_protein)
        valid_sampler = EpisodeSampler(
            valid_residual, valid_scaffold, K_PRIMARY, QUERY_SIZE)

    for step in range(1, int(steps) + 1):
        model.train()
        episodic, coordinates = _episode_loss(
            model, ligand_tensor, protein_tensor, residual, sampler, rng)
        zero, dd = _auxiliary_losses(
            model, ligand_tensor, protein_tensor, residual, rng)
        coordinate_penalty = torch.mean(coordinates[:, :3] ** 2)
        loss = episodic + 0.5 * zero + 0.25 * dd + 1e-3 * coordinate_penalty
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        optimizer.step()

        if validation is not None and step % VALIDATE_EVERY == 0:
            valid = _validation_loss(
                model, valid_ligand_tensor, valid_protein_tensor,
                valid_residual, valid_sampler)
            item = {"step": step, "training_loss": float(loss.detach()),
                    "episodic_loss": float(episodic.detach()),
                    "support_free_loss": float(zero.detach()),
                    "double_difference_loss": float(dd.detach()),
                    "validation_loss": valid}
            history.append(item)
            print(f"step={step} train={item['training_loss']:.6f} "
                  f"valid={valid:.6f}", flush=True)
            if valid < best["loss"] - 1e-6:
                best = {
                    "step": step, "loss": valid,
                    "state": {key: value.detach().clone()
                              for key, value in model.state_dict().items()},
                }
                stale = 0
            else:
                stale += 1
            if stale >= PATIENCE:
                break

    if validation is not None:
        if best["state"] is None:
            raise RuntimeError("validation did not create a checkpoint")
        model.load_state_dict(best["state"])
        return model, {"best_step": best["step"],
                       "best_validation_loss": best["loss"], "history": history}
    return model, {"steps": int(steps)}


@torch.no_grad()
def _all_surfaces(model, ligand_views, protein_views, chunk=8192):
    model.eval()
    n_ligand, n_target = len(ligand_views[0]), len(protein_views[0])
    ligand_tensor, protein_tensor = _tensor_views(ligand_views), _tensor_views(protein_views)
    flat_ligand = np.repeat(np.arange(n_ligand), n_target)
    flat_target = np.tile(np.arange(n_target), n_ligand)
    values = []
    for start in range(0, len(flat_ligand), chunk):
        stop = min(len(flat_ligand), start + chunk)
        values.append(_pair(
            model, ligand_tensor, protein_tensor,
            flat_ligand[start:stop], flat_target[start:stop],
        ).cpu().numpy())
    return np.concatenate(values).reshape(n_ligand, n_target, 3)


@torch.no_grad()
def _null_surfaces(model, ligand_views, protein_views):
    model.eval()
    n_ligand, n_target = len(ligand_views[0]), len(protein_views[0])
    ligand_tensor = _tensor_views(ligand_views)
    zero_protein = [torch.zeros((n_target, value.shape[1]), dtype=torch.float32)
                    for value in protein_views]
    flat_ligand = np.repeat(np.arange(n_ligand), n_target)
    flat_target = np.tile(np.arange(n_target), n_ligand)
    values = _pair(
        model, ligand_tensor, zero_protein, flat_ligand, flat_target).cpu().numpy()
    return values.reshape(n_ligand, n_target, 3)


def _nearest_nonself(protein_views: list[np.ndarray]) -> np.ndarray:
    feature = np.concatenate(protein_views, axis=1).astype(np.float64)
    squared = np.sum((feature[:, None, :] - feature[None, :, :]) ** 2, axis=-1)
    np.fill_diagonal(squared, np.inf)
    return np.argmin(squared, axis=1).astype(np.int64)


def _d_optimal(scaffolds, factor, finite, k, policy, seed):
    finite = np.asarray(finite, dtype=np.int64)
    by_scaffold = {}
    for index in finite:
        by_scaffold.setdefault(str(scaffolds[index]), []).append(int(index))
    if len(by_scaffold) < k:
        return None
    rng = np.random.default_rng(seed)
    if policy == "random":
        selected_scaffolds = rng.choice(sorted(by_scaffold), size=k, replace=False)
        return np.sort(np.asarray([
            rng.choice(by_scaffold[str(key)]) for key in selected_scaffolds
        ], dtype=np.int64))
    if policy != "d_optimal":
        raise ValueError(policy)
    representatives = np.asarray([
        min(indices) for _, indices in sorted(by_scaffold.items())
    ], dtype=np.int64)
    design = np.concatenate([
        np.ones((len(representatives), 1)), factor[representatives]
    ], axis=1)
    gram = np.eye(design.shape[1]) * 1e-3
    selected, remaining = [], list(range(len(representatives)))
    for _ in range(k):
        scores = []
        for candidate in remaining:
            sign, value = np.linalg.slogdet(
                gram + np.outer(design[candidate], design[candidate]))
            scores.append(value if sign > 0 else -np.inf)
        best = max(scores)
        tied = [remaining[index] for index, value in enumerate(scores)
                if np.isclose(value, best, rtol=1e-12, atol=1e-14)]
        chosen = min(tied)
        selected.append(chosen)
        gram += np.outer(design[chosen], design[chosen])
        remaining.remove(chosen)
    return np.sort(representatives[np.asarray(selected, dtype=np.int64)])


def _query(y, target, support, scaffolds):
    forbidden = {str(scaffolds[index]) for index in support}
    return np.asarray([
        ligand for ligand in range(y.shape[0])
        if np.isfinite(y[ligand, target])
        and ligand not in set(support)
        and str(scaffolds[ligand]) not in forbidden
    ], dtype=np.int64)


@torch.no_grad()
def _apply_numpy(model, support_surface, support_residual, query_surface):
    support_surface = torch.as_tensor(support_surface[None], dtype=torch.float32)
    support_residual = torch.as_tensor(support_residual[None], dtype=torch.float32)
    query_surface = torch.as_tensor(query_surface[None], dtype=torch.float32)
    coordinate = model.adapt(support_surface, support_residual)
    prediction = model.apply_section(query_surface, coordinate)
    return prediction[0].cpu().numpy(), coordinate[0].cpu().numpy()


def _safe_spearman(y, prediction):
    if len(y) < 3 or np.std(y) < 1e-12 or np.std(prediction) < 1e-12:
        return 0.0
    return float(spearmanr(y, prediction).statistic)


def _evaluate(model, y, additive, surfaces, null_surfaces, protein_views,
              scaffolds, k, policy, seeds):
    names = (
        "support_free", "location_only", "conditional", "zero_protein",
        "nearest_protein", "wrong_support", "permuted_support",
    )
    n_target = y.shape[1]
    per_target = {name: {metric: np.full((n_target, len(seeds)), np.nan)
                         for metric in ("mse", "mae", "spearman", "interaction_mse")}
                  for name in names}
    coordinates = np.full((n_target, len(seeds), 4), np.nan)
    condition = np.full((n_target, len(seeds)), np.nan)
    nearest = _nearest_nonself(protein_views)
    interaction_y = _masked_double_center(y)

    for target in range(n_target):
        wrong_target = int(nearest[target])
        jointly_finite = np.flatnonzero(
            np.isfinite(y[:, target]) & np.isfinite(y[:, wrong_target]))
        for seed_index, seed in enumerate(seeds):
            support = _d_optimal(
                scaffolds, surfaces[:, target, 1:3], jointly_finite, k,
                policy, seed + 104729 * target)
            if support is None:
                continue
            query = _query(y, target, support, scaffolds)
            if len(query) < 3:
                continue
            residual = y[support, target] - additive[support, target]
            wrong_residual = (
                y[support, wrong_target] - additive[support, wrong_target])
            rng = np.random.default_rng(seed + 65537 * target)

            correct, coordinate = _apply_numpy(
                model, surfaces[support, target], residual,
                surfaces[query, target])
            no_protein, _ = _apply_numpy(
                model, null_surfaces[support, target], residual,
                null_surfaces[query, target])
            nearest_prediction, _ = _apply_numpy(
                model, surfaces[support, wrong_target], residual,
                surfaces[query, wrong_target])
            wrong_support, _ = _apply_numpy(
                model, surfaces[support, target], wrong_residual,
                surfaces[query, target])
            permuted, _ = _apply_numpy(
                model, surfaces[support, target], residual[rng.permutation(k)],
                surfaces[query, target])
            location = float(np.sum(residual) / (k + LOCATION_RIDGE))
            residual_predictions = {
                "support_free": surfaces[query, target, 0],
                "location_only": np.full(len(query), location),
                "conditional": correct,
                "zero_protein": no_protein,
                "nearest_protein": nearest_prediction,
                "wrong_support": wrong_support,
                "permuted_support": permuted,
            }
            raw_predictions = {
                name: np.clip(additive[query, target] + value, 0.0, 1.0)
                for name, value in residual_predictions.items()
            }
            outcome = y[query, target]
            interaction_outcome = interaction_y[query, target]
            for name, prediction in raw_predictions.items():
                error = outcome - prediction
                per_target[name]["mse"][target, seed_index] = float(np.mean(error ** 2))
                per_target[name]["mae"][target, seed_index] = float(np.mean(np.abs(error)))
                per_target[name]["spearman"][target, seed_index] = _safe_spearman(
                    outcome, prediction)
                interaction_error = interaction_outcome - residual_predictions[name]
                per_target[name]["interaction_mse"][target, seed_index] = float(
                    np.mean(interaction_error ** 2))
            coordinates[target, seed_index] = coordinate
            design = np.concatenate([
                np.ones((k, 1)), surfaces[support, target, 1:3]
            ], axis=1)
            singular = np.linalg.svd(design, compute_uv=False)
            condition[target, seed_index] = float(singular[-1] / max(singular[0], 1e-12))
    return per_target, coordinates, condition, nearest


def _masked_double_center(value):
    value = np.asarray(value, dtype=np.float64)
    row = np.nanmean(value, axis=1, keepdims=True)
    column = np.nanmean(value, axis=0, keepdims=True)
    return value - row - column + float(np.nanmean(value))


def _bootstrap(comparator, candidate, metric="mse", seed=SEED):
    difference = (np.nanmean(comparator[metric], axis=1)
                  - np.nanmean(candidate[metric], axis=1))
    difference = difference[np.isfinite(difference)]
    if not len(difference):
        return {"estimate": None, "ci95": [None, None], "n_targets": 0}
    rng = np.random.default_rng(seed)
    draws = difference[rng.integers(
        0, len(difference), size=(BOOTSTRAP_DRAWS, len(difference)))].mean(axis=1)
    return {
        "estimate": float(np.mean(difference)),
        "ci95": [float(np.quantile(draws, 0.025)),
                 float(np.quantile(draws, 0.975))],
        "n_targets": int(len(difference)), "draws": BOOTSTRAP_DRAWS, "seed": seed,
    }


def _summarize(per_target, coordinates, condition, nearest):
    metrics = {
        name: {metric: float(np.nanmean(array)) for metric, array in values.items()}
        for name, values in per_target.items()
    }
    contrasts = {
        comparator: _bootstrap(per_target[comparator], per_target["conditional"])
        for comparator in per_target if comparator != "conditional"
    }
    interaction = {
        comparator: _bootstrap(
            per_target[comparator], per_target["conditional"], metric="interaction_mse")
        for comparator in ("support_free", "nearest_protein")
    }
    finite_condition = condition[np.isfinite(condition)]
    coordinate_summary = {}
    for index, name in enumerate(("tau", "u1", "u2", "confidence")):
        values = coordinates[..., index]
        coordinate_summary[name] = {
            "mean": float(np.nanmean(values)), "std": float(np.nanstd(values)),
            "q05": float(np.nanquantile(values, 0.05)),
            "q95": float(np.nanquantile(values, 0.95)),
        }
    return {
        "metrics": metrics,
        "mse_reduction_conditional_vs": contrasts,
        "interaction_mse_reduction_conditional_vs": interaction,
        "coordinates": coordinate_summary,
        "condition_ratio": {
            "median": float(np.median(finite_condition)),
            "q10": float(np.quantile(finite_condition, 0.1)),
            "q90": float(np.quantile(finite_condition, 0.9)),
        },
        "nearest_nonself_target_index": nearest.tolist(),
    }


def _source_split(y, ligand_raw, protein_raw, scaffolds, protein_groups):
    ligand_folds = np.asarray([stable_fold(value, 3) for value in scaffolds])
    protein_folds = np.asarray([stable_fold(value, 3) for value in protein_groups])
    ligand_train, ligand_valid = np.flatnonzero(ligand_folds != 0), np.flatnonzero(ligand_folds == 0)
    protein_train, protein_valid = np.flatnonzero(protein_folds != 0), np.flatnonzero(protein_folds == 0)
    transform = _fit_view_transforms(
        ligand_raw, protein_raw, ligand_train, protein_train)
    all_ligand = transform.ligand_apply(ligand_raw)
    all_protein = transform.protein_apply(protein_raw)
    training = (
        [value[ligand_train] for value in all_ligand],
        [value[protein_train] for value in all_protein],
        double_center(y[np.ix_(ligand_train, protein_train)]),
        np.asarray(scaffolds, dtype=object)[ligand_train],
    )
    validation = (
        [value[ligand_valid] for value in all_ligand],
        [value[protein_valid] for value in all_protein],
        double_center(y[np.ix_(ligand_valid, protein_valid)]),
        np.asarray(scaffolds, dtype=object)[ligand_valid],
    )
    detail = {
        "fold": 0, "train_ligands": len(ligand_train),
        "valid_ligands": len(ligand_valid), "train_targets": len(protein_train),
        "valid_targets": len(protein_valid),
    }
    return training, validation, detail


def run(args):
    started = time.time()
    torch.set_num_threads(int(args.threads))
    data = _load(args)
    y_source = data["y_source"]
    source_ligands = data["source_ligands"]
    source_targets = data["source_targets"]
    source_scaffolds = np.asarray([
        record.generic_scaffold for record in source_ligands], dtype=object)
    source_groups = [item["record"].group for item in source_targets]
    weights = _subpocket_weights(args.kissim_distances)
    source_ligand_raw = _ligand_raw_views(source_ligands)
    source_protein_raw = _protein_raw_views(source_targets, weights)

    print("selecting training horizon on source-only dual-cold split", flush=True)
    training, validation, split_detail = _source_split(
        y_source, source_ligand_raw, source_protein_raw,
        source_scaffolds, source_groups)
    _, selection = _train(*training, steps=MAX_STEPS,
                          validation=validation, seed=SEED)
    selected_steps = int(selection["best_step"])
    print(f"selected_steps={selected_steps}; fitting full PKIS1", flush=True)

    full_transform = _fit_view_transforms(source_ligand_raw, source_protein_raw)
    full_ligand = full_transform.ligand_apply(source_ligand_raw)
    full_protein = full_transform.protein_apply(source_protein_raw)
    model, full_training = _train(
        full_ligand, full_protein, double_center(y_source), source_scaffolds,
        steps=selected_steps, validation=None, seed=SEED)

    source_x_ligand = np.stack([
        record.nuisance_features for record in source_ligands])
    source_x_protein = np.stack([
        protein_model_features(item["record"].pocket) for item in source_targets])
    additive_model = _fit_additive(y_source, source_x_ligand, source_x_protein)

    transfers = {}
    for panel in data["panels"]:
        print(f"evaluating {panel['name']}", flush=True)
        ligand_raw = _ligand_raw_views(panel["ligands"])
        protein_raw = _protein_raw_views(panel["targets"], weights)
        ligand_views = full_transform.ligand_apply(ligand_raw)
        protein_views = full_transform.protein_apply(protein_raw)
        surfaces = _all_surfaces(model, ligand_views, protein_views)
        null = _null_surfaces(model, ligand_views, protein_views)
        additive = additive_model.predict(
            np.stack([record.nuisance_features for record in panel["ligands"]]),
            np.stack([protein_model_features(item["record"].pocket)
                      for item in panel["targets"]]),
        )
        scaffolds = np.asarray([
            record.generic_scaffold for record in panel["ligands"]], dtype=object)
        panel_result = {
            "role": panel["role"], "n_ligands": len(panel["ligands"]),
            "n_targets": len(panel["targets"]),
            "n_finite_cells": int(np.isfinite(panel["y"]).sum()),
            "support_sizes": {},
        }
        for k in (5, 20):
            panel_result["support_sizes"][str(k)] = {}
            for policy in ("random", "d_optimal"):
                seeds = SEEDS if policy == "random" else (SEED,)
                evaluated = _evaluate(
                    model, panel["y"], additive, surfaces, null, protein_views,
                    scaffolds, k, policy, seeds)
                panel_result["support_sizes"][str(k)][policy] = _summarize(*evaluated)
        transfers[panel["name"]] = panel_result

    primary = transfers["PKIS2"]["support_sizes"]["5"]["d_optimal"]
    external = transfers["Anastassiadis2011"]["support_sizes"]["5"]["d_optimal"]
    required = (
        "support_free", "location_only", "zero_protein", "nearest_protein",
        "wrong_support", "permuted_support",
    )
    pkis_raw = all(
        primary["mse_reduction_conditional_vs"][name]["ci95"][0] > 0.0
        for name in required)
    external_raw = all(
        external["mse_reduction_conditional_vs"][name]["estimate"] > 0.0
        for name in required)
    pkis_interaction = all(
        primary["interaction_mse_reduction_conditional_vs"][name]["ci95"][0] > 0.0
        for name in ("support_free", "nearest_protein"))
    passed = bool(pkis_raw and external_raw and pkis_interaction)

    result = {
        "schema": "MetaSieve.SectionConditionedBilinear.F1B.v1",
        "selected": {
            "d_adapt": 4, "k_primary": K_PRIMARY,
            "selected_training_steps": selected_steps,
            "ligand_pca_dim": LIGAND_PCA_DIM,
            "protein_pca_dim": PROTEIN_PCA_DIM,
            "bilinear_width": BILINEAR_WIDTH,
            "view_names": VIEW_NAMES,
            "location_ridge": LOCATION_RIDGE,
        },
        "source": {
            "panel": "PKIS1", "n_ligands": len(source_ligands),
            "n_targets": len(source_targets),
            "interaction_variance_fraction": float(
                np.var(double_center(y_source)) / np.var(y_source)),
            "split": split_detail, "selection": selection,
            "full_training": full_training,
        },
        "transfers": transfers,
        "gate": {
            "pkis2_raw_all_controls": pkis_raw,
            "anastassiadis_raw_point_estimates": external_raw,
            "pkis2_interaction_controls": pkis_interaction,
            "passed": passed,
            "verdict": ("F1B_BIOLOGICAL_SECTION_ADMISSIBLE" if passed
                        else "F1B_BIOLOGICAL_SECTION_NOT_ADMISSIBLE"),
        },
        "read_firewall": {
            "pkis1": "source", "pkis2": "consumed_development",
            "anastassiadis": "consumed_development",
            "kcgs_numeric_outcomes": "NOT_READ", "davis_labels": "NOT_READ",
            "recipient_labels": "NOT_READ",
        },
        "limitations": [
            "The affinity decoder is a diagnostic bridge, not the frozen law-valued operator.",
            "PKIS2 and Anastassiadis2011 are consumed development panels.",
            "Random and D-optimal raw MSE are not directly contrasted because scaffold exclusion changes their query sets.",
        ],
        "elapsed_seconds": time.time() - started,
        "software": {"torch": torch.__version__, "numpy": np.__version__},
    }
    output = Path(args.output).resolve()
    _write_json(output / "result.json", result)
    output.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), output / "model_state.pt")
    _write_json(output / "manifest.json", {
        "sha256": {
            "preregistration": _sha256(Path(__file__).with_name("F1B_PREREGISTRATION.md")),
            "script": _sha256(Path(__file__)),
            "kissim_distances": _sha256(Path(args.kissim_distances).resolve()),
            "pkis1_labels": _sha256(
                Path(args.informers_root).resolve() / "data/pkis1_continuous_labels.csv"),
        },
        "parameters": vars(args),
    })
    print(json.dumps(result["gate"], indent=2), flush=True)


def parser():
    item = argparse.ArgumentParser()
    item.add_argument("--informers-root", default="../external/informers")
    item.add_argument("--klifs-json", default="../external/klifs/kinase_information_human.json")
    item.add_argument("--kissim-distances",
                      default="../external/kissim/kissim/data/min_max_distances_fine.csv")
    item.add_argument("--anastassiadis-workbook",
                      default="external/anastassiadis/NIHMS328213-supplement-3.xls")
    item.add_argument("--anastassiadis-identities",
                      default="external/anastassiadis/compound_identities.json")
    item.add_argument("--threads", default="4")
    item.add_argument("--output", default="research/section_operator_pilot/artifacts/f1b")
    return item


if __name__ == "__main__":
    run(parser().parse_args())
