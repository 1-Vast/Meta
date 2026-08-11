"""Small, testable primitives for the preregistered R0 distance Gate.

This module contains no data loading and no affinity code.  It keeps scoring,
label-only ceiling calculation, and identity interventions independent from the
GPU runner so their contracts can be tested before any model result exists.
"""
from __future__ import annotations

from collections import defaultdict
import hashlib
from typing import Hashable, Sequence

import numpy as np

from contracts.mechanism import DISTANCE_BINS_ANGSTROM


def distance_bin_labels(distances: np.ndarray) -> np.ndarray:
    """Map Angstrom distances to the frozen P1B ordered-bin contract."""
    values = np.asarray(distances)
    if values.ndim != 2 or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("distances must be a finite nonnegative [atom,residue] matrix")
    edges = np.asarray(DISTANCE_BINS_ANGSTROM, dtype=np.float64)
    if float(values.max(initial=0.0)) >= edges[-1]:
        raise ValueError("distance exceeds the frozen overflow edge")
    return np.digitize(values, edges[1:-1]).astype(np.int64)


def ranked_probability_score(probability: np.ndarray,
                             labels: np.ndarray) -> np.ndarray:
    """Return pair-level RPS for an ordered categorical posterior."""
    probability = np.asarray(probability, dtype=np.float64)
    labels = np.asarray(labels)
    if probability.ndim < 2 or labels.shape != probability.shape[:-1]:
        raise ValueError("probability and label shapes disagree")
    bins = probability.shape[-1]
    if bins < 2 or not np.isfinite(probability).all() or (probability < 0).any():
        raise ValueError("probability must be finite, nonnegative, and have >=2 bins")
    total = probability.sum(axis=-1)
    if not np.allclose(total, 1.0, atol=1e-6):
        raise ValueError("active posterior rows must sum to one")
    if labels.dtype.kind not in "iu" or ((labels < 0) | (labels >= bins)).any():
        raise ValueError("labels must be integer bin indices")
    predicted_cdf = probability.cumsum(axis=-1)[..., :-1]
    thresholds = np.arange(bins - 1)
    observed_cdf = labels[..., None] <= thresholds
    return np.square(predicted_cdf - observed_cdf).mean(axis=-1)


def negative_log_likelihood(probability: np.ndarray,
                            labels: np.ndarray) -> np.ndarray:
    """Return pair-level full-bin NLL with finite clipping for reports."""
    probability = np.asarray(probability, dtype=np.float64)
    labels = np.asarray(labels)
    # Reuse all posterior and label validation in the primary score.
    ranked_probability_score(probability, labels)
    selected = np.take_along_axis(probability, labels[..., None], axis=-1)[..., 0]
    return -np.log(np.clip(selected, np.finfo(np.float64).tiny, 1.0))


def component_macro(per_system: dict[Hashable, float],
                    component_of: dict[Hashable, Hashable]) -> dict[Hashable, float]:
    """Average system scores inside components without pair-count weighting."""
    if set(per_system) - set(component_of):
        raise ValueError("one or more scored systems lack a closure component")
    grouped: dict[Hashable, list[float]] = defaultdict(list)
    for system, score in per_system.items():
        if not np.isfinite(score):
            raise ValueError("system scores must be finite")
        grouped[component_of[system]].append(float(score))
    if not grouped:
        raise ValueError("component macro requires at least one system")
    return {component: float(np.mean(values))
            for component, values in grouped.items()}


def _group_mean(values: np.ndarray, mask: np.ndarray, groups: np.ndarray,
                group_count: int) -> np.ndarray:
    """Masked group means used by the additive checkerboard backfit."""
    total = np.bincount(groups[mask], weights=values[mask], minlength=group_count)
    count = np.bincount(groups[mask], minlength=group_count)
    return np.divide(total, count, out=np.zeros(group_count, dtype=np.float64),
                     where=count > 0)


def _additive_threshold_prediction(target: np.ndarray, slot_of: np.ndarray, *,
                                   iterations: int = 30) -> np.ndarray:
    """Cross-fit atom + slot + within-slot-residue additive probabilities.

    Four parity checkerboards are held out in turn.  The residue effect is
    centered inside its slot after every update, which makes the nested slot and
    residue terms identifiable without reading a held-out cell.
    """
    atoms, residues = target.shape
    atom_index = np.broadcast_to(np.arange(atoms)[:, None], target.shape).ravel()
    residue_index = np.broadcast_to(np.arange(residues)[None, :], target.shape).ravel()
    slot_index = np.broadcast_to(slot_of[None, :], target.shape).ravel()
    y = target.astype(np.float64).ravel()
    folds = 2 * (atom_index % 2) + (residue_index % 2)
    prediction = np.empty_like(y)

    for fold in range(4):
        train = folds != fold
        test = ~train
        if not train.any() or not test.any():
            continue
        mu = float(y[train].mean())
        atom_effect = np.zeros(atoms, dtype=np.float64)
        slot_effect = np.zeros(int(slot_of.max()) + 1, dtype=np.float64)
        residue_effect = np.zeros(residues, dtype=np.float64)
        for _ in range(iterations):
            residual = y - mu - slot_effect[slot_index] - residue_effect[residue_index]
            atom_effect = _group_mean(residual, train, atom_index, atoms)
            residual = y - mu - atom_effect[atom_index] - residue_effect[residue_index]
            slot_effect = _group_mean(
                residual, train, slot_index, len(slot_effect))
            residual = y - mu - atom_effect[atom_index] - slot_effect[slot_index]
            residue_effect = _group_mean(residual, train, residue_index, residues)
            # Enforce sum-to-zero of exact-residue deviations inside each slot.
            for slot in np.unique(slot_of):
                members = np.flatnonzero(slot_of == slot)
                residue_effect[members] -= residue_effect[members].mean()
            fitted_without_mu = (atom_effect[atom_index] + slot_effect[slot_index]
                                 + residue_effect[residue_index])
            mu = float((y[train] - fitted_without_mu[train]).mean())
        prediction[test] = np.clip(
            mu + atom_effect[atom_index[test]] + slot_effect[slot_index[test]]
            + residue_effect[residue_index[test]], 0.0, 1.0)
    if not np.isfinite(prediction).all():
        raise RuntimeError("additive checkerboard did not score every cell")
    return prediction.reshape(target.shape)


def additive_checkerboard_rps(labels: np.ndarray, slot_of: np.ndarray, *,
                              iterations: int = 30) -> float:
    """Label-only cross-fit additive oracle under the R0 ordered RPS."""
    labels = np.asarray(labels)
    slot_of = np.asarray(slot_of)
    if labels.ndim != 2 or slot_of.shape != (labels.shape[1],):
        raise ValueError("labels need [atom,residue] and one slot per residue")
    if min(labels.shape) < 2:
        raise ValueError("checkerboard additive scoring needs >=2 atoms and residues")
    if labels.dtype.kind not in "iu" or ((labels < 0) | (labels >= 5)).any():
        raise ValueError("labels are outside the frozen five-bin contract")
    if slot_of.dtype.kind not in "iu" or (slot_of < 0).any():
        raise ValueError("slot indices must be nonnegative integers")
    cumulative = np.stack([
        _additive_threshold_prediction(labels <= threshold, slot_of,
                                       iterations=iterations)
        for threshold in range(4)
    ], axis=-1)
    observed = labels[..., None] <= np.arange(4)
    return float(np.square(cumulative - observed).mean())


def deterministic_derangement(groups: Sequence[Hashable], *,
                              namespace: str) -> np.ndarray:
    """Derange indices inside each non-singleton group with no fixed points."""
    if not namespace:
        raise ValueError("derangement namespace cannot be empty")
    grouped: dict[Hashable, list[int]] = defaultdict(list)
    for index, group in enumerate(groups):
        grouped[group].append(index)
    result = np.arange(len(groups), dtype=np.int64)
    for group, members in grouped.items():
        if len(members) < 2:
            continue
        ordered = sorted(members, key=lambda index: hashlib.sha256(
            f"{namespace}|{group}|{index}".encode("utf-8")).hexdigest())
        for position, index in enumerate(ordered):
            result[index] = ordered[(position + 1) % len(ordered)]
    return result


def paired_component_bootstrap(left: dict[Hashable, float],
                               right: dict[Hashable, float], *,
                               seed: int, draws: int = 10_000) -> dict[str, float | int]:
    """Bootstrap the paired improvement `left - right` by closure component."""
    keys = sorted(set(left) & set(right), key=str)
    if len(keys) < 2 or draws < 1:
        raise ValueError("paired bootstrap needs >=2 components and positive draws")
    delta = np.asarray([left[key] - right[key] for key in keys], dtype=np.float64)
    if not np.isfinite(delta).all():
        raise ValueError("component contrasts must be finite")
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(keys), size=(draws, len(keys)))
    samples = delta[indices].mean(axis=1)
    return {
        "delta": float(delta.mean()),
        "lcb95_one_sided": float(np.percentile(samples, 5.0)),
        "ucb95_one_sided": float(np.percentile(samples, 95.0)),
        "components": len(keys),
    }
