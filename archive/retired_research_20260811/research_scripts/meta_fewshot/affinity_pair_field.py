"""Research-only exact residue--atom field for biological admission tests.

The field keeps every local value attached to an explicit ligand atom and
protein residue.  The frozen P1B contact/distance output is an input prior, not
an exact-residue label: a small bilinear residual can refine that slot-level
prior with exact residue states before typed interaction energies are pooled.

This module deliberately does not import the meta-section or the production
operator.  It must pass structural and measured-affinity Gates before either
integration is allowed.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Hashable, Sequence

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor


INTERACTION_CHANNELS = (
    "hbond",
    "ionic",
    "hydrophobic",
    "aromatic",
    "vdw",
    "other",
)


@dataclass(frozen=True)
class PairFieldPrediction:
    """Auditable local field and its bounded-dimensional summary."""

    typed_pair_energy: Tensor
    typed_summary: Tensor
    score: Tensor
    distance_logits: Tensor
    distance_prob: Tensor
    contact_prob: Tensor
    geometry_gate: Tensor
    pair_mask: Tensor
    atom_indices: Tensor
    residue_indices: Tensor


@dataclass(frozen=True)
class PairDifference:
    """A set of measured two-cell differences: score[left]-score[right]."""

    left: Tensor
    right: Tensor
    target: Tensor


@dataclass(frozen=True)
class RectangleDifference:
    """Measured 2x2 differences ordered as (P1L1, P1L2, P2L1, P2L2)."""

    indices: Tensor
    target: Tensor


def _validate_exact_indices(indices: Tensor, mask: Tensor, name: str) -> None:
    if indices.shape != mask.shape:
        raise ValueError(f"{name} indices and mask must have the same shape")
    if indices.dtype not in (torch.int32, torch.int64):
        raise TypeError(f"{name} indices must be integer tensors")
    active = mask.bool()
    if bool((indices[active] < 0).any().item()):
        raise ValueError(f"active {name} indices must be nonnegative")
    if bool((indices[~active] != -1).any().item()):
        raise ValueError(f"masked {name} indices must use -1")
    for row_indices, row_mask in zip(indices, active):
        values = row_indices[row_mask]
        if values.numel() != torch.unique(values).numel():
            raise ValueError(f"active {name} indices must be unique per sample")


def coarse_interaction_compatibility(atom_classes: Tensor,
                                     residue_classes: Tensor) -> Tensor:
    """Return fixed, overlapping typed compatibility channels.

    Atom classes follow the retained T-BASIS order:
    hydrophobe/aromatic/donor/acceptor/positive/negative/halogen/other.
    Residue classes follow its six coarse chemistry groups:
    aliphatic/aromatic/polar/basic/acidic/other. Padding is encoded by -1.
    """
    if atom_classes.ndim != 2 or residue_classes.ndim != 2:
        raise ValueError("chemistry class tensors must have shape [B,A] and [B,R]")
    if atom_classes.shape[0] != residue_classes.shape[0]:
        raise ValueError("atom and residue chemistry batches differ")
    if atom_classes.dtype not in (torch.int32, torch.int64) or \
            residue_classes.dtype not in (torch.int32, torch.int64):
        raise TypeError("chemistry classes must be integer tensors")
    if bool(((atom_classes < -1) | (atom_classes > 7)).any().item()):
        raise ValueError("atom chemistry class is outside [-1, 7]")
    if bool(((residue_classes < -1) | (residue_classes > 5)).any().item()):
        raise ValueError("residue chemistry class is outside [-1, 5]")

    atom = atom_classes.unsqueeze(-1)
    residue = residue_classes.unsqueeze(-2)
    valid = (atom >= 0) & (residue >= 0)
    result = torch.zeros(
        *valid.shape, len(INTERACTION_CHANNELS),
        dtype=torch.float32, device=atom_classes.device,
    )
    result[..., 0] = (
        ((atom == 2) & ((residue == 2) | (residue == 4))) |
        ((atom == 3) & ((residue == 2) | (residue == 3)))
    ).float()
    result[..., 1] = (
        ((atom == 4) & (residue == 4)) |
        ((atom == 5) & (residue == 3))
    ).float()
    result[..., 2] = (
        ((atom == 0) | (atom == 6)) &
        ((residue == 0) | (residue == 1))
    ).float()
    result[..., 3] = (
        (atom == 1) & ((residue == 1) | (residue == 3))
    ).float()
    result[..., 4] = valid.float()
    result[..., 5] = (valid & ((atom == 7) | (residue == 5))).float()
    return result


def expand_slot_geometry_prior(contact_by_slot: Tensor, distance_by_slot: Tensor,
                               residue_slot: Tensor,
                               residue_mask: Tensor) -> tuple[Tensor, Tensor]:
    """Lift the frozen P1B slot prior onto explicit residue identities.

    Residues in one slot intentionally receive the same prior.  They remain
    separate rows and can only be distinguished by their exact residue state;
    this function never presents the lifted prior as exact structural truth.
    """
    if contact_by_slot.ndim != 3 or distance_by_slot.ndim != 4:
        raise ValueError("slot priors must have shape [B,A,S] and [B,A,S,D]")
    if contact_by_slot.shape != distance_by_slot.shape[:3]:
        raise ValueError("contact and distance slot priors disagree")
    if residue_slot.ndim != 2 or residue_mask.shape != residue_slot.shape:
        raise ValueError("residue slot map and mask must have shape [B,R]")
    if residue_slot.shape[0] != contact_by_slot.shape[0]:
        raise ValueError("slot priors and exact residues have different batches")
    if residue_slot.dtype not in (torch.int32, torch.int64):
        raise TypeError("residue slot map must be an integer tensor")
    active = residue_mask.bool()
    slot_count = contact_by_slot.shape[2]
    if bool(((residue_slot[active] < 0) |
             (residue_slot[active] >= slot_count)).any().item()):
        raise ValueError("active exact residue points outside the P1B slot axis")
    if bool((residue_slot[~active] != -1).any().item()):
        raise ValueError("masked exact residues must use slot -1")
    safe_slot = residue_slot.clamp_min(0).long()
    atoms = contact_by_slot.shape[1]
    contact_index = safe_slot.unsqueeze(1).expand(-1, atoms, -1)
    contact = torch.gather(contact_by_slot, 2, contact_index)
    distance_index = contact_index.unsqueeze(-1).expand(
        -1, -1, -1, distance_by_slot.shape[-1])
    distance = torch.gather(distance_by_slot, 2, distance_index)
    mask = active.unsqueeze(1).to(contact.dtype)
    return contact * mask, distance * mask.unsqueeze(-1)


class AffinityDirectedPairField(nn.Module):
    """Low-rank typed field with an exact-residue distance posterior.

    The architecture is intentionally narrow: one bilinear signal per typed
    channel, one linear residual over a frozen distance prior, and a six-value
    sum of local contributions.  There is no pooled protein/ligand shortcut and
    no free pair embedding.
    """

    channels = INTERACTION_CHANNELS

    def __init__(self, d_atom: int, d_residue: int, *, rank: int = 8,
                 n_distance_bins: int = 5, contact_bins: int = 2,
                 dtype: torch.dtype = torch.float32,
                 interaction_mode: str = "bilinear") -> None:
        super().__init__()
        if min(d_atom, d_residue, rank) < 1:
            raise ValueError("state dimensions and rank must be positive")
        if n_distance_bins < 2 or not 1 <= contact_bins < n_distance_bins:
            raise ValueError("contact bins must be a nonempty prefix of distance bins")
        self.d_atom = int(d_atom)
        self.d_residue = int(d_residue)
        self.rank = int(rank)
        self.n_distance_bins = int(n_distance_bins)
        self.contact_bins = int(contact_bins)
        if interaction_mode not in {"bilinear", "additive"}:
            raise ValueError("interaction_mode must be bilinear or additive")
        self.interaction_mode = interaction_mode
        channel_count = len(self.channels)
        self.atom_projection = nn.Linear(
            d_atom, channel_count * rank, bias=False, dtype=dtype)
        self.residue_projection = nn.Linear(
            d_residue, channel_count * rank, bias=False, dtype=dtype)
        self.distance_residual = nn.Linear(
            channel_count, n_distance_bins, bias=False, dtype=dtype)
        self.geometry_mix_logits = nn.Parameter(
            self._initial_geometry_mix(n_distance_bins, dtype))
        self.channel_coefficients = nn.Parameter(
            torch.full((channel_count,), 1.0 / channel_count, dtype=dtype))

    @staticmethod
    def _initial_geometry_mix(n_bins: int, dtype: torch.dtype) -> Tensor:
        # Columns are frozen contact prior, refined contact posterior, then bins.
        # The saturated binary teacher is available but is not initially dominant.
        result = torch.zeros(len(INTERACTION_CHANNELS), n_bins + 2, dtype=dtype)
        near = min(2, n_bins)
        result[:3, 1] = 2.0
        result[:3, 2:2 + near] = 1.5
        result[3, 1] = 1.0
        result[3, 2:2 + min(3, n_bins)] = 1.0
        result[4, 1:2 + min(3, n_bins)] = 0.75
        return result

    def forward(self, atom_states: Tensor, atom_mask: Tensor,
                residue_states: Tensor, residue_mask: Tensor,
                contact_prior: Tensor, distance_prior: Tensor,
                compatibility: Tensor, atom_indices: Tensor,
                residue_indices: Tensor) -> PairFieldPrediction:
        if atom_states.ndim != 3 or residue_states.ndim != 3:
            raise ValueError("states must have shape [B,A,D] and [B,R,D]")
        if atom_states.shape[-1] != self.d_atom or \
                residue_states.shape[-1] != self.d_residue:
            raise ValueError("state dimensions differ from the field contract")
        if atom_mask.shape != atom_states.shape[:2] or \
                residue_mask.shape != residue_states.shape[:2]:
            raise ValueError("state masks have incompatible shapes")
        if bool((atom_mask.sum(1) <= 0).any().item()) or \
                bool((residue_mask.sum(1) <= 0).any().item()):
            raise ValueError("every sample needs an atom and an exact residue")
        _validate_exact_indices(atom_indices, atom_mask, "atom")
        _validate_exact_indices(residue_indices, residue_mask, "residue")

        batch, atoms = atom_states.shape[:2]
        residues = residue_states.shape[1]
        pair_shape = (batch, atoms, residues)
        if contact_prior.shape != pair_shape:
            raise ValueError("contact prior must have shape [B,A,R]")
        if distance_prior.shape != (*pair_shape, self.n_distance_bins):
            raise ValueError("distance prior has the wrong shape")
        if compatibility.shape != (*pair_shape, len(self.channels)):
            raise ValueError("typed compatibility has the wrong shape")
        tensors = (contact_prior, distance_prior, compatibility)
        if any(not bool(torch.isfinite(value).all().item()) for value in tensors):
            raise ValueError("pair priors and compatibility must be finite")
        if bool(((contact_prior < 0) | (contact_prior > 1)).any().item()):
            raise ValueError("contact prior must lie in [0, 1]")
        if bool((distance_prior < 0).any().item()):
            raise ValueError("distance prior must be nonnegative")
        if bool(((compatibility < 0) | (compatibility > 1)).any().item()):
            raise ValueError("typed compatibility must lie in [0, 1]")

        pair_mask = atom_mask.bool().unsqueeze(-1) & residue_mask.bool().unsqueeze(-2)
        prior_total = distance_prior.sum(-1)
        if bool((prior_total[pair_mask] <= 0).any().item()):
            raise ValueError("active distance priors must carry positive mass")
        normalized_prior = distance_prior / prior_total.clamp_min(
            torch.finfo(distance_prior.dtype).tiny).unsqueeze(-1)

        channel_count = len(self.channels)
        atom = self.atom_projection(atom_states).reshape(
            batch, atoms, channel_count, self.rank)
        residue = self.residue_projection(residue_states).reshape(
            batch, residues, channel_count, self.rank)
        if self.interaction_mode == "bilinear":
            pair_signal = torch.einsum(
                "bacd,brcd->barc", atom, residue) / math.sqrt(self.rank)
            relation = pair_signal * compatibility.to(pair_signal.dtype)
        else:
            atom_marginal = atom.mean(dim=-1).unsqueeze(2)
            residue_marginal = residue.mean(dim=-1).unsqueeze(1)
            relation = (atom_marginal + residue_marginal) / math.sqrt(2.0)

        log_prior = normalized_prior.clamp_min(
            torch.finfo(normalized_prior.dtype).tiny).log()
        distance_logits = log_prior + self.distance_residual(relation)
        distance_prob = torch.softmax(distance_logits, dim=-1)
        mask_value = pair_mask.to(distance_prob.dtype)
        distance_prob = distance_prob * mask_value.unsqueeze(-1)
        contact_prob = distance_prob[..., :self.contact_bins].sum(-1)

        geometry_basis = torch.cat((
            contact_prior.to(distance_prob.dtype).unsqueeze(-1),
            contact_prob.unsqueeze(-1),
            distance_prob,
        ), dim=-1)
        geometry_weights = torch.softmax(self.geometry_mix_logits, dim=-1)
        geometry_gate = torch.einsum("barq,cq->barc", geometry_basis, geometry_weights)
        geometry_gate = geometry_gate * mask_value.unsqueeze(-1)
        typed_pair_energy = relation * geometry_gate * mask_value.unsqueeze(-1)
        typed_summary = typed_pair_energy.sum(dim=(1, 2))
        score = typed_summary @ self.channel_coefficients
        return PairFieldPrediction(
            typed_pair_energy=typed_pair_energy,
            typed_summary=typed_summary,
            score=score,
            distance_logits=distance_logits,
            distance_prob=distance_prob,
            contact_prob=contact_prob,
            geometry_gate=geometry_gate,
            pair_mask=pair_mask,
            atom_indices=atom_indices,
            residue_indices=residue_indices,
        )


def exact_distance_loss(prediction: PairFieldPrediction,
                        distance_bin: Tensor, *,
                        ordinal_weight: float = 0.25) -> Tensor:
    """NLL plus discrete CRPS on ordered exact residue--atom distances."""
    if ordinal_weight < 0:
        raise ValueError("ordinal distance weight must be nonnegative")
    if distance_bin.shape != prediction.pair_mask.shape:
        raise ValueError("exact distance labels must have shape [B,A,R]")
    if distance_bin.dtype not in (torch.int32, torch.int64):
        raise TypeError("exact distance bins must be integer tensors")
    labels = distance_bin[prediction.pair_mask]
    if labels.numel() == 0:
        raise ValueError("exact distance loss has no active pairs")
    bins = prediction.distance_logits.shape[-1]
    if bool(((labels < 0) | (labels >= bins)).any().item()):
        raise ValueError("exact distance label is outside the configured bins")
    logits = prediction.distance_logits[prediction.pair_mask]
    nll = F.cross_entropy(logits, labels.long())
    if ordinal_weight == 0:
        return nll
    probability = torch.softmax(logits, dim=-1)
    observed = F.one_hot(labels.long(), num_classes=bins).to(probability.dtype)
    # Discrete CRPS respects adjacency: moving mass by one bin is penalized
    # less than moving it across the full distance range.
    crps = torch.square(
        probability.cumsum(-1) - observed.cumsum(-1)
    ).mean()
    return nll + ordinal_weight * crps


def exact_distance_loss_per_system(
        prediction: PairFieldPrediction, distance_bin: Tensor, *,
        ordinal_weight: float = 0.25) -> Tensor:
    """Return one equal-weight NLL+CRPS loss per complex."""
    if ordinal_weight < 0:
        raise ValueError("ordinal distance weight must be nonnegative")
    if distance_bin.shape != prediction.pair_mask.shape:
        raise ValueError("exact distance labels must have shape [B,A,R]")
    if distance_bin.dtype not in (torch.int32, torch.int64):
        raise TypeError("exact distance bins must be integer tensors")
    mask = prediction.pair_mask
    bins = prediction.distance_logits.shape[-1]
    safe = distance_bin.clamp(0, bins - 1).long()
    if bool((((distance_bin < 0) | (distance_bin >= bins)) & mask).any().item()):
        raise ValueError("exact distance label is outside the configured bins")
    log_probability = torch.log_softmax(prediction.distance_logits, dim=-1)
    nll = -torch.gather(log_probability, -1, safe.unsqueeze(-1)).squeeze(-1)
    probability = torch.softmax(prediction.distance_logits, dim=-1)
    observed = F.one_hot(safe, num_classes=bins).to(probability.dtype)
    crps = torch.square(
        probability.cumsum(-1) - observed.cumsum(-1)
    ).mean(dim=-1)
    pair_loss = nll + ordinal_weight * crps
    weight = mask.to(pair_loss.dtype)
    count = weight.sum(dim=(1, 2))
    if bool((count <= 0).any().item()):
        raise ValueError("every complex needs at least one active pair")
    return (pair_loss * weight).sum(dim=(1, 2)) / count


def rectangle_values(scores: Tensor, rectangles: RectangleDifference) -> Tensor:
    indices = rectangles.indices
    if indices.ndim != 2 or indices.shape[1] != 4:
        raise ValueError("rectangle indices must have shape [N,4]")
    if indices.dtype not in (torch.int32, torch.int64):
        raise TypeError("rectangle indices must be integer tensors")
    if indices.numel() and bool(((indices < 0) | (indices >= len(scores))).any().item()):
        raise IndexError("rectangle index is outside the score vector")
    if any(torch.unique(row).numel() != 4 for row in indices):
        raise ValueError("every measured rectangle must contain four distinct cells")
    values = scores[indices.long()]
    return values[:, 0] - values[:, 1] - values[:, 2] + values[:, 3]


class AffinityContrastLoss(nn.Module):
    """Huber objective over measured differences; absolute affinity is opt-in."""

    def __init__(self, *, within_target_weight: float = 1.0,
                 cross_protein_weight: float = 1.0,
                 rectangle_weight: float = 1.0,
                 absolute_weight: float = 0.0, delta: float = 1.0) -> None:
        super().__init__()
        weights = (within_target_weight, cross_protein_weight,
                   rectangle_weight, absolute_weight)
        if any(value < 0 for value in weights) or delta <= 0:
            raise ValueError("loss weights must be nonnegative and delta positive")
        if not any(value > 0 for value in weights[:3]):
            raise ValueError("at least one measured contrast must have positive weight")
        self.weights = weights
        self.delta = float(delta)

    @staticmethod
    def _pair_values(scores: Tensor, pairs: PairDifference, name: str) -> Tensor:
        if pairs.left.ndim != 1 or pairs.right.ndim != 1 or pairs.target.ndim != 1:
            raise ValueError(f"{name} contrast tensors must be one-dimensional")
        if not (len(pairs.left) == len(pairs.right) == len(pairs.target)):
            raise ValueError(f"{name} contrast tensors must have equal lengths")
        if pairs.left.dtype not in (torch.int32, torch.int64) or \
                pairs.right.dtype not in (torch.int32, torch.int64):
            raise TypeError(f"{name} contrast indices must be integers")
        for indices in (pairs.left, pairs.right):
            if indices.numel() and bool(((indices < 0) | (indices >= len(scores))).any().item()):
                raise IndexError(f"{name} contrast index is outside the score vector")
        if bool((pairs.left == pairs.right).any().item()):
            raise ValueError(f"{name} measured differences need two distinct cells")
        return scores[pairs.left.long()] - scores[pairs.right.long()]

    def forward(self, scores: Tensor, *,
                within_target: PairDifference | None = None,
                cross_protein: PairDifference | None = None,
                rectangles: RectangleDifference | None = None,
                absolute_target: Tensor | None = None) -> dict[str, Tensor]:
        if scores.ndim != 1:
            raise ValueError("affinity scores must be one-dimensional")
        total = scores.sum() * 0.0
        output: dict[str, Tensor] = {}
        entries = (
            ("within_target", within_target, self.weights[0]),
            ("cross_protein", cross_protein, self.weights[1]),
        )
        for name, values, weight in entries:
            if values is None or weight == 0:
                continue
            predicted = self._pair_values(scores, values, name)
            loss = F.huber_loss(predicted, values.target.to(scores), delta=self.delta)
            output[name] = loss
            total = total + weight * loss
        if rectangles is not None and self.weights[2] > 0:
            predicted = rectangle_values(scores, rectangles)
            if rectangles.target.ndim != 1 or len(rectangles.target) != len(predicted):
                raise ValueError("rectangle targets must have shape [N]")
            loss = F.huber_loss(
                predicted, rectangles.target.to(scores), delta=self.delta)
            output["rectangle"] = loss
            total = total + self.weights[2] * loss
        if absolute_target is not None and self.weights[3] > 0:
            if absolute_target.shape != scores.shape:
                raise ValueError("absolute targets must match the score vector")
            loss = F.huber_loss(scores, absolute_target.to(scores), delta=self.delta)
            output["absolute"] = loss
            total = total + self.weights[3] * loss
        if not output:
            raise ValueError("the batch contains no enabled measured objective")
        output["total"] = total
        return output


def cluster_partner_necessity(correct_prediction: Tensor | np.ndarray,
                              wrong_prediction: Tensor | np.ndarray,
                              target: Tensor | np.ndarray,
                              clusters: Sequence[Hashable], *,
                              n_bootstrap: int = 5000,
                              seed: int = 17) -> dict[str, float | int]:
    """Cluster-bootstrap N_P = MSE(wrong protein)-MSE(correct protein)."""
    if n_bootstrap < 1:
        raise ValueError("n_bootstrap must be positive")

    def as_numpy(value) -> np.ndarray:
        if isinstance(value, Tensor):
            value = value.detach().cpu().numpy()
        return np.asarray(value, dtype=np.float64)

    correct, wrong, truth = map(as_numpy, (correct_prediction, wrong_prediction, target))
    if correct.ndim != 1 or correct.shape != wrong.shape or correct.shape != truth.shape:
        raise ValueError("partner predictions and targets must be equal one-dimensional arrays")
    if len(clusters) != len(correct):
        raise ValueError("one cluster label is required per prediction")
    if not np.isfinite(correct).all() or not np.isfinite(wrong).all() or \
            not np.isfinite(truth).all():
        raise ValueError("partner Gate inputs must be finite")
    labels = np.asarray(list(clusters), dtype=object)
    unique = list(dict.fromkeys(labels.tolist()))
    if len(unique) < 2:
        raise ValueError("cluster bootstrap needs at least two dependency clusters")
    row_delta = np.square(wrong - truth) - np.square(correct - truth)
    cluster_delta = np.asarray(
        [row_delta[labels == label].mean() for label in unique], dtype=np.float64)
    rng = np.random.default_rng(seed)
    draws = cluster_delta[
        rng.integers(0, len(cluster_delta), size=(n_bootstrap, len(cluster_delta)))
    ].mean(axis=1)
    return {
        "partner_necessity": float(cluster_delta.mean()),
        "lcb95_one_sided": float(np.percentile(draws, 5.0)),
        "clusters": len(unique),
        "rows": len(correct),
    }
