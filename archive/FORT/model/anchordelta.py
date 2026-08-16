"""Support-anchored, protein-conditioned relative affinity operator.

The encoder is intentionally kept outside the trainable head.  This makes the
P0 experiment a representation kill test: only the small comparator head is
fit while the existing protein-ligand representation stays frozen.
"""

from __future__ import annotations

from typing import Literal

import torch
from torch import nn
from torch.nn import functional as F


def _validate_features(left: torch.Tensor, right: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    left = left.float()
    right = right.float()
    if left.ndim != 2 or right.ndim != 2 or left.shape != right.shape:
        raise ValueError("pair features must have matching shape [rows, width]")
    if not bool(torch.isfinite(left).all()) or not bool(torch.isfinite(right).all()):
        raise ValueError("pair features must be finite")
    return left, right


class ComparatorHead(nn.Module):
    """Small ordered pair scorer used only inside the antisymmetric operator."""

    def __init__(self, feature_dim: int, hidden_dim: int = 64) -> None:
        super().__init__()
        if feature_dim <= 0 or hidden_dim <= 0:
            raise ValueError("feature_dim and hidden_dim must be positive")
        self.feature_dim = feature_dim
        self.network = nn.Sequential(
            nn.Linear(4 * feature_dim, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, left: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
        left, right = _validate_features(left, right)
        values = torch.cat(
            (left, right, left * right, (left - right).abs()), dim=-1
        )
        return self.network(values).squeeze(-1)


class AnchorDelta(nn.Module):
    """Antisymmetric comparator over frozen per-ligand interaction features.

    ``delta(q, i)`` is exactly the negative of ``delta(i, q)`` and is exactly
    zero on the diagonal, independent of the head parameters.
    """

    def __init__(self, feature_dim: int, hidden_dim: int = 64) -> None:
        super().__init__()
        self.head = ComparatorHead(feature_dim, hidden_dim)

    def delta(self, query: torch.Tensor, anchor: torch.Tensor) -> torch.Tensor:
        query, anchor = _validate_features(query, anchor)
        return 0.5 * (self.head(query, anchor) - self.head(anchor, query))

    def forward(self, query: torch.Tensor, anchor: torch.Tensor) -> torch.Tensor:
        return self.delta(query, anchor)

    def matrix(self, query: torch.Tensor, anchor: torch.Tensor) -> torch.Tensor:
        """Return ``[queries, anchors]`` deltas without requiring equal counts."""

        if query.ndim != 2 or anchor.ndim != 2 or query.shape[1] != anchor.shape[1]:
            raise ValueError("query and anchor features must be [rows, width]")
        queries = query[:, None, :].expand(-1, anchor.shape[0], -1).reshape(-1, query.shape[1])
        anchors = anchor[None, :, :].expand(query.shape[0], -1, -1).reshape(-1, anchor.shape[1])
        return self.delta(queries, anchors).reshape(query.shape[0], anchor.shape[0])


class EncodedAnchorDelta(AnchorDelta):
    """AnchorDelta wrapper around an existing ``InteractionEncoder``.

    The interaction encoder is frozen by default.  It must expose the existing
    ``encodeprotein`` and ``pairfromprotein`` methods, so no new protein path is
    introduced by the P0 experiment.
    """

    def __init__(
        self,
        interaction: nn.Module,
        feature_dim: int,
        hidden_dim: int = 64,
        *,
        freeze_encoder: bool = True,
    ) -> None:
        super().__init__(feature_dim, hidden_dim)
        if not hasattr(interaction, "encodeprotein") or not hasattr(interaction, "pairfromprotein"):
            raise TypeError("interaction must expose encodeprotein and pairfromprotein")
        self.interaction = interaction
        if freeze_encoder:
            self.interaction.requires_grad_(False)
            self.interaction.eval()
        self.freeze_encoder = freeze_encoder

    def train(self, mode: bool = True) -> "EncodedAnchorDelta":
        super().train(mode)
        if self.freeze_encoder:
            self.interaction.eval()
        return self

    def encode(self, protein_tokens: torch.Tensor, ligands: torch.Tensor) -> torch.Tensor:
        if self.freeze_encoder:
            with torch.no_grad():
                protein = self.interaction.encodeprotein(protein_tokens)
                return self.interaction.pairfromprotein(protein, ligands)
        protein = self.interaction.encodeprotein(protein_tokens)
        return self.interaction.pairfromprotein(protein, ligands)

    def forward(
        self,
        protein_tokens: torch.Tensor,
        query_ligand: torch.Tensor,
        anchor_ligand: torch.Tensor,
    ) -> torch.Tensor:
        query = self.encode(protein_tokens, query_ligand)
        anchor = self.encode(protein_tokens, anchor_ligand)
        return self.matrix(query, anchor)


def anchordeltaloss(
    prediction: torch.Tensor,
    target: torch.Tensor,
    *,
    order_weight: float = 0.25,
    order_temperature: float = 1.0,
    noise_floor: float = 0.0,
) -> dict[str, torch.Tensor]:
    """Huber difference loss plus optional order loss with noise-floor masking."""

    prediction, target = _validate_features(prediction.reshape(-1, 1), target.reshape(-1, 1))
    prediction = prediction[:, 0]
    target = target[:, 0]
    if order_weight < 0 or order_temperature <= 0 or noise_floor < 0:
        raise ValueError("loss weights and thresholds must be non-negative")
    regression = F.huber_loss(prediction, target)
    usable = target.abs() > noise_floor
    if bool(usable.any()):
        labels = (target[usable] > 0).float()
        ordering = F.binary_cross_entropy_with_logits(
            prediction[usable] / order_temperature, labels
        )
    else:
        ordering = prediction.sum() * 0.0
    return {"loss": regression + order_weight * ordering, "difference": regression, "order": ordering}


def aggregateanchors(
    anchor_values: torch.Tensor,
    *,
    method: Literal["uniform", "median", "weighted"] = "uniform",
    weights: torch.Tensor | None = None,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Aggregate ``[queries, anchors]`` absolute predictions label-free.

    Weighted aggregation accepts only externally supplied, label-free weights;
    support labels never enter this function.  It returns the estimate and the
    across-anchor disagreement variance.
    """

    values = anchor_values.float()
    if values.ndim != 2 or values.shape[1] == 0 or not bool(torch.isfinite(values).all()):
        raise ValueError("anchor_values must be finite [queries, anchors]")
    if method == "uniform":
        estimate = values.mean(dim=1)
        used = values.new_full(values.shape, 1.0 / values.shape[1])
    elif method == "median":
        estimate = values.median(dim=1).values
        used = values.new_full(values.shape, 1.0 / values.shape[1])
    elif method == "weighted":
        if weights is None or weights.shape != values.shape:
            raise ValueError("weighted aggregation requires matching label-free weights")
        weights = weights.float()
        if bool((weights < 0).any()) or not bool(torch.isfinite(weights).all()):
            raise ValueError("weights must be finite and non-negative")
        normalizer = weights.sum(dim=1, keepdim=True)
        if bool((normalizer <= 0).any()):
            raise ValueError("each query needs positive aggregate weight")
        used = weights / normalizer
        estimate = (used * values).sum(dim=1)
    else:
        raise ValueError("unknown anchor aggregation method")
    variance = (used * (values - estimate[:, None]).square()).sum(dim=1)
    return estimate, variance


def anchorabsolute(
    support_labels: torch.Tensor,
    deltas: torch.Tensor,
    *,
    method: Literal["uniform", "median", "weighted"] = "uniform",
    weights: torch.Tensor | None = None,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Convert pairwise deltas to absolute query predictions via support anchors."""

    labels = support_labels.to(device=deltas.device, dtype=deltas.dtype).reshape(-1)
    if deltas.ndim != 2 or deltas.shape[1] != labels.numel():
        raise ValueError("deltas must have shape [queries, support] and match labels")
    return aggregateanchors(labels[None, :] + deltas, method=method, weights=weights)


__all__ = [
    "AnchorDelta",
    "ComparatorHead",
    "EncodedAnchorDelta",
    "anchorabsolute",
    "aggregateanchors",
    "anchordeltaloss",
]
