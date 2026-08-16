"""Endpoint- and source-aware observation likelihoods."""
from __future__ import annotations

import math

import torch
from torch import nn
from torch.nn import functional as F


def _inverse_softplus(value: float) -> float:
    return math.log(math.expm1(value))


def _validate_rows(*values: torch.Tensor) -> None:
    rows = {value.shape[0] for value in values}
    if len(rows) != 1:
        raise ValueError("all observation tensors must have matching rows")


class ObservationHeads(nn.Module):
    """Endpoint- and source-aware likelihoods applied after latent scoring."""

    def __init__(
        self,
        endpoints: int,
        sources: int,
        ordinal_classes: int = 4,
        minimum_scale: float = 1e-4,
    ) -> None:
        super().__init__()
        if endpoints < 1 or sources < 1:
            raise ValueError("endpoints and sources must be positive")
        if ordinal_classes < 3:
            raise ValueError("ordinal_classes must be at least three")
        if minimum_scale <= 0:
            raise ValueError("minimum_scale must be positive")

        self.endpoints = endpoints
        self.sources = sources
        self.ordinal_classes = ordinal_classes
        self.minimum_scale = minimum_scale

        self.exact_endpoint_bias = nn.Parameter(torch.zeros(endpoints))
        self.exact_source_bias = nn.Parameter(torch.zeros(endpoints, sources))
        initial_raw_scale = _inverse_softplus(0.5 - minimum_scale)
        self.exact_raw_scale = nn.Parameter(
            torch.full((endpoints,), initial_raw_scale)
        )

        self.binary_bias = nn.Parameter(torch.zeros(()))
        self.binary_source_bias = nn.Parameter(torch.zeros(sources))

        self.ordinal_source_bias = nn.Parameter(torch.zeros(sources))
        self.ordinal_first_cut = nn.Parameter(torch.tensor(-0.75))
        initial_gap = _inverse_softplus(0.75)
        self.ordinal_raw_gaps = nn.Parameter(
            torch.full((ordinal_classes - 2,), initial_gap)
        )

    def latent(self, score: torch.Tensor) -> torch.Tensor:
        """Return the deployable score without endpoint or source shortcuts."""

        if score.ndim != 1:
            raise ValueError("latent scores must be a vector")
        return score

    def neutral_parameters(
        self,
        score: torch.Tensor,
        endpoint: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """Observation parameters with endpoint calibration but neutral source."""

        endpoint = endpoint.to(dtype=torch.long)
        _validate_rows(score, endpoint)
        if endpoint.ndim != 1:
            raise ValueError("endpoint indices must be a vector")
        if bool((endpoint < 0).any()) or bool((endpoint >= self.endpoints).any()):
            raise ValueError("endpoint index is out of range")
        mean = score + self.exact_endpoint_bias[endpoint]
        scale = F.softplus(self.exact_raw_scale[endpoint]) + self.minimum_scale
        return mean, scale

    def _indices(
        self,
        endpoint: torch.Tensor,
        source: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        endpoint = endpoint.to(dtype=torch.long)
        source = source.to(dtype=torch.long)
        _validate_rows(endpoint, source)
        if endpoint.ndim != 1 or source.ndim != 1:
            raise ValueError("endpoint and source indices must be vectors")
        if bool((endpoint < 0).any()) or bool((endpoint >= self.endpoints).any()):
            raise ValueError("endpoint index is out of range")
        if bool((source < 0).any()) or bool((source >= self.sources).any()):
            raise ValueError("source index is out of range")
        return endpoint, source

    def exact_parameters(
        self,
        score: torch.Tensor,
        endpoint: torch.Tensor,
        source: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        endpoint, source = self._indices(endpoint, source)
        _validate_rows(score, endpoint)
        mean = (
            score
            + self.exact_endpoint_bias[endpoint]
            + self.exact_source_bias[endpoint, source]
        )
        scale = F.softplus(self.exact_raw_scale[endpoint]) + self.minimum_scale
        return mean, scale

    def exact_nll(
        self,
        score: torch.Tensor,
        value: torch.Tensor,
        endpoint: torch.Tensor,
        source: torch.Tensor,
    ) -> torch.Tensor:
        _validate_rows(score, value)
        mean, scale = self.exact_parameters(score, endpoint, source)
        standardized = (value - mean) / scale
        return (
            0.5 * standardized.square()
            + torch.log(scale)
            + 0.5 * math.log(2.0 * math.pi)
        ).mean()

    def censored_nll(
        self,
        score: torch.Tensor,
        bound: torch.Tensor,
        direction: torch.Tensor,
        endpoint: torch.Tensor,
        source: torch.Tensor,
    ) -> torch.Tensor:
        """Stable one-sided Gaussian likelihood.

        ``direction`` is +1 for ``Y >= bound`` and -1 for ``Y <= bound``.
        """
        _validate_rows(score, bound, direction)
        if not bool(torch.all((direction == 1) | (direction == -1))):
            raise ValueError("censor direction must contain only -1 or +1")
        mean, scale = self.exact_parameters(score, endpoint, source)
        signed_tail = direction.to(score.dtype) * (mean - bound) / scale
        return -torch.special.log_ndtr(signed_tail).mean()

    def binary_logits(
        self,
        score: torch.Tensor,
        source: torch.Tensor,
    ) -> torch.Tensor:
        source = source.to(dtype=torch.long)
        _validate_rows(score, source)
        if source.ndim != 1:
            raise ValueError("source indices must be a vector")
        if bool((source < 0).any()) or bool((source >= self.sources).any()):
            raise ValueError("source index is out of range")
        return score + self.binary_bias + self.binary_source_bias[source]

    def binary_nll(
        self,
        score: torch.Tensor,
        label: torch.Tensor,
        source: torch.Tensor,
    ) -> torch.Tensor:
        _validate_rows(score, label)
        logits = self.binary_logits(score, source)
        return F.binary_cross_entropy_with_logits(logits, label.to(score.dtype))

    def ordinal_cutpoints(self) -> torch.Tensor:
        gaps = F.softplus(self.ordinal_raw_gaps) + 1e-4
        return torch.cat(
            (
                self.ordinal_first_cut.reshape(1),
                self.ordinal_first_cut + torch.cumsum(gaps, dim=0),
            )
        )

    def ordinal_log_probabilities(
        self,
        score: torch.Tensor,
        source: torch.Tensor,
    ) -> torch.Tensor:
        source = source.to(dtype=torch.long)
        _validate_rows(score, source)
        if source.ndim != 1:
            raise ValueError("source indices must be a vector")
        if bool((source < 0).any()) or bool((source >= self.sources).any()):
            raise ValueError("source index is out of range")

        eta = score + self.ordinal_source_bias[source]
        cut = self.ordinal_cutpoints()
        cumulative = torch.sigmoid(cut.unsqueeze(0) - eta.unsqueeze(1))
        probabilities = torch.cat(
            (
                cumulative[:, :1],
                cumulative[:, 1:] - cumulative[:, :-1],
                1.0 - cumulative[:, -1:],
            ),
            dim=1,
        )
        return torch.log(probabilities.clamp_min(torch.finfo(score.dtype).tiny))

    def ordinal_nll(
        self,
        score: torch.Tensor,
        label: torch.Tensor,
        source: torch.Tensor,
    ) -> torch.Tensor:
        label = label.to(dtype=torch.long)
        _validate_rows(score, label)
        if bool((label < 0).any()) or bool((label >= self.ordinal_classes).any()):
            raise ValueError("ordinal label is out of range")
        log_probabilities = self.ordinal_log_probabilities(score, source)
        return -log_probabilities.gather(1, label.unsqueeze(1)).mean()

__all__ = ["ObservationHeads"]
