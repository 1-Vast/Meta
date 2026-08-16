"""Protein sequence blocks for frozen language-model residue tokens."""

from __future__ import annotations

import math

import torch
from torch import nn
from torch.nn import functional as F

try:
    from mamba_ssm.modules.mamba_simple import Mamba
except ImportError as error:  # pragma: no cover
    Mamba = None
    MAMBA_ERROR = error
else:
    MAMBA_ERROR = None


def makemamba(dmodel: int, state: int = 16, conv: int = 4, expand: int = 2) -> nn.Module:
    if Mamba is None:
        raise RuntimeError("mamba-ssm is required in the drug environment") from MAMBA_ERROR
    return Mamba(d_model=dmodel, d_state=state, d_conv=conv, expand=expand)


class BidirectionalMamba(nn.Module):
    """Position-aligned forward and reverse selective scans."""

    def __init__(self, dmodel: int, state: int = 16, conv: int = 4, expand: int = 2) -> None:
        super().__init__()
        self.forwardscan = makemamba(dmodel, state, conv, expand)
        self.backwardscan = makemamba(dmodel, state, conv, expand)
        self.output = nn.Linear(2 * dmodel, dmodel)

    def forward(self, tokens: torch.Tensor) -> torch.Tensor:
        forward = self.forwardscan(tokens)
        backward = torch.flip(self.backwardscan(torch.flip(tokens, dims=(1,))), dims=(1,))
        return self.output(torch.cat((forward, backward), dim=-1))


class LandmarkAttention(nn.Module):
    """Residue attention against fixed-count segment summaries."""

    def __init__(self, dmodel: int, landmarks: int) -> None:
        super().__init__()
        self.landmarks = landmarks
        self.attention = nn.MultiheadAttention(dmodel, num_heads=8, batch_first=True)
        self.norm = nn.LayerNorm(dmodel)

    def summarize(self, tokens: torch.Tensor) -> torch.Tensor:
        if tokens.ndim != 3 or tokens.shape[1] == 0:
            raise ValueError("protein tokens must have shape [batch, residues, width]")
        length = tokens.shape[1]
        stride = max(1, math.ceil(length / self.landmarks))
        padding = (-length) % stride
        padded = F.pad(tokens, (0, 0, 0, padding)) if padding else tokens
        mask = torch.ones(tokens.shape[:2], device=tokens.device, dtype=tokens.dtype)
        mask = F.pad(mask, (0, padding)) if padding else mask
        groups = padded.reshape(tokens.shape[0], -1, stride, tokens.shape[2])
        weights = mask.reshape(tokens.shape[0], -1, stride, 1)
        return (groups * weights).sum(dim=2) / weights.sum(dim=2).clamp_min(1.0)

    def forward(self, tokens: torch.Tensor) -> torch.Tensor:
        landmarks = self.summarize(tokens)
        attended, _ = self.attention(tokens, landmarks, landmarks, need_weights=False)
        return self.norm(tokens + attended)


class ProteinStage(nn.Module):
    """Sequence stage with optional support conditioning for the gradient baseline."""

    def __init__(
        self,
        dmodel: int,
        taskdim: int,
        landmarks: int,
        backbone: str,
        conditioned: bool = True,
    ) -> None:
        super().__init__()
        self.conditioned = conditioned
        self.scan = BidirectionalMamba(dmodel) if backbone != "transformer" else None
        self.attention = LandmarkAttention(dmodel, landmarks) if backbone != "mamba" else None
        self.film = nn.Linear(taskdim, 2 * dmodel) if conditioned else None
        self.recap = nn.Linear(taskdim, dmodel, bias=False) if conditioned else None
        self.norm = nn.LayerNorm(dmodel)

    def forward(self, tokens: torch.Tensor, code: torch.Tensor | None = None) -> torch.Tensor:
        values = tokens
        if self.scan is not None:
            values = values + self.scan(values)
        if self.attention is not None:
            values = self.attention(values)
        if self.conditioned:
            if code is None:
                raise ValueError("conditioned protein stages require a support code")
            positions = torch.arange(values.shape[1], device=values.device)
            boundary = (positions % 64 == 0).to(values.dtype)
            values = values + self.recap(code)[None, None, :] * boundary[None, :, None]
            scale, shift = self.film(code).chunk(2)
            values = values * (1.0 + 0.1 * torch.tanh(scale)[None, None, :])
            values = values + 0.1 * torch.tanh(shift)[None, None, :]
        elif code is not None:
            raise ValueError("Bayesian protein stages do not accept support-label codes")
        return self.norm(values)


__all__ = ["BidirectionalMamba", "LandmarkAttention", "ProteinStage"]
