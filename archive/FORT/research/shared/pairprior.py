"""Sequence-conditioned global affinity prior for a zero-support research gate."""

from __future__ import annotations

import math

import torch
from torch import nn
from torch.nn import functional as F


class SegmentBlock(nn.Module):
    """Pre-norm attention block over frozen protein sequence segments."""

    def __init__(self, width: int, heads: int) -> None:
        super().__init__()
        self.norm = nn.LayerNorm(width)
        self.attention = nn.MultiheadAttention(
            width, heads, dropout=0.0, batch_first=True
        )
        self.feednorm = nn.LayerNorm(width)
        self.feed = nn.Sequential(
            nn.Linear(width, 4 * width),
            nn.GELU(),
            nn.Linear(4 * width, width),
        )

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        normalized = self.norm(values)
        attended, _ = self.attention(
            normalized, normalized, normalized, need_weights=False
        )
        values = values + attended
        return values + self.feed(self.feednorm(values))


class PairPrior(nn.Module):
    """Predict a Gaussian affinity prior from protein segments and a ligand."""

    def __init__(
        self,
        proteindim: int = 1280,
        liganddim: int = 1034,
        width: int = 192,
        heads: int = 8,
        blocks: int = 2,
    ) -> None:
        super().__init__()
        if width % heads:
            raise ValueError("width must be divisible by attention heads")
        self.proteindim = proteindim
        self.liganddim = liganddim
        self.proteinprojection = nn.Linear(proteindim, width)
        self.blocks = nn.ModuleList(SegmentBlock(width, heads) for _ in range(blocks))
        self.poolscore = nn.Linear(width, 1, bias=False)
        self.proteinnorm = nn.LayerNorm(width)
        self.ligand = nn.Sequential(
            nn.Linear(liganddim, 2 * width),
            nn.GELU(),
            nn.LayerNorm(2 * width),
            nn.Linear(2 * width, width),
            nn.GELU(),
        )
        self.film = nn.Linear(width, 2 * width)
        self.head = nn.Sequential(
            nn.Linear(5 * width, 2 * width),
            nn.GELU(),
            nn.LayerNorm(2 * width),
            nn.Linear(2 * width, width),
            nn.GELU(),
            nn.Linear(width, 2),
        )

    def encodeprotein(self, segments: torch.Tensor) -> torch.Tensor:
        if segments.ndim != 3 or segments.shape[-1] != self.proteindim:
            raise ValueError(
                "protein segments must have shape [batch, segments, protein_dim]"
            )
        values = self.proteinprojection(segments.float())
        for block in self.blocks:
            values = block(values)
        weight = torch.softmax(self.poolscore(values), dim=1)
        return self.proteinnorm((weight * values).sum(dim=1))

    def forward(
        self, segments: torch.Tensor, ligand: torch.Tensor
    ) -> dict[str, torch.Tensor]:
        if ligand.ndim != 2 or ligand.shape[-1] != self.liganddim:
            raise ValueError("ligands must have shape [batch, ligand_dim]")
        if segments.shape[0] != ligand.shape[0]:
            raise ValueError("protein and ligand batches must align")
        protein = self.encodeprotein(segments)
        chemical = self.ligand(ligand.float())
        scale, shift = self.film(protein).chunk(2, dim=1)
        conditioned = chemical * (1.0 + torch.tanh(scale)) + shift
        pair = torch.cat(
            (
                chemical,
                protein,
                chemical * protein,
                (chemical - protein).abs(),
                conditioned,
            ),
            dim=1,
        )
        raw = self.head(pair)
        return {
            "prediction": raw[:, 0],
            "logvariance": raw[:, 1].clamp(-4.0, 4.0),
            "protein": protein,
        }

    def predict(
        self,
        segments: torch.Tensor,
        ligand: torch.Tensor,
        supportlabel: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        """Return the k=0 prior; nonempty support is outside this research gate."""

        if supportlabel is not None and supportlabel.numel() != 0:
            raise ValueError("the global pair prior accepts k=0 only")
        return self(segments, ligand)


def gaussiannll(
    output: dict[str, torch.Tensor], label: torch.Tensor
) -> torch.Tensor:
    """Per-row Gaussian negative log likelihood including its constant."""

    prediction = output["prediction"]
    logvariance = output["logvariance"]
    label = label.float().reshape(-1)
    if prediction.shape != label.shape or logvariance.shape != label.shape:
        raise ValueError("prediction, variance, and label must align")
    error = prediction - label
    return 0.5 * (
        math.log(2.0 * math.pi) + logvariance + error.square() * torch.exp(-logvariance)
    )


def evidencecontrast(
    correct: dict[str, torch.Tensor],
    wrong: dict[str, torch.Tensor],
    label: torch.Tensor,
    margin: float = 0.05,
    weight: torch.Tensor | None = None,
) -> torch.Tensor:
    """Prefer correct-protein evidence without treating unknown pairs as negatives."""

    correctnll = gaussiannll(correct, label)
    wrongnll = gaussiannll(wrong, label)
    value = F.softplus(margin + correctnll - wrongnll)
    if weight is None:
        return value.mean()
    weight = weight.to(device=value.device, dtype=value.dtype).reshape(-1)
    if weight.shape != value.shape or not bool((weight > 0).all()):
        raise ValueError("contrast weights must be positive and aligned")
    return (value * weight).sum() / weight.sum()


__all__ = ["PairPrior", "SegmentBlock", "evidencecontrast", "gaussiannll"]
