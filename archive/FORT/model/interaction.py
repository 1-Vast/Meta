"""Low-rank protein-ligand interaction feature encoder."""

from __future__ import annotations

import math

import torch
from torch import nn

from .protein import ProteinStage


class InteractionEncoder(nn.Module):
    """Encode one target and a batch of ligands into fixed-width pair features."""

    def __init__(
        self,
        proteindim: int,
        liganddim: int,
        dmodel: int = 256,
        taskdim: int = 8,
        stages: int = 2,
        landmarks: int = 32,
        backbone: str = "hybrid",
        conditioned: bool = True,
        interactiononly: bool = False,
    ) -> None:
        super().__init__()
        if backbone not in {"transformer", "mamba", "hybrid"}:
            raise ValueError(f"unsupported backbone: {backbone}")
        self.proteindim = proteindim
        self.liganddim = liganddim
        self.conditioned = conditioned
        self.interactiononly = interactiononly
        self.proteinprojection = nn.Linear(
            proteindim, dmodel, bias=not interactiononly
        )
        self.ligandprojection = nn.Linear(
            liganddim, dmodel, bias=not interactiononly
        )
        self.stages = nn.ModuleList(
            ProteinStage(dmodel, taskdim, landmarks, backbone, conditioned)
            for _ in range(stages)
        )
        self.poolkey = nn.Linear(dmodel, dmodel, bias=False)
        self.poolvalue = nn.Linear(dmodel, dmodel, bias=False)
        pairwidth = dmodel if interactiononly else 4 * dmodel
        self.pairfeature = nn.Sequential(
            nn.Linear(pairwidth, dmodel, bias=not interactiononly),
            nn.GELU(),
            nn.Linear(dmodel, taskdim, bias=not interactiononly),
        )
        self.pairnorm = nn.LayerNorm(taskdim, elementwise_affine=False)

    def encodeprotein(
        self, tokens: torch.Tensor, code: torch.Tensor | None = None
    ) -> torch.Tensor:
        if tokens.ndim != 2 or tokens.shape[1] != self.proteindim:
            raise ValueError("protein tokens must have shape [residues, protein_dim]")
        values = self.proteinprojection(tokens.float()).unsqueeze(0)
        for stage in self.stages:
            values = stage(values, code)
        return values.squeeze(0)

    def pairfromprotein(
        self, protein: torch.Tensor, ligand: torch.Tensor
    ) -> torch.Tensor:
        if ligand.ndim != 2 or ligand.shape[1] != self.liganddim:
            raise ValueError("ligand features must have shape [rows, ligand_dim]")
        ligandhidden = self.ligandprojection(ligand.float())
        logits = ligandhidden @ self.poolkey(protein).T / math.sqrt(protein.shape[1])
        pooled = torch.softmax(logits, dim=1) @ self.poolvalue(protein)
        if self.interactiononly:
            combined = pooled * ligandhidden
        else:
            combined = torch.cat(
                (
                    pooled,
                    ligandhidden,
                    pooled * ligandhidden,
                    (pooled - ligandhidden).abs(),
                ),
                dim=1,
            )
        return self.pairnorm(self.pairfeature(combined))

    def forward(
        self,
        tokens: torch.Tensor,
        ligand: torch.Tensor,
        code: torch.Tensor | None = None,
    ) -> torch.Tensor:
        protein = self.encodeprotein(tokens, code)
        return self.pairfromprotein(protein, ligand)


__all__ = ["InteractionEncoder"]
