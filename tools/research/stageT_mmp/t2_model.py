"""The T2 discriminator: a small, ordinary gradient-trained response model.

    R(tau, p)                      per-(transformation, protein) response
    D_hat(tau, p1, p2) = R(tau,p1) - R(tau,p2)

Because the prediction is a difference of a per-protein response, identity
(`p1 = p2` gives exactly 0), antisymmetry and protein-cycle consistency hold for
every parameter setting, before and after training. They are properties of the
construction, not of the fit.

Deliberately small. The inputs are a **frozen** protein embedding and a
**structured, hand-specified** transformation descriptor -- no learned ligand
encoder, no whole-molecule representation, so this cannot quietly become another
DTA model. Forbidden as input and absent from the signature: target ID,
document, assay, panel, target index, component ID.
"""
from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn as nn


@dataclass(frozen=True)
class DiscriminatorConfig:
    descriptor_dim: int
    protein_dim: int
    width: int = 64
    hidden: int = 128
    depth: int = 2
    # "zero"            -> R is a learned constant; D_hat is identically 0
    # "transformation"  -> R = f(tau); D_hat is identically 0 (mu_tau cancels)
    # "protein"         -> R = f(tau, p); the only arm that can express D
    mode: str = "protein"


class Response(nn.Module):
    """R(tau, p)."""

    def __init__(self, config: DiscriminatorConfig) -> None:
        super().__init__()
        self.config = config
        self.mode = config.mode
        if self.mode not in {"zero", "transformation", "protein"}:
            raise ValueError(f"unknown mode: {self.mode}")
        self.constant = nn.Parameter(torch.zeros(()))
        self.transformation = nn.Sequential(
            nn.Linear(config.descriptor_dim, config.width), nn.SiLU(),
            nn.Linear(config.width, config.width))
        self.protein = nn.Sequential(
            nn.Linear(config.protein_dim, config.width), nn.SiLU(),
            nn.Linear(config.width, config.width))
        layers: list[nn.Module] = [nn.LayerNorm(2 * config.width),
                                   nn.Linear(2 * config.width, config.hidden),
                                   nn.SiLU()]
        for _ in range(config.depth - 1):
            layers += [nn.Linear(config.hidden, config.hidden), nn.SiLU()]
        # No output bias: an additive constant on R cancels exactly in
        # R(tau,p1) - R(tau,p2), so it can never receive gradient and can never
        # change a prediction. Keeping it would inflate the parameter count with
        # a provably dead weight; the no-dead-parameter test pins its absence.
        layers.append(nn.Linear(config.hidden, 1, bias=False))
        self.head = nn.Sequential(*layers)

    def forward(self, descriptor: torch.Tensor,
                protein: torch.Tensor) -> torch.Tensor:
        if self.mode == "zero":
            return self.constant.expand(descriptor.shape[0])
        edit = self.transformation(descriptor)
        if self.mode == "transformation":
            # Protein branch is deliberately not consulted. D_hat is then
            # identically zero, which is the point: mu_tau cancels in the
            # double difference and a protein-free model cannot express D.
            state = torch.cat((edit, torch.zeros_like(edit)), dim=-1)
        else:
            state = torch.cat((edit, self.protein(protein)), dim=-1)
        return self.head(state).squeeze(-1)


class DoubleDifferenceModel(nn.Module):
    """D_hat(tau, p1, p2) = R(tau, p1) - R(tau, p2)."""

    def __init__(self, config: DiscriminatorConfig) -> None:
        super().__init__()
        self.config = config
        self.response = Response(config)

    def forward(self, descriptor: torch.Tensor, protein_left: torch.Tensor,
                protein_right: torch.Tensor) -> torch.Tensor:
        return (self.response(descriptor, protein_left)
                - self.response(descriptor, protein_right))


def parameter_report(model: DoubleDifferenceModel) -> dict:
    return {
        "total": sum(p.numel() for p in model.parameters()),
        "transformation_branch": sum(
            p.numel() for p in model.response.transformation.parameters()),
        "protein_branch": sum(
            p.numel() for p in model.response.protein.parameters()),
        "head": sum(p.numel() for p in model.response.head.parameters()),
        "mode": model.config.mode,
    }
