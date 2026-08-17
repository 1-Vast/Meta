"""The cross-target protein-conditioned SAR field.

The construction under test:

    phi(L)                 target-independent shared chemical coordinate.
                           No protein tensor is on its input path.
    u_ab = phi(L_b)-phi(L_a)   the directed chemical change (protein-free).
    alpha(P)               protein response coefficients.
    V(P,L)                 a nonlinear SAR potential, FiLM-modulated by alpha.
    dy_hat(P,a,b) = V(P,L_b) - V(P,L_a)

Why the prediction is a *potential difference* and not a learned pair function:
a difference of a scalar field is a conservative (curl-free) field, so

    dy_hat(P,a,b) = -dy_hat(P,b,a)                 antisymmetry
    dy_hat(P,a,a) = 0                              identity
    dy_hat_ab + dy_hat_bc + dy_hat_ca = 0          cycle consistency

hold identically in exact arithmetic, for every parameter setting, before and
after training.  They are properties of the construction, not of the fit.  In
IEEE-754 the first two are bitwise exact (negation and self-subtraction are
exact operations); the third telescopes to a rounding residual at machine
epsilon, which `tests/test_structural.py` measures rather than assumes.

Two constructions are deliberately *not* used:

* the proposed explicit quadratic `e^T H e` term.  It is even under direction
  reversal (`(-e)^T H (-e) = e^T H e`), so adding it to a signed prediction
  destroys antisymmetry.  There is no such term anywhere in this module;
* double protein conditioning.  `phi` never sees a protein and `alpha` never
  sees a ligand; the two meet exactly once, inside the potential.

No 3D coordinate, pose, docking or protein-ligand contact is constructed or
claimed here.  `phi` is a 2D-graph + fingerprint coordinate and `alpha` is a
frozen-PLM sequence readout.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

import torch
import torch.nn as nn

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from model.encoders import LigandEncoder


@dataclass(frozen=True)
class FieldConfig:
    hidden: int = 96
    graph_layers: int = 3
    fingerprint_bits: int = 1024
    coordinate: int = 64
    response: int = 64
    potential_width: int = 128
    potential_depth: int = 3
    protein_dim: int = 640
    protein_conditioned: bool = True
    dtype: torch.dtype = torch.float32


class LigandCoordinate(nn.Module):
    """phi(L). Protein-free by construction: `forward` has no protein argument."""

    def __init__(self, config: FieldConfig) -> None:
        super().__init__()
        self.graph = LigandEncoder(
            config.hidden, n_layers=config.graph_layers, dtype=config.dtype)
        self.fingerprint = nn.Sequential(
            nn.Linear(config.fingerprint_bits, config.hidden, dtype=config.dtype),
            nn.SiLU(),
            nn.Linear(config.hidden, config.hidden, dtype=config.dtype))
        self.mix = nn.Sequential(
            nn.LayerNorm(2 * config.hidden, dtype=config.dtype),
            nn.Linear(2 * config.hidden, config.coordinate, dtype=config.dtype))

    def forward(self, atoms: torch.Tensor, bonds: torch.Tensor,
                mask: torch.Tensor, fingerprint: torch.Tensor) -> torch.Tensor:
        pooled, _states = self.graph(atoms, bonds, mask)
        return self.mix(torch.cat((pooled, self.fingerprint(fingerprint)), dim=-1))


class ProteinResponse(nn.Module):
    """alpha(P). Ligand-free by construction: `forward` has no ligand argument."""

    def __init__(self, config: FieldConfig) -> None:
        super().__init__()
        self.pooled = nn.Sequential(
            nn.Linear(config.protein_dim, config.hidden, dtype=config.dtype),
            nn.SiLU(),
            nn.Linear(config.hidden, config.hidden, dtype=config.dtype))
        self.residue = nn.Linear(config.protein_dim, config.hidden,
                                 dtype=config.dtype)
        self.head = nn.Sequential(
            nn.LayerNorm(2 * config.hidden, dtype=config.dtype),
            nn.Linear(2 * config.hidden, config.hidden, dtype=config.dtype),
            nn.SiLU(),
            nn.Linear(config.hidden, config.response, dtype=config.dtype))

    def forward(self, pooled: torch.Tensor, residues: torch.Tensor,
                mask: torch.Tensor) -> torch.Tensor:
        gate = mask.to(residues.dtype).unsqueeze(-1)
        summary = (self.residue(residues) * gate).sum(1) / gate.sum(1).clamp_min(1.0)
        return self.head(torch.cat((self.pooled(pooled), summary), dim=-1))


class SARPotential(nn.Module):
    """V(P,L): a nonlinear scalar potential over phi, FiLM-modulated by alpha.

    The FiLM projections start small (std 0.02) rather than at exactly zero, so
    a gradient reaches `alpha` at the very first step; an exact-zero init makes
    dV/dalpha identically zero at step 0 and the protein-path gradient test
    would then be trivially satisfied only after the first update.
    """

    def __init__(self, config: FieldConfig) -> None:
        super().__init__()
        width, depth = config.potential_width, config.potential_depth
        self.inp = nn.Linear(config.coordinate, width, dtype=config.dtype)
        self.blocks = nn.ModuleList(
            nn.Linear(width, width, dtype=config.dtype) for _ in range(depth))
        self.scale = nn.ModuleList(
            nn.Linear(config.response, width, dtype=config.dtype)
            for _ in range(depth))
        self.shift = nn.ModuleList(
            nn.Linear(config.response, width, dtype=config.dtype)
            for _ in range(depth))
        self.out = nn.Linear(width, 1, dtype=config.dtype)
        for module in (*self.scale, *self.shift):
            nn.init.normal_(module.weight, std=0.02)
            nn.init.zeros_(module.bias)

    def forward(self, coordinate: torch.Tensor,
                response: torch.Tensor) -> torch.Tensor:
        state = torch.nn.functional.silu(self.inp(coordinate))
        for block, scale, shift in zip(self.blocks, self.scale, self.shift):
            modulated = block(state) * (1.0 + scale(response)) + shift(response)
            state = state + torch.nn.functional.silu(modulated)
        return self.out(state).squeeze(-1)


class SARField(nn.Module):
    """The full field. `protein_conditioned=False` is arm A (ligand-only)."""

    def __init__(self, config: FieldConfig) -> None:
        super().__init__()
        self.config = config
        self.coordinate = LigandCoordinate(config)
        self.potential_module = SARPotential(config)
        self.protein_conditioned = bool(config.protein_conditioned)
        if self.protein_conditioned:
            self.response_module = ProteinResponse(config)
        else:
            # A learned constant response: the same potential, the same FiLM
            # path, the same gradient route -- with the protein identity
            # removed and nothing else changed.
            self.constant_response = nn.Parameter(
                torch.zeros(config.response, dtype=config.dtype))

    # -- the three parts, exposed separately so the tests can isolate them ---

    def phi(self, atoms, bonds, mask, fingerprint) -> torch.Tensor:
        return self.coordinate(atoms, bonds, mask, fingerprint)

    def alpha(self, pooled, residues, protein_mask) -> torch.Tensor:
        if self.protein_conditioned:
            return self.response_module(pooled, residues, protein_mask)
        return self.constant_response.unsqueeze(0).expand(pooled.shape[0], -1)

    @staticmethod
    def direction(phi_a: torch.Tensor, phi_b: torch.Tensor) -> torch.Tensor:
        """u_ab = phi(L_b) - phi(L_a). Protein-free, exactly antisymmetric."""
        return phi_b - phi_a

    def potential(self, coordinate: torch.Tensor,
                  response: torch.Tensor) -> torch.Tensor:
        return self.potential_module(coordinate, response)

    def forward(self, phi_a: torch.Tensor, phi_b: torch.Tensor,
                response: torch.Tensor) -> torch.Tensor:
        """dy_hat = V(P, L_b) - V(P, L_a)."""
        return (self.potential_module(phi_b, response)
                - self.potential_module(phi_a, response))


def build_field(config: FieldConfig) -> SARField:
    return SARField(config)


def parameter_report(field: SARField) -> dict:
    """Parameter counts by path, so arm budgets can be compared honestly."""
    def total(module: nn.Module) -> int:
        return sum(p.numel() for p in module.parameters())

    counts = {
        "ligand_coordinate": total(field.coordinate),
        "potential": total(field.potential_module),
        "protein_response": (total(field.response_module)
                             if field.protein_conditioned else 0),
        "constant_response": (0 if field.protein_conditioned
                              else field.constant_response.numel()),
    }
    counts["total"] = sum(p.numel() for p in field.parameters())
    return counts
