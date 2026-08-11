"""Support-identifiable low-rank adaptation for unseen-target regression."""
from __future__ import annotations

from dataclasses import dataclass

import torch
from torch import Tensor, nn


@dataclass(frozen=True)
class SectionState:
    coefficients: Tensor
    projector: Tensor
    adaptation_map: Tensor
    covariance: Tensor | None
    rank: int
    condition_number: float


class IdentifiableMetaSection(nn.Module):
    """Learn a shared basis while restricting task adaptation to support span."""

    def __init__(self, input_dim: int, section_dim: int, ridge: float = 0.1,
                 eps: float = 1e-12) -> None:
        super().__init__()
        if input_dim < 1 or not 1 <= section_dim <= 5:
            raise ValueError("input_dim must be positive and section_dim must be in [1, 5]")
        if section_dim > input_dim:
            raise ValueError("section_dim cannot exceed input_dim")
        if ridge <= 0:
            raise ValueError("ridge must be strictly positive")
        self.input_dim = input_dim
        self.section_dim = section_dim
        self.ridge = float(ridge)
        self.eps = float(eps)
        initial = torch.randn(input_dim, section_dim, dtype=torch.float64)
        self.raw_basis = nn.Parameter(initial)

    def basis(self) -> Tensor:
        return torch.linalg.qr(self.raw_basis, mode="reduced").Q

    def coordinates(self, phi: Tensor) -> Tensor:
        if phi.shape[-1] != self.input_dim:
            raise ValueError("phi has the wrong final dimension")
        return phi @ self.basis()

    def adapt(self, support_phi: Tensor, support_residual: Tensor,
              measurement_covariance: Tensor | None = None) -> SectionState:
        if support_phi.ndim != 2 or support_residual.ndim != 1:
            raise ValueError("support_phi must be [k,D] and residual must be [k]")
        if support_phi.shape[0] != support_residual.shape[0] or support_phi.shape[0] < 1:
            raise ValueError("support features and residuals must have the same nonzero k")
        matrix = self.coordinates(support_phi)
        gram = matrix @ matrix.T
        identity = torch.eye(len(matrix), dtype=matrix.dtype, device=matrix.device)
        inverse_residual = torch.linalg.solve(
            gram + self.ridge * identity, support_residual
        )
        coefficients = matrix.T @ inverse_residual
        adaptation_map = matrix.T @ torch.linalg.solve(
            gram + self.ridge * identity, identity
        )
        projector = matrix.T @ torch.linalg.pinv(gram) @ matrix

        singular = torch.linalg.svdvals(matrix)
        tolerance = max(matrix.shape) * torch.finfo(matrix.dtype).eps * singular.max()
        retained = singular[singular > tolerance]
        rank = int(retained.numel())
        condition = float((retained.max() / retained.min()).detach()) if rank else float("inf")

        covariance = None
        if measurement_covariance is not None:
            expected = (len(matrix), len(matrix))
            if measurement_covariance.shape != expected:
                raise ValueError(f"measurement covariance must have shape {expected}")
            covariance = adaptation_map @ measurement_covariance @ adaptation_map.T
        return SectionState(
            coefficients=coefficients,
            projector=projector,
            adaptation_map=adaptation_map,
            covariance=covariance,
            rank=rank,
            condition_number=condition,
        )

    def query(self, query_phi: Tensor, state: SectionState) -> tuple[Tensor, Tensor, Tensor | None]:
        coordinates = self.coordinates(query_phi)
        observable = coordinates @ state.projector
        correction = observable @ state.coefficients
        coverage = observable.square().sum(dim=-1) / (
            coordinates.square().sum(dim=-1) + self.eps
        )
        uncertainty = None
        if state.covariance is not None:
            uncertainty = torch.sqrt(torch.clamp(
                torch.einsum("...i,ij,...j->...", coordinates, state.covariance, coordinates),
                min=0.0,
            ))
        return correction, coverage, uncertainty

    def bounded_noise_radius(self, query_phi: Tensor, state: SectionState,
                             support_error_bound: Tensor) -> Tensor:
        """Worst-case correction change for elementwise bounded support errors."""
        if support_error_bound.ndim != 1 or support_error_bound.shape[0] != state.adaptation_map.shape[1]:
            raise ValueError("support_error_bound must have one value per support row")
        coordinates = self.coordinates(query_phi)
        leverage = coordinates @ state.adaptation_map
        return leverage.abs() @ support_error_bound
