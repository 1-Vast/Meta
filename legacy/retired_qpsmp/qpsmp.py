"""Legacy analytic centered-ridge diagnostic for flat-feature experiments.

The learnable biological meta-operator lives in :mod:`model.qpsmp_meta`. This
legacy module has independent heads and is not an identical-family comparator
or coverage certificate for that neural checkpoint.
"""
from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn as nn
from torch import Tensor


@dataclass(frozen=True)
class QPSMPOutput:
    ridge_prediction: Tensor
    section_midpoint: Tensor
    section_radius: Tensor
    center_radius: Tensor
    diagnostic_total_radius: Tensor
    task_state: Tensor
    level_shift: Tensor
    query_basis: Tensor
    support_baseline: Tensor


class QPSMPCore(nn.Module):
    """Analytic centered-ridge comparator over precomputed features."""

    def __init__(
            self, feature_dim: int, task_dim: int,
            ridge: float = 1.0, section_radius_bound: float = 1.0,
            baseline_hidden_dim: int = 0, baseline_dim: int | None = None,
            dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if feature_dim < 1 or task_dim < 0 or ridge <= 0 or section_radius_bound < 0:
            raise ValueError("invalid QPSMP dimensions, ridge, or section bound")
        if baseline_hidden_dim < 0:
            raise ValueError("baseline hidden dimension cannot be negative")
        self.feature_dim = int(feature_dim)
        self.baseline_dim = int(feature_dim if baseline_dim is None else baseline_dim)
        if self.baseline_dim < 1:
            raise ValueError("baseline feature dimension must be positive")
        self.task_dim = int(task_dim)
        self.ridge = float(ridge)
        self.section_radius_bound = float(section_radius_bound)
        self.baseline = (
            nn.Sequential(
                nn.Linear(self.baseline_dim, baseline_hidden_dim, dtype=dtype),
                nn.SiLU(),
                nn.Linear(baseline_hidden_dim, 1, dtype=dtype),
            )
            if baseline_hidden_dim else nn.Linear(self.baseline_dim, 1, dtype=dtype)
        )
        self.zero_shot = nn.Linear(feature_dim, 1, bias=False, dtype=dtype)
        self.section = (
            nn.Linear(feature_dim, task_dim, bias=False, dtype=dtype)
            if task_dim else None
        )

    def _check_features(self, value: Tensor, name: str) -> None:
        if value.ndim != 2 or value.shape[1] != self.feature_dim:
            raise ValueError(f"{name} must have shape [N,{self.feature_dim}]")

    def _check_baseline_features(self, value: Tensor, name: str) -> None:
        if value.ndim != 2 or value.shape[1] != self.baseline_dim:
            raise ValueError(f"{name} must have shape [N,{self.baseline_dim}]")

    def scalar_components(
            self, features: Tensor, baseline_features: Tensor | None = None) -> tuple[Tensor, Tensor]:
        self._check_features(features, "features")
        baseline_input = features if baseline_features is None else baseline_features
        self._check_baseline_features(baseline_input, "baseline_features")
        baseline = self.baseline(baseline_input).squeeze(-1)
        zero_shot = self.zero_shot(features).squeeze(-1)
        return baseline, zero_shot

    def basis(self, features: Tensor) -> Tensor:
        self._check_features(features, "features")
        if self.task_dim == 0:
            return torch.empty(
                features.shape[0], 0, device=features.device, dtype=features.dtype)
        assert self.section is not None
        return self.section(features)

    def centered_state(
            self, support_features: Tensor, support_y: Tensor,
            support_baseline_features: Tensor | None = None,
            *, adapt: bool = True) -> tuple[Tensor, Tensor, Tensor, Tensor]:
        self._check_features(support_features, "support_features")
        if support_y.ndim != 1 or support_y.shape[0] != support_features.shape[0]:
            raise ValueError("support_y must have one value per support feature")
        baseline, zero_shot = self.scalar_components(
            support_features, support_baseline_features)
        support_baseline = baseline + zero_shot
        zero_level = torch.zeros((), device=support_features.device, dtype=support_features.dtype)
        if not adapt or support_y.numel() == 0:
            state = torch.zeros(
                self.task_dim, device=support_features.device, dtype=support_features.dtype)
            return state, zero_level, support_baseline, self.basis(support_features)
        # This is the fixed T_level(S) channel. It is invariant to support order
        # and to ligand-label rebinding because it uses only separate means.
        level_shift = support_y.mean() - support_baseline.mean()
        if self.task_dim == 0:
            state = torch.empty(0, device=support_features.device, dtype=support_features.dtype)
            return state, level_shift, support_baseline, self.basis(support_features)
        phi = self.basis(support_features)
        residual = support_y - support_baseline
        k = support_y.shape[0]
        centered_phi = phi - phi.mean(dim=0, keepdim=True)
        centered_residual = residual - residual.mean()
        identity = torch.eye(self.task_dim, device=phi.device, dtype=phi.dtype)
        gram = (centered_phi.T @ centered_phi) / k + self.ridge * identity
        rhs = (centered_phi.T @ centered_residual) / k
        state = torch.linalg.solve(gram, rhs)
        return state, level_shift, support_baseline, phi

    def forward(
            self, support_features: Tensor, support_y: Tensor,
            query_features: Tensor, *,
            support_baseline_features: Tensor | None = None,
            query_baseline_features: Tensor | None = None,
            adapt: bool = True,
            level_shift_override: float | Tensor | None = None,
            task_state_override: Tensor | None = None,
            repr_radius: float | Tensor = 0.0,
            trans_radius: float | Tensor = 0.0,
            obs_radius: float | Tensor = 0.0) -> QPSMPOutput:
        state, level_shift, support_baseline, _ = self.centered_state(
            support_features, support_y, support_baseline_features, adapt=adapt)
        if level_shift_override is not None:
            level_shift = torch.as_tensor(
                level_shift_override,
                device=query_features.device,
                dtype=query_features.dtype)
        if task_state_override is not None:
            if task_state_override.shape != state.shape:
                raise ValueError("task_state_override has the wrong shape")
            state = task_state_override.to(device=state.device, dtype=state.dtype)
        query_baseline, query_zero_shot = self.scalar_components(
            query_features, query_baseline_features)
        query_base = query_baseline + query_zero_shot + level_shift
        query_phi = self.basis(query_features)
        if self.task_dim:
            ridge_prediction = query_base + query_phi @ state
        else:
            ridge_prediction = query_base
        section_contribution_midpoint, section_radius = self.point_section(
            query_phi, support_features)
        section_midpoint = query_base + section_contribution_midpoint
        center_radius = torch.abs(ridge_prediction - section_midpoint)
        diagnostic_total_radius = (
            center_radius + section_radius
            + torch.as_tensor(repr_radius, device=query_features.device, dtype=query_features.dtype)
            + torch.as_tensor(trans_radius, device=query_features.device, dtype=query_features.dtype)
            + torch.as_tensor(obs_radius, device=query_features.device, dtype=query_features.dtype)
        )
        return QPSMPOutput(
            ridge_prediction=ridge_prediction,
            section_midpoint=section_midpoint,
            section_radius=section_radius,
            center_radius=center_radius,
            diagnostic_total_radius=diagnostic_total_radius,
            task_state=state,
            level_shift=level_shift,
            query_basis=query_phi,
            support_baseline=support_baseline,
        )

    def point_section(self, query_phi: Tensor, support_features: Tensor) -> tuple[Tensor, Tensor]:
        """Return a row-space ambiguity diagnostic, not a coverage certificate."""
        if query_phi.ndim != 2 or query_phi.shape[1] != self.task_dim:
            raise ValueError("query_phi must have shape [Q,task_dim]")
        if self.task_dim == 0:
            zeros = torch.zeros(query_phi.shape[0], device=query_phi.device, dtype=query_phi.dtype)
            return zeros, zeros
        support_phi = self.basis(support_features)
        centered_support = support_phi - support_phi.mean(dim=0, keepdim=True)
        if centered_support.numel() == 0:
            null_projection = query_phi
        else:
            row_projector = torch.linalg.pinv(centered_support) @ centered_support
            null_projection = query_phi - query_phi @ row_projector
        radius = self.section_radius_bound * torch.linalg.vector_norm(null_projection, dim=1)
        midpoint = torch.zeros_like(radius)
        return midpoint, radius

    @staticmethod
    def delta(predictions: Tensor, left: Tensor, right: Tensor) -> Tensor:
        return predictions.index_select(0, right) - predictions.index_select(0, left)

    @staticmethod
    def rectangle(delta_a: Tensor, delta_b: Tensor) -> Tensor:
        return delta_a - delta_b
