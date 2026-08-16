"""Trainable MetaSieve V1 core for cold-target few-shot adaptation.

The inference path keeps the retained uncentered positive dual-ridge section.
Task scheduling and label perturbation are training interventions and are not
required to load a V1 predictor at test time.
"""
from __future__ import annotations

import math

import torch
import torch.nn as nn
from torch import Tensor


class MetaSieveV1(nn.Module):
    """Ligand population predictor plus identifiable target section."""

    def __init__(self, input_dim: int = 288, section_dim: int = 2,
                 ridge: float = 1.0, *, support_only_section: bool = False,
                 population_hidden_dim: int = 0,
                 pair_hidden_dim: int = 0,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if (input_dim < 1 or not 0 <= section_dim <= 5 or ridge <= 0
                or population_hidden_dim < 0):
            raise ValueError("invalid MetaSieve V1 dimensions or ridge")
        self.input_dim = int(input_dim)
        self.section_dim = int(section_dim)
        self.ridge = float(ridge)
        self.support_only_section = bool(support_only_section)
        self.population_hidden_dim = int(population_hidden_dim)
        self.pair_hidden_dim = int(pair_hidden_dim)
        if pair_hidden_dim < 0:
            raise ValueError("pair hidden dimension cannot be negative")
        self.population = (
            nn.Sequential(
                nn.Linear(input_dim, population_hidden_dim, dtype=dtype),
                nn.SiLU(),
                nn.Linear(population_hidden_dim, 1, dtype=dtype),
            )
            if population_hidden_dim else
            nn.Linear(input_dim, 1, dtype=dtype)
        )
        if section_dim:
            self.pair_encoder = (
                nn.Sequential(
                    nn.Linear(input_dim, pair_hidden_dim, dtype=dtype),
                    nn.SiLU(),
                )
                if pair_hidden_dim else nn.Identity()
            )
            self.raw_basis = nn.Parameter(
                torch.randn(
                    pair_hidden_dim or input_dim, section_dim, dtype=dtype)
                / math.sqrt(pair_hidden_dim or input_dim))
            if not self.support_only_section:
                self.population_coordinate = nn.Parameter(
                    torch.zeros(section_dim, dtype=dtype))

    def basis(self) -> Tensor:
        if not self.section_dim:
            raise RuntimeError("population-only V1 has no section basis")
        return torch.linalg.qr(self.raw_basis, mode="reduced").Q

    def _components_with_basis(
            self, ligand: Tensor, pair: Tensor,
            basis: Tensor | None) -> tuple[Tensor, Tensor | None]:
        population = self.population(ligand).squeeze(-1)
        if not self.section_dim:
            return population, None
        coordinates = self.pair_encoder(pair) @ basis
        if not self.support_only_section:
            population = population + coordinates @ self.population_coordinate
        return population, coordinates

    def components(self, ligand: Tensor, pair: Tensor) -> tuple[Tensor, Tensor | None]:
        if ligand.ndim != 2 or pair.ndim != 2 or ligand.shape != pair.shape:
            raise ValueError("ligand and pair features need the same [N,D] shape")
        if ligand.shape[1] != self.input_dim:
            raise ValueError("V1 feature dimension differs from the model contract")
        return self._components_with_basis(
            ligand, pair, self.basis() if self.section_dim else None)

    def batched_components(
            self, ligand: Tensor, pair: Tensor) -> tuple[Tensor, Tensor | None]:
        if (ligand.ndim != 3 or pair.shape != ligand.shape
                or ligand.shape[2] != self.input_dim):
            raise ValueError("batched components need matching [B,N,D] tensors")
        return self._components_with_basis(
            ligand, pair, self.basis() if self.section_dim else None)

    def episode(self, support_ligand: Tensor, support_pair: Tensor,
                support_y: Tensor, query_ligand: Tensor,
                query_pair: Tensor) -> Tensor:
        """Predict a query after an uncentered positive dual-ridge support solve."""
        if support_y.ndim != 1 or len(support_y) != len(support_ligand):
            raise ValueError("support labels must have one value per support row")
        if len(support_y) < 1:
            raise ValueError("V1 adaptation requires nonempty support")
        basis = self.basis() if self.section_dim else None
        support_population, support_coordinates = self._components_with_basis(
            support_ligand, support_pair, basis)
        query_population, query_coordinates = self._components_with_basis(
            query_ligand, query_pair, basis)
        if not self.section_dim:
            return query_population
        residual = support_y - support_population
        identity = torch.eye(
            len(support_y), device=support_y.device, dtype=support_y.dtype)
        dual = torch.linalg.solve(
            support_coordinates @ support_coordinates.T + self.ridge * identity,
            residual,
        )
        coefficient = support_coordinates.T @ dual
        return query_population + query_coordinates @ coefficient

    def batched_episode(
            self, support_ligand: Tensor, support_pair: Tensor,
            support_y: Tensor, query_ligand: Tensor,
            query_pair: Tensor) -> tuple[Tensor, Tensor]:
        """Vectorized episodes with shapes [B,K,D], [B,K], and [B,Q,D]."""
        if (support_ligand.ndim != 3 or support_pair.shape != support_ligand.shape
                or query_ligand.ndim != 3 or query_pair.shape != query_ligand.shape
                or support_y.shape != support_ligand.shape[:2]
                or support_ligand.shape[0] != query_ligand.shape[0]
                or support_ligand.shape[2] != self.input_dim
                or query_ligand.shape[2] != self.input_dim):
            raise ValueError("invalid batched MetaSieve V1 episode shapes")
        basis = self.basis() if self.section_dim else None
        support_population, support_coordinates = self._components_with_basis(
            support_ligand, support_pair, basis)
        query_population, query_coordinates = self._components_with_basis(
            query_ligand, query_pair, basis)
        if not self.section_dim:
            return query_population, support_population
        residual = support_y - support_population
        batch, support_size = support_y.shape
        identity = torch.eye(
            support_size, device=support_y.device, dtype=support_y.dtype)
        gram = support_coordinates @ support_coordinates.transpose(-1, -2)
        dual = torch.linalg.solve(
            gram + self.ridge * identity.expand(batch, -1, -1),
            residual.unsqueeze(-1),
        )
        coefficient = support_coordinates.transpose(-1, -2) @ dual
        prediction = query_population + (
            query_coordinates @ coefficient).squeeze(-1)
        return prediction, support_population


class TaskScheduler(nn.Module):
    """Small source-only scheduler over loss, gradient agreement and progress."""

    def __init__(self, hidden_dim: int = 8,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if hidden_dim < 1:
            raise ValueError("scheduler hidden dimension must be positive")
        self.network = nn.Sequential(
            nn.LayerNorm(3, dtype=dtype),
            nn.Linear(3, hidden_dim, dtype=dtype),
            nn.Tanh(),
            nn.Linear(hidden_dim, 1, dtype=dtype),
        )

    def forward(self, statistics: Tensor) -> Tensor:
        if statistics.ndim != 2 or statistics.shape[1] != 3:
            raise ValueError("task statistics must have shape [tasks,3]")
        torch._assert_async(
            torch.isfinite(statistics).all(), "task statistics must be finite")
        return self.network(statistics).squeeze(-1)

    def probabilities(self, statistics: Tensor) -> Tensor:
        return torch.softmax(self(statistics), dim=0)

    def sample(self, statistics: Tensor, count: int,
               generator: torch.Generator) -> Tensor:
        if not 1 <= count <= len(statistics):
            raise ValueError("selected task count is outside the candidate pool")
        return torch.multinomial(
            self.probabilities(statistics), count, replacement=False,
            generator=generator)


def uniform_label_noise(labels: Tensor, standard_deviation: float, *,
                        generator: torch.Generator) -> Tensor:
    """Add zero-mean uniform noise with a declared standard deviation."""
    if standard_deviation < 0:
        raise ValueError("label-noise standard deviation cannot be negative")
    if standard_deviation == 0:
        return labels
    bound = math.sqrt(3.0) * standard_deviation
    noise = torch.empty_like(labels).uniform_(-bound, bound, generator=generator)
    return labels + noise
