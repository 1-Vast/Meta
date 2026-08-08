"""Pair-local E0 Mechanistic Affinity Potential; production pipeline is untouched."""
from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn as nn
import torch.nn.functional as F

from research.e0_identifiability.mechanistic_affinity_contract import (
    ATOM_CHEMISTRY_DIM, ATOM_STATE_DIM, DISTANCE_PROB_DIM, MAP_HIDDEN_DIM,
    RESIDUE_CHEMISTRY_DIM, RESIDUE_STATE_DIM,
)


@dataclass
class MechanisticAffinityOutput:
    potential: torch.Tensor
    contact_mass: torch.Tensor
    pair_contribution: torch.Tensor


class LocalMechanisticAffinityPotential(nn.Module):
    """Apply nonlinear chemistry interpretation before atom-residue pooling."""

    def __init__(self, hidden_dim: int = MAP_HIDDEN_DIM):
        super().__init__()
        self.atom = nn.Sequential(
            nn.Linear(ATOM_STATE_DIM + ATOM_CHEMISTRY_DIM, hidden_dim), nn.SiLU())
        self.residue = nn.Sequential(
            nn.Linear(RESIDUE_STATE_DIM + RESIDUE_CHEMISTRY_DIM, hidden_dim), nn.SiLU())
        self.geometry = nn.Sequential(
            nn.Linear(1 + DISTANCE_PROB_DIM, hidden_dim // 2), nn.SiLU())
        pair_dim = hidden_dim * 3 + hidden_dim // 2
        self.pair = nn.Sequential(
            nn.Linear(pair_dim, hidden_dim), nn.SiLU(), nn.Linear(hidden_dim, 1))
        self.mass_weight = nn.Parameter(torch.zeros(()))

    def forward(self, atom_state: torch.Tensor, atom_chemistry: torch.Tensor,
                atom_mask: torch.Tensor, residue_state: torch.Tensor,
                residue_chemistry: torch.Tensor, residue_mask: torch.Tensor,
                contact_prob: torch.Tensor, distance_prob: torch.Tensor, *,
                use_geometry: bool = True) -> MechanisticAffinityOutput:
        if atom_state.shape[-1] != ATOM_STATE_DIM or \
                atom_chemistry.shape[-1] != ATOM_CHEMISTRY_DIM:
            raise ValueError("atom-local state does not match E0 contract")
        if residue_state.shape[-1] != RESIDUE_STATE_DIM or \
                residue_chemistry.shape[-1] != RESIDUE_CHEMISTRY_DIM:
            raise ValueError("residue-local state does not match E0 contract")
        pair_shape = (atom_state.shape[0], atom_state.shape[1], residue_state.shape[1])
        if contact_prob.shape != pair_shape or distance_prob.shape != (*pair_shape, 5):
            raise ValueError("geometry does not match E0 pair axes")
        pair_mask = atom_mask.float().unsqueeze(2) * residue_mask.float().unsqueeze(1)
        atom = self.atom(torch.cat((atom_state.float(), atom_chemistry.float()), dim=-1))
        residue = self.residue(torch.cat(
            (residue_state.float(), residue_chemistry.float()), dim=-1))
        if use_geometry:
            contact = contact_prob.float() * pair_mask
            geometry_values = torch.cat((contact.unsqueeze(-1), distance_prob.float()), dim=-1)
            weight = contact
        else:
            geometry_values = torch.zeros((*pair_shape, 6), device=atom.device)
            weight = pair_mask
        geometry = self.geometry(geometry_values)
        atom_pair = atom.unsqueeze(2).expand(-1, -1, pair_shape[2], -1)
        residue_pair = residue.unsqueeze(1).expand(-1, pair_shape[1], -1, -1)
        pair_input = torch.cat((atom_pair, residue_pair, atom_pair * residue_pair, geometry), -1)
        contribution = self.pair(pair_input).squeeze(-1) * pair_mask
        normalizer = weight.sum((1, 2)).clamp_min(1e-6)
        mean_energy = (contribution * weight).sum((1, 2)) / normalizer
        contact_mass = (contact_prob.float() * pair_mask).sum((1, 2)) / \
            pair_mask.sum((1, 2)).clamp_min(1.0)
        potential = mean_energy + self.mass_weight * contact_mass
        return MechanisticAffinityOutput(potential, contact_mass, contribution)


class EndpointCalibration(nn.Module):
    """Shared potential with monotone Ki/Kd scale and endpoint offset."""

    def __init__(self):
        super().__init__()
        self.raw_scale = nn.Parameter(torch.zeros(2))
        self.offset = nn.Parameter(torch.zeros(2))

    def forward(self, potential: torch.Tensor, endpoint_index: torch.Tensor) -> torch.Tensor:
        scale = F.softplus(self.raw_scale) + 1e-6
        return scale[endpoint_index] * potential + self.offset[endpoint_index]


def pairwise_rank_loss(scores: torch.Tensor, labels: torch.Tensor,
                       task_index: torch.Tensor) -> torch.Tensor:
    """Logistic ranking over all non-tied pairs within each assay task."""
    losses = []
    for task in torch.unique(task_index):
        indices = torch.where(task_index == task)[0]
        if len(indices) < 2:
            continue
        left, right = torch.triu_indices(len(indices), len(indices), offset=1,
                                         device=scores.device)
        signs = torch.sign(labels[indices[left]] - labels[indices[right]])
        keep = signs != 0
        if keep.any():
            delta = scores[indices[left[keep]]] - scores[indices[right[keep]]]
            losses.append(F.softplus(-signs[keep] * delta))
    if not losses:
        return scores.sum() * 0.0
    return torch.cat(losses).mean()


def e0_loss(scores: torch.Tensor, residual_targets: torch.Tensor,
            labels: torch.Tensor, task_index: torch.Tensor,
            rank_weight: float = 1.0) -> torch.Tensor:
    return F.huber_loss(scores, residual_targets) + rank_weight * pairwise_rank_loss(
        scores, labels, task_index)
