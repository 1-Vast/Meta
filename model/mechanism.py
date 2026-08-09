"""Low-rank atom-residue geometry heads for the mechanistic bridge."""
from __future__ import annotations

from dataclasses import dataclass
import math

import torch
import torch.nn as nn

from contracts.mechanism import DISTANCE_BINS_ANGSTROM, MECHANISM_SCHEMA


@dataclass
class DenseMechanismPrediction:
    contact_logits: torch.Tensor
    distance_logits: torch.Tensor
    pair_mask: torch.Tensor

    @property
    def contact_prob(self) -> torch.Tensor:
        return torch.sigmoid(self.contact_logits) * self.pair_mask

    @property
    def distance_prob(self) -> torch.Tensor:
        return torch.softmax(self.distance_logits, dim=-1) * self.pair_mask.unsqueeze(-1)


class MechanisticInteractionBridge(nn.Module):
    """Predict contact and distance bins without a [B,N,L,H] pair tensor."""

    schema = MECHANISM_SCHEMA

    def __init__(self, d_atom: int, d_residue: int, rank: int = 32,
                 distance_bins=DISTANCE_BINS_ANGSTROM, dtype=torch.float32):
        super().__init__()
        if min(d_atom, d_residue, rank) < 1:
            raise ValueError("bridge dimensions must be positive")
        bins = tuple(float(value) for value in distance_bins)
        if len(bins) < 3 or bins[0] != 0.0 or any(b <= a for a, b in zip(bins, bins[1:])):
            raise ValueError("distance bin edges must start at zero and increase")
        self.distance_bins = bins
        self.n_distance_bins = len(bins) - 1
        self.rank = int(rank)
        self.contact_atom = nn.Linear(d_atom, rank, bias=False, dtype=dtype)
        self.contact_residue = nn.Linear(d_residue, rank, bias=False, dtype=dtype)
        self.distance_atom = nn.Linear(d_atom, rank * self.n_distance_bins,
                                       bias=False, dtype=dtype)
        self.distance_residue = nn.Linear(d_residue, rank * self.n_distance_bins,
                                          bias=False, dtype=dtype)

    def forward(self, ligand_atoms: torch.Tensor, ligand_mask: torch.Tensor,
                protein_residues: torch.Tensor, residue_mask: torch.Tensor
                ) -> DenseMechanismPrediction:
        if ligand_atoms.ndim != 3 or protein_residues.ndim != 3:
            raise ValueError("bridge states must have shape [B,N,D] and [B,L,D]")
        if ligand_mask.shape != ligand_atoms.shape[:2] or \
                residue_mask.shape != protein_residues.shape[:2]:
            raise ValueError("bridge masks do not match state shapes")
        if bool((ligand_mask.sum(dim=1) <= 0).any().item()):
            raise ValueError("bridge received a zero-atom ligand")
        if bool((residue_mask.sum(dim=1) <= 0).any().item()):
            raise ValueError("bridge received a zero-residue protein")
        pair_mask = ligand_mask.to(ligand_atoms.dtype).unsqueeze(-1) * \
            residue_mask.to(ligand_atoms.dtype).unsqueeze(-2)
        scale = math.sqrt(self.rank)
        contact = torch.einsum(
            "bnr,blr->bnl", self.contact_atom(ligand_atoms),
            self.contact_residue(protein_residues)) / scale
        batch, atoms = ligand_atoms.shape[:2]
        residues = protein_residues.shape[1]
        atom_distance = self.distance_atom(ligand_atoms).reshape(
            batch, atoms, self.n_distance_bins, self.rank)
        residue_distance = self.distance_residue(protein_residues).reshape(
            batch, residues, self.n_distance_bins, self.rank)
        distance = torch.einsum("bndr,bldr->bnld", atom_distance,
                                residue_distance) / scale
        return DenseMechanismPrediction(contact, distance, pair_mask)
