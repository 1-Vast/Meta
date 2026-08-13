"""Low-rank atom-residue geometry heads for the mechanistic bridge."""
from __future__ import annotations

from dataclasses import dataclass
import math

import torch
import torch.nn as nn

from .bpsf import BipartitePairSectionFormer

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


class GeometrySupervisionHead(nn.Module):
    """Source-only geometry teacher attached to the deployment pair field."""

    def __init__(self, pair_dim: int, distance_bins=DISTANCE_BINS_ANGSTROM,
                 dtype=torch.float32):
        super().__init__()
        self.n_distance_bins = len(tuple(distance_bins)) - 1
        self.contact = nn.Linear(pair_dim, 1, dtype=dtype)
        self.distance = nn.Linear(pair_dim, self.n_distance_bins, dtype=dtype)

    def forward(self, pair: torch.Tensor, pair_mask: torch.Tensor
                ) -> DenseMechanismPrediction:
        if pair.ndim != 4 or pair_mask.shape != pair.shape[:3]:
            raise ValueError("geometry head requires [B,A,R,H] pair states and mask")
        return DenseMechanismPrediction(
            self.contact(pair).squeeze(-1), self.distance(pair), pair_mask)


class PairGeometryTeacher(nn.Module):
    """Train BPSF pair states from holo labels without requiring holo input."""

    def __init__(self, hidden_dim: int, section_dim: int = 16,
                 blocks: int = 2, latents: int = 8, heads: int = 4,
                 chunk_size: int = 32, dtype=torch.float32):
        super().__init__()
        self.trunk = BipartitePairSectionFormer(
            hidden_dim, section_dim, blocks, latents, heads, chunk_size, dtype)
        self.geometry = GeometrySupervisionHead(hidden_dim, dtype=dtype)

    def forward(self, ligand_atoms: torch.Tensor, ligand_mask: torch.Tensor,
                protein_residues: torch.Tensor, residue_mask: torch.Tensor
                ) -> DenseMechanismPrediction:
        encoded = self.trunk(
            ligand_atoms, protein_residues, ligand_mask, residue_mask,
            return_pair=True)
        return self.geometry(encoded.pair, encoded.pair_mask)
