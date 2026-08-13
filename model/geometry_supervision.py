"""Source-only geometry supervision for the active BPSF pair field."""
from __future__ import annotations

from dataclasses import dataclass
import torch
import torch.nn as nn

from contracts.mechanism import DISTANCE_BINS_ANGSTROM, MECHANISM_SCHEMA
from .bpsf import BipartitePairSectionFormer


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
        return torch.softmax(self.distance_logits, -1) * self.pair_mask.unsqueeze(-1)


class GeometrySupervisionHead(nn.Module):
    schema = MECHANISM_SCHEMA

    def __init__(self, pair_dim: int, distance_bins=DISTANCE_BINS_ANGSTROM,
                 dtype=torch.float32):
        super().__init__()
        self.n_distance_bins = len(tuple(distance_bins)) - 1
        self.contact = nn.Linear(pair_dim, 1, dtype=dtype)
        self.distance = nn.Linear(pair_dim, self.n_distance_bins, dtype=dtype)

    def forward(self, pair: torch.Tensor,
                pair_mask: torch.Tensor) -> DenseMechanismPrediction:
        if pair.ndim != 4 or pair_mask.shape != pair.shape[:3]:
            raise ValueError("geometry head requires [B,A,R,Dp] pair states and mask")
        return DenseMechanismPrediction(
            self.contact(pair).squeeze(-1), self.distance(pair), pair_mask)


class PairGeometryTeacher(nn.Module):
    def __init__(self, hidden_dim: int = 128, section_dim: int = 32,
                 pair_dim: int = 64, blocks: int = 3, latents: int = 16,
                 heads: int = 8, chunk_size: int = 16,
                 dtype=torch.float32):
        super().__init__()
        self.trunk = BipartitePairSectionFormer(
            hidden_dim, section_dim, pair_dim, blocks, latents, heads,
            chunk_size, dtype)
        self.geometry = GeometrySupervisionHead(pair_dim, dtype=dtype)

    def forward(self, ligand_atoms: torch.Tensor, ligand_mask: torch.Tensor,
                protein_residues: torch.Tensor, residue_mask: torch.Tensor,
                adjacency: torch.Tensor | None = None) -> DenseMechanismPrediction:
        encoded = self.trunk(
            ligand_atoms, protein_residues, ligand_mask, residue_mask,
            adjacency, return_pair=True)
        return self.geometry(encoded.pair, encoded.pair_mask)
