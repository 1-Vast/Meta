"""Local protein and ligand encoders validated by the P1B geometry gate."""
from __future__ import annotations

import torch
import torch.nn as nn

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM


class ProteinEncoder(nn.Module):
    """Project cached protein-LM pooled and residue-local states."""

    def __init__(self, protein_dim: int, d_h: int, dtype=torch.float64):
        super().__init__()
        self.proj = nn.Sequential(
            nn.Linear(protein_dim, d_h, dtype=dtype),
            nn.SiLU(),
            nn.Linear(d_h, d_h, dtype=dtype),
        )
        self.bank_proj = nn.Linear(protein_dim, d_h, dtype=dtype)
        self.norm = nn.LayerNorm(d_h, dtype=dtype)

    def forward(self, pooled: torch.Tensor, bank: torch.Tensor):
        return self.norm(self.proj(pooled)), self.bank_proj(bank)


class GINELayer(nn.Module):
    """Dense deterministic GINE update used by the ligand tower."""

    def __init__(self, d_h: int, bond_dim: int, dtype=torch.float64):
        super().__init__()
        self.eps = nn.Parameter(torch.zeros((), dtype=dtype))
        self.edge = nn.Linear(bond_dim, d_h, dtype=dtype)
        self.mlp = nn.Sequential(
            nn.Linear(d_h, 2 * d_h, dtype=dtype),
            nn.SiLU(),
            nn.Linear(2 * d_h, d_h, dtype=dtype),
        )
        self.norm = nn.LayerNorm(d_h, dtype=dtype)

    def forward(self, states, bonds, adjacency, mask):
        edge = self.edge(bonds)
        messages = torch.relu(states.unsqueeze(1) + edge) * adjacency.unsqueeze(-1)
        output = self.mlp((1.0 + self.eps) * states + messages.sum(dim=2))
        return self.norm(output) * mask.unsqueeze(-1)


class LigandEncoder(nn.Module):
    """Map a padded molecular graph to pooled and atom-local states."""

    def __init__(self, d_h: int, n_layers: int = 4,
                 atom_dim: int = ATOM_FEAT_DIM,
                 bond_dim: int = BOND_FEAT_DIM, dtype=torch.float64):
        super().__init__()
        self.inp = nn.Linear(atom_dim, d_h, dtype=dtype)
        self.layers = nn.ModuleList(
            GINELayer(d_h, bond_dim, dtype) for _ in range(n_layers)
        )
        self.out = nn.Linear(2 * d_h, d_h, dtype=dtype)
        self.norm = nn.LayerNorm(d_h, dtype=dtype)

    def forward(self, atoms, bonds, mask):
        adjacency = (bonds.abs().sum(-1) > 0).to(atoms.dtype)
        states = self.inp(atoms) * mask.unsqueeze(-1)
        for layer in self.layers:
            states = states + layer(states, bonds, adjacency, mask)
        denominator = mask.sum(1, keepdim=True).clamp(min=1.0)
        mean = (states * mask.unsqueeze(-1)).sum(1) / denominator
        maximum = states.masked_fill(
            mask.unsqueeze(-1) == 0, torch.finfo(states.dtype).min
        ).amax(1)
        maximum = torch.nan_to_num(maximum, neginf=0.0)
        pooled = self.norm(self.out(torch.cat((mean, maximum), dim=-1)))
        return pooled, states
