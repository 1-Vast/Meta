"""Standalone retired QPSMP representation and analytic-section comparators."""
from __future__ import annotations

import math
import torch
import torch.nn as nn
from torch import Tensor


class PooledInteraction(nn.Module):
    def __init__(self, hidden_dim: int, dtype=torch.float32):
        super().__init__()
        self.key = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.query = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.value = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)

    def forward(self, residues: Tensor, residue_mask: Tensor, ligand: Tensor) -> Tensor:
        logits = torch.einsum("bmh,brh->bmr", self.query(ligand), self.key(residues))
        logits = logits / math.sqrt(residues.shape[-1])
        logits = logits.masked_fill(~residue_mask[:, None].bool(), -torch.inf)
        localized = torch.einsum(
            "bmr,brh->bmh", torch.softmax(logits, -1), self.value(residues))
        return localized * ligand


class AtomResiduePool(nn.Module):
    def __init__(self, hidden_dim: int, dtype=torch.float32):
        super().__init__()
        self.atom = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.residue = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.value = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)

    def forward(self, residues: Tensor, residue_mask: Tensor,
                atoms: Tensor, atom_mask: Tensor) -> Tensor:
        logits = torch.einsum(
            "bmnh,brh->bmnr", self.atom(atoms), self.residue(residues))
        logits = logits / math.sqrt(residues.shape[-1])
        logits = logits.masked_fill(~residue_mask[:, None, None].bool(), -torch.inf)
        localized = torch.einsum(
            "bmnr,brh->bmnh", torch.softmax(logits, -1), self.value(residues))
        pair = atoms * localized * atom_mask.unsqueeze(-1)
        return pair.sum(-2) / atom_mask.sum(-1, keepdim=True).clamp_min(1)


class CenteredRidgeSection(nn.Module):
    def __init__(self, ridge: float = 0.1):
        super().__init__()
        if ridge <= 0:
            raise ValueError("ridge must be positive")
        self.ridge = float(ridge)

    def forward(self, support: Tensor, query: Tensor, residual: Tensor):
        basis = support - support.mean(-2, keepdim=True)
        residual = residual - residual.mean(-1, keepdim=True)
        gram = basis @ basis.transpose(-1, -2)
        identity = torch.eye(gram.shape[-1], device=gram.device, dtype=gram.dtype)
        alpha = torch.linalg.solve(
            gram + self.ridge * gram.shape[-1] * identity, residual.unsqueeze(-1))
        state = (basis.transpose(-1, -2) @ alpha).squeeze(-1)
        return state, (query @ state.unsqueeze(-1)).squeeze(-1)


class LowRankGeometryBridge(nn.Module):
    def __init__(self, atom_dim: int, residue_dim: int, rank: int = 32,
                 dtype=torch.float32):
        super().__init__()
        self.atom = nn.Linear(atom_dim, rank, bias=False, dtype=dtype)
        self.residue = nn.Linear(residue_dim, rank, bias=False, dtype=dtype)

    def forward(self, atoms: Tensor, residues: Tensor) -> Tensor:
        return torch.einsum("bnr,blr->bnl", self.atom(atoms), self.residue(residues)) \
            / math.sqrt(self.atom.out_features)
