"""Bipartite pair and support-identifiable section modules for QPSMP."""
from __future__ import annotations

from dataclasses import dataclass
import math

import torch
import torch.nn as nn
from torch import Tensor


@dataclass(frozen=True)
class PairSectionEncoding:
    endpoint: Tensor
    section: Tensor
    pair: Tensor | None = None
    pair_mask: Tensor | None = None


class BipartitePairBlock(nn.Module):
    """Linear-complexity rectangular update over an atom-residue pair field."""

    def __init__(self, hidden_dim: int, dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.atom_feedback = nn.Linear(2 * hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.residue_feedback = nn.Linear(2 * hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.outer_atom = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.outer_residue = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.update = nn.Sequential(
            nn.LayerNorm(hidden_dim, dtype=dtype),
            nn.Linear(hidden_dim, 2 * hidden_dim, dtype=dtype),
            nn.SiLU(),
            nn.Linear(2 * hidden_dim, hidden_dim, bias=False, dtype=dtype),
        )

    def forward(self, pair: Tensor, atoms: Tensor, residues: Tensor,
                atom_mask: Tensor, residue_mask: Tensor) -> Tensor:
        mask = atom_mask[:, :, None] * residue_mask[:, None, :]
        weighted = pair * mask.unsqueeze(-1)
        atom_context = weighted.sum(2) / residue_mask.sum(1).clamp_min(1).view(-1, 1, 1)
        residue_context = weighted.sum(1) / atom_mask.sum(1).clamp_min(1).view(-1, 1, 1)
        atom_update = self.atom_feedback(torch.cat((atoms, atom_context), dim=-1))
        residue_update = self.residue_feedback(torch.cat((residues, residue_context), dim=-1))
        outer = self.outer_atom(atom_update)[:, :, None, :] * \
            self.outer_residue(residue_update)[:, None, :, :]
        value = pair + atom_update[:, :, None, :] + residue_update[:, None, :, :] + outer
        return (pair + self.update(value)) * mask.unsqueeze(-1)


class SectionLatentEncoder(nn.Module):
    """Compress a masked pair field with a small bank of learned latent queries."""

    def __init__(self, hidden_dim: int, section_dim: int, latent_count: int = 8,
                 heads: int = 4, dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if latent_count < 1 or hidden_dim % heads:
            raise ValueError("latent_count must be positive and hidden_dim divisible by heads")
        self.latents = nn.Parameter(torch.empty(latent_count, hidden_dim, dtype=dtype))
        nn.init.normal_(self.latents, std=hidden_dim ** -0.5)
        self.cross = nn.MultiheadAttention(
            hidden_dim, heads, batch_first=True, dtype=dtype)
        self.norm1 = nn.LayerNorm(hidden_dim, dtype=dtype)
        self.ffn = nn.Sequential(
            nn.LayerNorm(hidden_dim, dtype=dtype),
            nn.Linear(hidden_dim, 2 * hidden_dim, dtype=dtype),
            nn.SiLU(),
            nn.Linear(2 * hidden_dim, hidden_dim, dtype=dtype),
        )
        self.endpoint = nn.Sequential(
            nn.LayerNorm(hidden_dim, dtype=dtype),
            nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype),
        )
        self.section = nn.Linear(hidden_dim, section_dim, bias=False, dtype=dtype)

    def forward(self, pair: Tensor, pair_mask: Tensor) -> tuple[Tensor, Tensor]:
        batch, _, _, hidden = pair.shape
        values = pair.reshape(batch, -1, hidden)
        valid = pair_mask.reshape(batch, -1).bool()
        latent = self.latents.unsqueeze(0).expand(batch, -1, -1)
        attended, _ = self.cross(
            latent, values, values, key_padding_mask=~valid, need_weights=False)
        latent = self.norm1(latent + attended)
        latent = latent + self.ffn(latent)
        weights = valid.to(values.dtype).unsqueeze(-1)
        pair_pooled = (values * weights).sum(1) / weights.sum(1).clamp_min(1.0)
        return self.endpoint(pair_pooled), self.section(latent.mean(1))


class BipartitePairSectionFormer(nn.Module):
    """Persistent atom-residue field followed by pair-to-section attention."""

    def __init__(self, hidden_dim: int, section_dim: int, blocks: int = 2,
                 latent_count: int = 8, heads: int = 4, chunk_size: int = 64,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if min(hidden_dim, section_dim, blocks, chunk_size) < 1:
            raise ValueError("pair-section dimensions must be positive")
        self.chunk_size = int(chunk_size)
        self.atom = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.residue = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.pair_init = nn.Sequential(
            nn.Linear(3 * hidden_dim, hidden_dim, bias=False, dtype=dtype),
            nn.SiLU(),
            nn.LayerNorm(hidden_dim, dtype=dtype),
        )
        self.blocks = nn.ModuleList(
            BipartitePairBlock(hidden_dim, dtype) for _ in range(blocks))
        self.latent = SectionLatentEncoder(
            hidden_dim, section_dim, latent_count, heads, dtype)

    def _chunk(self, atoms: Tensor, residues: Tensor, atom_mask: Tensor,
               residue_mask: Tensor, return_pair: bool) -> PairSectionEncoding:
        atom = self.atom(atoms)
        residue = self.residue(residues)
        pair = self.pair_init(torch.cat((
            atom[:, :, None, :].expand(-1, -1, residue.shape[1], -1),
            residue[:, None, :, :].expand(-1, atom.shape[1], -1, -1),
            atom[:, :, None, :] * residue[:, None, :, :]), dim=-1))
        pair_mask = atom_mask[:, :, None] * residue_mask[:, None, :]
        pair = pair * pair_mask.unsqueeze(-1)
        for block in self.blocks:
            pair = block(pair, atom, residue, atom_mask, residue_mask)
        endpoint, section = self.latent(pair, pair_mask)
        return PairSectionEncoding(
            endpoint, section, pair if return_pair else None,
            pair_mask if return_pair else None)

    def forward(self, atoms: Tensor, residues: Tensor, atom_mask: Tensor,
                residue_mask: Tensor, *, return_pair: bool = False) -> PairSectionEncoding:
        if atoms.ndim != 3 or residues.ndim != 3:
            raise ValueError("pair-section states must have rank three")
        if atom_mask.shape != atoms.shape[:2] or residue_mask.shape != residues.shape[:2]:
            raise ValueError("pair-section masks do not match states")
        if bool((atom_mask.sum(1) == 0).any()) or bool((residue_mask.sum(1) == 0).any()):
            raise ValueError("pair-section inputs must contain valid atoms and residues")
        outputs = []
        for start in range(0, atoms.shape[0], self.chunk_size):
            stop = min(start + self.chunk_size, atoms.shape[0])
            outputs.append(self._chunk(
                atoms[start:stop], residues[start:stop], atom_mask[start:stop],
                residue_mask[start:stop], return_pair))
        return PairSectionEncoding(
            torch.cat([item.endpoint for item in outputs]),
            torch.cat([item.section for item in outputs]),
            torch.cat([item.pair for item in outputs]) if return_pair else None,
            torch.cat([item.pair_mask for item in outputs]) if return_pair else None)


class QuotientSupportSetOperator(nn.Module):
    """Learned support-to-state map with exact quotient and row-span invariants."""

    def __init__(self, section_dim: int, hidden_dim: int = 64, heads: int = 4,
                 state_bound: float = 1.0,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if hidden_dim % heads:
            raise ValueError("hidden_dim must be divisible by heads")
        self.state_bound = float(state_bound)
        self.token = nn.Sequential(
            nn.Linear(section_dim + 3, hidden_dim, dtype=dtype),
            nn.SiLU(),
            nn.LayerNorm(hidden_dim, dtype=dtype),
        )
        self.attention = nn.MultiheadAttention(
            hidden_dim, heads, batch_first=True, dtype=dtype)
        self.ffn = nn.Sequential(
            nn.LayerNorm(hidden_dim, dtype=dtype),
            nn.Linear(hidden_dim, 2 * hidden_dim, dtype=dtype),
            nn.SiLU(),
            nn.Linear(2 * hidden_dim, hidden_dim, dtype=dtype),
        )
        self.weight = nn.Linear(hidden_dim, 1, bias=False, dtype=dtype)

    def forward(self, support_basis: Tensor, query_basis: Tensor,
                centered_residual: Tensor) -> tuple[Tensor, Tensor, Tensor]:
        support_dim = support_basis.ndim - 2
        centered_residual = centered_residual - centered_residual.mean(
            dim=-1, keepdim=True)
        center = support_basis.mean(dim=support_dim, keepdim=True)
        basis = support_basis - center
        query = query_basis - center
        gram = basis @ basis.transpose(-1, -2)
        token = self.token(torch.cat((
            basis, centered_residual.unsqueeze(-1),
            gram.diagonal(dim1=-2, dim2=-1).unsqueeze(-1),
            gram.mean(-1, keepdim=True)), dim=-1))
        update, _ = self.attention(token, token, token, need_weights=False)
        token = token + update
        token = token + self.ffn(token)
        alpha = centered_residual * torch.tanh(self.weight(token).squeeze(-1))
        alpha = alpha - alpha.mean(dim=-1, keepdim=True)
        state = (basis.transpose(-1, -2) @ alpha.unsqueeze(-1)).squeeze(-1)
        state = state / max(support_basis.shape[support_dim], 1)
        state = self.state_bound * state / (
            self.state_bound + torch.linalg.vector_norm(state, dim=-1, keepdim=True))
        sar = (query @ state.unsqueeze(-1)).squeeze(-1)
        return state, query, sar
