"""Scalable bipartite pair-section modules for the active QPSMP model."""
from __future__ import annotations

from dataclasses import dataclass

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
    """Progressively update atom, residue, and rectangular pair states."""

    def __init__(self, token_dim: int, pair_dim: int,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.atom_neighbor = nn.Sequential(
            nn.Linear(token_dim, 2 * token_dim, bias=False, dtype=dtype),
            nn.SiLU(), nn.Linear(2 * token_dim, token_dim, bias=False, dtype=dtype))
        self.residue_local = nn.Conv1d(
            token_dim, token_dim, kernel_size=5, padding=2,
            groups=token_dim, bias=False, dtype=dtype)
        self.pair_to_atom = nn.Linear(pair_dim, token_dim, bias=False, dtype=dtype)
        self.pair_to_residue = nn.Linear(pair_dim, token_dim, bias=False, dtype=dtype)
        self.atom_norm = nn.LayerNorm(token_dim, dtype=dtype)
        self.residue_norm = nn.LayerNorm(token_dim, dtype=dtype)
        self.atom_pair = nn.Linear(token_dim, pair_dim, bias=False, dtype=dtype)
        self.residue_pair = nn.Linear(token_dim, pair_dim, bias=False, dtype=dtype)
        self.pair_update = nn.Sequential(
            nn.LayerNorm(4 * pair_dim, dtype=dtype),
            nn.Linear(4 * pair_dim, 4 * pair_dim, dtype=dtype),
            nn.SiLU(), nn.Linear(4 * pair_dim, pair_dim, bias=False, dtype=dtype))

    def forward(self, pair: Tensor, atoms: Tensor, residues: Tensor,
                atom_mask: Tensor, residue_mask: Tensor,
                adjacency: Tensor) -> tuple[Tensor, Tensor, Tensor]:
        pair_mask = atom_mask[:, :, None] * residue_mask[:, None, :]
        weighted = pair * pair_mask.unsqueeze(-1)
        atom_context = weighted.sum(2) / residue_mask.sum(1).clamp_min(1).view(-1, 1, 1)
        residue_context = weighted.sum(1) / atom_mask.sum(1).clamp_min(1).view(-1, 1, 1)
        degree = adjacency.sum(-1, keepdim=True).clamp_min(1)
        neighbor = adjacency @ atoms / degree
        atoms = self.atom_norm(
            atoms + self.atom_neighbor(neighbor) + self.pair_to_atom(atom_context))
        local = self.residue_local(residues.transpose(1, 2)).transpose(1, 2)
        residues = self.residue_norm(
            residues + local + self.pair_to_residue(residue_context))
        atoms = atoms * atom_mask.unsqueeze(-1)
        residues = residues * residue_mask.unsqueeze(-1)
        atom_pair = self.atom_pair(atoms)[:, :, None, :]
        residue_pair = self.residue_pair(residues)[:, None, :, :]
        refresh = self.pair_update(torch.cat((
            pair, atom_pair.expand_as(pair), residue_pair.expand_as(pair),
            atom_pair * residue_pair), dim=-1))
        pair = (pair + refresh) * pair_mask.unsqueeze(-1)
        return pair, atoms, residues


class _LatentReadout(nn.Module):
    def __init__(self, pair_dim: int, output_dim: int, latent_count: int,
                 heads: int, dtype: torch.dtype) -> None:
        super().__init__()
        self.latents = nn.Parameter(torch.empty(latent_count, pair_dim, dtype=dtype))
        nn.init.normal_(self.latents, std=pair_dim ** -0.5)
        self.cross = nn.MultiheadAttention(pair_dim, heads, batch_first=True, dtype=dtype)
        self.blocks = nn.ModuleList([
            nn.Sequential(nn.LayerNorm(pair_dim, dtype=dtype),
                          nn.Linear(pair_dim, 4 * pair_dim, dtype=dtype), nn.GELU(),
                          nn.Linear(4 * pair_dim, pair_dim, dtype=dtype))
            for _ in range(2)])
        self.norm = nn.LayerNorm(pair_dim, dtype=dtype)
        self.output = nn.Linear(pair_dim, output_dim, bias=False, dtype=dtype)

    def forward(self, pair: Tensor, pair_mask: Tensor) -> Tensor:
        batch, _, _, width = pair.shape
        values = pair.reshape(batch, -1, width)
        valid = pair_mask.reshape(batch, -1).bool()
        latent = self.latents.unsqueeze(0).expand(batch, -1, -1)
        attended, _ = self.cross(
            latent, values, values, key_padding_mask=~valid, need_weights=False)
        latent = latent + attended
        for block in self.blocks:
            latent = latent + block(latent)
        pooled = self.norm(latent).mean(1)
        return self.output(pooled)


class SectionLatentEncoder(nn.Module):
    """Independent endpoint and meta-section readouts over one pair field."""

    def __init__(self, pair_dim: int, endpoint_dim: int, section_dim: int,
                 latent_count: int, heads: int,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.endpoint = _LatentReadout(
            pair_dim, endpoint_dim, latent_count, heads, dtype)
        self.section = _LatentReadout(
            pair_dim, section_dim, latent_count, heads, dtype)

    def forward(self, pair: Tensor, pair_mask: Tensor) -> tuple[Tensor, Tensor]:
        return self.endpoint(pair, pair_mask), self.section(pair, pair_mask)


class BipartitePairSectionFormer(nn.Module):
    def __init__(self, token_dim: int, section_dim: int, pair_dim: int = 48,
                 blocks: int = 3, latent_count: int = 12, heads: int = 4,
                 chunk_size: int = 16,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if min(token_dim, pair_dim, section_dim, blocks, latent_count, chunk_size) < 1:
            raise ValueError("BPSF dimensions must be positive")
        if pair_dim % heads:
            raise ValueError("pair_dim must be divisible by attention heads")
        self.chunk_size = int(chunk_size)
        self.atom_pair = nn.Linear(token_dim, pair_dim, bias=False, dtype=dtype)
        self.residue_pair = nn.Linear(token_dim, pair_dim, bias=False, dtype=dtype)
        self.pair_init = nn.Sequential(
            nn.Linear(3 * pair_dim, pair_dim, bias=False, dtype=dtype),
            nn.SiLU(), nn.LayerNorm(pair_dim, dtype=dtype))
        self.blocks = nn.ModuleList(
            BipartitePairBlock(token_dim, pair_dim, dtype) for _ in range(blocks))
        self.latent = SectionLatentEncoder(
            pair_dim, token_dim, section_dim, latent_count, heads, dtype)

    def _forward_chunk(self, atoms: Tensor, residues: Tensor, atom_mask: Tensor,
                       residue_mask: Tensor, adjacency: Tensor,
                       return_pair: bool) -> PairSectionEncoding:
        atom_pair = self.atom_pair(atoms)[:, :, None, :]
        residue_pair = self.residue_pair(residues)[:, None, :, :]
        pair = self.pair_init(torch.cat((
            atom_pair.expand(-1, -1, residues.shape[1], -1),
            residue_pair.expand(-1, atoms.shape[1], -1, -1),
            atom_pair * residue_pair), dim=-1))
        mask = atom_mask[:, :, None] * residue_mask[:, None, :]
        pair = pair * mask.unsqueeze(-1)
        for block in self.blocks:
            pair, atoms, residues = block(
                pair, atoms, residues, atom_mask, residue_mask, adjacency)
        endpoint, section = self.latent(pair, mask)
        return PairSectionEncoding(endpoint, section,
                                   pair if return_pair else None,
                                   mask if return_pair else None)

    def forward(self, atoms: Tensor, residues: Tensor, atom_mask: Tensor,
                residue_mask: Tensor, adjacency: Tensor | None = None, *,
                return_pair: bool = False) -> PairSectionEncoding:
        if atoms.ndim != 3 or residues.ndim != 3:
            raise ValueError("BPSF states must have rank three")
        if atom_mask.shape != atoms.shape[:2] or residue_mask.shape != residues.shape[:2]:
            raise ValueError("BPSF masks do not match states")
        if bool((atom_mask.sum(1) == 0).any()):
            raise ValueError("BPSF received a zero-atom sample")
        if bool((residue_mask.sum(1) == 0).any()):
            raise ValueError("BPSF received a zero-residue sample")
        if adjacency is None:
            adjacency = torch.zeros(
                atoms.shape[0], atoms.shape[1], atoms.shape[1],
                device=atoms.device, dtype=atoms.dtype)
        outputs = []
        for start in range(0, atoms.shape[0], self.chunk_size):
            stop = min(start + self.chunk_size, atoms.shape[0])
            outputs.append(self._forward_chunk(
                atoms[start:stop], residues[start:stop], atom_mask[start:stop],
                residue_mask[start:stop], adjacency[start:stop], return_pair))
        return PairSectionEncoding(
            torch.cat([x.endpoint for x in outputs]),
            torch.cat([x.section for x in outputs]),
            torch.cat([x.pair for x in outputs]) if return_pair else None,
            torch.cat([x.pair_mask for x in outputs]) if return_pair else None)


class QuotientSupportSetOperator(nn.Module):
    def __init__(self, section_dim: int, hidden_dim: int = 96, heads: int = 4,
                 blocks: int = 2, state_bound: float = 1.0,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if hidden_dim % heads:
            raise ValueError("support hidden_dim must be divisible by heads")
        self.state_bound = float(state_bound)
        self.token = nn.Sequential(
            nn.Linear(section_dim + 3, hidden_dim, dtype=dtype), nn.GELU(),
            nn.LayerNorm(hidden_dim, dtype=dtype))
        self.attention = nn.ModuleList(
            nn.MultiheadAttention(hidden_dim, heads, batch_first=True, dtype=dtype)
            for _ in range(blocks))
        self.ffn = nn.ModuleList(
            nn.Sequential(nn.LayerNorm(hidden_dim, dtype=dtype),
                          nn.Linear(hidden_dim, 4 * hidden_dim, dtype=dtype), nn.GELU(),
                          nn.Linear(4 * hidden_dim, hidden_dim, dtype=dtype))
            for _ in range(blocks))
        self.weight = nn.Linear(hidden_dim, 1, bias=False, dtype=dtype)

    def forward(self, support: Tensor, query: Tensor,
                residual: Tensor) -> tuple[Tensor, Tensor, Tensor]:
        support_dim = support.ndim - 2
        residual = residual - residual.mean(-1, keepdim=True)
        center = support.mean(support_dim, keepdim=True)
        basis = support - center
        query_basis = query - center
        gram = basis @ basis.transpose(-1, -2)
        token = self.token(torch.cat((
            basis, residual.unsqueeze(-1),
            gram.diagonal(dim1=-2, dim2=-1).unsqueeze(-1),
            gram.mean(-1, keepdim=True)), dim=-1))
        for attention, ffn in zip(self.attention, self.ffn):
            update, _ = attention(token, token, token, need_weights=False)
            token = token + update
            token = token + ffn(token)
        alpha = residual * torch.tanh(self.weight(token).squeeze(-1))
        alpha = alpha - alpha.mean(-1, keepdim=True)
        state = (basis.transpose(-1, -2) @ alpha.unsqueeze(-1)).squeeze(-1)
        state = state / max(support.shape[support_dim], 1)
        state = self.state_bound * state / (
            self.state_bound + torch.linalg.vector_norm(state, dim=-1, keepdim=True))
        sar = (query_basis @ state.unsqueeze(-1)).squeeze(-1)
        return state, query_basis, sar
