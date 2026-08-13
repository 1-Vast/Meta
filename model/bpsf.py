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
    mechanism_slots: Tensor
    pair: Tensor | None = None
    pair_mask: Tensor | None = None


class BipartitePairBlock(nn.Module):
    """Progressively update atom, residue, and rectangular pair states."""

    def __init__(self, token_dim: int, pair_dim: int, task_dim: int = 0,
                 adapter_rank: int = 0, adapter_scale: float = 0.25,
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
        self.atom_localizer = nn.Linear(pair_dim, 1, bias=False, dtype=dtype)
        self.residue_localizer = nn.Linear(pair_dim, 1, bias=False, dtype=dtype)
        self.atom_norm = nn.LayerNorm(token_dim, dtype=dtype)
        self.residue_norm = nn.LayerNorm(token_dim, dtype=dtype)
        self.atom_pair = nn.Linear(token_dim, pair_dim, bias=False, dtype=dtype)
        self.residue_pair = nn.Linear(token_dim, pair_dim, bias=False, dtype=dtype)
        self.pair_update = nn.Sequential(
            nn.LayerNorm(4 * pair_dim, dtype=dtype),
            nn.Linear(4 * pair_dim, 4 * pair_dim, dtype=dtype),
            nn.SiLU(), nn.Linear(4 * pair_dim, pair_dim, bias=False, dtype=dtype))
        self.adapter_scale = float(adapter_scale)
        if adapter_rank:
            if task_dim < 1:
                raise ValueError("task_dim is required for an adaptive pair block")
            self.adapter_gate = nn.Linear(
                task_dim, adapter_rank, bias=False, dtype=dtype)
            self.adapter_down = nn.Linear(
                pair_dim, adapter_rank, bias=False, dtype=dtype)
            self.adapter_up = nn.Linear(
                adapter_rank, pair_dim, bias=False, dtype=dtype)
        else:
            self.adapter_gate = None

    def forward(self, pair: Tensor, atoms: Tensor, residues: Tensor,
                atom_mask: Tensor, residue_mask: Tensor,
                adjacency: Tensor,
                task_code: Tensor | None = None) -> tuple[Tensor, Tensor, Tensor]:
        pair_mask = atom_mask[:, :, None] * residue_mask[:, None, :]
        atom_gate = atom_mask.unsqueeze(-1)
        residue_gate = residue_mask.unsqueeze(-1)
        atoms = atoms * atom_gate
        residues = residues * residue_gate
        adjacency = adjacency * atom_mask[:, :, None] * atom_mask[:, None, :]
        atom_score = self.atom_localizer(pair).squeeze(-1)
        residue_score = self.residue_localizer(pair).squeeze(-1)
        atom_score = atom_score.masked_fill(
            ~pair_mask.bool(), torch.finfo(atom_score.dtype).min)
        residue_score = residue_score.masked_fill(
            ~pair_mask.bool(), torch.finfo(residue_score.dtype).min)
        atom_weight = torch.softmax(atom_score, dim=2) * pair_mask
        residue_weight = torch.softmax(residue_score, dim=1) * pair_mask
        atom_context = (pair * atom_weight.unsqueeze(-1)).sum(2)
        residue_context = (pair * residue_weight.unsqueeze(-1)).sum(1)
        degree = adjacency.sum(-1, keepdim=True).clamp_min(1)
        neighbor = adjacency @ atoms / degree
        atoms = self.atom_norm(
            atoms + self.atom_neighbor(neighbor) + self.pair_to_atom(atom_context))
        local = self.residue_local(residues.transpose(1, 2)).transpose(1, 2)
        local = local * residue_gate
        residues = self.residue_norm(
            residues + local + self.pair_to_residue(residue_context))
        atoms = atoms * atom_gate
        residues = residues * residue_gate
        atom_pair = self.atom_pair(atoms)[:, :, None, :]
        residue_pair = self.residue_pair(residues)[:, None, :, :]
        refresh = self.pair_update(torch.cat((
            pair, atom_pair.expand_as(pair), residue_pair.expand_as(pair),
            atom_pair * residue_pair), dim=-1))
        if task_code is not None and self.adapter_gate is not None:
            gate = torch.tanh(self.adapter_gate(task_code))[:, None, None, :]
            refresh = refresh + self.adapter_scale * self.adapter_up(
                self.adapter_down(pair) * gate)
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

    def forward(self, pair: Tensor, pair_mask: Tensor) -> tuple[Tensor, Tensor]:
        batch, _, _, width = pair.shape
        values = pair.reshape(batch, -1, width)
        valid = pair_mask.reshape(batch, -1).bool()
        latent = self.latents.unsqueeze(0).expand(batch, -1, -1)
        attended, _ = self.cross(
            latent, values, values, key_padding_mask=~valid, need_weights=False)
        latent = latent + attended
        for block in self.blocks:
            latent = latent + block(latent)
        slots = self.norm(latent)
        return self.output(slots.mean(1)), slots

    def project_slots(self, slots: Tensor) -> Tensor:
        """Project retained aligned slots through the original pooled readout."""
        if slots.ndim < 3 or slots.shape[-2:] != self.latents.shape:
            raise ValueError(
                "mechanism slots must end in latent_count x pair_dim")
        return self.output(slots.mean(-2))


class SectionLatentEncoder(nn.Module):
    """Shared interaction latent with endpoint and meta-section projections."""

    def __init__(self, pair_dim: int, endpoint_dim: int, section_dim: int,
                 latent_count: int, heads: int,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.interaction = _LatentReadout(
            pair_dim, pair_dim, latent_count, heads, dtype)
        self.endpoint = nn.Linear(pair_dim, endpoint_dim, bias=False, dtype=dtype)
        self.section = nn.Linear(pair_dim, section_dim, bias=False, dtype=dtype)
        self.section_norm = nn.LayerNorm(
            section_dim, elementwise_affine=False, dtype=dtype)

    def project_slots(self, slots: Tensor) -> tuple[Tensor, Tensor]:
        shared = self.interaction.project_slots(slots)
        return self.endpoint(shared), self.section_norm(self.section(shared))

    def forward(self, pair: Tensor, pair_mask: Tensor) -> tuple[Tensor, Tensor, Tensor]:
        shared, slots = self.interaction(pair, pair_mask)
        return (self.endpoint(shared), self.section_norm(self.section(shared)),
                slots)


class BipartitePairSectionFormer(nn.Module):
    def __init__(self, token_dim: int, section_dim: int, pair_dim: int = 48,
                 blocks: int = 3, latent_count: int = 12, heads: int = 4,
                 chunk_size: int = 16, task_dim: int = 0,
                 adapter_rank: int = 0, adaptive_blocks: int = 0,
                 adapter_scale: float = 0.25,
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
        if adaptive_blocks < 0 or adaptive_blocks > blocks:
            raise ValueError("adaptive_blocks must be within the pair trunk")
        self.blocks = nn.ModuleList(
            BipartitePairBlock(
                token_dim, pair_dim,
                task_dim if index >= blocks - adaptive_blocks else 0,
                adapter_rank if index >= blocks - adaptive_blocks else 0,
                adapter_scale, dtype)
            for index in range(blocks))
        self.latent = SectionLatentEncoder(
            pair_dim, token_dim, section_dim, latent_count, heads, dtype)

    def _forward_chunk(self, atoms: Tensor, residues: Tensor, atom_mask: Tensor,
                       residue_mask: Tensor, adjacency: Tensor,
                       return_pair: bool,
                       task_code: Tensor | None = None) -> PairSectionEncoding:
        atoms = atoms * atom_mask.unsqueeze(-1)
        residues = residues * residue_mask.unsqueeze(-1)
        adjacency = adjacency * atom_mask[:, :, None] * atom_mask[:, None, :]
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
                pair, atoms, residues, atom_mask, residue_mask, adjacency,
                task_code)
        endpoint, section, slots = self.latent(pair, mask)
        return PairSectionEncoding(endpoint, section, slots,
                                   pair if return_pair else None,
                                   mask if return_pair else None)

    def forward(self, atoms: Tensor, residues: Tensor, atom_mask: Tensor,
                residue_mask: Tensor, adjacency: Tensor | None = None, *,
                task_code: Tensor | None = None,
                return_pair: bool = False) -> PairSectionEncoding:
        if atoms.ndim != 3 or residues.ndim != 3:
            raise ValueError("BPSF states must have rank three")
        if atom_mask.shape != atoms.shape[:2] or residue_mask.shape != residues.shape[:2]:
            raise ValueError("BPSF masks do not match states")
        if bool((atom_mask.sum(1) == 0).any()):
            raise ValueError("BPSF received a zero-atom sample")
        if bool((residue_mask.sum(1) == 0).any()):
            raise ValueError("BPSF received a zero-residue sample")
        if task_code is not None and task_code.shape[0] != atoms.shape[0]:
            raise ValueError("task code batch does not match BPSF states")
        if adjacency is None:
            adjacency = torch.zeros(
                atoms.shape[0], atoms.shape[1], atoms.shape[1],
                device=atoms.device, dtype=atoms.dtype)
        outputs = []
        for start in range(0, atoms.shape[0], self.chunk_size):
            stop = min(start + self.chunk_size, atoms.shape[0])
            outputs.append(self._forward_chunk(
                atoms[start:stop], residues[start:stop], atom_mask[start:stop],
                residue_mask[start:stop], adjacency[start:stop], return_pair,
                task_code[start:stop] if task_code is not None else None))
        return PairSectionEncoding(
            torch.cat([x.endpoint for x in outputs]),
            torch.cat([x.section for x in outputs]),
            torch.cat([x.mechanism_slots for x in outputs]),
            torch.cat([x.pair for x in outputs]) if return_pair else None,
            torch.cat([x.pair_mask for x in outputs]) if return_pair else None)
