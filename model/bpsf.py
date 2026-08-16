"""Scalable bipartite pair-section modules for the active QPSMP model."""
from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn as nn
from torch import Tensor

from contracts.mechanism import DISTANCE_BINS_ANGSTROM, MECHANISM_SCHEMA


@dataclass(frozen=True)
class PairSectionEncoding:
    endpoint: Tensor
    section: Tensor
    mechanism_slots: Tensor
    mechanism_response: Tensor
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
        # Localized residues are an attention-selected set, not a contiguous
        # sequence. A shared token update preserves permutation equivariance.
        self.residue_local = nn.Sequential(
            nn.Linear(token_dim, 2 * token_dim, bias=False, dtype=dtype),
            nn.SiLU(),
            nn.Linear(2 * token_dim, token_dim, bias=False, dtype=dtype))
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
        local = self.residue_local(residues)
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
        # Each index is a globally shared interaction primitive, rather than
        # an exchangeable episode-specific latent.
        self.latents = nn.Parameter(torch.empty(latent_count, pair_dim, dtype=dtype))
        nn.init.normal_(self.latents, std=pair_dim ** -0.5)
        del heads
        self.value = nn.Linear(pair_dim, pair_dim, bias=False, dtype=dtype)
        self.response_weight = nn.Parameter(
            torch.empty(latent_count, pair_dim, dtype=dtype))
        nn.init.normal_(self.response_weight, std=pair_dim ** -0.5)
        self.blocks = nn.ModuleList([
            nn.Sequential(nn.LayerNorm(pair_dim, dtype=dtype),
                          nn.Linear(pair_dim, 4 * pair_dim, dtype=dtype), nn.GELU(),
                          nn.Linear(4 * pair_dim, pair_dim, dtype=dtype))
            for _ in range(2)])
        self.norm = nn.LayerNorm(pair_dim, dtype=dtype)
        self.output = nn.Linear(pair_dim, output_dim, bias=False, dtype=dtype)

    def primitive_response(self, slots: Tensor) -> Tensor:
        if slots.shape[-2:] != self.latents.shape:
            raise ValueError("mechanism slots do not match primitive dictionary")
        raw = torch.einsum("...md,md->...m", slots, self.response_weight) \
            / self.response_weight.shape[-1] ** 0.5
        return torch.sigmoid(raw)

    def forward(self, pair: Tensor, pair_mask: Tensor,
                primitive_bias: Tensor | None = None) -> tuple[Tensor, Tensor, Tensor]:
        batch, _, _, width = pair.shape
        values = pair.reshape(batch, -1, width)
        valid = pair_mask.reshape(batch, -1).bool()
        score = torch.einsum("bnd,md->bmn", values, self.latents) / width ** 0.5
        if primitive_bias is not None:
            expected = (*pair.shape[:3], self.latents.shape[0])
            if primitive_bias.shape != expected:
                raise ValueError(f"primitive bias must have shape {expected}")
            score = score + primitive_bias.flatten(1, 2).transpose(1, 2)
        score = score.masked_fill(~valid[:, None], torch.finfo(score.dtype).min)
        weight = torch.softmax(score, -1)
        attended = torch.einsum("bmn,bnd->bmd", weight, self.value(values))
        latent = self.latents.unsqueeze(0).expand(batch, -1, -1)
        latent = latent + attended
        for block in self.blocks:
            latent = latent + block(latent)
        slots = self.norm(latent)
        # Responses exclude the dictionary-only identity baseline, preventing
        # pair-independent primitive constants from acting as level experts.
        response = self.primitive_response(self.norm(attended))
        return self.output(slots.mean(1)), slots, response

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

    def forward(self, pair: Tensor, pair_mask: Tensor,
                primitive_bias: Tensor | None = None) -> tuple[Tensor, Tensor, Tensor, Tensor]:
        shared, slots, response = self.interaction(
            pair, pair_mask, primitive_bias)
        return (self.endpoint(shared), self.section_norm(self.section(shared)),
                slots, response)


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
        # Sequence/2D-compatible physicochemical compatibility bias. Atom
        # features are the governed RDKit vector; protein region types are
        # inferred from residue embeddings without requiring coordinates.
        self.atom_primitive = nn.Linear(32, latent_count, bias=False, dtype=dtype)
        self.residue_primitive = nn.Linear(
            token_dim, latent_count, bias=False, dtype=dtype)
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
                       atom_features: Tensor | None = None,
                       atom_chemistry: Tensor | None = None,
                       residue_chemistry: Tensor | None = None,
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
        primitive_bias = None
        if atom_features is not None:
            if atom_features.shape[:2] != atoms.shape[:2] or atom_features.shape[-1] != 32:
                raise ValueError("raw atom features do not match the primitive field")
            atom_bias = self.atom_primitive(atom_features)[:, :, None]
            residue_bias = self.residue_primitive(residues)[:, None]
            primitive_bias = torch.tanh(atom_bias + residue_bias)
        if atom_chemistry is not None or residue_chemistry is not None:
            if atom_chemistry is None or residue_chemistry is None:
                raise ValueError("anchored primitives need atom and residue chemistry")
            if atom_chemistry.shape != (*atoms.shape[:2], 4) or \
                    residue_chemistry.shape != (*residues.shape[:2], 4):
                raise ValueError("anchored primitive chemistry has wrong shape")
            positive_a, negative_a, aromatic_a, hydrophobic_a = atom_chemistry.unbind(-1)
            positive_r, negative_r, aromatic_r, hydrophobic_r = residue_chemistry.unbind(-1)
            anchored = torch.stack((
                positive_a[:, :, None] * negative_r[:, None, :],
                negative_a[:, :, None] * positive_r[:, None, :],
                aromatic_a[:, :, None] * aromatic_r[:, None, :],
                hydrophobic_a[:, :, None] * hydrophobic_r[:, None, :]), -1)
            missing = self.latent.interaction.latents.shape[0] - anchored.shape[-1]
            if missing < 0:
                raise ValueError("primitive dictionary cannot hold anchored channels")
            anchored = torch.cat((anchored, anchored.new_zeros(
                *anchored.shape[:-1], missing)), -1)
            primitive_bias = anchored if primitive_bias is None else primitive_bias + anchored
        endpoint, section, slots, response = self.latent(
            pair, mask, primitive_bias)
        return PairSectionEncoding(endpoint, section, slots, response,
                                   pair if return_pair else None,
                                   mask if return_pair else None)

    def forward(self, atoms: Tensor, residues: Tensor, atom_mask: Tensor,
                residue_mask: Tensor, adjacency: Tensor | None = None, *,
                atom_features: Tensor | None = None,
                atom_chemistry: Tensor | None = None,
                residue_chemistry: Tensor | None = None,
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
                atom_features[start:stop] if atom_features is not None else None,
                atom_chemistry[start:stop] if atom_chemistry is not None else None,
                residue_chemistry[start:stop] if residue_chemistry is not None else None,
                task_code[start:stop] if task_code is not None else None))
        return PairSectionEncoding(
            torch.cat([x.endpoint for x in outputs]),
            torch.cat([x.section for x in outputs]),
            torch.cat([x.mechanism_slots for x in outputs]),
            torch.cat([x.mechanism_response for x in outputs]),
            torch.cat([x.pair for x in outputs]) if return_pair else None,
            torch.cat([x.pair_mask for x in outputs]) if return_pair else None)


@dataclass(frozen=True)
class DenseMechanismPrediction:
    """Source-only contact and distance predictions from a BPSF pair field."""

    contact_logits: Tensor
    distance_logits: Tensor
    pair_mask: Tensor

    @property
    def contact_prob(self) -> Tensor:
        return torch.sigmoid(self.contact_logits) * self.pair_mask

    @property
    def distance_prob(self) -> Tensor:
        return torch.softmax(self.distance_logits, -1) * self.pair_mask.unsqueeze(-1)


class GeometrySupervisionHead(nn.Module):
    """Attach source-only geometry supervision to an existing BPSF field."""

    schema = MECHANISM_SCHEMA

    def __init__(self, pair_dim: int, distance_bins=DISTANCE_BINS_ANGSTROM,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.n_distance_bins = len(tuple(distance_bins)) - 1
        self.contact = nn.Linear(pair_dim, 1, dtype=dtype)
        self.distance = nn.Linear(pair_dim, self.n_distance_bins, dtype=dtype)

    def forward(self, pair: Tensor, pair_mask: Tensor) -> DenseMechanismPrediction:
        if pair.ndim != 4 or pair_mask.shape != pair.shape[:3]:
            raise ValueError("geometry head requires [B,A,R,Dp] pair states and mask")
        return DenseMechanismPrediction(
            self.contact(pair).squeeze(-1), self.distance(pair), pair_mask)


class PairGeometryTeacher(nn.Module):
    """BPSF plus its sole source-only geometry supervision consumer."""

    def __init__(self, hidden_dim: int = 128, section_dim: int = 32,
                 pair_dim: int = 64, blocks: int = 3, latents: int = 16,
                 heads: int = 8, chunk_size: int = 16,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.trunk = BipartitePairSectionFormer(
            hidden_dim, section_dim, pair_dim, blocks, latents, heads,
            chunk_size, dtype=dtype)
        self.geometry = GeometrySupervisionHead(pair_dim, dtype=dtype)

    def forward(self, ligand_atoms: Tensor, ligand_mask: Tensor,
                protein_residues: Tensor, residue_mask: Tensor,
                adjacency: Tensor | None = None) -> DenseMechanismPrediction:
        encoded = self.trunk(
            ligand_atoms, protein_residues, ligand_mask, residue_mask,
            adjacency, return_pair=True)
        return self.geometry(encoded.pair, encoded.pair_mask)
