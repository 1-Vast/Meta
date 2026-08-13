"""Active trainable QPSMP-BPSF model for cold-target few-shot DTA."""
from __future__ import annotations

from dataclasses import dataclass
import torch
import torch.nn as nn
from torch import Tensor

from .bpsf import BipartitePairSectionFormer, QuotientSupportSetOperator
from .encoders import LigandEncoder, ProteinEncoder


@dataclass(frozen=True)
class QPSMPMetaOutput:
    prediction: Tensor
    additive: Tensor
    ligand_only: Tensor
    cross_zero_shot: Tensor
    level_baseline: Tensor
    level_adjustment: Tensor
    sar_adaptation: Tensor
    adaptation: Tensor
    zero_shot: Tensor
    task_state: Tensor
    level_shift: Tensor
    query_basis: Tensor
    support_residual_quotient: Tensor
    support_evidence: Tensor
    evidence_score: Tensor
    level_shrinkage: Tensor
    shape_scale: Tensor
    sar_scale: Tensor


class QPSMPMetaLearner(nn.Module):
    """Learned scalar potential and quotient-preserving episodic operator."""

    def __init__(self, hidden_dim: int, section_dim: int,
                 support_hidden_dim: int = 96, support_blocks: int = 2,
                 state_bound: float = 1.0,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        heads = next(h for h in (8, 4, 2, 1) if support_hidden_dim % h == 0)
        self.ligand_baseline = nn.Sequential(
            nn.Linear(hidden_dim, 2 * hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(2 * hidden_dim, 1, dtype=dtype))
        self.protein_level = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, 1, dtype=dtype))
        self.cross_head = nn.Sequential(
            nn.LayerNorm(hidden_dim, dtype=dtype),
            nn.Linear(hidden_dim, 2 * hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(2 * hidden_dim, 1, bias=False, dtype=dtype))
        self.level_prior_raw = nn.Parameter(torch.tensor(0.5413249, dtype=dtype))
        self.level_noise_raw = nn.Parameter(torch.tensor(-0.2981850, dtype=dtype))
        self.section_operator = QuotientSupportSetOperator(
            section_dim, support_hidden_dim, heads, support_blocks,
            state_bound, dtype)

    def forward(self, protein: Tensor, support_ligand: Tensor,
                support_endpoint: Tensor, support_section: Tensor,
                support_y: Tensor, query_ligand: Tensor,
                query_endpoint: Tensor, query_section: Tensor, *,
                adapt: bool = True) -> QPSMPMetaOutput:
        batched = support_ligand.ndim == 3
        support_count = support_ligand.shape[1 if batched else 0]
        protein_level = self.protein_level(protein).squeeze(-1)
        support_ligand_value = self.ligand_baseline(support_ligand).squeeze(-1)
        query_ligand_value = self.ligand_baseline(query_ligand).squeeze(-1)
        if batched:
            protein_level = protein_level.unsqueeze(-1)
        support_add = support_ligand_value + protein_level
        query_add = query_ligand_value + protein_level
        support_cross = self.cross_head(support_endpoint).squeeze(-1)
        query_cross = self.cross_head(query_endpoint).squeeze(-1)
        support_zero = support_add + support_cross
        zero_shot = query_add + query_cross
        if not adapt or support_count == 0:
            state = torch.zeros(zero_shot.shape[:-1] + (query_section.shape[-1],),
                                device=zero_shot.device, dtype=zero_shot.dtype)
            centered = support_y.new_zeros(support_y.shape)
            level_shift = torch.zeros(zero_shot.shape[:-1], device=zero_shot.device)
            shrinkage = torch.zeros_like(level_shift)
            level_adjustment = torch.zeros_like(zero_shot)
            sar = torch.zeros_like(zero_shot)
        else:
            residual = support_y - support_zero
            centered = residual - residual.mean(-1, keepdim=True)
            tolerance = (32 * torch.finfo(centered.dtype).eps
                         * residual.abs().amax(-1, keepdim=True).clamp_min(1))
            centered = torch.where(centered.abs() <= tolerance,
                                   torch.zeros_like(centered), centered)
            level_shift = residual.mean(-1)
            prior = torch.nn.functional.softplus(self.level_prior_raw)
            noise = torch.nn.functional.softplus(self.level_noise_raw)
            shrinkage = support_count * prior / (noise + support_count * prior)
            level_adjustment = shrinkage * level_shift
            if batched:
                level_adjustment = level_adjustment.unsqueeze(-1)
            state, query_section, sar = self.section_operator(
                support_section, query_section, centered)
        level_baseline = zero_shot + level_adjustment
        prediction = level_baseline + sar
        evidence = torch.linalg.vector_norm(centered, dim=-1) / max(support_count, 1) ** 0.5
        evidence_score = torch.tanh(torch.linalg.vector_norm(state, dim=-1))
        one = torch.ones((), device=prediction.device, dtype=prediction.dtype)
        return QPSMPMetaOutput(
            prediction, query_add, query_ligand_value, query_cross,
            level_baseline, level_adjustment, sar, prediction - zero_shot,
            zero_shot, state, level_shift, query_section, centered, evidence,
            evidence_score, shrinkage, one, one)

    @staticmethod
    def delta(predictions: Tensor, left: Tensor, right: Tensor) -> Tensor:
        return predictions.index_select(0, right) - predictions.index_select(0, left)

    @staticmethod
    def rectangle(delta_a: Tensor, delta_b: Tensor) -> Tensor:
        return delta_a - delta_b


class QPSMPBioModel(nn.Module):
    """Small-to-medium BPSF model; the only active QPSMP architecture."""

    def __init__(self, protein_dim: int, hidden_dim: int = 128,
                 task_dim: int = 32, ligand_layers: int = 3,
                 pair_dim: int = 64, pair_blocks: int = 3,
                 pair_latents: int = 16, pair_heads: int = 8,
                 pair_chunk_size: int = 16, support_hidden_dim: int = 128,
                 support_blocks: int = 2, state_bound: float = 1.0,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.protein_encoder = ProteinEncoder(protein_dim, hidden_dim, dtype)
        self.ligand_encoder = LigandEncoder(hidden_dim, ligand_layers, dtype=dtype)
        self.pair_section = BipartitePairSectionFormer(
            hidden_dim, task_dim, pair_dim, pair_blocks, pair_latents,
            pair_heads, pair_chunk_size, dtype)
        self.meta = QPSMPMetaLearner(
            hidden_dim, task_dim, support_hidden_dim, support_blocks,
            state_bound, dtype)

    def forward(self, protein_pooled: Tensor, protein_tokens: Tensor,
                protein_mask: Tensor, support_atoms: Tensor, support_bonds: Tensor,
                support_mask: Tensor, support_y: Tensor, query_atoms: Tensor,
                query_bonds: Tensor, query_mask: Tensor, *,
                adapt: bool = True) -> QPSMPMetaOutput:
        parameter = next(self.parameters())
        device, dtype = parameter.device, parameter.dtype
        batched = protein_pooled.ndim == 2
        sequence_dim = 1 if batched else 0
        pooled, residues = self.protein_encoder(
            (protein_pooled if batched else protein_pooled.unsqueeze(0)).to(device, dtype),
            (protein_tokens if batched else protein_tokens.unsqueeze(0)).to(device, dtype))
        atoms = torch.cat((support_atoms, query_atoms), sequence_dim).to(device, dtype)
        bonds = torch.cat((support_bonds, query_bonds), sequence_dim).to(device, dtype)
        mask = torch.cat((support_mask, query_mask), sequence_dim).to(device, dtype)
        if batched:
            batch, count, atom_count = atoms.shape[:3]
            ligand, atom_states = self.ligand_encoder(
                atoms.flatten(0, 1), bonds.flatten(0, 1), mask.flatten(0, 1))
            ligand = ligand.reshape(batch, count, -1)
            atom_states = atom_states.reshape(batch, count, atom_count, -1)
            flat_residues = residues[:, None].expand(-1, count, -1, -1).flatten(0, 1)
            flat_residue_mask = protein_mask[:, None].expand(-1, count, -1).flatten(0, 1)
            adjacency = (bonds.abs().sum(-1) > 0).to(dtype).flatten(0, 1)
            encoded = self.pair_section(
                atom_states.flatten(0, 1), flat_residues, mask.flatten(0, 1),
                flat_residue_mask.to(device), adjacency)
            endpoint = encoded.endpoint.reshape(batch, count, -1)
            section = encoded.section.reshape(batch, count, -1)
        else:
            ligand, atom_states = self.ligand_encoder(atoms, bonds, mask)
            count = atoms.shape[0]
            adjacency = (bonds.abs().sum(-1) > 0).to(dtype)
            encoded = self.pair_section(
                atom_states, residues.squeeze(0).unsqueeze(0).expand(count, -1, -1),
                mask, protein_mask.to(device).unsqueeze(0).expand(count, -1), adjacency)
            endpoint, section = encoded.endpoint, encoded.section
        support_count = support_atoms.shape[sequence_dim]
        support_ligand, query_ligand = torch.split(
            ligand, (support_count, query_atoms.shape[sequence_dim]), sequence_dim)
        support_endpoint, query_endpoint = torch.split(
            endpoint, (support_count, query_atoms.shape[sequence_dim]), sequence_dim)
        support_section, query_section = torch.split(
            section, (support_count, query_atoms.shape[sequence_dim]), sequence_dim)
        return self.meta(
            pooled if batched else pooled.squeeze(0), support_ligand,
            support_endpoint, support_section, support_y.to(device, dtype),
            query_ligand, query_endpoint, query_section, adapt=adapt)
