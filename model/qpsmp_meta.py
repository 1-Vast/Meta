"""QPSMP-HyperSAR for cold-target few-shot drug-target affinity prediction."""
from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn as nn
from torch import Tensor

from .bpsf import BipartitePairSectionFormer
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
    support_match_loss: Tensor


@dataclass(frozen=True)
class HyperSARState:
    support_zero: Tensor
    query_add: Tensor
    query_ligand_value: Tensor
    query_cross: Tensor
    zero_shot: Tensor
    task_code: Tensor
    centered_residual: Tensor
    level_shift: Tensor
    level_adjustment: Tensor
    level_gate: Tensor
    support_evidence: Tensor


class AmortizedTargetConditioner(nn.Module):
    """Permutation-invariant support-to-code inference without a linear solve."""

    def __init__(self, interaction_dim: int, ligand_dim: int, task_dim: int = 32,
                 hidden_dim: int = 128, blocks: int = 2,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        heads = next(h for h in (8, 4, 2, 1) if hidden_dim % h == 0)
        self.token = nn.Sequential(
            nn.Linear(interaction_dim + ligand_dim + 2, hidden_dim, dtype=dtype),
            nn.GELU(), nn.LayerNorm(hidden_dim, dtype=dtype))
        self.attention = nn.ModuleList(
            nn.MultiheadAttention(hidden_dim, heads, batch_first=True, dtype=dtype)
            for _ in range(blocks))
        self.ffn = nn.ModuleList(
            nn.Sequential(
                nn.LayerNorm(hidden_dim, dtype=dtype),
                nn.Linear(hidden_dim, 4 * hidden_dim, dtype=dtype), nn.GELU(),
                nn.Linear(4 * hidden_dim, hidden_dim, dtype=dtype))
            for _ in range(blocks))
        self.pool_seed = nn.Parameter(torch.empty(1, hidden_dim, dtype=dtype))
        nn.init.normal_(self.pool_seed, std=hidden_dim ** -0.5)
        self.pool = nn.MultiheadAttention(
            hidden_dim, heads, batch_first=True, dtype=dtype)
        self.output = nn.Sequential(
            nn.LayerNorm(hidden_dim, dtype=dtype),
            nn.Linear(hidden_dim, task_dim, dtype=dtype), nn.Tanh())
        self.binding_token = nn.Sequential(
            nn.Linear(interaction_dim + ligand_dim, hidden_dim, dtype=dtype),
            nn.GELU(), nn.LayerNorm(hidden_dim, dtype=dtype))
        self.binding_output = nn.Linear(
            hidden_dim, task_dim, bias=False, dtype=dtype)

    def _encode(self, interaction: Tensor, ligand: Tensor,
                centered_residual: Tensor) -> Tensor:
        token = self.token(torch.cat((
            interaction, ligand, centered_residual.unsqueeze(-1),
            centered_residual.abs().unsqueeze(-1)), dim=-1))
        for attention, ffn in zip(self.attention, self.ffn):
            update, _ = attention(token, token, token, need_weights=False)
            token = token + update
            token = token + ffn(token)
        seed = self.pool_seed.unsqueeze(0).expand(token.shape[0], -1, -1)
        pooled, _ = self.pool(seed, token, token, need_weights=False)
        return self.output(pooled.squeeze(1))

    def forward(self, interaction: Tensor, ligand: Tensor,
                centered_residual: Tensor) -> Tensor:
        unbatched = interaction.ndim == 2
        if unbatched:
            interaction = interaction.unsqueeze(0)
            ligand = ligand.unsqueeze(0)
            centered_residual = centered_residual.unsqueeze(0)
        pooled_code = self._encode(interaction, ligand, centered_residual)
        binding_token = self.binding_token(torch.cat((interaction, ligand), -1))
        denominator = centered_residual.abs().sum(-1, keepdim=True).clamp_min(
            torch.finfo(centered_residual.dtype).eps)
        binding = (binding_token * centered_residual.unsqueeze(-1)).sum(-2)
        binding = binding / denominator
        # This moment structurally binds every residual to its support ligand
        # and interaction. Joint support permutations preserve it; label-only
        # permutations generally change it.
        code = torch.tanh(pooled_code + self.binding_output(binding))
        evidence = torch.linalg.vector_norm(centered_residual, dim=-1, keepdim=True)
        code = code * torch.tanh(evidence)
        return code.squeeze(0) if unbatched else code


class SiameseRelativeConditioner(nn.Module):
    """Turn one support reference code into query-specific modulation codes."""

    def __init__(self, interaction_dim: int, ligand_dim: int, task_dim: int,
                 hidden_dim: int, dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.query = nn.Sequential(
            nn.Linear(interaction_dim + ligand_dim, hidden_dim, dtype=dtype),
            nn.GELU(), nn.LayerNorm(hidden_dim, dtype=dtype))
        self.reference = nn.Linear(task_dim, hidden_dim, dtype=dtype)
        self.relative = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, task_dim, bias=False, dtype=dtype), nn.Tanh())

    def forward(self, reference_code: Tensor, query_interaction: Tensor,
                query_ligand: Tensor) -> Tensor:
        if query_interaction.ndim == 3:
            reference = reference_code.unsqueeze(1).expand(
                -1, query_interaction.shape[1], -1)
        else:
            reference = reference_code.unsqueeze(0).expand(
                query_interaction.shape[0], -1)
        query = self.query(torch.cat((query_interaction, query_ligand), dim=-1))
        anchor = self.reference(reference)
        gate = self.relative(torch.cat((query - anchor, query * anchor), dim=-1))
        # Exact zero preservation is structural, not a learned convention.
        return reference * (1.0 + gate)


class QPSMPMetaLearner(nn.Module):
    """Shared scalar heads and amortized target conditioner for HyperSAR."""

    def __init__(self, hidden_dim: int, task_dim: int = 32,
                 support_hidden_dim: int = 128, support_blocks: int = 2,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.task_dim = int(task_dim)
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
        self.conditioner = AmortizedTargetConditioner(
            hidden_dim, hidden_dim, task_dim, support_hidden_dim,
            support_blocks, dtype)
        self.relative_conditioner = SiameseRelativeConditioner(
            hidden_dim, hidden_dim, task_dim, support_hidden_dim, dtype)
        self.code_to_endpoint = nn.Linear(
            task_dim, hidden_dim, bias=False, dtype=dtype)
        self.match_query = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim, dtype=dtype),
            nn.GELU(), nn.LayerNorm(hidden_dim, dtype=dtype))
        self.match_support = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim, dtype=dtype),
            nn.GELU(), nn.LayerNorm(hidden_dim, dtype=dtype))
        self.match_temperature = hidden_dim ** -0.5
        self.match_logit_gain = nn.Parameter(
            torch.tensor(-2.0, dtype=dtype))
        self.reliability_scale = nn.Parameter(
            torch.tensor(2.0, dtype=dtype))
        self.reliability_bias = nn.Parameter(
            torch.tensor(-1.5, dtype=dtype))
        self.level_gate = nn.Sequential(
            nn.Linear(3, 16, dtype=dtype), nn.GELU(),
            nn.Linear(16, 1, dtype=dtype), nn.Sigmoid())

    def condition_queries(self, task_code: Tensor, query_ligand: Tensor,
                          query_endpoint: Tensor) -> Tensor:
        return self.relative_conditioner(
            task_code, query_endpoint, query_ligand)

    def inject_query_code(self, endpoint: Tensor, query_code: Tensor) -> Tensor:
        modulation = self.code_to_endpoint(query_code)
        return endpoint + 0.25 * modulation * (1.0 + torch.tanh(endpoint))

    def relative_residual_match(
            self, support_endpoint: Tensor, support_ligand: Tensor,
            centered_residual: Tensor, query_endpoint: Tensor,
            query_ligand: Tensor, reliability: Tensor | None = None) -> Tensor:
        """Neural reference-query matching over label-bound support residuals."""
        support = self.match_support(torch.cat((support_endpoint, support_ligand), -1))
        query = self.match_query(torch.cat((query_endpoint, query_ligand), -1))
        scores = torch.matmul(query, support.transpose(-1, -2))
        weights = torch.softmax(scores * self.match_temperature, dim=-1)
        gain = torch.sigmoid(self.match_logit_gain)
        correction = gain * torch.matmul(
            weights, centered_residual.unsqueeze(-1)).squeeze(-1)
        if reliability is None:
            reliability = self.support_match_reliability(
                support_endpoint, support_ligand, centered_residual)
        if correction.ndim > reliability.ndim:
            reliability = reliability.unsqueeze(-1)
        return reliability * correction

    def support_match_reliability(
            self, support_endpoint: Tensor, support_ligand: Tensor,
            centered_residual: Tensor) -> Tensor:
        """Support-only leave-one-out gate for the relative SAR hypothesis."""
        count = centered_residual.shape[-1]
        if count < 2:
            return centered_residual.new_zeros(centered_residual.shape[:-1])
        support = self.match_support(torch.cat((support_endpoint, support_ligand), -1))
        scores = torch.matmul(support, support.transpose(-1, -2))
        diagonal = torch.eye(count, device=scores.device, dtype=torch.bool)
        scores = scores.masked_fill(diagonal, torch.finfo(scores.dtype).min)
        weights = torch.softmax(scores * self.match_temperature, dim=-1)
        prediction = torch.matmul(
            weights, centered_residual.unsqueeze(-1)).squeeze(-1)
        null_error = centered_residual.square().mean(-1)
        match_error = (prediction - centered_residual).square().mean(-1)
        improvement = (null_error - match_error) / null_error.clamp_min(
            32 * torch.finfo(null_error.dtype).eps)
        logit = (torch.nn.functional.softplus(self.reliability_scale)
                 * improvement + self.reliability_bias)
        return torch.sigmoid(logit)

    def support_match_loss(
            self, support_endpoint: Tensor, support_ligand: Tensor,
            centered_residual: Tensor) -> Tensor:
        """Leave-one-out support objective for the shared matching metric."""
        count = centered_residual.shape[-1]
        if count < 2:
            return centered_residual.new_zeros(())
        support = self.match_support(torch.cat((support_endpoint, support_ligand), -1))
        scores = torch.matmul(support, support.transpose(-1, -2))
        diagonal = torch.eye(count, device=scores.device, dtype=torch.bool)
        scores = scores.masked_fill(diagonal, torch.finfo(scores.dtype).min)
        prediction = torch.matmul(
            torch.softmax(scores * self.match_temperature, dim=-1),
            centered_residual.unsqueeze(-1)).squeeze(-1)
        return (prediction - centered_residual).square().mean()

    def infer(self, protein: Tensor, support_ligand: Tensor,
              support_endpoint: Tensor, support_y: Tensor,
              query_ligand: Tensor, query_endpoint: Tensor, *,
              adapt: bool = True,
              task_state_override: Tensor | None = None) -> HyperSARState:
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
        state_shape = zero_shot.shape[:-1] + (self.task_dim,)
        if not adapt or support_count == 0:
            centered = support_y.new_zeros(support_y.shape)
            level_shift = torch.zeros(zero_shot.shape[:-1], device=zero_shot.device)
            gate = torch.zeros_like(level_shift)
            level_adjustment = torch.zeros_like(zero_shot)
            task_code = torch.zeros(
                state_shape, device=zero_shot.device, dtype=zero_shot.dtype)
        else:
            residual = support_y - support_zero
            centered = residual - residual.mean(-1, keepdim=True)
            tolerance = (32 * torch.finfo(centered.dtype).eps
                         * residual.abs().amax(-1, keepdim=True).clamp_min(1))
            centered = torch.where(
                centered.abs() <= tolerance, torch.zeros_like(centered), centered)
            level_shift = residual.mean(-1)
            mad = (residual - residual.mean(-1, keepdim=True)).abs().mean(-1)
            count = torch.full_like(level_shift, float(support_count)).log1p()
            gate = self.level_gate(torch.stack((level_shift, mad, count), -1)).squeeze(-1)
            support_mean = support_y.mean(-1)
            support_zero_mean = support_zero.mean(-1)
            if batched:
                support_mean = support_mean.unsqueeze(-1)
                support_zero_mean = support_zero_mean.unsqueeze(-1)
                shape_gate = gate.unsqueeze(-1)
            else:
                shape_gate = gate
            # Anchor the absolute level to observed support. The learned gate
            # controls only transfer of zero-shot relative query shape, so a
            # zero gate is exactly the robust support-mean baseline.
            level_baseline = support_mean + shape_gate * (
                zero_shot - support_zero_mean)
            level_adjustment = level_baseline - zero_shot
            task_code = self.conditioner(
                support_endpoint, support_ligand, centered)
            if task_state_override is not None:
                if task_state_override.shape != task_code.shape:
                    raise ValueError(
                        "task_state_override must match the target-level "
                        f"reference code shape {tuple(task_code.shape)}, got "
                        f"{tuple(task_state_override.shape)}")
                task_code = task_state_override.to(task_code)
        evidence = torch.linalg.vector_norm(centered, dim=-1) / max(support_count, 1) ** 0.5
        return HyperSARState(
            support_zero, query_add, query_ligand_value, query_cross, zero_shot,
            task_code, centered, level_shift, level_adjustment, gate, evidence)

    @staticmethod
    def delta(predictions: Tensor, left: Tensor, right: Tensor) -> Tensor:
        return predictions.index_select(0, right) - predictions.index_select(0, left)

    @staticmethod
    def rectangle(delta_a: Tensor, delta_b: Tensor) -> Tensor:
        return delta_a - delta_b


class QPSMPBioModel(nn.Module):
    """Localized pair trunk with support-generated low-rank interaction adapters."""

    def __init__(self, protein_dim: int, hidden_dim: int = 128,
                 task_dim: int = 32, ligand_layers: int = 3,
                 pair_dim: int = 64, pair_blocks: int = 3,
                 pair_latents: int = 16, pair_heads: int = 8,
                 pair_chunk_size: int = 16, support_hidden_dim: int = 128,
                 support_blocks: int = 2, adapter_rank: int = 4,
                 adaptive_blocks: int = 2, adapter_scale: float = 0.25,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.protein_encoder = ProteinEncoder(protein_dim, hidden_dim, dtype)
        self.ligand_encoder = LigandEncoder(hidden_dim, ligand_layers, dtype=dtype)
        self.pair_section = BipartitePairSectionFormer(
            hidden_dim, task_dim, pair_dim, pair_blocks, pair_latents,
            pair_heads, pair_chunk_size, task_dim, adapter_rank,
            min(adaptive_blocks, pair_blocks), adapter_scale, dtype)
        self.meta = QPSMPMetaLearner(
            hidden_dim, task_dim, support_hidden_dim, support_blocks, dtype)

    def _pair_encode(self, atom_states: Tensor, residues: Tensor, mask: Tensor,
                     protein_mask: Tensor, bonds: Tensor,
                     task_code: Tensor | None = None):
        adjacency = (bonds.abs().sum(-1) > 0).to(atom_states.dtype)
        return self.pair_section(
            atom_states, residues, mask, protein_mask, adjacency,
            task_code=task_code)

    def forward(self, protein_pooled: Tensor, protein_tokens: Tensor,
                protein_mask: Tensor, support_atoms: Tensor, support_bonds: Tensor,
                support_mask: Tensor, support_y: Tensor, query_atoms: Tensor,
                query_bonds: Tensor, query_mask: Tensor, *,
                adapt: bool = True,
                task_state_override: Tensor | None = None) -> QPSMPMetaOutput:
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
            flat_protein_mask = protein_mask[:, None].expand(-1, count, -1).flatten(0, 1)
            encoded = self._pair_encode(
                atom_states.flatten(0, 1), flat_residues, mask.flatten(0, 1),
                flat_protein_mask.to(device), bonds.flatten(0, 1))
            endpoint = encoded.endpoint.reshape(batch, count, -1)
            section = encoded.section.reshape(batch, count, -1)
        else:
            ligand, atom_states = self.ligand_encoder(atoms, bonds, mask)
            count = atoms.shape[0]
            expanded_residues = residues.squeeze(0).unsqueeze(0).expand(count, -1, -1)
            expanded_protein_mask = protein_mask.to(device).unsqueeze(0).expand(count, -1)
            encoded = self._pair_encode(
                atom_states, expanded_residues, mask, expanded_protein_mask, bonds)
            endpoint, section = encoded.endpoint, encoded.section
        support_count = support_atoms.shape[sequence_dim]
        query_count = query_atoms.shape[sequence_dim]
        support_ligand, query_ligand = torch.split(
            ligand, (support_count, query_count), sequence_dim)
        support_endpoint, query_endpoint = torch.split(
            endpoint, (support_count, query_count), sequence_dim)
        _, query_section = torch.split(section, (support_count, query_count), sequence_dim)
        state = self.meta.infer(
            pooled if batched else pooled.squeeze(0), support_ligand,
            support_endpoint, support_y.to(device, dtype), query_ligand,
            query_endpoint, adapt=adapt, task_state_override=task_state_override)

        if adapt and support_count:
            sar_reliability = self.meta.support_match_reliability(
                support_endpoint, support_ligand, state.centered_residual)
            query_code = self.meta.condition_queries(
                state.task_code, query_ligand, query_endpoint)
            if batched:
                query_code = query_code * sar_reliability[:, None, None]
            else:
                query_code = query_code * sar_reliability
            relative_code = query_code
            if batched:
                query_atom_states = atom_states[:, support_count:].flatten(0, 1)
                query_residues = residues[:, None].expand(
                    -1, query_count, -1, -1).flatten(0, 1)
                query_protein_mask = protein_mask[:, None].expand(
                    -1, query_count, -1).flatten(0, 1).to(device)
                query_code = query_code.flatten(0, 1)
                adapted = self._pair_encode(
                    query_atom_states, query_residues,
                    mask[:, support_count:].flatten(0, 1), query_protein_mask,
                    bonds[:, support_count:].flatten(0, 1), query_code)
                adapted_endpoint = adapted.endpoint.reshape(batch, query_count, -1)
                query_section = adapted.section.reshape(batch, query_count, -1)
            else:
                adapted = self._pair_encode(
                    atom_states[support_count:], expanded_residues[support_count:],
                    mask[support_count:], expanded_protein_mask[support_count:],
                    bonds[support_count:], query_code)
                adapted_endpoint, query_section = adapted.endpoint, adapted.section
            adapted_endpoint = self.meta.inject_query_code(
                adapted_endpoint, relative_code)
            adapted_cross = self.meta.cross_head(adapted_endpoint).squeeze(-1)
            active = torch.linalg.vector_norm(state.task_code, dim=-1) > 0
            if batched:
                active = active.unsqueeze(-1)
            adapted_cross = torch.where(active, adapted_cross, state.query_cross)
        else:
            adapted_cross = state.query_cross
            sar_reliability = state.zero_shot.new_zeros(
                state.zero_shot.shape[:-1])
        adapted_zero = state.query_add + adapted_cross
        matched_sar = self.meta.relative_residual_match(
            support_endpoint, support_ligand, state.centered_residual,
            query_endpoint, query_ligand, sar_reliability
        ) if adapt and support_count else 0.0
        sar = adapted_zero - state.zero_shot + matched_sar
        level_baseline = state.zero_shot + state.level_adjustment
        prediction = state.zero_shot + state.level_adjustment + sar
        evidence_score = torch.tanh(torch.linalg.vector_norm(state.task_code, dim=-1))
        support_match_loss = self.meta.support_match_loss(
            support_endpoint, support_ligand, state.centered_residual)
        one = torch.ones((), device=prediction.device, dtype=prediction.dtype)
        return QPSMPMetaOutput(
            prediction, state.query_add, state.query_ligand_value,
            state.query_cross, level_baseline, state.level_adjustment, sar,
            prediction - state.zero_shot, state.zero_shot, state.task_code,
            state.level_shift, query_section, state.centered_residual,
            state.support_evidence, evidence_score, state.level_gate, one,
            sar_reliability.mean(), support_match_loss)
