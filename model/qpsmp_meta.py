"""Mechanism-evidence meta-learning for cold-target affinity prediction.

The active path is solver-free and single stage.  It retains aligned
protein--ligand interaction slots, binds support residuals to slot
sensitivities, and transports that evidence to each query through a strict
support--query difference path.
"""
from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn as nn
from torch import Tensor

from .bpsf import BipartitePairSectionFormer
from .cartesian import SparseCartesianMechanismEncoder
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
class MechanismMetaState:
    support_zero: Tensor
    query_add: Tensor
    query_ligand_value: Tensor
    query_cross: Tensor
    zero_shot: Tensor
    task_prompts: Tensor
    centered_residual: Tensor
    level_shift: Tensor
    level_adjustment: Tensor
    level_gate: Tensor
    support_evidence: Tensor
    query_gate: Tensor
    reference_correction: Tensor
    support_auxiliary_loss: Tensor


class MechanismEvidenceMetaTransformer(nn.Module):
    """Convert label-bound slot evidence into query-specific scalar gates.

    Support order has no positional encoding.  Query transport scores receive
    only aligned slot differences; absolute support slots are confined to the
    task-evidence path.  The returned scalar gates can modulate scalar, vector,
    or tensor channels without changing their O(3) transformation law.
    """

    def __init__(self, slot_dim: int, prompt_dim: int, slot_count: int,
                 blocks: int = 2, dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        heads = next(h for h in (8, 4, 2, 1) if prompt_dim % h == 0)
        self.slot_count = int(slot_count)
        self.prompt_dim = int(prompt_dim)
        self.sensitivity = nn.Sequential(
            nn.LayerNorm(slot_dim, dtype=dtype),
            nn.Linear(slot_dim, prompt_dim, dtype=dtype), nn.SiLU(),
            nn.Linear(prompt_dim, 1, bias=False, dtype=dtype))
        self.evidence = nn.Sequential(
            nn.LayerNorm(slot_dim + 4, dtype=dtype),
            nn.Linear(slot_dim + 4, prompt_dim, dtype=dtype), nn.GELU())
        self.slot_identity = nn.Parameter(
            torch.empty(slot_count, prompt_dim, dtype=dtype))
        nn.init.normal_(self.slot_identity, std=prompt_dim ** -0.5)
        self.attention = nn.ModuleList(
            nn.MultiheadAttention(
                prompt_dim, heads, batch_first=True, dtype=dtype)
            for _ in range(blocks))
        self.feedforward = nn.ModuleList(
            nn.Sequential(
                nn.LayerNorm(prompt_dim, dtype=dtype),
                nn.Linear(prompt_dim, 4 * prompt_dim, dtype=dtype), nn.GELU(),
                nn.Linear(4 * prompt_dim, prompt_dim, dtype=dtype))
            for _ in range(blocks))
        self.prompt_norm = nn.LayerNorm(prompt_dim, dtype=dtype)
        self.difference = nn.Sequential(
            nn.LayerNorm(slot_dim, dtype=dtype),
            nn.Linear(slot_dim, prompt_dim, bias=False, dtype=dtype), nn.GELU(),
            nn.Linear(prompt_dim, prompt_dim, bias=False, dtype=dtype))
        self.difference_score = nn.Linear(
            prompt_dim, 1, bias=False, dtype=dtype)
        self.gate = nn.Sequential(
            nn.LayerNorm(prompt_dim, dtype=dtype),
            nn.Linear(prompt_dim, 1, bias=False, dtype=dtype))
        nn.init.normal_(self.gate[-1].weight, std=0.02)

    def forward(self, support_slots: Tensor, support_residual: Tensor,
                query_slots: Tensor, *,
                task_state_override: Tensor | None = None,
                adapt: bool = True) -> tuple[Tensor, Tensor, Tensor, Tensor]:
        if support_slots.ndim != 4 or query_slots.ndim != 4:
            raise ValueError("mechanism slot batches must have rank four")
        batch, support_count, slot_count, slot_dim = support_slots.shape
        if slot_count != self.slot_count or query_slots.shape[0] != batch \
                or query_slots.shape[2:] != (slot_count, slot_dim):
            raise ValueError("support/query mechanism slot contracts disagree")
        if support_residual.shape != (batch, support_count):
            raise ValueError("support residuals do not match mechanism slots")
        query_count = query_slots.shape[1]
        empty_state = query_slots.new_zeros(
            batch, slot_count, self.prompt_dim)
        empty_gate = query_slots.new_zeros(batch, query_count, slot_count)
        if not adapt or support_count == 0:
            if task_state_override is not None:
                raise ValueError("cannot override mechanism state without support")
            return empty_state, empty_gate, empty_gate, query_slots.new_zeros(())

        centered = support_residual - support_residual.mean(-1, keepdim=True)
        sensitivity = self.sensitivity(support_slots).squeeze(-1)
        pseudo_gradient = -support_residual.unsqueeze(-1) * sensitivity
        broadcast = lambda value: value[:, :, None, None].expand(
            -1, -1, slot_count, 1)
        evidence = self.evidence(torch.cat((
            support_slots, pseudo_gradient.unsqueeze(-1),
            broadcast(support_residual), broadcast(centered),
            broadcast(support_residual.abs())), dim=-1))
        evidence = evidence + self.slot_identity[None, None]
        token = evidence.flatten(1, 2)
        for attention, feedforward in zip(self.attention, self.feedforward):
            update, _ = attention(token, token, token, need_weights=False)
            token = token + update
            token = token + feedforward(token)
        token = self.prompt_norm(token).reshape(
            batch, support_count, slot_count, self.prompt_dim)
        task_prompts = token.mean(1)
        if task_state_override is not None:
            if task_state_override.shape != task_prompts.shape:
                raise ValueError(
                    "task_state_override must match target mechanism prompts "
                    f"{tuple(task_prompts.shape)}, got "
                    f"{tuple(task_state_override.shape)}")
            task_prompts = task_state_override.to(task_prompts)

        # Aligned slots permit a strict difference-only reference path.  No
        # absolute support/query representation enters the transport scores.
        difference = query_slots[:, :, None] - support_slots[:, None]
        relative = self.difference(difference)
        weights = torch.softmax(
            self.difference_score(relative).squeeze(-1), dim=2)
        evidence_context = torch.einsum("bqkm,bkmh->bqmh", weights, token)
        relative_context = (weights.unsqueeze(-1) * relative).sum(2)
        query_prompts = (task_prompts[:, None] + evidence_context
                         + relative_context)
        gate = torch.tanh(self.gate(query_prompts).squeeze(-1))
        # Evidence magnitude is a diagnostic, not a minimization target: a
        # direct penalty would teach the model to erase its own sensitivity.
        auxiliary = pseudo_gradient.new_zeros(())
        return task_prompts, gate, pseudo_gradient, auxiliary


class QPSMPMetaLearner(nn.Module):
    """Shared scalar heads plus difference-constrained mechanism adaptation."""

    def __init__(self, hidden_dim: int, slot_dim: int, slot_count: int,
                 task_dim: int = 32, support_blocks: int = 2,
                 mechanism_scale: float = 0.25,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.mechanism_scale = float(mechanism_scale)
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
        self.mechanism = MechanismEvidenceMetaTransformer(
            slot_dim, task_dim, slot_count, support_blocks, dtype)
        self.reference_delta = nn.Sequential(
            nn.LayerNorm(slot_dim, dtype=dtype),
            nn.Linear(slot_dim, task_dim, bias=False, dtype=dtype), nn.GELU(),
            nn.Linear(task_dim, 1, bias=False, dtype=dtype))
        nn.init.zeros_(self.reference_delta[-1].weight)
        self.reference_score = nn.Sequential(
            nn.LayerNorm(slot_dim, dtype=dtype),
            nn.Linear(slot_dim, 1, bias=False, dtype=dtype))
        self.reference_blend = nn.Sequential(
            nn.Linear(3, 16, dtype=dtype), nn.GELU(),
            nn.Linear(16, 1, dtype=dtype), nn.Sigmoid())
        nn.init.constant_(self.reference_blend[-2].bias, -2.0)
        self.level_gate = nn.Sequential(
            nn.Linear(3, 16, dtype=dtype), nn.GELU(),
            nn.Linear(16, 1, dtype=dtype), nn.Sigmoid())

    def infer(self, protein: Tensor, support_ligand: Tensor,
              support_endpoint: Tensor, support_slots: Tensor,
              support_y: Tensor, query_ligand: Tensor,
              query_endpoint: Tensor, query_slots: Tensor, *,
              adapt: bool = True,
              task_state_override: Tensor | None = None) -> MechanismMetaState:
        support_count = support_ligand.shape[1]
        protein_level = self.protein_level(protein).squeeze(-1).unsqueeze(-1)
        support_ligand_value = self.ligand_baseline(support_ligand).squeeze(-1)
        query_ligand_value = self.ligand_baseline(query_ligand).squeeze(-1)
        support_add = support_ligand_value + protein_level
        query_add = query_ligand_value + protein_level
        support_cross = self.cross_head(support_endpoint).squeeze(-1)
        query_cross = self.cross_head(query_endpoint).squeeze(-1)
        support_zero = support_add + support_cross
        zero_shot = query_add + query_cross
        if not adapt or support_count == 0:
            residual = support_y.new_zeros(support_y.shape)
            centered = residual
            level_shift = zero_shot.new_zeros(zero_shot.shape[0])
            level_gate = level_shift
            level_adjustment = torch.zeros_like(zero_shot)
        else:
            residual = support_y - support_zero
            centered = residual - residual.mean(-1, keepdim=True)
            level_shift = residual.mean(-1)
            mad = centered.abs().mean(-1)
            count = torch.full_like(level_shift, float(support_count)).log1p()
            level_gate = self.level_gate(
                torch.stack((level_shift, mad, count), -1)).squeeze(-1)
            support_mean = support_y.mean(-1, keepdim=True)
            support_zero_mean = support_zero.mean(-1, keepdim=True)
            level_baseline = support_mean + level_gate.unsqueeze(-1) * (
                zero_shot - support_zero_mean)
            level_adjustment = level_baseline - zero_shot
        task_prompts, query_gate, pseudo_gradient, auxiliary = self.mechanism(
            support_slots, residual, query_slots,
            task_state_override=task_state_override, adapt=adapt)
        evidence = (pseudo_gradient.square().mean((-2, -1)).sqrt()
                    if support_count else zero_shot.new_zeros(zero_shot.shape[0]))
        if adapt and support_count:
            difference = query_slots[:, :, None] - support_slots[:, None]
            gated_difference = difference * (
                1.0 + query_gate[:, :, None, :, None])
            # With the query-conditioned gate held fixed, the odd construction
            # makes an identity reference zero and reverses sign with the slot
            # difference.  It does not claim global support/query exchange
            # antisymmetry because the gate itself is directional.
            delta = 0.5 * (
                self.reference_delta(gated_difference)
                - self.reference_delta(-gated_difference)).squeeze(-1).mean(-1)
            reference_weight = torch.softmax(
                self.reference_score(difference).squeeze(-1).mean(-1), dim=-1)
            anchors = (support_y[:, None] + zero_shot[:, :, None]
                       - support_zero[:, None] + delta)
            reference_prediction = (reference_weight * anchors).sum(-1)
            mean_distance = torch.linalg.vector_norm(
                difference, dim=-1).mean((-2, -1))
            blend = self.reference_blend(torch.stack((
                evidence[:, None].expand_as(mean_distance), mean_distance,
                torch.full_like(mean_distance, float(support_count)).log1p()),
                -1)).squeeze(-1)
            reference_correction = blend * (
                reference_prediction - (zero_shot + level_adjustment))
            if support_count >= 2:
                support_difference = (
                    support_slots[:, :, None] - support_slots[:, None])
                task_gate = torch.tanh(
                    self.mechanism.gate(task_prompts).squeeze(-1))
                gated_support_difference = support_difference * (
                    1.0 + task_gate[:, None, None, :, None])
                support_delta = 0.5 * (
                    self.reference_delta(gated_support_difference)
                    - self.reference_delta(-gated_support_difference)
                ).squeeze(-1).mean(-1)
                support_score = self.reference_score(
                    support_difference).squeeze(-1).mean(-1)
                diagonal = torch.eye(
                    support_count, dtype=torch.bool, device=support_score.device)
                support_score = support_score.masked_fill(
                    diagonal, torch.finfo(support_score.dtype).min)
                support_weight = torch.softmax(support_score, -1)
                support_anchor = (
                    support_y[:, None] + support_zero[:, :, None]
                    - support_zero[:, None] + support_delta)
                reconstructed = (support_weight * support_anchor).sum(-1)
                auxiliary = torch.nn.functional.smooth_l1_loss(
                    reconstructed, support_y)
        else:
            reference_correction = torch.zeros_like(zero_shot)
        return MechanismMetaState(
            support_zero, query_add, query_ligand_value, query_cross, zero_shot,
            task_prompts, centered, level_shift, level_adjustment, level_gate,
            evidence, query_gate, reference_correction, auxiliary)

    @staticmethod
    def delta(predictions: Tensor, left: Tensor, right: Tensor) -> Tensor:
        return predictions.index_select(0, right) - predictions.index_select(0, left)

    @staticmethod
    def rectangle(delta_a: Tensor, delta_b: Tensor) -> Tensor:
        return delta_a - delta_b


class QPSMPBioModel(nn.Module):
    """BPSF mechanism slots with single-stage episodic meta-adaptation."""

    def __init__(self, protein_dim: int, hidden_dim: int = 128,
                 task_dim: int = 32, ligand_layers: int = 3,
                 pair_dim: int = 64, pair_blocks: int = 3,
                 pair_latents: int = 16, pair_heads: int = 8,
                 pair_chunk_size: int = 16, support_hidden_dim: int = 128,
                 support_blocks: int = 2, adapter_rank: int = 4,
                 adaptive_blocks: int = 2, adapter_scale: float = 0.25,
                 use_cartesian: bool = False,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        del support_hidden_dim, adapter_rank, adaptive_blocks
        self.protein_encoder = ProteinEncoder(protein_dim, hidden_dim, dtype)
        self.ligand_encoder = LigandEncoder(hidden_dim, ligand_layers, dtype=dtype)
        self.pair_section = BipartitePairSectionFormer(
            hidden_dim, task_dim, pair_dim, pair_blocks, pair_latents,
            pair_heads, pair_chunk_size, dtype=dtype)
        self.cartesian_encoder = (SparseCartesianMechanismEncoder(
            hidden_dim, scalar_dim=max(pair_dim // 2, 8),
            vector_channels=max(pair_dim // 4, 4),
            tensor_channels=max(pair_dim // 8, 2), layers=2,
            mechanism_slots=pair_latents, slot_dim=pair_dim, dtype=dtype)
            if use_cartesian else None)
        self.geometry_scale = (nn.Parameter(torch.tensor(0.1, dtype=dtype))
                               if use_cartesian else None)
        self.meta = QPSMPMetaLearner(
            hidden_dim, pair_dim, pair_latents, task_dim,
            support_blocks, adapter_scale, dtype)

    def _pair_encode(self, atom_states: Tensor, residues: Tensor, mask: Tensor,
                     protein_mask: Tensor, bonds: Tensor):
        adjacency = (bonds.abs().sum(-1) > 0).to(atom_states.dtype)
        return self.pair_section(
            atom_states, residues, mask, protein_mask, adjacency)

    def forward(self, protein_pooled: Tensor, protein_tokens: Tensor,
                protein_mask: Tensor, support_atoms: Tensor, support_bonds: Tensor,
                support_mask: Tensor, support_y: Tensor, query_atoms: Tensor,
                query_bonds: Tensor, query_mask: Tensor, *,
                adapt: bool = True,
                task_state_override: Tensor | None = None,
                geometry_coordinates: Tensor | None = None,
                geometry_edge_index: Tensor | None = None,
                geometry_available: Tensor | None = None,
                geometry_common_frame: Tensor | None = None) -> QPSMPMetaOutput:
        parameter = next(self.parameters())
        device, dtype = parameter.device, parameter.dtype
        unbatched = protein_pooled.ndim == 1
        if unbatched:
            protein_pooled = protein_pooled.unsqueeze(0)
            protein_tokens = protein_tokens.unsqueeze(0)
            protein_mask = protein_mask.unsqueeze(0)
            support_atoms = support_atoms.unsqueeze(0)
            support_bonds = support_bonds.unsqueeze(0)
            support_mask = support_mask.unsqueeze(0)
            support_y = support_y.unsqueeze(0)
            query_atoms = query_atoms.unsqueeze(0)
            query_bonds = query_bonds.unsqueeze(0)
            query_mask = query_mask.unsqueeze(0)
            if task_state_override is not None:
                task_state_override = task_state_override.unsqueeze(0)
            if geometry_coordinates is not None:
                geometry_coordinates = geometry_coordinates.unsqueeze(0)
            if geometry_available is not None:
                geometry_available = geometry_available.unsqueeze(0)
            if geometry_common_frame is not None:
                geometry_common_frame = geometry_common_frame.unsqueeze(0)
        pooled, residues = self.protein_encoder(
            protein_pooled.to(device, dtype), protein_tokens.to(device, dtype))
        atoms = torch.cat((support_atoms, query_atoms), 1).to(device, dtype)
        bonds = torch.cat((support_bonds, query_bonds), 1).to(device, dtype)
        mask = torch.cat((support_mask, query_mask), 1).to(device, dtype)
        batch, count, atom_count = atoms.shape[:3]
        ligand, atom_states = self.ligand_encoder(
            atoms.flatten(0, 1), bonds.flatten(0, 1), mask.flatten(0, 1))
        ligand = ligand.reshape(batch, count, -1)
        atom_states = atom_states.reshape(batch, count, atom_count, -1)
        encoded = self._pair_encode(
            atom_states.flatten(0, 1),
            residues[:, None].expand(-1, count, -1, -1).flatten(0, 1),
            mask.flatten(0, 1),
            protein_mask[:, None].expand(-1, count, -1).flatten(0, 1).to(device),
            bonds.flatten(0, 1))
        endpoint = encoded.endpoint.reshape(batch, count, -1)
        section = encoded.section.reshape(batch, count, -1)
        slots = encoded.mechanism_slots.reshape(batch, count, *encoded.mechanism_slots.shape[1:])
        if geometry_coordinates is not None or geometry_edge_index is not None:
            if self.cartesian_encoder is None:
                raise ValueError("Cartesian inputs require use_cartesian=True")
            if geometry_coordinates is None or geometry_edge_index is None:
                raise ValueError("coordinates and Cartesian edges must be provided together")
            residue_count = residues.shape[1]
            expected = (batch, count, atom_count + residue_count, 3)
            if geometry_coordinates.shape != expected:
                raise ValueError(
                    f"geometry coordinates must have shape {expected}")
            if geometry_available is None:
                geometry_available = torch.ones(
                    batch, count, dtype=torch.bool, device=device)
            if geometry_available.shape != (batch, count):
                raise ValueError("geometry availability mask has wrong shape")
            if geometry_common_frame is None \
                    or geometry_common_frame.shape != (batch, count):
                raise ValueError(
                    "Cartesian interaction requires an explicit common-frame mask")
            if bool((geometry_available.bool()
                     & ~geometry_common_frame.bool()).any()):
                raise ValueError(
                    "available Cartesian interaction coordinates must share a common frame")
            geometry_nodes = torch.cat((
                atom_states,
                residues[:, None].expand(-1, count, -1, -1)), dim=2)
            geometry_mask = torch.cat((mask, protein_mask[:, None].expand(
                -1, count, -1).to(device)), dim=2).bool()
            geometry_mask = geometry_mask & geometry_available.to(device).unsqueeze(-1)
            geometry = self.cartesian_encoder(
                geometry_nodes.flatten(0, 1), geometry_mask.flatten(0, 1),
                geometry_coordinates.to(device, dtype).flatten(0, 1),
                geometry_edge_index.to(device))
            geometry_slots = geometry.mechanism_slots.reshape(
                batch, count, *geometry.mechanism_slots.shape[1:])
            slots = slots + torch.tanh(self.geometry_scale) * geometry_slots
        elif self.cartesian_encoder is not None and geometry_available is not None \
                and bool(geometry_available.any()):
            raise ValueError("available geometry requires coordinates and edges")
        support_count = support_atoms.shape[1]
        support_ligand, query_ligand = torch.split(
            ligand, (support_count, count - support_count), 1)
        support_endpoint, query_endpoint = torch.split(
            endpoint, (support_count, count - support_count), 1)
        support_slots, query_slots = torch.split(
            slots, (support_count, count - support_count), 1)
        _, query_section = torch.split(
            section, (support_count, count - support_count), 1)
        state = self.meta.infer(
            pooled, support_ligand, support_endpoint, support_slots,
            support_y.to(device, dtype), query_ligand, query_endpoint,
            query_slots, adapt=adapt,
            task_state_override=task_state_override)
        if adapt and support_count:
            adapted_slots = query_slots * (
                1.0 + self.meta.mechanism_scale * state.query_gate.unsqueeze(-1))
            adapted_endpoint, query_section = self.pair_section.latent.project_slots(
                adapted_slots)
            adapted_cross = self.meta.cross_head(adapted_endpoint).squeeze(-1)
            sar = adapted_cross - state.query_cross
        else:
            sar = torch.zeros_like(state.query_cross)
        level_baseline = state.zero_shot + state.level_adjustment
        sar = sar + state.reference_correction
        prediction = level_baseline + sar
        evidence_score = torch.tanh(state.support_evidence)
        sar_scale = state.query_gate.abs().mean((-2, -1))
        one = torch.ones_like(state.level_gate)
        output = QPSMPMetaOutput(
            prediction, state.query_add, state.query_ligand_value,
            state.query_cross, level_baseline, state.level_adjustment, sar,
            prediction - state.zero_shot, state.zero_shot, state.task_prompts,
            state.level_shift, query_section, state.centered_residual,
            state.support_evidence, evidence_score, state.level_gate, one,
            sar_scale, state.support_auxiliary_loss)
        if not unbatched:
            return output
        return QPSMPMetaOutput(*(
            value.squeeze(0) if isinstance(value, Tensor) and value.ndim else value
            for value in output.__dict__.values()))
