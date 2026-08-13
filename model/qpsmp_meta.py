"""CIPF + TERM for cold-target zero/few-shot affinity prediction.

CIPF exposes globally indexed protein--ligand interaction primitives from
sequence/residue embeddings and a 2D ligand graph.  TERM uses the exact loss
gradient with respect to virtual primitive coefficients; it performs no
closed-form solve, inner-loop update, or support-label copying.
"""
from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn as nn
import torch.nn.functional as F
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
class TERMState:
    support_zero: Tensor
    query_add: Tensor
    query_ligand_value: Tensor
    query_cross: Tensor
    zero_shot: Tensor
    task_evidence: Tensor
    residual: Tensor
    level_shift: Tensor
    level_adjustment: Tensor
    level_gate: Tensor
    exact_gradient: Tensor
    coefficients: Tensor
    reliability: Tensor
    primitive_regularizer: Tensor


class TriadicEvidenceRouter(nn.Module):
    """Route exact primitive-gradient evidence through P--Li--Lq triangles."""

    def __init__(self, ligand_dim: int, primitive_count: int,
                 hidden_dim: int = 32,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.primitive_count = int(primitive_count)
        self.hidden_dim = int(hidden_dim)
        self.primitive_identity = nn.Parameter(
            torch.empty(primitive_count, hidden_dim, dtype=dtype))
        nn.init.normal_(self.primitive_identity, std=hidden_dim ** -0.5)
        self.protein_prior = nn.Linear(
            ligand_dim, primitive_count, dtype=dtype)
        self.ligand_change = nn.Sequential(
            nn.LayerNorm(4 * ligand_dim, dtype=dtype),
            nn.Linear(4 * ligand_dim, hidden_dim, dtype=dtype), nn.GELU())
        self.support_evidence = nn.Sequential(
            nn.Linear(4, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim, dtype=dtype))
        self.triad = nn.Sequential(
            nn.LayerNorm(4 * hidden_dim + 2, dtype=dtype),
            nn.Linear(4 * hidden_dim + 2, 2 * hidden_dim, dtype=dtype),
            nn.GELU(), nn.Linear(2 * hidden_dim, 1, bias=False, dtype=dtype))
        self.confidence = nn.Sequential(
            nn.LayerNorm(2 * hidden_dim + 2, dtype=dtype),
            nn.Linear(2 * hidden_dim + 2, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, 1, bias=False, dtype=dtype))
        self.log_temperature = nn.Parameter(torch.zeros((), dtype=dtype))
        self.reliability_bias = nn.Parameter(torch.tensor(0.0, dtype=dtype))
        self.reliability_entropy = nn.Parameter(torch.tensor(1.0, dtype=dtype))
        self.reliability_count = nn.Parameter(torch.tensor(0.5, dtype=dtype))

    @staticmethod
    def exact_coefficient_gradient(residual: Tensor,
                                   primitive: Tensor) -> Tensor:
        """d[.5*(y-(level+a*phi))^2]/da at a=0."""
        if residual.shape != primitive.shape[:-1]:
            raise ValueError("residual and primitive responses disagree")
        return -residual.unsqueeze(-1) * primitive

    def forward(self, protein: Tensor, support_ligand: Tensor,
                support_primitive: Tensor, support_residual: Tensor,
                query_ligand: Tensor, query_primitive: Tensor, *,
                adapt: bool = True,
                task_state_override: Tensor | None = None
                ) -> tuple[Tensor, Tensor, Tensor, Tensor]:
        batch, support_count, primitive_count = support_primitive.shape
        if primitive_count != self.primitive_count:
            raise ValueError("primitive dictionary size changed across the episode")
        if query_primitive.shape[:1] + query_primitive.shape[2:] != (
                batch, primitive_count):
            raise ValueError("support/query primitive dictionaries disagree")
        query_count = query_primitive.shape[1]
        zero_coeff = query_primitive.new_zeros(batch, query_count, primitive_count)
        zero_state = query_primitive.new_zeros(
            batch, primitive_count, self.hidden_dim)
        zero_reliability = query_primitive.new_zeros(batch, query_count)
        if not adapt or support_count == 0:
            if task_state_override is not None:
                raise ValueError("cannot override TERM evidence without support")
            return zero_state, zero_coeff, zero_reliability, \
                query_primitive.new_zeros(batch, support_count, primitive_count)

        gradient = self.exact_coefficient_gradient(
            support_residual, support_primitive)
        centered = support_residual - support_residual.mean(-1, keepdim=True)
        evidence_input = torch.stack((
            gradient, support_primitive,
            support_residual.unsqueeze(-1).expand_as(gradient),
            centered.unsqueeze(-1).expand_as(gradient)), -1)
        evidence = self.support_evidence(evidence_input)
        task_evidence = evidence.mean(1) + self.primitive_identity[None]
        if task_state_override is not None:
            if task_state_override.shape != task_evidence.shape:
                raise ValueError(
                    f"task_state_override must match TERM evidence "
                    f"{tuple(task_evidence.shape)}, got "
                    f"{tuple(task_state_override.shape)}")
            task_evidence = task_state_override.to(task_evidence)

        support = support_ligand[:, None].expand(-1, query_count, -1, -1)
        query = query_ligand[:, :, None].expand(-1, -1, support_count, -1)
        change = self.ligand_change(torch.cat((
            support, query, query - support, query * support), -1))
        support_evidence = evidence[:, None].expand(-1, query_count, -1, -1, -1)
        task = task_evidence[:, None, None].expand(
            -1, query_count, support_count, -1, -1)
        change = change[:, :, :, None].expand(-1, -1, -1, primitive_count, -1)
        identity = self.primitive_identity[None, None, None].expand(
            batch, query_count, support_count, -1, -1)
        scalar = torch.stack((
            support_primitive[:, None].expand(-1, query_count, -1, -1),
            query_primitive[:, :, None].expand(-1, -1, support_count, -1)), -1)
        routed = self.triad(torch.cat((
            support_evidence, task, change, identity, scalar), -1)).squeeze(-1)
        # Set aggregation with stable evidence scale: sum/sqrt(k), not
        # mean/sqrt(k), which would attenuate repeated evidence as k grows.
        coefficients = routed.sum(2) / support_count ** 0.5

        prior = self.protein_prior(protein)
        task_for_confidence = task_evidence[:, None].expand(
            -1, query_count, -1, -1)
        identity_for_confidence = self.primitive_identity[None, None].expand(
            batch, query_count, -1, -1)
        confidence_logits = self.confidence(torch.cat((
            task_for_confidence, identity_for_confidence,
            prior[:, None, :, None].expand(-1, query_count, -1, -1),
            coefficients.unsqueeze(-1)), -1)).squeeze(-1)
        temperature = self.log_temperature.exp().clamp(0.25, 4.0)
        probability = torch.softmax(confidence_logits / temperature, -1)
        entropy = -(probability * probability.clamp_min(1e-8).log()).sum(-1)
        entropy = entropy / torch.tensor(
            float(primitive_count), device=entropy.device,
            dtype=entropy.dtype).log()
        reliability = torch.sigmoid(
            self.reliability_bias
            - F.softplus(self.reliability_entropy) * entropy.detach()
            + F.softplus(self.reliability_count)
            * torch.full_like(entropy, float(support_count)).log1p())
        # Identity and ligand-change paths may route evidence but may never
        # create an adaptation when the label-bound coefficient score is zero.
        evidence_norm = gradient.square().mean((1, 2)).sqrt()
        evidence_gate = evidence_norm / (evidence_norm + 1.0)
        reliability = reliability * evidence_gate[:, None]
        return task_evidence, torch.tanh(coefficients), reliability, gradient


# Explicit alias for downstream code and reports.
TERM = TriadicEvidenceRouter


class QPSMPMetaLearner(nn.Module):
    """Scalar baseline, support level calibration, and TERM composition."""

    def __init__(self, hidden_dim: int, primitive_count: int,
                 task_dim: int = 32, mechanism_scale: float = 0.25,
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
        self.level_gate = nn.Sequential(
            nn.Linear(3, 16, dtype=dtype), nn.GELU(),
            nn.Linear(16, 1, dtype=dtype), nn.Sigmoid())
        self.term = TriadicEvidenceRouter(
            hidden_dim, primitive_count, task_dim, dtype)

    @staticmethod
    def primitive_regularizer(query_primitive: Tensor) -> Tensor:
        centered = query_primitive - query_primitive.mean(1, keepdim=True)
        mean_loss = query_primitive.mean(1).square().mean()
        count = max(query_primitive.shape[1] - 1, 1)
        covariance = torch.einsum("bqm,bqn->bmn", centered, centered) / count
        scale = covariance.diagonal(dim1=-2, dim2=-1).clamp_min(1e-6).sqrt()
        correlation = covariance / (scale[:, :, None] * scale[:, None, :])
        identity = torch.eye(
            correlation.shape[-1], device=correlation.device,
            dtype=correlation.dtype)
        return mean_loss + (correlation - identity).square().mean()

    def infer(self, protein: Tensor, support_ligand: Tensor,
              support_endpoint: Tensor, support_primitive: Tensor,
              support_y: Tensor, query_ligand: Tensor,
              query_endpoint: Tensor, query_primitive: Tensor, *,
              adapt: bool = True,
              task_state_override: Tensor | None = None) -> TERMState:
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
            level_shift = zero_shot.new_zeros(zero_shot.shape[0])
            level_gate = level_shift
            level_adjustment = torch.zeros_like(zero_shot)
        else:
            # Raw zero-shot residual preserves one-shot evidence. The separate
            # primitive mean penalty discourages TERM from becoming a level head.
            residual = support_y - support_zero
            centered = residual - residual.mean(-1, keepdim=True)
            level_shift = residual.mean(-1)
            count = torch.full_like(level_shift, float(support_count)).log1p()
            level_gate = self.level_gate(torch.stack((
                level_shift, centered.abs().mean(-1), count), -1)).squeeze(-1)
            # A scalar shift cannot alter within-target ligand ordering. A
            # learned shrinkage leaves non-zero k=1 evidence for TERM.
            level_adjustment = level_gate.unsqueeze(-1) * level_shift.unsqueeze(-1)
        task, coefficients, reliability, gradient = self.term(
            protein, support_ligand, support_primitive, residual.detach(),
            query_ligand, query_primitive, adapt=adapt,
            task_state_override=task_state_override)
        return TERMState(
            support_zero, query_add, query_ligand_value, query_cross, zero_shot,
            task, residual, level_shift, level_adjustment, level_gate, gradient,
            coefficients, reliability,
            self.primitive_regularizer(query_primitive))

    @staticmethod
    def delta(predictions: Tensor, left: Tensor, right: Tensor) -> Tensor:
        return predictions.index_select(0, right) - predictions.index_select(0, left)

    @staticmethod
    def rectangle(delta_a: Tensor, delta_b: Tensor) -> Tensor:
        return delta_a - delta_b


class QPSMPBioModel(nn.Module):
    """Sequence/2D CIPF with solver-free TERM few-shot composition."""

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
        del support_hidden_dim, support_blocks, adapter_rank, adaptive_blocks
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
            hidden_dim, pair_latents, task_dim, adapter_scale, dtype)

    def _pair_encode(self, atom_states: Tensor, residues: Tensor, mask: Tensor,
                     protein_mask: Tensor, bonds: Tensor,
                     atom_features: Tensor):
        adjacency = (bonds.abs().sum(-1) > 0).to(atom_states.dtype)
        return self.pair_section(
            atom_states, residues, mask, protein_mask, adjacency,
            atom_features=atom_features)

    def forward(self, protein_pooled: Tensor, protein_tokens: Tensor,
                protein_mask: Tensor, support_atoms: Tensor, support_bonds: Tensor,
                support_mask: Tensor, support_y: Tensor, query_atoms: Tensor,
                query_bonds: Tensor, query_mask: Tensor, *, adapt: bool = True,
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
        raw_atoms = torch.cat((support_atoms, query_atoms), 1).to(device, dtype)
        bonds = torch.cat((support_bonds, query_bonds), 1).to(device, dtype)
        mask = torch.cat((support_mask, query_mask), 1).to(device, dtype)
        batch, count, atom_count = raw_atoms.shape[:3]
        ligand, atom_states = self.ligand_encoder(
            raw_atoms.flatten(0, 1), bonds.flatten(0, 1), mask.flatten(0, 1))
        ligand = ligand.reshape(batch, count, -1)
        atom_states = atom_states.reshape(batch, count, atom_count, -1)
        encoded = self._pair_encode(
            atom_states.flatten(0, 1),
            residues[:, None].expand(-1, count, -1, -1).flatten(0, 1),
            mask.flatten(0, 1),
            protein_mask[:, None].expand(-1, count, -1).flatten(0, 1).to(device),
            bonds.flatten(0, 1), raw_atoms.flatten(0, 1))
        endpoint = encoded.endpoint.reshape(batch, count, -1)
        section = encoded.section.reshape(batch, count, -1)
        primitive = encoded.mechanism_response.reshape(batch, count, -1)
        if geometry_coordinates is not None or geometry_edge_index is not None:
            if self.cartesian_encoder is None:
                raise ValueError("Cartesian inputs require use_cartesian=True")
            if geometry_coordinates is None or geometry_edge_index is None:
                raise ValueError("coordinates and Cartesian edges must be provided together")
            residue_count = residues.shape[1]
            expected = (batch, count, atom_count + residue_count, 3)
            if geometry_coordinates.shape != expected:
                raise ValueError(f"geometry coordinates must have shape {expected}")
            if geometry_available is None:
                geometry_available = torch.ones(
                    batch, count, dtype=torch.bool, device=device)
            if geometry_available.shape != (batch, count):
                raise ValueError("geometry availability mask has wrong shape")
            if geometry_common_frame is None or geometry_common_frame.shape != (batch, count):
                raise ValueError("Cartesian interaction requires an explicit common-frame mask")
            if bool((geometry_available.bool() & ~geometry_common_frame.bool()).any()):
                raise ValueError("available Cartesian interaction coordinates must share a common frame")
            geometry_nodes = torch.cat((
                atom_states,
                residues[:, None].expand(-1, count, -1, -1)), 2)
            geometry_mask = torch.cat((
                mask, protein_mask[:, None].expand(-1, count, -1).to(device)), 2
                ).bool() & geometry_available.to(device).unsqueeze(-1)
            geometry = self.cartesian_encoder(
                geometry_nodes.flatten(0, 1), geometry_mask.flatten(0, 1),
                geometry_coordinates.to(device, dtype).flatten(0, 1),
                geometry_edge_index.to(device))
            geometry_response = self.pair_section.latent.interaction.primitive_response(
                geometry.mechanism_slots).reshape(batch, count, -1)
            primitive = primitive + torch.tanh(self.geometry_scale) * geometry_response
        elif self.cartesian_encoder is not None and geometry_available is not None \
                and bool(geometry_available.any()):
            raise ValueError("available geometry requires coordinates and edges")

        support_count = support_atoms.shape[1]
        support_ligand, query_ligand = torch.split(
            ligand, (support_count, count - support_count), 1)
        support_endpoint, query_endpoint = torch.split(
            endpoint, (support_count, count - support_count), 1)
        support_primitive, query_primitive = torch.split(
            primitive, (support_count, count - support_count), 1)
        _, query_section = torch.split(
            section, (support_count, count - support_count), 1)
        state = self.meta.infer(
            pooled, support_ligand, support_endpoint, support_primitive,
            support_y.to(device, dtype), query_ligand, query_endpoint,
            query_primitive, adapt=adapt,
            task_state_override=task_state_override)
        if adapt and support_count:
            sar = (state.reliability.unsqueeze(-1) * state.coefficients
                   * query_primitive).sum(-1)
            sar = self.meta.mechanism_scale * sar
        else:
            sar = torch.zeros_like(state.zero_shot)
        level_baseline = state.zero_shot + state.level_adjustment
        prediction = level_baseline + sar
        evidence = (state.exact_gradient.square().mean((-2, -1)).sqrt()
                    if support_count else state.zero_shot.new_zeros(batch))
        one = torch.ones_like(state.level_gate)
        output = QPSMPMetaOutput(
            prediction, state.query_add, state.query_ligand_value,
            state.query_cross, level_baseline, state.level_adjustment, sar,
            prediction - state.zero_shot, state.zero_shot, state.task_evidence,
            state.level_shift, query_primitive, state.residual,
            state.exact_gradient, evidence, state.level_gate, one,
            state.reliability.mean(-1), state.primitive_regularizer)
        if not unbatched:
            return output
        return QPSMPMetaOutput(*(value.squeeze(0)
            if isinstance(value, Tensor) and value.ndim else value
            for value in output.__dict__.values()))
