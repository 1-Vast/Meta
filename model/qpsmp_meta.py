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


class EvidenceLockedMetaTransport(nn.Module):
    """Label-blind keys transport label-bound residual values to queries."""

    def __init__(self, ligand_dim: int, primitive_count: int,
                 hidden_dim: int = 32,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.primitive_count = int(primitive_count)
        self.hidden_dim = int(hidden_dim)
        self.primitive_identity = nn.Parameter(
            torch.empty(primitive_count, hidden_dim, dtype=dtype))
        nn.init.normal_(self.primitive_identity, std=hidden_dim ** -0.5)
        self.protein_key = nn.Linear(ligand_dim, hidden_dim, dtype=dtype)
        self.ligand_change = nn.Sequential(
            nn.LayerNorm(4 * ligand_dim, dtype=dtype),
            nn.Linear(4 * ligand_dim, hidden_dim, dtype=dtype), nn.GELU())
        width = 3 * hidden_dim + 2
        self.score = nn.Sequential(
            nn.LayerNorm(width, dtype=dtype),
            nn.Linear(width, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, 1, bias=False, dtype=dtype))
        self.direction = nn.Sequential(
            nn.LayerNorm(width, dtype=dtype),
            nn.Linear(width, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, 1, bias=False, dtype=dtype))
        self.reliability = nn.Parameter(
            torch.tensor([0.0, 1.0, 1.0, 0.5], dtype=dtype))

    @staticmethod
    def label_value(residual: Tensor, primitive: Tensor) -> Tensor:
        if residual.shape != primitive.shape[:-1]:
            raise ValueError("residual and primitive responses disagree")
        return residual.unsqueeze(-1) * primitive

    def forward(self, protein: Tensor, support_ligand: Tensor,
                support_primitive: Tensor, support_residual: Tensor,
                query_ligand: Tensor, query_primitive: Tensor, *,
                adapt: bool = True,
                task_state_override: Tensor | None = None
                ) -> tuple[Tensor, Tensor, Tensor, Tensor]:
        batch, support_count, primitive_count = support_primitive.shape
        query_count = query_primitive.shape[1]
        if primitive_count != self.primitive_count or \
                query_primitive.shape != (batch, query_count, primitive_count):
            raise ValueError("support/query primitive dictionaries disagree")
        if task_state_override is not None:
            raise ValueError("ELMT does not permit transplanted task states")
        zero_coeff = query_primitive.new_zeros(batch, query_count, primitive_count)
        zero_state = query_primitive.new_zeros(batch, primitive_count, self.hidden_dim)
        zero_reliability = query_primitive.new_zeros(batch, query_count)
        value = self.label_value(support_residual, support_primitive)
        if not adapt or support_count == 0:
            return zero_state, zero_coeff, zero_reliability, value

        support = support_ligand[:, None].expand(-1, query_count, -1, -1)
        query = query_ligand[:, :, None].expand(-1, -1, support_count, -1)
        change = self.ligand_change(torch.cat((
            support, query, query - support, query * support), -1))
        change = change[:, :, :, None].expand(-1, -1, -1, primitive_count, -1)
        identity = self.primitive_identity[None, None, None].expand(
            batch, query_count, support_count, -1, -1)
        protein_key = self.protein_key(protein)[:, None, None, None].expand(
            -1, query_count, support_count, primitive_count, -1)
        scalar = torch.stack((
            support_primitive[:, None].expand(-1, query_count, -1, -1),
            query_primitive[:, :, None].expand(-1, -1, support_count, -1)), -1)
        key = torch.cat((change, identity, protein_key, scalar), -1)
        weights = torch.softmax(self.score(key).squeeze(-1), dim=2)
        direction = torch.tanh(self.direction(key).squeeze(-1))
        coefficients = (weights * direction * value[:, None]).sum(2)

        coherence = value.sum(1).abs() / value.abs().sum(1).clamp_min(1e-8)
        coherence = coherence.mean(-1, keepdim=True).expand(-1, query_count)
        support_norm = F.normalize(support_ligand, dim=-1)
        query_norm = F.normalize(query_ligand, dim=-1)
        coverage = 1.0 - torch.einsum(
            "bqh,bkh->bqk", query_norm, support_norm).amax(-1)
        b0, b1, b2, b3 = self.reliability
        reliability = torch.sigmoid(
            b0 + F.softplus(b1) * coherence - F.softplus(b2) * coverage
            + F.softplus(b3) * query_primitive.new_tensor(float(support_count)).log1p())
        task_state = value.mean(1).unsqueeze(-1) * self.primitive_identity[None]
        return task_state, coefficients, reliability, value


ELMT = EvidenceLockedMetaTransport
# Compatibility names remain importable, but the active implementation is ELMT.
TERM = EvidenceLockedMetaTransport
TriadicEvidenceRouter = EvidenceLockedMetaTransport


class QPSMPMetaLearner(nn.Module):
    """Scalar baseline, conserved level residual, and ELMT composition."""

    def __init__(self, hidden_dim: int, primitive_count: int,
                 task_dim: int = 32, mechanism_scale: float = 1.0,
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
        self.term = EvidenceLockedMetaTransport(
            hidden_dim, primitive_count, task_dim, dtype)

    @staticmethod
    def primitive_regularizer(query_primitive: Tensor) -> Tensor:
        centered = query_primitive - query_primitive.mean(1, keepdim=True)
        mean_loss = (query_primitive.mean(1) - 0.5).square().mean()
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
            raw_residual = support_y - support_zero
            centered = raw_residual - raw_residual.mean(-1, keepdim=True)
            level_shift = raw_residual.mean(-1)
            count = torch.full_like(level_shift, float(support_count)).log1p()
            level_gate = self.level_gate(torch.stack((
                level_shift, centered.abs().mean(-1), count), -1)).squeeze(-1)
            # A scalar shift cannot alter within-target ligand ordering. A
            # learned shrinkage leaves non-zero k=1 evidence for ELMT.
            level_adjustment = level_gate.unsqueeze(-1) * level_shift.unsqueeze(-1)
            # Conservation: the same residual cannot be fully spent by both
            # level calibration and structural transport.
            residual = raw_residual - level_adjustment
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
    """Localized sequence/2D CIPF with evidence-locked meta transport."""

    def __init__(self, protein_dim: int, hidden_dim: int = 128,
                 task_dim: int = 32, ligand_layers: int = 3,
                 pair_dim: int = 64, pair_blocks: int = 3,
                 pair_latents: int = 16, pair_heads: int = 8,
                 pair_chunk_size: int = 16, support_hidden_dim: int = 128,
                 support_blocks: int = 2, adapter_rank: int = 4,
                 adaptive_blocks: int = 2, adapter_scale: float = 0.25,
                 locality_topk: int = 32,
                 use_cartesian: bool = False,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        del support_hidden_dim, support_blocks, adapter_rank, adaptive_blocks
        self.protein_encoder = ProteinEncoder(protein_dim, hidden_dim, dtype)
        self.ligand_encoder = LigandEncoder(hidden_dim, ligand_layers, dtype=dtype)
        self.locality_topk = int(locality_topk)
        self.residue_query = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.residue_key = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
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
            hidden_dim, pair_latents, task_dim, 1.0, dtype)

    def _pair_encode(self, atom_states: Tensor, residues: Tensor, mask: Tensor,
                     protein_mask: Tensor, bonds: Tensor,
                     atom_features: Tensor, atom_chemistry: Tensor,
                     residue_chemistry: Tensor):
        adjacency = (bonds.abs().sum(-1) > 0).to(atom_states.dtype)
        return self.pair_section(
            atom_states, residues, mask, protein_mask, adjacency,
            atom_features=atom_features, atom_chemistry=atom_chemistry,
            residue_chemistry=residue_chemistry)

    @staticmethod
    def _atom_chemistry(raw_atoms: Tensor) -> Tensor:
        positive = raw_atoms[..., 21]
        negative = raw_atoms[..., 18] + raw_atoms[..., 19]
        aromatic = raw_atoms[..., 27]
        hydrophobic = raw_atoms[..., [1, 4, 7, 8, 9]].sum(-1).clamp_max(1)
        return torch.stack((positive, negative, aromatic, hydrophobic), -1)

    def _localize(self, ligand: Tensor, residues: Tensor, protein_mask: Tensor,
                  residue_chemistry: Tensor) -> tuple[Tensor, Tensor, Tensor, Tensor]:
        scores = torch.einsum(
            "bch,brh->bcr", self.residue_query(ligand),
            self.residue_key(residues)) / residues.shape[-1] ** 0.5
        scores = scores.masked_fill(
            ~protein_mask[:, None].bool(), torch.finfo(scores.dtype).min)
        keep = min(self.locality_topk, residues.shape[1])
        selected_score, selected = scores.topk(keep, dim=-1)
        weights = torch.softmax(selected_score, -1) * keep
        batch, count = ligand.shape[:2]
        expanded = residues[:, None].expand(-1, count, -1, -1)
        local = torch.gather(expanded, 2, selected[..., None].expand(
            -1, -1, -1, residues.shape[-1]))
        local = local * weights[..., None]
        chemistry = torch.gather(
            residue_chemistry[:, None].expand(-1, count, -1, -1), 2,
            selected[..., None].expand(-1, -1, -1, residue_chemistry.shape[-1]))
        mask = torch.gather(
            protein_mask[:, None].expand(-1, count, -1), 2, selected)
        probability = torch.softmax(scores, -1)
        positions = torch.linspace(0, 1, residues.shape[1], device=residues.device,
                                   dtype=residues.dtype)
        center = (probability * positions).sum(-1, keepdim=True)
        locality = (probability * (positions - center).abs()).sum(-1).mean()
        return local, mask, chemistry, locality

    def forward(self, protein_pooled: Tensor, protein_tokens: Tensor,
                protein_mask: Tensor, support_atoms: Tensor, support_bonds: Tensor,
                support_mask: Tensor, support_y: Tensor, query_atoms: Tensor,
                query_bonds: Tensor, query_mask: Tensor, *, adapt: bool = True,
                protein_chemistry: Tensor | None = None,
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
            if protein_chemistry is not None:
                protein_chemistry = protein_chemistry.unsqueeze(0)
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
        if protein_chemistry is None:
            protein_chemistry = residues.new_zeros(
                residues.shape[0], residues.shape[1], 4)
        else:
            protein_chemistry = protein_chemistry.to(device, dtype)
        raw_atoms = torch.cat((support_atoms, query_atoms), 1).to(device, dtype)
        bonds = torch.cat((support_bonds, query_bonds), 1).to(device, dtype)
        mask = torch.cat((support_mask, query_mask), 1).to(device, dtype)
        batch, count, atom_count = raw_atoms.shape[:3]
        ligand, atom_states = self.ligand_encoder(
            raw_atoms.flatten(0, 1), bonds.flatten(0, 1), mask.flatten(0, 1))
        ligand = ligand.reshape(batch, count, -1)
        atom_states = atom_states.reshape(batch, count, atom_count, -1)
        local_residues, local_mask, local_chemistry, locality_loss = self._localize(
            ligand, residues, protein_mask.to(device), protein_chemistry)
        atom_chemistry = self._atom_chemistry(raw_atoms)
        encoded = self._pair_encode(
            atom_states.flatten(0, 1),
            local_residues.flatten(0, 1),
            mask.flatten(0, 1),
            local_mask.flatten(0, 1), bonds.flatten(0, 1),
            raw_atoms.flatten(0, 1), atom_chemistry.flatten(0, 1),
            local_chemistry.flatten(0, 1))
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
            state.reliability.mean(-1),
            state.primitive_regularizer + locality_loss)
        if not unbatched:
            return output
        return QPSMPMetaOutput(*(value.squeeze(0)
            if isinstance(value, Tensor) and value.ndim else value
            for value in output.__dict__.values()))
