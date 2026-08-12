"""Learned quotient-preserving meta-potential for cold-target DTA.

The analytic centered section in :mod:`model.qpsmp` is a diagnostic helper.
This module is the trainable meta-learning path: source query losses train the
protein localizer, scalar potential, section basis, and support-set adapter.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import torch
import torch.nn as nn
from torch import Tensor

from .encoders import LigandEncoder, ProteinEncoder


@dataclass(frozen=True)
class QPSMPMetaOutput:
    prediction: Tensor
    additive: Tensor
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


class LigandConditionedProteinLocalizer(nn.Module):
    """Pool residue tokens separately for each ligand representation."""

    def __init__(self, hidden_dim: int, dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if hidden_dim < 1:
            raise ValueError("hidden_dim must be positive")
        self.key = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.query = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.value = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.norm = nn.LayerNorm(hidden_dim, elementwise_affine=False, dtype=dtype)

    def forward(
            self, protein_tokens: Tensor, protein_mask: Tensor,
            ligand_states: Tensor) -> Tensor:
        if protein_tokens.ndim != 2 or ligand_states.ndim != 2:
            raise ValueError("tokens and ligand states must be rank-two")
        if protein_mask.ndim != 1 or protein_mask.shape[0] != protein_tokens.shape[0]:
            raise ValueError("protein mask does not match protein tokens")
        if protein_tokens.shape[1] != ligand_states.shape[1]:
            raise ValueError("protein and ligand hidden dimensions differ")
        if not bool(protein_mask.any().item()):
            raise ValueError("protein mask contains no valid residue token")
        logits = self.query(ligand_states) @ self.key(protein_tokens).T
        logits = logits / math.sqrt(protein_tokens.shape[1])
        logits = logits.masked_fill(~protein_mask.bool().unsqueeze(0), -torch.inf)
        weights = torch.softmax(logits, dim=-1)
        localized = weights @ self.value(protein_tokens)
        return self.norm(localized)


class CenteredNeuralAdapter(nn.Module):
    """Permutation-invariant learned task update with an exact quotient null."""

    def __init__(
            self, interaction_dim: int, task_dim: int,
            hidden_dim: int, state_bound: float = 1.0,
            dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if min(interaction_dim, task_dim, hidden_dim) < 1 or state_bound <= 0:
            raise ValueError("adapter dimensions and state_bound must be positive")
        self.state_bound = float(state_bound)
        self.value = nn.Linear(interaction_dim, hidden_dim, bias=False, dtype=dtype)
        # Bias-free maps make c(0)=0 an architectural identity.
        self.update = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype),
            nn.Tanh(),
            nn.Linear(hidden_dim, task_dim, bias=False, dtype=dtype),
        )

    def forward(self, support_interaction: Tensor, centered_residual: Tensor) -> Tensor:
        if support_interaction.ndim != 2 or centered_residual.ndim != 1:
            raise ValueError("adapter inputs have invalid ranks")
        if support_interaction.shape[0] != centered_residual.shape[0]:
            raise ValueError("adapter inputs have different support sizes")
        evidence = (
            centered_residual.unsqueeze(-1) * self.value(support_interaction)
        ).mean(dim=0)
        raw = self.update(evidence)
        return self.state_bound * raw / (1.0 + torch.linalg.vector_norm(raw))


class QPSMPMetaLearner(nn.Module):
    """Shared scalar potential plus a learned centered episodic update."""

    def __init__(
            self, hidden_dim: int, task_dim: int,
            adapter_hidden_dim: int | None = None,
            state_bound: float = 1.0,
            dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if hidden_dim < 1 or task_dim < 1:
            raise ValueError("hidden_dim and task_dim must be positive")
        adapter_hidden_dim = hidden_dim if adapter_hidden_dim is None else adapter_hidden_dim
        self.localizer = LigandConditionedProteinLocalizer(hidden_dim, dtype)
        self.interaction = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype),
            nn.SiLU(),
            nn.LayerNorm(hidden_dim, elementwise_affine=False, dtype=dtype),
        )
        self.ligand_baseline = nn.Linear(hidden_dim, 1, dtype=dtype)
        self.protein_level = nn.Linear(hidden_dim, 1, dtype=dtype)
        self.zero_shot_head = nn.Linear(hidden_dim, 1, bias=False, dtype=dtype)
        self.section_head = nn.Linear(hidden_dim, task_dim, bias=False, dtype=dtype)
        # Start near the strong support-level baseline while retaining a
        # nonzero query-shape path and gradient.
        # Positive random-intercept and observation variances induce a
        # support-size-aware empirical-Bayes level shrinkage.
        self.level_prior_raw = nn.Parameter(torch.tensor(0.5413249, dtype=dtype))
        self.level_noise_raw = nn.Parameter(torch.tensor(-0.2981850, dtype=dtype))
        self.shape_logit = nn.Parameter(torch.tensor(-2.1972246, dtype=dtype))
        self.sar_logit = nn.Parameter(torch.tensor(-2.1972246, dtype=dtype))
        self.adapter = CenteredNeuralAdapter(
            hidden_dim, task_dim, adapter_hidden_dim, state_bound, dtype)

    def interaction_features(
            self, protein_pooled: Tensor, protein_tokens: Tensor,
            protein_mask: Tensor, ligand_states: Tensor) -> Tensor:
        parameter = next(self.parameters())
        target_device, target_dtype = parameter.device, parameter.dtype
        protein_pooled = protein_pooled.to(device=target_device, dtype=target_dtype)
        protein_tokens = protein_tokens.to(device=target_device, dtype=target_dtype)
        protein_mask = protein_mask.to(device=target_device)
        ligand_states = ligand_states.to(device=target_device, dtype=target_dtype)
        if protein_pooled.ndim != 1 or ligand_states.ndim != 2:
            raise ValueError("pooled protein and ligand states have invalid ranks")
        if protein_pooled.shape[0] != ligand_states.shape[1]:
            raise ValueError("protein and ligand hidden dimensions differ")
        localized = self.localizer(protein_tokens, protein_mask, ligand_states)
        # Only a crossed term enters the interaction heads. Protein-only and
        # ligand-only information have separate additive baseline channels.
        return self.interaction(localized * ligand_states)

    def scalar_components(
            self, protein_pooled: Tensor, ligand_states: Tensor,
            interaction: Tensor) -> tuple[Tensor, Tensor, Tensor]:
        ligand = self.ligand_baseline(ligand_states).squeeze(-1)
        protein = self.protein_level(protein_pooled).squeeze(-1)
        zero_shot = self.zero_shot_head(interaction).squeeze(-1)
        return ligand + protein, zero_shot, self.section_head(interaction)

    def forward(
            self, protein_pooled: Tensor, protein_tokens: Tensor,
            protein_mask: Tensor, support_ligand: Tensor, support_y: Tensor,
            query_ligand: Tensor, *, adapt: bool = True,
            ) -> QPSMPMetaOutput:
        parameter = next(self.parameters())
        device, dtype = parameter.device, parameter.dtype
        protein_pooled = protein_pooled.to(device=device, dtype=dtype)
        protein_tokens = protein_tokens.to(device=device, dtype=dtype)
        protein_mask = protein_mask.to(device=device)
        support_y = support_y.to(device=device, dtype=dtype)
        support_ligand = support_ligand.to(device=device, dtype=dtype)
        query_ligand = query_ligand.to(device=device, dtype=dtype)
        if support_y.ndim != 1 or support_y.shape[0] != support_ligand.shape[0]:
            raise ValueError("support_y must have one value per support ligand")
        all_ligand = torch.cat((support_ligand, query_ligand), dim=0)
        all_interaction = self.interaction_features(
            protein_pooled, protein_tokens, protein_mask, all_ligand)
        support_count = support_ligand.shape[0]
        support_interaction = all_interaction[:support_count]
        query_interaction = all_interaction[support_count:]
        support_add, support_zero, _ = self.scalar_components(
            protein_pooled, support_ligand, support_interaction)
        query_add, query_zero, query_basis = self.scalar_components(
            protein_pooled, query_ligand, query_interaction)
        support_scalar = support_add + support_zero
        additive = query_add
        cross_zero_shot = query_zero
        zero_shot = additive + cross_zero_shot

        if not adapt or support_count == 0:
            state = torch.zeros(
                query_basis.shape[1], device=query_basis.device, dtype=query_basis.dtype)
            level_shift = torch.zeros((), device=query_basis.device, dtype=query_basis.dtype)
            centered_residual = support_y.new_zeros(support_y.shape)
            prediction = zero_shot
            level_baseline = zero_shot
            level_adjustment = torch.zeros_like(prediction)
            sar_adaptation = torch.zeros_like(prediction)
            adaptation = torch.zeros_like(prediction)
            evidence_score = torch.zeros((), device=query_basis.device, dtype=query_basis.dtype)
            level_shrinkage = torch.zeros((), device=query_basis.device, dtype=query_basis.dtype)
            shape_scale = torch.sigmoid(self.shape_logit)
            sar_scale = torch.sigmoid(self.sar_logit)
        else:
            residual = support_y - support_scalar
            centered_residual = residual - residual.mean()
            level_shift = residual.mean()
            state = self.adapter(support_interaction, centered_residual)
            evidence_score = torch.tanh(torch.linalg.vector_norm(state))
            prior_variance = torch.nn.functional.softplus(self.level_prior_raw)
            noise_variance = torch.nn.functional.softplus(self.level_noise_raw)
            level_shrinkage = (
                support_count * prior_variance
                / (noise_variance + support_count * prior_variance)
            )
            shape_scale = torch.sigmoid(self.shape_logit)
            sar_scale = torch.sigmoid(self.sar_logit)
            support_center = support_scalar.mean()
            calibrated_level = support_center + level_shrinkage * level_shift
            level_baseline = calibrated_level + shape_scale * (zero_shot - support_center)
            level_adjustment = level_baseline - zero_shot
            sar_adaptation = sar_scale * (query_basis @ state)
            prediction = level_baseline + sar_adaptation
            adaptation = prediction - zero_shot
        evidence = torch.linalg.vector_norm(centered_residual) / max(support_count, 1) ** 0.5
        return QPSMPMetaOutput(
            prediction=prediction,
            additive=additive,
            cross_zero_shot=cross_zero_shot,
            level_baseline=level_baseline,
            level_adjustment=level_adjustment,
            sar_adaptation=sar_adaptation,
            adaptation=adaptation,
            zero_shot=zero_shot,
            task_state=state,
            level_shift=level_shift,
            query_basis=query_basis,
            support_residual_quotient=centered_residual,
            support_evidence=evidence,
            evidence_score=evidence_score,
            level_shrinkage=level_shrinkage,
            shape_scale=shape_scale,
            sar_scale=sar_scale,
        )

    @staticmethod
    def delta(predictions: Tensor, left: Tensor, right: Tensor) -> Tensor:
        return predictions.index_select(0, right) - predictions.index_select(0, left)

    @staticmethod
    def rectangle(delta_a: Tensor, delta_b: Tensor) -> Tensor:
        return delta_a - delta_b


class QPSMPBioModel(nn.Module):
    """End-to-end wrapper from cached biological tensors to QPSMP outputs."""

    def __init__(
            self, protein_dim: int, hidden_dim: int, task_dim: int,
            ligand_layers: int = 2, adapter_hidden_dim: int | None = None,
            state_bound: float = 1.0,
            dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.protein_encoder = ProteinEncoder(protein_dim, hidden_dim, dtype)
        self.ligand_encoder = LigandEncoder(hidden_dim, ligand_layers, dtype=dtype)
        self.meta = QPSMPMetaLearner(
            hidden_dim, task_dim, adapter_hidden_dim, state_bound, dtype)

    def forward(
            self, protein_pooled: Tensor, protein_tokens: Tensor,
            protein_mask: Tensor, support_atoms: Tensor, support_bonds: Tensor,
            support_mask: Tensor, support_y: Tensor, query_atoms: Tensor,
            query_bonds: Tensor, query_mask: Tensor, *, adapt: bool = True) -> QPSMPMetaOutput:
        parameter = next(self.parameters())
        device, dtype = parameter.device, parameter.dtype
        pooled, tokens = self.protein_encoder(
            protein_pooled.to(device=device, dtype=dtype).unsqueeze(0),
            protein_tokens.to(device=device, dtype=dtype).unsqueeze(0))
        all_atoms = torch.cat((support_atoms, query_atoms), dim=0).to(device=device, dtype=dtype)
        all_bonds = torch.cat((support_bonds, query_bonds), dim=0).to(device=device, dtype=dtype)
        all_mask = torch.cat((support_mask, query_mask), dim=0).to(device=device, dtype=dtype)
        ligand, _ = self.ligand_encoder(all_atoms, all_bonds, all_mask)
        support_count = support_atoms.shape[0]
        return self.meta(
            pooled.squeeze(0), tokens.squeeze(0), protein_mask.to(device=device),
            ligand[:support_count], support_y,
            ligand[support_count:], adapt=adapt)
