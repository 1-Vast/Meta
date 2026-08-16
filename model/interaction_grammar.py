"""Shared interaction grammar with transferability-gated residual transport.

Motivation (Stage 0 audit of the retained BPSF baseline):

* the zero-shot endpoint collapsed to a near constant (within-episode spread
  0.065 pK), so k=0 carried no protein- or ligand-specific information;
* swapping the protein for a cross-component donor moved the zero-shot output
  by 0.009 pK, i.e. the trunk was protein-blind;
* the pooled-protein, ligand-baseline, section and task-state branches received
  no gradient from the query loss at any k;
* at k=1 the residual kernel reduced to `mean(r)`, so no query-specific
  adaptation existed by construction.

This module replaces the deep bipartite pair trunk with an atom-resolved cross
attention onto residue slots, a globally shared contact-type dictionary, and a
few-shot transport whose per-support coefficient depends on the query.  Every
trainable tensor is on the path from inputs to `prediction`.

Contracts preserved from the active model:

* `adapt=False` or `k=0` returns exactly the zero-shot endpoint;
* support labels enter only as residual values (label locked);
* support ordering does not change the output;
* queries are processed independently (permutation equivariance);
* query labels never enter the model.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor

from .encoders import LigandEncoder
from .qpsmp_meta import QPSMPMetaOutput


class ResidueEncoder(nn.Module):
    """Project cached ESM residue slots and a pooled protein summary."""

    def __init__(self, protein_dim: int, hidden_dim: int,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.residue = nn.Linear(protein_dim, hidden_dim, dtype=dtype)
        self.residue_mlp = nn.Sequential(
            nn.LayerNorm(hidden_dim, dtype=dtype),
            nn.Linear(hidden_dim, 2 * hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(2 * hidden_dim, hidden_dim, dtype=dtype))
        self.residue_norm = nn.LayerNorm(hidden_dim, dtype=dtype)
        self.chemistry = nn.Linear(4, hidden_dim, bias=False, dtype=dtype)
        self.pooled = nn.Sequential(
            nn.Linear(protein_dim, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim, dtype=dtype))
        self.pooled_norm = nn.LayerNorm(hidden_dim, dtype=dtype)

    def forward(self, pooled: Tensor, tokens: Tensor, mask: Tensor,
                chemistry: Tensor | None) -> tuple[Tensor, Tensor]:
        gate = mask.unsqueeze(-1)
        residues = self.residue(tokens)
        if chemistry is not None:
            residues = residues + self.chemistry(chemistry)
        residues = self.residue_norm(residues + self.residue_mlp(residues)) * gate
        summary = self.pooled_norm(
            self.pooled(pooled) + (residues.sum(1) / gate.sum(1).clamp_min(1.0)))
        return residues, summary


class ContactGrammar(nn.Module):
    """Atom-to-residue cross attention with a globally shared type dictionary.

    The dictionary index is a source-shared interaction primitive: the same
    contact types are used for every target, so few-shot evidence about one
    ligand is expressed in coordinates that transfer to unseen proteins.
    """

    def __init__(self, hidden_dim: int, contact_types: int, heads: int,
                 embed_dim: int, dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if hidden_dim % heads:
            raise ValueError("hidden_dim must be divisible by contact heads")
        self.heads = int(heads)
        self.head_dim = hidden_dim // heads
        self.contact_types = int(contact_types)
        self.atom_query = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.residue_key = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.residue_value = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.atom_chemistry = nn.Linear(4, hidden_dim, bias=False, dtype=dtype)
        self.atom_context = nn.Sequential(
            nn.LayerNorm(3 * hidden_dim, dtype=dtype),
            nn.Linear(3 * hidden_dim, 2 * hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(2 * hidden_dim, hidden_dim, dtype=dtype))
        self.type_logit = nn.Linear(hidden_dim, contact_types, dtype=dtype)
        self.type_strength = nn.Linear(hidden_dim, contact_types, dtype=dtype)

    def forward(self, atoms: Tensor, atom_mask: Tensor, atom_chemistry: Tensor,
                residues: Tensor, residue_mask: Tensor
                ) -> tuple[Tensor, Tensor, Tensor]:
        pairs, atom_count, hidden = atoms.shape
        residue_count = residues.shape[1]
        atoms = atoms + self.atom_chemistry(atom_chemistry)
        query = self.atom_query(atoms).reshape(
            pairs, atom_count, self.heads, self.head_dim).transpose(1, 2)
        key = self.residue_key(residues).reshape(
            pairs, residue_count, self.heads, self.head_dim).transpose(1, 2)
        value = self.residue_value(residues).reshape(
            pairs, residue_count, self.heads, self.head_dim).transpose(1, 2)
        score = query @ key.transpose(-1, -2) / self.head_dim ** 0.5
        score = score.masked_fill(
            ~residue_mask[:, None, None, :].bool(), torch.finfo(score.dtype).min)
        weight = torch.softmax(score, dim=-1)
        context = (weight @ value).transpose(1, 2).reshape(
            pairs, atom_count, hidden)
        state = self.atom_context(
            torch.cat((atoms, context, atoms * context), dim=-1))
        state = state * atom_mask.unsqueeze(-1)
        # Occupancy of each shared contact type, saturating in the atom count so
        # that large ligands cannot win by size alone.
        occupancy = (torch.sigmoid(self.type_logit(state))
                     * F.softplus(self.type_strength(state)))
        occupancy = (occupancy * atom_mask.unsqueeze(-1)).sum(1)
        occupancy = occupancy / (1.0 + occupancy)
        denominator = atom_mask.sum(1, keepdim=True).clamp_min(1.0)
        mean_state = (state * atom_mask.unsqueeze(-1)).sum(1) / denominator
        max_state = state.masked_fill(
            atom_mask.unsqueeze(-1) == 0, torch.finfo(state.dtype).min).amax(1)
        max_state = torch.nan_to_num(max_state, neginf=0.0)
        return occupancy, mean_state, max_state


class TransferabilityTransport(nn.Module):
    """Label-locked residual transport with a query-specific support gate.

    ``f = f0(q) + shrink(n) * sum_k softmax_k(sim) * rho(q,k) * r_k``

    ``rho`` is a bounded per-(query, support) transferability coefficient, so a
    single support observation already produces a query-specific correction.
    Setting ``rho == 1`` and a flat ``sim`` recovers the shrunken support mean,
    which keeps level-only abstention inside the hypothesis class.
    """

    def __init__(self, embed_dim: int, hidden_dim: int = 64,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.key = nn.Linear(embed_dim, hidden_dim, bias=False, dtype=dtype)
        self.gate = nn.Sequential(
            nn.Linear(4 * embed_dim, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, 1, dtype=dtype))
        self.log_temperature = nn.Parameter(torch.tensor(1.443, dtype=dtype))
        self.log_shrinkage = nn.Parameter(torch.tensor(1.8546, dtype=dtype))
        # Start close to `rho == 1` (the shrunken support mean) but not exactly
        # on it: a zero readout would leave the gate hidden layer without
        # gradient and would make the k=1 correction a pure level shift.
        nn.init.normal_(self.gate[-1].weight, std=5e-2)
        nn.init.zeros_(self.gate[-1].bias)

    def shrinkage(self, support_count: int, reference: Tensor) -> Tensor:
        strength = F.softplus(self.log_shrinkage)
        count = reference.new_tensor(float(support_count))
        return count / (count + strength)

    def forward(self, support_embed: Tensor, query_embed: Tensor,
                support_residual: Tensor) -> tuple[Tensor, Tensor, Tensor]:
        support_count = support_embed.shape[1]
        query_count = query_embed.shape[1]
        support_key = F.normalize(self.key(support_embed), dim=-1)
        query_key = F.normalize(self.key(query_embed), dim=-1)
        temperature = F.softplus(self.log_temperature) + 1.0
        weight = torch.softmax(
            temperature * torch.einsum("bqh,bkh->bqk", query_key, support_key), -1)
        left = query_embed[:, :, None, :].expand(-1, -1, support_count, -1)
        right = support_embed[:, None, :, :].expand(-1, query_count, -1, -1)
        gate = 2.0 * torch.sigmoid(self.gate(torch.cat(
            (left, right, left * right, (left - right).abs()), -1)).squeeze(-1))
        transport = torch.einsum(
            "bqk,bqk,bk->bq", weight, gate, support_residual)
        # Parameter-free evidence summary: a reported diagnostic must never
        # carry trainable weights that the query loss cannot reach.
        summary = (support_embed * support_residual.unsqueeze(-1)).mean(1)
        return transport, gate, summary


class InteractionGrammarModel(nn.Module):
    """Protein-conditioned zero-shot trunk plus query-specific few-shot transport."""

    def __init__(self, protein_dim: int, hidden_dim: int = 192,
                 task_dim: int = 48, ligand_layers: int = 4,
                 pair_dim: int = 96, pair_blocks: int = 4,
                 pair_latents: int = 24, pair_heads: int = 8,
                 pair_chunk_size: int = 8, support_hidden_dim: int = 192,
                 support_blocks: int = 2, adapter_rank: int = 4,
                 adaptive_blocks: int = 2, adapter_scale: float = 0.25,
                 locality_topk: int = 32, use_cartesian: bool = False,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        del pair_blocks, support_blocks, adapter_rank, adaptive_blocks
        del adapter_scale, locality_topk
        if use_cartesian:
            raise ValueError(
                "the interaction-grammar trunk has no common-frame coordinate path")
        embed_dim = int(pair_dim)
        contact_types = int(pair_latents)
        self.contact_types = contact_types
        self.embed_dim = embed_dim
        self.protein_encoder = ResidueEncoder(protein_dim, hidden_dim, dtype)
        self.ligand_encoder = LigandEncoder(hidden_dim, ligand_layers, dtype=dtype)
        self.grammar = ContactGrammar(
            hidden_dim, contact_types, pair_heads, embed_dim, dtype)
        self.embed = nn.Sequential(
            nn.Linear(4 * hidden_dim + contact_types, 2 * embed_dim, dtype=dtype),
            nn.GELU(),
            nn.Linear(2 * embed_dim, embed_dim, dtype=dtype))
        self.embed_norm = nn.LayerNorm(embed_dim, dtype=dtype)
        self.section = nn.Linear(embed_dim, task_dim, bias=False, dtype=dtype)
        self.section_norm = nn.LayerNorm(
            task_dim, elementwise_affine=False, dtype=dtype)
        self.ligand_head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, 1, dtype=dtype))
        self.protein_head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, 1, dtype=dtype))
        self.interaction_head = nn.Sequential(
            nn.Linear(embed_dim + task_dim, embed_dim, dtype=dtype), nn.GELU(),
            nn.Linear(embed_dim, 1, dtype=dtype))
        self.contact_weight = nn.Linear(contact_types, 1, bias=False, dtype=dtype)
        self.transport = TransferabilityTransport(embed_dim, task_dim, dtype)

    @staticmethod
    def atom_chemistry(raw_atoms: Tensor) -> Tensor:
        positive = raw_atoms[..., 21]
        negative = raw_atoms[..., 18] + raw_atoms[..., 19]
        aromatic = raw_atoms[..., 27]
        hydrophobic = raw_atoms[..., [1, 4, 7, 8, 9]].sum(-1).clamp_max(1)
        return torch.stack((positive, negative, aromatic, hydrophobic), -1)

    def dictionary_regularizer(self, occupancy: Tensor) -> Tensor:
        """Keep shared contact types distinct and used."""
        centered = occupancy - occupancy.mean(1, keepdim=True)
        count = max(occupancy.shape[1] - 1, 1)
        covariance = torch.einsum("bqm,bqn->bmn", centered, centered) / count
        scale = covariance.diagonal(dim1=-2, dim2=-1).clamp_min(1e-6).sqrt()
        correlation = covariance / (scale[:, :, None] * scale[:, None, :])
        identity = torch.eye(correlation.shape[-1], device=occupancy.device,
                             dtype=occupancy.dtype)
        usage = (occupancy.mean((0, 1)) - 0.5).square().mean()
        return usage + (correlation - identity).square().mean()

    def refine_slots(self, residues: Tensor, protein_mask: Tensor) -> Tensor:
        """Hook for protein-slot refinement. Identity in this model.

        Subclasses may override to refine the ordered protein slots before the
        interaction trunk consumes them. The default is an exact identity so
        that every existing architecture is numerically unchanged.
        """
        return residues

    def encode(self, protein_pooled: Tensor, protein_tokens: Tensor,
               protein_mask: Tensor, raw_atoms: Tensor, bonds: Tensor,
               mask: Tensor, protein_chemistry: Tensor | None):
        batch, count = raw_atoms.shape[:2]
        residues, summary = self.protein_encoder(
            protein_pooled, protein_tokens, protein_mask, protein_chemistry)
        residues = self.refine_slots(residues, protein_mask)
        ligand, atom_states = self.ligand_encoder(
            raw_atoms.flatten(0, 1), bonds.flatten(0, 1), mask.flatten(0, 1))
        residue_count = residues.shape[1]
        wide_residues = residues[:, None].expand(-1, count, -1, -1).reshape(
            batch * count, residue_count, -1)
        wide_mask = protein_mask[:, None].expand(-1, count, -1).reshape(
            batch * count, residue_count)
        occupancy, mean_state, max_state = self.grammar(
            atom_states, mask.flatten(0, 1),
            self.atom_chemistry(raw_atoms).flatten(0, 1),
            wide_residues, wide_mask)
        wide_summary = summary[:, None].expand(-1, count, -1).reshape(
            batch * count, -1)
        embed = self.embed_norm(self.embed(torch.cat(
            (ligand, mean_state, max_state, wide_summary, occupancy), -1)))
        section = self.section_norm(self.section(embed))
        ligand_value = self.ligand_head(ligand).squeeze(-1)
        protein_value = self.protein_head(wide_summary).squeeze(-1)
        interaction = (self.interaction_head(
            torch.cat((embed, section), -1)).squeeze(-1)
            + self.contact_weight(occupancy).squeeze(-1))
        endpoint = ligand_value + protein_value + interaction
        shape = (batch, count)
        return (endpoint.reshape(shape), ligand_value.reshape(shape),
                protein_value.reshape(shape), embed.reshape(*shape, -1),
                occupancy.reshape(*shape, -1))

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
        if task_state_override is not None:
            raise ValueError("the grammar trunk does not accept transplanted states")
        if geometry_coordinates is not None or geometry_edge_index is not None or (
                geometry_available is not None and bool(geometry_available.any())):
            raise ValueError("the grammar trunk has no legal coordinate path")
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
        raw_atoms = torch.cat((support_atoms, query_atoms), 1).to(device, dtype)
        bonds = torch.cat((support_bonds, query_bonds), 1).to(device, dtype)
        mask = torch.cat((support_mask, query_mask), 1).to(device, dtype)
        endpoint, ligand_value, protein_value, embed, occupancy = self.encode(
            protein_pooled.to(device, dtype), protein_tokens.to(device, dtype),
            protein_mask.to(device, dtype), raw_atoms, bonds, mask,
            None if protein_chemistry is None
            else protein_chemistry.to(device, dtype))
        support_count = support_atoms.shape[1]
        query_count = endpoint.shape[1] - support_count
        support_zero, zero_shot = torch.split(
            endpoint, (support_count, query_count), 1)
        support_embed, query_embed = torch.split(
            embed, (support_count, query_count), 1)
        _, query_occupancy = torch.split(
            occupancy, (support_count, query_count), 1)
        _, query_ligand_value = torch.split(
            ligand_value, (support_count, query_count), 1)
        _, query_protein_value = torch.split(
            protein_value, (support_count, query_count), 1)
        batch = endpoint.shape[0]
        if not adapt or support_count == 0:
            level_shift = zero_shot.new_zeros(batch)
            level_gate = zero_shot.new_zeros(batch)
            level_adjustment = zero_shot.new_zeros(batch, 1)
            residual = support_y.to(device, dtype).new_zeros(batch, support_count)
            transport = torch.zeros_like(zero_shot)
            gate = zero_shot.new_zeros(batch, query_count, max(support_count, 1))
            summary = query_embed.new_zeros(batch, self.embed_dim)
            evidence = zero_shot.new_zeros(batch)
        else:
            raw_residual = support_y.to(device, dtype) - support_zero
            locked = raw_residual.detach()
            level_shift = locked.mean(-1)
            shrink = self.transport.shrinkage(support_count, locked)
            level_gate = shrink.expand_as(level_shift)
            level_adjustment = shrink * level_shift.unsqueeze(-1)
            residual = locked - level_adjustment
            transport, gate, summary = self.transport(
                support_embed, query_embed, locked)
            transport = shrink * transport
            evidence = (locked.unsqueeze(-1) * occupancy[:, :support_count]
                        ).square().mean((-2, -1)).sqrt()
        level_baseline = zero_shot + level_adjustment
        sar = transport - level_adjustment
        prediction = zero_shot + transport
        task_state = summary[:, None].expand(-1, self.contact_types, -1)
        output = QPSMPMetaOutput(
            prediction=prediction,
            additive=query_ligand_value + query_protein_value,
            ligand_only=query_ligand_value,
            cross_zero_shot=zero_shot,
            level_baseline=level_baseline,
            level_adjustment=level_adjustment,
            sar_adaptation=sar,
            adaptation=prediction - zero_shot,
            zero_shot=zero_shot,
            task_state=task_state,
            level_shift=level_shift,
            query_basis=query_occupancy,
            support_residual_quotient=residual,
            support_evidence=locked.unsqueeze(-1) * occupancy[:, :support_count]
            if (adapt and support_count) else residual.unsqueeze(-1),
            evidence_score=evidence,
            level_shrinkage=level_gate,
            shape_scale=torch.ones_like(level_gate),
            sar_scale=gate.mean(-1).mean(-1) if gate.numel() else level_gate,
            support_match_loss=self.dictionary_regularizer(query_occupancy))
        if not unbatched:
            return output
        return QPSMPMetaOutput(*(value.squeeze(0)
            if isinstance(value, Tensor) and value.ndim else value
            for value in output.__dict__.values()))
