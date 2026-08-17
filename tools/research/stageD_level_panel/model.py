"""Stage E candidate: panel-context level head + orthogonal level/shape routing.

Framework innovation (I1): the zero-shot protein level is replaced by a
panel-set level readout. For the query panel (the model's legal inputs at every
k), order-invariant mean/max pooling of the ligand encoder outputs plus the
protein summary is mapped by a small MLP to a per-episode scalar level. The
hypothesis, backed by the D0 measurements, is that a BindingDB target's mean
affinity is partly a property of which ligands were tested against it (assay
history / library composition), so the tested-panel composition carries level
information that the protein sequence alone does not (D0_LEVEL_IDENTIFIABILITY:
panel features level MSE 1.887 vs 2.155 constant, shuffle control 5.075;
D0_LEVEL_ANATOMY: panel composition held-out R^2 +0.239 vs protein +0.119).

The interaction/contact path and the ligand baseline stay exactly as in the
incumbent; only the additive protein_value(P) branch is replaced, so the level
head is the single architectural difference.

The level head consumes query ligand ENCODINGS and the protein summary only.
It never consumes support labels, query labels, or the transport correction.
"""
from __future__ import annotations

import torch
import torch.nn as nn
from torch import Tensor

from model.interaction_grammar import InteractionGrammarModel
from model.similarity_grammar import SimilarityTransport
from model.qpsmp_meta import QPSMPMetaOutput


class PanelLevelShapeModel(InteractionGrammarModel):
    """Incumbent trunk + Tanimoto transport + panel-set level head."""

    def __init__(self, *args, level_hidden: int = 64,
                 use_learned_key: bool = False, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        del self.protein_head          # replaced by the panel level head
        del self.transport             # replaced with the fixed Tanimoto transport
        self.transport = SimilarityTransport(
            self.embed_dim, kwargs.get("task_dim", 64),
            similarity_scale=8.0, use_learned_key=use_learned_key,
            dtype=next(self.parameters()).dtype)
        hidden_dim = int(kwargs.get("hidden_dim", 192))
        # level input: protein summary [H] + panel mean [H] + panel max [H]
        self.panel_level = nn.Sequential(
            nn.Linear(3 * hidden_dim, level_hidden,
                      dtype=next(self.parameters()).dtype),
            nn.GELU(),
            nn.Linear(level_hidden, 1, dtype=next(self.parameters()).dtype))

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
        interaction = (self.interaction_head(
            torch.cat((embed, section), -1)).squeeze(-1)
            + self.contact_weight(occupancy).squeeze(-1))
        endpoint = ligand_value + interaction
        shape = (batch, count)
        return (endpoint.reshape(shape), ligand_value.reshape(shape),
                ligand.reshape(*shape, -1), summary, embed.reshape(*shape, -1),
                occupancy.reshape(*shape, -1))

    def forward(self, protein_pooled: Tensor, protein_tokens: Tensor,
                protein_mask: Tensor, support_atoms: Tensor, support_bonds: Tensor,
                support_mask: Tensor, support_y: Tensor, query_atoms: Tensor,
                query_bonds: Tensor, query_mask: Tensor, *, adapt: bool = True,
                protein_chemistry: Tensor | None = None,
                support_fingerprint: Tensor | None = None,
                query_fingerprint: Tensor | None = None,
                task_state_override: Tensor | None = None,
                geometry_coordinates: Tensor | None = None,
                geometry_edge_index: Tensor | None = None,
                geometry_available: Tensor | None = None,
                geometry_common_frame: Tensor | None = None) -> QPSMPMetaOutput:
        if task_state_override is not None:
            raise ValueError("the panel-level trunk does not accept transplanted states")
        if geometry_coordinates is not None or geometry_edge_index is not None or (
                geometry_available is not None and bool(geometry_available.any())):
            raise ValueError("no BindingDB deployment pair has a common-frame complex")
        parameter = next(self.parameters())
        device, dtype = parameter.device, parameter.dtype
        unbatched = protein_pooled.ndim == 1
        if unbatched:
            protein_pooled = protein_pooled.unsqueeze(0)
            protein_tokens = protein_tokens.unsqueeze(0)
            protein_mask = protein_mask.unsqueeze(0)
            if protein_chemistry is not None:
                protein_chemistry = protein_chemistry.unsqueeze(0)
            if support_fingerprint is not None:
                support_fingerprint = support_fingerprint.unsqueeze(0)
            if query_fingerprint is not None:
                query_fingerprint = query_fingerprint.unsqueeze(0)
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
        endpoint, ligand_value, ligand, summary, embed, occupancy = self.encode(
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
        _, query_ligand = torch.split(ligand, (support_count, query_count), 1)
        _, query_occupancy = torch.split(
            occupancy, (support_count, query_count), 1)
        _, query_ligand_value = torch.split(
            ligand_value, (support_count, query_count), 1)
        batch = endpoint.shape[0]
        # Panel-set level: order-invariant mean/max pooling over the query
        # ligands plus the protein summary. Constant across a target's queries.
        panel_mean = query_ligand.mean(1)
        panel_max = query_ligand.amax(1)
        panel_level = self.panel_level(
            torch.cat((summary, panel_mean, panel_max), -1))
        zero_shot = zero_shot + panel_level
        if not adapt or support_count == 0:
            level_shift = zero_shot.new_zeros(batch)
            level_gate = zero_shot.new_zeros(batch)
            level_adjustment = zero_shot.new_zeros(batch, 1)
            residual = support_y.to(device, dtype).new_zeros(batch, support_count)
            transport = torch.zeros_like(zero_shot)
            weight = zero_shot.new_zeros(
                batch, query_count, max(support_count, 1))
            similarity = weight
            evidence = zero_shot.new_zeros(batch)
        else:
            if support_fingerprint is None or query_fingerprint is None:
                raise ValueError(
                    "chemistry-grounded weighting requires support and query "
                    "fingerprints")
            similarity = self.transport_similarity(
                query_fingerprint.to(device, dtype),
                support_fingerprint.to(device, dtype))
            locked = (support_y.to(device, dtype) - support_zero).detach()
            level_shift = locked.mean(-1)
            shrink = self.transport.shrinkage(support_count, locked)
            level_gate = shrink.expand_as(level_shift)
            level_adjustment = shrink * level_shift.unsqueeze(-1)
            residual = locked - level_adjustment
            transport, weight = self.transport(
                support_embed, query_embed, locked, similarity)
            transport = shrink * transport
            evidence = similarity.amax(-1).mean(-1)
        level_baseline = zero_shot + level_adjustment
        prediction = zero_shot + transport
        summary_row = weight.mean(1) if weight.numel() else zero_shot.new_zeros(batch, 1)
        output = QPSMPMetaOutput(
            prediction=prediction,
            additive=query_ligand_value + panel_level,
            ligand_only=query_ligand_value,
            cross_zero_shot=zero_shot,
            level_baseline=level_baseline,
            level_adjustment=level_adjustment,
            sar_adaptation=transport - level_adjustment,
            adaptation=prediction - zero_shot,
            zero_shot=zero_shot,
            task_state=summary_row[:, None].expand(-1, self.contact_types, -1),
            level_shift=level_shift,
            query_basis=query_occupancy,
            support_residual_quotient=residual,
            support_evidence=similarity,
            evidence_score=evidence,
            level_shrinkage=level_gate,
            shape_scale=torch.ones_like(level_gate),
            sar_scale=(weight.amax(-1).mean(-1) if weight.numel() else level_gate),
            support_match_loss=self.dictionary_regularizer(query_occupancy))
        if not unbatched:
            return output
        return QPSMPMetaOutput(*(value.squeeze(0)
            if isinstance(value, Tensor) and value.ndim else value
            for value in output.__dict__.values()))

    @staticmethod
    def transport_similarity(query: Tensor, support: Tensor) -> Tensor:
        intersection = torch.einsum("bqf,bkf->bqk", query, support)
        query_sum = query.sum(-1)[:, :, None]
        support_sum = support.sum(-1)[:, None, :]
        return intersection / (query_sum + support_sum - intersection).clamp_min(1e-6)

    def panel_level_value(self, summary: Tensor, query_ligand: Tensor) -> Tensor:
        panel_mean = query_ligand.mean(1)
        panel_max = query_ligand.amax(1)
        return self.panel_level(torch.cat((summary, panel_mean, panel_max), -1))
