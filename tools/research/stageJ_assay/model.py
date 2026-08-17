""""Stage J candidate: assay-aware level head with paired level alignment.

Framework innovation (I1j): the level head consumes THREE legal covariate
families - the protein summary, order-invariant mean/max pooling over the
query ligand encodings (panel composition), and a learned embedding of the
query panel's journal/publisher codes (assay provenance, parsed from
panel_ids). D0c measured the journal family at level MSE 1.619 (vs 2.155
constant, shuffled control 2.522) and 100% of meta_val episodes share a
journal code with meta_train, so the feature is usable at inference.

Training innovation (I2j): paired cross-target level alignment. Within each
optimization step, the predicted level gap between two episodes is regressed
against the true (transport-residual) level gap - a BatchDTA-inspired
implicit alignment that pins a consistent global level scale without using
any query-label statistic at inference. The full-prediction smooth_l1 term
is retained so the head learns its RESIDUAL role alongside the transport
(the Stage E double-fit failure mode is thereby removed by construction).

Journal embeddings are trained in the same single stage (ordinary
parameters). Query labels remain loss-only.
"""
from __future__ import annotations

import torch
import torch.nn as nn
from torch import Tensor

from model.similarity_grammar import SimilarityGrammarModel, SimilarityTransport
from model.qpsmp_meta import QPSMPMetaOutput


class AssayLevelModel(SimilarityGrammarModel):
    """Incumbent trunk + Tanimoto transport + assay-aware level head."""

    def __init__(self, *args, journal_vocab: int = 0, level_hidden: int = 64,
                 journal_dim: int = 32, use_learned_key: bool = False,
                 **kwargs) -> None:
        super().__init__(*args, **kwargs)
        del self.protein_head
        del self.transport
        self.transport = SimilarityTransport(
            self.embed_dim, kwargs.get("task_dim", 64),
            similarity_scale=8.0, use_learned_key=use_learned_key,
            dtype=next(self.parameters()).dtype)
        hidden_dim = int(kwargs.get("hidden_dim", 192))
        self.journal_vocab = int(journal_vocab)
        self.journal_dim = int(journal_dim)
        self.journal_embed = nn.Embedding(
            max(self.journal_vocab, 1), self.journal_dim,
            dtype=next(self.parameters()).dtype)
        nn.init.normal_(self.journal_embed.weight, std=0.05)
        # level input: summary [H] + panel mean [H] + panel max [H] + journal [D]
        self.panel_level = nn.Sequential(
            nn.Linear(3 * hidden_dim + self.journal_dim, level_hidden,
                      dtype=next(self.parameters()).dtype),
            nn.GELU(),
            nn.Linear(level_hidden, 1, dtype=next(self.parameters()).dtype))

    def encode(self, protein_pooled, protein_tokens, protein_mask, raw_atoms,
               bonds, mask, protein_chemistry):
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

    def level_value(self, summary, query_ligand, journal_ids):
        panel_mean = query_ligand.mean(1)
        panel_max = query_ligand.amax(1)
        if journal_ids is None or self.journal_vocab == 0:
            journal = panel_mean.new_zeros(summary.shape[0], self.journal_dim)
        else:
            ids = journal_ids.clamp(min=0)
            mask = (journal_ids >= 0).float().unsqueeze(-1)
            embedded = self.journal_embed(ids) * mask
            count = mask.sum(1).clamp_min(1.0)
            journal = embedded.sum(1) / count
        return self.panel_level(torch.cat(
            (summary, panel_mean, panel_max, journal), -1))

    def forward(self, protein_pooled, protein_tokens, protein_mask,
                support_atoms, support_bonds, support_mask, support_y,
                query_atoms, query_bonds, query_mask, *, adapt=True,
                protein_chemistry=None, support_fingerprint=None,
                query_fingerprint=None, journal_ids=None, level_gate=None,
                task_state_override=None, geometry_coordinates=None,
                geometry_edge_index=None, geometry_available=None,
                geometry_common_frame=None):
        if task_state_override is not None:
            raise ValueError("the assay-level trunk rejects transplanted states")
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
            if journal_ids is not None:
                journal_ids = journal_ids.unsqueeze(0)
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
        panel_level = self.level_value(summary, query_ligand, journal_ids)
        if level_gate is not None:
            panel_level = panel_level * float(level_gate)
        zero_shot = zero_shot + panel_level
        if not adapt or support_count == 0:
            level_shift = zero_shot.new_zeros(batch)
            level_gate = zero_shot.new_zeros(batch)
            level_adjustment = zero_shot.new_zeros(batch, 1)
            residual = support_y.to(device, dtype).new_zeros(batch, support_count)
            transport = torch.zeros_like(zero_shot)
            weight = zero_shot.new_zeros(batch, query_count, max(support_count, 1))
            similarity = weight
            evidence = zero_shot.new_zeros(batch)
        else:
            if support_fingerprint is None or query_fingerprint is None:
                raise ValueError("transport requires support and query fingerprints")
            from model.similarity_grammar import tanimoto
            similarity = tanimoto(query_fingerprint.to(device, dtype),
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
