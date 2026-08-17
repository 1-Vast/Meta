"""Stage Q model: decoupled frozen-feature level head on the incumbent trunk."""
from __future__ import annotations

import torch
import torch.nn as nn
from torch import Tensor

from model.similarity_grammar import SimilarityGrammarModel, SimilarityTransport
from model.qpsmp_meta import QPSMPMetaOutput


class FrozenHeadModel(SimilarityGrammarModel):
    """Incumbent trunk + Tanimoto transport + decoupled frozen-feature level head.

    The head consumes ONLY frozen features: the frozen ESM pooled vector
    (the bank vector already carried by the episode), handcrafted panel
    statistics and a trainable journal-embedding table. No trunk-derived
    representation enters the head, so its gradients cannot reshape the
    trunk (the measured Stage L failure mode).
    """

    def __init__(self, *args, journal_vocab: int = 0, panel_dim: int = 35,
                 level_hidden: int = 64, journal_dim: int = 32,
                 use_learned_key: bool = False, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        del self.transport
        self.transport = SimilarityTransport(
            self.embed_dim, kwargs.get("task_dim", 64),
            similarity_scale=8.0, use_learned_key=use_learned_key,
            dtype=next(self.parameters()).dtype)
        protein_dim = int(args[0]) if args else int(kwargs["protein_dim"])
        self.journal_vocab = int(journal_vocab)
        self.journal_dim = int(journal_dim)
        self.panel_dim = int(panel_dim)
        self.journal_embed = nn.Embedding(
            max(self.journal_vocab, 1), self.journal_dim,
            dtype=next(self.parameters()).dtype)
        nn.init.normal_(self.journal_embed.weight, std=0.05)
        # inputs: frozen ESM pooled [protein_dim] + panel stats [panel_dim]
        # + journal [D]
        input_dim = int(protein_dim) + self.panel_dim + self.journal_dim
        self.level_head = nn.Sequential(
            nn.Linear(input_dim, level_hidden,
                      dtype=next(self.parameters()).dtype),
            nn.GELU(),
            nn.Linear(level_hidden, 1, dtype=next(self.parameters()).dtype))

    def level_value(self, esm_pooled, panel_stats, journal_ids):
        if journal_ids is None or self.journal_vocab == 0:
            journal = esm_pooled.new_zeros(esm_pooled.shape[0], self.journal_dim)
        else:
            ids = journal_ids.clamp(min=0)
            mask = (journal_ids >= 0).float().unsqueeze(-1)
            embedded = self.journal_embed(ids) * mask
            count = mask.sum(1).clamp_min(1.0)
            journal = embedded.sum(1) / count
        x = torch.cat((esm_pooled, panel_stats, journal), -1)
        return self.level_head(x)

    def forward(self, protein_pooled, protein_tokens, protein_mask,
                support_atoms, support_bonds, support_mask, support_y,
                query_atoms, query_bonds, query_mask, *, adapt=True,
                protein_chemistry=None, support_fingerprint=None,
                query_fingerprint=None, journal_ids=None, panel_stats=None,
                level_gate=None, task_state_override=None,
                geometry_coordinates=None, geometry_edge_index=None,
                geometry_available=None, geometry_common_frame=None):
        if task_state_override is not None:
            raise ValueError("the frozen-head trunk rejects transplanted states")
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
            if panel_stats is not None:
                panel_stats = panel_stats.unsqueeze(0)
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
        if panel_stats is None:
            panel_stats = query_ligand_value.new_zeros(batch, self.panel_dim)
        panel_stats = panel_stats.to(device, dtype)
        if journal_ids is not None:
            journal_ids = journal_ids.to(device)
        head_level = self.level_value(
            protein_pooled.to(device, dtype), panel_stats, journal_ids)
        if level_gate is not None:
            head_level = head_level * float(level_gate)
        zero_shot = zero_shot + head_level
        if not adapt or support_count == 0:
            level_shift = zero_shot.new_zeros(batch)
            level_gate_out = zero_shot.new_zeros(batch)
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
            level_gate_out = shrink.expand_as(level_shift)
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
            additive=query_ligand_value + query_protein_value,
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
            level_shrinkage=level_gate_out,
            shape_scale=torch.ones_like(level_gate_out),
            sar_scale=(weight.amax(-1).mean(-1) if weight.numel() else level_gate_out),
            support_match_loss=self.dictionary_regularizer(query_occupancy))
        if not unbatched:
            return output
        return QPSMPMetaOutput(*(value.squeeze(0)
            if isinstance(value, Tensor) and value.ndim else value
            for value in output.__dict__.values()))
