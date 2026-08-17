""""Stage K candidate: contrastive coembedding branch on the incumbent trunk.

The zero-shot path, Tanimoto transport and level branch are byte-identical
to the similarity_only incumbent. A coembedding branch projects the protein
summary and each ligand encoding to 128 dims; the training objective (added
in the trainer) is episodic InfoNCE with other query ligands as negatives.
"""
from __future__ import annotations

import torch
import torch.nn as nn
from torch import Tensor

from model.similarity_grammar import SimilarityGrammarModel, SimilarityTransport
from model.qpsmp_meta import QPSMPMetaOutput


class CoembedModel(SimilarityGrammarModel):
    def __init__(self, *args, coembed_dim: int = 128, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        del self.transport
        self.transport = SimilarityTransport(
            self.embed_dim, kwargs.get("task_dim", 64),
            similarity_scale=8.0, use_learned_key=False,
            dtype=next(self.parameters()).dtype)
        hidden_dim = int(kwargs.get("hidden_dim", 192))
        self.coembed_dim = int(coembed_dim)
        self.zp = nn.Linear(hidden_dim, self.coembed_dim,
                            dtype=next(self.parameters()).dtype)
        self.zl = nn.Linear(hidden_dim, self.coembed_dim,
                            dtype=next(self.parameters()).dtype)
        self._stash = None

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
        protein_value = self.protein_head(wide_summary).squeeze(-1)
        interaction = (self.interaction_head(
            torch.cat((embed, section), -1)).squeeze(-1)
            + self.contact_weight(occupancy).squeeze(-1))
        endpoint = ligand_value + protein_value + interaction
        shape = (batch, count)
        return (endpoint.reshape(shape), ligand_value.reshape(shape),
                protein_value.reshape(shape), ligand.reshape(*shape, -1),
                summary, embed.reshape(*shape, -1),
                occupancy.reshape(*shape, -1))

    def forward(self, protein_pooled, protein_tokens, protein_mask,
                support_atoms, support_bonds, support_mask, support_y,
                query_atoms, query_bonds, query_mask, *, adapt=True,
                protein_chemistry=None, support_fingerprint=None,
                query_fingerprint=None, task_state_override=None,
                geometry_coordinates=None, geometry_edge_index=None,
                geometry_available=None, geometry_common_frame=None):
        if task_state_override is not None:
            raise ValueError("the coembed trunk rejects transplanted states")
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
        (endpoint, ligand_value, protein_value, ligand, summary, embed,
         occupancy) = self.encode(
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
        _, query_protein_value = torch.split(
            protein_value, (support_count, query_count), 1)
        batch = endpoint.shape[0]
        # coembedding branch: stashed for the trainer's contrastive loss
        zp = self.zp(summary)                       # [B, D]
        zl = self.zl(query_ligand)                  # [B, Q, D]
        self._stash = (zp, zl)
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
            level_shrinkage=level_gate,
            shape_scale=torch.ones_like(level_gate),
            sar_scale=(weight.amax(-1).mean(-1) if weight.numel() else level_gate),
            support_match_loss=self.dictionary_regularizer(query_occupancy))
        if not unbatched:
            return output
        return QPSMPMetaOutput(*(value.squeeze(0)
            if isinstance(value, Tensor) and value.ndim else value
            for value in output.__dict__.values()))
