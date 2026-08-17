"""Stage F candidate: pairwise learned interaction transport.

Framework innovation (I1f): the transport kernel is a learned pairwise
operator over (query, support) pairs. For each pair, a small edge MLP consumes
the query embedding, the support embedding and the support residual value and
emits a logit; the fixed Morgan/Tanimoto kernel remains as an additive anchor
(learned scale), so the collapse-to-uniform failure mode of every prior
learned kernel keeps a chemically grounded floor. The transport stays
label-locked: only support residual VALUES are transported, never query
labels.

Motivation (measured): the Stage L audit found a pairwise signed-gap direction
in `embed` (r +0.270 [+0.128, +0.418]) that is orthogonal to Tanimoto
(correlation +0.026) and to every moment-form adaptation tested in A2. A
pairwise operator is the smallest mechanism that can consume it.

Training innovation (I2f): pairwise signed-gap supervision. Within each
episode the predicted affinity gaps p(q) - p(k) are regressed against the
true signed gaps y(q) - y(k) over the query x support grid. Query labels
remain loss-only; inference consumes only support labels.

The zero-shot trunk, protein level branch, contact dictionary and dictionary
regularizer are byte-identical to the incumbent similarity_only model: the
transport and its supervision are the only differences.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor

from model.similarity_grammar import SimilarityGrammarModel, tanimoto
from model.qpsmp_meta import QPSMPMetaOutput


class PairwiseTransport(nn.Module):
    """Label-locked transport with a learned pairwise edge logit."""

    def __init__(self, embed_dim: int, hidden_dim: int = 64, dtype=None) -> None:
        super().__init__()
        self.query_proj = nn.Linear(embed_dim, hidden_dim, bias=False, dtype=dtype)
        self.support_proj = nn.Linear(embed_dim, hidden_dim, bias=False, dtype=dtype)
        self.edge = nn.Sequential(
            nn.Linear(2 * hidden_dim + 1, hidden_dim, dtype=dtype),
            nn.GELU(),
            nn.Linear(hidden_dim, 1, dtype=dtype))
        self.log_temperature = nn.Parameter(torch.tensor(1.443, dtype=dtype))
        self.log_shrinkage = nn.Parameter(torch.tensor(1.8546, dtype=dtype))
        self.similarity_scale = nn.Parameter(torch.tensor(8.0, dtype=dtype))
        # start with the edge nearly silent so training begins inside the
        # regime where the fixed Tanimoto kernel already works
        nn.init.normal_(self.edge[-1].weight, std=1e-2)
        nn.init.zeros_(self.edge[-1].bias)

    def shrinkage(self, support_count: int, reference: Tensor) -> Tensor:
        strength = F.softplus(self.log_shrinkage)
        count = reference.new_tensor(float(support_count))
        return count / (count + strength)

    def forward(self, support_embed: Tensor, query_embed: Tensor,
                support_residual: Tensor, similarity: Tensor
                ) -> tuple[Tensor, Tensor]:
        q = self.query_proj(query_embed)[:, :, None, :]      # B,Q,1,H
        s = self.support_proj(support_embed)[:, None, :, :]  # B,1,K,H
        r = support_residual[:, None, :, None]               # B,1,K,1
        edge = self.edge(torch.cat((q.expand(-1, -1, s.shape[2], -1),
                                    s.expand(-1, q.shape[1], -1, -1),
                                    r.expand(-1, q.shape[1], -1, -1)),
                                   dim=-1)).squeeze(-1)      # B,Q,K
        logits = self.similarity_scale * similarity + edge
        weight = torch.softmax(logits, -1)
        transport = torch.einsum("bqk,bk->bq", weight, support_residual)
        return transport, weight


class PairwiseTransportModel(SimilarityGrammarModel):
    """Incumbent trunk with the pairwise learned transport."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        del self.transport
        self.transport = PairwiseTransport(
            self.embed_dim, kwargs.get("task_dim", 64),
            dtype=next(self.parameters()).dtype)

    def forward(self, protein_pooled, protein_tokens, protein_mask,
                support_atoms, support_bonds, support_mask, support_y,
                query_atoms, query_bonds, query_mask, *, adapt: bool = True,
                protein_chemistry=None, support_fingerprint=None,
                query_fingerprint=None, task_state_override=None,
                geometry_coordinates=None, geometry_edge_index=None,
                geometry_available=None, geometry_common_frame=None):
        if task_state_override is not None:
            raise ValueError("the pairwise trunk does not accept transplanted states")
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
        endpoint, ligand_value, protein_value, embed, occupancy = self.encode(
            protein_pooled.to(device, dtype), protein_tokens.to(device, dtype),
            protein_mask.to(device, dtype), raw_atoms, bonds, mask,
            None if protein_chemistry is None else protein_chemistry.to(device, dtype))
        support_count = support_atoms.shape[1]
        query_count = endpoint.shape[1] - support_count
        support_zero, zero_shot = torch.split(endpoint, (support_count, query_count), 1)
        support_embed, query_embed = torch.split(embed, (support_count, query_count), 1)
        _, query_occupancy = torch.split(occupancy, (support_count, query_count), 1)
        _, query_ligand_value = torch.split(ligand_value, (support_count, query_count), 1)
        _, query_protein_value = torch.split(protein_value, (support_count, query_count), 1)
        batch = endpoint.shape[0]
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
            similarity = tanimoto(
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
        summary = weight.mean(1) if weight.numel() else zero_shot.new_zeros(batch, 1)
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
            task_state=summary[:, None].expand(-1, self.contact_types, -1),
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
