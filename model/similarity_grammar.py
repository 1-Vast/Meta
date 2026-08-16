"""Chemically grounded support weighting on the interaction-grammar trunk.

Evidence that motivates this exact change
(`report/meta_fewshot/stage6_bottleneck_20260815/SIGNAL_meta_val.json`):
with no learning at all, a Tanimoto-similarity-weighted average of the support
labels beats their plain mean on `meta_val` by

| k | support mean | Tanimoto-softmax | oracle best single support |
|---|---:|---:|---:|
| 2 | 1.163 | 0.976 | 0.570 |
| 3 | 1.096 | 0.887 | 0.357 |
| 5 | 1.002 | 0.747 | 0.186 |

and the within-target correlation between Tanimoto similarity and absolute
affinity gap is -0.35. Every learned transport in this project instead collapsed
to uniform weights. The bottleneck was therefore the *similarity representation*
used for support weighting, not the training objective: a 1024-bit Morgan
fingerprint captures what the learned kernel keys did not find.

This is a **fixed Morgan/Tanimoto support kernel with learned scalar
calibration**, not a non-learned estimator: `similarity_scale` and
`log_shrinkage` are trained. In practice `similarity_scale` moves very little
(7.98-7.99 against an initialisation of 8.0), so almost all of the behaviour is
carried by the fixed kernel.

Transport:

```text
f = f0(P,Lq) + s(n) * sum_k w_qk * r_k
w_qk = softmax_k( tau <key_q,key_k> + gamma * Tanimoto(Lq, L_k) )
```

`gamma -> 0` and `tau -> 0` recover the shrunken support mean exactly, so the
previous behaviour remains the floor. No bounded rescaling gate is used: Stage 4
showed that channel trades ranking for squared error. At k=1 the weighting is
degenerate by construction and this mechanism makes no k=1 claim.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor

from .interaction_grammar import InteractionGrammarModel
from .qpsmp_meta import QPSMPMetaOutput


def tanimoto(query: Tensor, support: Tensor) -> Tensor:
    """Tanimoto similarity between binary fingerprint rows: [B,Q,F] x [B,K,F]."""
    intersection = torch.einsum("bqf,bkf->bqk", query, support)
    query_sum = query.sum(-1)[:, :, None]
    support_sum = support.sum(-1)[:, None, :]
    return intersection / (query_sum + support_sum - intersection).clamp_min(1e-6)


class SimilarityTransport(nn.Module):
    """Label-locked residual transport weighted by learned + chemical similarity."""

    def __init__(self, embed_dim: int, hidden_dim: int = 64,
                 similarity_scale: float = 8.0, use_learned_key: bool = True,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.use_learned_key = bool(use_learned_key)
        self.key = nn.Linear(embed_dim, hidden_dim, bias=False, dtype=dtype)
        self.log_temperature = nn.Parameter(torch.tensor(1.443, dtype=dtype))
        self.log_shrinkage = nn.Parameter(torch.tensor(1.8546, dtype=dtype))
        # Initialised at the value the label-only audit used, so training starts
        # inside the regime the evidence says is useful instead of searching for
        # it from scratch.
        self.similarity_scale = nn.Parameter(
            torch.tensor(float(similarity_scale), dtype=dtype))
        if not self.use_learned_key:
            # `key` and `log_temperature` are constructed unconditionally so
            # that state-dict keys stay compatible across both variants, but
            # they take no part in the forward pass here. Freezing them keeps
            # them out of the optimizer and out of the trainable-parameter
            # census instead of leaving silently dead trainable tensors.
            self.key.weight.requires_grad_(False)
            self.log_temperature.requires_grad_(False)

    def shrinkage(self, support_count: int, reference: Tensor) -> Tensor:
        strength = F.softplus(self.log_shrinkage)
        count = reference.new_tensor(float(support_count))
        return count / (count + strength)

    def forward(self, support_embed: Tensor, query_embed: Tensor,
                support_residual: Tensor, similarity: Tensor
                ) -> tuple[Tensor, Tensor]:
        logits = self.similarity_scale * similarity
        if self.use_learned_key:
            support_key = F.normalize(self.key(support_embed), dim=-1)
            query_key = F.normalize(self.key(query_embed), dim=-1)
            temperature = F.softplus(self.log_temperature) + 1.0
            logits = logits + temperature * torch.einsum(
                "bqh,bkh->bqk", query_key, support_key)
        weight = torch.softmax(logits, -1)
        transport = torch.einsum("bqk,bk->bq", weight, support_residual)
        return transport, weight


class SimilarityGrammarModel(InteractionGrammarModel):
    """Interaction-grammar trunk with chemistry-grounded support weighting.

    The trunk, encoders and zero-shot endpoint are inherited unchanged, so
    `--arch grammar` remains an exact comparator for everything but the
    transport.
    """

    def __init__(self, *args, similarity_scale: float = 8.0,
                 use_learned_key: bool = True, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        del self.transport
        self.transport = SimilarityTransport(
            self.embed_dim, kwargs.get("task_dim", 64),
            similarity_scale=similarity_scale, use_learned_key=use_learned_key,
            dtype=next(self.parameters()).dtype)

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
            raise ValueError("the similarity trunk does not accept transplanted states")
        if geometry_coordinates is not None or geometry_edge_index is not None or (
                geometry_available is not None and bool(geometry_available.any())):
            raise ValueError(
                "no BindingDB deployment pair has a common-frame complex; "
                "see stage5 GEOMETRY_COVERAGE_AUDIT.json")
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
            weight = zero_shot.new_zeros(
                batch, query_count, max(support_count, 1))
            similarity = weight
            evidence = zero_shot.new_zeros(batch)
        else:
            if support_fingerprint is None or query_fingerprint is None:
                raise ValueError(
                    "chemistry-grounded weighting requires support and query "
                    "fingerprints; rebuild episodes with QPSMPData.materialize")
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
            # Reported mechanism signal: how far the weighting is from uniform.
            sar_scale=(weight.amax(-1).mean(-1) if weight.numel() else level_gate),
            support_match_loss=self.dictionary_regularizer(query_occupancy))
        if not unbatched:
            return output
        return QPSMPMetaOutput(*(value.squeeze(0)
            if isinstance(value, Tensor) and value.ndim else value
            for value in output.__dict__.values()))
