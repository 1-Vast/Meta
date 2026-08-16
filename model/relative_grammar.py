"""Signed reference-query relative transport on the interaction-grammar trunk.

Why this shape, from the Stage 4 admission failure:

* the previous transport multiplied each support residual by a bounded
  query-dependent scalar `rho(q,k)`.  A bounded rescaling can only shrink the
  transported residual toward the target level, which is exactly what the
  evidence showed: MSE fell slightly, no component-bootstrap lower bound was
  positive, and the concordance index fell from 0.647 to 0.571-0.610;
* permuting the support labels leaves `mean(r)` exactly invariant, so the
  permutation contrast isolates the query-specific channel.  It was negative at
  k=2 and k=3, i.e. support identity was not being used at all.

The borrowed idea is PBCNet2.0's Siamese *relative* formulation -- predict the
signed difference between a reference and a query rather than an absolute
correction -- with its geometry removed, because the coverage audit
(`GEOMETRY_COVERAGE_AUDIT.json`) found **zero** BindingDB deployment pairs with
a common-frame protein-ligand complex.  The support-reliability weighting is the
task-difficulty idea from AdaMBind, expressed without MAML, without an inner
loop and without test-time gradients.

Prediction:

```text
f(P,Lq,S) = f0(P,Lq) + s(n) * sum_k w_qk * [ r_k + delta(P, L_k -> Lq) ]
r_k       = y_k - f0(P,L_k)
s(n)      = n / (n + lambda),  s(0) = 0
delta(a->b) = m(e_a, e_b) - m(e_b, e_a)          exactly antisymmetric
```

Nesting and floors:

* `n = 0` returns exactly `f0`;
* `delta == 0` with flat weights returns exactly the shrunken support-mean
  calibration, i.e. the previous safe behaviour is the floor, not a special
  case that must be relearned;
* `delta(a->a) == 0` and `delta(a->b) == -delta(b->a)` hold algebraically.

`delta` is the only new degree of freedom: it is the deviation of the learned
relative operator from the zero-shot endpoint difference.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor

from .interaction_grammar import InteractionGrammarModel
from .qpsmp_meta import QPSMPMetaOutput


class RelativeDifferenceTransport(nn.Module):
    """Antisymmetric reference-query difference operator with label locking."""

    def __init__(self, embed_dim: int, hidden_dim: int = 64,
                 use_reliability: bool = True,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.use_reliability = bool(use_reliability)
        self.key = nn.Linear(embed_dim, hidden_dim, bias=False, dtype=dtype)
        # The readout carries no bias: `delta` is built as `m(a,b) - m(b,a)`,
        # so an output bias cancels exactly and would be an unidentifiable
        # parameter that can never receive gradient.
        self.pair = nn.Sequential(
            nn.Linear(4 * embed_dim, 2 * hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(2 * hidden_dim, 1, bias=False, dtype=dtype))
        self.log_temperature = nn.Parameter(torch.tensor(1.443, dtype=dtype))
        self.log_shrinkage = nn.Parameter(torch.tensor(1.8546, dtype=dtype))
        self.log_reliability = nn.Parameter(torch.tensor(0.5413, dtype=dtype))
        # Start near the safe floor (delta == 0) without zeroing the readout,
        # which would leave the hidden layer without gradient.
        nn.init.normal_(self.pair[-1].weight, std=5e-2)

    def shrinkage(self, support_count: int, reference: Tensor) -> Tensor:
        strength = F.softplus(self.log_shrinkage)
        count = reference.new_tensor(float(support_count))
        return count / (count + strength)

    def delta(self, target_embed: Tensor, source_embed: Tensor) -> Tensor:
        """`delta[b, i, j]` is the signed shift from `source[j]` to `target[i]`."""
        rows = target_embed.shape[1]
        columns = source_embed.shape[1]
        left = target_embed[:, :, None, :].expand(-1, -1, columns, -1)
        right = source_embed[:, None, :, :].expand(-1, rows, -1, -1)

        def score(a: Tensor, b: Tensor) -> Tensor:
            return self.pair(torch.cat((a, b, a * b, a - b), -1)).squeeze(-1)

        return score(left, right) - score(right, left)

    def reliability(self, support_embed: Tensor, support_residual: Tensor,
                    logits: Tensor) -> Tensor:
        """Leave-one-out agreement of each support label with the others.

        Under a label permutation the residuals stop agreeing with the
        structural difference operator, so this term falls -- which is what
        makes the permuted-support control structurally worse rather than
        accidentally worse.
        """
        support_count = support_embed.shape[1]
        if support_count < 2:
            return torch.zeros_like(support_residual)
        support_delta = self.delta(support_embed, support_embed)
        mask = torch.eye(support_count, dtype=torch.bool,
                         device=support_embed.device)
        loo = logits.masked_fill(mask, torch.finfo(logits.dtype).min)
        weight = torch.softmax(loo, -1)
        predicted = torch.einsum(
            "bkj,bkj->bk", weight,
            support_residual[:, None, :].expand(-1, support_count, -1)
            + support_delta)
        inconsistency = (support_residual - predicted).abs()
        return -inconsistency / F.softplus(self.log_reliability)

    def forward(self, support_embed: Tensor, query_embed: Tensor,
                support_residual: Tensor) -> tuple[Tensor, Tensor, Tensor, Tensor]:
        support_key = F.normalize(self.key(support_embed), dim=-1)
        query_key = F.normalize(self.key(query_embed), dim=-1)
        temperature = F.softplus(self.log_temperature) + 1.0
        logits = temperature * torch.einsum(
            "bqh,bkh->bqk", query_key, support_key)
        credit = torch.zeros_like(support_residual)
        if self.use_reliability:
            support_logits = temperature * torch.einsum(
                "bkh,bjh->bkj", support_key, support_key)
            credit = self.reliability(
                support_embed, support_residual, support_logits)
            logits = logits + credit[:, None, :]
        weight = torch.softmax(logits, -1)
        query_delta = self.delta(query_embed, support_embed)
        transport = torch.einsum(
            "bqk,bqk->bq", weight,
            support_residual[:, None, :].expand_as(query_delta) + query_delta)
        return transport, query_delta, weight, credit


class RelativeGrammarModel(InteractionGrammarModel):
    """Interaction-grammar trunk with signed relative-difference transport.

    The trunk, encoders and zero-shot endpoint are inherited unchanged from
    `InteractionGrammarModel`, so `--arch grammar` remains an exact comparator
    for everything except the transport.
    """

    def __init__(self, *args, use_reliability: bool = True, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        del self.transport
        self.transport = RelativeDifferenceTransport(
            self.embed_dim, kwargs.get("task_dim", 64),
            use_reliability=use_reliability,
            dtype=next(self.parameters()).dtype)

    def query_difference(self, protein_pooled: Tensor, protein_tokens: Tensor,
                         protein_mask: Tensor, query_atoms: Tensor,
                         query_bonds: Tensor, query_mask: Tensor,
                         protein_chemistry: Tensor | None = None
                         ) -> tuple[Tensor, Tensor]:
        """Signed within-episode query differences and their endpoint baseline.

        Used only as a training signal on query labels, never as an input.
        Returns `(delta[b,i,j], endpoint[b,i])`.
        """
        parameter = next(self.parameters())
        device, dtype = parameter.device, parameter.dtype
        endpoint, _, _, embed, _ = self.encode(
            protein_pooled.to(device, dtype), protein_tokens.to(device, dtype),
            protein_mask.to(device, dtype), query_atoms.to(device, dtype),
            query_bonds.to(device, dtype), query_mask.to(device, dtype),
            None if protein_chemistry is None
            else protein_chemistry.to(device, dtype))
        return self.transport.delta(embed, embed), endpoint

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
            raise ValueError("the relative trunk does not accept transplanted states")
        if geometry_coordinates is not None or geometry_edge_index is not None or (
                geometry_available is not None and bool(geometry_available.any())):
            raise ValueError(
                "no BindingDB deployment pair has a common-frame complex; "
                "see report/.../GEOMETRY_COVERAGE_AUDIT.json")
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
            query_delta = zero_shot.new_zeros(
                batch, query_count, max(support_count, 1))
            credit = residual
            evidence = zero_shot.new_zeros(batch)
        else:
            locked = (support_y.to(device, dtype) - support_zero).detach()
            level_shift = locked.mean(-1)
            shrink = self.transport.shrinkage(support_count, locked)
            level_gate = shrink.expand_as(level_shift)
            level_adjustment = shrink * level_shift.unsqueeze(-1)
            residual = locked - level_adjustment
            transport, query_delta, _, credit = self.transport(
                support_embed, query_embed, locked)
            transport = shrink * transport
            evidence = query_delta.square().mean((-2, -1)).sqrt()
        level_baseline = zero_shot + level_adjustment
        prediction = zero_shot + transport
        summary = (query_delta.mean(1) if query_delta.numel()
                   else zero_shot.new_zeros(batch, 1))
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
            support_evidence=query_delta,
            evidence_score=evidence,
            level_shrinkage=level_gate,
            shape_scale=torch.ones_like(level_gate),
            sar_scale=(credit.mean(-1) if credit.numel() else level_gate),
            support_match_loss=self.dictionary_regularizer(query_occupancy))
        if not unbatched:
            return output
        return QPSMPMetaOutput(*(value.squeeze(0)
            if isinstance(value, Tensor) and value.ndim else value
            for value in output.__dict__.values()))
