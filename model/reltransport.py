"""Core Innovation A: a protein-conditioned relative interaction potential
shared by the zero-shot shape and the few-shot transport.

Endpoint (k=0):

    f0(P, L) = ligand_prior(L) + target_level(P) + shape(P, L)
    shape(P, L) = mean_m delta(P, e(P, L), anchor_m)

Transport (k >= 1):

    t(q) = shrink(n) * sum_k a(q, k) * rho(q, k) * r_k
    a(q, k) = softmax_k(gamma * Tanimoto(q, k)
                        + tau * cos(key(e_q), key(e_k)))
    rho(q, k) = 1 + u_g(P)^T [g(e_q) * h(e_k) + g(e_k) * h(e_q)]
    r_k = y_k - f0(P, L_k)                (label locked: detached)

`rho` is a linear query-specific transferability gate, zero-initialised so
the model starts exactly at the shrunken support mean and learns per-pair
residual scaling — active at k=1, where the softmax is degenerate. The
relative potential `delta` provides the zero-shot shape and is supervised
on within-target label differences; the gate's pair features live in the
same interaction-embedding space.

`delta` is an antisymmetric relative interaction potential:

    delta(P, i, j) = u(P)^T [psi(e_i) * phi(e_j) - psi(e_j) * phi(e_i)]

built from the per-ligand interaction embedding `e(P, L)` of the grammar
trunk (atom-to-residue cross attention). `psi`, `phi` are learned linear maps,
`u(P)` is a protein-conditioned direction vector, and `anchor_1..M` are learned
points in interaction-embedding space — model parameters, never episode data.

Why this exact structure
------------------------
Stage R0 measured the binding failure on activity cliffs: the trained
zero-shot endpoint orders at chance (0.519) where a parameter-free chemical
prior reaches 0.716, and the Stage 9/R3 decomposition shows the interaction
branch contributes essentially nothing to within-target shape. The root cause
is not capacity (Stage 0: the trunk *can* express the interaction) — it is
that nothing forces the endpoint to be a *relative* predictor, and a scalar
head trained on absolute affinity buys nearly all of its squared error from
the target level, which is constant across the queries of one target.

This module makes relative structure the parameterization itself:

* the endpoint's within-target differences ARE the delta predictions:
  shape(P,L_i) - shape(P,L_j) = delta(P,L_i,L_j) exactly, so zero-shot
  ordering ability is mechanically the quality of `delta` — it cannot hide in
  a target-level constant;
* `delta` is antisymmetric by construction (delta(i,j) = -delta(j,i), and
  delta(i,i) = 0), the PBCNet-style relative-recognition inductive bias,
  without any 3D input;
* the anchor centering pins the shape branch: `mean_m shape(P, anchor_m) = 0`
  exactly for every protein, so the branch has no constant component and the
  factorization level/shape is structural, not cosmetic;
* the transport is a query-specific residual gate, so one-support episodes
  already produce a query-specific correction:
  correction(q,s,P) = shrink * rho(P,L_q,L_s) * r_s, which
  depends on (protein, support ligand, support label, query ligand) jointly;
* `gamma -> 0, tau -> 0, delta -> 0` recovers the shrunken support mean
  (level-only abstention stays inside the hypothesis class), and
  `tau -> 0, delta -> 0` recovers the Stage 6/7 fixed Morgan/Tanimoto
  transport, retained here strictly as a baseline special case.

Contracts preserved
-------------------
* one scalar endpoint; `adapt=False` or k=0 returns exactly the zero-shot
  endpoint;
* support labels enter only as residuals against that endpoint (label locked);
* support ordering does not change the output; queries are equivariant;
* query labels are never a model input;
* no common-frame coordinate path (none exists in this corpus).
"""
from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch import Tensor

from .encoders import LigandEncoder
from .interaction_grammar import ContactGrammar, ResidueEncoder
from .similarity_grammar import tanimoto


@dataclass
class RelTransportOutput:
    """Every field is a scalar-per-query tensor unless noted."""
    prediction: Tensor           # deployed output: endpoint + transport
    endpoint: Tensor             # zero-shot endpoint f(P, L)
    ligand_prior: Tensor
    target_level: Tensor         # constant across the queries of one protein
    shape: Tensor                # anchor-centered relative interaction
    transport: Tensor            # few-shot correction, exactly 0 at k=0
    support_endpoint: Tensor     # endpoint on the support ligands
    support_residual: Tensor     # y_k - endpoint(L_k), label locked
    weight: Tensor               # [B, Q, K] transport weights
    rho: Tensor                  # [B, Q, K] transferability gate, 0 at k=0
    similarity: Tensor           # [B, Q, K] Tanimoto
    delta: Tensor                # [B, Q, K] relative potential (0 at k=0)
    anchor_shape_mean: Tensor    # [B] mean_m shape(P, anchor_m) == 0 exactly


class RelativePotential(nn.Module):
    """Antisymmetric protein-conditioned relative interaction potential.

    ``delta(P, i, j) = u(P)^T [psi(e_i) * phi(e_j) - psi(e_j) * phi(e_i)]``
    — a bilinear antisymmetric form over the interaction embeddings, so
    antisymmetry is exact by construction and, crucially, the zero-shot
    shape's within-target differences ARE the delta predictions:

        shape(P, L_i) - shape(P, L_j) = delta(P, L_i, L_j)   (exactly)

    The R6 screening falsified the nonlinear joint pair function on this
    identity's absence: its anchor-shape was invisible to the relative
    supervision (anchors do not appear in delta(i,j)), the gate measured
    inert (nogate gap 0.000) and the endpoint spread collapsed to 0.053 pK.
    The bilinear form makes the relative supervision a *direct* supervision
    of the endpoint's own pairwise differences.

    Consequence recorded honestly: under the bilinear form the anchor set
    acts through its two summary vectors ``psi(mean anchor)`` and
    ``phi(mean anchor)``, so the M anchors parameterize two learned reference
    directions rather than M independent reference ligands. Every anchor
    still receives gradient (through the mean), the exact anchor-mean-zero
    property is preserved, and no constant can enter the shape branch.
    """

    def __init__(self, embed_dim: int, rank: int, hidden_dim: int,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.psi = nn.Linear(embed_dim, rank, bias=False, dtype=dtype)
        self.phi = nn.Linear(embed_dim, rank, bias=False, dtype=dtype)
        self.direction = nn.Sequential(
            nn.LayerNorm(2 * hidden_dim, dtype=dtype),
            nn.Linear(2 * hidden_dim, rank, bias=False, dtype=dtype))

    def direction_vector(self, summary: Tensor,
                         residue_mean: Tensor) -> Tensor:
        """u(P) in R^rank: which relative direction this protein cares about."""
        return self.direction(torch.cat((summary, residue_mean), -1))

    def delta_matrix(self, u: Tensor, embed_i: Tensor,
                     embed_j: Tensor) -> Tensor:
        """[B,R] x [B,Q,e] x [B,K,e] -> [B,Q,K]."""
        forward = torch.einsum(
            "br,bqr,bkr->bqk", u, self.psi(embed_i), self.phi(embed_j))
        reverse = torch.einsum(
            "br,bkr,bqr->bqk", u, self.psi(embed_j), self.phi(embed_i))
        return forward - reverse


class SymmetricGate(nn.Module):
    """Linear query-specific transferability gate for the residual transport.

    ``rho(q, k) = 1 + u(P)^T [g(e_q) * h(e_k) + g(e_k) * h(e_q)]`` with
    zero-initialised `g`/`h`, so the model starts exactly at the shrunken
    support mean (rho == 1) and learns per-(query, support) scaling of the
    label-locked residual. The R6a screening falsified the saturating form
    `1 + tanh(delta)`: a bounded gate cannot express the optimal per-query
    residual scaling. A linear gate can, and its zero initialisation keeps
    level-only abstention at the starting point of optimization, not a limit.
    """

    def __init__(self, embed_dim: int, rank: int, hidden_dim: int,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.g = nn.Linear(embed_dim, rank, bias=False, dtype=dtype)
        self.h = nn.Linear(embed_dim, rank, bias=False, dtype=dtype)
        self.direction = nn.Sequential(
            nn.LayerNorm(2 * hidden_dim, dtype=dtype),
            nn.Linear(2 * hidden_dim, rank, bias=False, dtype=dtype))
        nn.init.normal_(self.g.weight, std=1e-3)
        nn.init.normal_(self.h.weight, std=1e-3)

    def direction_vector(self, summary: Tensor,
                         residue_mean: Tensor) -> Tensor:
        return self.direction(torch.cat((summary, residue_mean), -1))

    def gate_matrix(self, u: Tensor, query_embed: Tensor,
                    support_embed: Tensor) -> Tensor:
        """[B,R] x [B,Q,e] x [B,K,e] -> [B,Q,K], 1 at initialisation."""
        forward = torch.einsum(
            "br,bqr,bkr->bqk", u, self.g(query_embed), self.h(support_embed))
        reverse = torch.einsum(
            "br,bkr,bqr->bqk", u, self.g(support_embed), self.h(query_embed))
        return 1.0 + forward + reverse


class RelativeTransport(nn.Module):
    """Similarity-weighted, query-specific residual transport.

    ``t(q) = shrink(n) * sum_k a(q,k) * rho(q,k) * r_k``: the fixed
    Tanimoto kernel plus a learned key weight the supports; the linear gate
    rescales each label-locked residual per (query, support) pair — active at
    k=1, where the softmax is degenerate by construction. `rho == 1` with a
    flat weight recovers the shrunken support mean (level-only abstention),
    and `rho == 1` alone recovers the Stage 6/7 Tanimoto transport, retained
    strictly as a baseline special case. The R6b additive correction was
    falsified here: `r_k + delta_hat - delta_f0` self-cancels as the relative
    potential converges to the endpoint's implied relative, which the
    measured inertness (nogate gap 0.001-0.000 in R6b/R6c) confirmed.
    """

    def __init__(self, embed_dim: int, similarity_scale: float = 8.0,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.similarity_scale = nn.Parameter(
            torch.tensor(float(similarity_scale), dtype=dtype))
        self.log_temperature = nn.Parameter(torch.tensor(0.0, dtype=dtype))
        self.log_shrinkage = nn.Parameter(torch.tensor(1.8546, dtype=dtype))
        self.key = nn.Linear(embed_dim, embed_dim, bias=False, dtype=dtype)

    def shrinkage(self, support_count: int, reference: Tensor) -> Tensor:
        strength = F.softplus(self.log_shrinkage)
        count = reference.new_tensor(float(support_count))
        return count / (count + strength)

    def forward(self, support_embed: Tensor, query_embed: Tensor,
                support_residual: Tensor, similarity: Tensor,
                rho: Tensor) -> tuple[Tensor, Tensor]:
        logits = self.similarity_scale * similarity
        if support_embed.shape[1] > 1:
            support_key = F.normalize(self.key(support_embed), dim=-1)
            query_key = F.normalize(self.key(query_embed), dim=-1)
            temperature = F.softplus(self.log_temperature) + 1.0
            logits = logits + temperature * torch.einsum(
                "bqh,bkh->bqk", query_key, support_key)
        weight = torch.softmax(logits, -1)
        transport = torch.einsum(
            "bqk,bqk,bk->bq", weight, rho, support_residual)
        return transport, weight


class RelTransportModel(nn.Module):
    """Grammar trunk + relative interaction potential + gated residual transport."""

    def __init__(self, protein_dim: int, hidden_dim: int = 192,
                 task_dim: int = 48, ligand_layers: int = 4,
                 pair_dim: int = 96, pair_blocks: int = 4,
                 pair_latents: int = 24, pair_heads: int = 8,
                 pair_chunk_size: int = 8, support_hidden_dim: int = 192,
                 support_blocks: int = 2, adapter_rank: int = 4,
                 adaptive_blocks: int = 2, adapter_scale: float = 0.25,
                 anchors: int = 16, rank: int = 96,
                 similarity_scale: float = 8.0, use_cartesian: bool = False,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        del pair_blocks, pair_chunk_size, support_hidden_dim
        del support_blocks, adapter_rank, adaptive_blocks, adapter_scale
        if use_cartesian:
            raise ValueError(
                "no BindingDB deployment pair has a common-frame complex; "
                "see stage5 GEOMETRY_COVERAGE_AUDIT.json")
        self.hidden_dim = int(hidden_dim)
        self.embed_dim = int(pair_dim)
        self.rank = int(rank)
        self.anchor_count = int(anchors)
        self.contact_types = int(pair_latents)
        self.protein_encoder = ResidueEncoder(protein_dim, hidden_dim, dtype)
        self.ligand_encoder = LigandEncoder(hidden_dim, ligand_layers, dtype=dtype)
        self.grammar = ContactGrammar(
            hidden_dim, self.contact_types, pair_heads, self.embed_dim, dtype)
        self.embed = nn.Sequential(
            nn.Linear(4 * hidden_dim + self.contact_types,
                      2 * self.embed_dim, dtype=dtype), nn.GELU(),
            nn.Linear(2 * self.embed_dim, self.embed_dim, dtype=dtype))
        self.embed_norm = nn.LayerNorm(self.embed_dim, dtype=dtype)
        self.relative = RelativePotential(
            self.embed_dim, rank, hidden_dim, dtype)
        self.gate = SymmetricGate(
            self.embed_dim, rank, hidden_dim, dtype)
        # Learned ligand-side reference points in interaction-embedding space.
        self.anchor = nn.Parameter(
            torch.randn(self.anchor_count, self.embed_dim, dtype=dtype) * 0.02)
        self.ligand_head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, 1, dtype=dtype))
        self.level_head = nn.Sequential(
            nn.LayerNorm(2 * hidden_dim, dtype=dtype),
            nn.Linear(2 * hidden_dim, task_dim, dtype=dtype), nn.GELU(),
            nn.Linear(task_dim, 1, dtype=dtype))
        self.transport = RelativeTransport(
            self.embed_dim, similarity_scale, dtype)

    def target_level(self, residues: Tensor, summary: Tensor,
                     protein_mask: Tensor) -> Tensor:
        """Per-protein level with attention pooling over the residue slots.

        The R6a/b screenings both eliminated the routed design at 300 steps
        on k=0 calibration (A2 1.66 vs A0 1.29): a two-vector
        (summary, mean) input cannot carry the calibration the incumbent's
        whole trunk provides. Attention pooling gives the level readout the
        full residue set, with `summary` as the query — one changed module,
        no routing or shape changes.
        """
        scores = (summary.unsqueeze(1)
                  @ residues.transpose(-1, -2)) / self.hidden_dim ** 0.5
        scores = scores.masked_fill(
            ~protein_mask[:, None, :].bool(), torch.finfo(scores.dtype).min)
        context = torch.einsum("bsd,bsh->bh", torch.softmax(scores, -1),
                               residues)
        return self.level_head(torch.cat((summary, context), -1)).squeeze(-1)

    @staticmethod
    def atom_chemistry(raw_atoms: Tensor) -> Tensor:
        positive = raw_atoms[..., 21]
        negative = raw_atoms[..., 18] + raw_atoms[..., 19]
        aromatic = raw_atoms[..., 27]
        hydrophobic = raw_atoms[..., [1, 4, 7, 8, 9]].sum(-1).clamp_max(1)
        return torch.stack((positive, negative, aromatic, hydrophobic), -1)

    # ---------------- branches ---------------------------------------------
    def encode_protein(self, protein_pooled: Tensor, protein_tokens: Tensor,
                       protein_mask: Tensor, protein_chemistry: Tensor | None):
        return self.protein_encoder(
            protein_pooled, protein_tokens, protein_mask, protein_chemistry)

    def encode_ligand(self, raw_atoms: Tensor, bonds: Tensor,
                      mask: Tensor) -> tuple[Tensor, Tensor]:
        """[B,N,A,F] -> ligand mean [B,N,h] and atom states [B,N,A,h]."""
        batch, count = raw_atoms.shape[:2]
        ligand, atom_states = self.ligand_encoder(
            raw_atoms.flatten(0, 1), bonds.flatten(0, 1), mask.flatten(0, 1))
        return (ligand.reshape(batch, count, -1),
                atom_states.reshape(batch, count, *atom_states.shape[1:]))

    def interaction_embed(self, atom_states: Tensor, atom_mask: Tensor,
                          raw_atoms: Tensor, residues: Tensor,
                          residue_mask: Tensor,
                          summary: Tensor) -> Tensor:
        """Atom-to-residue grammar -> per-ligand interaction embedding [B,N,e]."""
        batch, count, atom_count, _ = atom_states.shape
        residue_count = residues.shape[1]
        wide_residues = residues[:, None].expand(-1, count, -1, -1).reshape(
            batch * count, residue_count, -1)
        wide_mask = residue_mask[:, None].expand(-1, count, -1).reshape(
            batch * count, residue_count)
        occupancy, mean_state, max_state = self.grammar(
            atom_states.flatten(0, 1), atom_mask.flatten(0, 1),
            self.atom_chemistry(raw_atoms).flatten(0, 1),
            wide_residues, wide_mask)
        flat_mask = atom_mask.flatten(0, 1)
        denominator = flat_mask.sum(1, keepdim=True).clamp_min(1.0)
        ligand_mean = ((atom_states.flatten(0, 1) * flat_mask.unsqueeze(-1))
                       .sum(1) / denominator)
        wide_summary = summary[:, None].expand(-1, count, -1).reshape(
            batch * count, -1)
        embed = self.embed_norm(self.embed(torch.cat(
            (ligand_mean, mean_state, max_state, wide_summary, occupancy), -1)))
        return embed.reshape(batch, count, -1)

    def anchor_shape(self, u: Tensor, embed: Tensor) -> tuple[Tensor, Tensor]:
        """shape(P,L) = mean_m delta(P, L, anchor_m); exact zero anchor mean."""
        anchors = self.anchor.unsqueeze(0).expand(embed.shape[0], -1, -1)
        delta = self.relative.delta_matrix(u, embed, anchors)
        # mean_m shape(P, anchor_m): the anchor-set mean of the anchor shapes,
        # exactly zero by antisymmetry (diagnostic, never subtracted).
        anchor_delta = self.relative.delta_matrix(u, anchors, anchors)
        anchor_mean = anchor_delta.mean(-1).mean(-1)
        return delta.mean(-1), anchor_mean

    def forward_parts(self, protein_pooled: Tensor, protein_tokens: Tensor,
                      protein_mask: Tensor, raw_atoms: Tensor, bonds: Tensor,
                      mask: Tensor, protein_chemistry: Tensor | None):
        """One pass over all ligands of one protein: every intermediate."""
        residues, summary = self.encode_protein(
            protein_pooled, protein_tokens, protein_mask, protein_chemistry)
        ligand, atom_states = self.encode_ligand(raw_atoms, bonds, mask)
        embed = self.interaction_embed(
            atom_states, mask, raw_atoms, residues, protein_mask, summary)
        gate = protein_mask.unsqueeze(-1)
        residue_mean = (residues * gate).sum(1) / gate.sum(1).clamp_min(1.0)
        u = self.relative.direction_vector(summary, residue_mean)
        u_gate = self.gate.direction_vector(summary, residue_mean)
        prior = self.ligand_head(ligand).squeeze(-1)
        level = self.target_level(residues, summary, protein_mask)
        shape, anchor_mean = self.anchor_shape(u, embed)
        level = level.unsqueeze(-1).expand_as(prior)
        endpoint = prior + level + shape
        return endpoint, prior, level, shape, anchor_mean, embed, u, u_gate

    # ---------------- full forward -----------------------------------------
    def forward(self, protein_pooled: Tensor, protein_tokens: Tensor,
                protein_mask: Tensor, support_atoms: Tensor, support_bonds: Tensor,
                support_mask: Tensor, support_y: Tensor, query_atoms: Tensor,
                query_bonds: Tensor, query_mask: Tensor, *, adapt: bool = True,
                protein_chemistry: Tensor | None = None,
                support_fingerprint: Tensor | None = None,
                query_fingerprint: Tensor | None = None,
                geometry_coordinates: Tensor | None = None,
                geometry_edge_index: Tensor | None = None,
                geometry_available: Tensor | None = None,
                geometry_common_frame: Tensor | None = None) -> RelTransportOutput:
        if geometry_coordinates is not None or geometry_edge_index is not None or (
                geometry_available is not None and bool(geometry_available.any())):
            raise ValueError(
                "no BindingDB deployment pair has a common-frame complex; "
                "see stage5 GEOMETRY_COVERAGE_AUDIT.json")
        parameter = next(self.parameters())
        device, dtype = parameter.device, parameter.dtype
        unbatched = protein_pooled.ndim == 1
        if unbatched:
            def add(value):
                return None if value is None else value.unsqueeze(0)
            protein_pooled, protein_tokens = add(protein_pooled), add(protein_tokens)
            protein_mask, protein_chemistry = add(protein_mask), add(protein_chemistry)
            support_atoms, support_bonds = add(support_atoms), add(support_bonds)
            support_mask, support_y = add(support_mask), add(support_y)
            query_atoms, query_bonds = add(query_atoms), add(query_bonds)
            query_mask = add(query_mask)
            support_fingerprint = add(support_fingerprint)
            query_fingerprint = add(query_fingerprint)
        support_count = support_atoms.shape[1]
        # Pad both groups to one common atom width before concatenation, so
        # the model accepts episodes whose support and query panels were
        # materialized with different padded widths.
        width = max(support_atoms.shape[-2], query_atoms.shape[-2], 1)

        def pad_atoms(values: Tensor, target: int) -> Tensor:
            return F.pad(values, (0, 0, 0, target - values.shape[-2]))

        def pad_bonds(values: Tensor, target: int) -> Tensor:
            return F.pad(values, (0, 0, 0, target - values.shape[-2],
                                  0, target - values.shape[-2]))

        def pad_mask(values: Tensor, target: int) -> Tensor:
            return F.pad(values, (0, target - values.shape[-1]))

        support_atoms, support_bonds, support_mask = (
            pad_atoms(support_atoms, width), pad_bonds(support_bonds, width),
            pad_mask(support_mask, width))
        query_atoms, query_bonds, query_mask = (
            pad_atoms(query_atoms, width), pad_bonds(query_bonds, width),
            pad_mask(query_mask, width))
        raw_atoms = torch.cat((support_atoms, query_atoms), 1).to(device, dtype)
        bonds = torch.cat((support_bonds, query_bonds), 1).to(device, dtype)
        mask = torch.cat((support_mask, query_mask), 1).to(device, dtype)
        endpoint, prior, level, shape, anchor_mean, embed, u, u_gate = \
            self.forward_parts(
                protein_pooled.to(device, dtype), protein_tokens.to(device, dtype),
                protein_mask.to(device, dtype), raw_atoms, bonds, mask,
                None if protein_chemistry is None
                else protein_chemistry.to(device, dtype))
        query_count = endpoint.shape[1] - support_count
        support_endpoint, query_endpoint = torch.split(
            endpoint, (support_count, query_count), 1)
        support_embed, query_embed = torch.split(
            embed, (support_count, query_count), 1)

        batch = endpoint.shape[0]
        if not adapt or support_count == 0:
            transport = torch.zeros_like(query_endpoint)
            residual = query_endpoint.new_zeros(batch, support_count)
            weight = query_endpoint.new_zeros(
                batch, query_count, max(support_count, 1))
            rho = weight
            similarity = weight
            delta = weight
        else:
            if support_fingerprint is None or query_fingerprint is None:
                raise ValueError("the residual transport requires fingerprints")
            similarity = tanimoto(query_fingerprint.to(device, dtype),
                                  support_fingerprint.to(device, dtype))
            residual = (support_y.to(device, dtype) - support_endpoint).detach()
            delta = self.relative.delta_matrix(u, query_embed, support_embed)
            rho = self.gate.gate_matrix(u_gate, query_embed, support_embed)
            shrink = self.transport.shrinkage(support_count, residual)
            transport, weight = self.transport(
                support_embed, query_embed, residual, similarity, rho)
            transport = shrink * transport

        output = RelTransportOutput(
            prediction=query_endpoint + transport, endpoint=query_endpoint,
            ligand_prior=prior[:, support_count:],
            target_level=level[:, support_count:],
            shape=shape[:, support_count:], transport=transport,
            support_endpoint=support_endpoint, support_residual=residual,
            weight=weight, rho=rho, similarity=similarity,
            delta=delta, anchor_shape_mean=anchor_mean)
        if not unbatched:
            return output
        return RelTransportOutput(*(value.squeeze(0)
            if isinstance(value, Tensor) and value.ndim else value
            for value in output.__dict__.values()))
