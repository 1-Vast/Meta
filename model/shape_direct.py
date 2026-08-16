"""Direct interaction-head shape with anchor centering and direct
difference supervision (R13 candidate, a fresh model family).

The R9-R12 ladder established the failure precisely: the relative-transport
family supervises a bilinear potential ``delta(P,i,j) = u^T[psi(e_i) phi(e_j)
- psi(e_j) phi(e_i)]`` while the deployed zero-shot ordering uses a DIFFERENT
quantity — ``shape(L_i) - shape(L_j)`` with ``shape(L) = mean_m delta(L,a_m)``
— because anchor means ``phi_bar``/``psi_bar`` replace the partner's
``phi(e_j)``/``psi(e_i)``. The supervised quantity and the deployed quantity
coincide only in the degenerate case. This module removes the leak: the
shape function itself is a direct readout ``s(e_L)`` of the interaction
embedding, anchor-centered for the no-constant guarantee, and the relative
supervision targets ``s(e_i) - s(e_j)`` — the exact quantity that decides
within-target ordering.

    f0(P, L)   = ligand_prior(L) + target_level(P) + shape(P, L)
    shape(P,L) = s(e(P,L)) - mean_m s(anchor_m)      (anchors: learned
                                                      interaction embeddings)
    t(q)       = shrink(n) * sum_k a(q,k) * r_k      (retained Tanimoto+key
                                                      baseline, no gate)

Contracts preserved: k=0 returns the endpoint exactly; support labels enter
only as residuals; support permutation invariance; query equivariance; query
labels never an input; no coordinate path; level-only abstention inside the
transport class. The closed gate family is not revived: there is no
query-specific transport mechanism here.
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
class ShapeDirectOutput:
    prediction: Tensor           # endpoint + transport
    endpoint: Tensor             # zero-shot endpoint
    ligand_prior: Tensor
    target_level: Tensor         # constant across the queries of one protein
    shape: Tensor                # anchor-centered direct interaction shape
    transport: Tensor            # few-shot correction, 0 at k=0
    support_endpoint: Tensor
    support_residual: Tensor     # label locked
    weight: Tensor               # [B, Q, K]
    similarity: Tensor           # [B, Q, K] Tanimoto
    anchor_centering: Tensor     # [B] mean_m s(anchor_m); shape's anchor-mean
                                 # is exactly zero by construction


class TanimotoTransport(nn.Module):
    """Retained Stage 6/7 Tanimoto+key residual weighting (baseline)."""

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
                support_residual: Tensor, similarity: Tensor
                ) -> tuple[Tensor, Tensor]:
        logits = self.similarity_scale * similarity
        if support_embed.shape[1] > 1:
            support_key = F.normalize(self.key(support_embed), dim=-1)
            query_key = F.normalize(self.key(query_embed), dim=-1)
            temperature = F.softplus(self.log_temperature) + 1.0
            logits = logits + temperature * torch.einsum(
                "bqh,bkh->bqk", query_key, support_key)
        weight = torch.softmax(logits, -1)
        transport = torch.einsum("bqk,bk->bq", weight, support_residual)
        return transport, weight


class ShapeHead(nn.Module):
    """Protein-conditioned direct interaction readout s([e; u(P)]) with
    learned anchor centering. The anchor set is centered exactly, and the
    supervised difference `s(e_i;P) - s(e_j;P)` is the deployed ordering
    quantity — the difference supervision is direct by construction."""

    def __init__(self, embed_dim: int, hidden_dim: int, protein_dim: int,
                 anchors: int, dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.protein = nn.Sequential(
            nn.LayerNorm(protein_dim, dtype=dtype),
            nn.Linear(protein_dim, hidden_dim, dtype=dtype))
        self.s = nn.Sequential(
            nn.LayerNorm(embed_dim + hidden_dim, dtype=dtype),
            nn.Linear(embed_dim + hidden_dim, hidden_dim, dtype=dtype),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, 1, bias=False, dtype=dtype))
        self.anchor = nn.Parameter(
            torch.randn(anchors, embed_dim, dtype=dtype) * 0.02)

    def forward(self, embed: Tensor, protein_vector: Tensor
                ) -> tuple[Tensor, Tensor]:
        """shape(L) = s([e_L; u(P)]) - mean_m s([anchor_m; u(P)])."""
        u = self.protein(protein_vector)                       # [B, h]
        ligand_in = torch.cat(
            (embed, u[:, None, :].expand(-1, embed.shape[1], -1)), -1)
        anchor_in = torch.cat(
            (self.anchor.unsqueeze(0).expand(embed.shape[0], -1, -1),
             u[:, None, :].expand(-1, self.anchor.shape[0], -1)), -1)
        values = self.s(ligand_in).squeeze(-1)
        anchor_values = self.s(anchor_in).squeeze(-1)
        anchor_mean = anchor_values.mean(1)
        return values - anchor_mean, anchor_mean


class ShapeDirectModel(nn.Module):
    """Grammar trunk + direct interaction-head shape + Tanimoto transport."""

    def __init__(self, protein_dim: int, hidden_dim: int = 192,
                 task_dim: int = 48, ligand_layers: int = 4,
                 pair_dim: int = 96, pair_blocks: int = 4,
                 pair_latents: int = 24, pair_heads: int = 8,
                 pair_chunk_size: int = 8, support_hidden_dim: int = 192,
                 support_blocks: int = 2, adapter_rank: int = 4,
                 adaptive_blocks: int = 2, adapter_scale: float = 0.25,
                 anchors: int = 16, shape_hidden: int = 96,
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
        self.shape_head = ShapeHead(self.embed_dim, int(shape_hidden),
                                    2 * hidden_dim, int(anchors), dtype)
        self.ligand_head = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, 1, dtype=dtype))
        self.level_head = nn.Sequential(
            nn.LayerNorm(2 * hidden_dim, dtype=dtype),
            nn.Linear(2 * hidden_dim, task_dim, dtype=dtype), nn.GELU(),
            nn.Linear(task_dim, 1, dtype=dtype))
        self.transport = TanimotoTransport(
            self.embed_dim, similarity_scale, dtype)

    @staticmethod
    def atom_chemistry(raw_atoms: Tensor) -> Tensor:
        positive = raw_atoms[..., 21]
        negative = raw_atoms[..., 18] + raw_atoms[..., 19]
        aromatic = raw_atoms[..., 27]
        hydrophobic = raw_atoms[..., [1, 4, 7, 8, 9]].sum(-1).clamp_max(1)
        return torch.stack((positive, negative, aromatic, hydrophobic), -1)

    def encode_protein(self, protein_pooled: Tensor, protein_tokens: Tensor,
                       protein_mask: Tensor, protein_chemistry: Tensor | None):
        return self.protein_encoder(
            protein_pooled, protein_tokens, protein_mask, protein_chemistry)

    def target_level(self, residues: Tensor, summary: Tensor,
                     protein_mask: Tensor) -> Tensor:
        scores = (summary.unsqueeze(1)
                  @ residues.transpose(-1, -2)) / self.hidden_dim ** 0.5
        scores = scores.masked_fill(
            ~protein_mask[:, None, :].bool(), torch.finfo(scores.dtype).min)
        context = torch.einsum("bsd,bsh->bh", torch.softmax(scores, -1),
                               residues)
        return self.level_head(torch.cat((summary, context), -1)).squeeze(-1)

    def forward_parts(self, protein_pooled: Tensor, protein_tokens: Tensor,
                      protein_mask: Tensor, raw_atoms: Tensor, bonds: Tensor,
                      mask: Tensor, protein_chemistry: Tensor | None):
        residues, summary = self.encode_protein(
            protein_pooled, protein_tokens, protein_mask, protein_chemistry)
        batch, count, atom_count, _ = raw_atoms.shape
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
        flat_mask = mask.flatten(0, 1)
        denominator = flat_mask.sum(1, keepdim=True).clamp_min(1.0)
        ligand_mean = ((atom_states * flat_mask.unsqueeze(-1)).sum(1)
                       / denominator)
        wide_summary = summary[:, None].expand(-1, count, -1).reshape(
            batch * count, -1)
        embed = self.embed_norm(self.embed(torch.cat(
            (ligand_mean, mean_state, max_state, wide_summary, occupancy), -1)))
        embed = embed.reshape(batch, count, -1)
        ligand = ligand.reshape(batch, count, -1)
        gate = protein_mask.unsqueeze(-1)
        residue_mean = (residues * gate).sum(1) / gate.sum(1).clamp_min(1.0)
        protein_vector = torch.cat((summary, residue_mean), -1)
        shape, anchor_mean = self.shape_head(embed, protein_vector)
        prior = self.ligand_head(ligand).squeeze(-1)
        level = self.target_level(residues, summary, protein_mask)
        level = level.unsqueeze(-1).expand_as(prior)
        endpoint = prior + level + shape
        return endpoint, prior, level, shape, anchor_mean, embed, residue_mean

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
                geometry_common_frame: Tensor | None = None) -> ShapeDirectOutput:
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
        endpoint, prior, level, shape, anchor_mean, embed, _ = self.forward_parts(
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
            similarity = weight
        else:
            if support_fingerprint is None or query_fingerprint is None:
                raise ValueError("the residual transport requires fingerprints")
            similarity = tanimoto(query_fingerprint.to(device, dtype),
                                  support_fingerprint.to(device, dtype))
            residual = (support_y.to(device, dtype) - support_endpoint).detach()
            shrink = self.transport.shrinkage(support_count, residual)
            transport, weight = self.transport(
                support_embed, query_embed, residual, similarity)
            transport = shrink * transport
        output = ShapeDirectOutput(
            prediction=query_endpoint + transport, endpoint=query_endpoint,
            ligand_prior=prior[:, support_count:],
            target_level=level[:, support_count:],
            shape=shape[:, support_count:], transport=transport,
            support_endpoint=support_endpoint, support_residual=residual,
            weight=weight, similarity=similarity,
            anchor_centering=anchor_mean)
        if not unbatched:
            return output
        return ShapeDirectOutput(*(value.squeeze(0)
            if isinstance(value, Tensor) and value.ndim else value
            for value in output.__dict__.values()))
