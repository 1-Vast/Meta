"""Core Innovation A: a level-shape factorized cold-target predictor.

    f(P, L) = ligand_prior(L) + target_level(P) + centered_interaction(P, L)

Why this decomposition
----------------------
Stage 9 split cold-target k=0 error exactly into `(mean(p-y))^2 + var(p-y)` and
found 59% of it is target-level calibration, with the trained endpoint
contributing essentially nothing to within-target shape (concordance 0.525
against a 0.500 coin flip). Stage R0 sharpened that: on activity-cliff pairs —
within-target pairs with Tanimoto >= 0.6 and a >= 1.0 pK gap, exactly where
chemistry has to be read rather than recalled — the trained endpoint orders at
**chance** (0.519) while a parameter-free Morgan kernel reaches 0.716.

Four separate interventions (capacity, budget, a learned relative operator, a
sequence-locality prior) failed to move within-target shape. The common cause is
that a single scalar head trained on absolute affinity can buy almost all of its
squared error from the target level, and a level shift is constant across the
queries of a target, so it cannot change ordering. This module removes that
escape route architecturally.

The anchor-centering device
---------------------------
`centered_interaction` must not be able to relearn a constant target offset —
otherwise the factorization is cosmetic. It is therefore defined as

    centered_interaction(P, L) = s(P, L) - mean_m s(P, anchor_m)

where `anchor_1..M` are **learned ligand-side embeddings owned by the model**,
not episode data. The subtracted term depends only on the protein and on
parameters, so:

* `mean_m centered_interaction(P, anchor_m) = 0` exactly, for every protein —
  the branch has no constant component in the anchor basis;
* prediction for a query never depends on the other queries in the batch, so
  this is an ordinary inductive predictor. No query-panel mean, no transductive
  statistic. That distinction is what invalidated Stage 9's 25.5% figure and it
  is not repeated here.

Contracts preserved
-------------------
* one scalar endpoint; each affinity contribution has exactly one source;
* `adapt=False` or `k=0` returns exactly the zero-shot endpoint;
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
from .interaction_grammar import ResidueEncoder
from .similarity_grammar import tanimoto

# positive, negative, aromatic, hydrophobic, and an untyped global type
TYPES = 5
# Learned query slots per pharmacophore type. One mean-pooled vector per type
# proved too coarse: the incumbent trunk attends **per atom** to residue slots,
# and collapsing a ligand to five vectors before any protein contact cost most
# of the within-target ordering (Stage R3 arm A1). Several slots per type
# restore atom-level resolution while keeping every token in a space a learned
# anchor can also occupy, which is what the centering requires.
SLOTS_PER_TYPE = 3
TOKENS = TYPES * SLOTS_PER_TYPE + 1        # + the encoder's own pooled readout
CHANNELS = TOKENS                          # name kept for the public contract


@dataclass
class LevelShapeOutput:
    """Every field is a scalar-per-query tensor unless noted."""
    prediction: Tensor           # deployed output: endpoint + transport
    endpoint: Tensor             # zero-shot endpoint f(P, L)
    ligand_prior: Tensor
    target_level: Tensor         # constant across the queries of one protein
    centered: Tensor             # anchor-centered interaction
    transport: Tensor            # few-shot correction, exactly 0 at k=0
    support_endpoint: Tensor     # endpoint on the support ligands
    support_residual: Tensor     # y_k - endpoint(L_k), label locked
    weight: Tensor               # [B, Q, K] transport weights
    similarity: Tensor           # [B, Q, K] Tanimoto
    anchor_centering: Tensor     # [B] the subtracted per-protein constant


class TypedLigandChannels(nn.Module):
    """Pharmacophore-typed pooling of GINE atom states into `CHANNELS` vectors.

    Typed pooling is what makes an anchor meaningful: an anchor is a learned
    point in the same typed ligand space, so it can be pushed through the exact
    interaction module that real ligands use.

    The last channel is the ligand encoder's own pooled readout rather than
    another mask-weighted mean. That is deliberate: it keeps `LigandEncoder.out`
    and `LigandEncoder.norm` on the path from inputs to the prediction, so the
    model has no dead trainable branch.
    """

    def __init__(self, hidden_dim: int, dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.slot = nn.Parameter(
            torch.randn(TYPES, SLOTS_PER_TYPE, hidden_dim, dtype=dtype) * 0.02)
        self.key = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.value = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.project = nn.Linear(hidden_dim, hidden_dim, dtype=dtype)
        self.norm = nn.LayerNorm(hidden_dim, dtype=dtype)

    @staticmethod
    def chemistry(raw_atoms: Tensor) -> Tensor:
        positive = raw_atoms[..., 21]
        negative = raw_atoms[..., 18] + raw_atoms[..., 19]
        aromatic = raw_atoms[..., 27]
        hydrophobic = raw_atoms[..., [1, 4, 7, 8, 9]].sum(-1).clamp_max(1)
        return torch.stack((positive, negative, aromatic, hydrophobic), -1)

    def forward(self, atom_states: Tensor, atom_mask: Tensor, raw_atoms: Tensor,
                pooled: Tensor) -> Tensor:
        """[N,A,d], [N,A], [N,A,F], [N,d] -> [N, TOKENS, d].

        Each pharmacophore type owns `SLOTS_PER_TYPE` learned queries that
        attend over the atom states. The type membership enters as an additive
        log-bias rather than a hard mask, so a slot concentrates on atoms of its
        type without being unable to see the rest of the molecule.
        """
        typed = self.chemistry(raw_atoms) * atom_mask.unsqueeze(-1)
        weights = torch.cat((typed, atom_mask.unsqueeze(-1)), -1)      # [N,A,TYPES]
        hidden = atom_states.shape[-1]
        score = torch.einsum("tsd,nad->ntsa", self.slot,
                             self.key(atom_states)) / hidden ** 0.5
        bias = weights.clamp_min(1e-4).log().permute(0, 2, 1)[:, :, None, :]
        score = score + bias
        score = score.masked_fill(
            ~atom_mask.bool()[:, None, None, :], torch.finfo(score.dtype).min)
        tokens = torch.einsum("ntsa,nad->ntsd", torch.softmax(score, -1),
                              self.value(atom_states))
        stacked = torch.cat(
            (tokens.reshape(-1, TYPES * SLOTS_PER_TYPE, hidden),
             pooled.unsqueeze(1)), 1)
        return self.norm(self.project(stacked))


class TanimotoTransport(nn.Module):
    """Fixed Morgan/Tanimoto support kernel with two learned scalars.

    Retained unchanged from Stage 6/7, where it was the only mechanism in this
    project to improve squared error and ranking together on both splits. It is
    a **baseline** here, not a contribution: Innovation A and B are evaluated
    against it, and it is inactive at k=0.
    """

    def __init__(self, similarity_scale: float = 8.0,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        self.similarity_scale = nn.Parameter(
            torch.tensor(float(similarity_scale), dtype=dtype))
        self.log_shrinkage = nn.Parameter(torch.tensor(1.8546, dtype=dtype))

    def shrinkage(self, support_count: int, reference: Tensor) -> Tensor:
        strength = F.softplus(self.log_shrinkage)
        count = reference.new_tensor(float(support_count))
        return count / (count + strength)

    def forward(self, support_residual: Tensor,
                similarity: Tensor) -> tuple[Tensor, Tensor]:
        weight = torch.softmax(self.similarity_scale * similarity, -1)
        return torch.einsum("bqk,bk->bq", weight, support_residual), weight


class ChannelInteraction(nn.Module):
    """Typed ligand channels attend to protein residue slots; outputs a scalar.

    Shared across real ligands and learned anchors, which is what makes the
    anchor centering exact rather than approximate.

    The final `mix` layer and the readout carry **no bias**. Any additive
    constant introduced after the last nonlinearity appears identically in
    `s(P, L)` and in `mean_m s(P, anchor_m)` and therefore cancels in the
    centered output: such a bias is structurally unidentifiable and would sit in
    the optimizer as a dead trainable tensor. Removing it is required by the
    centering, not a tuning choice.
    """

    def __init__(self, hidden_dim: int, heads: int = 4,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        if hidden_dim % heads:
            raise ValueError("hidden_dim must be divisible by interaction heads")
        self.heads, self.head_dim = int(heads), hidden_dim // int(heads)
        self.query = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.key = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.value = nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype)
        self.channel = nn.Parameter(torch.zeros(TOKENS, hidden_dim, dtype=dtype))
        self.mix = nn.Sequential(
            nn.LayerNorm(3 * hidden_dim, dtype=dtype),
            nn.Linear(3 * hidden_dim, hidden_dim, dtype=dtype), nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim, bias=False, dtype=dtype))
        self.readout = nn.Linear(TOKENS * hidden_dim, 1, bias=False, dtype=dtype)
        nn.init.normal_(self.channel, std=0.02)

    def forward(self, channels: Tensor, residues: Tensor,
                residue_mask: Tensor) -> Tensor:
        """[B,N,C,d] x [B,S,d] x [B,S] -> [B,N] scalar interaction."""
        batch, count = channels.shape[:2]
        hidden = channels.shape[-1]
        flat = (channels + self.channel).reshape(batch, count * TOKENS, hidden)
        query = self.query(flat).reshape(
            batch, count * TOKENS, self.heads, self.head_dim).transpose(1, 2)
        key = self.key(residues).reshape(
            batch, -1, self.heads, self.head_dim).transpose(1, 2)
        value = self.value(residues).reshape(
            batch, -1, self.heads, self.head_dim).transpose(1, 2)
        score = query @ key.transpose(-1, -2) / self.head_dim ** 0.5
        score = score.masked_fill(
            ~residue_mask[:, None, None, :].bool(), torch.finfo(score.dtype).min)
        context = (torch.softmax(score, -1) @ value).transpose(1, 2).reshape(
            batch, count * TOKENS, hidden)
        state = self.mix(torch.cat((flat, context, flat * context), -1))
        return self.readout(state.reshape(batch, count, CHANNELS * hidden)).squeeze(-1)


class LevelShapeModel(nn.Module):
    """Level-shape factorized predictor with a fixed Tanimoto residual transport."""

    def __init__(self, protein_dim: int, hidden_dim: int = 192,
                 task_dim: int = 48, ligand_layers: int = 4,
                 pair_dim: int = 96, pair_blocks: int = 4,
                 pair_latents: int = 24, pair_heads: int = 4,
                 pair_chunk_size: int = 8, support_hidden_dim: int = 192,
                 support_blocks: int = 2, adapter_rank: int = 4,
                 adaptive_blocks: int = 2, adapter_scale: float = 0.25,
                 anchors: int = 32, similarity_scale: float = 8.0,
                 use_cartesian: bool = False,
                 dtype: torch.dtype = torch.float32) -> None:
        super().__init__()
        del pair_blocks, pair_latents, pair_chunk_size, support_hidden_dim
        del support_blocks, adapter_rank, adaptive_blocks, adapter_scale
        if use_cartesian:
            raise ValueError(
                "no BindingDB deployment pair has a common-frame complex; "
                "see stage5 GEOMETRY_COVERAGE_AUDIT.json")
        self.hidden_dim = int(hidden_dim)
        self.anchor_count = int(anchors)
        self.protein_encoder = ResidueEncoder(protein_dim, hidden_dim, dtype)
        self.ligand_encoder = LigandEncoder(hidden_dim, ligand_layers, dtype=dtype)
        self.channels = TypedLigandChannels(hidden_dim, dtype)
        self.interaction = ChannelInteraction(hidden_dim, pair_heads, dtype)
        # Learned ligand-side reference points. They are model parameters, never
        # episode data, so the centering they induce is inductive.
        self.anchor = nn.Parameter(
            torch.randn(self.anchor_count, TOKENS, hidden_dim, dtype=dtype) * 0.02)
        self.ligand_head = nn.Sequential(
            nn.LayerNorm(TOKENS * hidden_dim, dtype=dtype),
            nn.Linear(TOKENS * hidden_dim, pair_dim, dtype=dtype), nn.GELU(),
            nn.Linear(pair_dim, 1, dtype=dtype))
        self.level_head = nn.Sequential(
            nn.LayerNorm(2 * hidden_dim, dtype=dtype),
            nn.Linear(2 * hidden_dim, task_dim, dtype=dtype), nn.GELU(),
            nn.Linear(task_dim, 1, dtype=dtype))
        self.transport = TanimotoTransport(similarity_scale, dtype)

    # ---------------- branches ---------------------------------------------
    def encode_protein(self, protein_pooled: Tensor, protein_tokens: Tensor,
                       protein_mask: Tensor, protein_chemistry: Tensor | None):
        return self.protein_encoder(
            protein_pooled, protein_tokens, protein_mask, protein_chemistry)

    def encode_ligand(self, raw_atoms: Tensor, bonds: Tensor,
                      mask: Tensor) -> Tensor:
        """[B,N,A,F] -> typed channels [B,N,CHANNELS,d]."""
        batch, count = raw_atoms.shape[:2]
        pooled, atom_states = self.ligand_encoder(
            raw_atoms.flatten(0, 1), bonds.flatten(0, 1), mask.flatten(0, 1))
        channels = self.channels(
            atom_states, mask.flatten(0, 1), raw_atoms.flatten(0, 1), pooled)
        return channels.reshape(batch, count, TOKENS, self.hidden_dim)

    def ligand_prior(self, channels: Tensor) -> Tensor:
        return self.ligand_head(channels.flatten(-2)).squeeze(-1)

    def target_level(self, residues: Tensor, summary: Tensor,
                     protein_mask: Tensor) -> Tensor:
        gate = protein_mask.unsqueeze(-1)
        mean = (residues * gate).sum(1) / gate.sum(1).clamp_min(1.0)
        return self.level_head(torch.cat((summary, mean), -1)).squeeze(-1)

    def centered_interaction(self, channels: Tensor, residues: Tensor,
                             protein_mask: Tensor) -> tuple[Tensor, Tensor]:
        raw = self.interaction(channels, residues, protein_mask)
        anchors = self.anchor.unsqueeze(0).expand(channels.shape[0], -1, -1, -1)
        centering = self.interaction(anchors, residues, protein_mask).mean(-1)
        return raw - centering.unsqueeze(-1), centering

    # ---------------- endpoint ---------------------------------------------
    def endpoint_with_channels(self, channels: Tensor, protein_pooled: Tensor,
                               protein_tokens: Tensor, protein_mask: Tensor,
                               protein_chemistry: Tensor | None):
        """Endpoint for ligand channels that have already been encoded.

        `encode_ligand` is protein-blind, so a counterfactual protein only needs
        the protein encoder, the level head and the interaction re-run. That is
        a direct consequence of the factorization and makes the wrong-protein
        objective cheap enough to compute in every training step.
        """
        residues, summary = self.encode_protein(
            protein_pooled, protein_tokens, protein_mask, protein_chemistry)
        prior = self.ligand_prior(channels)
        level = self.target_level(residues, summary, protein_mask)
        centered, centering = self.centered_interaction(
            channels, residues, protein_mask)
        level = level.unsqueeze(-1).expand_as(prior)
        return prior + level + centered, prior, level, centered, centering

    def endpoint(self, protein_pooled: Tensor, protein_tokens: Tensor,
                 protein_mask: Tensor, raw_atoms: Tensor, bonds: Tensor,
                 mask: Tensor, protein_chemistry: Tensor | None):
        residues, summary = self.encode_protein(
            protein_pooled, protein_tokens, protein_mask, protein_chemistry)
        channels = self.encode_ligand(raw_atoms, bonds, mask)
        prior = self.ligand_prior(channels)
        level = self.target_level(residues, summary, protein_mask)
        centered, centering = self.centered_interaction(
            channels, residues, protein_mask)
        level = level.unsqueeze(-1).expand_as(prior)
        return prior + level + centered, prior, level, centered, centering

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
                geometry_common_frame: Tensor | None = None) -> LevelShapeOutput:
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
        raw_atoms = torch.cat((support_atoms, query_atoms), 1).to(device, dtype)
        bonds = torch.cat((support_bonds, query_bonds), 1).to(device, dtype)
        mask = torch.cat((support_mask, query_mask), 1).to(device, dtype)
        endpoint, prior, level, centered, centering = self.endpoint(
            protein_pooled.to(device, dtype), protein_tokens.to(device, dtype),
            protein_mask.to(device, dtype), raw_atoms, bonds, mask,
            None if protein_chemistry is None
            else protein_chemistry.to(device, dtype))
        query_count = endpoint.shape[1] - support_count
        support_endpoint, query_endpoint = torch.split(
            endpoint, (support_count, query_count), 1)

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
            transport, weight = self.transport(residual, similarity)
            transport = shrink * transport

        output = LevelShapeOutput(
            prediction=query_endpoint + transport, endpoint=query_endpoint,
            ligand_prior=prior[:, support_count:],
            target_level=level[:, support_count:],
            centered=centered[:, support_count:], transport=transport,
            support_endpoint=support_endpoint, support_residual=residual,
            weight=weight, similarity=similarity, anchor_centering=centering)
        if not unbatched:
            return output
        return LevelShapeOutput(*(value.squeeze(0)
            if isinstance(value, Tensor) and value.ndim else value
            for value in output.__dict__.values()))
