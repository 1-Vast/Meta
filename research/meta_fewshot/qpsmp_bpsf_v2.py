"""Research-only BPSF variants; none are authorized for the active model."""
from __future__ import annotations

import torch
import torch.nn as nn
from torch import Tensor

from model.bpsf import BipartitePairSectionFormer, PairSectionEncoding
from model.qpsmp_meta import QPSMPBioModel


class SharedPairLatent(nn.Module):
    """One pair-to-latent representation with light endpoint/section heads."""

    def __init__(self, independent: nn.Module) -> None:
        super().__init__()
        endpoint, section = independent.endpoint, independent.section
        if hasattr(independent, "interaction"):
            interaction = independent.interaction
            self.latents = interaction.latents
            self.value = interaction.value
            self.response_weight = interaction.response_weight
            self.blocks = interaction.blocks
            self.norm = interaction.norm
            self.endpoint = endpoint
            self.section = section
            return
        if endpoint.latents.shape != section.latents.shape:
            raise ValueError("endpoint and section latent shapes differ")
        self.latents = nn.Parameter(endpoint.latents.detach().clone())
        self.cross = endpoint.cross
        self.blocks = endpoint.blocks
        self.norm = endpoint.norm
        pair_dim = self.latents.shape[-1]
        self.endpoint = nn.Linear(
            pair_dim, endpoint.output.out_features, bias=False,
            dtype=endpoint.output.weight.dtype)
        self.section = nn.Linear(
            pair_dim, section.output.out_features, bias=False,
            dtype=section.output.weight.dtype)
        with torch.no_grad():
            self.endpoint.weight.copy_(endpoint.output.weight)
            self.section.weight.copy_(section.output.weight)

    def _representation(self, pair: Tensor, pair_mask: Tensor) -> Tensor:
        batch, _, _, width = pair.shape
        values = pair.reshape(batch, -1, width)
        valid = pair_mask.reshape(batch, -1).bool()
        latent = self.latents.unsqueeze(0).expand(batch, -1, -1)
        score = torch.einsum("bnd,md->bmn", values, self.latents) / width ** 0.5
        score = score.masked_fill(~valid[:, None], torch.finfo(score.dtype).min)
        attended = torch.einsum(
            "bmn,bnd->bmd", torch.softmax(score, -1), self.value(values))
        latent = latent + attended
        for block in self.blocks:
            latent = latent + block(latent)
        return self.norm(latent)

    def forward(self, pair: Tensor, pair_mask: Tensor, primitive_bias=None):
        del primitive_bias
        slots = self._representation(pair, pair_mask)
        representation = slots.mean(1)
        response = torch.einsum(
            "bmd,md->bm", slots, self.response_weight) \
            / self.response_weight.shape[-1] ** 0.5
        return self.endpoint(representation), self.section(representation), slots, response


class RelevanceWeightedBPSF(BipartitePairSectionFormer):
    """Masked pair relevance before latent pooling; no sparsity is claimed."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        pair_dim = self.atom_pair.out_features
        dtype = self.atom_pair.weight.dtype
        self.relevance = nn.Sequential(
            nn.LayerNorm(pair_dim, dtype=dtype),
            nn.Linear(pair_dim, pair_dim, dtype=dtype), nn.SiLU(),
            nn.Linear(pair_dim, 1, dtype=dtype))

    def _forward_chunk(self, atoms: Tensor, residues: Tensor, atom_mask: Tensor,
                       residue_mask: Tensor, adjacency: Tensor,
                       return_pair: bool,
                       atom_features: Tensor | None = None,
                       atom_chemistry: Tensor | None = None,
                       residue_chemistry: Tensor | None = None,
                       task_code: Tensor | None = None) -> PairSectionEncoding:
        atoms = atoms * atom_mask.unsqueeze(-1)
        residues = residues * residue_mask.unsqueeze(-1)
        adjacency = adjacency * atom_mask[:, :, None] * atom_mask[:, None, :]
        atom_pair = self.atom_pair(atoms)[:, :, None, :]
        residue_pair = self.residue_pair(residues)[:, None, :, :]
        pair = self.pair_init(torch.cat((
            atom_pair.expand(-1, -1, residues.shape[1], -1),
            residue_pair.expand(-1, atoms.shape[1], -1, -1),
            atom_pair * residue_pair), dim=-1))
        mask = atom_mask[:, :, None] * residue_mask[:, None, :]
        pair = pair * mask.unsqueeze(-1)
        for block in self.blocks:
            pair, atoms, residues = block(
                pair, atoms, residues, atom_mask, residue_mask, adjacency,
                task_code)
        gate = torch.sigmoid(self.relevance(pair)) * mask.unsqueeze(-1)
        weighted = pair * gate
        primitive_bias = None
        if atom_features is not None:
            atom_bias = self.atom_primitive(atom_features)[:, :, None]
            residue_bias = self.residue_primitive(residues)[:, None]
            primitive_bias = torch.tanh(atom_bias + residue_bias)
        if atom_chemistry is not None or residue_chemistry is not None:
            if atom_chemistry is None or residue_chemistry is None:
                raise ValueError("anchored primitives need atom and residue chemistry")
            positive_a, negative_a, aromatic_a, hydrophobic_a = atom_chemistry.unbind(-1)
            positive_r, negative_r, aromatic_r, hydrophobic_r = residue_chemistry.unbind(-1)
            anchored = torch.stack((
                positive_a[:, :, None] * negative_r[:, None, :],
                negative_a[:, :, None] * positive_r[:, None, :],
                aromatic_a[:, :, None] * aromatic_r[:, None, :],
                hydrophobic_a[:, :, None] * hydrophobic_r[:, None, :]), -1)
            dictionary = (self.latent.interaction.latents
                          if hasattr(self.latent, "interaction")
                          else self.latent.latents)
            missing = dictionary.shape[0] - anchored.shape[-1]
            if missing < 0:
                raise ValueError("primitive dictionary cannot hold anchored channels")
            anchored = torch.cat((anchored, anchored.new_zeros(
                *anchored.shape[:-1], missing)), -1)
            primitive_bias = anchored if primitive_bias is None else primitive_bias + anchored
        endpoint, section, slots, response = self.latent(
            weighted, mask, primitive_bias)
        return PairSectionEncoding(
            endpoint, section, slots, response, weighted if return_pair else None,
            mask if return_pair else None)


class ResearchQPSMP(QPSMPBioModel):
    """Factory-compatible research model with isolated representation switches."""

    shared_latent: bool = False
    relevance_weighted: bool = False

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        if self.relevance_weighted:
            base = self.pair_section
            replacement = RelevanceWeightedBPSF(
                token_dim=base.atom_pair.in_features,
                section_dim=base.latent.section.out_features,
                pair_dim=base.atom_pair.out_features,
                blocks=len(base.blocks),
                latent_count=base.latent.interaction.latents.shape[0],
                heads=base.latent.interaction.cross.num_heads,
                chunk_size=base.chunk_size,
                dtype=base.atom_pair.weight.dtype)
            replacement.load_state_dict(base.state_dict(), strict=False)
            self.pair_section = replacement
        if self.shared_latent:
            self.pair_section.latent = SharedPairLatent(self.pair_section.latent)

    def freeze_for_section_training(self) -> None:
        for parameter in self.parameters():
            parameter.requires_grad_(False)
        if self.shared_latent:
            for parameter in self.pair_section.latent.section.parameters():
                parameter.requires_grad_(True)
        else:
            for parameter in self.pair_section.latent.section.parameters():
                parameter.requires_grad_(True)
        for parameter in self.meta.section_operator.parameters():
            parameter.requires_grad_(True)


def research_model_factory(shared_latent: bool, relevance_weighted: bool):
    class Variant(ResearchQPSMP):
        pass
    Variant.shared_latent = shared_latent
    Variant.relevance_weighted = relevance_weighted
    Variant.__name__ = (
        f"ResearchQPSMP_shared{int(shared_latent)}_relevance{int(relevance_weighted)}")
    return Variant
