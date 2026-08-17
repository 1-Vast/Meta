"""Frozen-representation capture, including the tensors `encode` never returns.

The families that propose learned pooling (4), basis reallocation (1) and
bilinear fusion (5) all make claims about tensors *inside* `ContactGrammar`.
`tools/research/a2_readiness_v2/_features.py` captures the pooled outputs, which
is one step too late: it cannot distinguish "the cross attention never produced
useful protein-ligand information" from "it did and the pooling discarded it".
That distinction is the whole content of the fusion hypothesis, so it needs the
pre-pooling tensors.

Two hooks recover them without editing the model:

* a pre-hook on `grammar` captures `atom_mask`, its second positional argument;
* a hook on `grammar.atom_context` captures its input and output. The input is
  `cat((atoms, context, atoms * context), -1)`, so slicing the middle third
  recovers `context` — the raw attention read-out over residue slots, before
  any fusion MLP touches it. The output is `state`, after fusion but before the
  mask and before pooling.

Nothing here trains, and nothing here modifies the model.
"""
from __future__ import annotations

import numpy as np
import torch

# Per-atom tensors, captured pre-pooling, each reduced by the poolings below.
ATOM_LEVEL = ("context", "state")
POOLINGS = ("mean", "max", "rms")

# Pooled tensors the model itself computes.
POOLED = ("occupancy", "mean_state", "max_state", "embed", "section", "ligand")


def pool(values: torch.Tensor, mask: torch.Tensor, how: str) -> torch.Tensor:
    """Reduce `[pairs, atoms, width]` over the masked atom axis."""
    gate = mask.unsqueeze(-1).to(values.dtype)
    denominator = gate.sum(1).clamp_min(1.0)
    if how == "mean":
        return (values * gate).sum(1) / denominator
    if how == "rms":
        return ((values * gate).square().sum(1) / denominator).sqrt()
    if how == "max":
        filled = values.masked_fill(gate == 0, torch.finfo(values.dtype).min)
        return torch.nan_to_num(filled.amax(1), neginf=0.0)
    raise ValueError(f"unknown pooling {how!r}")


class Capture:
    """Forward hooks that expose the pre-fusion and pre-pooling tensors."""

    def __init__(self, model):
        self.values: dict[str, torch.Tensor] = {}
        self.atom_mask: torch.Tensor | None = None
        grammar = model.grammar
        self.hidden = grammar.atom_query.out_features
        self.handles = [
            grammar.register_forward_pre_hook(self._grammar_inputs),
            grammar.register_forward_hook(self._grammar_output),
            grammar.atom_context.register_forward_hook(self._fusion),
            model.section_norm.register_forward_hook(self._store("section")),
            model.embed_norm.register_forward_hook(self._store("embed")),
            model.ligand_encoder.register_forward_hook(self._ligand),
        ]

    # `ContactGrammar.forward(atoms, atom_mask, atom_chemistry, residues, ...)`
    def _grammar_inputs(self, module, inputs):
        self.atom_mask = inputs[1].detach()

    def _grammar_output(self, module, inputs, output):
        occupancy, mean_state, max_state = output
        self.values["occupancy"] = occupancy.detach()
        self.values["mean_state"] = mean_state.detach()
        self.values["max_state"] = max_state.detach()

    def _fusion(self, module, inputs, output):
        fused = inputs[0].detach()
        # cat((atoms, context, atoms * context), -1) — the middle third is the
        # raw attention read-out over residue slots.
        self.values["context"] = fused[..., self.hidden:2 * self.hidden]
        self.values["state"] = output.detach()

    def _store(self, name: str):
        def hook(module, inputs, output):
            self.values[name] = output.detach()
        return hook

    def _ligand(self, module, inputs, output):
        self.values["ligand"] = output[0].detach()

    def close(self) -> None:
        for handle in self.handles:
            handle.remove()


def extract(model, protein_parts, query_atoms, query_bonds, query_mask,
            query_fingerprint) -> dict[str, np.ndarray]:
    """Every representation for one target's query panel, one row per ligand.

    `adapt=False` with an empty support: the k=0 path exactly. No label is read
    and the transport contributes nothing.
    """
    pooled, tokens, mask, chemistry = protein_parts
    capture = Capture(model)
    try:
        with torch.no_grad():
            output = model(
                pooled, tokens, mask,
                query_atoms[:, :0], query_bonds[:, :0], query_mask[:, :0],
                torch.zeros(1, 0, device=pooled.device, dtype=pooled.dtype),
                query_atoms, query_bonds, query_mask,
                adapt=False, protein_chemistry=chemistry,
                support_fingerprint=query_fingerprint[:, :0],
                query_fingerprint=query_fingerprint)
    finally:
        values = dict(capture.values)
        atom_mask = capture.atom_mask
        capture.close()

    features: dict[str, np.ndarray] = {}
    for name in POOLED:
        features[name] = values[name].float().cpu().numpy()
    for name in ATOM_LEVEL:
        for how in POOLINGS:
            features[f"{name}_{how}"] = (
                pool(values[name], atom_mask, how).float().cpu().numpy())

    endpoint = output.zero_shot.squeeze(0).float().cpu().numpy()
    ligand_value = output.ligand_only.squeeze(0).float().cpu().numpy()
    additive = output.additive.squeeze(0).float().cpu().numpy()
    features["interaction"] = (endpoint - additive)[:, None]
    features["_endpoint"] = endpoint[:, None]
    features["_ligand_value"] = ligand_value[:, None]
    features["_protein_value"] = (additive - ligand_value)[:, None]
    return features


# The representations a family may name as its information source, in pipeline
# order. `ligand` is the protein-blind reference, not a candidate.
REPRESENTATIONS = (
    "ligand",                                        # protein-blind reference
    "context_mean", "context_max", "context_rms",    # pre-fusion attention
    "state_mean", "state_max", "state_rms",          # post-fusion, pre-pool
    "occupancy", "mean_state", "max_state",          # the model's own pooling
    "embed", "section", "interaction",               # readout path
)
