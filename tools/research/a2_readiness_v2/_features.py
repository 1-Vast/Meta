"""Frozen internal representations of the A0 trunk, extracted without training.

A2-min does not consume the endpoint scalar — it consumes `e0(P, L)`, an
internal representation. The v1 probe measured only the scalar, so it could not
distinguish "the trunk has no protein-conditioned ligand information" from "it
has some and the readout discards it". These are different diagnoses with
different remedies, and only the first would reject A2.

Six representations are extracted, in pipeline order:

| name | width | where |
|---|---|---|
| `occupancy` | `pair_latents` = 24 | `ContactGrammar`, contact-type occupancy |
| `mean_state` | `hidden_dim` = 192 | `ContactGrammar`, mean over atoms |
| `max_state` | `hidden_dim` = 192 | `ContactGrammar`, max over atoms |
| `embed` | `pair_dim` = 96 | **the `e0` of the A2 plan** |
| `section` | `task_dim` = 48 | normalised projection of `embed` |
| `interaction` | 1 | the scalar the v1 probe measured |

`ligand` (the protein-blind encoder output, 192) is extracted as the reference
for what a representation looks like when it provably carries no protein.
"""
from __future__ import annotations

import numpy as np
import torch

REPRESENTATIONS = ("occupancy", "mean_state", "max_state", "embed",
                   "section", "ligand", "interaction")

# Representations that are protein-conditioned by construction. `ligand` is not
# (it is the protein-blind encoder output), and it is kept as the negative
# reference rather than as a candidate.
PROTEIN_CONDITIONED = ("occupancy", "mean_state", "max_state", "embed",
                       "section", "interaction")


class Capture:
    """Forward hooks on the modules whose outputs `encode` does not return."""

    def __init__(self, model):
        self.values: dict[str, torch.Tensor] = {}
        self.handles = [
            model.grammar.register_forward_hook(self._grammar),
            model.section_norm.register_forward_hook(self._store("section")),
            model.embed_norm.register_forward_hook(self._store("embed")),
            model.ligand_encoder.register_forward_hook(self._ligand),
        ]

    def _grammar(self, module, inputs, output):
        occupancy, mean_state, max_state = output
        self.values["occupancy"] = occupancy.detach()
        self.values["mean_state"] = mean_state.detach()
        self.values["max_state"] = max_state.detach()

    def _store(self, name: str):
        def hook(module, inputs, output):
            self.values[name] = output.detach()
        return hook

    def _ligand(self, module, inputs, output):
        # LigandEncoder returns (pooled, atom_states); the pooled vector is the
        # protein-blind ligand representation.
        self.values["ligand"] = output[0].detach()

    def close(self) -> None:
        for handle in self.handles:
            handle.remove()


def extract(model, protein_parts: list, query_atoms, query_bonds, query_mask,
            query_fingerprint) -> dict[str, np.ndarray]:
    """Every representation for one target's query panel, one row per query.

    `adapt=False` with an empty support, so no label is read and the transport
    contributes nothing: this is the k=0 path exactly.
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
        capture.close()

    endpoint = output.zero_shot.squeeze(0).detach().float().cpu().numpy()
    ligand_value = output.ligand_only.squeeze(0).detach().float().cpu().numpy()
    additive = output.additive.squeeze(0).detach().float().cpu().numpy()

    features = {name: values[name].float().cpu().numpy()
                for name in ("occupancy", "mean_state", "max_state",
                             "section", "embed", "ligand")}
    features["interaction"] = (endpoint - additive)[:, None]
    features["_endpoint"] = endpoint[:, None]
    features["_protein_value"] = (additive - ligand_value)[:, None]
    features["_ligand_value"] = ligand_value[:, None]
    return features


def differential_geometry(correct: np.ndarray, wrong: np.ndarray) -> dict:
    """How much of a representation's *ligand-differential* the protein sets.

    The within-target mean is the level; removing it leaves exactly the part
    that can order ligands. Four numbers answering different questions:

    * `level_relative_shift` — does the protein move the representation at all?
      A positive control: if this is ~0 the substitution did not happen and no
      other number here is interpretable.
    * `differential_cosine` — is the *arrangement* of ligands relative to one
      another the same under both proteins? 1.0 means the protein changes
      nothing about how this representation orders ligands.
    * `differential_relative_shift` — the size of the change in that
      arrangement, as a fraction of the arrangement itself.
    * `differential_norm_ratio` — whether the wrong protein merely shrinks or
      inflates the arrangement rather than rotating it.
    """
    correct_level, wrong_level = correct.mean(0), wrong.mean(0)
    correct_diff = correct - correct_level
    wrong_diff = wrong - wrong_level
    correct_norm = float(np.linalg.norm(correct_diff))
    wrong_norm = float(np.linalg.norm(wrong_diff))
    level_norm = float(np.linalg.norm(correct_level))
    inner = float((correct_diff * wrong_diff).sum())
    return {
        "level_relative_shift": (
            float(np.linalg.norm(correct_level - wrong_level) / level_norm)
            if level_norm > 1e-12 else float("nan")),
        "differential_cosine": (
            inner / (correct_norm * wrong_norm)
            if correct_norm > 1e-12 and wrong_norm > 1e-12 else float("nan")),
        "differential_relative_shift": (
            float(np.linalg.norm(correct_diff - wrong_diff) / correct_norm)
            if correct_norm > 1e-12 else float("nan")),
        "differential_norm_ratio": (
            wrong_norm / correct_norm if correct_norm > 1e-12 else float("nan")),
    }
