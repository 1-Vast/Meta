"""Shared arm construction for the A2-readiness probes.

A `randinit` arm is the same architecture with re-randomised weights. It
separates two diagnoses of a protein-inert interaction branch: an architecture
that cannot express protein-conditioned ligand ordering, versus training that
removed it. Only the second is a training-innovation target.
"""
from __future__ import annotations

from pathlib import Path

import torch
import torch.nn as nn

from scripts.stageR6_compare_arms import load_arm

RANDOM_INIT_SEED = 20260816


def build_arm(path: str | Path, tag: str, data, device: str):
    """Load a checkpoint; re-randomise every weight when `tag` says so."""
    model, kind, seed = load_arm(Path(path), data, device)
    if not tag.startswith("randinit"):
        return model, kind, seed
    generator = torch.Generator(device="cpu").manual_seed(RANDOM_INIT_SEED)

    def reinitialise(module: nn.Module) -> None:
        if isinstance(module, (nn.Linear, nn.Embedding)):
            weight = torch.empty(module.weight.shape, device="cpu",
                                 dtype=torch.float32)
            nn.init.xavier_uniform_(weight, generator=generator)
            with torch.no_grad():
                module.weight.copy_(weight.to(module.weight.device,
                                              module.weight.dtype))
                if getattr(module, "bias", None) is not None:
                    module.bias.zero_()

    model.apply(reinitialise)
    return model, kind, -1
