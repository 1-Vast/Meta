"""Arm construction for Phase 1: trained A0 seeds and independent random inits.

The v1 probe used a single random-initialisation arm and reported that its
protein sensitivity was 110x the trained model's. One draw cannot distinguish
a structural property of the architecture from the noise of one initialisation,
so this version builds ten, each with its own seed, and the analysis measures
across-seed agreement rather than magnitude alone.
"""
from __future__ import annotations

from pathlib import Path

import torch
import torch.nn as nn

from scripts.stageR6_compare_arms import load_arm


def randomise(model: nn.Module, seed: int) -> nn.Module:
    """Re-initialise every Linear/Embedding weight from a seeded generator.

    Xavier-uniform with zeroed biases, matching the v1 procedure exactly so
    that the two cycles' random arms are comparable. LayerNorm affine
    parameters are left at their loaded values, which for a trained checkpoint
    means the normalisation scales are *not* random — recorded here because it
    makes this arm a conservative test of "the architecture can express it":
    a fully random model would be further from the trained one, not closer.
    """
    generator = torch.Generator(device="cpu").manual_seed(int(seed))

    def apply(module: nn.Module) -> None:
        if isinstance(module, (nn.Linear, nn.Embedding)):
            weight = torch.empty(module.weight.shape, device="cpu",
                                 dtype=torch.float32)
            nn.init.xavier_uniform_(weight, generator=generator)
            with torch.no_grad():
                module.weight.copy_(weight.to(module.weight.device,
                                              module.weight.dtype))
                if getattr(module, "bias", None) is not None:
                    module.bias.zero_()

    model.apply(apply)
    return model


def trained_arm(path: Path, data, device: str):
    model, kind, seed = load_arm(Path(path), data, device)
    return model, kind, int(seed)


def random_arm(path: Path, data, device: str, seed: int):
    """Same architecture and shapes as the checkpoint, weights re-randomised."""
    model, kind, _ = load_arm(Path(path), data, device)
    return randomise(model, seed), kind, int(seed)
