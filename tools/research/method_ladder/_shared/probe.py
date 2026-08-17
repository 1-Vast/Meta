"""The low-cost discriminator: a tiny head trained by ordinary SGD.

Deliberately *not* a ridge fit. The cycle's constraints prohibit closed-form
solvers, and the prohibition applies to diagnostics as well as to candidates —
a probe that is fitted analytically measures the information a linear solve can
extract, which is not the quantity a gradient-trained model will actually get.

The probe predicts the **within-panel centered** label. Centering removes the
per-target level exactly, so the probe cannot score by knowing a target's
average affinity: it is measured only on how it orders that target's ligands.

Feature standardization always uses `meta_train` statistics. Selection (the
weight decay, the only knob) always uses `meta_train` component folds.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
import torch.nn as nn

from ._frozen import (
    FOLD_SEED, PROBE_BATCH_TARGETS, PROBE_FOLDS, PROBE_LR, PROBE_SEEDS,
    PROBE_STEPS, PROBE_WEIGHT_DECAYS,
)
from .panels import within_target_r


@dataclass(frozen=True)
class Standardizer:
    mean: np.ndarray
    scale: np.ndarray

    def apply(self, values: np.ndarray) -> np.ndarray:
        return (values - self.mean) / self.scale


def fit_standardizer(blocks: list[np.ndarray]) -> Standardizer:
    stacked = np.concatenate(blocks, axis=0).astype(np.float64)
    mean = stacked.mean(0)
    scale = stacked.std(0)
    scale[scale < 1e-8] = 1.0
    return Standardizer(mean=mean, scale=scale)


class LinearProbe(nn.Module):
    """`y_hat = <w, x>`. One linear functional, stated without a fake rank.

    The zero init makes the fit reproducible and removes initialisation as a
    source of variance, but it also means the probe seed enters only through
    `train_probe`'s minibatch sampler: with at most `PROBE_BATCH_TARGETS`
    panels the run is deterministic full-batch descent and the seed does
    nothing. The ladder's `meta_train` has 346 panels against a batch of 32, so
    seeds do separate runs there — but a small-panel diagnostic must not be
    reported as if three seeds gave three independent fits.
    """

    def __init__(self, width: int):
        super().__init__()
        self.weight = nn.Linear(width, 1, bias=False)
        nn.init.zeros_(self.weight.weight)

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return self.weight(values).squeeze(-1)


class MLPProbe(nn.Module):
    """The capacity-matched control: same inputs, strictly more capacity.

    Its role is to show that a null is not an artifact of linearity. If a
    representation carries ordering information that only a nonlinear readout
    can reach, this probe finds it and the linear one does not.
    """

    def __init__(self, width: int, hidden: int = 64):
        super().__init__()
        self.body = nn.Sequential(
            nn.Linear(width, hidden), nn.GELU(), nn.Linear(hidden, 1))

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return self.body(values).squeeze(-1)


def centered_batch_loss(prediction: torch.Tensor, truth: torch.Tensor
                        ) -> torch.Tensor:
    """Squared error on the centered panel. One panel, one term."""
    if prediction.numel() < 2:
        return prediction.new_zeros(())
    p = prediction - prediction.mean()
    t = truth - truth.mean()
    return (p - t).square().mean()


def train_probe(train_blocks: list[tuple[np.ndarray, np.ndarray]],
                width: int, weight_decay: float, seed: int,
                device: str = "cpu", kind: str = "linear") -> nn.Module:
    """Fit on a list of `(features, labels)` panels. Ordinary AdamW, no solve."""
    torch.manual_seed(seed)
    probe = (LinearProbe(width) if kind == "linear"
             else MLPProbe(width)).to(device).double()
    optimizer = torch.optim.AdamW(probe.parameters(), lr=PROBE_LR,
                                  weight_decay=weight_decay)
    tensors = [(torch.as_tensor(x, dtype=torch.float64, device=device),
                torch.as_tensor(y, dtype=torch.float64, device=device))
               for x, y in train_blocks]
    rng = np.random.default_rng(seed)
    for _ in range(PROBE_STEPS):
        picks = rng.choice(len(tensors),
                           size=min(PROBE_BATCH_TARGETS, len(tensors)),
                           replace=False)
        loss = torch.zeros((), dtype=torch.float64, device=device)
        for index in picks:
            features, labels = tensors[int(index)]
            loss = loss + centered_batch_loss(probe(features), labels)
        loss = loss / len(picks)
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
    probe.eval()
    return probe


def evaluate_probe(probe: nn.Module, blocks: list[tuple[np.ndarray, np.ndarray]],
                   device: str = "cpu") -> list[float]:
    """Within-target `r` for each panel."""
    out = []
    with torch.no_grad():
        for features, labels in blocks:
            prediction = probe(torch.as_tensor(
                features, dtype=torch.float64, device=device)).cpu().numpy()
            out.append(within_target_r(prediction, labels))
    return out


def component_folds(components: list[str], folds: int = PROBE_FOLDS,
                    seed: int = FOLD_SEED) -> list[list[str]]:
    """Partition distinct components into `folds` groups.

    Grouping by component, never by target, is what makes the held-out fold a
    genuine cold-protein estimate: two targets of one homology component are
    not independent evidence about an unseen component.
    """
    unique = sorted(set(components))
    order = np.random.default_rng(seed).permutation(len(unique))
    buckets: list[list[str]] = [[] for _ in range(folds)]
    for rank, index in enumerate(order):
        buckets[rank % folds].append(unique[int(index)])
    return buckets


def select_weight_decay(blocks: list[tuple[np.ndarray, np.ndarray]],
                        components: list[str], width: int, device: str = "cpu",
                        kind: str = "linear") -> tuple[float, dict]:
    """Choose the single knob on `meta_train` component folds only.

    Returns the chosen decay and the full fold table, so the selection is
    auditable rather than asserted.
    """
    buckets = component_folds(components)
    table: dict[str, float] = {}
    for decay in PROBE_WEIGHT_DECAYS:
        scores: list[float] = []
        for held_out in buckets:
            held = set(held_out)
            train = [b for b, c in zip(blocks, components) if c not in held]
            test = [b for b, c in zip(blocks, components) if c in held]
            if not train or not test:
                continue
            for seed in PROBE_SEEDS[:1]:      # one seed for selection, three for the report
                probe = train_probe(train, width, decay, seed, device, kind)
                scores.extend(evaluate_probe(probe, test, device))
        table[f"{decay:g}"] = float(np.mean(scores)) if scores else float("nan")
    best = max(table, key=lambda key: table[key])
    return float(best), {"fold_mean_r": table, "selected": best,
                         "folds": len(buckets), "kind": kind}
