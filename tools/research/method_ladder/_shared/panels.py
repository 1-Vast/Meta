"""Deterministic per-target ligand panels and the protein conditions.

One panel per target: unique ligands in a `stable_seed` order, capped at
`PANEL_MAX`, kept only when at least `PANEL_MIN` ligands survive. The cap keeps
large targets from dominating an equal-target statistic; the floor keeps a
within-target Pearson `r` from being reported on a handful of points.

Every panel is built from `QPSMPData.cells` alone, so it inherits the
fail-closed `meta_test` seal: a sealed cell is not in `cells` and therefore
cannot enter a panel.

Panels are label-blind in construction. Labels are attached afterwards and are
used only as the probe's target and as the reported truth.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch

from scripts.qpsmp_data import EpisodeSpec, stable_seed

from ._frozen import PANEL_MAX, PANEL_MIN, PANEL_SEED


@dataclass(frozen=True)
class Panel:
    split: str
    component: str
    target: str
    cells: tuple[int, ...]
    ligands: tuple[str, ...]
    labels: np.ndarray          # pK, one per ligand
    donor: str                  # matched-wrong protein, different component


def build_panels(data, split: str, donors: dict[str, dict[str, tuple[str, float]]],
                 stratum: str = "nearest") -> tuple[Panel, ...]:
    """Every eligible target of `split`, as a deterministic ligand panel.

    The donor is the *nearest* cross-component protein by default: the hardest
    wrong protein, and therefore the one a specificity claim must beat. A donor
    that is trivially far away would make any substitution look decisive.
    """
    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    panels: list[Panel] = []
    for target in sorted(data.tasks[split]):
        indices = data.tasks[split][target]
        rng = np.random.default_rng(stable_seed("ladder-panel", PANEL_SEED,
                                                split, target))
        ordered = data._unique_ligand_order(indices, rng)
        if len(ordered) < PANEL_MIN:
            continue
        chosen = tuple(int(index) for index in ordered[:PANEL_MAX])
        panels.append(Panel(
            split=split,
            component=component_of[target],
            target=target,
            cells=chosen,
            ligands=tuple(data.cells[i]["ligand_id"] for i in chosen),
            labels=np.asarray([data.cells[i]["pK"] for i in chosen],
                              dtype=np.float64),
            donor=donors[target][stratum][0]))
    return tuple(panels)


def panel_inputs(data, panel: Panel, device: str, dtype=torch.float32):
    """Batched ligand tensors for one panel, with a leading batch axis of 1."""
    spec = EpisodeSpec(panel.split, panel.component, panel.target,
                       (), panel.cells, panel.donor)
    batch = data.materialize(spec)
    return (batch.query_atoms.unsqueeze(0).to(device, dtype),
            batch.query_bonds.unsqueeze(0).to(device, dtype),
            batch.query_mask.unsqueeze(0).to(device, dtype),
            batch.query_fingerprint.unsqueeze(0).to(device, dtype))


def protein_parts(data, target: str, device: str, dtype=torch.float32):
    """The four protein tensors for `target`, batched, ready for `extract`."""
    pooled, tokens, mask = data.protein_for_target(target)
    chemistry = data.protein_chemistry_for_target(target)
    return (pooled.unsqueeze(0).to(device, dtype),
            tokens.unsqueeze(0).to(device, dtype),
            mask.unsqueeze(0).to(device, dtype),
            chemistry.unsqueeze(0).to(device, dtype))


def permuted_protein_assignment(panels: tuple[Panel, ...], seed: int
                                ) -> dict[str, str]:
    """A derangement of target -> protein across panels.

    Distinct from the matched-wrong donor: the donor is *chosen* to be the
    nearest legal protein, while this control destroys the target-protein
    correspondence globally. A representation that survives the donor but not
    the permutation is responding to protein similarity structure rather than
    to the identity of the recipient.
    """
    targets = [panel.target for panel in panels]
    rng = np.random.default_rng(seed)
    for _ in range(64):
        shuffled = list(rng.permutation(targets))
        if all(a != b for a, b in zip(targets, shuffled)):
            return dict(zip(targets, shuffled))
    # A derangement always exists for n >= 2; the loop only guards the sampler.
    rotated = targets[1:] + targets[:1]
    return dict(zip(targets, rotated))


def centered(values: np.ndarray) -> np.ndarray:
    return values - values.mean()


def within_target_r(prediction: np.ndarray, truth: np.ndarray) -> float:
    """Pearson `r` after removing each panel's own mean.

    Centering removes the per-target level exactly, so this number cannot be
    improved by predicting the right average affinity for the target. It speaks
    to ordering and nothing else.
    """
    p, t = centered(np.asarray(prediction, dtype=np.float64)), centered(
        np.asarray(truth, dtype=np.float64))
    denominator = float(np.sqrt((p ** 2).sum()) * np.sqrt((t ** 2).sum()))
    return float((p * t).sum() / denominator) if denominator > 1e-12 else 0.0
