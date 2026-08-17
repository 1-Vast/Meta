"""Crossed double-difference banks: D(tau, t1, t2) = dy(t1,tau) - dy(t2,tau).

Why this estimand. `delta_y(t, tau) = mu_tau + delta(t, tau) + noise`. The
double difference cancels the target-level affinity offset **and** the generic
chemical effect `mu_tau` exactly, so a model that scores on `D` cannot be
scoring on target level, on generic medicinal chemistry, or on a target-identity
key. Those are the three shortcuts that explained every earlier positive in this
repository.

Construction rules, all label-blind except the aggregation of `delta_y` itself:

* `delta_y(t, tau)` is the **median** of that target's observations of `tau`
  (a target can realise one transformation through more than one core);
* only observations in the same-panel strata (S1/S2) enter the primary bank;
  S3 cross-panel observations are a separate weak bank and are never pooled;
* a `D` row exists for every unordered target pair sharing `tau`, in canonical
  target-id sort order. The model is exactly antisymmetric, so one orientation
  carries all the information and a test verifies the other is its negation;
* `D` rows are **never** formed across different transformation keys.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.qpsmp_data import stable_seed
from tools.research.stageT_mmp.mmp import (
    DESCRIPTOR_WIDTH, Transformation, descriptor,
)
from tools.research.stageT_mmp.observations import MMPObservation

SAME_PANEL_STRATA = ("S1_same_panel_single", "S2_same_panel_multi")


@dataclass(frozen=True)
class TargetEffect:
    """delta_y(t, tau) after robust within-target aggregation."""

    key: str
    target: str
    component: str
    delta_y: float
    observations: int
    stratum: str
    activity_cliff: bool
    descriptor: tuple[float, ...]


@dataclass(frozen=True)
class DoubleDifference:
    key: str
    target_left: str
    target_right: str
    component_left: str
    component_right: str
    value: float
    observations: int
    cross_component: bool
    activity_cliff: bool
    descriptor: tuple[float, ...]

    @property
    def row_id(self) -> str:
        return f"{self.key}|{self.target_left}|{self.target_right}"


def _descriptor_for(item: MMPObservation) -> tuple[float, ...]:
    """Rebuild the structured descriptor from the transformation key parts."""
    context = eval(item.exact_key.split("|", 1)[0])  # noqa: S307 - our own repr
    edit = item.exact_key.split("|", 1)[1]
    r_a, r_b = edit.split(">>", 1)
    built = Transformation(
        core=item.core, r_a=r_a, r_b=r_b, context=context,
        coarse_context=(context[0], context[1]),
        charge_change=item.charge_change, stereo_edit=item.stereo_edit)
    return tuple(descriptor(built))


def target_effects(observations: list[MMPObservation],
                   same_panel_only: bool = True) -> list[TargetEffect]:
    grouped: dict[tuple[str, str], list[MMPObservation]] = defaultdict(list)
    for item in observations:
        if same_panel_only and item.stratum not in SAME_PANEL_STRATA:
            continue
        grouped[(item.exact_key, item.target)].append(item)
    out: list[TargetEffect] = []
    for (key, target), rows in sorted(grouped.items()):
        values = [row.delta_y for row in rows]
        out.append(TargetEffect(
            key=key, target=target, component=rows[0].component,
            delta_y=float(np.median(values)), observations=len(rows),
            stratum=min(row.stratum for row in rows),
            activity_cliff=any(row.activity_cliff for row in rows),
            descriptor=_descriptor_for(rows[0])))
    return out


def double_differences(effects: list[TargetEffect]) -> list[DoubleDifference]:
    by_key: dict[str, list[TargetEffect]] = defaultdict(list)
    for effect in effects:
        by_key[effect.key].append(effect)
    out: list[DoubleDifference] = []
    for key in sorted(by_key):
        rows = sorted(by_key[key], key=lambda value: value.target)
        for left, right in combinations(rows, 2):
            out.append(DoubleDifference(
                key=key, target_left=left.target, target_right=right.target,
                component_left=left.component, component_right=right.component,
                value=left.delta_y - right.delta_y,
                observations=left.observations + right.observations,
                cross_component=left.component != right.component,
                activity_cliff=left.activity_cliff or right.activity_cliff,
                descriptor=left.descriptor))
    return out


def shuffle_within_key(rows: list[DoubleDifference], seed: int
                       ) -> list[DoubleDifference]:
    """Arm F: permute D **inside** each transformation key.

    The marginal distribution of `D` for every key is preserved exactly; what is
    destroyed is which protein pair got which value. A model that still scores
    after this is reading something other than protein x transformation.
    """
    by_key: dict[str, list[int]] = defaultdict(list)
    for position, row in enumerate(rows):
        by_key[row.key].append(position)
    out = list(rows)
    for key in sorted(by_key):
        positions = by_key[key]
        if len(positions) < 2:
            continue
        rng = np.random.default_rng(stable_seed("stageT-label-shuffle", seed, key))
        values = [rows[p].value for p in positions]
        order = rng.permutation(len(values))
        for slot, position in enumerate(positions):
            row = rows[position]
            out[position] = DoubleDifference(
                row.key, row.target_left, row.target_right, row.component_left,
                row.component_right, float(values[int(order[slot])]),
                row.observations, row.cross_component, row.activity_cliff,
                row.descriptor)
    return out


def split_by_key_overlap(rows: list[DoubleDifference], known_keys: set[str]
                         ) -> tuple[list[DoubleDifference], list[DoubleDifference]]:
    """Surface 1 (repeated keys) and surface 2 (transformation-disjoint)."""
    repeated = [row for row in rows if row.key in known_keys]
    disjoint = [row for row in rows if row.key not in known_keys]
    return repeated, disjoint


def effective_independent_units(rows: list[DoubleDifference]) -> dict:
    components = {c for row in rows
                  for c in (row.component_left, row.component_right)}
    keys = {row.key for row in rows}
    return {
        "rows": len(rows),
        "protein_components": len(components),
        "transformation_keys": len(keys),
        "effective_independent_units": min(len(components), len(keys)),
        "note": ("rows sharing a component or a key are correlated; the binding "
                 "constraint on power is the smaller of the two cluster counts"),
    }


DESCRIPTOR_DIM = DESCRIPTOR_WIDTH
