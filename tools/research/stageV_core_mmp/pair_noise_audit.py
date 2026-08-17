"""Pair-level same-panel delta-pK reliability audit (post-hoc, descriptive).

The preregistered V1 used the Stage T/T0 cell-level noise reference
`sigma2_same = 2 * L1 residual variance = 0.858 pK^2`. That reference assumes
two cell-level measurement errors add independently in a delta; it cannot see
a panel offset shared by both cells of one MMP pair, and it is inflated by
single-cell outliers.

The raw provenance supports a **direct** estimate: for an MMP pair whose two
cells were both measured in more than one shared panel, the per-panel deltas
are repeated realisations of exactly the estimand's supervision. This module
measures that subset, reports its strong selection and its curation-duplication
bias, and recomputes V1 against the alternative noise references.

**Post-hoc:** V0/V1 gates already fired. This file cannot change the verdict
and is excluded from every gate. It exists to make the closure's measurement
basis complete rather than convenient.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.internal_validation import partition_components
from tools.research.stageU_mmp_interaction.observation_cache import load_observations
from tools.research.stageU_mmp_interaction.observations import load_governed
from tools.research.stageU_mmp_interaction.provenance import (
    group_by_cell, load_cache,
)
from tools.research.stageV_core_mmp.core_mmp import Observation, target_effects
from tools.research.stageV_core_mmp.v0_census import interaction_variance

HERE = Path(__file__).resolve().parent
DRAWS = 1000
SEED = 20260820
T0_SIGMA2_SAME = 0.8576301151359423


def _pair_groups(data, rows):
    by_cell = group_by_cell(rows)
    cell_panel: dict[str, dict[str, float]] = defaultdict(dict)
    for cid, cell_rows in by_cell.items():
        by_panel: dict[str, list[float]] = defaultdict(list)
        for row in cell_rows:
            by_panel[row["panel_id"]].append(row["pK"])
        cell_panel[cid] = {p: float(np.median(v)) for p, v in by_panel.items()}
    cell_id = [cell["cell_id"] for cell in data.cells]
    groups = []
    for item in load_observations():
        if not item.same_panel:
            continue
        left = cell_panel.get(cell_id[item.cell_a], {})
        right = cell_panel.get(cell_id[item.cell_b], {})
        values = np.asarray([right[p] - left[p] for p in set(left) & set(right)])
        if values.size >= 2:
            groups.append(values)
    return groups


def _pair_cluster_bootstrap(groups, draws=DRAWS, seed=SEED):
    if not groups:
        return {"groups": 0, "identifiable": False}
    rng = np.random.default_rng(seed)
    samples = []
    for _ in range(draws):
        picked = [groups[int(i)] for i in rng.integers(
            0, len(groups), size=len(groups))]
        residuals = np.concatenate(
            [values - np.median(values) for values in picked])
        samples.append(float(residuals.var(ddof=1)) if residuals.size > 1
                       else float("nan"))
    values = np.asarray(samples)
    finite = values[np.isfinite(values)]
    return {
        "groups": len(groups),
        "identifiable": bool(finite.size),
        "draws": int(finite.size),
        "lo": float(np.quantile(finite, 0.025)),
        "median": float(np.quantile(finite, 0.5)),
        "hi": float(np.quantile(finite, 0.975)),
        "unit": "repeated MMP pair (both cells in a shared panel)",
    }


def _point_variance(groups):
    residuals = np.concatenate([g - np.median(g) for g in groups])
    return float(residuals.var(ddof=1)) if residuals.size > 1 else float("nan")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "PAIR_LEVEL_NOISE_AUDIT.json")
    parser.add_argument("--draws", type=int, default=DRAWS)
    args = parser.parse_args()

    data, seal = load_governed()
    rows = load_cache()
    fit, internal = partition_components(data)
    fit_set = set(fit)

    all_groups = _pair_groups(data, rows)
    zero_groups = [g for g in all_groups if g.max() - g.min() == 0.0]
    disagreeing = [g for g in all_groups if g.max() - g.min() > 0.0]
    report: dict = {
        "schema": "MetaSieve.StageV.PairLevelNoiseAudit.v1",
        "stage": "stageV_core_mmp",
        "disclosure": (
            "post-hoc descriptive audit computed after V0/V1 fired; excluded "
            "from every gate; cannot change the Stage V verdict"),
        "meta_test": seal,
        "population": {
            "same_panel_mmp_pairs": int(sum(
                1 for o in load_observations() if o.same_panel)),
            "repeated_shared_panel_pairs": len(all_groups),
            "zero_range_curation_duplicate_groups": len(zero_groups),
            "disagreeing_groups": len(disagreeing),
            "selection_warning": (
                "only pairs whose two cells were both measured in more than "
                "one shared panel can appear here; they are reference-heavy "
                "and are NOT a random sample of all same-panel MMP pairs"),
            "duplication_warning": (
                "zero-range groups are consistent with one physical "
                "measurement curated twice, not with independent replication; "
                "pooling them deflates the noise estimate"),
        },
        "pair_level_variance": {
            "all_repeated_pairs": {
                "point": _point_variance(all_groups),
                "bootstrap": _pair_cluster_bootstrap(
                    all_groups, args.draws, SEED),
            },
            "disagreeing_pairs_only": {
                "point": _point_variance(disagreeing),
                "bootstrap": _pair_cluster_bootstrap(
                    disagreeing, args.draws, SEED),
            },
        },
    }

    # Recompute V1 against every defensible noise reference, on all keys and
    # on keys spanning distinct components.
    def _stage_v_observation(item):
        return Observation(
            target=item.target, component=item.component, core=item.core,
            exact_key=item.exact_key, coarse_key=item.coarse_key,
            cell_a=item.cell_a, cell_b=item.cell_b, delta_y=item.delta_y,
            same_panel=item.same_panel, stratum=item.stratum,
            tanimoto=item.tanimoto, activity_cliff=item.activity_cliff,
            stereo_edit=item.stereo_edit, charge_change=item.charge_change,
            edit=(0.0,))

    fit_obs = [_stage_v_observation(o) for o in load_observations()
               if o.component in fit_set and o.same_panel]
    effects = target_effects(fit_obs)
    by_key = defaultdict(list)
    for effect in effects:
        by_key[effect.key].append(effect)
    cross_effects = [e for rows in by_key.values()
                     if len({r.component for r in rows}) >= 2 for e in rows]

    # L1 groups for the preregistered T0 reference bootstrap.
    l1_groups: list[list[float]] = []
    for _cid, cell_rows in group_by_cell(rows).items():
        by_panel: dict[str, list[dict]] = defaultdict(list)
        for row in cell_rows:
            by_panel[row["panel_id"]].append(row)
        for panel_rows in by_panel.values():
            by_protocol: dict[tuple, list[float]] = defaultdict(list)
            for row in panel_rows:
                by_protocol[tuple(row.get("assay_protocols") or ())].append(
                    row["pK"])
            for values in by_protocol.values():
                if len(values) >= 2:
                    l1_groups.append(values)

    references = {
        "T0_cell_level_preregistered": (T0_SIGMA2_SAME, l1_groups),
        "pair_level_all_repeated": (
            report["pair_level_variance"]["all_repeated_pairs"]["point"],
            [list(g) for g in all_groups]),
        "pair_level_disagreeing_only": (
            report["pair_level_variance"]["disagreeing_pairs_only"]["point"],
            [list(g) for g in disagreeing]),
    }
    v1 = {}
    for name, (sigma2, groups) in references.items():
        v1[name] = {
            "all_keys": interaction_variance(
                effects, sigma2, groups, draws=args.draws, seed=SEED),
            "cross_component_keys_only": interaction_variance(
                cross_effects, sigma2, groups, draws=args.draws, seed=SEED),
        }
    report["v1_against_alternative_noise_references"] = v1
    report["reading"] = (
        "Under the preregistered cell-level reference the negative is "
        "resolved. The direct pair-level references are lower but heavily "
        "selected and duplicated; under the conservative disagreeing-only "
        "reference no V1 contrast has a positive lower bound, and only the "
        "downward-biased all-group reference would pass cross-component keys. "
        "The interaction variance is therefore not identifiable above the "
        "defensible noise envelope on this corpus.")

    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
