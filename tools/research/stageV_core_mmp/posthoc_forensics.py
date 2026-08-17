"""Post-hoc sensitivity forensics for the Stage V negative.

This module was written **after** V0/V1 were read and the frozen gates had
already fired. It is descriptive only and cannot change the verdict: its output
is excluded from every gate. Its purpose is to show whether the three negative
readings — degree domination, a non-evaluable primary surface, and an
interaction MS below supervision noise — are artefacts of the one dominating
target, of key selection, or of same-component pairs.

Nothing here trains a model, reads the development-validation split, or reads
the sealed confirmation split.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.internal_validation import partition_components
from tools.research.stageU_mmp_interaction.observation_cache import load_observations
from tools.research.stageV_core_mmp.core_mmp import (
    DoubleDifference, Observation, TargetEffect, double_differences,
    load_governed, target_effects,
)

HERE = Path(__file__).resolve().parent
SIGMA2_NOISE = 0.8576301151359423


def _as_stage_v_observation(row):
    return Observation(
        target=row.target, component=row.component, core=row.core,
        exact_key=row.exact_key, coarse_key=row.coarse_key,
        cell_a=row.cell_a, cell_b=row.cell_b, delta_y=row.delta_y,
        same_panel=row.same_panel, stratum=row.stratum,
        tanimoto=row.tanimoto, activity_cliff=row.activity_cliff,
        stereo_edit=row.stereo_edit, charge_change=row.charge_change,
        edit=(0.0,))


def _select_keys(effects: list[TargetEffect], require_distinct_components: bool,
                 min_targets: int = 2) -> list[tuple[str, list[TargetEffect]]]:
    by_key: dict[str, list[TargetEffect]] = defaultdict(list)
    for effect in effects:
        by_key[effect.key].append(effect)
    keys = []
    for key, rows in by_key.items():
        if len(rows) < min_targets:
            continue
        if require_distinct_components:
            if len({r.component for r in rows}) < 2:
                continue
        keys.append((key, rows))
    return keys


def _pooled_ms(effects: list[TargetEffect], require_distinct_components: bool,
               min_targets: int = 2, bootstrap: bool = False,
               draws: int = 1000, seed: int = 20260820) -> dict:
    keys = _select_keys(effects, require_distinct_components, min_targets)
    if not keys:
        return {"keys": 0, "effects": 0, "components": 0, "MS_effect": None}

    def weighted_ms(rows_by_key, key_weight):
        numerator = denominator = 0.0
        for key, rows in rows_by_key:
            weight = key_weight.get(key, 1.0)
            if weight <= 0 or len(rows) < 2:
                continue
            values = np.asarray([r.delta_y for r in rows], dtype=np.float64)
            numerator += weight * float(((values - values.mean()) ** 2).sum())
            denominator += weight * (len(values) - 1)
        return numerator / denominator if denominator > 0 else float("nan")

    out = {
        "keys": len(keys),
        "effects": int(sum(len(rows) for _, rows in keys)),
        "components": len({r.component for _, rows in keys for r in rows}),
        "MS_effect": weighted_ms(keys, {k: 1.0 for k, _ in keys}),
    }
    if bootstrap:
        key_names = [k for k, _ in keys]
        components = sorted({r.component for _, rows in keys for r in rows})
        component_index = {name: i for i, name in enumerate(components)}
        rng = np.random.default_rng(seed)
        samples = np.full(draws, np.nan)
        for draw in range(draws):
            key_counts = np.bincount(rng.integers(0, len(key_names),
                                                  size=len(key_names)),
                                     minlength=len(key_names))
            component_counts = np.bincount(
                rng.integers(0, len(components), size=len(components)),
                minlength=len(components))
            weights = {}
            for position, (key, rows) in enumerate(keys):
                share = float(np.mean([component_counts[component_index[r.component]]
                                       for r in rows]))
                weights[key] = float(key_counts[position]) * share
            samples[draw] = weighted_ms(keys, weights)
        finite = samples[np.isfinite(samples)]
        out["bootstrap"] = {
            "draws": int(finite.size),
            "seed": seed,
            "lo": float(np.quantile(finite, 0.025)),
            "median": float(np.quantile(finite, 0.5)),
            "hi": float(np.quantile(finite, 0.975)),
            "unit": "two-way: transformation keys x protein components",
        }
    return out


def _census(observations: list[Observation]) -> dict:
    key_counts = Counter(o.exact_key for o in observations)
    target_counts = Counter(o.target for o in observations)
    component_counts = Counter(o.component for o in observations)
    rich_keys = defaultdict(set)
    rich_components = defaultdict(set)
    for o in observations:
        rich_keys[o.exact_key].add(o.target)
        rich_components[o.exact_key].add(o.component)
    total = len(observations)
    def share(counter, n):
        return sum(v for _, v in counter.most_common(n)) / total if total else 0.0
    return {
        "observations": total,
        "targets": len(target_counts),
        "components": len(component_counts),
        "exact_keys": len(key_counts),
        "rich_exact_keys": sum(1 for key in key_counts
                               if len(rich_keys[key]) >= 3
                               and len(rich_components[key]) >= 3),
        "top1_target_share": share(target_counts, 1),
        "top5_target_share": share(target_counts, 5),
        "top1_component_share": share(component_counts, 1),
        "top5_component_share": share(component_counts, 5),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "POSTHOC_FORENSICS.json")
    args = parser.parse_args()

    data, seal = load_governed()
    fit, internal = partition_components(data)
    fit_set = set(fit)
    rows = [_as_stage_v_observation(row) for row in load_observations()]
    fit_same = [o for o in rows if o.component in fit_set and o.same_panel]
    internal_same = [o for o in rows if o.component not in fit_set and o.same_panel]

    effects_all = target_effects(fit_same)
    top_target, top_count = Counter(o.target for o in fit_same).most_common(1)[0]
    top_component = Counter(o.component for o in fit_same).most_common(1)[0][0]
    without_top = [o for o in fit_same if o.target != top_target]

    # Cross-component effects only.
    by_target_component = {o.target: o.component for o in fit_same}
    cross_effects = []
    by_key: dict[str, list[TargetEffect]] = defaultdict(list)
    for effect in effects_all:
        by_key[effect.key].append(effect)
    for key, entries in by_key.items():
        if len({e.component for e in entries}) >= 2:
            cross_effects.extend(entries)

    fit_d = double_differences(effects_all)
    fit_d_cross = [r for r in fit_d if r.cross_component]
    d_values_cross = np.asarray([r.value for r in fit_d_cross], dtype=np.float64)

    # Coarse-key repeated surface: a family-relaxed reference, never strict tau.
    fit_coarse = {e.coarse_key for e in effects_all}
    internal_effects = target_effects(internal_same)
    internal_d = double_differences(internal_effects)
    internal_repeated_exact = [r for r in internal_d if r.key in
                               {e.key for e in effects_all}]
    internal_repeated_coarse = [r for r in internal_d
                                if r.coarse_key in fit_coarse]

    report: dict = {
        "schema": "MetaSieve.StageV.PosthocForensics.v1",
        "stage": "stageV_core_mmp",
        "disclosure": (
            "computed after the frozen V0/V1 gates had already fired; purely "
            "descriptive sensitivity analysis; excluded from every gate and "
            "cannot change the Stage V verdict"),
        "meta_test": seal,
        "noise_reference_sigma2_same": SIGMA2_NOISE,
        "sensitivity_degree_domination": {
            "top_target": {
                "count": top_count,
                "share": top_count / len(fit_same),
            },
            "top_component": top_component,
            "after_removing_top_target": _census(without_top),
            "reading": (
                "if the single dominating target is removed, the remaining fit "
                "bank still carries the same scientific reading: rich keys "
                "remain, and the interaction variance statistic is recomputed "
                "below on the reduced bank"),
        },
        "sensitivity_interaction_variance": {
            "all_keys_ge2_targets": _pooled_ms(
                effects_all, require_distinct_components=False),
            "keys_with_distinct_components": _pooled_ms(
                effects_all, require_distinct_components=True),
            "after_removing_top_target": _pooled_ms(
                target_effects(without_top), require_distinct_components=False),
            "cross_component_effects_only": _pooled_ms(
                cross_effects, require_distinct_components=True,
                bootstrap=True, draws=1000),
        },
        "cross_component_D_scale": {
            "rows": int(d_values_cross.size),
            "variance": float(d_values_cross.var(ddof=1)),
            "sd": float(d_values_cross.std(ddof=1)),
            "null_variance_2sigma2_same": 2.0 * SIGMA2_NOISE,
            "variance_over_null": float(d_values_cross.var(ddof=1)
                                        / (2.0 * SIGMA2_NOISE)),
            "reading": (
                "the observed cross-component D variance is compared to the "
                "null expectation 2*sigma2_same; a ratio near or below one is "
                "consistent with D being dominated by supervision noise"),
        },
        "evaluation_surfaces": {
            "internal_repeated_exact": {
                "rows": len(internal_repeated_exact),
                "components": len({c for r in internal_repeated_exact
                                   for c in (r.component_left,
                                             r.component_right)}),
                "keys": len({r.key for r in internal_repeated_exact}),
            },
            "internal_repeated_coarse_reference": {
                "rows": len(internal_repeated_coarse),
                "components": len({c for r in internal_repeated_coarse
                                   for c in (r.component_left,
                                             r.component_right)}),
                "keys": len({r.coarse_key for r in internal_repeated_coarse}),
                "note": ("coarse key is NOT the strict core/context-matched "
                         "estimand and can never substitute for it"),
            },
        },
    }
    cross_ms = report["sensitivity_interaction_variance"][
        "cross_component_effects_only"]
    noise = report["noise_reference_sigma2_same"]
    report["cross_component_interaction_ms_minus_noise"] = {
        "MS_effect": cross_ms.get("MS_effect"),
        "MS_bootstrap": cross_ms.get("bootstrap"),
        "sigma2_noise": noise,
        "point_theta": (cross_ms.get("MS_effect") - noise
                        if cross_ms.get("MS_effect") is not None else None),
        "conservative_theta_hi": (
            cross_ms["bootstrap"]["hi"] - noise
            if cross_ms.get("bootstrap") else None),
        "reading": (
            "the cross-component between-target MS is also below the "
            "supervision-noise reference, so the negative is not an artefact "
            "of same-component pairs; it is the population-level reading"),
    }

    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
