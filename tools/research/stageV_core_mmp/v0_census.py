"""V0 — core-inclusive census, and V1 — interaction variance vs supervision noise.

Both gates are frozen in `PREREGISTRATION.md` and both are inherited verbatim
from Stage U, which froze them before any core-inclusive number existed.

V1 is the decisive cheap measurement: if the between-target variance of
`delta_y` *within one complete transformation* does not exceed the supervision
noise, there is nothing for a protein model to explain and no network is a
rescue.

Run:
    python -m tools.research.stageV_core_mmp.v0_census
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
from tools.research.stageV_core_mmp.core_mmp import (
    Observation, build_observations, double_differences,
    effective_independent_units, load_governed, target_effects,
)

HERE = Path(__file__).resolve().parent
STAGE_T = HERE.parent / "stageT_mmp"

# Frozen, inherited verbatim from Stage U PREREGISTRATION.md section 2.5.
V0_THRESHOLDS = {
    "same_panel_fit_observations": 2000,
    "fit_targets": 50,
    "rich_exact_keys": 30,
    "internal_same_panel_observations": 300,
    "internal_components": 10,
}
DOMINATION_CAPS = {
    "top1_key_share": 0.05, "top10_key_share": 0.20,
    "top1_target_share": 0.25, "top5_target_share": 0.75,
    "top1_component_share": 0.25, "top5_component_share": 0.75,
}
MIN_TARGETS_PER_KEY = 3
MIN_COMPONENTS_PER_KEY = 3
# Frozen, inherited from Stage U section 4.3 / gate 10.
MIN_EVALUABLE_ROWS = 100
BOOTSTRAP_DRAWS = 2000
BOOTSTRAP_SEED = 20260820


def _share(counter: Counter, top: int) -> float:
    total = sum(counter.values())
    if not total:
        return 0.0
    return sum(value for _, value in counter.most_common(top)) / total


def census(observations: list[Observation]) -> dict:
    keys: dict[str, set[str]] = defaultdict(set)
    components: dict[str, set[str]] = defaultdict(set)
    key_counts: Counter = Counter()
    target_counts: Counter = Counter()
    component_counts: Counter = Counter()
    for item in observations:
        keys[item.exact_key].add(item.target)
        components[item.exact_key].add(item.component)
        key_counts[item.exact_key] += 1
        target_counts[item.target] += 1
        component_counts[item.component] += 1
    rich = [k for k in key_counts
            if len(keys[k]) >= MIN_TARGETS_PER_KEY
            and len(components[k]) >= MIN_COMPONENTS_PER_KEY]
    return {
        "observations": len(observations),
        "targets": len(target_counts),
        "components": len(component_counts),
        "exact_keys": len(key_counts),
        "coarse_keys": len({o.coarse_key for o in observations}),
        "keys_with_ge2_targets": sum(1 for v in keys.values() if len(v) >= 2),
        "rich_exact_keys": len(rich),
        "activity_cliffs": sum(1 for o in observations if o.activity_cliff),
        "stereo_edits": sum(1 for o in observations if o.stereo_edit),
        "charge_changing": sum(1 for o in observations if o.charge_change != 0),
        "domination": {
            "top1_key_share": _share(key_counts, 1),
            "top10_key_share": _share(key_counts, 10),
            "top1_target_share": _share(target_counts, 1),
            "top5_target_share": _share(target_counts, 5),
            "top1_component_share": _share(component_counts, 1),
            "top5_component_share": _share(component_counts, 5),
        },
    }


def interaction_variance(effects, sigma2_noise: float, noise_groups: list[list[float]],
                         draws: int = BOOTSTRAP_DRAWS,
                         seed: int = BOOTSTRAP_SEED) -> dict:
    """V1: pooled between-target mean square within a key, minus supervision noise."""
    by_key: dict[str, list] = defaultdict(list)
    for effect in effects:
        by_key[effect.key].append(effect)
    usable = {k: v for k, v in by_key.items() if len(v) >= 2}
    if not usable:
        return {"identifiable": False, "reason": "no exact key has >= 2 targets"}

    def pooled_ms(selection: dict[str, list], weights: dict[str, float]) -> float:
        numerator = denominator = 0.0
        for key, rows in selection.items():
            weight = weights.get(key, 1.0)
            if weight <= 0 or len(rows) < 2:
                continue
            values = np.asarray([r.delta_y for r in rows], dtype=np.float64)
            numerator += weight * float(((values - values.mean()) ** 2).sum())
            denominator += weight * (len(values) - 1)
        return numerator / denominator if denominator > 0 else float("nan")

    ms_effect = pooled_ms(usable, {k: 1.0 for k in usable})
    theta = ms_effect - sigma2_noise

    keys = sorted(usable)
    components = sorted({e.component for rows in usable.values() for e in rows})
    component_index = {name: i for i, name in enumerate(components)}
    rng = np.random.default_rng(seed)
    samples = np.full(draws, np.nan)
    ms_samples = np.full(draws, np.nan)
    noise_samples = np.full(draws, np.nan)
    for draw in range(draws):
        key_counts = np.bincount(rng.integers(0, len(keys), size=len(keys)),
                                 minlength=len(keys))
        component_counts = np.bincount(
            rng.integers(0, len(components), size=len(components)),
            minlength=len(components))
        weights = {}
        for position, key in enumerate(keys):
            rows = usable[key]
            share = float(np.mean([component_counts[component_index[r.component]]
                                   for r in rows]))
            weights[key] = float(key_counts[position]) * share
        ms_b = pooled_ms(usable, weights)
        picked = rng.integers(0, len(noise_groups), size=len(noise_groups))
        residuals: list[float] = []
        for index in picked:
            values = noise_groups[int(index)]
            centre = float(np.median(values))
            residuals.extend(float(v - centre) for v in values)
        noise_b = 2.0 * float(np.var(residuals, ddof=1)) if len(residuals) > 1 \
            else float("nan")
        ms_samples[draw] = ms_b
        noise_samples[draw] = noise_b
        samples[draw] = ms_b - noise_b

    finite = samples[np.isfinite(samples)]
    lo = float(np.quantile(finite, 0.025)) if finite.size else float("nan")
    hi = float(np.quantile(finite, 0.975)) if finite.size else float("nan")
    return {
        "identifiable": True,
        "keys_with_ge2_targets": len(usable),
        "target_effects_used": int(sum(len(v) for v in usable.values())),
        "components": len(components),
        "MS_effect": ms_effect,
        "sigma2_noise": sigma2_noise,
        "theta": theta,
        "ratio_MS_over_noise": ms_effect / sigma2_noise if sigma2_noise else float("nan"),
        "theta_lo": lo,
        "theta_hi": hi,
        "pass": bool(np.isfinite(lo) and lo > 0.0),
        "MS_effect_interval": [
            float(np.quantile(ms_samples[np.isfinite(ms_samples)], 0.025)),
            float(np.quantile(ms_samples[np.isfinite(ms_samples)], 0.975))],
        "sigma2_noise_interval": [
            float(np.quantile(noise_samples[np.isfinite(noise_samples)], 0.025)),
            float(np.quantile(noise_samples[np.isfinite(noise_samples)], 0.975))],
        "draws": int(finite.size),
    }


def noise_groups_from_t0() -> tuple[float, list[list[float]]]:
    """L1 repeated-measure groups from Stage T's T0 provenance cache."""
    from tools.research.stageT_mmp.provenance import group_by_cell, load_cache

    rows = load_cache()
    groups: list[list[float]] = []
    for _cell, cell_rows in group_by_cell(rows).items():
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
                    groups.append(values)
    residuals: list[float] = []
    for values in groups:
        centre = float(np.median(values))
        residuals.extend(float(v - centre) for v in values)
    sigma2_same = 2.0 * float(np.var(residuals, ddof=1))
    return sigma2_same, groups


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "V0_V1_RESULT.json")
    args = parser.parse_args()

    data, seal = load_governed()
    if not seal["isolation"]["physically_isolated"]:
        raise SystemExit("refusing to run without the physical split view")
    fit, internal = partition_components(data)
    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    fit_targets = sorted(t for t, c in component_of.items() if c in set(fit))
    internal_targets = sorted(t for t, c in component_of.items()
                              if c in set(internal))

    built = build_observations(data, fit_targets + internal_targets)
    same_panel = [o for o in built["observations"] if o.same_panel]
    fit_set = set(fit_targets)
    fit_rows = [o for o in same_panel if o.target in fit_set]
    internal_rows = [o for o in same_panel if o.target not in fit_set]

    fit_census = census(fit_rows)
    internal_census = census(internal_rows)

    fit_effects = target_effects(fit_rows)
    internal_effects = target_effects(internal_rows)
    fit_d = double_differences(fit_effects)
    internal_d = double_differences(internal_effects)
    fit_keys = {e.key for e in fit_effects}
    fit_coarse = {e.coarse_key for e in fit_effects}
    repeated = [r for r in internal_d if r.key in fit_keys]
    disjoint = [r for r in internal_d if r.coarse_key not in fit_coarse]

    report: dict = {
        "schema": "MetaSieve.StageV.V0V1.v1",
        "stage": "stageV_core_mmp",
        "meta_test": seal,
        "construction": built["construction"],
        "supersedes": {
            "stage_t_defect": "tools/research/stageT_mmp/CORRECTION_20260817_CORE_KEY.md",
            "stage_u_decision": "STAGE_U_GOVERNANCE_AUDIT.md",
        },
        "census": {"fit_same_panel": fit_census,
                   "internal_same_panel": internal_census},
    }

    checks = {
        "same_panel_fit_observations": fit_census["observations"],
        "fit_targets": fit_census["targets"],
        "rich_exact_keys": fit_census["rich_exact_keys"],
        "internal_same_panel_observations": internal_census["observations"],
        "internal_components": internal_census["components"],
    }
    v0 = {name: {"measured": checks[name], "threshold": value,
                 "pass": bool(checks[name] >= value)}
          for name, value in V0_THRESHOLDS.items()}
    for name, cap in DOMINATION_CAPS.items():
        measured = fit_census["domination"][name]
        v0[name] = {"measured": measured, "cap": cap,
                    "pass": bool(measured <= cap)}
    report["v0_gate"] = {
        "thresholds_inherited_from": "Stage U PREREGISTRATION.md section 2.5",
        "checks": v0,
        "all_pass": bool(all(item["pass"] for item in v0.values())),
    }

    surfaces = {
        "fit_all": effective_independent_units(fit_d),
        "internal_all": effective_independent_units(internal_d),
        "internal_repeated": effective_independent_units(repeated),
        "internal_disjoint": effective_independent_units(disjoint),
    }
    for name, value in surfaces.items():
        value["evaluable"] = bool(value["rows"] >= MIN_EVALUABLE_ROWS)
    report["double_difference_surfaces"] = surfaces
    report["v0b_evaluability"] = {
        "rule_inherited_from": "Stage U PREREGISTRATION.md section 4.3 / gate 10",
        "min_rows": MIN_EVALUABLE_ROWS,
        "evaluable_internal_surfaces": [
            name for name in ("internal_repeated", "internal_disjoint",
                              "internal_all") if surfaces[name]["evaluable"]],
        "primary_surface_internal_repeated_evaluable":
            surfaces["internal_repeated"]["evaluable"],
    }

    sigma2_noise, groups = noise_groups_from_t0()
    report["v1_gate"] = {
        "noise_reference": {
            "sigma2_same_recomputed": sigma2_noise,
            "source": "tools/research/stageT_mmp T0 provenance cache, L1 groups",
            "groups": len(groups),
            "disclosure": ("supervision reliability on a small selected subset; "
                           "NOT an MSE floor"),
        },
        "fit": interaction_variance(fit_effects, sigma2_noise, groups),
        "internal_descriptive_only": interaction_variance(
            internal_effects, sigma2_noise, groups, draws=500),
    }
    report["verdict"] = {
        "v0_pass": report["v0_gate"]["all_pass"],
        "v0b_primary_surface_evaluable":
            surfaces["internal_repeated"]["evaluable"],
        "v1_pass": bool(report["v1_gate"]["fit"].get("pass")),
    }

    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "v0": {k: v["pass"] for k, v in v0.items()},
        "surfaces": {k: {"rows": v["rows"], "EIU": v["effective_independent_units"],
                         "evaluable": v["evaluable"]} for k, v in surfaces.items()},
        "v1": {k: report["v1_gate"]["fit"].get(k) for k in
               ("MS_effect", "sigma2_noise", "theta", "theta_lo", "theta_hi",
                "ratio_MS_over_noise", "keys_with_ge2_targets", "pass")},
        "verdict": report["verdict"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
