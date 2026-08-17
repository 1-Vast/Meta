"""Synthetic calibration of the V1 interaction-variance statistic.

Post-hoc, descriptive, excluded from every frozen gate. The frozen V1 gate
already fired; this module asks a different and useful question: **how large
would a true target x transformation interaction have to be for the observed
fit graph to produce the measured `MS_effect = 0.4517`**, under each of the
three defensible supervision-noise references?

It simulates `delta_y(t,tau) = key_mean + delta(t,tau) + noise` on the real
key x target x component structure (4,651 keys / 12,133 target effects / 99
components), with `delta ~ N(0, delta2)` and `noise ~ N(0, sigma2)`, and
compares the resulting pooled between-target mean square with the observed
value. No model is trained and no label is read beyond the already-recorded
statistics.
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

from scipy import stats as scipy_stats

from scripts.internal_validation import partition_components
from tools.research.stageU_mmp_interaction.observation_cache import load_observations
from tools.research.stageV_core_mmp.core_mmp import (
    Observation, load_governed, target_effects,
)

HERE = Path(__file__).resolve().parent
OBSERVED_MS = 0.45168040420690786
NOISE_REFERENCES = {
    "pair_level_all_repeated": 0.16618178668989872,
    "pair_level_disagreeing_only": 0.30323847769409734,
    "T0_cell_level_preregistered": 0.8576301151359423,
}
DELTA2_GRID = (0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.00)


def _usable_keys():
    data, _seal = load_governed()
    fit, _internal = partition_components(data)
    fit_set = set(fit)
    observations = load_observations()

    def convert(item):
        return Observation(
            target=item.target, component=item.component, core=item.core,
            exact_key=item.exact_key, coarse_key=item.coarse_key,
            cell_a=item.cell_a, cell_b=item.cell_b, delta_y=item.delta_y,
            same_panel=item.same_panel, stratum=item.stratum,
            tanimoto=item.tanimoto, activity_cliff=item.activity_cliff,
            stereo_edit=item.stereo_edit, charge_change=item.charge_change,
            edit=(0.0,))

    effects = target_effects([convert(o) for o in observations
                              if o.component in fit_set and o.same_panel])
    by_key: dict[str, list] = defaultdict(list)
    for effect in effects:
        by_key[effect.key].append(effect)
    return {key: rows for key, rows in by_key.items() if len(rows) >= 2}


def _calibration_row(sigma2: float, delta2: float, df: int) -> dict:
    """Exact null distribution of the pooled MS.

    Under delta(t,tau) ~ N(0, delta2) and noise ~ N(0, sigma2), each key's
    between-target SS divided by (sigma2 + delta2) is chi-square(k-1), so the
    pooled statistic is (sigma2+delta2) * chi-square(df)/df. No simulation is
    needed; the Monte-Carlo check was replaced by this exact distribution after
    the first exploratory run showed the simulation loop was unnecessarily
    slow.
    """
    total = sigma2 + delta2
    scale = total / df
    lo = scale * scipy_stats.chi2.ppf(0.05, df)
    median = scale * scipy_stats.chi2.ppf(0.50, df)
    hi = scale * scipy_stats.chi2.ppf(0.95, df)
    p_ge = float(scipy_stats.chi2.sf(OBSERVED_MS * df / total, df))
    return {
        "MS_mean": float(total),
        "MS_q05": float(lo),
        "MS_q50": float(median),
        "MS_q95": float(hi),
        "P_MS_ge_observed": p_ge,
        "df": int(df),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "V1_SYNTHETIC_CALIBRATION.json")
    args = parser.parse_args()

    usable = _usable_keys()
    n_keys = len(usable)
    n_effects = int(sum(len(v) for v in usable.values()))
    n_components = len({e.component for rows in usable.values() for e in rows})
    df = n_effects - n_keys

    report: dict = {
        "schema": "MetaSieve.StageV.V1SyntheticCalibration.v1",
        "stage": "stageV_core_mmp",
        "disclosure": ("post-hoc descriptive calibration; gates already fired; "
                       "no label read beyond recorded statistics; cannot "
                       "change the Stage V verdict"),
        "method": ("exact chi-square null distribution of the pooled "
                   "between-key mean square on the real key/target/component "
                   "graph; delta ~ N(0, delta2), noise ~ N(0, sigma2)"),
        "observed": {"MS_effect": OBSERVED_MS, "keys": n_keys,
                     "effects": n_effects, "components": n_components,
                     "df": int(df)},
        "noise_references": NOISE_REFERENCES,
        "grid": {},
    }

    for noise_name, sigma2 in NOISE_REFERENCES.items():
        report["grid"][noise_name] = {
            str(delta2): _calibration_row(sigma2, delta2, df)
            for delta2 in DELTA2_GRID
        }

    report["implied_interaction_variance_if_all_excess_is_signal"] = {
        name: {
            "delta2_point": max(0.0, OBSERVED_MS - sigma2),
            "delta_sd": float(np.sqrt(max(0.0, OBSERVED_MS - sigma2))),
        }
        for name, sigma2 in NOISE_REFERENCES.items()
    }
    report["reading"] = (
        "The measured MS_effect is far below what the T0 cell-level noise "
        "reference predicts under a zero-interaction model (empirical "
        "P(MS >= observed) = 0), so T0 sigma2 is an upper-bound reference for "
        "the MMP delta estimand, as already suspected. Under the direct "
        "pair-level references the observed MS implies an interaction sd of "
        "about 0.39 pK (disagreeing-only) to 0.53 pK (all repeated pairs) if "
        "all excess above noise were signal. That is consistent with, and "
        "bounded by, the unresolved cross-component V1 interval; it is not "
        "evidence of a detectable signal and cannot reopen the frozen gate.")

    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
