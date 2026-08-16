"""Paired component-level contrasts between R14 arms on the ordering metrics.

`scripts/r14_dispersion_audit.py` writes one row per (arm, seed, target, k).
This turns those rows into the preregistered gate statistics: paired
per-target differences between two arms, aggregated equal-component and
resampled by component.

The primary quantity is the **ordering floor** `Var(y)(1-r²)` rather than the
mean of `r`, because the floor is what enters the MSE and the two can
disagree — targets with a large label spread dominate the floor while
contributing one value to a mean of `r`. Both are reported.

Usage::

    python -m scripts.r14_contrast rows.jsonl --reference A0repro --k 0
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.stageR0_retrieval_falsification import (
    component_bootstrap, component_target_mean,
)

FIELDS = ("shape_ordering_floor", "pearson_r", "shape_pk", "calibration_pk",
          "ci", "sd_pred")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rows", type=Path)
    parser.add_argument("--reference", default="A0repro")
    parser.add_argument("--k", type=int, default=0)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--output", type=Path, default=None)
    arguments = parser.parse_args()

    rows = [json.loads(line) for line in
            arguments.rows.read_text(encoding="utf-8").splitlines() if line.strip()]
    rows = [r for r in rows if int(r["k"]) == arguments.k]
    arms = sorted({r["arm"] for r in rows})
    if arguments.reference not in arms:
        raise SystemExit(f"reference {arguments.reference!r} not among {arms}")

    # Average over seeds per (arm, component, target) before differencing, so
    # the pairing is by target and seed variance is not resampled — the same
    # convention the R6-R12 comparisons use, and its limitation is the same:
    # intervals are conditional on the trained seeds.
    by_key: dict[tuple[str, str, str], dict[str, list[float]]] = defaultdict(
        lambda: defaultdict(list))
    for row in rows:
        key = (row["arm"], row["component"], row["target"])
        for field in FIELDS:
            value = row.get(field)
            if value is not None and np.isfinite(value):
                by_key[key][field].append(float(value))

    seeds = sorted({r["seed"] for r in rows})
    report: dict = {
        "schema": "MetaSieve.R14Contrast.v1",
        "k": arguments.k,
        "reference": arguments.reference,
        "seeds": seeds,
        "note": "shape_ordering_floor is primary; a positive contrast means the "
                "arm has a LOWER floor than the reference, i.e. better ordering",
        "levels": {}, "contrasts": {},
    }

    for arm in arms:
        report["levels"][arm] = {
            field: component_target_mean(
                (component, target, float(np.mean(values[field])))
                for (a, component, target), values in by_key.items()
                if a == arm and values[field])
            for field in FIELDS}

    for arm in arms:
        if arm == arguments.reference:
            continue
        cell = {}
        for field in FIELDS:
            # sign convention: positive = arm is better
            better_when_lower = field in ("shape_ordering_floor", "shape_pk",
                                          "calibration_pk")
            paired = []
            for (a, component, target), values in by_key.items():
                if a != arm or not values[field]:
                    continue
                reference = by_key.get((arguments.reference, component, target))
                if not reference or not reference[field]:
                    continue
                arm_value = float(np.mean(values[field]))
                reference_value = float(np.mean(reference[field]))
                delta = (reference_value - arm_value if better_when_lower
                         else arm_value - reference_value)
                paired.append((component, target, delta))
            cell[field] = component_bootstrap(paired, arguments.bootstrap_draws,
                                              20260816)
        report["contrasts"][f"{arm}_vs_{arguments.reference}"] = cell

    print(f"k={arguments.k}, seeds {seeds}, reference {arguments.reference}\n")
    print(f"{'arm':<12}{'ord floor':>11}{'r':>8}{'shape':>8}{'calib':>8}"
          f"{'CI':>8}{'sd_pred':>9}")
    for arm, level in report["levels"].items():
        print(f"{arm:<12}{level['shape_ordering_floor']:>11.4f}"
              f"{level['pearson_r']:>8.3f}{level['shape_pk']:>8.4f}"
              f"{level['calibration_pk']:>8.4f}{level['ci']:>8.3f}"
              f"{level['sd_pred']:>9.3f}")

    print("\npaired component bootstrap (positive = better than reference):")
    for name, cell in report["contrasts"].items():
        floor, r = cell["shape_ordering_floor"], cell["pearson_r"]
        print(f"  {name}")
        print(f"    ordering floor {floor['mean']:+.4f} "
              f"[{floor['lo']:+.4f}, {floor['hi']:+.4f}]  "
              f"{'RESOLVED' if floor['lo'] > 0 else 'unresolved'}")
        print(f"    pearson r      {r['mean']:+.4f} "
              f"[{r['lo']:+.4f}, {r['hi']:+.4f}]  "
              f"{'RESOLVED' if r['lo'] > 0 else 'unresolved'}")

    if arguments.output is not None:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(json.dumps(report, indent=1), encoding="utf-8")
        print(f"\nwrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
