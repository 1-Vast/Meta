"""Paired bootstrap between two transports evaluated on the same frozen trunk.

Pairs within (checkpoint, component, target, k), so the only difference between
the two arms is the transport rule. Reports target-level and homology-component
intervals. Seeds are averaged per (component, target) before resampling, so the
intervals are conditional on the trained checkpoints and do not include
retraining variance.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def load(path: Path) -> list[dict]:
    return [json.loads(line) for line in
            path.read_text(encoding="utf-8").splitlines() if line.strip()]


def effects(rows: list[dict], k: int, metric: str, treatment: str,
            control: str) -> dict[tuple[str, str], float]:
    paired: dict[tuple[str, str, str], dict[str, float]] = {}
    for row in rows:
        if row["k"] != k or row["arm"] not in (treatment, control):
            continue
        value = row.get(metric)
        if value is None or not np.isfinite(value):
            continue
        key = (row["arm_name"], row["component"], row["target"])
        paired.setdefault(key, {})[row["arm"]] = float(value)
    pooled: dict[tuple[str, str], list[float]] = {}
    for (_, component, target), values in paired.items():
        if {treatment, control} <= set(values):
            delta = values[treatment] - values[control]
            if metric == "mse_pk":
                delta = -delta          # positive means the treatment helps
            pooled.setdefault((component, target), []).append(delta)
    return {key: float(np.mean(v)) for key, v in pooled.items()}


def bootstrap(values: dict[tuple[str, str], float], level: str,
              draws: int, seed: int) -> dict:
    if not values:
        return {"units": 0}
    if level == "component":
        grouped: dict[str, list[float]] = {}
        for (component, _), value in values.items():
            grouped.setdefault(component, []).append(value)
        sample = np.asarray([float(np.mean(v)) for v in grouped.values()])
    else:
        sample = np.asarray(list(values.values()))
    rng = np.random.default_rng(seed)
    draws_ = np.asarray([sample[rng.integers(sample.size, size=sample.size)].mean()
                         for _ in range(draws)])
    low, high = float(np.quantile(draws_, 0.025)), float(np.quantile(draws_, 0.975))
    return {"units": int(sample.size), "mean": float(sample.mean()),
            "ci95": [low, high], "favours_treatment": bool(low > 0),
            "favours_control": bool(high < 0)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("rows", type=Path)
    parser.add_argument("--contrast", action="append", required=True,
                        help="treatment:control (repeatable)")
    parser.add_argument("--draws", type=int, default=9999)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = load(args.rows)
    ks = sorted({row["k"] for row in rows})
    report = {"schema": "MetaSieve.TransportContrastBootstrap.v1",
              "interval_semantics": ("conditional on the trained checkpoints; "
                                     "seeds averaged before component resampling"),
              "contrasts": []}
    for pair in args.contrast:
        treatment, _, control = pair.partition(":")
        print(f"== {treatment} vs {control} (positive favours {treatment}) ==")
        for k in ks:
            for metric in ("mse_pk", "ci", "spearman"):
                values = effects(rows, k, metric, treatment, control)
                entry = {"treatment": treatment, "control": control, "k": k,
                         "metric": metric,
                         "target_level": bootstrap(values, "target", args.draws, 7 + k),
                         "component_level": bootstrap(values, "component",
                                                      args.draws, 77 + k)}
                report["contrasts"].append(entry)
                comp = entry["component_level"]
                print("  k=%d %-9s component %+.4f [%+.4f,%+.4f] LB>0=%s"
                      % (k, metric, comp["mean"], comp["ci95"][0], comp["ci95"][1],
                         comp["favours_treatment"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")


if __name__ == "__main__":
    main()
