"""Paired component bootstrap between two arms across model seeds.

Consumes the `*.rows.jsonl` dump written by `evaluate_arms_ranking.py`. Pairs on
(component, target, k, arm-role) so the contrast is within-episode, then
resamples homology components, which is the unit the cold-target split makes
independent.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def load(path: Path) -> list[dict]:
    return [json.loads(line) for line in
            path.read_text(encoding="utf-8").splitlines() if line.strip()]


def paired_component_effects(rows: list[dict], treatment_prefix: str,
                             control_prefix: str, k: int, arm: str
                             ) -> dict[str, list[float]]:
    """`control - treatment` per component; positive means treatment is better."""
    by_key: dict[tuple, dict[str, list[float]]] = {}
    for row in rows:
        if row["k"] != k or row["arm"] != arm:
            continue
        value = row["mse_pk"]
        if value is None or not np.isfinite(value):
            continue
        name = row["arm_name"]
        if name.endswith(treatment_prefix):
            role = "treatment"
        elif name.endswith(control_prefix):
            role = "control"
        else:
            continue
        seed = name.split("_")[0]
        key = (seed, row["component"], row["target"])
        by_key.setdefault(key, {}).setdefault(role, []).append(float(value))
    effects: dict[str, list[float]] = {}
    for (_, component, _), values in by_key.items():
        if set(values) != {"treatment", "control"}:
            continue
        effects.setdefault(component, []).append(
            float(np.mean(values["control"])) - float(np.mean(values["treatment"])))
    return effects


def bootstrap(effects: dict[str, list[float]], draws: int, seed: int) -> dict:
    components = sorted(effects)
    values = np.asarray([float(np.mean(effects[c])) for c in components])
    if values.size == 0:
        return {"components": 0}
    rng = np.random.default_rng(seed)
    samples = np.asarray([
        values[rng.integers(values.size, size=values.size)].mean()
        for _ in range(draws)])
    return {
        "components": int(values.size),
        "mean_mse_reduction_pk": float(values.mean()),
        "paired_component_95_ci": [float(np.quantile(samples, 0.025)),
                                   float(np.quantile(samples, 0.975))],
        "lower_bound_positive": bool(np.quantile(samples, 0.025) > 0),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("rows", type=Path)
    parser.add_argument("--treatment", required=True,
                        help="suffix identifying the treatment arm names")
    parser.add_argument("--control", required=True)
    parser.add_argument("--arms", nargs="+", default=["full"])
    parser.add_argument("--draws", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=20260815)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = load(args.rows)
    ks = sorted({row["k"] for row in rows})
    report = {"schema": "MetaSieve.PairedArmBootstrap.v1",
              "treatment": args.treatment, "control": args.control,
              "seeds": sorted({row["arm_name"].split("_")[0] for row in rows}),
              "contrasts": []}
    for arm in args.arms:
        for k in ks:
            effects = paired_component_effects(
                rows, args.treatment, args.control, k, arm)
            entry = {"arm": arm, "k": k,
                     **bootstrap(effects, args.draws, args.seed + k)}
            report["contrasts"].append(entry)
            if entry.get("components"):
                low, high = entry["paired_component_95_ci"]
                print(f"{arm} k={k}: {entry['mean_mse_reduction_pk']:+.4f} "
                      f"[{low:+.4f}, {high:+.4f}] "
                      f"components={entry['components']} "
                      f"LB>0={entry['lower_bound_positive']}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")


if __name__ == "__main__":
    main()
