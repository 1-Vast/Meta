"""Stage 6 admission analysis: within-checkpoint and cross-arm paired bootstrap.

Two contrasts, deliberately separated:

* **within-checkpoint** `full` against the *same* checkpoint's `level_only`.
  This is the causally informative one: identical trunk, identical weights,
  identical episodes, differing only by whether the support weighting is
  applied. It is the contrast that isolates the mechanism.
* **cross-arm** F `full` against A `full`. This compares independently trained
  trunks and therefore mixes the mechanism with training-dynamics variation.

Bootstrap is reported at both the target level and the homology-component level;
the component level is the unit the cold-target split makes independent and is
the binding one.

**Interval semantics.** Seeds are averaged per (component, target) *before*
resampling, so seed-to-seed variance is not resampled. The reported intervals
are therefore conditional on the trained checkpoints: they quantify uncertainty
over homology components given these seeds, not uncertainty over retraining.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


METRICS = ("mse_pk", "ci", "spearman")


def load_rows(path: Path) -> list[dict]:
    return [json.loads(line) for line in
            path.read_text(encoding="utf-8").splitlines() if line.strip()]


def _paired(rows: list[dict], k: int, metric: str,
            treat: tuple[str, str], control: tuple[str, str]) -> dict[tuple[str, str], float]:
    """Return {(component, target): treatment - control} for one seed pair."""
    picked: dict[tuple[str, str], dict[str, list[float]]] = {}
    for row in rows:
        if row["k"] != k:
            continue
        value = row.get(metric)
        if value is None or not np.isfinite(value):
            continue
        key = (row["component"], row["target"])
        if (row["arm_name"], row["arm"]) == treat:
            picked.setdefault(key, {}).setdefault("t", []).append(float(value))
        elif (row["arm_name"], row["arm"]) == control:
            picked.setdefault(key, {}).setdefault("c", []).append(float(value))
    return {key: float(np.mean(v["t"])) - float(np.mean(v["c"]))
            for key, v in picked.items() if {"t", "c"} <= set(v)}


def bootstrap(effects: dict[tuple[str, str], float], level: str,
              draws: int, seed: int) -> dict:
    if not effects:
        return {"units": 0}
    if level == "component":
        grouped: dict[str, list[float]] = {}
        for (component, _), value in effects.items():
            grouped.setdefault(component, []).append(value)
        values = np.asarray([float(np.mean(v)) for v in grouped.values()])
    else:
        values = np.asarray(list(effects.values()))
    rng = np.random.default_rng(seed)
    samples = np.asarray([values[rng.integers(values.size, size=values.size)].mean()
                          for _ in range(draws)])
    low, high = float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))
    return {"units": int(values.size), "mean": float(values.mean()),
            "ci95": [low, high], "excludes_zero": bool(low > 0 or high < 0),
            "favours_treatment": bool(low > 0)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("rows", type=Path)
    parser.add_argument("--seeds", nargs="+",
                        default=["20260812", "20260813", "20260814"])
    parser.add_argument("--treatment-suffix", default="F_similarity_only")
    parser.add_argument("--control-suffix", default="A_grammar")
    parser.add_argument("--draws", type=int, default=9999)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = load_rows(args.rows)
    ks = sorted({row["k"] for row in rows})
    report = {"schema": "MetaSieve.Stage6PairedAnalysis.v1",
              "seeds": args.seeds, "per_seed": [], "aggregate": []}

    # ---- per-seed within-checkpoint and cross-arm point estimates -----------
    for seed in args.seeds:
        for suffix in (args.treatment_suffix, args.control_suffix):
            name = f"{seed}_{suffix}"
            for k in ks:
                entry = {"seed": seed, "arm": suffix, "k": k}
                for metric in METRICS:
                    selected = [r for r in rows if r["arm_name"] == name
                                and r["k"] == k and r["arm"] in ("full", "level_only")
                                and r.get(metric) is not None
                                and np.isfinite(r[metric])]
                    for role in ("full", "level_only"):
                        values = [r[metric] for r in selected if r["arm"] == role]
                        entry[f"{role}_{metric}"] = (
                            float(np.mean(values)) if values else None)
                    if entry.get(f"full_{metric}") is not None:
                        entry[f"delta_{metric}"] = (
                            entry[f"full_{metric}"] - entry[f"level_only_{metric}"])
                perm = [r["mse_pk"] for r in rows if r["arm_name"] == name
                        and r["k"] == k and r["arm"] == "permuted_state"]
                full = [r["mse_pk"] for r in rows if r["arm_name"] == name
                        and r["k"] == k and r["arm"] == "full"]
                entry["permutation_gap"] = (
                    float(np.mean(perm)) - float(np.mean(full)) if perm and full
                    else None)
                report["per_seed"].append(entry)

    # ---- pooled paired bootstrap -------------------------------------------
    for contrast, treat_arm, control_arm, same_checkpoint in (
            ("within_checkpoint_full_vs_level", "full", "level_only", True),
            ("cross_arm_F_full_vs_A_full", "full", "full", False)):
        for k in ks:
            for metric in METRICS:
                pooled: dict[tuple[str, str], list[float]] = {}
                for seed in args.seeds:
                    treat_name = f"{seed}_{args.treatment_suffix}"
                    control_name = (treat_name if same_checkpoint
                                    else f"{seed}_{args.control_suffix}")
                    effects = _paired(rows, k, metric,
                                      (treat_name, treat_arm),
                                      (control_name, control_arm))
                    for key, value in effects.items():
                        pooled.setdefault(key, []).append(value)
                averaged = {key: float(np.mean(v)) for key, v in pooled.items()}
                # For MSE lower is better, so flip the sign to make "positive
                # means the treatment helps" hold for every metric.
                signed = ({key: -value for key, value in averaged.items()}
                          if metric == "mse_pk" else averaged)
                report["aggregate"].append({
                    "contrast": contrast, "k": k, "metric": metric,
                    "direction": ("mse reduction" if metric == "mse_pk"
                                  else f"{metric} increase"),
                    "target_level": bootstrap(signed, "target", args.draws, 11 + k),
                    "component_level": bootstrap(signed, "component", args.draws, 91 + k),
                })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")

    print("== within-checkpoint: F full vs the SAME checkpoint's level-only ==")
    for entry in report["aggregate"]:
        if entry["contrast"] != "within_checkpoint_full_vs_level" or entry["k"] < 2:
            continue
        comp, tgt = entry["component_level"], entry["target_level"]
        print(f"  k={entry['k']} {entry['metric']:9s} "
              f"component {comp['mean']:+.4f} [{comp['ci95'][0]:+.4f},{comp['ci95'][1]:+.4f}] "
              f"LB>0={comp['favours_treatment']} | "
              f"target {tgt['mean']:+.4f} LB>0={tgt['favours_treatment']}")
    print("== cross-arm: F full vs A full ==")
    for entry in report["aggregate"]:
        if entry["contrast"] != "cross_arm_F_full_vs_A_full" or entry["metric"] != "mse_pk":
            continue
        comp = entry["component_level"]
        print(f"  k={entry['k']} mse reduction component {comp['mean']:+.4f} "
              f"[{comp['ci95'][0]:+.4f},{comp['ci95'][1]:+.4f}] "
              f"LB>0={comp['favours_treatment']}")


if __name__ == "__main__":
    main()
