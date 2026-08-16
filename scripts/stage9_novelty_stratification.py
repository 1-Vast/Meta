"""Stratify the k=0 decomposition by ligand novelty and target distance.

The CD-HIT40 split is component-hard on **proteins only**. Ligands may recur
across splits, so a retrieval predictor can be partly memorising a ligand's
typical potency. This measures how the retrieval advantage behaves as the query
ligands become genuinely novel, which decides whether the advantage is a real
capability or a recall artefact.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def component_target_mean(rows: list[dict], field: str) -> float:
    by_target: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        value = row.get(field)
        if value is None or not np.isfinite(value):
            continue
        by_target.setdefault((row["component"], row["target"]), []).append(float(value))
    by_component: dict[str, list[float]] = {}
    for (component, _), values in by_target.items():
        by_component.setdefault(component, []).append(float(np.mean(values)))
    return (float(np.mean([np.mean(v) for v in by_component.values()]))
            if by_component else float("nan"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("rows", type=Path)
    parser.add_argument("--estimator", action="append", required=True)
    parser.add_argument("--variable", default="ligand_novelty")
    parser.add_argument("--edges", type=float, nargs="+",
                        default=[0.0, 0.4, 0.6, 0.8, 1.01])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    rows = [json.loads(line) for line in
            args.rows.read_text(encoding="utf-8").splitlines() if line.strip()]
    report = {"schema": "MetaSieve.K0NoveltyStratification.v1",
              "variable": args.variable, "edges": args.edges, "buckets": []}
    print("%-46s %6s %8s %8s %8s %7s %7s" % (
        "estimator", "n", args.variable, "MSE", "calib", "CI", "rho"))
    for low, high in zip(args.edges[:-1], args.edges[1:]):
        for estimator in args.estimator:
            selected = [r for r in rows if r["estimator"] == estimator
                        and low <= r[args.variable] < high]
            if not selected:
                continue
            entry = {
                "estimator": estimator, "low": low, "high": high,
                "targets": len({(r["component"], r["target"]) for r in selected}),
                "mean_variable": float(np.mean([r[args.variable] for r in selected])),
                "mse_pk": component_target_mean(selected, "mse_pk"),
                "calibration_pk": component_target_mean(selected, "calibration_pk"),
                "shape_pk": component_target_mean(selected, "shape_pk"),
                "ci": component_target_mean(selected, "ci"),
                "spearman": component_target_mean(selected, "spearman"),
            }
            report["buckets"].append(entry)
            print("%-46s %6d %8.3f %8.4f %8.4f %7.4f %7.4f" % (
                f"[{low:.2f},{high:.2f}) {estimator}", entry["targets"],
                entry["mean_variable"], entry["mse_pk"], entry["calibration_pk"],
                entry["ci"], entry["spearman"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")


if __name__ == "__main__":
    main()
