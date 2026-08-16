"""Render the admission table for a governed nested-k evaluation directory."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


ARMS = ("full", "zero_shot", "level_only", "sar_cut", "permuted_state",
        "foreign_code_state", "wrong_protein_state")


def component_target_metric(rows: list[dict], metric: str) -> float:
    target_values: dict[tuple[str, str], list[float]] = {}
    for row in rows:
        if not np.isfinite(row[metric]):
            continue
        target_values.setdefault((row["component"], row["target"]), []).append(
            float(row[metric]))
    component_values: dict[str, list[float]] = {}
    for (component, _), values in target_values.items():
        component_values.setdefault(component, []).append(float(np.mean(values)))
    if not component_values:
        return float("nan")
    return float(np.mean([np.mean(v) for v in component_values.values()]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    result = json.loads((args.directory / "RESULT.json").read_text(encoding="utf-8"))
    rows = [json.loads(line) for line in
            (args.directory / "PREDICTIONS.jsonl").read_text(
                encoding="utf-8").splitlines() if line.strip()]
    seeds = sorted({row["model_seed"] for row in rows})
    ks = sorted({row["k"] for row in rows})

    lines = [f"seeds: {seeds}", f"targets: {result['targets']}",
             f"components: {result['components']}",
             f"parameters: {result['training'][0].get('trainable_parameters')}"]

    lines.append("\n## Pooled metrics (equal component, then target, then draw)\n")
    lines.append("| k | arm | MSE | RMSE | CI | Spearman |")
    lines.append("|---|---|---:|---:|---:|---:|")
    for k in ks:
        for arm in ARMS:
            selected = [r for r in rows if r["k"] == k and r["arm"] == arm]
            if not selected:
                continue
            lines.append("| {} | {} | {:.4f} | {:.4f} | {:.4f} | {:.4f} |".format(
                k, arm,
                component_target_metric(selected, "mse_pk"),
                component_target_metric(selected, "rmse_pk"),
                component_target_metric(selected, "ci"),
                component_target_metric(selected, "spearman")))

    lines.append("\n## Per-seed full-arm MSE\n")
    lines.append("| k | " + " | ".join(f"seed {s}" for s in seeds) + " |")
    lines.append("|---" * (len(seeds) + 1) + "|")
    for k in ks:
        cells = []
        for seed in seeds:
            selected = [r for r in rows if r["k"] == k and r["arm"] == "full"
                        and r["model_seed"] == seed]
            cells.append(f"{component_target_metric(selected, 'mse_pk'):.4f}")
        lines.append(f"| {k} | " + " | ".join(cells) + " |")

    lines.append("\n## Paired component bootstrap, full versus control\n")
    lines.append("| k | control | mean MSE reduction | 95% CI | lower bound > 0 |")
    lines.append("|---|---|---:|---|---|")
    for contrast in result["contrasts"]:
        low, high = contrast["paired_component_95_ci"]
        lines.append("| {} | {} | {:+.4f} | [{:+.4f}, {:+.4f}] | {} |".format(
            contrast["k"], contrast["control"],
            contrast["mean_mse_reduction_pk"], low, high,
            "yes" if low > 0 else "no"))

    text = "\n".join(lines)
    print(text)
    if args.output:
        args.output.write_text(text + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
