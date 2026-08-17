"""Render Stage B tables and evaluate the preregistered gates from the artifact.

Gate verdicts are computed here, not transcribed, so the report cannot disagree
with the JSON it cites.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

SUPPORT_SIZES = ("0", "1", "2", "3", "5")
FEW_SHOT = ("1", "2", "3", "5")


def cell(value) -> str:
    if value is None:
        return "—"
    if isinstance(value, float) and value != value:
        return "nan"
    return f"{value:.4f}"


def interval(block: dict) -> str:
    if not block or not block.get("components"):
        return "—"
    return f"{block['mean']:+.4f} [{block['lo']:+.4f}, {block['hi']:+.4f}]"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--candidate", default="C")
    parser.add_argument("--baseline", default="T")
    arguments = parser.parse_args()
    payload = json.loads(arguments.input.read_text(encoding="utf-8"))
    arms = list(payload["arm_metrics"])
    metrics = payload["arm_metrics"]

    print("## Per-arm metrics on meta_val (condition = correct support)\n")
    print("| k | arm | MSE | RMSE | Pearson | Spearman | CI | R2 | centered MSE | cliff |")
    print("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for k in SUPPORT_SIZES:
        for arm in arms:
            block = metrics[arm][k].get("correct", {})
            print(f"| {k} | {arm} | " + " | ".join(
                cell(block.get(f)) for f in
                ("mse_pk", "rmse_pk", "pearson", "spearman", "ci", "r2",
                 "centered_mse_pk", "cliff_sign")) + " |")

    print("\n## Paired contrasts (component bootstrap; within-checkpoint only)\n")
    for name, block in payload["contrasts"].items():
        left, right = name.split("_vs_")
        if right != arguments.baseline and left != arguments.candidate:
            continue
        print(f"### {left} minus {right}\n")
        print("| k | MSE | Pearson | Spearman | CI | centered MSE |")
        print("|---|---|---|---|---|---|")
        for k in SUPPORT_SIZES:
            row = block[k]
            print(f"| {k} | " + " | ".join(
                interval(row.get(f, {})) for f in
                ("mse_pk", "pearson", "spearman", "ci", "centered_mse_pk"))
                + " |")
        print()

    print("## Counterfactuals: control minus correct MSE (positive = control worse)\n")
    controls = ("no_adaptation", "no_transport", "permuted_support",
                "matched_wrong_support", "wrong_protein", "level_only",
                "shape_only")
    for arm in arms:
        print(f"### {arm}\n")
        print("| k | " + " | ".join(controls) + " |")
        print("|---|" + "---|" * len(controls))
        for k in FEW_SHOT:
            row = payload["counterfactuals"][arm][k]
            print(f"| {k} | " + " | ".join(
                interval(row.get(c, {})) for c in controls) + " |")
        print()

    print("## Incremental label dependence (arm minus baseline, on each control)\n")
    print("Positive means the arm loses **more** than the baseline when support")
    print("labels are corrupted, i.e. it depends on correct binding more.\n")
    print("| k | control | " + " | ".join(
        a for a in arms if a != arguments.baseline) + " |")
    print("|---|---|" + "---|" * (len(arms) - 1))
    for k in FEW_SHOT:
        for control in ("permuted_support", "matched_wrong_support"):
            cells = []
            for arm in arms:
                if arm == arguments.baseline:
                    continue
                a = payload["counterfactuals"][arm][k].get(control, {})
                b = payload["counterfactuals"][arguments.baseline][k].get(control, {})
                if not a.get("components") or not b.get("components"):
                    cells.append("—")
                else:
                    cells.append(f"{a['mean'] - b['mean']:+.4f}")
            print(f"| {k} | {control} | " + " | ".join(cells) + " |")

    print("\n## Novelty stratum: low-recall only (max train Tanimoto < 0.4)\n")
    print("| k | " + " | ".join(payload["novelty_strata"]) + " |")
    print("|---|" + "---|" * len(payload["novelty_strata"]))
    for k in FEW_SHOT:
        cells = [interval(block[k]["low_recall_mse"])
                 for block in payload["novelty_strata"].values()]
        print(f"| {k} | " + " | ".join(cells) + " |")

    print("\n## Mechanism magnitudes (mean absolute pK per term)\n")
    print("| k | arm | meta | meta shape | transport | level | complementary |")
    print("|---|---|---:|---:|---:|---:|---:|")
    rows_path = arguments.input.with_suffix(".rows.jsonl")
    if rows_path.exists():
        import collections
        acc = collections.defaultdict(lambda: collections.defaultdict(list))
        with rows_path.open(encoding="utf-8") as handle:
            for line in handle:
                row = json.loads(line)
                if row["condition"] != "correct":
                    continue
                for field in ("meta_abs", "meta_shape_abs", "transport_abs",
                              "level_abs", "complementary_abs"):
                    if row.get(field) is not None:
                        acc[(row["arm"], row["k"])][field].append(row[field])
        for k in SUPPORT_SIZES:
            for arm in arms:
                block = acc.get((arm, int(k)), {})
                print(f"| {k} | {arm} | " + " | ".join(
                    cell(float(np.mean(block[f])) if block.get(f) else None)
                    for f in ("meta_abs", "meta_shape_abs", "transport_abs",
                              "level_abs", "complementary_abs")) + " |")

    # --- gates -------------------------------------------------------------
    print("\n## Preregistered Stage 2 gates\n")
    candidate, baseline = arguments.candidate, arguments.baseline
    if candidate not in metrics or baseline not in metrics:
        print(f"(candidate {candidate} or baseline {baseline} absent)")
        return 0
    few = [metrics[candidate][k]["correct"]["mse_pk"] for k in FEW_SHOT]
    base_few = [metrics[baseline][k]["correct"]["mse_pk"] for k in FEW_SHOT]
    relative = 1.0 - float(np.mean(few)) / float(np.mean(base_few))
    contrast = payload["contrasts"].get(f"{candidate}_vs_{baseline}", {})
    gains = [-contrast[k]["mse_pk"]["mean"] for k in FEW_SHOT]
    k0 = (metrics[candidate]["0"]["correct"]["mse_pk"]
          / metrics[baseline]["0"]["correct"]["mse_pk"] - 1.0)
    verdicts = {
        "G1_mean_mse_gain_ge_5pct": (relative >= 0.05, f"{relative:+.2%}"),
        "G2_non_decreasing_with_support": (
            gains[-1] >= gains[0], f"k1 {gains[0]:+.4f} -> k5 {gains[-1]:+.4f}"),
        "G3_k0_degradation_le_1pct": (k0 <= 0.01, f"{k0:+.2%}"),
        "G4_no_material_ranking_loss": (
            all(contrast[k]["spearman"]["mean"] > -0.02 for k in FEW_SHOT)
            and all(contrast[k]["ci"]["mean"] > -0.02 for k in FEW_SHOT),
            "spearman " + ", ".join(f"{contrast[k]['spearman']['mean']:+.4f}"
                                    for k in FEW_SHOT)),
    }
    if "M" in metrics:
        meta_only = payload["contrasts"].get(f"{candidate}_vs_M", {})
        verdicts["G5_beats_transport_only_and_meta_only"] = (
            all(contrast[k]["mse_pk"]["mean"] < 0 for k in FEW_SHOT)
            and all(meta_only[k]["mse_pk"]["mean"] < 0 for k in FEW_SHOT),
            "vs T " + ", ".join(f"{contrast[k]['mse_pk']['mean']:+.4f}" for k in FEW_SHOT)
            + " | vs M " + ", ".join(f"{meta_only[k]['mse_pk']['mean']:+.4f}"
                                     for k in FEW_SHOT))
    strat = payload["novelty_strata"].get(f"{candidate}_vs_{baseline}", {})
    if strat:
        verdicts["G7_survives_low_recall_stratum"] = (
            all(strat[k]["low_recall_mse"]["mean"] < 0 for k in FEW_SHOT
                if strat[k]["low_recall_mse"].get("components")),
            ", ".join(f"{strat[k]['low_recall_mse']['mean']:+.4f}" for k in FEW_SHOT))

    print("| gate | verdict | value |")
    print("|---|---|---|")
    for name, (passes, detail) in verdicts.items():
        print(f"| {name} | {'PASS' if passes else 'FAIL'} | {detail} |")
    print(f"\n**Overall: {'PASS' if all(v[0] for v in verdicts.values()) else 'FAIL'}**")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
