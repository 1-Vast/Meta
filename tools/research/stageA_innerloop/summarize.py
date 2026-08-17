"""Render the Stage A report tables from `STAGE_A_meta_val.json`.

Tables in the report are generated from the leaf artifact rather than
transcribed, so a number in the prose cannot drift from the number in the JSON.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

SUPPORT_SIZES = ("0", "1", "2", "3", "5")
ARMS = ("A0", "A1", "A2")


def cell(value) -> str:
    if value is None:
        return "—"
    if isinstance(value, float) and value != value:
        return "nan"
    return f"{value:.4f}"


def interval(block: dict) -> str:
    return (f"{block['mean']:+.4f} [{block['lo']:+.4f}, {block['hi']:+.4f}]"
            if block and block.get("components") else "—")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    arguments = parser.parse_args()
    payload = json.loads(arguments.input.read_text(encoding="utf-8"))

    def condition_of(arm: str, k: str) -> str:
        if k == "0":
            return "steps0"
        return "steps0" if payload["arm_configs"][arm]["inner_steps"] == 0 else "steps1"

    print("## Per-arm metrics on meta_val (each arm at its own operating condition)\n")
    print("| k | arm | MSE | RMSE | Pearson | Spearman | CI | R2 |")
    print("|---|---|---:|---:|---:|---:|---:|---:|")
    for k in SUPPORT_SIZES:
        for arm in ARMS:
            block = payload["arm_metrics"][arm][k].get(condition_of(arm, k), {})
            print(f"| {k} | {arm} | " + " | ".join(
                cell(block.get(f)) for f in
                ("mse_pk", "rmse_pk", "pearson", "spearman", "ci", "r2")) + " |")

    print("\n## Paired contrasts (component bootstrap, 9999 draws)\n")
    for name in ("A1_vs_A0", "A2_vs_A1", "A2_vs_A0"):
        if name not in payload["contrasts"]:
            continue
        print(f"### {name.replace('_vs_', ' minus ')}\n")
        print("| k | MSE | Pearson | Spearman | CI |")
        print("|---|---|---|---|---|")
        for k in SUPPORT_SIZES:
            row = payload["contrasts"][name][k]
            print(f"| {k} | " + " | ".join(
                interval(row.get(f, {})) for f in
                ("mse_pk", "pearson", "spearman", "ci")) + " |")
        print()

    print("## Counterfactuals (control minus correct MSE; positive = control worse)\n")
    for arm in ("A1", "A2"):
        if arm not in payload.get("counterfactuals", {}):
            continue
        print(f"### {arm}\n")
        controls = ("permuted_support", "matched_wrong_support",
                    "no_adaptation", "wrong_protein", "keep_bias", "keep_weight")
        print("| k | " + " | ".join(controls) + " |")
        print("|---|" + "---|" * len(controls))
        for k in ("1", "2", "3", "5"):
            row = payload["counterfactuals"][arm][k]
            print(f"| {k} | " + " | ".join(
                interval(row.get(c, {})) for c in controls) + " |")
        print()

    print("## Inner-step sweep (MSE at each evaluation-time inner step count)\n")
    print("| k | arm | steps0 | steps1 | steps2 | steps3 |")
    print("|---|---|---:|---:|---:|---:|")
    for k in SUPPORT_SIZES:
        for arm in ARMS:
            block = payload["arm_metrics"][arm][k]
            print(f"| {k} | {arm} | " + " | ".join(
                cell(block.get(f"steps{s}", {}).get("mse_pk")) for s in (0, 1, 2, 3))
                + " |")

    print("\n## Cost\n")
    print("| arm | steps | encoder forwards | wall time (s) | peak GPU (MB) | best step |")
    print("|---|---:|---:|---:|---:|---:|")
    for arm in ARMS:
        report = payload["arm_reports"][arm]
        print(f"| {arm} | {report['optimization_steps']} | "
              f"{report['encoder_forward_passes']} | "
              f"{report['wall_time_seconds']:.0f} | "
              f"{report['peak_cuda_memory_mb']:.0f} | {report['best_step']} |")

    selection = payload["arm_reports"]["A2"].get("selection_log", [])
    cosines = payload["arm_reports"]["A2"].get("gradient_cosine_samples", [])
    if selection:
        import numpy as np
        effective = [row["effective_tasks"] for row in selection]
        print(f"\n## A2 selection\n")
        print(f"- selection steps: {len(selection)}")
        print(f"- effective tasks per step: mean {np.mean(effective):.2f} "
              f"of {selection[0]['candidates']} candidates")
        print(f"- selected-candidate gradient cosine: mean "
              f"{np.mean([r['selected_cosine_mean'] for r in selection]):+.4f}")
        print(f"- all-candidate gradient cosine: mean "
              f"{np.mean([r['all_cosine_mean'] for r in selection]):+.4f}")
    if cosines:
        import numpy as np
        values = np.asarray(cosines)
        print(f"- cosine distribution over {len(values)} candidates: "
              f"mean {values.mean():+.4f}, sd {values.std():.4f}, "
              f"min {values.min():+.4f}, max {values.max():+.4f}, "
              f"fraction > 0: {float((values > 0).mean()):.3f}, "
              f"fraction exactly 0: {float((values == 0).mean()):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
