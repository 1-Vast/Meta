"""Stage 0 feasibility analysis: where the k=0 error actually lives, and what
MSE <= 1.00 pK^2 would require.

Every squared error decomposes exactly into a target-level term and a
within-target shape term:

    MSE = (mean(p) - mean(y))^2 + mean((p - mean(p)) - (y - mean(y)))^2
        = level^2 + centered_MSE

The two are attacked by completely different mechanisms — level by calibrating
an unseen protein's mean affinity, shape by ordering that protein's ligands — so
knowing the split determines which of them a candidate must move, and whether
the target is reachable at all.

Three oracles bound the problem without training anything:

* `oracle_level`   — replace the predicted target mean with the true one;
                     the remaining error is exactly `centered_MSE`.
* `oracle_shape`   — replace the centered prediction with the centered truth;
                     the remaining error is exactly `level^2`.
* `constant`       — predict each target's own mean; its `centered_MSE` is the
                     within-target label variance, the floor a level-only
                     mechanism can reach and the ceiling any ordering mechanism
                     must beat.

No training, no meta_test.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData                              # noqa: E402
from scripts.stageR0_retrieval_falsification import (                 # noqa: E402
    component_bootstrap, component_target_mean,
)
from scripts.stageR6_compare_arms import load_arm                     # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, normalized_episode, training_label_scale,
)
from tools.research.stageA_innerloop.train_meta import (              # noqa: E402
    align_atoms, encode_parts, episode_tensors,
)
from tools.research.stageB_complementary.arms import (                # noqa: E402
    StageBAdaptation, predict,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
SUPPORT_SIZES = (0, 1, 2, 3, 5)
EVALUATION_SEED = 73101
QUERY_SIZE = 16
DRAWS = 2
TARGET_MSE = 1.00
OUT = Path(__file__).resolve().parent / "FEASIBILITY.json"


def decompose(prediction: np.ndarray, truth: np.ndarray) -> dict:
    level = float((prediction.mean() - truth.mean()) ** 2)
    centered = float((((prediction - prediction.mean())
                       - (truth - truth.mean())) ** 2).mean())
    return {"mse": float(((prediction - truth) ** 2).mean()),
            "level_squared": level, "centered_mse": centered}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    arguments = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    specs = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, QUERY_SIZE, DRAWS, EVALUATION_SEED, None)
    banks = {k: tuple(compact_episode(normalized_episode(
        data.materialize(spec), label_scale)) for spec in items)
        for k, items in specs.items()}

    blob = torch.load(arguments.checkpoint, map_location="cpu",
                      weights_only=False)
    adaptation = StageBAdaptation.from_dict(blob["adaptation"])
    model, _, _ = load_arm(arguments.checkpoint, data, arguments.device)
    model.eval()

    rows: list[dict] = []
    scale, mean = label_scale.scale, label_scale.mean
    for k, bank in banks.items():
        for episode in bank:
            spec = episode.spec
            parts = align_atoms(episode_tensors(model, episode,
                                                arguments.device, torch.float32))
            truth = parts["query_y"].squeeze(0).cpu().numpy() * scale + mean
            with torch.no_grad():
                task = encode_parts(model, parts)
                prediction = predict(model, parts, task, adaptation)["prediction"]
            values = prediction.squeeze(0).cpu().numpy() * scale + mean
            actual = decompose(values, truth)
            # Oracles. Each replaces exactly one half of the prediction.
            oracle_level = values - values.mean() + truth.mean()
            oracle_shape = (truth - truth.mean()) + values.mean()
            constant = np.full_like(truth, values.mean())
            perfect_constant = np.full_like(truth, truth.mean())
            rows.append({
                "k": k, "component": spec.component, "target": spec.target,
                "n": int(len(truth)),
                **{f"actual_{key}": value for key, value in actual.items()},
                "oracle_level_mse": float(((oracle_level - truth) ** 2).mean()),
                "oracle_shape_mse": float(((oracle_shape - truth) ** 2).mean()),
                "constant_mse": float(((constant - truth) ** 2).mean()),
                "within_target_variance": float(((truth - truth.mean()) ** 2).mean()),
                "truth_mean": float(truth.mean()),
                "prediction_mean": float(values.mean()),
            })

    def summarize(field: str, k: int) -> float:
        return component_target_mean(
            [(r["component"], r["target"], r[field]) for r in rows if r["k"] == k])

    payload = {
        "schema": "MetaSieve.StageC.Feasibility.v1", "date": "2026-08-17",
        "checkpoint": str(arguments.checkpoint.resolve().relative_to(ROOT)),
        "arm_mode": adaptation.mode,
        "target_mse_pk": TARGET_MSE,
        "population": {"split": "meta_val", "episodes_per_k": len(banks[0]),
                       "targets": len({e.spec.target for e in banks[0]}),
                       "components": len({e.spec.component for e in banks[0]})},
        "by_k": {}, "meta_test": data.seal_record(),
    }
    print(f"{'k':>2} {'MSE':>8} {'level^2':>8} {'shape':>8} "
          f"{'oracleLvl':>10} {'oracleShp':>10} {'within-var':>11} {'const':>8}")
    for k in SUPPORT_SIZES:
        block = {field: summarize(field, k) for field in (
            "actual_mse", "actual_level_squared", "actual_centered_mse",
            "oracle_level_mse", "oracle_shape_mse", "constant_mse",
            "within_target_variance")}
        block["level_share"] = (block["actual_level_squared"]
                                / block["actual_mse"])
        block["reaches_target_with_perfect_level"] = bool(
            block["oracle_level_mse"] <= TARGET_MSE)
        block["reaches_target_with_perfect_shape"] = bool(
            block["oracle_shape_mse"] <= TARGET_MSE)
        block["level_budget_if_shape_unchanged"] = (
            TARGET_MSE - block["actual_centered_mse"])
        payload["by_k"][str(k)] = block
        print(f"{k:>2} {block['actual_mse']:>8.4f} "
              f"{block['actual_level_squared']:>8.4f} "
              f"{block['actual_centered_mse']:>8.4f} "
              f"{block['oracle_level_mse']:>10.4f} "
              f"{block['oracle_shape_mse']:>10.4f} "
              f"{block['within_target_variance']:>11.4f} "
              f"{block['constant_mse']:>8.4f}")

    # The decisive question: can the target be reached at all by fixing level?
    zero = payload["by_k"]["0"]
    payload["verdict"] = {
        "k0_level_share": zero["level_share"],
        "k0_perfect_level_mse": zero["oracle_level_mse"],
        "k0_target_reachable_by_level_alone": zero["oracle_level_mse"] <= TARGET_MSE,
        "k0_shape_budget_needed": (
            None if zero["oracle_level_mse"] <= TARGET_MSE
            else zero["actual_centered_mse"] - TARGET_MSE),
        "statement": (
            "with a PERFECT target-level predictor the k=0 MSE would be "
            f"{zero['oracle_level_mse']:.4f}; the target is "
            f"{'reachable' if zero['oracle_level_mse'] <= TARGET_MSE else 'NOT reachable'} "
            "by level calibration alone"),
    }
    print("\n" + payload["verdict"]["statement"])
    ratio = zero["actual_centered_mse"] / zero["within_target_variance"]
    payload["verdict"]["k0_centered_mse_over_within_variance"] = ratio
    print(f"k=0 centered MSE / within-target variance = {ratio:.4f} "
          f"(1.0 means no better than predicting the target's own mean)")

    OUT.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
