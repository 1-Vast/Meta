"""R14 diagnostic 1: is the within-target shape error ordering or amplitude?

No training. Loads frozen checkpoints, runs the identical fixed `meta_val`
bank, and decomposes each target's within-target shape term exactly.

For centered predictions `p` and centered labels `y` inside one target:

    shape = E[(p - y)^2]
          = Var(y) (1 - r^2)  +  (sd_p - r sd_y)^2
            \___ordering___/     \____amplitude____/

with `r = corr(p, y)`. The split matters because the two halves need
different mechanisms and have different costs:

* the **ordering** half cannot be reduced by any rescaling — it needs a
  genuinely better within-target readout;
* the **amplitude** half is removed by multiplying the centered prediction by
  a scalar. A positive scalar is a monotone transform, so it changes MSE and
  **cannot change CI, Spearman or any sign/ranking metric at all**.

If the amplitude half is large, then part of the k=0 MSE frontier is
reachable without paying anything on the ranking side — which is precisely
the Pareto conflict the R3R4-R13 ladder kept running into. If it is small,
Innovation A's protein-conditioned amplitude is not worth building and the
cycle must attack ordering directly.

Three scale estimators are reported, and they are **not** interchangeable:

``oracle_per_target``
    the per-target minimiser. Transductive — it reads the target's own query
    labels. A diagnostic upper bound, never deployable.
``global_in_sample``
    one scalar fitted on all of `meta_val` and applied to `meta_val`.
    Development-grade: selection and inference share the population.
``global_loco``
    one scalar fitted on 18 homology components and applied to the held-out
    19th, rotated. This is the honest, deployable-shaped estimate, and it is
    the number any claim must use.

`meta_test` is never read.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np
import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.qpsmp_data import QPSMPData
from scripts.stageR0_retrieval_falsification import (
    component_bootstrap, component_target_mean,
)
from scripts.stageR6_compare_arms import SUPPORT_SIZES, load_arm
from scripts.stageR9_pair_audit import predictions_for
from scripts.train_level_shape import matched_donors, normalized
from scripts.train_qpsmp import (
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, training_label_scale,
)


def concordance(prediction: np.ndarray, truth: np.ndarray) -> float:
    """Fraction of comparable pairs ordered correctly (ties count a half)."""
    rows, cols = np.triu_indices(len(truth), 1)
    delta_y = truth[rows] - truth[cols]
    comparable = delta_y != 0
    if not comparable.any():
        return float("nan")
    signed = np.sign(delta_y[comparable]) * (prediction[rows] - prediction[cols])[comparable]
    return float((signed > 0).mean() + 0.5 * (signed == 0).mean())


def decompose(prediction: np.ndarray, truth: np.ndarray) -> dict:
    """Exact ordering/amplitude split of one target's shape term."""
    centered_p = prediction - prediction.mean()
    centered_y = truth - truth.mean()
    var_p = float((centered_p ** 2).mean())
    var_y = float((centered_y ** 2).mean())
    covariance = float((centered_p * centered_y).mean())
    sd_p, sd_y = float(np.sqrt(var_p)), float(np.sqrt(var_y))
    r = covariance / (sd_p * sd_y) if sd_p > 0 and sd_y > 0 else 0.0
    shape = float(((centered_p - centered_y) ** 2).mean())
    ordering = var_y * (1.0 - r ** 2)
    amplitude = (sd_p - r * sd_y) ** 2
    return {
        "shape_pk": shape,
        "calibration_pk": float((prediction.mean() - truth.mean()) ** 2),
        "var_truth": var_y, "var_pred": var_p,
        "sd_truth": sd_y, "sd_pred": sd_p,
        "pearson_r": r,
        "shape_ordering_floor": ordering,
        "shape_amplitude_excess": amplitude,
        # the multiplier on the centered prediction that minimises shape
        "optimal_scale": covariance / var_p if var_p > 0 else float("nan"),
        # sufficient statistics for pooling a *single* scale over targets
        "cov_sum": covariance * len(truth), "var_pred_sum": var_p * len(truth),
        "n": int(len(truth)),
        "ci": concordance(prediction, truth),
    }


def rescaled_shape(row: dict, scale: float) -> float:
    """One target's shape when its centered prediction is multiplied by `scale`."""
    return (row["var_truth"]
            - 2.0 * scale * row["pearson_r"] * row["sd_pred"] * row["sd_truth"]
            + scale ** 2 * row["var_pred"])


def weighted(rows: list[dict], field: str) -> float:
    """Equal-component, then equal-target mean — the project's reporting protocol."""
    return component_target_mean(
        (row["component"], row["target"], row.get(field)) for row in rows)


def shape_at_scale(rows: list[dict], scale: float) -> float:
    """Equal-component mean shape after a single global rescale."""
    return component_target_mean(
        (row["component"], row["target"], rescaled_shape(row, scale))
        for row in rows)


def pooled_scale(rows: list[dict]) -> float:
    """The single scalar minimising total shape over a set of targets."""
    numerator = sum(row["cov_sum"] for row in rows)
    denominator = sum(row["var_pred_sum"] for row in rows)
    return numerator / denominator if denominator > 0 else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", action="append", required=True,
                        help="name=path/to/checkpoint.pt (repeat to add seeds)")
    parser.add_argument("--split-directory", type=Path, required=True)
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--query-size", type=int, default=16)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=args.split_directory,
                     include_meta_test=False)
    scale = training_label_scale(data)
    donors = matched_donors(data, "meta_val", donor_pool="meta_val",
                            whitening_pool="meta_train")
    banks = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, args.query_size, 1, args.evaluation_seed, None)

    models: dict[str, list] = defaultdict(list)
    for item in args.arm:
        name, _, path = item.partition("=")
        models[name].append(load_arm(Path(path), data, args.device))

    # rows[arm][k] = list of per-(seed, target) decompositions
    rows: dict[str, dict[int, list]] = defaultdict(lambda: defaultdict(list))
    for name, entries in models.items():
        for model, kind, seed in entries:
            dtype = next(model.parameters()).dtype
            with torch.no_grad():
                for k, specs in banks.items():
                    for spec in specs:
                        episode = compact_episode(
                            normalized(data.materialize(spec), scale))
                        prediction = predictions_for(
                            model, kind, data, spec, episode, donors[spec.target],
                            args.device, dtype) * scale.scale + scale.mean
                        truth = episode.query_y.numpy() * scale.scale + scale.mean
                        row = decompose(prediction, truth)
                        row.update(target=spec.target, component=spec.component,
                                   seed=seed)
                        rows[name][k].append(row)
            del model
            torch.cuda.empty_cache()

    audit: dict = {
        "schema": "MetaSieve.R14DispersionAudit.v1",
        "split": "meta_val",
        "split_directory": str(args.split_directory),
        "split_assignment_sha256": data.split_manifest["assignment_sha256"],
        "meta_test": data.seal_record(),
        "population": {"targets": len(banks[0]), "query_size": args.query_size,
                       "evaluation_seed": args.evaluation_seed},
        "arms": {},
    }

    for name in sorted(rows):
        arm: dict = {}
        for k in SUPPORT_SIZES:
            block = rows[name][k]
            if not block:
                continue
            components = sorted({row["component"] for row in block})

            for row in block:
                row["dispersion_ratio"] = (
                    row["sd_pred"] / row["sd_truth"] if row["sd_truth"] > 0 else None)

            # Leave-one-component-out: fit the single scale on the other
            # components, apply it to the held-out one. Honest, deployable.
            loco_rows, loco_scales, gains = [], [], []
            for held in components:
                scale_value = pooled_scale(
                    [row for row in block if row["component"] != held])
                loco_scales.append(scale_value)
                for row in block:
                    if row["component"] != held:
                        continue
                    after = rescaled_shape(row, scale_value)
                    loco_rows.append(dict(row, loco_shape=after))
                    gains.append((row["component"], row["target"],
                                  row["shape_pk"] - after))

            in_sample = pooled_scale(block)
            base_shape = weighted(block, "shape_pk")
            base_calibration = weighted(block, "calibration_pk")
            oracle_shape = weighted(block, "shape_ordering_floor")
            loco_shape = weighted(loco_rows, "loco_shape")

            arm[str(k)] = {
                "mse_pk": base_calibration + base_shape,
                "calibration_pk": base_calibration,
                "shape_pk": base_shape,
                "shape_ordering_floor": oracle_shape,
                "shape_amplitude_excess": base_shape - oracle_shape,
                "amplitude_share_of_shape": (
                    (base_shape - oracle_shape) / base_shape if base_shape else None),
                "mean_pearson_r": weighted(block, "pearson_r"),
                "mean_sd_pred": weighted(block, "sd_pred"),
                "mean_sd_truth": weighted(block, "sd_truth"),
                "mean_dispersion_ratio": weighted(block, "dispersion_ratio"),
                "mean_ci": weighted(block, "ci"),
                "scale": {
                    "oracle_per_target_shape": oracle_shape,
                    "global_in_sample": in_sample,
                    "global_in_sample_shape": shape_at_scale(block, in_sample),
                    "global_loco_mean": float(np.mean(loco_scales)),
                    "global_loco_shape": loco_shape,
                    "global_loco_mse_pk": base_calibration + loco_shape,
                },
                "loco_mse_gain_component_bootstrap": component_bootstrap(
                    gains, args.bootstrap_draws, 20260816),
                "note": "a positive rescale is monotone within a target, so it "
                        "changes MSE and cannot change CI/Spearman/sign accuracy",
            }
        audit["arms"][name] = arm

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, indent=1), encoding="utf-8")
    rows_path = args.output.with_suffix(".rows.jsonl")
    with rows_path.open("w", encoding="utf-8") as handle:
        for name in sorted(rows):
            for k, block in rows[name].items():
                for row in block:
                    handle.write(json.dumps(dict(row, arm=name, k=k)) + "\n")
    print(f"wrote {args.output} and {rows_path}")
    for name, arm in audit["arms"].items():
        for k, cell in arm.items():
            print(f"{name} k={k}: MSE {cell['mse_pk']:.4f} "
                  f"(calib {cell['calibration_pk']:.4f} + shape {cell['shape_pk']:.4f}) | "
                  f"r {cell['mean_pearson_r']:.3f} | "
                  f"sd_pred/sd_truth {cell['mean_dispersion_ratio']:.3f} | "
                  f"amplitude share {cell['amplitude_share_of_shape']:.3f} | "
                  f"LOCO scale {cell['scale']['global_loco_mean']:.3f} -> "
                  f"MSE {cell['scale']['global_loco_mse_pk']:.4f}")


if __name__ == "__main__":
    raise SystemExit(main())
