"""Stage T2 analysis: strata, two-way clustered intervals, and the frozen gate.

Rows sharing a protein component or a transformation key are correlated, so an
independent row bootstrap would be anticonservative by a large factor. The
preregistered estimator is a **two-way (multiway) cluster bootstrap**: each draw
samples protein components with replacement and transformation keys with
replacement, independently, and a row's multiplicity is the product of its two
draw counts. Both arms of a contrast are re-scored on the identical draw.

Effective independent units are `min(#components, #keys)` and are printed next
to every interval, because on this corpus that number is small and the intervals
must be read with it in view.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy import stats

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

HERE = Path(__file__).resolve().parent

# Frozen in PREREGISTRATION.md section 7.
GATES = {
    "pearson_gain_threshold": 0.05,
    "bootstrap_draws": 2000,
    "bootstrap_seed": 20260820,
    "max_single_transformation_share_of_effect": 0.5,
    "max_single_component_share_of_effect": 0.5,
}
CANDIDATE = "C_protein"
PRIMARY = "internal_all_correct"


class Rows:
    __slots__ = ("row_id", "truth", "prediction", "key", "component_left",
                 "component_right", "activity_cliff", "cross_component")

    def __init__(self, payload: dict) -> None:
        self.row_id = list(payload["row_id"])
        self.truth = np.asarray(payload["truth"], dtype=np.float64)
        self.prediction = np.asarray(payload["prediction"], dtype=np.float64)
        self.key = list(payload["key"])
        self.component_left = list(payload["component_left"])
        self.component_right = list(payload["component_right"])
        self.activity_cliff = np.asarray(payload["activity_cliff"], dtype=bool)
        self.cross_component = np.asarray(payload["cross_component"], dtype=bool)

    def __len__(self) -> int:
        return int(self.truth.size)

    def select(self, mask) -> "Rows":
        keep = np.flatnonzero(np.asarray(mask))
        return Rows({
            "row_id": [self.row_id[int(i)] for i in keep],
            "truth": self.truth[keep], "prediction": self.prediction[keep],
            "key": [self.key[int(i)] for i in keep],
            "component_left": [self.component_left[int(i)] for i in keep],
            "component_right": [self.component_right[int(i)] for i in keep],
            "activity_cliff": self.activity_cliff[keep],
            "cross_component": self.cross_component[keep],
        })


def _metrics(truth: np.ndarray, prediction: np.ndarray,
             weights: np.ndarray | None = None) -> dict:
    if truth.size < 3:
        return {name: float("nan") for name in
                ("mse", "mae", "pearson", "spearman", "ci", "sign_accuracy")}
    if weights is None:
        weights = np.ones_like(truth)
    total = float(weights.sum())
    error = (prediction - truth) ** 2
    mse = float((error * weights).sum() / total)
    mae = float((np.abs(prediction - truth) * weights).sum() / total)
    if np.std(truth) < 1e-12 or np.std(prediction) < 1e-12:
        pearson = float("nan")
    else:
        mean_t = float((truth * weights).sum() / total)
        mean_p = float((prediction * weights).sum() / total)
        cov = float((weights * (truth - mean_t) * (prediction - mean_p)).sum() / total)
        var_t = float((weights * (truth - mean_t) ** 2).sum() / total)
        var_p = float((weights * (prediction - mean_p) ** 2).sum() / total)
        pearson = cov / max(np.sqrt(var_t * var_p), 1e-12)
    expanded_t, expanded_p = truth, prediction
    if not np.allclose(weights, 1.0):
        counts = np.maximum(weights.astype(np.int64), 0)
        expanded_t = np.repeat(truth, counts)
        expanded_p = np.repeat(prediction, counts)
    spearman = stats.spearmanr(expanded_t, expanded_p).statistic \
        if expanded_t.size >= 3 else float("nan")
    tau = stats.kendalltau(expanded_t, expanded_p).statistic \
        if expanded_t.size >= 3 else float("nan")
    usable = np.abs(truth) > 0
    sign = (float((weights[usable] * (np.sign(truth[usable])
                                      == np.sign(prediction[usable]))).sum()
                  / max(float(weights[usable].sum()), 1e-12))
            if usable.any() else float("nan"))
    return {
        "mse": mse, "mae": mae,
        "pearson": float(pearson) if np.isfinite(pearson) else float("nan"),
        "spearman": float(spearman) if np.isfinite(spearman) else float("nan"),
        "ci": float((tau + 1.0) / 2.0) if np.isfinite(tau) else float("nan"),
        "sign_accuracy": sign,
    }


def metrics_of(rows: Rows) -> dict:
    out = _metrics(rows.truth, rows.prediction)
    out["n"] = len(rows)
    out["components"] = len(set(rows.component_left) | set(rows.component_right))
    out["keys"] = len(set(rows.key))
    out["effective_independent_units"] = min(out["components"], out["keys"])
    return out


STATISTICS = ("mse", "mae", "pearson", "spearman", "ci", "sign_accuracy")


def two_way_cluster_bootstrap(left: Rows, right: Rows, draws: int,
                              seed: int) -> dict:
    """Multiway cluster bootstrap over protein components and transformation keys."""
    if left.row_id != right.row_id:
        raise ValueError("paired bootstrap requires identical evaluation rows")
    components = sorted(set(left.component_left) | set(left.component_right))
    keys = sorted(set(left.key))
    component_index = {name: i for i, name in enumerate(components)}
    key_index = {name: i for i, name in enumerate(keys)}
    # A row belongs to two components; its component multiplicity is the mean of
    # the two draw counts, which keeps a cross-component row from being counted
    # under only one side.
    row_left = np.asarray([component_index[c] for c in left.component_left])
    row_right = np.asarray([component_index[c] for c in left.component_right])
    row_key = np.asarray([key_index[k] for k in left.key])

    rng = np.random.default_rng(seed)
    samples = {name: np.full(draws, np.nan) for name in STATISTICS}
    for draw in range(draws):
        component_counts = np.bincount(
            rng.integers(0, len(components), size=len(components)),
            minlength=len(components))
        key_counts = np.bincount(
            rng.integers(0, len(keys), size=len(keys)), minlength=len(keys))
        weight = (0.5 * (component_counts[row_left] + component_counts[row_right])
                  * key_counts[row_key]).astype(np.float64)
        if weight.sum() <= 0:
            continue
        keep = weight > 0
        if keep.sum() < 3:
            continue
        left_metrics = _metrics(left.truth[keep], left.prediction[keep],
                                weight[keep])
        right_metrics = _metrics(right.truth[keep], right.prediction[keep],
                                 weight[keep])
        for name in STATISTICS:
            samples[name][draw] = left_metrics[name] - right_metrics[name]

    point_left, point_right = metrics_of(left), metrics_of(right)
    out = {
        "n": len(left),
        "components": len(components),
        "keys": len(keys),
        "effective_independent_units": min(len(components), len(keys)),
        "draws": draws,
    }
    for name in STATISTICS:
        values = samples[name][np.isfinite(samples[name])]
        lo, hi = ((float(np.quantile(values, 0.025)),
                   float(np.quantile(values, 0.975))) if values.size
                  else (float("nan"), float("nan")))
        out[name] = {
            "left": point_left[name], "right": point_right[name],
            "delta": point_left[name] - point_right[name],
            "lo": lo, "hi": hi,
            "resolved": bool(np.isfinite(lo) and np.isfinite(hi)
                             and (lo > 0.0 or hi < 0.0)),
            "effective_draws": int(values.size),
        }
    return out


def leave_one_out_influence(rows: Rows, baseline: Rows, attribute: str) -> dict:
    """Gate 7: is the effect carried by one transformation or one component?"""
    groups: dict[str, list[int]] = defaultdict(list)
    for position in range(len(rows)):
        if attribute == "key":
            groups[rows.key[position]].append(position)
        else:
            groups[rows.component_left[position]].append(position)
            groups[rows.component_right[position]].append(position)
    full = metrics_of(rows)["pearson"] - metrics_of(baseline)["pearson"]
    worst_name, worst_value = None, 0.0
    for name, positions in groups.items():
        mask = np.ones(len(rows), dtype=bool)
        mask[np.asarray(sorted(set(positions)))] = False
        if mask.sum() < 10:
            continue
        without = (metrics_of(rows.select(mask))["pearson"]
                   - metrics_of(baseline.select(mask))["pearson"])
        drop = full - without
        if abs(drop) > abs(worst_value):
            worst_name, worst_value = name, drop
    return {
        "full_effect": full,
        "most_influential": worst_name,
        "effect_change_when_removed": worst_value,
        "share_of_effect": (abs(worst_value) / abs(full)) if full else float("nan"),
        "groups": len(groups),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, default=HERE / "runs")
    parser.add_argument("--output", type=Path, default=HERE / "T2_RESULT.json")
    parser.add_argument("--draws", type=int, default=GATES["bootstrap_draws"])
    args = parser.parse_args()

    arms = sorted(path.name for path in args.runs.iterdir() if path.is_dir())
    runs: dict[str, dict] = {}
    rows: dict[str, dict[str, Rows]] = {}
    for arm in arms:
        runs[arm] = json.loads((args.runs / arm / "RUN.json").read_text(
            encoding="utf-8"))
        rows[arm] = {}
        for path in sorted((args.runs / arm).glob("*.rows.json")):
            rows[arm][path.name.replace(".rows.json", "")] = Rows(
                json.loads(path.read_text(encoding="utf-8")))

    report: dict = {
        "schema": "MetaSieve.StageT.T2Result.v1",
        "stage": "stageT_mmp",
        "preregistration_sha256": hashlib.sha256(
            (HERE / "PREREGISTRATION.md").read_bytes()).hexdigest(),
        "gates_frozen": GATES,
        "meta_test": runs[arms[0]]["meta_test"],
        "selection": runs[arms[0]]["selection"],
        "banks": runs[arms[0]]["banks"],
        "arms": {arm: {"config": runs[arm]["arm"],
                       "parameters": runs[arm]["parameters"],
                       "steps": runs[arm]["steps"], "seed": runs[arm]["seed"],
                       "final_training_loss": runs[arm]["history"][-1]["loss"]}
                 for arm in arms},
        "metrics": {arm: {name: metrics_of(rows[arm][name])
                          for name in sorted(rows[arm])} for arm in arms},
    }

    def contrast(left: Rows, right: Rows) -> dict:
        return two_way_cluster_bootstrap(left, right, args.draws,
                                         GATES["bootstrap_seed"])

    candidate = rows[CANDIDATE]
    contrasts: dict[str, dict] = {}
    # Gates 1-2: correct protein vs the trained shuffled-protein arm.
    contrasts["C_vs_D_shuffled_arm"] = contrast(
        candidate[PRIMARY], rows["D_protein_shuffled"][PRIMARY])
    # Same comparison as a paired within-arm substitution: only the protein
    # input changes, the model is identical.
    contrasts["C_correct_vs_C_shuffled_input"] = contrast(
        candidate[PRIMARY], candidate["internal_all_shuffled"])
    # Gates 3-4: correct vs similarity-matched wrong protein, both forms.
    contrasts["C_correct_vs_C_matched_wrong_input"] = contrast(
        candidate[PRIMARY], candidate["internal_all_matched_wrong"])
    contrasts["C_vs_E_matched_wrong_arm"] = contrast(
        candidate[PRIMARY], rows["E_protein_matched_wrong"][PRIMARY])
    # Gate 6: label shuffle.
    contrasts["C_vs_F_label_shuffled"] = contrast(
        candidate[PRIMARY], rows["F_label_shuffled"][PRIMARY])
    # Gate 10: the two evaluation surfaces.
    contrasts["C_vs_D_repeated_keys"] = contrast(
        candidate["internal_repeated_correct"],
        rows["D_protein_shuffled"]["internal_repeated_correct"])
    contrasts["C_vs_D_transformation_disjoint"] = contrast(
        candidate["internal_disjoint_correct"],
        rows["D_protein_shuffled"]["internal_disjoint_correct"])
    # Gate 9: the fit-unsampled bank (target-key shortcut check).
    contrasts["C_vs_D_fit_unsampled"] = contrast(
        candidate["fit_unsampled_correct"],
        rows["D_protein_shuffled"]["fit_unsampled_correct"])
    report["contrasts"] = contrasts

    # -- strata -------------------------------------------------------------
    primary = candidate[PRIMARY]
    shuffled_arm = rows["D_protein_shuffled"][PRIMARY]
    strata = {
        "all": np.ones(len(primary), dtype=bool),
        "cross_component": primary.cross_component,
        "within_component": ~primary.cross_component,
        "activity_cliff": primary.activity_cliff,
        "non_cliff": ~primary.activity_cliff,
    }
    key_counts = defaultdict(int)
    for key in primary.key:
        key_counts[key] += 1
    frequent = np.asarray([key_counts[k] >= 3 for k in primary.key])
    strata["high_frequency_edits"] = frequent
    strata["low_frequency_edits"] = ~frequent
    report["strata"] = {
        name: {
            "n": int(mask.sum()),
            "C_protein": metrics_of(primary.select(mask)),
            "D_protein_shuffled": metrics_of(shuffled_arm.select(mask)),
        }
        for name, mask in strata.items() if int(mask.sum()) >= 10
    }

    # -- gate 7 influence ---------------------------------------------------
    report["influence"] = {
        "by_transformation_key": leave_one_out_influence(
            primary, shuffled_arm, "key"),
        "by_protein_component": leave_one_out_influence(
            primary, shuffled_arm, "component"),
    }

    # -- gate 8 alignment ---------------------------------------------------
    shift = (candidate[PRIMARY].prediction
             - candidate["internal_all_matched_wrong"].prediction)
    truth = candidate[PRIMARY].truth
    report["protein_shift"] = {
        "shift_sd": float(shift.std()),
        "truth_sd": float(truth.std()),
        "shift_over_truth_sd": float(shift.std() / max(truth.std(), 1e-12)),
        "alignment_with_truth": (float(np.corrcoef(shift, truth)[0, 1])
                                 if shift.std() > 1e-12 else float("nan")),
        "interpretation": ("large shift with near-zero alignment is the Stage P "
                           "/ Stage S failure mode: the protein moves the "
                           "prediction without carrying information"),
    }

    # -- the frozen gate ----------------------------------------------------
    threshold = GATES["pearson_gain_threshold"]
    shuffled = contrasts["C_vs_D_shuffled_arm"]["pearson"]
    wrong = contrasts["C_correct_vs_C_matched_wrong_input"]["pearson"]
    label = contrasts["C_vs_F_label_shuffled"]["pearson"]
    mse = contrasts["C_vs_D_shuffled_arm"]["mse"]
    ranking = {name: contrasts["C_vs_D_shuffled_arm"][name]
               for name in ("spearman", "ci", "sign_accuracy")}
    repeated = contrasts["C_vs_D_repeated_keys"]["pearson"]["delta"]
    disjoint = contrasts["C_vs_D_transformation_disjoint"]["pearson"]["delta"]
    fit_shortcut = contrasts["C_vs_D_fit_unsampled"]["pearson"]["delta"]
    influence = report["influence"]
    alignment = report["protein_shift"]["alignment_with_truth"]

    gate = {
        "1_correct_minus_shuffled_pearson": {
            "value": shuffled["delta"], "threshold": threshold,
            "pass": bool(shuffled["delta"] >= threshold)},
        "2_lower_bound_above_zero": {
            "lo": shuffled["lo"], "pass": bool(shuffled["lo"] > 0.0)},
        "3_correct_minus_matched_wrong_pearson": {
            "value": wrong["delta"], "threshold": threshold,
            "pass": bool(wrong["delta"] >= threshold)},
        "4_lower_bound_above_zero": {
            "lo": wrong["lo"], "pass": bool(wrong["lo"] > 0.0)},
        "5_error_and_ranking_both_improve": {
            "mse_delta": mse["delta"],
            "ranking": {name: value["delta"] for name, value in ranking.items()},
            "pass": bool(mse["delta"] < 0.0
                         and all(value["delta"] > 0.0
                                 for value in ranking.values()))},
        "6_label_shuffle_destroys_the_effect": {
            "delta": label["delta"], "lo": label["lo"],
            "pass": bool(label["delta"] >= threshold and label["lo"] > 0.0)},
        "7_not_confined_to_one_transformation_or_component": {
            "key_share": influence["by_transformation_key"]["share_of_effect"],
            "component_share": influence["by_protein_component"]["share_of_effect"],
            "high_frequency_only": bool(
                report["strata"].get("low_frequency_edits", {})
                .get("C_protein", {}).get("pearson", 0) is None),
            "pass": bool(
                abs(influence["by_transformation_key"]["share_of_effect"])
                <= GATES["max_single_transformation_share_of_effect"]
                and abs(influence["by_protein_component"]["share_of_effect"])
                <= GATES["max_single_component_share_of_effect"])},
        "8_protein_shift_aligned_with_truth": {
            "alignment": alignment,
            "shift_over_truth_sd": report["protein_shift"]["shift_over_truth_sd"],
            "pass": bool(np.isfinite(alignment) and alignment >= threshold)},
        "9_no_target_key_shortcut_on_fit_unsampled": {
            "fit_unsampled_gain": fit_shortcut,
            "internal_gain": shuffled["delta"],
            "pass": bool(fit_shortcut <= max(2.0 * shuffled["delta"], threshold))},
        "10_transformation_cold_does_not_reverse": {
            "repeated_keys_gain": repeated,
            "transformation_disjoint_gain": disjoint,
            "pass": bool(np.sign(repeated) == np.sign(disjoint)
                         or disjoint >= 0.0)},
    }
    report["gate"] = gate
    report["verdict"] = {
        "route_passes": bool(all(item["pass"] for item in gate.values())),
        "gates_passed": sum(1 for item in gate.values() if item["pass"]),
        "gates_total": len(gate),
    }

    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "verdict": report["verdict"],
        "gate": {name: item["pass"] for name, item in gate.items()},
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
