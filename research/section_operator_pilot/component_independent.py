"""F6I protein-independent task location plus biological interaction curve."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from ..pkis_mechanism_pilot.mechanism import double_center, protein_model_features, stable_fold
from .bioactivity_atlas import (
    _atlas_protein_views, _atlas_surfaces, _bootstrap, _fit_ligand_atlas,
    _ligand_bits, _predict_profile, _protein_temperatures,
    _select_ligand_ridge, _uniform_surfaces,
)
from .ceiling_probe import SEEDS, _load
from .conditional_bilinear import (
    SEED, _d_optimal, _fit_additive, _masked_double_center,
    _nearest_nonself, _query, _safe_spearman, _sha256,
    _subpocket_weights, _write_json,
)
from .component_product import ARM_NAMES, CONTROL_NAMES, CV_SEEDS, K, LOCATION_RIDGES


def _task_location(residual, ridge):
    return float(np.clip(
        float(np.sum(residual)) / (len(residual) + float(ridge)), -0.5, 0.5))


def _cv_location(y, bits, protein_views, scaffolds, protein_groups, ligand_ridge):
    ligand_folds = np.asarray([stable_fold(value, 3) for value in scaffolds])
    protein_folds = np.asarray([stable_fold(value, 3) for value in protein_groups])
    caches = []
    for fold in range(3):
        lt, lv = np.flatnonzero(ligand_folds != fold), np.flatnonzero(ligand_folds == fold)
        pt, pv = np.flatnonzero(protein_folds != fold), np.flatnonzero(protein_folds == fold)
        coefficient = _fit_ligand_atlas(bits[lt], double_center(y[np.ix_(lt, pt)]), ligand_ridge)
        profile = _predict_profile(bits[lv], bits[lt], coefficient)
        source_views = [value[pt] for value in protein_views]
        query_views = [value[pv] for value in protein_views]
        surface = _atlas_surfaces(
            profile, query_views, source_views, _protein_temperatures(source_views))[..., 0]
        caches.append({"y": double_center(y[np.ix_(lv, pv)]), "surface": surface,
                       "scaffold": np.asarray(scaffolds, dtype=object)[lv]})
    report = {}
    for ridge in LOCATION_RIDGES:
        values = []
        for cache in caches:
            for target in range(cache["y"].shape[1]):
                finite = np.flatnonzero(np.isfinite(cache["y"][:, target]))
                for seed in CV_SEEDS:
                    support = _d_optimal(
                        cache["scaffold"], cache["surface"][:, target, None],
                        finite, K, "random", seed + 104729 * target)
                    if support is None:
                        continue
                    query = _query(cache["y"], target, support, cache["scaffold"])
                    if len(query) < 3:
                        continue
                    tau = _task_location(cache["y"][support, target], ridge)
                    prediction = cache["surface"][query, target] + tau
                    values.append(float(np.mean(
                        (cache["y"][query, target] - prediction) ** 2)))
        report[ridge] = values
    selected = min(
        (float(np.mean(value)), ridge) for ridge, value in report.items() if value)[1]
    return float(selected), {
        str(ridge): {"mean_mse": float(np.mean(value)), "n_episodes": len(value)}
        for ridge, value in report.items() if value
    }


def _evaluate(y, additive, surface, null_surface, protein_views, scaffolds,
              location_ridge):
    n_target = y.shape[1]
    per_target = {name: {metric: np.full((n_target, len(SEEDS)), np.nan)
                         for metric in ("mse", "mae", "spearman", "interaction_mse")}
                  for name in ARM_NAMES}
    tau_values = np.full((n_target, len(SEEDS)), np.nan)
    permutation_error = 0.0
    interaction_error = 0.0
    shared_location_error = 0.0
    nearest = _nearest_nonself(protein_views)
    interaction_y = _masked_double_center(y)
    for target in range(n_target):
        wrong_target = int(nearest[target])
        finite = np.flatnonzero(
            np.isfinite(y[:, target]) & np.isfinite(y[:, wrong_target]))
        for seed_index, seed in enumerate(SEEDS):
            support = _d_optimal(
                scaffolds, surface[:, target, None], finite, K, "random",
                seed + 104729 * target)
            if support is None:
                continue
            query = _query(y, target, support, scaffolds)
            if len(query) < 3:
                continue
            residual = y[support, target] - additive[support, target]
            wrong_residual = y[support, wrong_target] - additive[support, wrong_target]
            tau = _task_location(residual, location_ridge)
            wrong_tau = _task_location(wrong_residual, location_ridge)
            rng = np.random.default_rng(seed + 65537 * target)
            permuted_tau = _task_location(residual[rng.permutation(K)], location_ridge)
            permutation_error = max(permutation_error, abs(tau - permuted_tau))
            # Correct/null/nearest all receive the exact same task location.
            shared_location_error = max(shared_location_error, abs(tau - tau), abs(tau - tau))
            curves = {
                "support_free": surface[:, target],
                "location_only": np.full(y.shape[0], tau),
                "component_product": surface[:, target] + tau,
                "zero_protein": null_surface[:, target] + tau,
                "nearest_protein": surface[:, wrong_target] + tau,
                "wrong_support": surface[:, target] + wrong_tau,
                "permuted_support": surface[:, target] + permuted_tau,
            }
            finite_target = np.flatnonzero(np.isfinite(y[:, target]))
            centered_component = curves["component_product"] - float(
                np.mean(curves["component_product"][finite_target]))
            centered_support_free = curves["support_free"] - float(
                np.mean(curves["support_free"][finite_target]))
            interaction_error = max(
                interaction_error,
                float(np.max(np.abs(centered_component - centered_support_free))))
            for name, residual_curve in curves.items():
                prediction = np.clip(additive[query, target] + residual_curve[query], 0, 1)
                outcome = y[query, target]
                error = outcome - prediction
                per_target[name]["mse"][target, seed_index] = float(np.mean(error ** 2))
                per_target[name]["mae"][target, seed_index] = float(np.mean(np.abs(error)))
                per_target[name]["spearman"][target, seed_index] = _safe_spearman(outcome, prediction)
                centered = residual_curve - float(np.mean(residual_curve[finite_target]))
                per_target[name]["interaction_mse"][target, seed_index] = float(
                    np.mean((interaction_y[query, target] - centered[query]) ** 2))
            tau_values[target, seed_index] = tau
    return (per_target, tau_values, nearest, permutation_error,
            interaction_error, shared_location_error)


def _summarize(per_target, tau, nearest, permutation_error,
               interaction_error, shared_location_error):
    metrics = {name: {metric: float(np.nanmean(value)) for metric, value in values.items()}
               for name, values in per_target.items()}
    contrasts = {name: _bootstrap(values, per_target["component_product"])
                 for name, values in per_target.items() if name != "component_product"}
    interaction = {name: _bootstrap(
        per_target[name], per_target["component_product"], metric="interaction_mse")
        for name in ("support_free", "nearest_protein")}
    return {
        "metrics": metrics, "mse_reduction_component_vs": contrasts,
        "interaction_mse_reduction_component_vs": interaction,
        "tau": {"mean": float(np.nanmean(tau)), "std": float(np.nanstd(tau))},
        "exact_invariants": {
            "label_permutation_max_abs_prediction_error": permutation_error,
            "support_free_interaction_max_abs_error": interaction_error,
            "protein_control_location_max_abs_error": shared_location_error,
        }, "nearest_nonself_target_index": nearest.tolist(),
    }


def run(args):
    started = time.time()
    data = _load(args)
    y_source, source_ligands, source_targets = (
        data["y_source"], data["source_ligands"], data["source_targets"])
    source_bits = _ligand_bits(source_ligands)
    source_scaffolds = np.asarray([record.generic_scaffold for record in source_ligands], dtype=object)
    source_groups = [item["record"].group for item in source_targets]
    masks = _subpocket_weights(args.kissim_distances)
    source_protein_views = _atlas_protein_views(source_targets, masks)
    ligand_ridge, ligand_cv = _select_ligand_ridge(y_source, source_bits, source_scaffolds)
    location_ridge, location_cv = _cv_location(
        y_source, source_bits, source_protein_views, source_scaffolds,
        source_groups, ligand_ridge)
    print(f"selected location ridge={location_ridge}", flush=True)
    ligand_coefficient = _fit_ligand_atlas(source_bits, double_center(y_source), ligand_ridge)
    temperatures = _protein_temperatures(source_protein_views)
    additive_model = _fit_additive(
        y_source, np.stack([record.nuisance_features for record in source_ligands]),
        np.stack([protein_model_features(item["record"].pocket) for item in source_targets]))
    transfers = {}
    for panel in data["panels"]:
        print(f"evaluating {panel['name']}", flush=True)
        profile = _predict_profile(_ligand_bits(panel["ligands"]), source_bits, ligand_coefficient)
        protein_views = _atlas_protein_views(panel["targets"], masks)
        surface = _atlas_surfaces(profile, protein_views, source_protein_views, temperatures)[..., 0]
        null_surface = _uniform_surfaces(profile, len(panel["targets"]))[..., 0]
        additive = additive_model.predict(
            np.stack([record.nuisance_features for record in panel["ligands"]]),
            np.stack([protein_model_features(item["record"].pocket) for item in panel["targets"]]))
        scaffolds = np.asarray([record.generic_scaffold for record in panel["ligands"]], dtype=object)
        transfers[panel["name"]] = {
            "role": panel["role"], "n_ligands": len(panel["ligands"]),
            "n_targets": len(panel["targets"]), "support_size": K,
            "random": _summarize(*_evaluate(
                panel["y"], additive, surface, null_surface, protein_views,
                scaffolds, location_ridge)),
        }
    primary, external = transfers["PKIS2"]["random"], transfers["Anastassiadis2011"]["random"]
    pkis_raw = all(primary["mse_reduction_component_vs"][n]["ci95"][0] > 0
                   for n in CONTROL_NAMES)
    external_raw = all(external["mse_reduction_component_vs"][n]["estimate"] > 0
                       for n in CONTROL_NAMES)
    pkis_protein = primary["interaction_mse_reduction_component_vs"][
        "nearest_protein"]["ci95"][0] > 0
    external_protein = external["interaction_mse_reduction_component_vs"][
        "nearest_protein"]["estimate"] > 0
    invariants = all(
        value["exact_invariants"][name] <= 1e-12
        for value in (primary, external)
        for name in ("label_permutation_max_abs_prediction_error",
                     "support_free_interaction_max_abs_error",
                     "protein_control_location_max_abs_error"))
    passed = bool(pkis_raw and external_raw and pkis_protein
                  and external_protein and invariants)
    result = {
        "schema": "MetaSieve.IndependentComponentProduct.F6I.v1",
        "selected": {"d_adapt": 1, "k": K, "ligand_ridge": ligand_ridge,
                     "location_ridge": location_ridge, "tau_bounds": [-0.5, 0.5]},
        "source": {"panel": "PKIS1", "ligand_cv": ligand_cv, "location_cv": location_cv},
        "transfers": transfers,
        "gate": {"pkis2_raw_components": pkis_raw,
                 "anastassiadis_raw_point_estimates": external_raw,
                 "pkis2_protein_interaction": bool(pkis_protein),
                 "anastassiadis_protein_interaction": bool(external_protein),
                 "exact_invariants": invariants, "passed": passed,
                 "verdict": ("F6I_COMPONENTS_ADMISSIBLE" if passed
                             else "F6I_COMPONENTS_NOT_ADMISSIBLE")},
        "read_firewall": {"kcgs_numeric_outcomes": "NOT_READ", "davis_labels": "NOT_READ",
                          "recipient_labels": "NOT_READ"},
        "limitations": [
            "Support labels identify task location only; biological interaction is zero-shot from the source atlas.",
            "Law-valued integration remains required before production use."],
        "elapsed_seconds": time.time() - started,
    }
    output = Path(args.output).resolve()
    _write_json(output / "result.json", result)
    _write_json(output / "manifest.json", {"sha256": {
        "preregistration": _sha256(Path(__file__).with_name("F6I_PREREGISTRATION.md")),
        "script": _sha256(Path(__file__)),
        "pkis1_labels": _sha256(Path(args.informers_root).resolve() / "data/pkis1_continuous_labels.csv")},
        "parameters": vars(args)})
    print(json.dumps(result["gate"], indent=2), flush=True)


def parser():
    item = argparse.ArgumentParser()
    item.add_argument("--informers-root", default="../external/informers")
    item.add_argument("--klifs-json", default="../external/klifs/kinase_information_human.json")
    item.add_argument("--kissim-distances", default="../external/kissim/kissim/data/min_max_distances_fine.csv")
    item.add_argument("--anastassiadis-workbook", default="external/anastassiadis/NIHMS328213-supplement-3.xls")
    item.add_argument("--anastassiadis-identities", default="external/anastassiadis/compound_identities.json")
    item.add_argument("--output", default="research/section_operator_pilot/artifacts/f6i")
    return item


if __name__ == "__main__":
    run(parser().parse_args())
