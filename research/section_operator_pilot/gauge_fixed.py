"""F4G gauge-fixed location/amplitude section over a biological atlas curve."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from ..pkis_mechanism_pilot.mechanism import (
    double_center,
    protein_model_features,
    stable_fold,
)
from .bioactivity_atlas import (
    _atlas_protein_views,
    _atlas_surfaces,
    _bootstrap,
    _fit_ligand_atlas,
    _ligand_bits,
    _predict_profile,
    _protein_temperatures,
    _select_ligand_ridge,
    _uniform_surfaces,
)
from .ceiling_probe import SEEDS, _load
from .conditional_bilinear import (
    SEED,
    _d_optimal,
    _fit_additive,
    _masked_double_center,
    _nearest_nonself,
    _query,
    _safe_spearman,
    _sha256,
    _subpocket_weights,
    _write_json,
)


K = 5
LOCATION_RIDGES = (1.0, 10.0, 100.0)
AMPLITUDE_RIDGES = (0.1, 1.0, 10.0, 100.0)
CV_SEEDS = SEEDS[:5]
ARM_NAMES = (
    "support_free", "location_only", "gauge_fixed", "zero_protein",
    "nearest_protein", "wrong_support", "permuted_support",
)
CONTROL_NAMES = tuple(name for name in ARM_NAMES if name != "gauge_fixed")


def _solve_gain(surface, residual, location_ridge, amplitude_ridge):
    surface = np.asarray(surface, dtype=np.float64)
    residual = np.asarray(residual, dtype=np.float64)
    design = np.stack([np.ones(len(surface)), surface], axis=1)
    penalty = np.diag([float(location_ridge), float(amplitude_ridge)])
    prior = np.asarray([0.0, 1.0])
    return np.linalg.solve(
        design.T @ design + penalty,
        design.T @ residual + penalty @ prior,
    )


def _curve(surface, coefficient):
    return coefficient[0] + coefficient[1] * np.asarray(surface)


def _cv_penalties(y, bits, protein_views, scaffolds, protein_groups, ligand_ridge):
    ligand_folds = np.asarray([stable_fold(value, 3) for value in scaffolds])
    protein_folds = np.asarray([stable_fold(value, 3) for value in protein_groups])
    caches = []
    for fold in range(3):
        lt, lv = np.flatnonzero(ligand_folds != fold), np.flatnonzero(ligand_folds == fold)
        pt, pv = np.flatnonzero(protein_folds != fold), np.flatnonzero(protein_folds == fold)
        coefficient = _fit_ligand_atlas(
            bits[lt], double_center(y[np.ix_(lt, pt)]), ligand_ridge)
        profile = _predict_profile(bits[lv], bits[lt], coefficient)
        source_views = [value[pt] for value in protein_views]
        query_views = [value[pv] for value in protein_views]
        surface = _atlas_surfaces(
            profile, query_views, source_views,
            _protein_temperatures(source_views))[..., 0]
        caches.append({
            "y": double_center(y[np.ix_(lv, pv)]), "surface": surface,
            "scaffold": np.asarray(scaffolds, dtype=object)[lv],
        })

    report = {}
    for location_ridge in LOCATION_RIDGES:
        for amplitude_ridge in AMPLITUDE_RIDGES:
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
                        query = _query(
                            cache["y"], target, support, cache["scaffold"])
                        if len(query) < 3:
                            continue
                        coefficient = _solve_gain(
                            cache["surface"][support, target],
                            cache["y"][support, target],
                            location_ridge, amplitude_ridge)
                        prediction = _curve(
                            cache["surface"][query, target], coefficient)
                        values.append(float(np.mean(
                            (cache["y"][query, target] - prediction) ** 2)))
            report[(location_ridge, amplitude_ridge)] = values
    selected = min(
        (float(np.mean(value)), pair) for pair, value in report.items() if value)[1]
    return selected, {
        f"location={pair[0]},amplitude={pair[1]}": {
            "mean_mse": float(np.mean(value)), "n_episodes": len(value)
        } for pair, value in report.items() if value
    }


def _evaluate(y, additive, surface, null_surface, protein_views, scaffolds,
              location_ridge, amplitude_ridge, seeds=SEEDS):
    n_target = y.shape[1]
    per_target = {name: {
        metric: np.full((n_target, len(seeds)), np.nan)
        for metric in ("mse", "mae", "spearman", "interaction_mse")
    } for name in ARM_NAMES}
    coefficients = np.full((n_target, len(seeds), 2), np.nan)
    nearest = _nearest_nonself(protein_views)
    interaction_y = _masked_double_center(y)
    for target in range(n_target):
        wrong_target = int(nearest[target])
        finite = np.flatnonzero(
            np.isfinite(y[:, target]) & np.isfinite(y[:, wrong_target]))
        for seed_index, seed in enumerate(seeds):
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
            correct_coefficient = _solve_gain(
                surface[support, target], residual,
                location_ridge, amplitude_ridge)
            null_coefficient = _solve_gain(
                null_surface[support, target], residual,
                location_ridge, amplitude_ridge)
            nearest_coefficient = _solve_gain(
                surface[support, wrong_target], residual,
                location_ridge, amplitude_ridge)
            wrong_coefficient = _solve_gain(
                surface[support, target], wrong_residual,
                location_ridge, amplitude_ridge)
            rng = np.random.default_rng(seed + 65537 * target)
            permuted_coefficient = _solve_gain(
                surface[support, target], residual[rng.permutation(K)],
                location_ridge, amplitude_ridge)
            location = float(np.sum(residual) / (K + location_ridge))
            residual_curves = {
                "support_free": surface[:, target],
                "location_only": np.full(y.shape[0], location),
                "gauge_fixed": _curve(surface[:, target], correct_coefficient),
                "zero_protein": _curve(null_surface[:, target], null_coefficient),
                "nearest_protein": _curve(
                    surface[:, wrong_target], nearest_coefficient),
                "wrong_support": _curve(surface[:, target], wrong_coefficient),
                "permuted_support": _curve(
                    surface[:, target], permuted_coefficient),
            }
            finite_target = np.flatnonzero(np.isfinite(y[:, target]))
            for name, residual_curve in residual_curves.items():
                prediction = np.clip(
                    additive[query, target] + residual_curve[query], 0.0, 1.0)
                outcome = y[query, target]
                error = outcome - prediction
                per_target[name]["mse"][target, seed_index] = float(np.mean(error ** 2))
                per_target[name]["mae"][target, seed_index] = float(np.mean(np.abs(error)))
                per_target[name]["spearman"][target, seed_index] = _safe_spearman(
                    outcome, prediction)
                centered = residual_curve - float(np.mean(residual_curve[finite_target]))
                per_target[name]["interaction_mse"][target, seed_index] = float(
                    np.mean((interaction_y[query, target] - centered[query]) ** 2))
            coefficients[target, seed_index] = correct_coefficient
    return per_target, coefficients, nearest


def _summarize(per_target, coefficients, nearest):
    metrics = {
        name: {metric: float(np.nanmean(value)) for metric, value in values.items()}
        for name, values in per_target.items()
    }
    contrasts = {
        name: _bootstrap(values, per_target["gauge_fixed"])
        for name, values in per_target.items() if name != "gauge_fixed"
    }
    interaction = {
        name: _bootstrap(
            per_target[name], per_target["gauge_fixed"], metric="interaction_mse")
        for name in ("support_free", "nearest_protein")
    }
    return {
        "metrics": metrics,
        "mse_reduction_gauge_fixed_vs": contrasts,
        "interaction_mse_reduction_gauge_fixed_vs": interaction,
        "coordinates": {
            "tau": {"mean": float(np.nanmean(coefficients[..., 0])),
                    "std": float(np.nanstd(coefficients[..., 0]))},
            "amplitude": {"mean": float(np.nanmean(coefficients[..., 1])),
                          "std": float(np.nanstd(coefficients[..., 1]))},
        },
        "nearest_nonself_target_index": nearest.tolist(),
    }


def run(args):
    started = time.time()
    data = _load(args)
    y_source = data["y_source"]
    source_ligands = data["source_ligands"]
    source_targets = data["source_targets"]
    source_bits = _ligand_bits(source_ligands)
    source_scaffolds = np.asarray([
        record.generic_scaffold for record in source_ligands], dtype=object)
    source_groups = [item["record"].group for item in source_targets]
    masks = _subpocket_weights(args.kissim_distances)
    source_protein_views = _atlas_protein_views(source_targets, masks)

    print("selecting source-only ligand ridge", flush=True)
    ligand_ridge, ligand_cv = _select_ligand_ridge(
        y_source, source_bits, source_scaffolds)
    print("selecting gauge-fixed penalties on source-only dual-cold episodes", flush=True)
    (location_ridge, amplitude_ridge), section_cv = _cv_penalties(
        y_source, source_bits, source_protein_views, source_scaffolds,
        source_groups, ligand_ridge)
    print(f"selected location={location_ridge} amplitude={amplitude_ridge}", flush=True)

    ligand_coefficient = _fit_ligand_atlas(
        source_bits, double_center(y_source), ligand_ridge)
    temperatures = _protein_temperatures(source_protein_views)
    additive_model = _fit_additive(
        y_source,
        np.stack([record.nuisance_features for record in source_ligands]),
        np.stack([protein_model_features(item["record"].pocket)
                  for item in source_targets]),
    )
    transfers = {}
    for panel in data["panels"]:
        print(f"evaluating {panel['name']}", flush=True)
        profile = _predict_profile(
            _ligand_bits(panel["ligands"]), source_bits, ligand_coefficient)
        protein_views = _atlas_protein_views(panel["targets"], masks)
        all_surface = _atlas_surfaces(
            profile, protein_views, source_protein_views, temperatures)
        surface = all_surface[..., 0]
        null_surface = _uniform_surfaces(profile, len(panel["targets"]))[..., 0]
        additive = additive_model.predict(
            np.stack([record.nuisance_features for record in panel["ligands"]]),
            np.stack([protein_model_features(item["record"].pocket)
                      for item in panel["targets"]]),
        )
        scaffolds = np.asarray([
            record.generic_scaffold for record in panel["ligands"]], dtype=object)
        evaluated = _evaluate(
            panel["y"], additive, surface, null_surface, protein_views,
            scaffolds, location_ridge, amplitude_ridge)
        transfers[panel["name"]] = {
            "role": panel["role"], "n_ligands": len(panel["ligands"]),
            "n_targets": len(panel["targets"]), "support_size": K,
            "random": _summarize(*evaluated),
        }

    primary = transfers["PKIS2"]["random"]
    external = transfers["Anastassiadis2011"]["random"]
    pkis_raw = all(
        primary["mse_reduction_gauge_fixed_vs"][name]["ci95"][0] > 0.0
        for name in CONTROL_NAMES)
    external_raw = all(
        external["mse_reduction_gauge_fixed_vs"][name]["estimate"] > 0.0
        for name in CONTROL_NAMES)
    pkis_interaction = all(
        primary["interaction_mse_reduction_gauge_fixed_vs"][name]["ci95"][0] > 0.0
        for name in ("support_free", "nearest_protein"))
    external_interaction = all(
        external["interaction_mse_reduction_gauge_fixed_vs"][name]["estimate"] > 0.0
        for name in ("support_free", "nearest_protein"))
    passed = bool(pkis_raw and external_raw and pkis_interaction and external_interaction)
    result = {
        "schema": "MetaSieve.GaugeFixedAtlasSection.F4G.v1",
        "selected": {
            "d_adapt": 2, "k": K, "ligand_ridge": ligand_ridge,
            "location_ridge": location_ridge,
            "amplitude_ridge": amplitude_ridge,
            "protein_temperatures": temperatures.tolist(),
        },
        "source": {"panel": "PKIS1", "ligand_cv": ligand_cv,
                   "section_cv": section_cv},
        "transfers": transfers,
        "gate": {
            "pkis2_raw_controls": pkis_raw,
            "anastassiadis_raw_point_estimates": external_raw,
            "pkis2_interaction_controls": pkis_interaction,
            "anastassiadis_interaction_point_estimates": external_interaction,
            "passed": passed,
            "verdict": ("F4G_GAUGE_FIXED_SECTION_ADMISSIBLE" if passed
                        else "F4G_GAUGE_FIXED_SECTION_NOT_ADMISSIBLE"),
        },
        "read_firewall": {
            "pkis1": "source", "pkis2": "consumed_development",
            "anastassiadis": "consumed_development",
            "kcgs_numeric_outcomes": "NOT_READ", "davis_labels": "NOT_READ",
            "recipient_labels": "NOT_READ",
        },
        "limitations": [
            "The scalar affinity decoder is a feasibility bridge, not the frozen law-valued operator.",
            "The gauge fix admits only one biological interaction shape per task.",
        ],
        "elapsed_seconds": time.time() - started,
    }
    output = Path(args.output).resolve()
    _write_json(output / "result.json", result)
    _write_json(output / "manifest.json", {
        "sha256": {
            "preregistration": _sha256(Path(__file__).with_name("F4G_PREREGISTRATION.md")),
            "script": _sha256(Path(__file__)),
            "pkis1_labels": _sha256(
                Path(args.informers_root).resolve() / "data/pkis1_continuous_labels.csv"),
        }, "parameters": vars(args),
    })
    print(json.dumps(result["gate"], indent=2), flush=True)


def parser():
    item = argparse.ArgumentParser()
    item.add_argument("--informers-root", default="../external/informers")
    item.add_argument("--klifs-json", default="../external/klifs/kinase_information_human.json")
    item.add_argument("--kissim-distances",
                      default="../external/kissim/kissim/data/min_max_distances_fine.csv")
    item.add_argument("--anastassiadis-workbook",
                      default="external/anastassiadis/NIHMS328213-supplement-3.xls")
    item.add_argument("--anastassiadis-identities",
                      default="external/anastassiadis/compound_identities.json")
    item.add_argument("--output", default="research/section_operator_pilot/artifacts/f4g")
    return item


if __name__ == "__main__":
    run(parser().parse_args())
