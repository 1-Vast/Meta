"""F2A protein-anchored bioactivity-atlas section pilot."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np

from ..pkis_mechanism_pilot.mechanism import (
    double_center,
    protein_model_features,
    protein_pair_properties,
    stable_fold,
)
from .ceiling_probe import BOOTSTRAP_DRAWS, SEEDS, _load
from .conditional_bilinear import (
    ADDITIVE_RIDGE,
    K_PRIMARY,
    LOCATION_RIDGE,
    PROTEIN_PROPERTIES,
    SEED,
    SUBPOCKET_NAMES,
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


LIGAND_RIDGES = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
LOCATION_RIDGES = (1.0, 10.0, 100.0)
TANGENT_RIDGES = (0.01, 0.1, 1.0, 10.0)
ANCHORS_PER_VIEW = 8


def _ligand_bits(records):
    return np.stack([record.nuisance_features[:1024] for record in records]).astype(np.float64)


def _tanimoto(left, right):
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    intersection = left @ right.T
    union = left.sum(axis=1, keepdims=True) + right.sum(axis=1)[None, :] - intersection
    return np.divide(intersection, union, out=np.zeros_like(intersection), where=union > 0)


def _atlas_protein_views(target_items, masks):
    normalizer = np.asarray([3.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0], dtype=np.float64)
    views = []
    for name, properties in zip(SUBPOCKET_NAMES, PROTEIN_PROPERTIES):
        rows = []
        for item in target_items:
            value = protein_pair_properties(item["record"].pocket).astype(np.float64)
            value = value / normalizer[:, None]
            rows.append((value[list(properties)] * masks[name][None, :]).reshape(-1))
        views.append(np.asarray(rows, dtype=np.float64))
    return views


def _squared_distance(left, right):
    left = np.asarray(left, dtype=np.float64)
    right = np.asarray(right, dtype=np.float64)
    value = (
        np.square(left).mean(axis=1, keepdims=True)
        + np.square(right).mean(axis=1)[None, :]
        - 2.0 * (left @ right.T) / left.shape[1]
    )
    return np.clip(value, 0.0, None)


def _protein_temperatures(source_views):
    values = []
    for view in source_views:
        distance = _squared_distance(view, view)
        upper = distance[np.triu_indices(len(distance), 1)]
        nonzero = upper[upper > 1e-12]
        values.append(float(np.median(nonzero)) if len(nonzero) else 1.0)
    return np.asarray(values, dtype=np.float64)


def _anchor_weights(query_views, source_views, temperatures,
                    n_anchor=ANCHORS_PER_VIEW):
    weights = []
    for query, source, temperature in zip(query_views, source_views, temperatures):
        distance = _squared_distance(query, source)
        value = np.zeros_like(distance)
        keep = min(int(n_anchor), source.shape[0])
        nearest = np.argpartition(distance, keep - 1, axis=1)[:, :keep]
        for row in range(len(query)):
            index = nearest[row]
            raw = np.exp(-distance[row, index] / max(float(temperature), 1e-12))
            value[row, index] = raw / max(float(raw.sum()), 1e-12)
        weights.append(value)
    return weights


def _fit_ligand_atlas(source_bits, residual, ridge):
    kernel = _tanimoto(source_bits, source_bits)
    coefficient = np.linalg.solve(
        kernel + float(ridge) * np.eye(len(kernel)), residual)
    return coefficient


def _predict_profile(query_bits, source_bits, coefficient):
    return _tanimoto(query_bits, source_bits) @ coefficient


def _atlas_surfaces(profile, query_protein_views, source_protein_views,
                    temperatures):
    weights = _anchor_weights(
        query_protein_views, source_protein_views, temperatures)
    view_surface = np.stack([
        profile @ value.T for value in weights
    ], axis=-1)
    surface = np.stack([
        np.mean(view_surface, axis=-1),
        view_surface[..., 0] - view_surface[..., 1],
        view_surface[..., 2] - view_surface[..., 1],
    ], axis=-1)
    # Query features are label-free, so task-wise centering is admissible.
    return surface - np.mean(surface, axis=0, keepdims=True)


def _uniform_surfaces(profile, n_target):
    value = np.mean(profile, axis=1)
    value = value - float(np.mean(value))
    out = np.zeros((len(profile), n_target, 3), dtype=np.float64)
    out[:, :, 0] = value[:, None]
    return out


def _solve_section(surface, residual, location_ridge, tangent_ridge):
    design = np.concatenate([
        np.ones((len(surface), 1)), surface[:, 1:3]
    ], axis=1)
    penalty = np.diag([
        float(location_ridge), float(tangent_ridge), float(tangent_ridge)])
    return np.linalg.solve(
        design.T @ design + penalty,
        design.T @ (residual - surface[:, 0]),
    )


def _apply_section(surface, coefficient):
    design = np.concatenate([
        np.ones((len(surface), 1)), surface[:, 1:3]
    ], axis=1)
    return surface[:, 0] + design @ coefficient


def _select_ligand_ridge(y, bits, scaffolds):
    folds = np.asarray([stable_fold(value, 3) for value in scaffolds])
    scores = {ridge: [] for ridge in LIGAND_RIDGES}
    for fold in range(3):
        train, valid = np.flatnonzero(folds != fold), np.flatnonzero(folds == fold)
        train_residual = double_center(y[train])
        valid_residual = double_center(y[valid])
        train_kernel = _tanimoto(bits[train], bits[train])
        cross = _tanimoto(bits[valid], bits[train])
        for ridge in LIGAND_RIDGES:
            coefficient = np.linalg.solve(
                train_kernel + ridge * np.eye(len(train)), train_residual)
            prediction = cross @ coefficient
            scores[ridge].append(float(np.mean((valid_residual - prediction) ** 2)))
    selected = min((float(np.mean(value)), ridge) for ridge, value in scores.items())[1]
    return float(selected), {
        str(ridge): {"fold_mse": value, "mean_mse": float(np.mean(value))}
        for ridge, value in scores.items()
    }


def _cv_penalties(y, bits, protein_views, scaffolds, protein_groups, ligand_ridge):
    ligand_folds = np.asarray([stable_fold(value, 3) for value in scaffolds])
    protein_folds = np.asarray([stable_fold(value, 3) for value in protein_groups])
    caches = []
    for fold in range(3):
        lt, lv = np.flatnonzero(ligand_folds != fold), np.flatnonzero(ligand_folds == fold)
        pt, pv = np.flatnonzero(protein_folds != fold), np.flatnonzero(protein_folds == fold)
        residual_train = double_center(y[np.ix_(lt, pt)])
        coefficient = _fit_ligand_atlas(bits[lt], residual_train, ligand_ridge)
        profile = _predict_profile(bits[lv], bits[lt], coefficient)
        source_views = [value[pt] for value in protein_views]
        query_views = [value[pv] for value in protein_views]
        temperatures = _protein_temperatures(source_views)
        surfaces = _atlas_surfaces(profile, query_views, source_views, temperatures)
        caches.append({
            "y": double_center(y[np.ix_(lv, pv)]),
            "surface": surfaces,
            "scaffold": np.asarray(scaffolds, dtype=object)[lv],
        })

    report = {}
    for location_ridge in LOCATION_RIDGES:
        for tangent_ridge in TANGENT_RIDGES:
            values = []
            for cache in caches:
                for target in range(cache["y"].shape[1]):
                    finite = np.flatnonzero(np.isfinite(cache["y"][:, target]))
                    support = _d_optimal(
                        cache["scaffold"], cache["surface"][:, target, 1:3],
                        finite, K_PRIMARY, "d_optimal", SEED + target)
                    if support is None:
                        continue
                    query = _query(
                        cache["y"], target, support, cache["scaffold"])
                    if len(query) < 3:
                        continue
                    coefficient = _solve_section(
                        cache["surface"][support, target],
                        cache["y"][support, target],
                        location_ridge, tangent_ridge)
                    prediction = _apply_section(
                        cache["surface"][query, target], coefficient)
                    values.append(float(np.mean(
                        (cache["y"][query, target] - prediction) ** 2)))
            report[(location_ridge, tangent_ridge)] = values
    selected = min(
        (float(np.mean(value)), pair) for pair, value in report.items() if value)[1]
    return selected, {
        f"location={pair[0]},tangent={pair[1]}": {
            "mean_mse": float(np.mean(value)), "n_episodes": len(value)
        } for pair, value in report.items() if value
    }


def _evaluate(y, additive, surface, null_surface, protein_views, scaffolds, k,
              policy, seeds, location_ridge, tangent_ridge):
    names = (
        "support_free", "location_only", "atlas_section", "zero_protein",
        "nearest_protein", "wrong_support", "permuted_support",
    )
    n_target = y.shape[1]
    per_target = {name: {metric: np.full((n_target, len(seeds)), np.nan)
                         for metric in ("mse", "mae", "spearman", "interaction_mse")}
                  for name in names}
    coefficients = np.full((n_target, len(seeds), 3), np.nan)
    condition = np.full((n_target, len(seeds)), np.nan)
    nearest = _nearest_nonself(protein_views)
    interaction_y = _masked_double_center(y)

    for target in range(n_target):
        wrong_target = int(nearest[target])
        finite = np.flatnonzero(
            np.isfinite(y[:, target]) & np.isfinite(y[:, wrong_target]))
        for seed_index, seed in enumerate(seeds):
            support = _d_optimal(
                scaffolds, surface[:, target, 1:3], finite, k, policy,
                seed + 104729 * target)
            if support is None:
                continue
            query = _query(y, target, support, scaffolds)
            if len(query) < 3:
                continue
            residual = y[support, target] - additive[support, target]
            wrong_residual = y[support, wrong_target] - additive[support, wrong_target]
            correct_coef = _solve_section(
                surface[support, target], residual, location_ridge, tangent_ridge)
            null_coef = _solve_section(
                null_surface[support, target], residual, location_ridge, tangent_ridge)
            nearest_coef = _solve_section(
                surface[support, wrong_target], residual, location_ridge, tangent_ridge)
            wrong_coef = _solve_section(
                surface[support, target], wrong_residual, location_ridge, tangent_ridge)
            rng = np.random.default_rng(seed + 65537 * target)
            permuted_coef = _solve_section(
                surface[support, target], residual[rng.permutation(k)],
                location_ridge, tangent_ridge)
            location = float(np.sum(residual) / (k + location_ridge))
            residual_prediction = {
                "support_free": surface[query, target, 0],
                "location_only": np.full(len(query), location),
                "atlas_section": _apply_section(surface[query, target], correct_coef),
                "zero_protein": _apply_section(
                    null_surface[query, target], null_coef),
                "nearest_protein": _apply_section(
                    surface[query, wrong_target], nearest_coef),
                "wrong_support": _apply_section(
                    surface[query, target], wrong_coef),
                "permuted_support": _apply_section(
                    surface[query, target], permuted_coef),
            }
            outcome = y[query, target]
            interaction_outcome = interaction_y[query, target]
            for name, residual_value in residual_prediction.items():
                prediction = np.clip(additive[query, target] + residual_value, 0.0, 1.0)
                error = outcome - prediction
                per_target[name]["mse"][target, seed_index] = float(np.mean(error ** 2))
                per_target[name]["mae"][target, seed_index] = float(np.mean(np.abs(error)))
                per_target[name]["spearman"][target, seed_index] = _safe_spearman(
                    outcome, prediction)
                per_target[name]["interaction_mse"][target, seed_index] = float(
                    np.mean((interaction_outcome - residual_value) ** 2))
            coefficients[target, seed_index] = correct_coef
            design = np.concatenate([
                np.ones((k, 1)), surface[support, target, 1:3]
            ], axis=1)
            singular = np.linalg.svd(design, compute_uv=False)
            condition[target, seed_index] = float(singular[-1] / max(singular[0], 1e-12))
    return per_target, coefficients, condition, nearest


def _bootstrap(comparator, candidate, metric="mse", seed=SEED):
    difference = (np.nanmean(comparator[metric], axis=1)
                  - np.nanmean(candidate[metric], axis=1))
    difference = difference[np.isfinite(difference)]
    if not len(difference):
        return {"estimate": None, "ci95": [None, None], "n_targets": 0}
    rng = np.random.default_rng(seed)
    draws = difference[rng.integers(
        0, len(difference), size=(BOOTSTRAP_DRAWS, len(difference)))].mean(axis=1)
    return {
        "estimate": float(np.mean(difference)),
        "ci95": [float(np.quantile(draws, 0.025)),
                 float(np.quantile(draws, 0.975))],
        "n_targets": int(len(difference)), "draws": BOOTSTRAP_DRAWS, "seed": seed,
    }


def _summarize(per_target, coefficients, condition, nearest):
    metrics = {
        name: {metric: float(np.nanmean(value)) for metric, value in values.items()}
        for name, values in per_target.items()
    }
    contrasts = {
        name: _bootstrap(values, per_target["atlas_section"])
        for name, values in per_target.items() if name != "atlas_section"
    }
    interaction = {
        name: _bootstrap(
            per_target[name], per_target["atlas_section"], metric="interaction_mse")
        for name in ("support_free", "nearest_protein")
    }
    finite_condition = condition[np.isfinite(condition)]
    coefficient_summary = {}
    for index, name in enumerate(("tau", "u_hinge_minus_dfg", "u_front_minus_dfg")):
        value = coefficients[..., index]
        coefficient_summary[name] = {
            "mean": float(np.nanmean(value)), "std": float(np.nanstd(value)),
            "q05": float(np.nanquantile(value, 0.05)),
            "q95": float(np.nanquantile(value, 0.95)),
        }
    return {
        "metrics": metrics,
        "mse_reduction_atlas_vs": contrasts,
        "interaction_mse_reduction_atlas_vs": interaction,
        "coordinates": coefficient_summary,
        "condition_ratio": {
            "median": float(np.median(finite_condition)),
            "q10": float(np.quantile(finite_condition, 0.1)),
            "q90": float(np.quantile(finite_condition, 0.9)),
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

    print("selecting ligand-atlas ridge on scaffold-cold PKIS1", flush=True)
    ligand_ridge, ligand_cv = _select_ligand_ridge(
        y_source, source_bits, source_scaffolds)
    print(f"selected ligand ridge={ligand_ridge}", flush=True)
    print("selecting section penalties on dual-cold PKIS1", flush=True)
    (location_ridge, tangent_ridge), section_cv = _cv_penalties(
        y_source, source_bits, source_protein_views, source_scaffolds,
        source_groups, ligand_ridge)
    print(f"selected location={location_ridge} tangent={tangent_ridge}", flush=True)

    source_residual = double_center(y_source)
    ligand_coefficient = _fit_ligand_atlas(
        source_bits, source_residual, ligand_ridge)
    temperatures = _protein_temperatures(source_protein_views)
    source_x_ligand = np.stack([
        record.nuisance_features for record in source_ligands])
    source_x_protein = np.stack([
        protein_model_features(item["record"].pocket) for item in source_targets])
    additive_model = _fit_additive(
        y_source, source_x_ligand, source_x_protein)

    transfers = {}
    for panel in data["panels"]:
        print(f"evaluating {panel['name']}", flush=True)
        bits = _ligand_bits(panel["ligands"])
        profile = _predict_profile(bits, source_bits, ligand_coefficient)
        protein_views = _atlas_protein_views(panel["targets"], masks)
        surface = _atlas_surfaces(
            profile, protein_views, source_protein_views, temperatures)
        null_surface = _uniform_surfaces(profile, len(panel["targets"]))
        additive = additive_model.predict(
            np.stack([record.nuisance_features for record in panel["ligands"]]),
            np.stack([protein_model_features(item["record"].pocket)
                      for item in panel["targets"]]),
        )
        scaffolds = np.asarray([
            record.generic_scaffold for record in panel["ligands"]], dtype=object)
        panel_result = {
            "role": panel["role"], "n_ligands": len(panel["ligands"]),
            "n_targets": len(panel["targets"]),
            "n_finite_cells": int(np.isfinite(panel["y"]).sum()),
            "support_sizes": {},
        }
        for k in (5, 20):
            panel_result["support_sizes"][str(k)] = {}
            for policy in ("random", "d_optimal"):
                seeds = SEEDS if policy == "random" else (SEED,)
                evaluated = _evaluate(
                    panel["y"], additive, surface, null_surface,
                    protein_views, scaffolds, k, policy, seeds,
                    location_ridge, tangent_ridge)
                panel_result["support_sizes"][str(k)][policy] = _summarize(*evaluated)
        transfers[panel["name"]] = panel_result

    primary = transfers["PKIS2"]["support_sizes"]["5"]["d_optimal"]
    external = transfers["Anastassiadis2011"]["support_sizes"]["5"]["d_optimal"]
    required = (
        "support_free", "location_only", "zero_protein", "nearest_protein",
        "wrong_support", "permuted_support",
    )
    pkis_raw = all(
        primary["mse_reduction_atlas_vs"][name]["ci95"][0] > 0.0
        for name in required)
    external_raw = all(
        external["mse_reduction_atlas_vs"][name]["estimate"] > 0.0
        for name in required)
    pkis_interaction = all(
        primary["interaction_mse_reduction_atlas_vs"][name]["ci95"][0] > 0.0
        for name in ("support_free", "nearest_protein"))
    passed = bool(pkis_raw and external_raw and pkis_interaction)
    result = {
        "schema": "MetaSieve.ProteinAnchoredBioactivityAtlas.F2A.v1",
        "selected": {
            "d_adapt": 3, "k_primary": K_PRIMARY,
            "ligand_ridge": ligand_ridge,
            "location_ridge": location_ridge,
            "tangent_ridge": tangent_ridge,
            "anchors_per_view": ANCHORS_PER_VIEW,
            "protein_temperatures": temperatures.tolist(),
            "additive_ridge": ADDITIVE_RIDGE,
        },
        "source": {
            "panel": "PKIS1", "n_ligands": len(source_ligands),
            "n_targets": len(source_targets),
            "interaction_variance_fraction": float(
                np.var(source_residual) / np.var(y_source)),
            "ligand_cv": ligand_cv, "section_cv": section_cv,
        },
        "transfers": transfers,
        "gate": {
            "pkis2_raw_all_controls": pkis_raw,
            "anastassiadis_raw_point_estimates": external_raw,
            "pkis2_interaction_controls": pkis_interaction,
            "passed": passed,
            "verdict": ("F2A_ATLAS_SECTION_ADMISSIBLE" if passed
                        else "F2A_ATLAS_SECTION_NOT_ADMISSIBLE"),
        },
        "read_firewall": {
            "pkis1": "source", "pkis2": "consumed_development",
            "anastassiadis": "consumed_development",
            "kcgs_numeric_outcomes": "NOT_READ", "davis_labels": "NOT_READ",
            "recipient_labels": "NOT_READ",
        },
        "limitations": [
            "Profile-QSAR is prior art; novelty is claimed only for the low-dimensional KLIFS-conditioned section construction.",
            "The scalar decoder is a feasibility diagnostic, not the frozen law-valued operator.",
            "PKIS2 and Anastassiadis2011 are consumed development panels.",
        ],
        "elapsed_seconds": time.time() - started,
    }
    output = Path(args.output).resolve()
    _write_json(output / "result.json", result)
    _write_json(output / "manifest.json", {
        "sha256": {
            "preregistration": _sha256(Path(__file__).with_name("F2A_PREREGISTRATION.md")),
            "script": _sha256(Path(__file__)),
            "kissim_distances": _sha256(Path(args.kissim_distances).resolve()),
            "pkis1_labels": _sha256(
                Path(args.informers_root).resolve() / "data/pkis1_continuous_labels.csv"),
        },
        "parameters": vars(args),
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
    item.add_argument("--output", default="research/section_operator_pilot/artifacts/f2a")
    return item


if __name__ == "__main__":
    run(parser().parse_args())
