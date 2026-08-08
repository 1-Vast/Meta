"""F5C counterfactual Fisher support design for protein identifiability."""
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
from .ceiling_probe import _load
from .conditional_bilinear import (
    SEED, _fit_additive, _masked_double_center, _nearest_nonself, _query,
    _safe_spearman, _sha256, _subpocket_weights, _write_json,
)
from .gauge_fixed import AMPLITUDE_RIDGES, K, LOCATION_RIDGES
from .gauge_orthogonal import _curve, _solve_orthogonal


NONINFERIORITY_MARGIN = 1e-4
ARM_NAMES = (
    "support_free", "location_only", "counterfactual", "zero_protein",
    "nearest_protein", "wrong_support", "permuted_support",
)
CONTROL_NAMES = tuple(name for name in ARM_NAMES if name != "counterfactual")


def _counterfactual_support(scaffolds, correct_surface, wrong_surface, finite, k=K):
    by_scaffold = {}
    for index in np.asarray(finite, dtype=np.int64):
        by_scaffold.setdefault(str(scaffolds[index]), []).append(int(index))
    if len(by_scaffold) < k:
        return None
    representatives = np.asarray([
        min(indices) for _, indices in sorted(by_scaffold.items())
    ], dtype=np.int64)
    delta = correct_surface[representatives] - wrong_surface[representatives]
    design = np.concatenate([
        np.ones((len(representatives), 1)),
        correct_surface[representatives, :1], delta,
    ], axis=1)
    scale = np.std(design[:, 1:], axis=0)
    scale = np.where(scale > 1e-8, scale, 1.0)
    design[:, 1:] /= scale[None, :]
    gram = np.eye(design.shape[1], dtype=np.float64) * 1e-3
    selected, remaining = [], list(range(len(representatives)))
    for _ in range(k):
        scores = []
        for candidate in remaining:
            sign, value = np.linalg.slogdet(
                gram + np.outer(design[candidate], design[candidate]))
            scores.append(value if sign > 0 else -np.inf)
        best = max(scores)
        tied = [remaining[index] for index, value in enumerate(scores)
                if np.isclose(value, best, rtol=1e-12, atol=1e-14)]
        chosen = min(tied)
        selected.append(chosen)
        gram += np.outer(design[chosen], design[chosen])
        remaining.remove(chosen)
    return np.sort(representatives[np.asarray(selected, dtype=np.int64)])


def _cv_penalties(y, bits, protein_views, scaffolds, protein_groups, ligand_ridge):
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
            profile, query_views, source_views, _protein_temperatures(source_views))
        caches.append({
            "y": double_center(y[np.ix_(lv, pv)]), "surface": surface,
            "nearest": _nearest_nonself(query_views),
            "scaffold": np.asarray(scaffolds, dtype=object)[lv],
        })
    report = {}
    for location_ridge in LOCATION_RIDGES:
        for amplitude_ridge in AMPLITUDE_RIDGES:
            values = []
            for cache in caches:
                for target in range(cache["y"].shape[1]):
                    wrong = int(cache["nearest"][target])
                    finite = np.flatnonzero(
                        np.isfinite(cache["y"][:, target])
                        & np.isfinite(cache["y"][:, wrong]))
                    support = _counterfactual_support(
                        cache["scaffold"], cache["surface"][:, target],
                        cache["surface"][:, wrong], finite)
                    if support is None:
                        continue
                    query = _query(cache["y"], target, support, cache["scaffold"])
                    if len(query) < 3:
                        continue
                    coordinate = _solve_orthogonal(
                        cache["surface"][support, target, 0],
                        cache["y"][support, target], location_ridge, amplitude_ridge)
                    prediction = _curve(cache["surface"][query, target, 0], coordinate)
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
              location_ridge, amplitude_ridge):
    n_target = y.shape[1]
    per_target = {name: {metric: np.full((n_target, 1), np.nan)
                         for metric in ("mse", "mae", "spearman", "interaction_mse")}
                  for name in ARM_NAMES}
    coordinates = np.full((n_target, 1, 3), np.nan)
    nearest = _nearest_nonself(protein_views)
    interaction_y = _masked_double_center(y)
    support_indices = {}
    for target in range(n_target):
        wrong_target = int(nearest[target])
        finite = np.flatnonzero(
            np.isfinite(y[:, target]) & np.isfinite(y[:, wrong_target]))
        support = _counterfactual_support(
            scaffolds, surface[:, target], surface[:, wrong_target], finite)
        if support is None:
            continue
        support_indices[str(target)] = support.tolist()
        query = _query(y, target, support, scaffolds)
        if len(query) < 3:
            continue
        residual = y[support, target] - additive[support, target]
        wrong_residual = y[support, wrong_target] - additive[support, wrong_target]
        correct = _solve_orthogonal(
            surface[support, target, 0], residual, location_ridge, amplitude_ridge)
        null = _solve_orthogonal(
            null_surface[support, target, 0], residual, location_ridge, amplitude_ridge)
        deranged = _solve_orthogonal(
            surface[support, wrong_target, 0], residual, location_ridge, amplitude_ridge)
        wrong = _solve_orthogonal(
            surface[support, target, 0], wrong_residual, location_ridge, amplitude_ridge)
        rng = np.random.default_rng(SEED + 65537 * target)
        permuted = _solve_orthogonal(
            surface[support, target, 0], residual[rng.permutation(K)],
            location_ridge, amplitude_ridge)
        location = float(np.sum(residual) / (K + location_ridge))
        curves = {
            "support_free": surface[:, target, 0],
            "location_only": np.full(y.shape[0], location),
            "counterfactual": _curve(surface[:, target, 0], correct),
            "zero_protein": _curve(null_surface[:, target, 0], null),
            "nearest_protein": _curve(surface[:, wrong_target, 0], deranged),
            "wrong_support": _curve(surface[:, target, 0], wrong),
            "permuted_support": _curve(surface[:, target, 0], permuted),
        }
        finite_target = np.flatnonzero(np.isfinite(y[:, target]))
        for name, residual_curve in curves.items():
            prediction = np.clip(additive[query, target] + residual_curve[query], 0, 1)
            outcome = y[query, target]
            error = outcome - prediction
            per_target[name]["mse"][target, 0] = float(np.mean(error ** 2))
            per_target[name]["mae"][target, 0] = float(np.mean(np.abs(error)))
            per_target[name]["spearman"][target, 0] = _safe_spearman(outcome, prediction)
            centered = residual_curve - float(np.mean(residual_curve[finite_target]))
            per_target[name]["interaction_mse"][target, 0] = float(
                np.mean((interaction_y[query, target] - centered[query]) ** 2))
        coordinates[target, 0] = correct
    return per_target, coordinates, nearest, support_indices


def _summarize(per_target, coordinates, nearest, support_indices):
    metrics = {name: {metric: float(np.nanmean(value)) for metric, value in values.items()}
               for name, values in per_target.items()}
    contrasts = {name: _bootstrap(values, per_target["counterfactual"])
                 for name, values in per_target.items() if name != "counterfactual"}
    interaction = {name: _bootstrap(
        per_target[name], per_target["counterfactual"], metric="interaction_mse")
        for name in ("support_free", "nearest_protein")}
    return {
        "metrics": metrics, "mse_reduction_counterfactual_vs": contrasts,
        "interaction_mse_reduction_counterfactual_vs": interaction,
        "coordinates": {
            "location": {"mean": float(np.nanmean(coordinates[..., 0])),
                         "std": float(np.nanstd(coordinates[..., 0]))},
            "amplitude": {"mean": float(np.nanmean(coordinates[..., 1])),
                          "std": float(np.nanstd(coordinates[..., 1]))}},
        "nearest_nonself_target_index": nearest.tolist(),
        "support_indices": support_indices,
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
    (location_ridge, amplitude_ridge), section_cv = _cv_penalties(
        y_source, source_bits, source_protein_views, source_scaffolds,
        source_groups, ligand_ridge)
    print(f"selected location={location_ridge} amplitude={amplitude_ridge}", flush=True)
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
        surface = _atlas_surfaces(profile, protein_views, source_protein_views, temperatures)
        null_surface = _uniform_surfaces(profile, len(panel["targets"]))
        additive = additive_model.predict(
            np.stack([record.nuisance_features for record in panel["ligands"]]),
            np.stack([protein_model_features(item["record"].pocket) for item in panel["targets"]]))
        scaffolds = np.asarray([record.generic_scaffold for record in panel["ligands"]], dtype=object)
        transfers[panel["name"]] = {
            "role": panel["role"], "n_ligands": len(panel["ligands"]),
            "n_targets": len(panel["targets"]), "support_size": K,
            "counterfactual_fisher": _summarize(*_evaluate(
                panel["y"], additive, surface, null_surface, protein_views,
                scaffolds, location_ridge, amplitude_ridge)),
        }
    primary = transfers["PKIS2"]["counterfactual_fisher"]
    external = transfers["Anastassiadis2011"]["counterfactual_fisher"]
    pkis_raw = all(primary["mse_reduction_counterfactual_vs"][n]["ci95"][0] > 0
                   for n in CONTROL_NAMES)
    external_raw = all(external["mse_reduction_counterfactual_vs"][n]["estimate"] > 0
                       for n in CONTROL_NAMES)
    pkis_nearest = primary["interaction_mse_reduction_counterfactual_vs"][
        "nearest_protein"]["ci95"][0] > 0
    external_nearest = external["interaction_mse_reduction_counterfactual_vs"][
        "nearest_protein"]["estimate"] > 0
    pkis_noninferior = primary["interaction_mse_reduction_counterfactual_vs"][
        "support_free"]["ci95"][0] > -NONINFERIORITY_MARGIN
    external_noninferior = external["interaction_mse_reduction_counterfactual_vs"][
        "support_free"]["estimate"] > -NONINFERIORITY_MARGIN
    passed = bool(pkis_raw and external_raw and pkis_nearest and external_nearest
                  and pkis_noninferior and external_noninferior)
    result = {
        "schema": "MetaSieve.CounterfactualFisherSection.F5C.v1",
        "selected": {"d_adapt": 2, "k": K, "ligand_ridge": ligand_ridge,
                     "location_ridge": location_ridge, "amplitude_ridge": amplitude_ridge,
                     "interaction_noninferiority_margin": NONINFERIORITY_MARGIN},
        "source": {"panel": "PKIS1", "ligand_cv": ligand_cv, "section_cv": section_cv},
        "transfers": transfers,
        "gate": {"pkis2_raw_controls": pkis_raw,
                 "anastassiadis_raw_point_estimates": external_raw,
                 "pkis2_nearest_interaction": bool(pkis_nearest),
                 "anastassiadis_nearest_interaction": bool(external_nearest),
                 "pkis2_support_free_noninferiority": bool(pkis_noninferior),
                 "anastassiadis_support_free_noninferiority": bool(external_noninferior),
                 "passed": passed,
                 "verdict": ("F5C_COUNTERFACTUAL_SECTION_ADMISSIBLE" if passed
                             else "F5C_COUNTERFACTUAL_SECTION_NOT_ADMISSIBLE")},
        "read_firewall": {"kcgs_numeric_outcomes": "NOT_READ", "davis_labels": "NOT_READ",
                          "recipient_labels": "NOT_READ"},
        "elapsed_seconds": time.time() - started,
    }
    output = Path(args.output).resolve()
    _write_json(output / "result.json", result)
    _write_json(output / "manifest.json", {"sha256": {
        "preregistration": _sha256(Path(__file__).with_name("F5C_PREREGISTRATION.md")),
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
    item.add_argument("--output", default="research/section_operator_pilot/artifacts/f5c")
    return item


if __name__ == "__main__":
    run(parser().parse_args())
