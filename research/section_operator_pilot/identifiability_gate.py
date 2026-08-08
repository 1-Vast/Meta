"""F2C conformal identifiability domain for the bioactivity-atlas section."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from ..pkis_mechanism_pilot.mechanism import (
    double_center,
    protein_model_features,
    stable_fold,
)
from .bioactivity_atlas import (
    _atlas_protein_views,
    _atlas_surfaces,
    _bootstrap,
    _cv_penalties,
    _fit_ligand_atlas,
    _ligand_bits,
    _predict_profile,
    _protein_temperatures,
    _select_ligand_ridge,
    _uniform_surfaces,
    _solve_section,
    _apply_section,
)
from .ceiling_probe import SEEDS, _load
from .conditional_bilinear import (
    K_PRIMARY,
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


CALIBRATION_SEEDS = SEEDS[:10]
CERTIFICATE_RIDGE = 10.0
CONFORMAL_LEVEL = 0.80

ARM_NAMES = (
    "support_free", "location_only", "atlas_section", "zero_protein",
    "nearest_protein", "wrong_support", "permuted_support",
)
CONTROL_NAMES = tuple(name for name in ARM_NAMES if name != "atlas_section")


def _section_curve(surface, residual, support, location_ridge, tangent_ridge):
    coefficient = _solve_section(
        surface[support], residual, location_ridge, tangent_ridge)
    return _apply_section(surface, coefficient), coefficient


def _loo_mse(surface, residual, location_ridge, tangent_ridge):
    prediction = np.empty(len(residual), dtype=np.float64)
    for held in range(len(residual)):
        train = np.asarray([index for index in range(len(residual)) if index != held])
        coefficient = _solve_section(
            surface[train], residual[train], location_ridge, tangent_ridge)
        prediction[held] = _apply_section(surface[[held]], coefficient)[0]
    return float(np.mean((residual - prediction) ** 2))


def _certificate_features(correct_surface, null_surface, nearest_surface,
                          support, residual, location_ridge, tangent_ridge):
    correct_support = correct_surface[support]
    null_support = null_surface[support]
    nearest_support = nearest_surface[support]
    correct_loo = _loo_mse(
        correct_support, residual, location_ridge, tangent_ridge)
    null_loo = _loo_mse(
        null_support, residual, location_ridge, tangent_ridge)
    nearest_loo = _loo_mse(
        nearest_support, residual, location_ridge, tangent_ridge)
    location_prediction = []
    for held in range(len(residual)):
        train = np.asarray([index for index in range(len(residual)) if index != held])
        value = float(np.sum(residual[train]) / (len(train) + location_ridge))
        location_prediction.append(value)
    location_loo = float(np.mean(
        (residual - np.asarray(location_prediction)) ** 2))
    support_free_loo = float(np.mean(
        (residual - correct_support[:, 0]) ** 2))
    coefficient = _solve_section(
        correct_support, residual, location_ridge, tangent_ridge)
    fitted = _apply_section(correct_support, coefficient)
    design = np.concatenate([
        np.ones((len(support), 1)), correct_support[:, 1:3]
    ], axis=1)
    singular = np.linalg.svd(design, compute_uv=False)
    condition = float(singular[-1] / max(singular[0], 1e-12))
    separation_all = float(np.mean((correct_surface - nearest_surface) ** 2))
    separation_support = float(np.mean(
        (correct_support - nearest_support) ** 2))
    return np.asarray([
        correct_loo,
        support_free_loo - correct_loo,
        location_loo - correct_loo,
        null_loo - correct_loo,
        nearest_loo - correct_loo,
        condition,
        separation_all,
        separation_support,
        float(np.std(residual)),
        float(np.ptp(residual)),
        float(np.mean((residual - fitted) ** 2)),
        float(np.linalg.norm(coefficient)),
        float(np.std(correct_surface[:, 1])),
        float(np.std(correct_surface[:, 2])),
    ], dtype=np.float64)


def _episode(y, additive, interaction_y, surface, null_surface, nearest_target,
             scaffolds, target, support, seed, location_ridge, tangent_ridge):
    query = _query(y, target, support, scaffolds)
    if len(query) < 3:
        return None
    wrong_target = int(nearest_target[target])
    residual = y[support, target] - additive[support, target]
    wrong_residual = y[support, wrong_target] - additive[support, wrong_target]
    correct_curve, correct_coefficient = _section_curve(
        surface[:, target], residual, support, location_ridge, tangent_ridge)
    null_curve, _ = _section_curve(
        null_surface[:, target], residual, support, location_ridge, tangent_ridge)
    nearest_curve, _ = _section_curve(
        surface[:, wrong_target], residual, support, location_ridge, tangent_ridge)
    wrong_curve, _ = _section_curve(
        surface[:, target], wrong_residual, support, location_ridge, tangent_ridge)
    rng = np.random.default_rng(seed + 65537 * target)
    permuted_curve, _ = _section_curve(
        surface[:, target], residual[rng.permutation(len(support))], support,
        location_ridge, tangent_ridge)
    location = float(np.sum(residual) / (len(support) + location_ridge))
    residual_curves = {
        "support_free": surface[:, target, 0],
        "location_only": np.full(y.shape[0], location),
        "atlas_section": correct_curve,
        "zero_protein": null_curve,
        "nearest_protein": nearest_curve,
        "wrong_support": wrong_curve,
        "permuted_support": permuted_curve,
    }
    raw = {}
    interaction = {}
    finite = np.flatnonzero(np.isfinite(y[:, target]))
    for name, residual_curve in residual_curves.items():
        prediction = np.clip(additive[query, target] + residual_curve[query], 0.0, 1.0)
        raw[name] = {
            "mse": float(np.mean((y[query, target] - prediction) ** 2)),
            "mae": float(np.mean(np.abs(y[query, target] - prediction))),
            "spearman": _safe_spearman(y[query, target], prediction),
        }
        centered = residual_curve - float(np.mean(residual_curve[finite]))
        interaction[name] = float(np.mean(
            (interaction_y[query, target] - centered[query]) ** 2))
    feature = _certificate_features(
        surface[:, target], null_surface[:, target], surface[:, wrong_target],
        support, residual, location_ridge, tangent_ridge)
    raw_gains = [raw[name]["mse"] - raw["atlas_section"]["mse"]
                 for name in CONTROL_NAMES]
    interaction_gains = [
        interaction[name] - interaction["atlas_section"]
        for name in ("support_free", "nearest_protein")
    ]
    return {
        "query": query, "feature": feature,
        "minimum_margin": float(min(raw_gains + interaction_gains)),
        "raw": raw, "interaction": interaction,
        "coefficient": correct_coefficient,
    }


def _source_calibration(y, bits, protein_views, scaffolds, protein_groups,
                        ligand_ridge, location_ridge, tangent_ridge):
    ligand_folds = np.asarray([stable_fold(value, 3) for value in scaffolds])
    protein_folds = np.asarray([stable_fold(value, 3) for value in protein_groups])
    features, margins, folds = [], [], []
    fold_report = []
    for fold in range(3):
        lt, lv = np.flatnonzero(ligand_folds != fold), np.flatnonzero(ligand_folds == fold)
        pt, pv = np.flatnonzero(protein_folds != fold), np.flatnonzero(protein_folds == fold)
        train_residual = double_center(y[np.ix_(lt, pt)])
        coefficient = _fit_ligand_atlas(bits[lt], train_residual, ligand_ridge)
        profile = _predict_profile(bits[lv], bits[lt], coefficient)
        source_views = [value[pt] for value in protein_views]
        query_views = [value[pv] for value in protein_views]
        temperatures = _protein_temperatures(source_views)
        surface = _atlas_surfaces(profile, query_views, source_views, temperatures)
        null_surface = _uniform_surfaces(profile, len(pv))
        outcome = double_center(y[np.ix_(lv, pv)])
        additive = np.zeros_like(outcome)
        interaction_y = double_center(outcome)
        nearest = _nearest_nonself(query_views)
        local_scaffold = np.asarray(scaffolds, dtype=object)[lv]
        count_before = len(features)
        for target in range(len(pv)):
            wrong_target = int(nearest[target])
            finite = np.flatnonzero(
                np.isfinite(outcome[:, target]) & np.isfinite(outcome[:, wrong_target]))
            for seed in CALIBRATION_SEEDS:
                support = _d_optimal(
                    local_scaffold, surface[:, target, 1:3], finite,
                    K_PRIMARY, "random", seed + 104729 * target)
                if support is None:
                    continue
                episode = _episode(
                    outcome, additive, interaction_y, surface, null_surface,
                    nearest, local_scaffold, target, support, seed,
                    location_ridge, tangent_ridge)
                if episode is None:
                    continue
                features.append(episode["feature"])
                margins.append(episode["minimum_margin"])
                folds.append(fold)
        fold_report.append({"fold": fold, "episodes": len(features) - count_before,
                            "valid_ligands": len(lv), "valid_targets": len(pv)})

    features = np.asarray(features, dtype=np.float64)
    margins = np.asarray(margins, dtype=np.float64)
    folds = np.asarray(folds, dtype=np.int64)
    crossfit = np.full(len(margins), np.nan)
    for fold in range(3):
        train, valid = folds != fold, folds == fold
        model = make_pipeline(
            StandardScaler(), Ridge(alpha=CERTIFICATE_RIDGE)
        ).fit(features[train], margins[train])
        crossfit[valid] = model.predict(features[valid])
    error = crossfit - margins
    quantile = float(np.quantile(error, CONFORMAL_LEVEL, method="higher"))
    lower = crossfit - quantile
    admitted = lower > 0.0
    final = make_pipeline(
        StandardScaler(), Ridge(alpha=CERTIFICATE_RIDGE)
    ).fit(features, margins)
    report = {
        "n_episodes": len(margins), "folds": fold_report,
        "conformal_level": CONFORMAL_LEVEL,
        "error_quantile": quantile,
        "crossfit_rmse": float(np.sqrt(np.mean((crossfit - margins) ** 2))),
        "crossfit_correlation": float(np.corrcoef(crossfit, margins)[0, 1]),
        "crossfit_admission_rate": float(np.mean(admitted)),
        "crossfit_admitted_margin_mean": (
            float(np.mean(margins[admitted])) if admitted.any() else None),
        "feature_dimension": int(features.shape[1]),
    }
    return final, quantile, report


def _evaluate_selective(y, additive, surface, null_surface, protein_views,
                        scaffolds, k, policy, seeds, location_ridge,
                        tangent_ridge, certificate_model, quantile):
    n_target = y.shape[1]
    per_target = {name: {
        metric: np.full((n_target, len(seeds)), np.nan)
        for metric in ("mse", "mae", "spearman", "interaction_mse")
    } for name in ARM_NAMES}
    certificate = np.full((n_target, len(seeds)), np.nan)
    admitted = np.zeros((n_target, len(seeds)), dtype=bool)
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
            episode = _episode(
                y, additive, interaction_y, surface, null_surface, nearest,
                scaffolds, target, support, seed, location_ridge, tangent_ridge)
            if episode is None:
                continue
            predicted = float(certificate_model.predict(
                episode["feature"][None])[0])
            value = predicted - quantile
            certificate[target, seed_index] = value
            if value <= 0.0:
                continue
            admitted[target, seed_index] = True
            for name in ARM_NAMES:
                for metric in ("mse", "mae", "spearman"):
                    per_target[name][metric][target, seed_index] = episode["raw"][name][metric]
                per_target[name]["interaction_mse"][target, seed_index] = episode["interaction"][name]
    return per_target, certificate, admitted, nearest


def _selective_summary(per_target, certificate, admitted, nearest):
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
    finite = np.isfinite(certificate)
    target_admitted = np.any(admitted, axis=1)
    return {
        "coverage": {
            "eligible_episodes": int(finite.sum()),
            "admitted_episodes": int(admitted.sum()),
            "episode_rate": float(admitted.sum() / max(int(finite.sum()), 1)),
            "admitted_target_clusters": int(target_admitted.sum()),
            "certificate_mean": float(np.nanmean(certificate)),
            "certificate_admitted_mean": (
                float(np.mean(certificate[admitted])) if admitted.any() else None),
        },
        "metrics": metrics,
        "mse_reduction_atlas_vs": contrasts,
        "interaction_mse_reduction_atlas_vs": interaction,
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

    print("recomputing source-only F2A parameters", flush=True)
    ligand_ridge, ligand_cv = _select_ligand_ridge(
        y_source, source_bits, source_scaffolds)
    (location_ridge, tangent_ridge), section_cv = _cv_penalties(
        y_source, source_bits, source_protein_views, source_scaffolds,
        source_groups, ligand_ridge)
    print("fitting cross-fold conformal identifiability certificate", flush=True)
    certificate_model, quantile, calibration = _source_calibration(
        y_source, source_bits, source_protein_views, source_scaffolds,
        source_groups, ligand_ridge, location_ridge, tangent_ridge)
    print(json.dumps(calibration, indent=2), flush=True)

    source_residual = double_center(y_source)
    ligand_coefficient = _fit_ligand_atlas(
        source_bits, source_residual, ligand_ridge)
    temperatures = _protein_temperatures(source_protein_views)
    additive_model = _fit_additive(
        y_source,
        np.stack([record.nuisance_features for record in source_ligands]),
        np.stack([protein_model_features(item["record"].pocket)
                  for item in source_targets]),
    )

    transfers = {}
    for panel in data["panels"]:
        print(f"selective evaluation {panel['name']}", flush=True)
        profile = _predict_profile(
            _ligand_bits(panel["ligands"]), source_bits, ligand_coefficient)
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
            "n_targets": len(panel["targets"]), "support_size": K_PRIMARY,
            "policies": {},
        }
        for policy in ("random", "d_optimal"):
            seeds = SEEDS if policy == "random" else (SEED,)
            evaluated = _evaluate_selective(
                panel["y"], additive, surface, null_surface, protein_views,
                scaffolds, K_PRIMARY, policy, seeds, location_ridge,
                tangent_ridge, certificate_model, quantile)
            panel_result["policies"][policy] = _selective_summary(*evaluated)
        transfers[panel["name"]] = panel_result

    primary = transfers["PKIS2"]["policies"]["random"]
    external = transfers["Anastassiadis2011"]["policies"]["random"]
    pkis_coverage = (primary["coverage"]["episode_rate"] >= 0.20
                     and primary["coverage"]["admitted_target_clusters"] >= 30)
    external_coverage = (external["coverage"]["episode_rate"] >= 0.10
                         and external["coverage"]["admitted_target_clusters"] >= 10)
    pkis_raw = all(
        primary["mse_reduction_atlas_vs"][name]["ci95"][0] > 0.0
        for name in CONTROL_NAMES)
    external_raw = all(
        external["mse_reduction_atlas_vs"][name]["estimate"] > 0.0
        for name in CONTROL_NAMES)
    pkis_interaction = all(
        primary["interaction_mse_reduction_atlas_vs"][name]["ci95"][0] > 0.0
        for name in ("support_free", "nearest_protein"))
    external_interaction = all(
        external["interaction_mse_reduction_atlas_vs"][name]["estimate"] > 0.0
        for name in ("support_free", "nearest_protein"))
    passed = bool(pkis_coverage and external_coverage and pkis_raw
                  and external_raw and pkis_interaction and external_interaction)
    result = {
        "schema": "MetaSieve.ConformalIdentifiabilitySection.F2C.v1",
        "selected": {
            "d_adapt": 4, "k_primary": K_PRIMARY,
            "ligand_ridge": ligand_ridge,
            "location_ridge": location_ridge,
            "tangent_ridge": tangent_ridge,
            "certificate_ridge": CERTIFICATE_RIDGE,
            "conformal_level": CONFORMAL_LEVEL,
            "conformal_error_quantile": quantile,
        },
        "source": {
            "panel": "PKIS1", "ligand_cv": ligand_cv,
            "section_cv": section_cv, "certificate_calibration": calibration,
        },
        "transfers": transfers,
        "gate": {
            "pkis2_coverage": pkis_coverage,
            "anastassiadis_coverage": external_coverage,
            "pkis2_raw_controls": pkis_raw,
            "anastassiadis_raw_point_estimates": external_raw,
            "pkis2_interaction_controls": pkis_interaction,
            "anastassiadis_interaction_point_estimates": external_interaction,
            "passed": passed,
            "verdict": ("F2C_PARTIAL_SECTION_ADMISSIBLE" if passed
                        else "F2C_PARTIAL_SECTION_NOT_ADMISSIBLE"),
        },
        "read_firewall": {
            "pkis1": "source", "pkis2": "consumed_development",
            "anastassiadis": "consumed_development",
            "kcgs_numeric_outcomes": "NOT_READ", "davis_labels": "NOT_READ",
            "recipient_labels": "NOT_READ",
        },
        "limitations": [
            "The 80% conformal certificate is marginal under source-episode exchangeability, not conditional under panel shift.",
            "Abstention is required outside the admitted partial task domain.",
            "The scalar decoder remains a feasibility diagnostic only.",
        ],
        "elapsed_seconds": time.time() - started,
    }
    output = Path(args.output).resolve()
    _write_json(output / "result.json", result)
    _write_json(output / "manifest.json", {
        "sha256": {
            "preregistration": _sha256(Path(__file__).with_name("F2C_PREREGISTRATION.md")),
            "script": _sha256(Path(__file__)),
            "f2a_script": _sha256(Path(__file__).with_name("bioactivity_atlas.py")),
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
    item.add_argument("--output", default="research/section_operator_pilot/artifacts/f2c")
    return item


if __name__ == "__main__":
    run(parser().parse_args())
