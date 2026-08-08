"""Run the single preregistered v2 product-kernel representation revision.

The script fits PKIS1 once, reports the already-consumed PKIS2 development
transfer, then evaluates the same frozen estimator on the external
Anastassiadis catalytic-activity panel.  It does not import or mutate model/.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import subprocess
import time
import warnings

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .kernel_revision import (
    LAMBDA_GRID,
    bounded_z,
    cross_rbf,
    fit_bounds,
    fit_convex_channel_weights,
    fit_separable_krr,
    kernel_ligand_from_smiles,
    predict_separable,
    protein_channel_features,
    source_rbf,
    tanimoto_kernel,
)
from .mechanism import (
    CHANNEL_NAMES,
    deterministic_derangement,
    double_center,
    load_klifs,
    load_smiles,
    map_target,
    normalize_name,
    protein_model_features,
    sha256_file,
    stable_fold,
)
from .run_pilot import _grouped_alpha, _target_mapping


SEED = 20260808
BOOTSTRAP_DRAWS = 10_000


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def _git_value(root: Path, *arguments: str) -> str:
    try:
        return subprocess.check_output(["git", *arguments], cwd=root, text=True).strip()
    except Exception:
        return "unavailable"


def _ligand_records(items: list[tuple[str, str]]):
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures

    factory = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    records = []
    exclusions = []
    for index, (molid, smiles) in enumerate(items):
        try:
            records.append(kernel_ligand_from_smiles(str(molid), str(smiles), factory=factory))
        except Exception as error:
            exclusions.append({"molid": str(molid), "reason": type(error).__name__, "detail": str(error)})
        if (index + 1) % 100 == 0:
            print(f"ligand features: {index + 1}/{len(items)}", flush=True)
    return records, exclusions


def _masked_double_center(matrix: np.ndarray) -> np.ndarray:
    value = np.asarray(matrix, dtype=np.float64)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        row = np.nanmean(value, axis=1, keepdims=True)
        column = np.nanmean(value, axis=0, keepdims=True)
        grand = float(np.nanmean(value))
    return value - row - column + grand


def _safe_pearson(y, prediction):
    mask = np.isfinite(y) & np.isfinite(prediction)
    if mask.sum() < 3 or np.std(y[mask]) < 1e-12 or np.std(prediction[mask]) < 1e-12:
        return 0.0
    return float(pearsonr(y[mask], prediction[mask]).statistic)


def _safe_spearman(y, prediction):
    mask = np.isfinite(y) & np.isfinite(prediction)
    if mask.sum() < 3 or np.std(y[mask]) < 1e-12 or np.std(prediction[mask]) < 1e-12:
        return 0.0
    return float(spearmanr(y[mask], prediction[mask]).statistic)


def _target_metrics(y: np.ndarray, arms: dict[str, np.ndarray], interaction: bool):
    summary = {}
    per_target = {}
    for name, prediction in arms.items():
        mse, secondary, counts = [], [], []
        for target in range(y.shape[1]):
            mask = np.isfinite(y[:, target]) & np.isfinite(prediction[:, target])
            counts.append(int(mask.sum()))
            if not mask.any():
                mse.append(np.nan)
                secondary.append(np.nan)
                continue
            error = y[mask, target] - prediction[mask, target]
            mse.append(float(np.mean(np.square(error))))
            secondary.append(
                _safe_pearson(y[:, target], prediction[:, target]) if interaction
                else _safe_spearman(y[:, target], prediction[:, target])
            )
        mse = np.asarray(mse, dtype=np.float64)
        secondary = np.asarray(secondary, dtype=np.float64)
        per_target[name] = {"mse": mse, "secondary": secondary, "n": np.asarray(counts)}
        summary[name] = {
            "target_macro_mse": float(np.nanmean(mse)),
            ("target_macro_pearson" if interaction else "target_macro_spearman"):
                float(np.nanmean(secondary)),
            "finite_target_clusters": int(np.isfinite(mse).sum()),
        }
    return summary, per_target


def _bootstrap(values, seed=SEED):
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return {"estimate": None, "ci95": [None, None], "n_target_clusters": 0,
                "draws": BOOTSTRAP_DRAWS, "seed": seed}
    rng = np.random.default_rng(seed)
    samples = values[rng.integers(0, len(values), size=(BOOTSTRAP_DRAWS, len(values)))].mean(axis=1)
    return {
        "estimate": float(values.mean()),
        "ci95": [float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))],
        "n_target_clusters": int(len(values)), "draws": BOOTSTRAP_DRAWS, "seed": seed,
    }


def _contrast(per_target, comparator, candidate="correct"):
    return _bootstrap(per_target[comparator]["mse"] - per_target[candidate]["mse"])


def _positive_lcb(value: dict) -> bool:
    return value["ci95"][0] is not None and value["ci95"][0] > 0.0


def _evaluate_stratum(name, y, raw_arms, interaction_arms, indices):
    indices = np.asarray(indices, dtype=np.int64)
    if len(indices) == 0:
        return {"name": name, "n_targets": 0, "status": "EMPTY"}
    target_y = y[:, indices]
    raw = {key: value[:, indices] for key, value in raw_arms.items()}
    residual = _masked_double_center(target_y)
    interaction = {
        key: _masked_double_center(value[:, indices]) if key != "zero" else value[:, indices]
        for key, value in interaction_arms.items()
    }
    raw_summary, raw_per_target = _target_metrics(target_y, raw, interaction=False)
    int_summary, int_per_target = _target_metrics(residual, interaction, interaction=True)
    raw_contrasts = {arm: _contrast(raw_per_target, arm)
                     for arm in ("population", "ligand", "protein", "additive", "deranged")}
    int_contrasts = {arm: _contrast(int_per_target, arm) for arm in ("zero", "deranged")}
    ligand_control = _bootstrap(raw_per_target["population"]["mse"]
                                - raw_per_target["ligand"]["mse"])
    correct_corr = _bootstrap(int_per_target["correct"]["secondary"])
    enough = bool(y.shape[0] >= 100 and len(indices) >= 20)
    interaction_pass = bool(
        enough and _positive_lcb(ligand_control)
        and _positive_lcb(int_contrasts["zero"])
        and _positive_lcb(int_contrasts["deranged"])
        and _positive_lcb(correct_corr)
    )
    location_pass = bool(enough and all(_positive_lcb(value) for value in raw_contrasts.values()))
    status = (
        "INTERACTION_AND_LOCATION_OBSERVED" if interaction_pass and location_pass
        else "INTERACTION_ONLY_OBSERVED" if interaction_pass
        else "LOCATION_ONLY_OBSERVED" if location_pass
        else "SIGNAL_NOT_OBSERVED"
    )
    return {
        "name": name, "n_targets": int(len(indices)), "n_ligands": int(y.shape[0]),
        "n_finite_cells": int(np.isfinite(target_y).sum()),
        "minimum_size_precondition": enough,
        "raw_metrics": raw_summary, "interaction_metrics": int_summary,
        "raw_mse_reduction_correct_vs": raw_contrasts,
        "interaction_mse_reduction_correct_vs": int_contrasts,
        "ligand_positive_control_mse_reduction": ligand_control,
        "correct_interaction_pearson": correct_corr,
        "interaction_signal_observed": interaction_pass,
        "raw_location_signal_observed": location_pass,
        "status": status,
    }


def _fit_source_kernel(y, ligands, target_items):
    residual = double_center(y)
    ligand_bits = np.stack([record.channel_fingerprints for record in ligands])
    protein_by_target = [protein_channel_features(item["record"].pocket) for item in target_items]
    protein_features = [
        np.stack([target[channel] for target in protein_by_target])
        for channel in range(len(CHANNEL_NAMES))
    ]
    ligand_kernels = [tanimoto_kernel(ligand_bits[:, channel], ligand_bits[:, channel])
                      for channel in range(len(CHANNEL_NAMES))]
    protein_kernels, protein_tau = [], []
    for channel in range(len(CHANNEL_NAMES)):
        kernel, tau = source_rbf(protein_features[channel])
        protein_kernels.append(kernel)
        protein_tau.append(tau)

    ligand_folds = np.asarray([stable_fold(record.generic_scaffold, 3) for record in ligands])
    target_folds = np.asarray([stable_fold(item["record"].group, 3) for item in target_items])
    prediction_cache = {}
    outcome_by_fold = {}
    fold_details = []
    score = {channel: {float(value): [] for value in LAMBDA_GRID}
             for channel in range(len(CHANNEL_NAMES))}
    for fold in range(3):
        ligand_train = np.flatnonzero(ligand_folds != fold)
        ligand_valid = np.flatnonzero(ligand_folds == fold)
        target_train = np.flatnonzero(target_folds != fold)
        target_valid = np.flatnonzero(target_folds == fold)
        if min(map(len, (ligand_train, ligand_valid, target_train, target_valid))) == 0:
            continue
        train_residual = double_center(y[np.ix_(ligand_train, target_train)])
        valid_residual = double_center(y[np.ix_(ligand_valid, target_valid)])
        outcome_by_fold[fold] = valid_residual.reshape(-1)
        fold_details.append({
            "fold": fold, "train_ligands": len(ligand_train), "valid_ligands": len(ligand_valid),
            "train_targets": len(target_train), "valid_targets": len(target_valid),
            "valid_cells": int(valid_residual.size),
        })
        print(f"kernel CV fold {fold + 1}/3", flush=True)
        for channel in range(len(CHANNEL_NAMES)):
            kl_train = ligand_kernels[channel][np.ix_(ligand_train, ligand_train)]
            kp_train = protein_kernels[channel][np.ix_(target_train, target_train)]
            kl_cross = ligand_kernels[channel][np.ix_(ligand_valid, ligand_train)]
            kp_cross = protein_kernels[channel][np.ix_(target_train, target_valid)]
            for regularization in LAMBDA_GRID:
                coefficient = fit_separable_krr(
                    kl_train, kp_train, train_residual, regularization
                )
                prediction = predict_separable(kl_cross, coefficient, kp_cross)
                prediction_cache[(fold, channel, float(regularization))] = prediction.reshape(-1)
                score[channel][float(regularization)].append(
                    float(np.mean(np.square(valid_residual - prediction)))
                )
    if len(outcome_by_fold) < 2:
        raise RuntimeError("dual-cold source CV produced fewer than two folds")

    selected = []
    score_report = {}
    for channel in range(len(CHANNEL_NAMES)):
        means = {value: float(np.mean(score[channel][value])) for value in score[channel]}
        best_value = min(means.values())
        chosen = max(value for value, item in means.items()
                     if np.isclose(item, best_value, rtol=1e-12, atol=1e-14))
        selected.append(chosen)
        score_report[CHANNEL_NAMES[channel]] = {
            "selected_lambda": chosen,
            "fold_mse": score[channel],
            "mean_mse": means,
        }

    oof_y, oof_prediction = [], []
    for fold in sorted(outcome_by_fold):
        oof_y.append(outcome_by_fold[fold])
        oof_prediction.append(np.stack([
            prediction_cache[(fold, channel, float(selected[channel]))]
            for channel in range(len(CHANNEL_NAMES))
        ], axis=1))
    oof_y = np.concatenate(oof_y)
    oof_prediction = np.concatenate(oof_prediction, axis=0)
    weight_fit = fit_convex_channel_weights(oof_prediction, oof_y)
    weights = weight_fit.pop("weights")

    coefficients = []
    source_contributions = []
    print("fitting full source product kernels", flush=True)
    for channel in range(len(CHANNEL_NAMES)):
        coefficient = fit_separable_krr(
            ligand_kernels[channel], protein_kernels[channel], residual, selected[channel]
        )
        coefficients.append(coefficient)
        fitted = predict_separable(
            ligand_kernels[channel], coefficient, protein_kernels[channel]
        )
        source_contributions.append(weights[channel] * fitted)
    source_contributions = np.stack(source_contributions, axis=-1)
    bounds = fit_bounds(source_contributions)
    return {
        "ligand_bits": ligand_bits,
        "protein_features": protein_features,
        "ligand_kernels": ligand_kernels,
        "protein_kernels": protein_kernels,
        "protein_tau": protein_tau,
        "selected_lambda": selected,
        "weights": weights,
        "coefficients": coefficients,
        "bounds": bounds,
        "report": {
            "lambda_grid": LAMBDA_GRID,
            "folds": fold_details,
            "channels": score_report,
            "convex_weight_fit": {**weight_fit, "weights": weights.tolist(),
                                  "channel_names": CHANNEL_NAMES},
            "protein_rbf_tau": dict(zip(CHANNEL_NAMES, protein_tau)),
            "source_channel_bounds": bounds,
        },
    }


def _panel_from_pkis2(frame, smiles_path, klifs_index, source_genes, source_smiles, source_scaffolds):
    targets_all, target_exclusions = _target_mapping(frame.columns, klifs_index)
    targets = [item for item in targets_all if item["record"].hgnc not in source_genes]
    smiles = load_smiles(smiles_path)
    all_ligands, ligand_exclusions = _ligand_records(
        [(str(molid), smiles[str(molid)]) for molid in frame.index.astype(str)]
    )
    ligands = [record for record in all_ligands
               if record.canonical_smiles not in source_smiles
               and record.generic_scaffold not in source_scaffolds]
    y = frame.loc[
        [record.molid for record in ligands],
        [item["assay_target"] for item in targets],
    ].to_numpy(dtype=np.float64)
    return {
        "name": "PKIS2", "role": "consumed_development_transfer",
        "ligands": ligands, "targets": targets, "y": y,
        "target_exclusions": target_exclusions,
        "ligand_exclusions": ligand_exclusions,
        "n_not_scaffold_cold": len(all_ligands) - len(ligands),
        "assay_transform": "identity_continuous_activity",
    }


def _panel_from_anastassiadis(
    workbook, identity_path, klifs_index, source_genes, source_smiles, source_scaffolds,
):
    identities = json.loads(Path(identity_path).read_text(encoding="utf-8"))
    resolved = [entry for entry in identities["entries"] if entry["status"] == "resolved"]
    all_ligands, ligand_exclusions = _ligand_records([
        (str(entry["workbook_column"]), entry["smiles"]) for entry in resolved
    ])
    ligands = [record for record in all_ligands
               if record.canonical_smiles not in source_smiles
               and record.generic_scaffold not in source_scaffolds]

    raw = pd.read_excel(workbook, sheet_name=0, header=None, skiprows=3)
    target_names = raw.iloc[:, 0].dropna().astype(str).tolist()
    targets_all, target_exclusions = _target_mapping(target_names, klifs_index)
    targets = [item for item in targets_all if item["record"].hgnc not in source_genes]
    row_by_target = {str(raw.iat[row, 0]): row for row in range(len(raw))
                     if not pd.isna(raw.iat[row, 0])}
    remaining = np.empty((len(ligands), len(targets)), dtype=np.float64)
    for target_index, item in enumerate(targets):
        row = row_by_target[item["assay_target"]]
        for ligand_index, ligand in enumerate(ligands):
            remaining[ligand_index, target_index] = pd.to_numeric(
                raw.iat[row, int(ligand.molid)], errors="coerce"
            )
    y = np.clip(1.0 - remaining / 100.0, 0.0, 1.0)
    unresolved = [entry for entry in identities["entries"] if entry["status"] != "resolved"]
    return {
        "name": "Anastassiadis2011", "role": "external_cross_assay_transfer",
        "ligands": ligands, "targets": targets, "y": y,
        "target_exclusions": target_exclusions,
        "ligand_exclusions": ligand_exclusions + [
            {"molid": str(entry["workbook_column"]), "reason": "pubchem_unresolved",
             "detail": entry.get("error", "")}
            for entry in unresolved
        ],
        "n_not_scaffold_cold": len(all_ligands) - len(ligands),
        "assay_transform": "clip(1 - percent_remaining_activity/100, 0, 1)",
        "identity_sidecar_sha256": sha256_file(identity_path),
    }


def _evaluate_panel(panel, source, source_targets, nuisance, output):
    ligands, targets, y = panel["ligands"], panel["targets"], panel["y"]
    if len(ligands) < 100 or len(targets) < 20:
        return {
            "schema": "MetaSieve.ProductKernelTransfer.v2", "panel": panel["name"],
            "verdict": "NOT_RUN_SAMPLE_PRECONDITION", "n_ligands": len(ligands),
            "n_targets": len(targets), "admission_to_biological_z_authorized": False,
        }
    ligand_bits = np.stack([record.channel_fingerprints for record in ligands])
    protein_by_target = [protein_channel_features(item["record"].pocket) for item in targets]
    protein_features = [
        np.stack([target[channel] for target in protein_by_target])
        for channel in range(len(CHANNEL_NAMES))
    ]
    contributions = []
    ligand_cross_kernels = []
    protein_cross_kernels = []
    for channel in range(len(CHANNEL_NAMES)):
        ligand_cross = tanimoto_kernel(ligand_bits[:, channel], source["ligand_bits"][:, channel])
        protein_cross = cross_rbf(
            source["protein_features"][channel], protein_features[channel],
            source["protein_tau"][channel],
        )
        predicted = predict_separable(
            ligand_cross, source["coefficients"][channel], protein_cross
        )
        contributions.append(source["weights"][channel] * predicted)
        ligand_cross_kernels.append(ligand_cross)
        protein_cross_kernels.append(protein_cross)
    contributions = np.stack(contributions, axis=-1)
    interaction_prediction = contributions.sum(axis=-1)
    candidate_z = bounded_z(contributions, source["bounds"])

    target_records = [item["record"] for item in targets]
    derangement = deterministic_derangement(target_records)
    deranged_contributions = contributions[:, derangement, :]
    deranged_interaction = deranged_contributions.sum(axis=-1)

    x_ligand = np.stack([record.nuisance_features for record in ligands])
    x_protein = np.stack([protein_model_features(item["record"].pocket) for item in targets])
    ligand_prediction = nuisance["ligand_model"].predict(x_ligand)
    protein_prediction = nuisance["protein_model"].predict(x_protein)
    mu = nuisance["mu"]
    additive = mu + ligand_prediction[:, None] + protein_prediction[None, :]
    raw_arms = {
        "population": np.full_like(y, mu),
        "ligand": np.broadcast_to(mu + ligand_prediction[:, None], y.shape),
        "protein": np.broadcast_to(mu + protein_prediction[None, :], y.shape),
        "additive": additive,
        "correct": additive + interaction_prediction,
        "deranged": additive + deranged_interaction,
    }
    raw_arms = {name: np.clip(value, 0.0, 1.0) for name, value in raw_arms.items()}
    interaction_arms = {
        "zero": np.zeros_like(y),
        "correct": interaction_prediction,
        "deranged": deranged_interaction,
    }

    source_families = {item["record"].family for item in source_targets}
    source_groups = {item["record"].group for item in source_targets}
    family_indices = np.asarray([
        index for index, item in enumerate(targets)
        if item["record"].family not in source_families
    ], dtype=np.int64)
    group_indices = np.asarray([
        index for index, item in enumerate(targets)
        if item["record"].group not in source_groups
    ], dtype=np.int64)
    main = _evaluate_stratum(
        "exact_target_cold_and_scaffold_cold", y, raw_arms, interaction_arms,
        np.arange(len(targets)),
    )
    family = _evaluate_stratum(
        "klifs_family_cold_and_scaffold_cold", y, raw_arms, interaction_arms, family_indices,
    )
    group = _evaluate_stratum(
        "klifs_group_cold_and_scaffold_cold", y, raw_arms, interaction_arms, group_indices,
    )

    steric_ligand_similarity = ligand_cross_kernels[-1].max(axis=1)
    mean_protein_kernel = np.mean(np.stack(protein_cross_kernels, axis=0), axis=0)
    protein_similarity = mean_protein_kernel.max(axis=0)
    confidence = np.sqrt(np.clip(
        steric_ligand_similarity[:, None] * protein_similarity[None, :], 0.0, 1.0
    ))
    mapping = []
    for index, item in enumerate(targets):
        mapping.append({
            "correct": item["record"].hgnc,
            "deranged": targets[int(derangement[index])]["record"].hgnc,
            "correct_group": item["record"].group,
            "deranged_group": targets[int(derangement[index])]["record"].group,
        })
    _write_json(output / f"{panel['name']}_derangement.json", {
        "mapping": mapping,
        "fixed_points": int(np.sum(derangement == np.arange(len(targets)))),
        "same_group_fraction": float(np.mean([
            targets[index]["record"].group == targets[int(derangement[index])]["record"].group
            for index in range(len(targets))
        ])),
    })
    with (output / f"{panel['name']}_target_mapping.jsonl").open("w", encoding="utf-8") as handle:
        for item in targets:
            handle.write(json.dumps({"assay_target": item["assay_target"],
                                     **asdict(item["record"])}, sort_keys=True) + "\n")
    _write_json(output / f"{panel['name']}_exclusions.json", {
        "target": panel["target_exclusions"], "ligand": panel["ligand_exclusions"],
        "ligands_not_strict_scaffold_cold": panel["n_not_scaffold_cold"],
    })
    verdict = (
        "CANDIDATE_INTERACTION_AND_LOCATION_VALIDATED" if
        main["interaction_signal_observed"] and main["raw_location_signal_observed"]
        else "CANDIDATE_INTERACTION_VALIDATED" if main["interaction_signal_observed"]
        else "CANDIDATE_LOCATION_VALIDATED" if main["raw_location_signal_observed"]
        else "REVISION_V2_NOT_VALIDATED"
    )
    return {
        "schema": "MetaSieve.ProductKernelTransfer.v2",
        "panel": panel["name"], "role": panel["role"], "verdict": verdict,
        "formal_gate": False, "assay_transform": panel["assay_transform"],
        "sample": {
            "n_ligands": len(ligands), "n_targets": len(targets),
            "n_family_cold_targets": int(len(family_indices)),
            "n_group_cold_targets": int(len(group_indices)),
            "n_finite_cells": int(np.isfinite(y).sum()),
        },
        "main_stratum": main, "family_cold_stratum": family, "group_cold_stratum": group,
        "candidate_biological_z": {
            "coordinate_names": CHANNEL_NAMES, "shape": list(candidate_z.shape),
            "minimum_by_coordinate": np.nanmin(candidate_z.reshape(-1, 5), axis=0).tolist(),
            "maximum_by_coordinate": np.nanmax(candidate_z.reshape(-1, 5), axis=0).tolist(),
            "mean_by_coordinate": np.nanmean(candidate_z.reshape(-1, 5), axis=0).tolist(),
            "bounded": bool(np.all((candidate_z >= 0.0) & (candidate_z <= 1.0))),
            "admitted_to_model_config": False,
        },
        "coverage": {
            "nearest_source_steric_tanimoto": {
                "min": float(steric_ligand_similarity.min()),
                "median": float(np.median(steric_ligand_similarity)),
                "max": float(steric_ligand_similarity.max()),
            },
            "nearest_source_mean_pocket_kernel": {
                "min": float(protein_similarity.min()),
                "median": float(np.median(protein_similarity)),
                "max": float(protein_similarity.max()),
            },
            "pair_confidence": {
                "min": float(confidence.min()), "median": float(np.median(confidence)),
                "max": float(confidence.max()),
            },
            "abstention_required": True,
        },
        "law_interface_contract": {
            "candidate_z_dimension": 5,
            "candidate_z_bounded": True,
            "frozen_operator_modified": False,
            "direct_scalar_output_used_for_deployment": False,
            "diagnostic_raw_activity_scored": True,
            "admission_to_biological_z_authorized": False,
        },
    }


def run(args):
    started = time.time()
    repo_root = Path(__file__).resolve().parents[2]
    informers = Path(args.informers_root).resolve()
    klifs_path = Path(args.klifs_json).resolve()
    workbook = Path(args.anastassiadis_workbook).resolve()
    identities = Path(args.anastassiadis_identities).resolve()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    paths = {
        "pkis1_labels": informers / "data/pkis1_continuous_labels.csv",
        "pkis2_labels": informers / "data/pkis2_continuous_labels.csv",
        "pkis1_smiles": informers / "data/compounds/pkis1_uniq_tested.can",
        "pkis2_smiles": informers / "data/compounds/pkis2_uniq_tested.can",
        "klifs": klifs_path, "anastassiadis_workbook": workbook,
        "anastassiadis_identities": identities,
        "preregistration": Path(__file__).with_name("REVISION_V2_PREREGISTRATION.md"),
    }
    for path in paths.values():
        if not path.exists():
            raise FileNotFoundError(path)

    pkis1 = pd.read_csv(paths["pkis1_labels"], dtype={"molid": str}).set_index("molid")
    pkis2 = pd.read_csv(paths["pkis2_labels"], dtype={"molid": str}).set_index("molid")
    klifs_index, _ = load_klifs(klifs_path)
    source_targets, source_target_exclusions = _target_mapping(pkis1.columns, klifs_index)
    source_smiles_map = load_smiles(paths["pkis1_smiles"])
    source_ligands, source_ligand_exclusions = _ligand_records([
        (str(molid), source_smiles_map[str(molid)]) for molid in pkis1.index.astype(str)
    ])
    source_ids = [record.molid for record in source_ligands]
    source_columns = [item["assay_target"] for item in source_targets]
    y_source = pkis1.loc[source_ids, source_columns].to_numpy(dtype=np.float64)
    source_genes = {item["record"].hgnc for item in source_targets}
    source_smiles = {record.canonical_smiles for record in source_ligands}
    source_scaffolds = {record.generic_scaffold for record in source_ligands}

    print("selecting the frozen five-channel source estimator", flush=True)
    source = _fit_source_kernel(y_source, source_ligands, source_targets)
    mu = float(y_source.mean())
    ligand_effect = y_source.mean(axis=1) - mu
    protein_effect = y_source.mean(axis=0) - mu
    x_ligand = np.stack([record.nuisance_features for record in source_ligands])
    x_protein = np.stack([protein_model_features(item["record"].pocket) for item in source_targets])
    ligand_alpha, ligand_cv = _grouped_alpha(
        x_ligand, ligand_effect, [record.generic_scaffold for record in source_ligands]
    )
    protein_alpha, protein_cv = _grouped_alpha(
        x_protein, protein_effect, [item["record"].group for item in source_targets]
    )
    ligand_model = make_pipeline(StandardScaler(), Ridge(alpha=ligand_alpha, solver="lsqr", tol=1e-7))
    protein_model = make_pipeline(StandardScaler(), Ridge(alpha=protein_alpha, solver="lsqr", tol=1e-7))
    ligand_model.fit(x_ligand, ligand_effect)
    protein_model.fit(x_protein, protein_effect)
    nuisance = {"mu": mu, "ligand_model": ligand_model, "protein_model": protein_model}

    panels = [
        _panel_from_pkis2(
            pkis2, paths["pkis2_smiles"], klifs_index, source_genes,
            source_smiles, source_scaffolds,
        ),
        _panel_from_anastassiadis(
            workbook, identities, klifs_index, source_genes, source_smiles, source_scaffolds,
        ),
    ]
    results = {}
    for panel in panels:
        print(f"evaluating frozen transfer: {panel['name']}", flush=True)
        results[panel["name"]] = _evaluate_panel(
            panel, source, source_targets, nuisance, output
        )

    external = results["Anastassiadis2011"]
    external_pass = bool(
        external.get("main_stratum", {}).get("interaction_signal_observed", False)
        or external.get("main_stratum", {}).get("raw_location_signal_observed", False)
    )
    overall = {
        "schema": "MetaSieve.ProductKernelRevision.v2",
        "formal_gate": False,
        "main_verdict": (
            "EXTERNAL_CANDIDATE_SIGNAL_OBSERVED_FORMAL_GATE_STILL_REQUIRED"
            if external_pass else "REVISION_V2_NOT_VALIDATED"
        ),
        "end_to_end_dta_validated": False,
        "admission_to_biological_z_authorized": False,
        "source": {
            "panel": "PKIS1", "n_ligands": len(source_ligands),
            "n_targets": len(source_targets), "mean_activity": mu,
            "model_selection": source["report"],
            "nuisance_model_selection": {
                "ligand": {"selected_alpha": ligand_alpha, **ligand_cv},
                "protein": {"selected_alpha": protein_alpha, **protein_cv},
            },
            "target_exclusions": source_target_exclusions,
            "ligand_exclusions": source_ligand_exclusions,
        },
        "transfers": results,
        "read_firewall": {
            "chembl_x1_labels": 0, "davis_labels": 0, "recipient_labels": 0,
            "pkis1_role": "source_training", "pkis2_role": "consumed_development_transfer",
            "anastassiadis_role": "external_cross_assay_transfer",
        },
        "revision_budget": {"permitted_representation_revisions": 1, "used": 1, "remaining": 0},
        "limitations": [
            "PKIS and Anastassiadis are single-concentration inhibition panels, not Ki/Kd affinity.",
            "The external validation remains kinase-focused; non-kinase protein generalization is untested.",
            "The five channels are transfer coordinates, not identified free-energy terms.",
            "The Anastassiadis workbook schema and a small header preview were inspected before freeze; no outcome statistic or fit was computed.",
            "This exploratory run cannot change X1 authorization or admit biological z.",
        ],
    }
    _write_json(output / "result.json", overall)
    _write_json(output / "source_model_selection.json", source["report"])

    manifest = {
        "schema": "MetaSieve.ProductKernelRevisionManifest.v2",
        "created_unix": time.time(), "runtime_seconds": time.time() - started,
        "seed": SEED, "bootstrap_draws": BOOTSTRAP_DRAWS,
        "git_head": _git_value(repo_root, "rev-parse", "HEAD"),
        "git_branch": _git_value(repo_root, "branch", "--show-current"),
        "input_sha256": {name: sha256_file(path) for name, path in paths.items()},
        "code_sha256": {
            "run_kernel_revision.py": sha256_file(__file__),
            "kernel_revision.py": sha256_file(Path(__file__).with_name("kernel_revision.py")),
            "mechanism.py": sha256_file(Path(__file__).with_name("mechanism.py")),
        },
        "outputs": {},
    }
    for path in sorted(output.iterdir()):
        if path.name != "manifest.json" and path.is_file():
            manifest["outputs"][path.name] = sha256_file(path)
    _write_json(output / "manifest.json", manifest)
    print(json.dumps({
        "verdict": overall["main_verdict"],
        "pkis2": results["PKIS2"].get("verdict"),
        "anastassiadis": external.get("verdict"),
        "runtime_seconds": manifest["runtime_seconds"],
    }, indent=2), flush=True)
    return overall


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--informers-root", required=True)
    parser.add_argument("--klifs-json", required=True)
    parser.add_argument("--anastassiadis-workbook", required=True)
    parser.add_argument("--anastassiadis-identities", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args(argv)


if __name__ == "__main__":
    run(parse_args())
