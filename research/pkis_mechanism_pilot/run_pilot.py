"""Run the preregistered exploratory PKIS1 -> PKIS2 mechanism pilot.

The script is CPU-only and writes a fully auditable result directory. It never
imports the GPU frontend, reads ChEMBL/DAVIS/recipient labels, or mutates
``model/``.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import gc
import json
from pathlib import Path
import subprocess
import time
import warnings

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from .mechanism import (
    CHANNEL_NAMES,
    SHELL_NAMES,
    bounded_biological_z,
    channel_contributions,
    deterministic_derangement,
    double_center,
    fit_channel_bounds,
    ligand_from_smiles,
    load_klifs,
    load_smiles,
    map_target,
    normalize_name,
    ordered_anchor_simplex,
    pair_feature_matrix,
    pocket_sitealign,
    protein_model_features,
    protein_pair_properties,
    sha256_file,
    stable_fold,
)


SEED = 20260808
ALPHAS = (0.01, 0.1, 1.0, 10.0, 100.0, 1000.0)
BOOTSTRAP_DRAWS = 10_000


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


def _git_value(root: Path, *arguments: str) -> str:
    try:
        return subprocess.check_output(["git", *arguments], cwd=root, text=True).strip()
    except Exception:
        return "unavailable"


def _target_mapping(columns, klifs_index):
    included = []
    excluded = []
    for column in columns:
        record, reason = map_target(column, klifs_index)
        if record is None:
            excluded.append({"assay_target": column, "reason": reason})
        else:
            included.append({"assay_target": column, "record": record})

    # A gene may occur under a modern HGNC name and an old synonym (for example
    # MYLK and MLCK). Prefer the literal HGNC/name match and exclude the alias so
    # one protein never contributes two assay contexts to this proxy panel.
    by_gene = {}
    for item in included:
        by_gene.setdefault(item["record"].hgnc, []).append(item)
    resolved = []
    for gene, items in by_gene.items():
        if len(items) == 1:
            resolved.append(items[0])
            continue
        direct = [item for item in items if normalize_name(item["assay_target"])
                  in {normalize_name(item["record"].hgnc), normalize_name(item["record"].name)}]
        if len(direct) == 1:
            resolved.append(direct[0])
            rejected = [item for item in items if item is not direct[0]]
        else:
            rejected = items
        for item in rejected:
            excluded.append({
                "assay_target": item["assay_target"], "mapped_hgnc": gene,
                "reason": "duplicate_hgnc_assay_context",
            })
    resolved.sort(key=lambda item: item["assay_target"])
    return resolved, excluded


def _ligand_records(frame: pd.DataFrame, smiles: dict[str, str]):
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures

    factory = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    records = []
    excluded = []
    for index, molid in enumerate(frame.index.astype(str)):
        try:
            records.append(ligand_from_smiles(molid, smiles[molid], factory=factory))
        except Exception as error:
            excluded.append({"molid": molid, "reason": type(error).__name__, "detail": str(error)})
        if (index + 1) % 100 == 0:
            print(f"ligands: {index + 1}/{len(frame)}", flush=True)
    return records, excluded


def _grouped_alpha(X, y, groups, alphas=ALPHAS):
    groups = np.asarray(groups)
    unique = np.unique(groups)
    splits = min(5, len(unique))
    if splits < 2:
        raise RuntimeError("source-only alpha selection needs at least two groups")
    cv = GroupKFold(n_splits=splits)
    scores = {float(alpha): [] for alpha in alphas}
    for train, valid in cv.split(X, y, groups):
        for alpha in alphas:
            model = make_pipeline(
                StandardScaler(), Ridge(alpha=float(alpha), solver="lsqr", tol=1e-7)
            )
            model.fit(X[train], y[train])
            prediction = model.predict(X[valid])
            scores[float(alpha)].append(float(np.mean(np.square(y[valid] - prediction))))
    means = {alpha: float(np.mean(values)) for alpha, values in scores.items()}
    best_value = min(means.values())
    best = max(alpha for alpha, value in means.items() if np.isclose(value, best_value, rtol=1e-12, atol=1e-14))
    return best, {"fold_mse": scores, "mean_mse": means}


def _pair_alpha(features, residual, ligand_folds, target_folds, alphas=ALPHAS):
    flat = features.reshape(-1, features.shape[-1])
    y = residual.reshape(-1)
    ligand_grid = np.repeat(np.asarray(ligand_folds), residual.shape[1])
    target_grid = np.tile(np.asarray(target_folds), residual.shape[0])
    scores = {float(alpha): [] for alpha in alphas}
    fold_counts = []
    for fold in range(3):
        train = (ligand_grid != fold) & (target_grid != fold)
        valid = (ligand_grid == fold) & (target_grid == fold)
        if train.sum() == 0 or valid.sum() == 0:
            continue
        fold_counts.append({"fold": fold, "train": int(train.sum()), "valid": int(valid.sum())})
        X_train = flat[train]
        X_valid = flat[valid]
        for alpha in alphas:
            model = Ridge(alpha=float(alpha), solver="lsqr", tol=1e-6, max_iter=5000)
            model.fit(X_train, y[train])
            prediction = model.predict(X_valid)
            scores[float(alpha)].append(float(np.mean(np.square(y[valid] - prediction))))
        del X_train, X_valid
        gc.collect()
    if len(fold_counts) < 2:
        raise RuntimeError("joint protein-group/ligand-scaffold source CV produced fewer than two folds")
    means = {alpha: float(np.mean(values)) for alpha, values in scores.items()}
    best_value = min(means.values())
    best = max(alpha for alpha, value in means.items() if np.isclose(value, best_value, rtol=1e-12, atol=1e-14))
    return best, {"fold_counts": fold_counts, "fold_mse": scores, "mean_mse": means}


def _safe_spearman(y, prediction):
    if np.std(y) < 1e-12 or np.std(prediction) < 1e-12:
        return 0.0
    return float(spearmanr(y, prediction).statistic)


def _safe_pearson(y, prediction):
    if np.std(y) < 1e-12 or np.std(prediction) < 1e-12:
        return 0.0
    return float(pearsonr(y, prediction).statistic)


def _raw_target_metrics(y, arms):
    out = {}
    per_target = {}
    for name, prediction in arms.items():
        error = y - prediction
        target_mse = np.mean(np.square(error), axis=0)
        target_mae = np.mean(np.abs(error), axis=0)
        target_spearman = np.asarray([
            _safe_spearman(y[:, index], prediction[:, index]) for index in range(y.shape[1])
        ])
        per_target[name] = {"mse": target_mse, "mae": target_mae, "spearman": target_spearman}
        out[name] = {
            "target_macro_mse": float(target_mse.mean()),
            "target_macro_mae": float(target_mae.mean()),
            "target_macro_spearman": float(target_spearman.mean()),
        }
    return out, per_target


def _interaction_target_metrics(residual, arms):
    out = {}
    per_target = {}
    for name, prediction in arms.items():
        target_mse = np.mean(np.square(residual - prediction), axis=0)
        target_corr = np.asarray([
            _safe_pearson(residual[:, index], prediction[:, index])
            for index in range(residual.shape[1])
        ])
        per_target[name] = {"mse": target_mse, "pearson": target_corr}
        out[name] = {
            "target_macro_mse": float(target_mse.mean()),
            "target_macro_pearson": float(target_corr.mean()),
        }
    return out, per_target


def _bootstrap(values, draws=BOOTSTRAP_DRAWS, seed=SEED):
    values = np.asarray(values, dtype=np.float64)
    rng = np.random.default_rng(seed)
    samples = values[rng.integers(0, len(values), size=(draws, len(values)))].mean(axis=1)
    return {
        "estimate": float(values.mean()),
        "ci95": [float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))],
        "n_target_clusters": int(len(values)),
        "draws": int(draws),
        "seed": int(seed),
    }


def _contrast(per_target, comparator, candidate="correct"):
    # Positive means candidate has lower MSE.
    return _bootstrap(per_target[comparator]["mse"] - per_target[candidate]["mse"])


def _spectral_summary(matrix):
    centered = double_center(matrix)
    singular = np.linalg.svd(centered, compute_uv=False)
    variance = np.square(singular)
    total = max(float(variance.sum()), 1e-15)
    return {
        "interaction_variance_fraction": float(np.square(centered).sum()
                                                / max(np.square(matrix - matrix.mean()).sum(), 1e-15)),
        "top_singular_values": singular[:10].tolist(),
        "top_1_variance_share": float(variance[:1].sum() / total),
        "top_5_variance_share": float(variance[:5].sum() / total),
        "top_10_variance_share": float(variance[:10].sum() / total),
    }


def _nearest_coverage(source_ligand, test_ligand, source_pockets, test_pockets):
    source_bits = source_ligand[:, :1024] > 0.5
    test_bits = test_ligand[:, :1024] > 0.5
    intersection = test_bits.astype(np.float32) @ source_bits.astype(np.float32).T
    union = (test_bits.sum(axis=1, keepdims=True) + source_bits.sum(axis=1)[None, :] - intersection)
    ligand_similarity = np.max(intersection / np.maximum(union, 1.0), axis=1)
    protein_identity = np.asarray([
        max(sum(a == b for a, b in zip(query, source)) / 85.0 for source in source_pockets)
        for query in test_pockets
    ])
    confidence = np.sqrt(np.clip(ligand_similarity[:, None] * protein_identity[None, :], 0.0, 1.0))
    return ligand_similarity, protein_identity, confidence


def _evaluate_stratum(name, y, residual, raw_arms, interaction_arms, target_indices):
    target_indices = np.asarray(target_indices, dtype=np.int64)
    if len(target_indices) == 0:
        return {"name": name, "n_targets": 0, "status": "EMPTY"}
    raw = {key: value[:, target_indices] for key, value in raw_arms.items()}
    interaction = {key: value[:, target_indices] for key, value in interaction_arms.items()}
    raw_summary, raw_per_target = _raw_target_metrics(y[:, target_indices], raw)
    int_summary, int_per_target = _interaction_target_metrics(residual[:, target_indices], interaction)
    raw_contrasts = {comparator: _contrast(raw_per_target, comparator)
                     for comparator in ("population", "ligand", "protein", "additive", "deranged")}
    int_contrasts = {comparator: _contrast(int_per_target, comparator)
                     for comparator in ("zero", "deranged")}
    ligand_control = _bootstrap(raw_per_target["population"]["mse"]
                                - raw_per_target["ligand"]["mse"])
    correct_corr = _bootstrap(int_per_target["correct"]["pearson"])
    enough = len(target_indices) >= 20 and y.shape[0] >= 100
    interaction_pass = bool(
        enough
        and int_contrasts["zero"]["ci95"][0] > 0.0
        and int_contrasts["deranged"]["ci95"][0] > 0.0
        and correct_corr["ci95"][0] > 0.0
        and ligand_control["ci95"][0] > 0.0
    )
    location_pass = bool(
        enough and all(contrast["ci95"][0] > 0.0 for contrast in raw_contrasts.values())
    )
    return {
        "name": name, "n_targets": int(len(target_indices)), "n_ligands": int(y.shape[0]),
        "raw_metrics": raw_summary, "interaction_metrics": int_summary,
        "raw_mse_reduction_correct_vs": raw_contrasts,
        "interaction_mse_reduction_correct_vs": int_contrasts,
        "ligand_positive_control_mse_reduction": ligand_control,
        "correct_interaction_pearson": correct_corr,
        "minimum_size_precondition": enough,
        "interaction_signal_observed": interaction_pass,
        "raw_location_signal_observed": location_pass,
        "status": "INTERACTION_AND_LOCATION_OBSERVED" if interaction_pass and location_pass
                  else "INTERACTION_ONLY_OBSERVED" if interaction_pass
                  else "LOCATION_ONLY_OBSERVED" if location_pass
                  else "SIGNAL_NOT_OBSERVED",
    }


def run(args):
    started = time.time()
    repo_root = Path(__file__).resolve().parents[2]
    informers = Path(args.informers_root).resolve()
    klifs_path = Path(args.klifs_json).resolve()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)

    paths = {
        "pkis1_labels": informers / "data/pkis1_continuous_labels.csv",
        "pkis2_labels": informers / "data/pkis2_continuous_labels.csv",
        "pkis1_smiles": informers / "data/compounds/pkis1_uniq_tested.can",
        "pkis2_smiles": informers / "data/compounds/pkis2_uniq_tested.can",
        "klifs": klifs_path,
        "preregistration": Path(__file__).with_name("PREREGISTRATION.md"),
    }
    for path in paths.values():
        if not path.exists():
            raise FileNotFoundError(path)

    pkis1 = pd.read_csv(paths["pkis1_labels"], dtype={"molid": str}).set_index("molid")
    pkis2 = pd.read_csv(paths["pkis2_labels"], dtype={"molid": str}).set_index("molid")
    klifs_index, _ = load_klifs(klifs_path)
    source_targets, source_target_exclusions = _target_mapping(pkis1.columns, klifs_index)
    transfer_targets, transfer_target_exclusions = _target_mapping(pkis2.columns, klifs_index)

    source_genes = {item["record"].hgnc for item in source_targets}
    cold_targets = [item for item in transfer_targets if item["record"].hgnc not in source_genes]
    source_families = {item["record"].family for item in source_targets}
    family_cold = np.asarray([
        index for index, item in enumerate(cold_targets)
        if item["record"].family not in source_families
    ], dtype=np.int64)
    group_cold = np.asarray([
        index for index, item in enumerate(cold_targets)
        if item["record"].group not in {source["record"].group for source in source_targets}
    ], dtype=np.int64)

    source_ligands, source_ligand_exclusions = _ligand_records(
        pkis1, load_smiles(paths["pkis1_smiles"])
    )
    transfer_ligands_all, transfer_ligand_exclusions = _ligand_records(
        pkis2, load_smiles(paths["pkis2_smiles"])
    )
    source_scaffolds = {record.generic_scaffold for record in source_ligands}
    source_smiles = {record.canonical_smiles for record in source_ligands}
    transfer_ligands = [record for record in transfer_ligands_all
                        if record.generic_scaffold not in source_scaffolds
                        and record.canonical_smiles not in source_smiles]

    if len(cold_targets) < 20 or len(transfer_ligands) < 100:
        raise RuntimeError("registered external sample-size precondition failed")

    source_columns = [item["assay_target"] for item in source_targets]
    transfer_columns = [item["assay_target"] for item in cold_targets]
    source_ids = [record.molid for record in source_ligands]
    transfer_ids = [record.molid for record in transfer_ligands]
    Y1 = pkis1.loc[source_ids, source_columns].to_numpy(dtype=np.float64)
    Y2 = pkis2.loc[transfer_ids, transfer_columns].to_numpy(dtype=np.float64)

    X_ligand_1 = np.stack([record.model_features for record in source_ligands])
    X_ligand_2 = np.stack([record.model_features for record in transfer_ligands])
    X_protein_1 = np.stack([protein_model_features(item["record"].pocket)
                            for item in source_targets])
    X_protein_2 = np.stack([protein_model_features(item["record"].pocket)
                            for item in cold_targets])
    P_pair_1 = np.stack([protein_pair_properties(item["record"].pocket)
                         for item in source_targets])
    P_pair_2 = np.stack([protein_pair_properties(item["record"].pocket)
                         for item in cold_targets])
    L_pair_1 = np.stack([record.pharmacophore_shells for record in source_ligands])
    L_pair_2 = np.stack([record.pharmacophore_shells for record in transfer_ligands])

    mu = float(Y1.mean())
    ligand_effect = Y1.mean(axis=1) - mu
    protein_effect = Y1.mean(axis=0) - mu
    source_residual = double_center(Y1)
    transfer_residual = double_center(Y2)

    ligand_alpha, ligand_cv = _grouped_alpha(
        X_ligand_1, ligand_effect,
        [record.generic_scaffold for record in source_ligands],
    )
    protein_alpha, protein_cv = _grouped_alpha(
        X_protein_1, protein_effect,
        [item["record"].group for item in source_targets],
    )
    ligand_model = make_pipeline(StandardScaler(), Ridge(alpha=ligand_alpha, solver="lsqr", tol=1e-7))
    protein_model = make_pipeline(StandardScaler(), Ridge(alpha=protein_alpha, solver="lsqr", tol=1e-7))
    ligand_model.fit(X_ligand_1, ligand_effect)
    protein_model.fit(X_protein_1, protein_effect)
    ligand_prediction = ligand_model.predict(X_ligand_2)
    protein_prediction = protein_model.predict(X_protein_2)

    print("building source typed tensor", flush=True)
    source_pair_features = pair_feature_matrix(L_pair_1, P_pair_1)
    ligand_folds = [stable_fold(record.generic_scaffold, 3) for record in source_ligands]
    target_folds = [stable_fold(item["record"].group, 3) for item in source_targets]
    pair_alpha, pair_cv = _pair_alpha(
        source_pair_features, source_residual, ligand_folds, target_folds
    )
    pair_model = Ridge(alpha=pair_alpha, solver="lsqr", tol=1e-6, max_iter=5000)
    pair_model.fit(source_pair_features.reshape(-1, source_pair_features.shape[-1]),
                   source_residual.reshape(-1))
    source_channel = channel_contributions(source_pair_features, pair_model.coef_)
    channel_bounds = fit_channel_bounds(source_channel)
    del source_pair_features
    gc.collect()

    print("building transfer typed tensor", flush=True)
    transfer_pair_features = pair_feature_matrix(L_pair_2, P_pair_2)
    correct_interaction_raw = pair_model.predict(
        transfer_pair_features.reshape(-1, transfer_pair_features.shape[-1])
    ).reshape(Y2.shape)
    transfer_channel = channel_contributions(transfer_pair_features, pair_model.coef_)
    biological_z = bounded_biological_z(transfer_channel, channel_bounds)
    del transfer_pair_features
    gc.collect()

    records = [item["record"] for item in cold_targets]
    derangement = deterministic_derangement(records)
    print("building deranged transfer tensor", flush=True)
    deranged_pair_features = pair_feature_matrix(L_pair_2, P_pair_2[derangement])
    deranged_interaction_raw = pair_model.predict(
        deranged_pair_features.reshape(-1, deranged_pair_features.shape[-1])
    ).reshape(Y2.shape)
    del deranged_pair_features
    gc.collect()

    additive = mu + ligand_prediction[:, None] + protein_prediction[None, :]
    raw_arms = {
        "population": np.full_like(Y2, mu),
        "ligand": np.broadcast_to(mu + ligand_prediction[:, None], Y2.shape),
        "protein": np.broadcast_to(mu + protein_prediction[None, :], Y2.shape),
        "additive": additive,
        "correct": additive + correct_interaction_raw,
        "deranged": additive + deranged_interaction_raw,
    }
    raw_arms = {name: np.clip(value, 0.0, 1.0) for name, value in raw_arms.items()}
    interaction_arms = {
        "zero": np.zeros_like(transfer_residual),
        "correct": double_center(correct_interaction_raw),
        "deranged": double_center(deranged_interaction_raw),
    }

    main = _evaluate_stratum(
        "exact_target_cold_and_scaffold_cold", Y2, transfer_residual,
        raw_arms, interaction_arms, np.arange(Y2.shape[1]),
    )
    family = _evaluate_stratum(
        "klifs_family_cold_and_scaffold_cold", Y2, transfer_residual,
        raw_arms, interaction_arms, family_cold,
    )
    group = _evaluate_stratum(
        "klifs_group_cold_and_scaffold_cold", Y2, transfer_residual,
        raw_arms, interaction_arms, group_cold,
    )

    group_effects = {}
    for group_name in sorted({record.group for record in records}):
        indices = np.asarray([index for index, record in enumerate(records) if record.group == group_name])
        summary, per_target = _interaction_target_metrics(
            transfer_residual[:, indices],
            {name: value[:, indices] for name, value in interaction_arms.items()},
        )
        group_effects[group_name] = {
            "n_targets": int(len(indices)), "metrics": summary,
            "correct_vs_zero_mse_reduction": _contrast(per_target, "zero"),
            "correct_vs_deranged_mse_reduction": _contrast(per_target, "deranged"),
        }

    source_pockets = [item["record"].pocket for item in source_targets]
    test_pockets = [item["record"].pocket for item in cold_targets]
    ligand_similarity, protein_identity, confidence = _nearest_coverage(
        X_ligand_1, X_ligand_2, source_pockets, test_pockets
    )
    simplex = ordered_anchor_simplex(raw_arms["correct"].reshape(-1),
                                     1.0 - confidence.reshape(-1))
    law_contract = {
        "shape": list(simplex.shape),
        "expected_columns": 8,
        "minimum_weight": float(simplex.min()),
        "maximum_row_sum_error": float(np.max(np.abs(simplex.sum(axis=1) - 1.0))),
        "population_column_mass": float(simplex[:, 0].mean()),
        "ordered_anchor_mass": float(simplex[:, 1:7].sum(axis=1).mean()),
        "uniform_abstention_mass": float(simplex[:, 7].mean()),
        "direct_scalar_output_used_for_deployment": False,
        "diagnostic_raw_activity_scored": True,
    }

    target_map_rows = []
    for dataset, rows in (("PKIS1", source_targets), ("PKIS2_COLD", cold_targets)):
        for item in rows:
            target_map_rows.append({"dataset": dataset, "assay_target": item["assay_target"],
                                    **asdict(item["record"])})
    with (output / "target_mapping.jsonl").open("w", encoding="utf-8") as handle:
        for row in target_map_rows:
            handle.write(json.dumps(row, sort_keys=True) + "\n")
    _write_json(output / "target_exclusions.json", {
        "pkis1": source_target_exclusions, "pkis2": transfer_target_exclusions,
    })
    _write_json(output / "ligand_exclusions.json", {
        "pkis1_parse": source_ligand_exclusions,
        "pkis2_parse": transfer_ligand_exclusions,
        "pkis2_not_scaffold_cold": int(len(transfer_ligands_all) - len(transfer_ligands)),
    })
    _write_json(output / "channel_bounds.json", {
        "channel_names": CHANNEL_NAMES, "shell_names": SHELL_NAMES,
        "bounds": channel_bounds,
        "coefficient_l2_by_channel": np.linalg.norm(
            pair_model.coef_.reshape(len(CHANNEL_NAMES), -1), axis=1
        ).tolist(),
        "bounded_z_min": biological_z.reshape(-1, len(CHANNEL_NAMES)).min(axis=0).tolist(),
        "bounded_z_max": biological_z.reshape(-1, len(CHANNEL_NAMES)).max(axis=0).tolist(),
    })
    _write_json(output / "derangement.json", {
        "mapping": [{"correct": records[index].hgnc,
                     "deranged": records[int(derangement[index])].hgnc,
                     "correct_group": records[index].group,
                     "deranged_group": records[int(derangement[index])].group}
                    for index in range(len(records))],
        "fixed_points": int(np.sum(derangement == np.arange(len(records)))),
        "same_group_fraction": float(np.mean([
            records[index].group == records[int(derangement[index])].group
            for index in range(len(records))
        ])),
    })

    result = {
        "schema": "MetaSieve.PKISTypedMechanismPilot.v1",
        "formal_gate": False,
        "stage_relation": "EXTERNAL_EXPLORATORY_NOT_E_AFF_X1",
        "main_verdict": main["status"],
        "end_to_end_dta_validated": bool(main["interaction_signal_observed"]
                                         and main["raw_location_signal_observed"]),
        "admission_to_biological_z_authorized": False,
        "sample": {
            "pkis1_ligands": len(source_ligands), "pkis1_targets": len(source_targets),
            "pkis2_scaffold_cold_ligands": len(transfer_ligands),
            "pkis2_exact_target_cold_targets": len(cold_targets),
            "pkis2_family_cold_targets": int(len(family_cold)),
            "pkis2_group_cold_targets": int(len(group_cold)),
            "exact_smiles_overlap": int(len(source_smiles & {r.canonical_smiles for r in transfer_ligands})),
            "generic_scaffold_overlap": int(len(source_scaffolds & {r.generic_scaffold for r in transfer_ligands})),
        },
        "source_label_summary": _spectral_summary(Y1),
        "transfer_label_summary": _spectral_summary(Y2),
        "source_only_model_selection": {
            "alphas": ALPHAS,
            "ligand": {"selected_alpha": ligand_alpha, **ligand_cv},
            "protein": {"selected_alpha": protein_alpha, **protein_cv},
            "typed_pair": {"selected_alpha": pair_alpha, **pair_cv},
        },
        "main_stratum": main,
        "family_cold_stratum": family,
        "group_cold_stratum": group,
        "per_klifs_group": group_effects,
        "coverage": {
            "ligand_nearest_tanimoto": {
                "min": float(ligand_similarity.min()), "median": float(np.median(ligand_similarity)),
                "max": float(ligand_similarity.max()),
            },
            "protein_nearest_pocket_identity": {
                "min": float(protein_identity.min()), "median": float(np.median(protein_identity)),
                "max": float(protein_identity.max()),
            },
            "pair_confidence": {
                "min": float(confidence.min()), "median": float(np.median(confidence)),
                "max": float(confidence.max()),
            },
        },
        "law_interface_contract": law_contract,
        "read_firewall": {
            "chembl_x1_labels": 0, "davis_labels": 0, "recipient_labels": 0,
            "pkis1_role": "source_development", "pkis2_role": "external_transfer",
        },
        "limitations": [
            "PKIS reports single-concentration inhibition, not equilibrium Ki/Kd.",
            "Both transfer panels are kinase-focused; non-kinase generalization remains untested.",
            "The five channels are mechanism proxies, not binding free-energy components.",
            "This exploratory result cannot change X1 authorization or admit biological z.",
        ],
    }
    _write_json(output / "result.json", result)

    manifest = {
        "schema": "MetaSieve.PKISTypedMechanismPilotManifest.v1",
        "created_unix": time.time(), "runtime_seconds": time.time() - started,
        "seed": SEED, "bootstrap_draws": BOOTSTRAP_DRAWS,
        "git_head": _git_value(repo_root, "rev-parse", "HEAD"),
        "git_branch": _git_value(repo_root, "branch", "--show-current"),
        "input_sha256": {name: sha256_file(path) for name, path in paths.items()},
        "code_sha256": {
            "run_pilot.py": sha256_file(__file__),
            "mechanism.py": sha256_file(Path(__file__).with_name("mechanism.py")),
        },
        "outputs": {},
    }
    for path in sorted(output.iterdir()):
        if path.name != "manifest.json" and path.is_file():
            manifest["outputs"][path.name] = sha256_file(path)
    _write_json(output / "manifest.json", manifest)
    print(json.dumps({
        "verdict": result["main_verdict"], "sample": result["sample"],
        "interaction": main["interaction_signal_observed"],
        "location": main["raw_location_signal_observed"],
        "runtime_seconds": manifest["runtime_seconds"],
    }, indent=2), flush=True)
    return result


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--informers-root", required=True)
    parser.add_argument("--klifs-json", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args(argv)


if __name__ == "__main__":
    warnings.filterwarnings("ignore", category=RuntimeWarning)
    run(parse_args())

