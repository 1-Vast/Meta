"""F0R: leakage-safe active support design plus a learned reliability gate."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.special import expit, logit
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from ..pkis_mechanism_pilot.mechanism import (
    protein_model_features,
    stable_fold,
)
from .ceiling_probe import (
    BOOTSTRAP_DRAWS,
    ENTITY_RIDGE,
    SEEDS,
    _fit_basis,
    _load,
    _predict_components,
    _sha256,
    _write_json,
)


RANK = 2  # frozen from F0 source-only selection
LOCATION_RIDGES = (0.1, 1.0, 10.0, 100.0)
INTERACTION_RIDGES = (0.1, 1.0, 10.0, 100.0)
GATE_RIDGE = 10.0


def _design(scaffold, factor, finite, k, policy, seed):
    finite = np.asarray(finite, dtype=np.int64)
    by_scaffold = {}
    for index in finite:
        by_scaffold.setdefault(str(scaffold[index]), []).append(int(index))
    if len(by_scaffold) < k:
        return None
    rng = np.random.default_rng(seed)
    if policy == "random":
        chosen_scaffolds = rng.choice(sorted(by_scaffold), size=k, replace=False)
        return np.sort(np.asarray([
            rng.choice(by_scaffold[str(value)]) for value in chosen_scaffolds
        ], dtype=np.int64))
    if policy != "d_optimal":
        raise ValueError(policy)

    # One deterministic representative per scaffold avoids scaffold duplication.
    representatives = np.asarray([
        min(indices) for _, indices in sorted(by_scaffold.items())
    ], dtype=np.int64)
    design = np.concatenate([
        np.ones((len(representatives), 1)), factor[representatives]
    ], axis=1)
    gram = np.eye(design.shape[1], dtype=np.float64) * 1e-3
    selected = []
    remaining = list(range(len(representatives)))
    for _ in range(k):
        scores = []
        for local in remaining:
            candidate = gram + np.outer(design[local], design[local])
            sign, value = np.linalg.slogdet(candidate)
            scores.append(value if sign > 0 else -np.inf)
        best_value = max(scores)
        tied = [remaining[i] for i, value in enumerate(scores)
                if np.isclose(value, best_value, rtol=1e-12, atol=1e-14)]
        chosen = min(tied)
        selected.append(chosen)
        gram += np.outer(design[chosen], design[chosen])
        remaining.remove(chosen)
    return np.sort(representatives[np.asarray(selected, dtype=np.int64)])


def _query_indices(y, target, support, scaffold):
    forbidden = {str(scaffold[index]) for index in support}
    return np.asarray([
        index for index in range(len(scaffold))
        if np.isfinite(y[index, target]) and index not in set(support)
        and str(scaffold[index]) not in forbidden
    ], dtype=np.int64)


def _solve(v, residual, prior, location_ridge, interaction_ridge):
    design = np.concatenate([np.ones((len(v), 1)), v], axis=1)
    prior_full = np.concatenate([[0.0], prior])
    diagonal = [float(location_ridge)] + [float(interaction_ridge)] * v.shape[1]
    penalty = np.diag(diagonal)
    return np.linalg.solve(
        design.T @ design + penalty,
        design.T @ residual + penalty @ prior_full,
    )


def _condition(v):
    design = np.concatenate([np.ones((len(v), 1)), v], axis=1)
    singular = np.linalg.svd(design, compute_uv=False)
    minimum = float(singular[-1]) if len(singular) else 0.0
    maximum = float(singular[0]) if len(singular) else 0.0
    ratio = minimum / maximum if maximum > 0 else 0.0
    return minimum, maximum, ratio


def _gate_features(v, residual, prior, coefficient):
    design = np.concatenate([np.ones((len(v), 1)), v], axis=1)
    prior_full = np.concatenate([[0.0], prior])
    pred_prior = design @ prior_full
    pred_adapt = design @ coefficient
    prior_rmse = float(np.sqrt(np.mean((residual - pred_prior) ** 2)))
    adapt_rmse = float(np.sqrt(np.mean((residual - pred_adapt) ** 2)))
    minimum, maximum, ratio = _condition(v)
    return np.asarray([
        np.log1p(len(v)),
        minimum,
        maximum,
        ratio,
        float(np.std(residual)),
        prior_rmse,
        adapt_rmse,
        (prior_rmse - adapt_rmse) / max(prior_rmse, 1e-8),
        float(np.linalg.norm(coefficient - prior_full)),
        float(np.linalg.norm(prior)),
        float(np.mean(np.linalg.norm(v, axis=1))),
    ], dtype=np.float64)


def _optimal_gate(outcome, prior_prediction, adapted_prediction):
    direction = adapted_prediction - prior_prediction
    denominator = float(direction @ direction)
    if denominator <= 1e-12:
        return 0.0
    return float(np.clip(direction @ (outcome - prior_prediction) / denominator, 0.0, 1.0))


def _safe_spearman(y, prediction):
    if len(y) < 3 or np.std(y) < 1e-12 or np.std(prediction) < 1e-12:
        return 0.0
    return float(spearmanr(y, prediction).statistic)


def _fold_cache(y, x_ligand, x_protein, ligand_groups, protein_groups):
    caches = []
    ligand_folds = np.asarray([stable_fold(value, 3) for value in ligand_groups])
    protein_folds = np.asarray([stable_fold(value, 3) for value in protein_groups])
    for fold in range(3):
        lt, lv = np.flatnonzero(ligand_folds != fold), np.flatnonzero(ligand_folds == fold)
        pt, pv = np.flatnonzero(protein_folds != fold), np.flatnonzero(protein_folds == fold)
        if min(len(lt), len(lv), len(pt), len(pv)) == 0:
            continue
        model = _fit_basis(y[np.ix_(lt, pt)], x_ligand[lt], x_protein[pt], RANK)
        additive, ligand_factor, protein_prior = _predict_components(
            model, x_ligand[lv], x_protein[pv])
        caches.append({
            "fold": fold,
            "y": y[np.ix_(lv, pv)],
            "additive": additive,
            "ligand_factor": ligand_factor,
            "protein_prior": protein_prior,
            "scaffold": np.asarray(ligand_groups, dtype=object)[lv],
        })
    return caches


def _cv_penalties(caches, k=5):
    scores = {}
    for location_ridge in LOCATION_RIDGES:
        for interaction_ridge in INTERACTION_RIDGES:
            values = []
            for cache in caches:
                y = cache["y"]
                for target in range(y.shape[1]):
                    finite = np.flatnonzero(np.isfinite(y[:, target]))
                    support = _design(cache["scaffold"], cache["ligand_factor"], finite,
                                      k, "d_optimal", 20260808 + target)
                    if support is None:
                        continue
                    query = _query_indices(y, target, support, cache["scaffold"])
                    if len(query) < 3:
                        continue
                    residual = y[support, target] - cache["additive"][support, target]
                    coefficient = _solve(
                        cache["ligand_factor"][support], residual,
                        cache["protein_prior"][target], location_ridge, interaction_ridge)
                    prediction = (cache["additive"][query, target] + coefficient[0]
                                  + cache["ligand_factor"][query] @ coefficient[1:])
                    values.append(float(np.mean((y[query, target]
                                                 - np.clip(prediction, 0, 1)) ** 2)))
            scores[(location_ridge, interaction_ridge)] = values
    selected = min(
        (float(np.mean(values)), pair) for pair, values in scores.items() if values
    )[1]
    return selected, {
        f"location={pair[0]},interaction={pair[1]}": {
            "mean_mse": float(np.mean(values)), "n_episodes": len(values)
        }
        for pair, values in scores.items() if values
    }


def _fit_reliability(caches, location_ridge, interaction_ridge, k=5):
    features, targets = [], []
    for cache in caches:
        y = cache["y"]
        for target in range(y.shape[1]):
            finite = np.flatnonzero(np.isfinite(y[:, target]))
            support = _design(cache["scaffold"], cache["ligand_factor"], finite,
                              k, "d_optimal", 20260808 + target)
            if support is None:
                continue
            query = _query_indices(y, target, support, cache["scaffold"])
            if len(query) < 3:
                continue
            v_support = cache["ligand_factor"][support]
            residual = y[support, target] - cache["additive"][support, target]
            prior = cache["protein_prior"][target]
            coefficient = _solve(v_support, residual, prior,
                                 location_ridge, interaction_ridge)
            p0 = cache["additive"][query, target] + cache["ligand_factor"][query] @ prior
            p1 = (cache["additive"][query, target] + coefficient[0]
                  + cache["ligand_factor"][query] @ coefficient[1:])
            features.append(_gate_features(v_support, residual, prior, coefficient))
            targets.append(_optimal_gate(y[query, target], p0, p1))
    features = np.asarray(features)
    targets = np.asarray(targets)
    transformed = logit(np.clip(targets, 0.02, 0.98))
    model = make_pipeline(StandardScaler(), Ridge(alpha=GATE_RIDGE)).fit(features, transformed)
    fitted = expit(model.predict(features))
    return model, {
        "n_episodes": len(targets),
        "target_mean": float(np.mean(targets)),
        "target_zero_fraction": float(np.mean(targets <= 1e-12)),
        "target_one_fraction": float(np.mean(targets >= 1.0 - 1e-12)),
        "training_mae": float(np.mean(np.abs(fitted - targets))),
        "fitted_mean": float(np.mean(fitted)),
    }


def _evaluate(y, additive, factor, prior, scaffold, k, policy,
              location_ridge, interaction_ridge, gate, seeds):
    names = (
        "support_free", "location_only", "adapted", "gated", "gated_no_protein",
        "gated_deranged_protein", "gated_wrong_support", "gated_permuted_support",
    )
    n_target = y.shape[1]
    per_target = {name: {metric: np.full((n_target, len(seeds)), np.nan)
                         for metric in ("mse", "mae", "spearman")}
                  for name in names}
    condition = np.full((n_target, len(seeds)), np.nan)
    gate_values = {name: np.full((n_target, len(seeds)), np.nan)
                   for name in names if name.startswith("gated")}
    deranged = np.roll(np.arange(n_target), 1)

    for target in range(n_target):
        wrong = int(deranged[target])
        finite = np.flatnonzero(np.isfinite(y[:, target]) & np.isfinite(y[:, wrong]))
        for seed_index, seed in enumerate(seeds):
            support = _design(scaffold, factor, finite, k, policy,
                              seed + 104729 * target)
            if support is None:
                continue
            query = _query_indices(y, target, support, scaffold)
            if len(query) < 3:
                continue
            v_support, v_query = factor[support], factor[query]
            residual = y[support, target] - additive[support, target]
            wrong_residual = y[support, wrong] - additive[support, wrong]
            correct_prior = prior[target]
            wrong_prior = prior[wrong]
            condition[target, seed_index] = _condition(v_support)[2]

            def adapt(this_residual, this_prior):
                coefficient = _solve(v_support, this_residual, this_prior,
                                     location_ridge, interaction_ridge)
                p0 = additive[query, target] + v_query @ this_prior
                p1 = additive[query, target] + coefficient[0] + v_query @ coefficient[1:]
                features = _gate_features(v_support, this_residual, this_prior, coefficient)
                weight = float(expit(gate.predict(features[None, :]))[0])
                return p0, p1, p0 + weight * (p1 - p0), weight

            p0, p1, pg, g = adapt(residual, correct_prior)
            _, _, png, ng = adapt(residual, np.zeros_like(correct_prior))
            _, _, pdg, dg = adapt(residual, wrong_prior)
            _, _, pwg, wg = adapt(wrong_residual, correct_prior)
            rng = np.random.default_rng(seed + 65537 * target)
            _, _, ppg, pgw = adapt(residual[rng.permutation(k)], correct_prior)
            location = float(np.sum(residual) / (k + location_ridge))
            predictions = {
                "support_free": p0,
                "location_only": additive[query, target] + location,
                "adapted": p1,
                "gated": pg,
                "gated_no_protein": png,
                "gated_deranged_protein": pdg,
                "gated_wrong_support": pwg,
                "gated_permuted_support": ppg,
            }
            gate_values["gated"][target, seed_index] = g
            gate_values["gated_no_protein"][target, seed_index] = ng
            gate_values["gated_deranged_protein"][target, seed_index] = dg
            gate_values["gated_wrong_support"][target, seed_index] = wg
            gate_values["gated_permuted_support"][target, seed_index] = pgw
            outcome = y[query, target]
            for name, prediction in predictions.items():
                prediction = np.clip(prediction, 0.0, 1.0)
                error = outcome - prediction
                per_target[name]["mse"][target, seed_index] = float(np.mean(error ** 2))
                per_target[name]["mae"][target, seed_index] = float(np.mean(np.abs(error)))
                per_target[name]["spearman"][target, seed_index] = _safe_spearman(
                    outcome, prediction)
    return per_target, condition, gate_values


def _bootstrap(comparator, candidate, seed=20260808):
    difference = np.nanmean(comparator, axis=1) - np.nanmean(candidate, axis=1)
    difference = difference[np.isfinite(difference)]
    rng = np.random.default_rng(seed)
    draws = difference[rng.integers(0, len(difference),
                                   size=(BOOTSTRAP_DRAWS, len(difference)))].mean(axis=1)
    return {
        "estimate": float(np.mean(difference)),
        "ci95": [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))],
        "n_targets": len(difference), "draws": BOOTSTRAP_DRAWS, "seed": seed,
    }


def _summary(per_target, condition, gates):
    metrics = {}
    for name, values in per_target.items():
        metrics[name] = {metric: float(np.nanmean(array))
                         for metric, array in values.items()}
    contrasts = {}
    for comparator in per_target:
        if comparator == "gated":
            continue
        contrasts[comparator] = _bootstrap(
            per_target[comparator]["mse"], per_target["gated"]["mse"])
    finite_condition = condition[np.isfinite(condition)]
    return {
        "metrics": metrics,
        "mse_reduction_gated_vs": contrasts,
        "condition_ratio": {
            "median": float(np.median(finite_condition)),
            "q10": float(np.quantile(finite_condition, 0.1)),
            "q90": float(np.quantile(finite_condition, 0.9)),
        },
        "gate": {name: float(np.nanmean(value)) for name, value in gates.items()},
    }


def run(args):
    started = time.time()
    data = _load(args)
    y_source = data["y_source"]
    source_ligands = data["source_ligands"]
    source_targets = data["source_targets"]
    x_ligand = np.stack([record.nuisance_features for record in source_ligands])
    x_protein = np.stack([protein_model_features(item["record"].pocket)
                          for item in source_targets])
    ligand_groups = [record.generic_scaffold for record in source_ligands]
    protein_groups = [item["record"].group for item in source_targets]

    print("building source-only dual-cold cache", flush=True)
    caches = _fold_cache(y_source, x_ligand, x_protein, ligand_groups, protein_groups)
    (location_ridge, interaction_ridge), cv = _cv_penalties(caches, k=5)
    print(f"selected location={location_ridge} interaction={interaction_ridge}", flush=True)
    gate, gate_report = _fit_reliability(
        caches, location_ridge, interaction_ridge, k=5)
    full = _fit_basis(y_source, x_ligand, x_protein, RANK)

    transfers = {}
    for panel in data["panels"]:
        print(f"evaluating {panel['name']}", flush=True)
        panel_x_ligand = np.stack([record.nuisance_features for record in panel["ligands"]])
        panel_x_protein = np.stack([protein_model_features(item["record"].pocket)
                                    for item in panel["targets"]])
        additive, factor, prior = _predict_components(full, panel_x_ligand, panel_x_protein)
        panel_result = {"support_sizes": {}, "n_ligands": len(panel["ligands"]),
                        "n_targets": len(panel["targets"]), "role": panel["role"]}
        scaffolds = np.asarray([record.generic_scaffold for record in panel["ligands"]],
                               dtype=object)
        for k in (5, 20):
            panel_result["support_sizes"][str(k)] = {}
            for policy in ("random", "d_optimal"):
                seeds = SEEDS if policy == "random" else (20260808,)
                episode = _evaluate(panel["y"], additive, factor, prior, scaffolds,
                                    k, policy, location_ridge, interaction_ridge,
                                    gate, seeds)
                panel_result["support_sizes"][str(k)][policy] = _summary(*episode)
        transfers[panel["name"]] = panel_result

    primary = transfers["PKIS2"]["support_sizes"]["5"]["d_optimal"]
    external = transfers["Anastassiadis2011"]["support_sizes"]["5"]["d_optimal"]
    required = ("support_free", "location_only", "gated_wrong_support",
                "gated_permuted_support")
    pkis_pass = all(primary["mse_reduction_gated_vs"][name]["ci95"][0] > 0
                    for name in required)
    external_pass = all(external["mse_reduction_gated_vs"][name]["estimate"] > 0
                        for name in required)
    protein_point = primary["mse_reduction_gated_vs"][
        "gated_deranged_protein"]["estimate"] > 0
    design_improvement = (
        transfers["PKIS2"]["support_sizes"]["5"]["random"]["metrics"]["gated"]["mse"]
        - primary["metrics"]["gated"]["mse"] > 0
    )
    condition_improvement = (
        primary["condition_ratio"]["median"]
        > transfers["PKIS2"]["support_sizes"]["5"]["random"]["condition_ratio"]["median"]
    )
    passed = bool(pkis_pass and external_pass and protein_point
                  and design_improvement and condition_improvement)
    result = {
        "schema": "MetaSieve.SectionReliability.F0R.v1",
        "selected": {"rank": RANK, "location_ridge": location_ridge,
                     "interaction_ridge": interaction_ridge,
                     "gate_ridge": GATE_RIDGE, "d_adapt": RANK + 1},
        "source_cv": cv,
        "reliability_gate": gate_report,
        "transfers": transfers,
        "gate": {
            "pkis2": pkis_pass,
            "anastassiadis_point_estimate": external_pass,
            "protein_point_estimate": protein_point,
            "d_optimal_beats_random": design_improvement,
            "condition_improved": condition_improvement,
            "passed": passed,
            "verdict": "F0R_CEILING_VIABLE" if passed else "F0R_CEILING_NOT_VIABLE",
        },
        "read_firewall": {"kcgs_numeric_outcomes": "NOT_READ",
                          "davis_labels": "NOT_READ", "recipient_labels": "NOT_READ"},
        "elapsed_seconds": time.time() - started,
    }
    output = Path(args.output).resolve()
    _write_json(output / "result.json", result)
    _write_json(output / "manifest.json", {
        "sha256": {
            "preregistration": _sha256(Path(__file__).with_name("F0R_PREREGISTRATION.md")),
            "script": _sha256(Path(__file__)),
        },
        "parameters": vars(args),
    })
    print(json.dumps(result["gate"], indent=2), flush=True)


def parser():
    item = argparse.ArgumentParser()
    item.add_argument("--informers-root", default="../external/informers")
    item.add_argument("--klifs-json", default="../external/klifs/kinase_information_human.json")
    item.add_argument("--anastassiadis-workbook",
                      default="external/anastassiadis/NIHMS328213-supplement-3.xls")
    item.add_argument("--anastassiadis-identities",
                      default="external/anastassiadis/compound_identities.json")
    item.add_argument("--output", default="research/section_operator_pilot/artifacts/f0r")
    return item


if __name__ == "__main__":
    run(parser().parse_args())
