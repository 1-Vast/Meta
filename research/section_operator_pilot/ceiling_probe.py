"""Source-only low-rank ceiling for a few-shot affinity section.

This is Stage F0 from PREREGISTRATION.md.  It deliberately uses a closed-form
positive-ridge section solve as an upper-bound/control.  It never imports or
modifies ``model/`` and it does not read DAVIS or recipient labels.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from ..pkis_mechanism_pilot.mechanism import (
    double_center,
    load_klifs,
    load_smiles,
    protein_model_features,
    stable_fold,
)
from ..pkis_mechanism_pilot.run_kernel_revision import (
    _ligand_records,
    _panel_from_anastassiadis,
    _panel_from_pkis2,
)
from ..pkis_mechanism_pilot.run_pilot import _target_mapping


SEEDS = tuple(range(20260808, 20260828))
RANKS = (1, 2, 3, 4)
ADAPT_RIDGES = (0.01, 0.1, 1.0, 10.0, 100.0)
ENTITY_RIDGE = 1000.0
LOCATION_RIDGE = 0.1
BOOTSTRAP_DRAWS = 10_000


def _write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, default=str) + "\n",
                    encoding="utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _ridge(alpha: float):
    return make_pipeline(
        StandardScaler(),
        Ridge(alpha=float(alpha), solver="lsqr", tol=1e-8),
    )


@dataclass
class FittedBasis:
    rank: int
    mu: float
    ligand_main: object
    protein_main: object
    ligand_factor: object
    protein_factor: object
    singular_values: np.ndarray


def _fit_basis(y: np.ndarray, x_ligand: np.ndarray, x_protein: np.ndarray,
               rank: int) -> FittedBasis:
    mu = float(np.mean(y))
    ligand_effect = np.mean(y, axis=1) - mu
    protein_effect = np.mean(y, axis=0) - mu
    ligand_main = _ridge(ENTITY_RIDGE).fit(x_ligand, ligand_effect)
    protein_main = _ridge(ENTITY_RIDGE).fit(x_protein, protein_effect)

    residual = double_center(y)
    left, singular, right_t = np.linalg.svd(residual, full_matrices=False)
    scale = np.sqrt(np.clip(singular[:rank], 0.0, None))
    ligand_factor_values = left[:, :rank] * scale[None, :]
    protein_factor_values = right_t[:rank].T * scale[None, :]
    ligand_factor = _ridge(ENTITY_RIDGE).fit(x_ligand, ligand_factor_values)
    protein_factor = _ridge(ENTITY_RIDGE).fit(x_protein, protein_factor_values)
    return FittedBasis(
        rank=rank,
        mu=mu,
        ligand_main=ligand_main,
        protein_main=protein_main,
        ligand_factor=ligand_factor,
        protein_factor=protein_factor,
        singular_values=singular,
    )


def _predict_components(model: FittedBasis, x_ligand: np.ndarray,
                        x_protein: np.ndarray):
    ligand_main = np.asarray(model.ligand_main.predict(x_ligand), dtype=np.float64)
    protein_main = np.asarray(model.protein_main.predict(x_protein), dtype=np.float64)
    ligand_factor = np.asarray(model.ligand_factor.predict(x_ligand), dtype=np.float64)
    protein_factor = np.asarray(model.protein_factor.predict(x_protein), dtype=np.float64)
    if ligand_factor.ndim == 1:
        ligand_factor = ligand_factor[:, None]
    if protein_factor.ndim == 1:
        protein_factor = protein_factor[:, None]
    additive = model.mu + ligand_main[:, None] + protein_main[None, :]
    return additive, ligand_factor, protein_factor


def _solve_section(v: np.ndarray, residual: np.ndarray, prior: np.ndarray,
                   interaction_ridge: float) -> np.ndarray:
    """Positive-ridge solve for `(location, interaction coordinates)`."""
    v = np.asarray(v, dtype=np.float64)
    residual = np.asarray(residual, dtype=np.float64)
    prior = np.asarray(prior, dtype=np.float64)
    design = np.concatenate([np.ones((len(v), 1)), v], axis=1)
    prior_full = np.concatenate([[0.0], prior])
    penalty = np.diag([LOCATION_RIDGE] + [float(interaction_ridge)] * v.shape[1])
    lhs = design.T @ design + penalty
    rhs = design.T @ residual + penalty @ prior_full
    return np.linalg.solve(lhs, rhs)


def _safe_spearman(y: np.ndarray, prediction: np.ndarray) -> float:
    mask = np.isfinite(y) & np.isfinite(prediction)
    if mask.sum() < 3 or np.std(y[mask]) < 1e-12 or np.std(prediction[mask]) < 1e-12:
        return 0.0
    return float(spearmanr(y[mask], prediction[mask]).statistic)


def _episode_predictions(y: np.ndarray, additive: np.ndarray,
                         ligand_factor: np.ndarray, protein_prior: np.ndarray,
                         k: int, interaction_ridge: float, seeds=SEEDS):
    n_ligand, n_target = y.shape
    deranged = np.roll(np.arange(n_target), 1)
    names = (
        "support_free", "additive", "location_only", "correct", "no_protein",
        "deranged_protein", "wrong_support", "permuted_support",
    )
    per_target = {name: {"mse": np.full((n_target, len(seeds)), np.nan),
                         "mae": np.full((n_target, len(seeds)), np.nan),
                         "spearman": np.full((n_target, len(seeds)), np.nan)}
                  for name in names}

    for target in range(n_target):
        wrong = int(deranged[target])
        jointly_finite = np.flatnonzero(np.isfinite(y[:, target]) & np.isfinite(y[:, wrong]))
        if len(jointly_finite) <= k + 2:
            continue
        for seed_index, seed in enumerate(seeds):
            rng = np.random.default_rng(seed + 104729 * target)
            support = np.sort(rng.choice(jointly_finite, size=k, replace=False))
            query_mask = np.isfinite(y[:, target])
            query_mask[support] = False
            query = np.flatnonzero(query_mask)
            if len(query) < 3:
                continue

            v_support = ligand_factor[support]
            v_query = ligand_factor[query]
            support_residual = y[support, target] - additive[support, target]
            wrong_residual = y[support, wrong] - additive[support, wrong]
            prior = protein_prior[target]
            wrong_prior = protein_prior[wrong]

            correct_coef = _solve_section(
                v_support, support_residual, prior, interaction_ridge)
            no_protein_coef = _solve_section(
                v_support, support_residual, np.zeros_like(prior), interaction_ridge)
            deranged_coef = _solve_section(
                v_support, support_residual, wrong_prior, interaction_ridge)
            wrong_support_coef = _solve_section(
                v_support, wrong_residual, prior, interaction_ridge)
            permutation = rng.permutation(k)
            permuted_coef = _solve_section(
                v_support, support_residual[permutation], prior, interaction_ridge)
            location = float(np.sum(support_residual) / (k + LOCATION_RIDGE))

            base_query = additive[query, target]
            predictions = {
                "support_free": base_query + v_query @ prior,
                "additive": base_query,
                "location_only": base_query + location,
                "correct": base_query + correct_coef[0] + v_query @ correct_coef[1:],
                "no_protein": base_query + no_protein_coef[0] + v_query @ no_protein_coef[1:],
                "deranged_protein": base_query + deranged_coef[0] + v_query @ deranged_coef[1:],
                "wrong_support": base_query + wrong_support_coef[0] + v_query @ wrong_support_coef[1:],
                "permuted_support": base_query + permuted_coef[0] + v_query @ permuted_coef[1:],
            }
            outcome = y[query, target]
            for name, prediction in predictions.items():
                prediction = np.clip(prediction, 0.0, 1.0)
                error = outcome - prediction
                per_target[name]["mse"][target, seed_index] = float(np.mean(error ** 2))
                per_target[name]["mae"][target, seed_index] = float(np.mean(np.abs(error)))
                per_target[name]["spearman"][target, seed_index] = _safe_spearman(
                    outcome, prediction)
    return per_target


def _paired_bootstrap(comparator: np.ndarray, candidate: np.ndarray,
                      seed: int = 20260808):
    difference = np.nanmean(comparator, axis=1) - np.nanmean(candidate, axis=1)
    difference = difference[np.isfinite(difference)]
    if not len(difference):
        return {"estimate": None, "ci95": [None, None], "n_targets": 0}
    rng = np.random.default_rng(seed)
    draws = difference[rng.integers(0, len(difference),
                                    size=(BOOTSTRAP_DRAWS, len(difference)))].mean(axis=1)
    return {
        "estimate": float(np.mean(difference)),
        "ci95": [float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))],
        "n_targets": int(len(difference)),
        "draws": BOOTSTRAP_DRAWS,
        "seed": seed,
    }


def _summarize(per_target):
    summary = {}
    for name, metrics in per_target.items():
        summary[name] = {
            key: float(np.nanmean(value)) for key, value in metrics.items()
        }
        summary[name]["finite_targets"] = int(
            np.isfinite(np.nanmean(metrics["mse"], axis=1)).sum())
    contrasts = {}
    for comparator in (
        "support_free", "additive", "location_only", "no_protein",
        "deranged_protein", "wrong_support", "permuted_support",
    ):
        contrasts[comparator] = _paired_bootstrap(
            per_target[comparator]["mse"], per_target["correct"]["mse"])
    return summary, contrasts


def _source_cv(y, x_ligand, x_protein, ligand_groups, protein_groups, k_values):
    report = {}
    for rank in RANKS:
        report[str(rank)] = {}
        for ridge in ADAPT_RIDGES:
            fold_mse = {k: [] for k in k_values}
            for fold in range(3):
                ligand_train = np.flatnonzero(
                    np.asarray([stable_fold(value, 3) for value in ligand_groups]) != fold)
                ligand_valid = np.flatnonzero(
                    np.asarray([stable_fold(value, 3) for value in ligand_groups]) == fold)
                protein_train = np.flatnonzero(
                    np.asarray([stable_fold(value, 3) for value in protein_groups]) != fold)
                protein_valid = np.flatnonzero(
                    np.asarray([stable_fold(value, 3) for value in protein_groups]) == fold)
                if min(len(ligand_train), len(ligand_valid),
                       len(protein_train), len(protein_valid)) == 0:
                    continue
                model = _fit_basis(
                    y[np.ix_(ligand_train, protein_train)],
                    x_ligand[ligand_train], x_protein[protein_train], rank)
                additive, ligand_factor, protein_prior = _predict_components(
                    model, x_ligand[ligand_valid], x_protein[protein_valid])
                y_valid = y[np.ix_(ligand_valid, protein_valid)]
                for k in k_values:
                    if len(ligand_valid) <= k + 2:
                        continue
                    episodes = _episode_predictions(
                        y_valid, additive, ligand_factor, protein_prior, k, ridge,
                        seeds=SEEDS[:5],
                    )
                    fold_mse[k].append(float(np.nanmean(episodes["correct"]["mse"])))
            report[str(rank)][str(ridge)] = {
                str(k): values for k, values in fold_mse.items()
            }
    candidates = []
    primary_k = min(k_values)
    for rank in RANKS:
        for ridge in ADAPT_RIDGES:
            values = report[str(rank)][str(ridge)][str(primary_k)]
            if values:
                candidates.append((float(np.mean(values)), rank, ridge))
    if not candidates:
        raise RuntimeError("source dual-cold CV produced no selectable candidate")
    _, selected_rank, selected_ridge = min(candidates)
    return selected_rank, selected_ridge, report


def _load(args):
    informers = Path(args.informers_root).resolve()
    klifs_path = Path(args.klifs_json).resolve()
    pkis1 = pd.read_csv(informers / "data/pkis1_continuous_labels.csv",
                        dtype={"molid": str}).set_index("molid")
    pkis2 = pd.read_csv(informers / "data/pkis2_continuous_labels.csv",
                        dtype={"molid": str}).set_index("molid")
    klifs_index, _ = load_klifs(klifs_path)
    source_targets, target_exclusions = _target_mapping(pkis1.columns, klifs_index)
    source_smiles_map = load_smiles(
        informers / "data/compounds/pkis1_uniq_tested.can")
    source_ligands, ligand_exclusions = _ligand_records([
        (str(molid), source_smiles_map[str(molid)])
        for molid in pkis1.index.astype(str)
    ])
    source_ids = [record.molid for record in source_ligands]
    source_columns = [item["assay_target"] for item in source_targets]
    y_source = pkis1.loc[source_ids, source_columns].to_numpy(dtype=np.float64)
    source_genes = {item["record"].hgnc for item in source_targets}
    source_smiles = {record.canonical_smiles for record in source_ligands}
    source_scaffolds = {record.generic_scaffold for record in source_ligands}
    panels = [
        _panel_from_pkis2(
            pkis2, informers / "data/compounds/pkis2_uniq_tested.can",
            klifs_index, source_genes, source_smiles, source_scaffolds),
        _panel_from_anastassiadis(
            Path(args.anastassiadis_workbook).resolve(),
            Path(args.anastassiadis_identities).resolve(), klifs_index,
            source_genes, source_smiles, source_scaffolds),
    ]
    return {
        "source_ligands": source_ligands,
        "source_targets": source_targets,
        "y_source": y_source,
        "panels": panels,
        "target_exclusions": target_exclusions,
        "ligand_exclusions": ligand_exclusions,
    }


def run(args):
    started = time.time()
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    data = _load(args)
    source_ligands = data["source_ligands"]
    source_targets = data["source_targets"]
    y_source = data["y_source"]
    x_ligand = np.stack([record.nuisance_features for record in source_ligands])
    x_protein = np.stack([
        protein_model_features(item["record"].pocket) for item in source_targets
    ])
    k_values = tuple(int(value) for value in args.k_values.split(","))

    print("selecting rank and section ridge on source-only dual-cold folds", flush=True)
    selected_rank, selected_ridge, cv = _source_cv(
        y_source, x_ligand, x_protein,
        [record.generic_scaffold for record in source_ligands],
        [item["record"].group for item in source_targets], k_values,
    )
    print(f"selected rank={selected_rank} ridge={selected_ridge}", flush=True)
    model = _fit_basis(y_source, x_ligand, x_protein, selected_rank)

    transfers = {}
    for panel in data["panels"]:
        print(f"evaluating {panel['name']}", flush=True)
        panel_x_ligand = np.stack([
            record.nuisance_features for record in panel["ligands"]
        ])
        panel_x_protein = np.stack([
            protein_model_features(item["record"].pocket) for item in panel["targets"]
        ])
        additive, ligand_factor, protein_prior = _predict_components(
            model, panel_x_ligand, panel_x_protein)
        panel_report = {
            "role": panel["role"],
            "n_ligands": len(panel["ligands"]),
            "n_targets": len(panel["targets"]),
            "n_finite_cells": int(np.isfinite(panel["y"]).sum()),
            "support_sizes": {},
        }
        for k in k_values:
            episodes = _episode_predictions(
                panel["y"], additive, ligand_factor, protein_prior,
                k, selected_ridge,
            )
            metrics, contrasts = _summarize(episodes)
            panel_report["support_sizes"][str(k)] = {
                "metrics": metrics,
                "mse_reduction_correct_vs": contrasts,
            }
        transfers[panel["name"]] = panel_report

    pkis_primary = transfers["PKIS2"]["support_sizes"][str(min(k_values))]
    ana_primary = transfers["Anastassiadis2011"]["support_sizes"][str(min(k_values))]
    required_pkis = ("support_free", "location_only", "wrong_support", "permuted_support")
    pkis_pass = all(
        pkis_primary["mse_reduction_correct_vs"][name]["ci95"][0] > 0.0
        for name in required_pkis
    )
    ana_pass = all(
        ana_primary["mse_reduction_correct_vs"][name]["estimate"] > 0.0
        for name in required_pkis
    )
    result = {
        "schema": "MetaSieve.SectionCeiling.F0.v1",
        "stage": "F0_CLOSED_FORM_CEILING_ONLY",
        "selected": {
            "rank": selected_rank,
            "adapt_ridge": selected_ridge,
            "entity_ridge": ENTITY_RIDGE,
            "location_ridge": LOCATION_RIDGE,
            "d_adapt": selected_rank + 1,
            "k_primary": min(k_values),
        },
        "source": {
            "panel": "PKIS1",
            "n_ligands": len(source_ligands),
            "n_targets": len(source_targets),
            "interaction_variance_fraction": float(
                np.var(double_center(y_source)) / np.var(y_source)),
            "singular_values": model.singular_values[:20].tolist(),
            "cv": cv,
        },
        "transfers": transfers,
        "f0_gate": {
            "pkis2": pkis_pass,
            "anastassiadis_point_estimate": ana_pass,
            "passed": bool(pkis_pass and ana_pass),
            "verdict": "F0_CEILING_VIABLE" if pkis_pass and ana_pass
                       else "F0_CEILING_NOT_VIABLE",
        },
        "read_firewall": {
            "pkis1": "source",
            "pkis2": "consumed_development",
            "anastassiadis": "consumed_development",
            "kcgs_numeric_outcomes": "NOT_READ",
            "davis_labels": "NOT_READ",
            "recipient_labels": "NOT_READ",
        },
        "limitations": [
            "Closed-form ridge is a ceiling/control, not the required learnable adapter.",
            "PKIS2 and Anastassiadis were already consumed development panels.",
            "This script does not establish law-valued operator validity or z admission.",
        ],
        "elapsed_seconds": time.time() - started,
    }
    _write_json(output / "result.json", result)
    manifest_paths = {
        "preregistration": Path(__file__).with_name("PREREGISTRATION.md"),
        "script": Path(__file__),
        "pkis1_labels": Path(args.informers_root).resolve() / "data/pkis1_continuous_labels.csv",
        "pkis2_labels": Path(args.informers_root).resolve() / "data/pkis2_continuous_labels.csv",
        "anastassiadis_workbook": Path(args.anastassiadis_workbook).resolve(),
    }
    _write_json(output / "manifest.json", {
        "sha256": {name: _sha256(path) for name, path in manifest_paths.items()},
        "parameters": vars(args),
    })
    print(json.dumps(result["f0_gate"], indent=2), flush=True)


def parser():
    item = argparse.ArgumentParser()
    item.add_argument("--informers-root", default="../external/informers")
    item.add_argument("--klifs-json", default="../external/klifs/kinase_information_human.json")
    item.add_argument("--anastassiadis-workbook",
                      default="external/anastassiadis/NIHMS328213-supplement-3.xls")
    item.add_argument("--anastassiadis-identities",
                      default="external/anastassiadis/compound_identities.json")
    item.add_argument("--k-values", default="5,20")
    item.add_argument("--output", default="research/section_operator_pilot/artifacts/f0")
    return item


if __name__ == "__main__":
    run(parser().parse_args())
