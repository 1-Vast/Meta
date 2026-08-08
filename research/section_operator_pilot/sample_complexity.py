"""F3K source-selected support budget for an identifiable partial section."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
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
from ..pkis_mechanism_pilot.run_kernel_revision import _ligand_records
from ..pkis_mechanism_pilot.run_pilot import _target_mapping
from .bioactivity_atlas import (
    _atlas_protein_views,
    _atlas_surfaces,
    _cv_penalties,
    _fit_ligand_atlas,
    _ligand_bits,
    _predict_profile,
    _protein_temperatures,
    _select_ligand_ridge,
    _uniform_surfaces,
)
from .ceiling_probe import BOOTSTRAP_DRAWS, SEEDS, _load
from .conditional_bilinear import (
    SEED,
    _d_optimal,
    _fit_additive,
    _nearest_nonself,
    _sha256,
    _subpocket_weights,
    _write_json,
)
from .identifiability_gate import (
    ARM_NAMES,
    CERTIFICATE_RIDGE,
    CONFORMAL_LEVEL,
    CONTROL_NAMES,
    _episode,
    _evaluate_selective,
    _selective_summary,
)


K_VALUES = (5, 10, 20, 40)
BUDGET_SEEDS = SEEDS[:5]


def _load_source(args):
    informers = Path(args.informers_root).resolve()
    pkis1 = pd.read_csv(
        informers / "data/pkis1_continuous_labels.csv",
        dtype={"molid": str},
    ).set_index("molid")
    klifs_index, _ = load_klifs(Path(args.klifs_json).resolve())
    targets, _ = _target_mapping(pkis1.columns, klifs_index)
    smiles = load_smiles(informers / "data/compounds/pkis1_uniq_tested.can")
    ligands, _ = _ligand_records([
        (str(molid), smiles[str(molid)]) for molid in pkis1.index.astype(str)
    ])
    ids = [record.molid for record in ligands]
    columns = [item["assay_target"] for item in targets]
    return {
        "ligands": ligands, "targets": targets,
        "y": pkis1.loc[ids, columns].to_numpy(dtype=np.float64),
    }


def _target_bootstrap_lcb(margins, admitted, targets, seed=SEED):
    target_values = []
    for target in sorted(set(targets.tolist())):
        mask = (targets == target) & admitted
        if mask.any():
            target_values.append(float(np.mean(margins[mask])))
    values = np.asarray(target_values, dtype=np.float64)
    if not len(values):
        return {"estimate": None, "ci95": [None, None], "n_targets": 0}
    rng = np.random.default_rng(seed)
    draws = values[rng.integers(
        0, len(values), size=(BOOTSTRAP_DRAWS, len(values)))].mean(axis=1)
    return {
        "estimate": float(np.mean(values)),
        "ci95": [float(np.quantile(draws, 0.025)),
                 float(np.quantile(draws, 0.975))],
        "n_targets": int(len(values)), "draws": BOOTSTRAP_DRAWS,
    }


def _fit_budget_certificate(features, margins, folds, targets):
    features = np.asarray(features, dtype=np.float64)
    margins = np.asarray(margins, dtype=np.float64)
    folds = np.asarray(folds, dtype=np.int64)
    targets = np.asarray(targets, dtype=np.int64)
    if not len(margins):
        return None, None, {
            "n_episodes": 0, "source_gate_passed": False,
            "failure": "no_valid_episodes",
        }
    crossfit = np.full(len(margins), np.nan)
    for fold in range(3):
        train, valid = folds != fold, folds == fold
        model = make_pipeline(
            StandardScaler(), Ridge(alpha=CERTIFICATE_RIDGE)
        ).fit(features[train], margins[train])
        crossfit[valid] = model.predict(features[valid])
    correlation = float(np.corrcoef(crossfit, margins)[0, 1])
    error = crossfit - margins
    quantile = float(np.quantile(error, CONFORMAL_LEVEL, method="higher"))
    certificate = crossfit - quantile
    admitted = certificate > 0.0
    interval = _target_bootstrap_lcb(margins, admitted, targets)
    rate = float(np.mean(admitted))
    passed = bool(
        np.isfinite(correlation) and correlation >= 0.20
        and rate >= 0.20 and interval["n_targets"] >= 30
        and interval["ci95"][0] is not None and interval["ci95"][0] > 0.0)
    final = make_pipeline(
        StandardScaler(), Ridge(alpha=CERTIFICATE_RIDGE)
    ).fit(features, margins)
    report = {
        "n_episodes": int(len(margins)),
        "crossfit_correlation": correlation,
        "crossfit_rmse": float(np.sqrt(np.mean((crossfit - margins) ** 2))),
        "conformal_error_quantile": quantile,
        "admission_rate": rate,
        "admitted_margin": interval,
        "source_gate_passed": passed,
    }
    return final, quantile, report


def _source_budget_ladder(y, bits, protein_views, scaffolds, protein_groups,
                          ligand_ridge, location_ridge, tangent_ridge):
    ligand_folds = np.asarray([stable_fold(value, 3) for value in scaffolds])
    protein_folds = np.asarray([stable_fold(value, 3) for value in protein_groups])
    records = {k: {"features": [], "margins": [], "folds": [], "targets": []}
               for k in K_VALUES}
    fold_details = []
    for fold in range(3):
        lt, lv = np.flatnonzero(ligand_folds != fold), np.flatnonzero(ligand_folds == fold)
        pt, pv = np.flatnonzero(protein_folds != fold), np.flatnonzero(protein_folds == fold)
        residual_train = double_center(y[np.ix_(lt, pt)])
        coefficient = _fit_ligand_atlas(bits[lt], residual_train, ligand_ridge)
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
        counts = {k: 0 for k in K_VALUES}
        for target in range(len(pv)):
            wrong_target = int(nearest[target])
            finite = np.flatnonzero(
                np.isfinite(outcome[:, target]) & np.isfinite(outcome[:, wrong_target]))
            for k in K_VALUES:
                for seed in BUDGET_SEEDS:
                    support = _d_optimal(
                        local_scaffold, surface[:, target, 1:3], finite, k,
                        "random", seed + 104729 * target)
                    if support is None:
                        continue
                    episode = _episode(
                        outcome, additive, interaction_y, surface, null_surface,
                        nearest, local_scaffold, target, support, seed,
                        location_ridge, tangent_ridge)
                    if episode is None:
                        continue
                    records[k]["features"].append(episode["feature"])
                    records[k]["margins"].append(episode["minimum_margin"])
                    records[k]["folds"].append(fold)
                    records[k]["targets"].append(int(pv[target]))
                    counts[k] += 1
        fold_details.append({
            "fold": fold, "valid_ligands": len(lv),
            "valid_targets": len(pv), "episodes": counts,
        })

    models, quantiles, reports = {}, {}, {}
    selected = None
    for k in K_VALUES:
        model, quantile, report = _fit_budget_certificate(
            records[k]["features"], records[k]["margins"],
            records[k]["folds"], records[k]["targets"])
        models[k], quantiles[k], reports[str(k)] = model, quantile, report
        print(f"k={k} {json.dumps(report)}", flush=True)
        if selected is None and report["source_gate_passed"]:
            selected = k
    return selected, models, quantiles, {
        "folds": fold_details, "budgets": reports,
    }


def _write_early_failure(args, started, source_report, selected_parameters):
    result = {
        "schema": "MetaSieve.IdentifiabilitySampleComplexity.F3K.v1",
        "selected": selected_parameters,
        "source": source_report,
        "transfers": {},
        "gate": {
            "source_budget_identified": False, "passed": False,
            "verdict": "F3K_NOT_IDENTIFIABLE_THROUGH_K40",
        },
        "read_firewall": {
            "pkis1": "source", "pkis2": "NOT_READ_BY_F3K",
            "anastassiadis": "NOT_READ_BY_F3K",
            "kcgs_numeric_outcomes": "NOT_READ", "davis_labels": "NOT_READ",
            "recipient_labels": "NOT_READ",
        },
        "elapsed_seconds": time.time() - started,
    }
    output = Path(args.output).resolve()
    _write_json(output / "result.json", result)
    _write_json(output / "manifest.json", {
        "sha256": {
            "preregistration": _sha256(Path(__file__).with_name("F3K_PREREGISTRATION.md")),
            "script": _sha256(Path(__file__)),
        }, "parameters": vars(args),
    })
    print(json.dumps(result["gate"], indent=2), flush=True)


def run(args):
    started = time.time()
    source = _load_source(args)
    y_source = source["y"]
    source_ligands = source["ligands"]
    source_targets = source["targets"]
    source_bits = _ligand_bits(source_ligands)
    source_scaffolds = np.asarray([
        record.generic_scaffold for record in source_ligands], dtype=object)
    source_groups = [item["record"].group for item in source_targets]
    masks = _subpocket_weights(args.kissim_distances)
    source_protein_views = _atlas_protein_views(source_targets, masks)

    print("selecting source-only atlas and section parameters", flush=True)
    ligand_ridge, ligand_cv = _select_ligand_ridge(
        y_source, source_bits, source_scaffolds)
    (location_ridge, tangent_ridge), section_cv = _cv_penalties(
        y_source, source_bits, source_protein_views, source_scaffolds,
        source_groups, ligand_ridge)
    print("evaluating source-only support-budget ladder", flush=True)
    selected_k, models, quantiles, ladder = _source_budget_ladder(
        y_source, source_bits, source_protein_views, source_scaffolds,
        source_groups, ligand_ridge, location_ridge, tangent_ridge)
    source_report = {
        "panel": "PKIS1", "ligand_cv": ligand_cv,
        "section_cv": section_cv, "budget_ladder": ladder,
    }
    selected_parameters = {
        "d_adapt": 4, "candidate_k": K_VALUES, "selected_k": selected_k,
        "ligand_ridge": ligand_ridge, "location_ridge": location_ridge,
        "tangent_ridge": tangent_ridge,
        "certificate_ridge": CERTIFICATE_RIDGE,
        "conformal_level": CONFORMAL_LEVEL,
    }
    if selected_k is None:
        _write_early_failure(
            args, started, source_report, selected_parameters)
        return

    print(f"source selected k={selected_k}; loading consumed development panels", flush=True)
    data = _load(args)
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
        print(f"evaluating selected budget on {panel['name']}", flush=True)
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
        evaluated = _evaluate_selective(
            panel["y"], additive, surface, null_surface, protein_views,
            scaffolds, selected_k, "random", SEEDS,
            location_ridge, tangent_ridge, models[selected_k], quantiles[selected_k])
        transfers[panel["name"]] = {
            "role": panel["role"], "n_ligands": len(panel["ligands"]),
            "n_targets": len(panel["targets"]), "support_size": selected_k,
            "random": _selective_summary(*evaluated),
        }

    primary = transfers["PKIS2"]["random"]
    external = transfers["Anastassiadis2011"]["random"]
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
        "schema": "MetaSieve.IdentifiabilitySampleComplexity.F3K.v1",
        "selected": selected_parameters,
        "source": source_report, "transfers": transfers,
        "gate": {
            "source_budget_identified": True,
            "pkis2_coverage": pkis_coverage,
            "anastassiadis_coverage": external_coverage,
            "pkis2_raw_controls": pkis_raw,
            "anastassiadis_raw_point_estimates": external_raw,
            "pkis2_interaction_controls": pkis_interaction,
            "anastassiadis_interaction_point_estimates": external_interaction,
            "passed": passed,
            "verdict": ("F3K_SAMPLE_COMPLEXITY_ADMISSIBLE" if passed
                        else "F3K_SELECTED_BUDGET_NOT_TRANSFERABLE"),
        },
        "read_firewall": {
            "pkis1": "source", "pkis2": "consumed_development",
            "anastassiadis": "consumed_development",
            "kcgs_numeric_outcomes": "NOT_READ", "davis_labels": "NOT_READ",
            "recipient_labels": "NOT_READ",
        },
        "elapsed_seconds": time.time() - started,
    }
    output = Path(args.output).resolve()
    _write_json(output / "result.json", result)
    _write_json(output / "manifest.json", {
        "sha256": {
            "preregistration": _sha256(Path(__file__).with_name("F3K_PREREGISTRATION.md")),
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
    item.add_argument("--output", default="research/section_operator_pilot/artifacts/f3k")
    return item


if __name__ == "__main__":
    run(parser().parse_args())
