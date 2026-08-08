"""Run the corrected research E-AFF-L0R protein affinity-location Gate.

Everything here is frozen by `EAFF_L0R_PREREGISTRATION.md`,
`EAFF_L0_PREREGISTRATION.md`,
`EAFF_L0_DATA_CONTRACT.md` and `l0_contract.py`. Ki only, per the Stage 4
identifiability verdict. Research-only: nothing enters `model/` or `scripts/`.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path
import platform

import numpy as np
import torch

from model.config import MetaSieveConfig
from model.mathematical import (
    band_loss,
    target_from_conditional_cdf,
)
from model import bands
from model.meta_operator import build_anchors
from research.e0_identifiability.l0_contract import (
    knn_conditional_cdf,
    Z_BIO_GAUGE_COORDINATES,
    Z_BIO_NAMES,
    assert_bounded_observable,
    band_mean_interval,
    bounded_projection,
    component_bootstrap,
    component_macro,
    pooled_replicate_sigma,
    sigma_confidence,
    z_bio,
)
from research.e0_identifiability.run_eaff_pilot import (
    _geometry_basis,
    _hash_key,
    _load_ligand_cache,
    _load_protein_cache,
    _read_jsonl,
    _write_json,
    _write_jsonl,
    build_derangement,
    load_affinity,
    make_observations,
)
from research.e0_identifiability.run_synthetic_pregate import _load_bridge
from scripts.source_affinity.common import sha256_file


STAGE = "P1R2B-E-AFF-L0R_CORRECTED_PROTEIN_AFFINITY_LOCATION_GATE"
ENDPOINT = "Ki"
MAX_TASKS = 195
TASKS_PER_COMPONENT = 3
POSITIVE_CONTROL_RATIO = 0.1
LIGANDS_PER_TASK = 20
FEATURE_DIM = 7
FOLDS = 5
MIN_REPLICATE_CELLS = 100
COVERAGE_TOLERANCE = 0.05
MARGIN_RATIO = 0.5
SEED = 20260808
DKW_EPS_MIN = 0.02
DKW_N_MIN = 30


def select_tasks(rows: list[dict], consumed: set[str]) -> tuple[list[dict], dict]:
    """Label-blind selection: Ki tasks never consumed by P0, H0A or H0C."""
    by_task: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    meta: dict[str, dict] = {}
    for row in rows:
        if row["endpoint_family"] != ENDPOINT or row["task_id"] in consumed:
            continue
        by_task[row["task_id"]][row["ligand_state_key"]].append(row)
        meta.setdefault(row["task_id"], row)
    eligible = [task for task, ligands in by_task.items() if len(ligands) >= LIGANDS_PER_TASK]
    eligible.sort(key=lambda task: _hash_key("L0-TASK", task))
    chosen, per_component = [], defaultdict(int)
    for task in eligible:
        component = meta[task]["closure_component_id"]
        if per_component[component] >= TASKS_PER_COMPONENT:
            continue
        per_component[component] += 1
        chosen.append(task)
        if len(chosen) >= MAX_TASKS:
            break
    used_components = set(per_component)
    selected = []
    for task in chosen:
        ligands = sorted(by_task[task], key=lambda value: _hash_key("L0-LIGAND", task, value))
        for ligand in ligands[:LIGANDS_PER_TASK]:
            selected.append(min(by_task[task][ligand],
                                key=lambda value: int(value["activity_id"])))
    audit = {
        "eligible_tasks": len(eligible),
        "selected_tasks": len(chosen),
        "closure_components": len(used_components),
        "consumed_tasks_excluded": len(consumed),
        "ligands_per_task": LIGANDS_PER_TASK,
    }
    return selected, audit


def population_band(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    if len(values) < DKW_N_MIN:
        return bands.join(np.zeros(len(grid)), np.ones(len(grid)))
    ecdf = (values[:, None] <= grid[None, :]).mean(axis=0)
    eps = max(float(np.sqrt(np.log(2.0 / 0.10) / (2.0 * len(values)))), DKW_EPS_MIN)
    return bands.band_from_ecdf(ecdf, eps)


def run(args) -> dict:
    output = Path(args.output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"output exists and is not empty: {output}")
    output.mkdir(parents=True, exist_ok=True)
    if not args.device.startswith("cuda") or not torch.cuda.is_available():
        return fail_closed(output, "L0R_NOT_RUN_NUMERICAL_PRECONDITION_FAILED",
                           "CUDA is required for the frozen P1B bridge")

    cfg = MetaSieveConfig()
    grid = cfg.grid()
    input_root = Path(args.input)
    estimand = json.loads(
        (Path(args.estimand) / "report.json").read_text(encoding="utf-8"))
    if ENDPOINT not in estimand["admitted_endpoints"]:
        return fail_closed(output, "L0R_NOT_RUN_LOCATION_ESTIMAND_NOT_IDENTIFIED",
                           f"{ENDPOINT} was not admitted by the estimand contract")
    contract = json.loads(
        (Path(args.contract) / "report.json").read_text(encoding="utf-8"))
    if contract["verdict"] != "L0_OPERATOR_AND_ANCHOR_CONTRACT_FROZEN":
        return fail_closed(output, "L0R_NOT_RUN_OPERATOR_OR_ANCHOR_CONTRACT_FAILED",
                           contract["verdict"])

    all_rows = list(_read_jsonl(input_root / "rows.label_blind.jsonl"))
    consumed = set()
    for name in args.consumed:
        path = Path(name)
        if path.is_file():
            consumed.update(row["task_id"] for row in _read_jsonl(path))
    selected, selection_audit = select_tasks(all_rows, consumed)
    if selection_audit["selected_tasks"] < 50:
        return fail_closed(output, "L0R_NOT_RUN_NUMERICAL_PRECONDITION_FAILED",
                           "fewer than 50 unconsumed Ki tasks are available")
    _write_jsonl(output / "selection.jsonl", selected)

    labels, label_audit = load_affinity(all_rows, Path(args.canonical_rows))
    replicate_groups = defaultdict(list)
    for row in all_rows:
        if row["endpoint_family"] != ENDPOINT:
            continue
        replicate_groups[(row["task_id"], row["ligand_state_key"])].append(
            labels[int(row["activity_id"])]["p_affinity"])
    sigma, replicate_cells, degrees = pooled_replicate_sigma(
        [np.asarray(values, dtype=np.float64) for values in replicate_groups.values()])
    if replicate_cells < MIN_REPLICATE_CELLS or not np.isfinite(sigma):
        return fail_closed(output, "L0R_NOT_RUN_NUMERICAL_PRECONDITION_FAILED",
                           f"only {replicate_cells} replicate cells available")
    sigma_low, sigma_high = sigma_confidence(sigma, degrees)
    margin = MARGIN_RATIO * sigma

    observations = make_observations(all_rows, labels)
    lookup = {(row["task_id"], row["ligand_state_key"]): index
              for index, row in enumerate(observations)}
    rows = [observations[lookup[(row["task_id"], row["ligand_state_key"])]]
            for row in selected]
    y_log = np.asarray([row["p_affinity"] for row in rows], dtype=np.float64)
    folds = np.asarray([int(row["outer_oof_fold"]) for row in rows])
    components = np.asarray([row["closure_component_id"] for row in rows])

    # frozen affine scaling of the affinity range onto V, from fitting data only
    low_q, high_q = np.quantile(y_log, [0.01, 0.99])
    span = float(high_q - low_q) if high_q > low_q else 1.0
    y_scaled = np.clip((y_log - low_q) / span, 0.0, 1.0)

    ligand_keys = {row["ligand_state_key"] for row in rows}
    pooled, local_ligands, _ = _load_ligand_cache(Path(args.cache), ligand_keys, ligand_keys)
    wrong_map, wrong_audit = build_derangement(rows, input_root, Path(args.governance))
    _write_jsonl(output / "derangement.jsonl", wrong_audit)
    protein_keys = ({row["protein_sequence_sha256"] for row in rows} | set(wrong_map.values()))
    proteins = _load_protein_cache(Path(args.cache), protein_keys)
    bridge = _load_bridge(Path(args.checkpoint), args.device)
    with np.load(args.tbasis_values, allow_pickle=False) as tbasis:
        correct_phi, deranged_phi = _geometry_basis(
            rows, bridge, proteins, local_ligands, wrong_map,
            tbasis["bin_rbf_expectation"].astype(np.float64),
            tbasis["calibration_coef"].astype(np.float64),
            tbasis["calibration_intercept"].astype(np.float64), args.device)

    heavy = np.asarray([len(local_ligands[row["ligand_state_key"]]["atoms"])
                        for row in rows], dtype=np.float64)
    z_correct = z_bio(correct_phi, heavy)
    z_deranged = z_bio(deranged_phi, heavy)
    assert_bounded_observable(z_correct, "z_bio(correct)")
    assert_bounded_observable(z_deranged, "z_bio(deranged)")

    ligand_state = np.stack([pooled[row["ligand_state_key"]] for row in rows]).astype(np.float64)
    protein_state = np.stack([
        proteins[row["protein_sequence_sha256"]]["residues"][
            proteins[row["protein_sequence_sha256"]]["mask"]].mean(axis=0)
        for row in rows]).astype(np.float64)

    anchors = build_anchors(cfg, device="cpu").numpy()
    arm_names = ("A0", "A1", "A2", "A3", "A4")
    predictions = {name: np.zeros((len(rows), cfg.band_dim)) for name in arm_names}
    p0_mass = {name: np.zeros(len(rows)) for name in arm_names}
    diagnostics = []

    for fold in range(FOLDS):
        train = folds != fold
        test = folds == fold
        if not test.any() or train.sum() < DKW_N_MIN:
            continue
        pop = population_band(y_scaled[train], grid)
        columns = np.concatenate([pop[None, :], anchors], axis=0).T  # (2G, m+1)

        features = {
            "A1": bounded_projection(ligand_state[train], ligand_state, FEATURE_DIM),
            "A2": bounded_projection(protein_state[train], protein_state, FEATURE_DIM),
            "A3": z_correct,
            "A4": z_deranged,
        }
        neighbours = int(np.ceil(np.sqrt(train.sum())))
        for name in arm_names:
            if name == "A0":
                p = np.zeros(cfg.m + 1)
                p[0] = 1.0
                predictions[name][test] = (columns @ p)[None, :]
                p0_mass[name][test] = 1.0
                continue
            estimator_features = features["A3"] if name == "A4" else features[name]
            query = knn_conditional_cdf(estimator_features[train], y_scaled[train],
                                        features[name][test], neighbours, grid)
            coefficients = np.stack([
                target_from_conditional_cdf(columns, row, cfg.lambda_w, cfg.mu)
                for row in query])
            predictions[name][test] = coefficients @ columns.T
            p0_mass[name][test] = coefficients[:, 0]
        diagnostics.append({"fold": fold, "train": int(train.sum()), "test": int(test.sum()),
                            "neighbours": neighbours,
                            "population_band_source_rows": int(train.sum())})

    metrics = {}
    for name in arm_names:
        beta = predictions[name]
        bands.assert_valid(beta, name=f"L0 emitted band {name}")
        loss = band_loss(beta, y_scaled, grid, cfg.lambda_w)
        lower, upper = bands.split(beta)
        interval = band_mean_interval(beta, grid, cfg.a_max)
        location = 0.5 * (interval[:, 0] + interval[:, 1])
        location_error = np.abs(location - y_scaled) * span
        coverage = float(np.mean((y_scaled >= interval[:, 0]) & (y_scaled <= interval[:, 1])))
        metrics[name] = {
            "band_loss": component_macro(loss, components),
            "band_loss_ci95": component_bootstrap(loss, components),
            "location_error_log_units": component_macro(location_error, components),
            "coverage": coverage,
            "mean_interval_width": float(np.mean(upper - lower)),
            "abstention_p0_mass": float(np.mean(p0_mass[name])),
            "_loss": loss,
            "_location_error": location_error,
        }

    control_delta = (metrics["A0"]["_location_error"] - metrics["A1"]["_location_error"])
    control_mean = component_macro(control_delta, components)
    control_ci = component_bootstrap(control_delta, components)
    control_required = POSITIVE_CONTROL_RATIO * sigma
    positive_control = {
        "location_error_gain_A0_minus_A1": control_mean,
        "ci95": control_ci,
        "required": control_required,
        "passed": bool(control_mean >= control_required and control_ci[0] > 0.0),
    }

    contrasts = {}
    passes = {}
    for control in ("A1", "A2", "A4"):
        loss_delta = metrics[control]["_loss"] - metrics["A3"]["_loss"]
        location_delta = metrics[control]["_location_error"] - metrics["A3"]["_location_error"]
        loss_ci = component_bootstrap(loss_delta, components)
        location_ci = component_bootstrap(location_delta, components)
        loss_mean = component_macro(loss_delta, components)
        location_mean = component_macro(location_delta, components)
        coverage_ok = metrics["A3"]["coverage"] >= metrics[control]["coverage"] - COVERAGE_TOLERANCE
        contrasts[f"A3_vs_{control}"] = {
            "band_loss_reduction": loss_mean,
            "band_loss_reduction_ci95": loss_ci,
            "location_error_reduction_log_units": location_mean,
            "location_error_reduction_ci95": location_ci,
            "margin_required": margin,
            "condition_1_band_loss": bool(loss_mean > 0.0 and loss_ci[0] > 0.0),
            "condition_2_location_margin": bool(location_mean >= margin and location_ci[0] > 0.0),
            "condition_3_coverage_guard": bool(coverage_ok),
        }
        passes[control] = all(contrasts[f"A3_vs_{control}"][key] for key in
                              ("condition_1_band_loss", "condition_2_location_margin",
                               "condition_3_coverage_guard"))

    if not positive_control["passed"]:
        verdict = "L0R_NOT_RUN_POSITIVE_CONTROL_ABSENT"
    elif all(passes.values()):
        verdict = "PROTEIN_SPECIFIC_AFFINITY_LOCATION_IDENTIFIED_IN_SOURCE"
    else:
        verdict = "PROTEIN_SPECIFIC_AFFINITY_LOCATION_NOT_IDENTIFIED"

    report = {
        "schema": "MetaSieve.EAffL0.v1",
        "stage": STAGE,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": verdict,
        "endpoint": ENDPOINT,
        "endpoint_not_identified": [name for name in ("Ki", "Kd")
                                    if name not in estimand["admitted_endpoints"]],
        "affinity_labels_read": True,
        "label_fields_accessed": ["p_affinity"],
        "davis_label_reads": 0,
        "recipient_label_reads": 0,
        "selection": selection_audit,
        "label_audit": label_audit,
        "sigma_assay": {
            "estimator": "pooled within-(task,ligand) replicate SD, weighted by n-1",
            "value": sigma, "ci95": [sigma_low, sigma_high],
            "replicate_cells": replicate_cells, "degrees_of_freedom": degrees,
            "margin_ratio": MARGIN_RATIO, "margin_L0": margin,
        },
        "affinity_scaling": {"q01": float(low_q), "q99": float(high_q), "span": span},
        "z_bio": {"names": list(Z_BIO_NAMES), "gauge_coordinates": list(Z_BIO_GAUGE_COORDINATES),
                  "bounds": [0.0, 1.0], "dimension": FEATURE_DIM},
        "arms": {name: {key: value for key, value in metrics[name].items()
                        if not key.startswith("_")} for name in arm_names},
        "positive_control": positive_control,
        "contrasts": contrasts,
        "arm_passes": passes,
        "folds": diagnostics,
        "rows": len(rows),
        "interpretation_limits": [
            "L0 addresses Claim A only and is not evidence about protein-by-ligand interaction",
            "a pass establishes PROTEIN_SPECIFIC_AFFINITY_LOCATION_IDENTIFIED_IN_SOURCE and nothing else",
            "no coordinate is a physical binding free energy",
            "ligand_size_gauge is a declared gauge coordinate",
        ],
    }
    _write_json(output / "report.json", report)
    _write_json(output / "manifest.json", {
        "stage": STAGE,
        "preregistration": sha256_file(
            Path(__file__).with_name("EAFF_L0R_PREREGISTRATION.md")),
        "l0_preregistration": sha256_file(
            Path(__file__).with_name("EAFF_L0_PREREGISTRATION.md")),
        "data_contract": sha256_file(Path(__file__).with_name("EAFF_L0_DATA_CONTRACT.md")),
        "l0_contract": sha256_file(Path(__file__).with_name("l0_contract.py")),
        "code": sha256_file(Path(__file__)),
        "inputs": {
            "label_blind_rows": sha256_file(input_root / "rows.label_blind.jsonl"),
            "canonical_rows": sha256_file(Path(args.canonical_rows)),
            "checkpoint": sha256_file(Path(args.checkpoint)),
            "tbasis_values": sha256_file(Path(args.tbasis_values)),
            "estimand_report": sha256_file(Path(args.estimand) / "report.json"),
            "contract_report": sha256_file(Path(args.contract) / "report.json"),
        },
        "outputs": {name: sha256_file(output / name) for name in
                    ("selection.jsonl", "derangement.jsonl", "report.json")},
        "seed": SEED,
        "environment": {
            "python": platform.python_version(), "torch": torch.__version__,
            "cuda": torch.version.cuda, "device": args.device,
            "numpy": np.__version__,
        },
        "label_reads": {"affinity_values": len(rows), "davis": 0, "recipient": 0},
    })
    return report


def fail_closed(output: Path, verdict: str, reason: str) -> dict:
    report = {"schema": "MetaSieve.EAffL0.v1", "stage": STAGE, "verdict": verdict,
              "reason": reason, "created_utc": datetime.now(timezone.utc).isoformat(),
              "davis_label_reads": 0, "recipient_label_reads": 0}
    _write_json(output / "report.json", report)
    return report


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="dataset/processed/source_affinity/e0_input_v1")
    parser.add_argument("--canonical-rows",
                        default="dataset/processed/source_affinity/energy_pilot_v1/canonical_rows.jsonl")
    parser.add_argument("--cache", default="dataset/processed/source_affinity/e0_local_states_v1")
    parser.add_argument("--governance",
                        default="dataset/processed/source_affinity/energy_pilot_v1_governance")
    parser.add_argument("--checkpoint",
                        default="report/mechanism_refactor/p1b_pilot20k_seed17_v1/best.pt")
    parser.add_argument("--tbasis-values",
                        default="research/e0_identifiability/artifacts/tbasis_r0_v1/basis_values.npz")
    parser.add_argument("--estimand",
                        default="research/e0_identifiability/artifacts/eaff_l0_estimand_v1")
    parser.add_argument("--contract",
                        default="research/e0_identifiability/artifacts/eaff_l0_contract_v1")
    parser.add_argument("--consumed", nargs="*", default=[
        "research/e0_identifiability/artifacts/eaff_p0_v1/selection.jsonl",
        "research/e0_identifiability/artifacts/eaff_h0a_v1/selection.jsonl",
        "research/e0_identifiability/artifacts/eaff_h0c_v1_run2/selection.jsonl",
        "research/e0_identifiability/artifacts/eaff_l0_v1/selection.jsonl",
    ])
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", default="research/e0_identifiability/artifacts/eaff_l0r_v1")
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), indent=2, sort_keys=True))
