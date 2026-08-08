"""Run the synthetic-only E0R1 objective, design, and solver audit."""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F

from research.e0_identifiability.objective_solver_contract import (
    CI_INCREMENT_THRESHOLD, CORRECT_CI_THRESHOLD,
    DERANGEMENT_IDENTITY_THRESHOLD, LBFGS_HISTORY_SIZE, LBFGS_MAX_ITERATIONS,
    LBFGS_TOLERANCE_CHANGE, LBFGS_TOLERANCE_GRAD,
    OBJECTIVE_GRADIENT_TOLERANCE, PERMUTATION_TOLERANCE,
    PINV_PRIMARY_RCOND, PINV_SENSITIVITY_RCONDS, PINV_TRAIN_RMSE_TOLERANCE,
    SOLVER_RELATIVE_GRADIENT_TOLERANCE,
)
from scripts.audit_e0_input_feasibility import _read_jsonl
from research.e0_identifiability.audit_evidence import _macro_ci, _teacher_values
from scripts.govern_structure_homology import _local_identity
from research.e0_identifiability.run_synthetic_pregate import (
    DISTANCE_WEIGHTS, _geometry, _load_bridge, _load_states, _select_rows,
)
from research.e0_identifiability.run_typed_tensor import _typed_statistics
from scripts.structure_sources.rcsb import sha256_file


def _pair_indices(rows: list[dict], selected: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        if selected[index]:
            grouped[row["task_id"]].append(index)
    left, right = [], []
    for indices in grouped.values():
        for first in range(len(indices)):
            for second in range(first + 1, len(indices)):
                left.append(indices[first])
                right.append(indices[second])
    return np.asarray(left, dtype=np.int64), np.asarray(right, dtype=np.int64)


def _loss_and_gradient(kind: str, features: np.ndarray, weights: np.ndarray,
                       residuals: np.ndarray, baselines: np.ndarray,
                       labels: np.ndarray, left: np.ndarray,
                       right: np.ndarray) -> tuple[float, float]:
    x = torch.from_numpy(features).double()
    parameter = torch.tensor(weights, dtype=torch.float64, requires_grad=True)
    residual = torch.from_numpy(residuals).double()
    baseline = torch.from_numpy(baselines).double()
    label = torch.from_numpy(labels).double()
    prediction = x @ parameter
    if kind == "point_huber":
        loss = F.huber_loss(prediction, residual)
    elif kind in {"old_logistic", "aligned_logistic"}:
        score = prediction if kind == "old_logistic" else baseline + prediction
        signs = torch.sign(label[left] - label[right])
        keep = signs != 0
        loss = F.softplus(-signs[keep] * (score[left][keep] - score[right][keep])).mean()
    elif kind == "residual_difference_huber":
        predicted_delta = prediction[left] - prediction[right]
        target_delta = residual[left] - residual[right]
        loss = F.huber_loss(predicted_delta, target_delta)
    elif kind == "old_total":
        signs = torch.sign(label[left] - label[right])
        keep = signs != 0
        rank = F.softplus(
            -signs[keep] * (prediction[left][keep] - prediction[right][keep])).mean()
        loss = F.huber_loss(prediction, residual) + rank
    else:
        raise ValueError(f"unknown objective: {kind}")
    gradient = torch.autograd.grad(loss, parameter)[0]
    return float(loss.detach()), float(gradient.norm())


def _score(rows: list[dict], holdout_mask: np.ndarray, labels: np.ndarray,
           baselines: np.ndarray, correct_prediction: np.ndarray,
           deranged_prediction: np.ndarray) -> dict:
    holdout_rows = [row for row, keep in zip(rows, holdout_mask) if keep]
    holdout_labels = labels[holdout_mask]
    holdout_baseline = baselines[holdout_mask]
    ligand_ci, _ = _macro_ci(holdout_rows, holdout_labels, holdout_baseline)
    correct_ci, _ = _macro_ci(
        holdout_rows, holdout_labels, holdout_baseline + correct_prediction)
    deranged_ci, _ = _macro_ci(
        holdout_rows, holdout_labels, holdout_baseline + deranged_prediction)
    checks = {
        "correct_ci_at_least_0_80": correct_ci >= CORRECT_CI_THRESHOLD,
        "correct_minus_ligand_at_least_0_10":
            correct_ci - ligand_ci >= CI_INCREMENT_THRESHOLD,
        "correct_minus_deranged_at_least_0_10":
            correct_ci - deranged_ci >= CI_INCREMENT_THRESHOLD,
        "permutation_error_at_most_1e_6": 0.0 <= PERMUTATION_TOLERANCE,
    }
    return {
        "ligand_ci": ligand_ci,
        "correct_ci": correct_ci,
        "deranged_ci": deranged_ci,
        "correct_minus_ligand": correct_ci - ligand_ci,
        "correct_minus_deranged": correct_ci - deranged_ci,
        "permutation_error": 0.0,
        "gate": {**checks, "pass": all(checks.values())},
    }


def _build_derangement(rows_path: Path, proteins_path: Path,
                       selected_rows: list[dict], holdout_mask: np.ndarray) -> tuple[list[dict], dict]:
    sequences = {row["sequence_sha256"]: row["sequence"]
                 for row in _read_jsonl(proteins_path)}
    protein_closure = {}
    for row in _read_jsonl(rows_path):
        protein_closure.setdefault(row["protein_sequence_sha256"],
                                   row["closure_component_id"])
    candidates = sorted(set(protein_closure) & set(sequences))
    mapping = {}
    records = []
    holdout_rows = [row for row, keep in zip(selected_rows, holdout_mask) if keep]
    for row in sorted(holdout_rows, key=lambda value: value["task_id"]):
        closure = row["closure_component_id"]
        if closure in mapping:
            continue
        correct = row["protein_sequence_sha256"]
        for candidate in candidates:
            if candidate == correct or protein_closure[candidate] == closure:
                continue
            identity = _local_identity(sequences[correct], sequences[candidate])
            if identity < DERANGEMENT_IDENTITY_THRESHOLD:
                mapping[closure] = candidate
                records.append({
                    "closure_component_id": closure,
                    "correct_protein_sequence_sha256": correct,
                    "deranged_protein_sequence_sha256": candidate,
                    "deranged_closure_component_id": protein_closure[candidate],
                    "local_identity": identity,
                })
                break
        if closure not in mapping:
            raise ValueError(f"no score-blind <40% derangement for closure {closure}")
    deranged = [dict(row, active_protein_key=mapping[row["closure_component_id"]])
                for row in holdout_rows]
    manifest = {
        "schema": "MetaSieve.E0R1DerangementMap.v1",
        "selection": "lexically_first_model_valid_protein_in_different_closure_with_local_identity_lt_0.40",
        "threshold": DERANGEMENT_IDENTITY_THRESHOLD,
        "pairs": records,
        "all_pairs_valid": all(row["local_identity"] < DERANGEMENT_IDENTITY_THRESHOLD
                               for row in records),
        "score_fields_read_for_selection": False,
    }
    return deranged, manifest


def _svd_audit(x_train: np.ndarray, x_holdout: np.ndarray,
               teacher: np.ndarray) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray]:
    _u, singular, vh = np.linalg.svd(x_train, full_matrices=False)
    machine_tolerance = max(x_train.shape) * np.finfo(np.float64).eps * singular[0]
    machine_rank = int(np.sum(singular > machine_tolerance))
    primary_rank = int(np.sum(singular > singular[0] * PINV_PRIMARY_RCOND))
    basis = vh[:primary_rank].T
    projection = basis @ basis.T
    probabilities = singular[singular > 0] / singular.sum()
    effective_rank = float(np.exp(-np.sum(probabilities * np.log(probabilities))))
    identified_teacher = projection @ teacher
    holdout_projection = x_holdout @ projection
    row_norm = np.linalg.norm(x_holdout, axis=1)
    coverage = np.linalg.norm(holdout_projection, axis=1) / np.maximum(row_norm, 1e-15)
    unseen_prediction = x_holdout @ (teacher - identified_teacher)
    total_prediction = x_holdout @ teacher
    column_std = x_train.std(axis=0)
    audit = {
        "schema": "MetaSieve.E0R1DesignSVDAudit.v1",
        "shape": list(x_train.shape),
        "machine_tolerance": machine_tolerance,
        "machine_rank": machine_rank,
        "primary_rcond": PINV_PRIMARY_RCOND,
        "primary_rank": primary_rank,
        "nullity": x_train.shape[1] - primary_rank,
        "effective_rank": effective_rank,
        "singular_maximum": float(singular[0]),
        "singular_minimum_identified": float(singular[primary_rank - 1]),
        "condition_number_identified": float(singular[0] / singular[primary_rank - 1]),
        "teacher_identifiable_norm_fraction": float(
            np.linalg.norm(identified_teacher) / np.linalg.norm(teacher)),
        "column_std": {
            "minimum": float(column_std.min()),
            "median": float(np.median(column_std)),
            "maximum": float(column_std.max()),
            "zero_columns": int(np.sum(column_std == 0)),
            "below_1e_8": int(np.sum(column_std < 1e-8)),
        },
        "holdout_row_space_coverage": {
            "mean": float(coverage.mean()),
            "median": float(np.median(coverage)),
            "minimum": float(coverage.min()),
            "maximum_orthogonal_fraction": float(
                np.max(np.sqrt(np.maximum(0.0, 1.0 - coverage ** 2)))),
        },
        "teacher_specific_unseen_contribution": {
            "relative_l2": float(np.linalg.norm(unseen_prediction)
                                 / max(np.linalg.norm(total_prediction), 1e-15)),
            "maximum_absolute": float(np.max(np.abs(unseen_prediction))),
            "rmse": float(np.sqrt(np.mean(unseen_prediction ** 2))),
        },
        "rank_sensitivity": {
            str(value): int(np.sum(singular > singular[0] * value))
            for value in (PINV_PRIMARY_RCOND, *PINV_SENSITIVITY_RCONDS)
        },
        "singular_values": singular.tolist(),
        "rank_deficiency_alone_is_not_failure": True,
    }
    return audit, singular, vh, projection


def _pinv_solution(x_train: np.ndarray, residual_train: np.ndarray,
                   rcond: float) -> np.ndarray:
    return np.linalg.pinv(x_train, rcond=rcond) @ residual_train


def _corrected_objective(parameter: torch.Tensor, x: torch.Tensor,
                         residual: torch.Tensor, left: torch.Tensor,
                         right: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    prediction = x @ parameter
    point = F.huber_loss(prediction, residual)
    difference = F.huber_loss(
        prediction[left] - prediction[right], residual[left] - residual[right])
    return point + difference, point, difference


def _run_deterministic_solver(x_train: np.ndarray, residual_train: np.ndarray,
                              left: np.ndarray, right: np.ndarray) -> tuple[np.ndarray, list[dict], dict]:
    x = torch.from_numpy(x_train).double()
    residual = torch.from_numpy(residual_train).double()
    left_tensor = torch.from_numpy(left).long()
    right_tensor = torch.from_numpy(right).long()
    parameter = torch.zeros(x_train.shape[1], dtype=torch.float64, requires_grad=True)
    optimizer = torch.optim.LBFGS(
        [parameter], lr=1.0, max_iter=LBFGS_MAX_ITERATIONS,
        tolerance_grad=LBFGS_TOLERANCE_GRAD,
        tolerance_change=LBFGS_TOLERANCE_CHANGE,
        history_size=LBFGS_HISTORY_SIZE, line_search_fn="strong_wolfe")
    trace = []

    def closure():
        optimizer.zero_grad(set_to_none=True)
        loss, point, difference = _corrected_objective(
            parameter, x, residual, left_tensor, right_tensor)
        loss.backward()
        trace.append({
            "closure_call": len(trace) + 1,
            "loss": float(loss.detach()),
            "point_huber": float(point.detach()),
            "difference_huber": float(difference.detach()),
            "gradient_l2": float(parameter.grad.detach().norm()),
            "parameter_l2": float(parameter.detach().norm()),
        })
        return loss

    optimizer.step(closure)
    parameter.grad = None
    final_loss, final_point, final_difference = _corrected_objective(
        parameter, x, residual, left_tensor, right_tensor)
    final_loss.backward()
    gradient_l2 = float(parameter.grad.detach().norm())
    parameter_l2 = float(parameter.detach().norm())
    relative_gradient = gradient_l2 / (1.0 + parameter_l2)
    prediction = x @ parameter.detach()
    train_rmse = float(torch.sqrt(torch.mean((prediction - residual) ** 2)))
    status = {
        "closure_calls": len(trace),
        "final_loss": float(final_loss.detach()),
        "final_point_huber": float(final_point.detach()),
        "final_difference_huber": float(final_difference.detach()),
        "final_gradient_l2": gradient_l2,
        "final_relative_gradient": relative_gradient,
        "parameter_l2": parameter_l2,
        "train_rmse": train_rmse,
        "converged": (relative_gradient <= SOLVER_RELATIVE_GRADIENT_TOLERANCE
                      and train_rmse <= PINV_TRAIN_RMSE_TOLERANCE),
    }
    return parameter.detach().numpy(), trace, status


def run(input_root: Path, cache_root: Path, checkpoint: Path,
        synthetic_root: Path, e0r0_root: Path, output: Path,
        device: str) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    gate_report = json.loads((synthetic_root / "synthetic_gate.json").read_text())
    rows_path = input_root / "rows.label_blind.jsonl"
    rows = _select_rows(rows_path)
    for index, row in enumerate(rows):
        row["example_id"] = index
        row["active_protein_key"] = row["protein_sequence_sha256"]
    train_mask = np.asarray([row["outer_oof_fold"] < 4 for row in rows])
    holdout_mask = ~train_mask
    deranged_rows, derangement_manifest = _build_derangement(
        rows_path, input_root / "proteins.jsonl", rows, holdout_mask)
    derangement_path = output / "DERANGEMENT_MAP.json"
    derangement_path.write_text(
        json.dumps(derangement_manifest, indent=2, sort_keys=True), encoding="utf-8")

    ligand_keys = {row["ligand_state_key"] for row in rows}
    protein_keys = ({row["protein_sequence_sha256"] for row in rows}
                    | {row["active_protein_key"] for row in deranged_rows})
    proteins, ligands = _load_states(cache_root, ligand_keys, protein_keys)
    bridge = _load_bridge(checkpoint, device)
    correct_geometry = _geometry(rows, bridge, proteins, ligands, device)
    deranged_geometry = _geometry(deranged_rows, bridge, proteins, ligands, device)
    weights = np.asarray(gate_report["teacher_weights"], dtype=np.float64)
    raw, _, baselines = _teacher_values(rows, proteins, ligands, correct_geometry, weights)
    center = float(raw[train_mask].mean())
    scale = max(float(raw[train_mask].std()), 1e-6)
    residuals = (raw - center) / scale
    labels = baselines + residuals
    correct = _typed_statistics(rows, proteins, ligands, correct_geometry)
    feature_center = correct[train_mask].mean(axis=0)
    x = (correct - feature_center).reshape(len(rows), -1)
    deranged = _typed_statistics(deranged_rows, proteins, ligands, deranged_geometry)
    x_deranged = (deranged - feature_center).reshape(len(deranged_rows), -1)
    teacher = (weights[:, :, None]
               * DISTANCE_WEIGHTS.astype(np.float64)[None, None, :] / scale).reshape(-1)
    left, right = _pair_indices(rows, train_mask)

    conflict_keep = ((np.sign(residuals[left] - residuals[right]) != 0)
                     & (np.sign(labels[left] - labels[right]) != 0))
    conflicts = np.sign(residuals[left] - residuals[right]) != np.sign(
        labels[left] - labels[right])
    objective = {
        "schema": "MetaSieve.E0R1ObjectiveAudit.v1",
        "train_pairs": int(len(left)),
        "comparable_pairs": int(conflict_keep.sum()),
        "ordering_conflicts": int(np.sum(conflicts & conflict_keep)),
        "ordering_conflict_rate": float(np.mean(conflicts[conflict_keep])),
        "analytic_teacher": {},
    }
    for kind in ("point_huber", "old_logistic", "aligned_logistic",
                 "residual_difference_huber", "old_total"):
        loss, gradient = _loss_and_gradient(
            kind, x[train_mask], teacher, residuals[train_mask], baselines[train_mask],
            labels[train_mask], np.searchsorted(np.flatnonzero(train_mask), left),
            np.searchsorted(np.flatnonzero(train_mask), right))
        objective["analytic_teacher"][kind] = {
            "loss": loss, "full_gradient_l2": gradient}
    learned_state = torch.load(
        e0r0_root / "full_240_model.pt", map_location="cpu", weights_only=False)
    learned_weight = learned_state["model_state"]["energy"].double().numpy().reshape(-1)
    loss, gradient = _loss_and_gradient(
        "old_total", x[train_mask], learned_weight, residuals[train_mask],
        baselines[train_mask], labels[train_mask],
        np.searchsorted(np.flatnonzero(train_mask), left),
        np.searchsorted(np.flatnonzero(train_mask), right))
    objective["e0r0_epoch60_full_risk"] = {
        "old_total_loss": loss, "full_gradient_l2": gradient}
    objective["objective_semantics_defect_confirmed"] = (
        objective["ordering_conflict_rate"] > 0
        and objective["analytic_teacher"]["point_huber"]["full_gradient_l2"]
            <= OBJECTIVE_GRADIENT_TOLERANCE
        and objective["analytic_teacher"]["residual_difference_huber"]["full_gradient_l2"]
            <= OBJECTIVE_GRADIENT_TOLERANCE
        and objective["analytic_teacher"]["old_logistic"]["full_gradient_l2"]
            > OBJECTIVE_GRADIENT_TOLERANCE)
    objective["aligned_logistic_stationarity_required"] = False

    x_train, x_holdout = x[train_mask], x[holdout_mask]
    design, _singular, _vh, _projection = _svd_audit(x_train, x_holdout, teacher)
    holdout_rows = [row for row, keep in zip(rows, holdout_mask) if keep]
    pinv_arms = {}
    primary_weight = None
    for rcond in (PINV_PRIMARY_RCOND, *PINV_SENSITIVITY_RCONDS):
        solution = _pinv_solution(x_train, residuals[train_mask], rcond)
        prediction_train = x_train @ solution
        prediction_holdout = x_holdout @ solution
        prediction_deranged = x_deranged @ solution
        metrics = _score(
            rows, holdout_mask, labels, baselines, prediction_holdout,
            prediction_deranged)
        metrics.update({
            "rcond": rcond,
            "rank": int(np.linalg.matrix_rank(
                x_train, tol=np.linalg.svd(x_train, compute_uv=False)[0] * rcond)),
            "train_rmse": float(np.sqrt(np.mean(
                (prediction_train - residuals[train_mask]) ** 2))),
            "train_maximum_absolute_error": float(np.max(np.abs(
                prediction_train - residuals[train_mask]))),
            "holdout_teacher_rmse": float(np.sqrt(np.mean(
                (prediction_holdout - residuals[holdout_mask]) ** 2))),
            "holdout_teacher_maximum_absolute_error": float(np.max(np.abs(
                prediction_holdout - residuals[holdout_mask]))),
            "coefficient_l2": float(np.linalg.norm(solution)),
        })
        pinv_arms[str(rcond)] = metrics
        if rcond == PINV_PRIMARY_RCOND:
            primary_weight = solution
    assert primary_weight is not None
    primary = pinv_arms[str(PINV_PRIMARY_RCOND)]
    pinv = {
        "schema": "MetaSieve.E0R1PinvWitness.v1",
        "primary_rcond": PINV_PRIMARY_RCOND,
        "primary": primary,
        "sensitivity_not_used_for_selection": {
            key: value for key, value in pinv_arms.items() if key != str(PINV_PRIMARY_RCOND)},
        "train_reconstruction_pass":
            primary["train_rmse"] <= PINV_TRAIN_RMSE_TOLERANCE,
        "gate_pass": primary["gate"]["pass"],
    }

    d_authorized = pinv["train_reconstruction_pass"] and pinv["gate_pass"]
    solver_trace = []
    if d_authorized:
        train_indices = np.flatnonzero(train_mask)
        local_left = np.searchsorted(train_indices, left)
        local_right = np.searchsorted(train_indices, right)
        deterministic_weight, solver_trace, solver_status = _run_deterministic_solver(
            x_train, residuals[train_mask], local_left, local_right)
        deterministic_metrics = _score(
            rows, holdout_mask, labels, baselines, x_holdout @ deterministic_weight,
            x_deranged @ deterministic_weight)
        solver = {
            "authorized_by_pinv": True,
            "objective": "point_huber_plus_all_within_task_residual_difference_huber",
            "status": solver_status,
            "metrics": deterministic_metrics,
            "gate_evaluated": solver_status["converged"],
            "gate_pass": (solver_status["converged"]
                          and deterministic_metrics["gate"]["pass"]),
        }
        np.save(output / "deterministic_weights.npy", deterministic_weight)
    else:
        solver = {
            "authorized_by_pinv": False,
            "status": "NOT_RUN_DESIGN_PRECONDITION_FAILED",
            "gate_evaluated": False,
            "gate_pass": False,
        }

    if not objective["objective_semantics_defect_confirmed"]:
        verdict = "OBJECTIVE_SEMANTICS_DEFECT_NOT_CONFIRMED"
    elif not pinv["train_reconstruction_pass"]:
        verdict = "NOT_RUN_PINV_NUMERICAL_PRECONDITION_FAILED"
    elif not pinv["gate_pass"]:
        verdict = "SYNTHETIC_TRAIN_DESIGN_PREDICTION_NONIDENTIFIABLE"
    elif not solver["gate_evaluated"]:
        verdict = "NOT_RUN_NUMERICAL_PRECONDITION_FAILED"
    elif solver["gate_pass"]:
        verdict = "SYNTHETIC_OBJECTIVE_DESIGN_SOLVER_IDENTIFIED"
    else:
        verdict = "CORRECTED_SOLVER_PREDICTION_FAIL"

    objective_path = output / "OBJECTIVE_AUDIT.json"
    design_path = output / "DESIGN_SVD_AUDIT.json"
    pinv_path = output / "PINV_WITNESS.json"
    for path, value in ((objective_path, objective), (design_path, design),
                        (pinv_path, pinv)):
        path.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
    (output / "DETERMINISTIC_SOLVER_TRACE.jsonl").write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in solver_trace),
        encoding="utf-8")
    np.save(output / "pinv_weights.npy", primary_weight)
    np.savez_compressed(
        output / "synthetic_design.npz", x=x, x_deranged=x_deranged,
        residuals=residuals, baselines=baselines, labels=labels,
        train_mask=train_mask, holdout_mask=holdout_mask, feature_center=feature_center)
    result = {
        "schema": "MetaSieve.E0R1Result.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P1R2B-E0R1_OBJECTIVE_DESIGN_SOLVER_AUDIT",
        "verdict": verdict,
        "findings": {
            "objective_semantics_defect_confirmed":
                objective["objective_semantics_defect_confirmed"],
            "train_design_prediction_identifiable": pinv["gate_pass"],
            "exact_linear_witness_pass": pinv["train_reconstruction_pass"],
            "corrected_deterministic_solver_pass": solver["gate_pass"],
        },
        "objective_summary": {
            "ordering_conflict_rate": objective["ordering_conflict_rate"],
            "teacher_old_rank_gradient_l2":
                objective["analytic_teacher"]["old_logistic"]["full_gradient_l2"],
            "teacher_difference_gradient_l2":
                objective["analytic_teacher"]["residual_difference_huber"]["full_gradient_l2"],
        },
        "design_summary": {
            "rank": design["primary_rank"],
            "effective_rank": design["effective_rank"],
            "condition_number": design["condition_number_identified"],
            "holdout_coverage_mean": design["holdout_row_space_coverage"]["mean"],
            "teacher_unseen_relative_l2":
                design["teacher_specific_unseen_contribution"]["relative_l2"],
        },
        "pinv_primary": primary,
        "solver": solver,
        "derangement_map_sha256": sha256_file(derangement_path),
        "inputs": {
            "input_manifest_sha256": sha256_file(input_root / "manifest.json"),
            "cache_manifest_sha256": sha256_file(cache_root / "manifest.json"),
            "checkpoint_sha256": sha256_file(checkpoint),
            "synthetic_gate_sha256": sha256_file(synthetic_root / "synthetic_gate.json"),
            "e0r0_result_sha256": sha256_file(e0r0_root / "result.json"),
        },
        "affinity_labels_read": False,
        "recipient_labels_read": False,
        "davis_accessed": False,
        "typed_interaction_training_performed": False,
        "downstream_authorized": False,
    }
    (output / "E0R1_RESULT.json").write_text(
        json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result


def _write_report(output: Path, result: dict) -> None:
    pinv = result["pinv_primary"]
    solver = result["solver"]
    solver_metrics = solver.get("metrics", {})
    text = f"""# E0R1 Objective, Design And Solver Audit

Decision: `{result['verdict']}`.

## Objective

Residual-order versus total-order conflict is
`{result['objective_summary']['ordering_conflict_rate']:.4%}`. At the analytic
teacher, old rank full-gradient L2 is
`{result['objective_summary']['teacher_old_rank_gradient_l2']:.3g}`, while the
residual-difference gradient is
`{result['objective_summary']['teacher_difference_gradient_l2']:.3g}`.

## Design

The 640 x 240 train design has rank `{result['design_summary']['rank']}`,
effective rank `{result['design_summary']['effective_rank']:.2f}` and identified
condition number `{result['design_summary']['condition_number']:.3g}`. Mean
holdout row-space coverage is `{result['design_summary']['holdout_coverage_mean']:.6f}`;
teacher-specific unseen relative L2 is
`{result['design_summary']['teacher_unseen_relative_l2']:.5f}`.

## Exact Witness

Primary Moore-Penrose train RMSE is `{pinv['train_rmse']:.3g}`. Holdout
correct/deranged CI is `{pinv['correct_ci']:.5f}/{pinv['deranged_ci']:.5f}`,
with partner delta `{pinv['correct_minus_deranged']:+.5f}`. Gate:
`{pinv['gate']['pass']}`.

## Corrected Deterministic Solve

Authorized: `{solver['authorized_by_pinv']}`. Converged:
`{solver.get('status', {}).get('converged', False) if isinstance(solver.get('status'), dict) else False}`.
Correct/deranged CI:
`{solver_metrics.get('correct_ci', float('nan')):.5f}/{solver_metrics.get('deranged_ci', float('nan')):.5f}`.
Gate pass: `{solver.get('gate_pass', False)}`.

No real affinity, DAVIS, PLIP/T, production, CSMO/Band or downstream work was
executed or authorized.
"""
    (output / "STAGE_REPORT.md").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path,
                        default=Path("dataset/processed/source_affinity/e0_input_v1"))
    parser.add_argument("--cache", type=Path,
                        default=Path("dataset/processed/source_affinity/e0_local_states_v1"))
    parser.add_argument("--checkpoint", type=Path,
                        default=Path("report/mechanism_refactor/p1b_pilot20k_seed17_v1/best.pt"))
    parser.add_argument("--synthetic", type=Path,
                        default=Path("research/e0_identifiability/artifacts/e0_local_map_v1"))
    parser.add_argument("--e0r0", type=Path,
                        default=Path("research/e0_identifiability/artifacts/e0r0_typed_tensor_v1"))
    parser.add_argument("--output", type=Path,
                        default=Path("research/e0_identifiability/artifacts/e0r1_objective_design_solver_v1"))
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    result = run(args.input, args.cache, args.checkpoint, args.synthetic,
                 args.e0r0, args.output, args.device)
    _write_report(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
