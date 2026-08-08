"""Run the synthetic-only E0R2 proposal and numerical-closure audit."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from research.e0_identifiability.proposal_triage_contract import (
    FULL_GRADIENT_TOLERANCE, HUBER_DELTA, OBJECTIVE_TOLERANCE, PRIMARY_RCOND,
    STAGE, TRAIN_RMSE_TOLERANCE,
)
from research.e0_identifiability.run_objective_design_solver import (
    _pair_indices, _score,
)
from research.e0_identifiability.run_synthetic_pregate import _select_rows
from scripts.structure_sources.rcsb import sha256_file


def _difference_rows(features: np.ndarray, targets: np.ndarray,
                     left: np.ndarray, right: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return features[left] - features[right], targets[left] - targets[right]


def _augmented_least_squares(features: np.ndarray, targets: np.ndarray,
                             left: np.ndarray, right: np.ndarray,
                             rcond: float = PRIMARY_RCOND) -> np.ndarray:
    differences, target_differences = _difference_rows(features, targets, left, right)
    augmented_x = np.vstack((
        features / np.sqrt(len(features)),
        differences / np.sqrt(len(differences)),
    ))
    augmented_y = np.concatenate((
        targets / np.sqrt(len(targets)),
        target_differences / np.sqrt(len(target_differences)),
    ))
    return np.linalg.lstsq(augmented_x, augmented_y, rcond=rcond)[0]


def _huber(values: np.ndarray, delta: float = HUBER_DELTA) -> np.ndarray:
    absolute = np.abs(values)
    return np.where(absolute <= delta, 0.5 * values ** 2,
                    delta * (absolute - 0.5 * delta))


def _corrected_objective(features: np.ndarray, targets: np.ndarray,
                         weights: np.ndarray, left: np.ndarray,
                         right: np.ndarray) -> dict:
    prediction = features @ weights
    point_error = prediction - targets
    differences, target_differences = _difference_rows(
        features, targets, left, right)
    difference_error = differences @ weights - target_differences
    clipped_point = np.clip(point_error, -HUBER_DELTA, HUBER_DELTA)
    clipped_difference = np.clip(
        difference_error, -HUBER_DELTA, HUBER_DELTA)
    gradient = (features.T @ clipped_point / len(features)
                + differences.T @ clipped_difference / len(differences))
    point = float(_huber(point_error).mean())
    difference = float(_huber(difference_error).mean())
    return {
        "loss": point + difference,
        "point_huber": point,
        "difference_huber": difference,
        "full_gradient_l2": float(np.linalg.norm(gradient)),
        "train_rmse": float(np.sqrt(np.mean(point_error ** 2))),
        "train_maximum_absolute_error": float(np.max(np.abs(point_error))),
        "all_errors_in_quadratic_region": bool(
            np.max(np.abs(point_error)) <= HUBER_DELTA
            and np.max(np.abs(difference_error)) <= HUBER_DELTA),
    }


def _viewpoint_audit() -> dict:
    return {
        "schema": "MetaSieve.E0R2ViewpointAudit.v1",
        "supported_by_existing_evidence": {
            "p1b_geometry_is_partner_specific": True,
            "global_mif_readout_did_not_identify_affinity_increment": True,
            "old_e0_rank_objective_had_residual_total_semantics_mismatch": True,
            "cross_fitted_ligand_residual_is_required_for_future_source_test": True,
            "frozen_operator_should_not_be_changed_to_hide_missing_biology": True,
        },
        "not_yet_established": {
            "p1b_has_lost_real_affinity_direction": (
                "Real affinity direction has not been tested; synthetic geometry is sufficient."),
            "orientation_or_typed_channels_are_required": (
                "Plausible hypothesis requiring a held-out structural identifiability gate."),
            "reference_state_improves_affinity_direction": (
                "Requires a frozen marginal-matched decoy distribution and later affinity gate."),
            "few_shot_energy_adapter_is_identified": (
                "Blocked until a population mechanism potential passes a source gate."),
            "residue_slot_pooling_erases_typed_energetics": (
                "Plausible but no pooled-versus-explicit residue experiment exists."),
        },
        "semantic_corrections": {
            "statistical_potential": (
                "-log(p_bound/p_ref) is a reference-dependent structural log-density ratio, "
                "not binding free energy."),
            "potential_sign": (
                "The proposal's negative-log-ratio convention conflicts with its later "
                "positive-is-favorable convention and must be fixed before use."),
            "uncertainty": (
                "Taking an expected potential under p(distance) propagates an expectation; "
                "it does not preserve the full uncertainty law."),
            "orthogonality": (
                "Cross-fitting prevents in-sample nuisance leakage but does not by itself "
                "establish Neyman orthogonality or causal interpretation."),
            "adapter_identifiability": (
                "d <= k is necessary, not sufficient; support design rank and query "
                "row-space coverage are also required."),
        },
        "frozen_followup_order": [
            "synthetic_numerical_closure",
            "structural_directional_identifiability",
            "reference_state_structural_log_odds",
            "source_oof_delta_affinity",
            "few_shot_mechanism_adapter",
            "production_z_and_operator_integration",
        ],
    }


def run(design_path: Path, rows_path: Path, e0r1_result_path: Path,
        output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    design = np.load(design_path)
    x = design["x"].astype(np.float64)
    x_deranged = design["x_deranged"].astype(np.float64)
    residuals = design["residuals"].astype(np.float64)
    baselines = design["baselines"].astype(np.float64)
    labels = design["labels"].astype(np.float64)
    train_mask = design["train_mask"].astype(bool)
    holdout_mask = design["holdout_mask"].astype(bool)
    rows = _select_rows(rows_path)
    if len(rows) != len(x):
        raise ValueError("frozen rows and E0R1 design have different lengths")

    global_left, global_right = _pair_indices(rows, train_mask)
    train_indices = np.flatnonzero(train_mask)
    left = np.searchsorted(train_indices, global_left)
    right = np.searchsorted(train_indices, global_right)
    x_train = x[train_mask]
    residual_train = residuals[train_mask]
    weights = _augmented_least_squares(
        x_train, residual_train, left, right)
    objective = _corrected_objective(
        x_train, residual_train, weights, left, right)
    metrics = _score(
        rows, holdout_mask, labels, baselines,
        x[holdout_mask] @ weights, x_deranged @ weights)
    checks = {
        "train_rmse": objective["train_rmse"] <= TRAIN_RMSE_TOLERANCE,
        "objective": objective["loss"] <= OBJECTIVE_TOLERANCE,
        "full_gradient": objective["full_gradient_l2"] <= FULL_GRADIENT_TOLERANCE,
        "quadratic_region": objective["all_errors_in_quadratic_region"],
        "historical_prediction_gate": metrics["gate"]["pass"],
    }
    passed = all(checks.values())
    result = {
        "schema": "MetaSieve.E0R2ProposalNumericalClosure.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": STAGE,
        "verdict": ("SYNTHETIC_OBJECTIVE_DESIGN_SOLVER_IDENTIFIED" if passed
                    else "NOT_RUN_NUMERICAL_PRECONDITION_FAILED"),
        "solver": {
            "method": "float64_svd_augmented_least_squares",
            "rcond": PRIMARY_RCOND,
            "selected_on_holdout": False,
            "objective": "point_huber_plus_all_within_task_residual_difference_huber",
            "train_rows": int(train_mask.sum()),
            "train_pairs": int(len(left)),
            "coefficient_l2": float(np.linalg.norm(weights)),
            **objective,
        },
        "historical_holdout_diagnostic": metrics,
        "checks": {**checks, "pass": passed},
        "control_limitations": {
            "holdout_is_untouched": False,
            "reason": "The same eight holdout tasks were used by E0R0/E0R1 diagnostics.",
            "derangement_reuse_bias_inherited": True,
            "generalization_claim_authorized": False,
        },
        "inputs": {
            "design_sha256": sha256_file(design_path),
            "rows_sha256": sha256_file(rows_path),
            "e0r1_result_sha256": sha256_file(e0r1_result_path),
        },
        "real_affinity_labels_read": False,
        "recipient_labels_read": False,
        "davis_accessed": False,
        "typed_interaction_training_performed": False,
        "downstream_authorized": False,
        "production_integration_authorized": False,
    }
    np.save(output / "weights.npy", weights)
    (output / "VIEWPOINT_AUDIT.json").write_text(
        json.dumps(_viewpoint_audit(), indent=2, sort_keys=True), encoding="utf-8")
    (output / "NUMERICAL_CLOSURE.json").write_text(
        json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result


def _write_report(output: Path, result: dict) -> None:
    solver = result["solver"]
    metrics = result["historical_holdout_diagnostic"]
    text = f"""# E0R2 Proposal And Numerical Closure

Decision: `{result['verdict']}`.

The corrected synthetic objective was solved once in float64 with a frozen SVD
augmented least-squares solve. Train RMSE is `{solver['train_rmse']:.3g}`, full
gradient L2 is `{solver['full_gradient_l2']:.3g}`, and corrected objective is
`{solver['loss']:.3g}`.

On the historical eight-task development diagnostic, correct/deranged CI is
`{metrics['correct_ci']:.5f}/{metrics['deranged_ci']:.5f}`; correct-minus-ligand
is `{metrics['correct_minus_ligand']:+.5f}` and correct-minus-deranged is
`{metrics['correct_minus_deranged']:+.5f}`.

This closes the synthetic objective/design/solver boundary only. The holdout is
not untouched and the inherited derangement control reuses wrong proteins. No
real affinity, directionality, reference-state potential, few-shot adapter, or
production biological statistic was tested. No code is authorized for promotion
to `model/` or normal `scripts/`.
"""
    (output / "STAGE_REPORT.md").write_text(text, encoding="utf-8")


def _write_manifest(output: Path, result: dict) -> None:
    root = Path(__file__).resolve().parent
    artifact_names = (
        "NUMERICAL_CLOSURE.json", "STAGE_REPORT.md", "VIEWPOINT_AUDIT.json",
        "weights.npy",
    )
    manifest = {
        "schema": "MetaSieve.E0R2ArtifactManifest.v1",
        "stage": STAGE,
        "verdict": result["verdict"],
        "inputs": result["inputs"],
        "code": {
            "runner_sha256": sha256_file(Path(__file__)),
            "contract_sha256": sha256_file(root / "proposal_triage_contract.py"),
            "preregistration_sha256": sha256_file(root / "E0R2_PREREGISTRATION.md"),
        },
        "artifacts": {
            name: sha256_file(output / name) for name in artifact_names
        },
    }
    (output / "ARTIFACT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--design", type=Path,
        default=Path("research/e0_identifiability/artifacts/"
                     "e0r1_objective_design_solver_v1/synthetic_design.npz"))
    parser.add_argument(
        "--rows", type=Path,
        default=Path("dataset/processed/source_affinity/e0_input_v1/"
                     "rows.label_blind.jsonl"))
    parser.add_argument(
        "--e0r1-result", type=Path,
        default=Path("research/e0_identifiability/artifacts/"
                     "e0r1_objective_design_solver_v1/E0R1_RESULT.json"))
    parser.add_argument(
        "--output", type=Path,
        default=Path("research/e0_identifiability/artifacts/"
                     "e0r2_proposal_numerical_closure_v1"))
    args = parser.parse_args()
    result = run(args.design, args.rows, args.e0r1_result, args.output)
    _write_report(args.output, result)
    _write_manifest(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
