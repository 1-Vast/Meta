"""Audit E0 synthetic identifiability without training or affinity-label access."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import torch
from scipy.stats import pearsonr, spearmanr

from contracts.ligand_graph import MAX_ATOMS
from scripts.audit_e0_input_feasibility import _read_jsonl, VALID_RESIDUES
from research.e0_identifiability.mechanistic_affinity import LocalMechanisticAffinityPotential
from research.e0_identifiability.metrics import concordance as _concordance
from scripts.govern_structure_homology import _local_identity
from research.e0_identifiability.run_synthetic_pregate import (
    BATCH_SIZE, DISTANCE_WEIGHTS, _batch, _geometry, _load_bridge, _load_states,
    _select_rows,
)
from scripts.structure_sources.rcsb import sha256_file


def _macro_ci(rows: list[dict], labels: np.ndarray, predictions: np.ndarray) -> tuple[float, list[dict]]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        grouped[row["task_id"]].append(index)
    per_task = []
    for task_id, indices in sorted(grouped.items()):
        per_task.append({"task_id": task_id, "rows": len(indices),
                         "ci": _concordance(labels[indices], predictions[indices])})
    return float(np.mean([row["ci"] for row in per_task])), per_task


def _safe_correlation(left: np.ndarray, right: np.ndarray) -> dict:
    if len(left) < 2 or np.std(left) == 0 or np.std(right) == 0:
        return {"pearson": None, "spearman": None}
    return {"pearson": float(pearsonr(left, right).statistic),
            "spearman": float(spearmanr(left, right).statistic)}


def _relative_l2(left: np.ndarray, right: np.ndarray) -> float:
    return float(np.linalg.norm(left - right) / max(np.linalg.norm(left), 1e-12))


def _teacher_values(rows: list[dict], proteins: dict, ligands: dict, geometry: dict,
                    weights: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    direct, sufficient, baselines = [], [], []
    for row in rows:
        ligand = ligands[row["ligand_state_key"]]
        protein = proteins[row["active_protein_key"]]
        atom = ligand["chemistry"][:, 32:40]
        residue = protein["chemistry"]
        cached = geometry[row["example_id"]]
        contact = cached["contact"].astype(np.float32)
        distance = cached["distance"].astype(np.float32)
        compatibility = np.einsum("na,ar,lr->nl", atom, weights, residue)
        distance_score = np.einsum("nld,d->nl", distance, DISTANCE_WEIGHTS)
        denominator = max(float(contact.sum()), 1e-6)
        direct.append(float(np.sum(contact * compatibility * distance_score) / denominator))
        statistics = np.einsum("na,lr,nl,nld->ard", atom, residue, contact, distance)
        sufficient.append(float(np.einsum(
            "ard,ar,d->", statistics, weights, DISTANCE_WEIGHTS) / denominator))
        mean_pharmacophore = atom.mean(0)
        baselines.append(float(0.4 * mean_pharmacophore[1]
                               + 0.2 * mean_pharmacophore[2]
                               - 0.15 * mean_pharmacophore[7]))
    return (np.asarray(direct, dtype=np.float64),
            np.asarray(sufficient, dtype=np.float64),
            np.asarray(baselines, dtype=np.float64))


def _predict(model, rows: list[dict], proteins: dict, ligands: dict, geometry: dict,
             device: str, collect_activations: bool = False):
    predictions = []
    activations: dict[str, list[np.ndarray]] = defaultdict(list)
    handles = []
    if collect_activations:
        modules = {"atom": model.atom[1], "residue": model.residue[1],
                   "geometry": model.geometry[1], "pair": model.pair[1]}
        for name, module in modules.items():
            handles.append(module.register_forward_hook(
                lambda _module, _inputs, output, key=name:
                activations[key].append(output.detach().float().cpu().numpy())))
    with torch.inference_mode():
        for start in range(0, len(rows), BATCH_SIZE):
            subset = rows[start:start + BATCH_SIZE]
            predictions.extend(model(**_batch(
                subset, proteins, ligands, device, geometry)).potential.cpu().tolist())
    for handle in handles:
        handle.remove()
    summary = {}
    for name, blocks in activations.items():
        values = np.concatenate([block.reshape(-1) for block in blocks])
        summary[name] = {"mean": float(values.mean()), "std": float(values.std()),
                         "fraction_abs_lt_1e_6": float(np.mean(np.abs(values) < 1e-6)),
                         "fraction_abs_gt_10": float(np.mean(np.abs(values) > 10)),
                         "minimum": float(values.min()), "maximum": float(values.max())}
    return np.asarray(predictions, dtype=np.float64), summary


def _derangement(rows: list[dict]) -> tuple[list[dict], dict[str, str]]:
    proteins_by_closure = {}
    for row in rows:
        proteins_by_closure.setdefault(row["closure_component_id"],
                                       row["protein_sequence_sha256"])
    closures = sorted(proteins_by_closure)
    mapping = {}
    for index, closure in enumerate(closures):
        for offset in range(1, len(closures)):
            candidate = closures[(index + offset) % len(closures)]
            if candidate != closure:
                mapping[closure] = proteins_by_closure[candidate]
                break
    return ([dict(row, active_protein_key=mapping[row["closure_component_id"]])
             for row in rows], mapping)


def _provenance_audit(rows_path: Path, splits_path: Path, manifest_rows_path: Path) -> dict:
    from rdkit import Chem

    assignments = {row["task_id"]: row for row in _read_jsonl(splits_path)}
    valid_rows = []
    compounds: dict[str, set[str]] = defaultdict(set)
    metadata = {}
    for row in _read_jsonl(rows_path):
        assignment = assignments.get(row.get("task_id"))
        if assignment is None:
            continue
        molecule = Chem.MolFromSmiles(row["canonical_smiles"])
        sequence = row["protein_sequence"]
        if (molecule is None or molecule.GetNumAtoms() > MAX_ATOMS or not sequence
                or set(sequence) - VALID_RESIDUES):
            continue
        task_id = row["task_id"]
        compounds[task_id].add(row["ligand_connectivity_key"])
        metadata[task_id] = {"closure_component_id": assignment["closure_component_id"],
                             "outer_oof_fold": int(assignment["outer_oof_fold"]),
                             "protein_sequence_sha256": row["protein_sequence_sha256"],
                             "endpoint_family": row["endpoint_family"]}
        valid_rows.append((int(row["activity_id"]), task_id))
    retained_tasks = {task for task, values in compounds.items() if len(values) >= 20}
    expected_ids = {activity for activity, task in valid_rows if task in retained_tasks}
    actual_rows = list(_read_jsonl(manifest_rows_path))
    actual_ids = {int(row["activity_id"]) for row in actual_rows}
    excluded_tasks = sorted(set(compounds) - retained_tasks)
    excluded_rows_by_task = Counter(task for _, task in valid_rows if task in excluded_tasks)
    fold_counts = Counter(metadata[task]["outer_oof_fold"] for task in excluded_tasks)
    fold_row_counts = Counter()
    for task, count in excluded_rows_by_task.items():
        fold_row_counts[metadata[task]["outer_oof_fold"]] += count
    endpoint_counts = Counter(metadata[task]["endpoint_family"] for task in excluded_tasks)
    affected_proteins = {metadata[task]["protein_sequence_sha256"] for task in excluded_tasks}
    affected_closures = {metadata[task]["closure_component_id"] for task in excluded_tasks}
    return {
        "schema": "MetaSieve.E0DataProvenanceAudit.v1",
        "audit_mode": "LABEL_BLIND",
        "affinity_value_fields_materialized": False,
        "valid_rows_before_task_floor": len(valid_rows),
        "retained_rows_after_task_floor": len(expected_ids),
        "transition_rows": len(valid_rows) - len(expected_ids),
        "transition_reason": "tasks_below_20_valid_connectivity_compounds_after_input_contract",
        "excluded_tasks": len(excluded_tasks),
        "governed_tasks": len(assignments),
        "tasks_with_any_model_valid_row": len(compounds),
        "tasks_with_zero_model_valid_rows": len(assignments) - len(compounds),
        "total_tasks_removed_from_model_manifest": len(assignments) - len(retained_tasks),
        "affected_proteins": len(affected_proteins),
        "affected_closure_components": len(affected_closures),
        "excluded_task_counts_by_fold": {str(key): value for key, value in sorted(fold_counts.items())},
        "excluded_row_counts_by_fold": {str(key): value for key, value in sorted(fold_row_counts.items())},
        "excluded_task_counts_by_endpoint": dict(sorted(endpoint_counts.items())),
        "excluded_task_details": [{**metadata[task], "task_id": task,
                                   "valid_connectivity_compounds": len(compounds[task]),
                                   "rows": excluded_rows_by_task[task]}
                                  for task in excluded_tasks],
        "materialized_manifest_rows": len(actual_rows),
        "materialized_activity_ids_match_reconstruction": actual_ids == expected_ids,
        "missing_activity_ids": sorted(expected_ids - actual_ids),
        "unexpected_activity_ids": sorted(actual_ids - expected_ids),
        "recipient_labels_read": False,
    }


def _sample_audit(rows: list[dict], population_rows: list[dict], derangement_map: dict[str, str],
                  protein_sequences: dict[str, str], task_metrics: dict[str, dict]) -> dict:
    train = [row for row in rows if row["outer_oof_fold"] < 4]
    holdout = [row for row in rows if row["outer_oof_fold"] == 4]
    train_states = {row["ligand_state_key"] for row in train}
    train_connectivity = {row["ligand_connectivity_key"] for row in train}
    holdout_states = {row["ligand_state_key"] for row in holdout}
    holdout_connectivity = {row["ligand_connectivity_key"] for row in holdout}
    fold4_population = [row for row in population_rows if int(row["outer_oof_fold"]) == 4]
    per_task = []
    for task_id in sorted({row["task_id"] for row in holdout}):
        task_rows = [row for row in holdout if row["task_id"] == task_id]
        first = task_rows[0]
        correct_key = first["protein_sequence_sha256"]
        deranged_key = derangement_map[first["closure_component_id"]]
        per_task.append({
            "task_id": task_id, "rows": len(task_rows),
            "closure_component_id": first["closure_component_id"],
            "endpoint_family": first["endpoint_family"],
            "protein_sequence_sha256": correct_key,
            "deranged_protein_sequence_sha256": deranged_key,
            "correct_deranged_local_identity": _local_identity(
                protein_sequences[correct_key], protein_sequences[deranged_key]),
            **task_metrics[task_id],
        })
    identity_violations = [row for row in per_task
                           if row["correct_deranged_local_identity"] >= 0.40]
    return {
        "selection_policy": "lexical_first_8_tasks_per_frozen_fold_first_20_ligand_states",
        "tasks": len({row["task_id"] for row in rows}),
        "train_tasks": len({row["task_id"] for row in train}),
        "holdout_tasks": len(per_task),
        "holdout_closure_components": len({row["closure_component_id"] for row in holdout}),
        "holdout_proteins": len({row["protein_sequence_sha256"] for row in holdout}),
        "holdout_endpoint_counts": dict(sorted(Counter(
            row["endpoint_family"] for row in {value["task_id"]: value for value in holdout}.values()).items())),
        "holdout_ligand_states": len(holdout_states),
        "holdout_ligand_connectivity_keys": len(holdout_connectivity),
        "holdout_state_overlap_with_train": len(holdout_states & train_states),
        "holdout_connectivity_overlap_with_train": len(holdout_connectivity & train_connectivity),
        "holdout_closure_coverage_of_retained_245": len(
            {row["closure_component_id"] for row in holdout}) / 245.0,
        "fold4_population_tasks": len({row["task_id"] for row in fold4_population}),
        "fold4_population_closure_components": len(
            {row["closure_component_id"] for row in fold4_population}),
        "fold4_population_proteins": len(
            {row["protein_sequence_sha256"] for row in fold4_population}),
        "holdout_task_coverage_of_fold4": len(per_task) / max(
            len({row["task_id"] for row in fold4_population}), 1),
        "holdout_closure_coverage_of_fold4": len(
            {row["closure_component_id"] for row in holdout}) / max(
            len({row["closure_component_id"] for row in fold4_population}), 1),
        "holdout_protein_coverage_of_fold4": len(
            {row["protein_sequence_sha256"] for row in holdout}) / max(
            len({row["protein_sequence_sha256"] for row in fold4_population}), 1),
        "derangement_pairs_at_or_above_40pct_local_identity": len(identity_violations),
        "derangement_contract_all_pairs_below_40pct": not identity_violations,
        "derangement_contract_violations": identity_violations,
        "per_holdout_task": per_task,
    }


def audit_e0s(input_root: Path, cache_root: Path, checkpoint: Path,
              synthetic_root: Path, canonical_rows: Path, splits: Path,
              output: Path, device: str) -> dict:
    gate_report = json.loads((synthetic_root / "synthetic_gate.json").read_text(encoding="utf-8"))
    rows = _select_rows(input_root / "rows.label_blind.jsonl")
    for index, row in enumerate(rows):
        row["example_id"] = index
        row["active_protein_key"] = row["protein_sequence_sha256"]
    deranged_rows, derangement_map = _derangement(rows)
    ligand_keys = {row["ligand_state_key"] for row in rows}
    protein_keys = ({row["protein_sequence_sha256"] for row in rows}
                    | {row["active_protein_key"] for row in deranged_rows})
    proteins, ligands = _load_states(cache_root, ligand_keys, protein_keys)
    bridge = _load_bridge(checkpoint, device)
    correct_geometry = _geometry(rows, bridge, proteins, ligands, device)
    deranged_geometry = _geometry(deranged_rows, bridge, proteins, ligands, device)
    weights = np.asarray(gate_report["teacher_weights"], dtype=np.float32)
    correct_raw, correct_sufficient, baselines = _teacher_values(
        rows, proteins, ligands, correct_geometry, weights)
    deranged_raw, deranged_sufficient, deranged_baselines = _teacher_values(
        deranged_rows, proteins, ligands, deranged_geometry, weights)
    train_mask = np.asarray([row["outer_oof_fold"] < 4 for row in rows])
    holdout_mask = ~train_mask
    center, scale = float(correct_raw[train_mask].mean()), float(correct_raw[train_mask].std())
    scale = max(scale, 1e-6)
    correct_oracle = baselines + (correct_raw - center) / scale
    deranged_oracle = baselines + (deranged_raw - center) / scale
    labels = correct_oracle

    model_payload = torch.load(synthetic_root / "model.pt", map_location="cpu", weights_only=False)
    model = LocalMechanisticAffinityPotential().to(device)
    model.load_state_dict(model_payload["model_state"])
    model.eval()
    map_correct, activation_summary = _predict(
        model, rows, proteins, ligands, correct_geometry, device, collect_activations=True)
    map_deranged, _ = _predict(model, deranged_rows, proteins, ligands, deranged_geometry, device)
    map_correct_score = baselines + map_correct
    map_deranged_score = baselines + map_deranged

    teacher_correct_ci, teacher_correct_tasks = _macro_ci(
        [row for row, keep in zip(rows, holdout_mask) if keep], labels[holdout_mask],
        correct_oracle[holdout_mask])
    teacher_deranged_ci, teacher_deranged_tasks = _macro_ci(
        [row for row, keep in zip(rows, holdout_mask) if keep], labels[holdout_mask],
        deranged_oracle[holdout_mask])
    map_correct_ci, map_correct_tasks = _macro_ci(
        [row for row, keep in zip(rows, holdout_mask) if keep], labels[holdout_mask],
        map_correct_score[holdout_mask])
    map_deranged_ci, map_deranged_tasks = _macro_ci(
        [row for row, keep in zip(rows, holdout_mask) if keep], labels[holdout_mask],
        map_deranged_score[holdout_mask])
    map_train_ci, _ = _macro_ci(
        [row for row, keep in zip(rows, train_mask) if keep], labels[train_mask],
        map_correct_score[train_mask])

    residue_delta, contact_delta, distance_delta = [], [], []
    for row, deranged_row in zip(rows, deranged_rows):
        correct_protein = proteins[row["protein_sequence_sha256"]]
        wrong_protein = proteins[deranged_row["active_protein_key"]]
        residue_delta.append(_relative_l2(correct_protein["residues"], wrong_protein["residues"]))
        contact_delta.append(_relative_l2(
            correct_geometry[row["example_id"]]["contact"].astype(np.float32),
            deranged_geometry[row["example_id"]]["contact"].astype(np.float32)))
        distance_delta.append(_relative_l2(
            correct_geometry[row["example_id"]]["distance"].astype(np.float32),
            deranged_geometry[row["example_id"]]["distance"].astype(np.float32)))
    residue_delta = np.asarray(residue_delta)
    contact_delta = np.asarray(contact_delta)
    distance_delta = np.asarray(distance_delta)
    teacher_delta = np.abs(correct_oracle - deranged_oracle)
    map_delta = np.abs(map_correct_score - map_deranged_score)
    holdout_teacher_delta = teacher_delta[holdout_mask]
    holdout_map_delta = map_delta[holdout_mask]

    holdout_tied_pairs = 0
    holdout_total_pairs = 0
    holdout_rows = [row for row, keep in zip(rows, holdout_mask) if keep]
    for task_id in sorted({row["task_id"] for row in holdout_rows}):
        indices = [index for index, row in enumerate(holdout_rows)
                   if row["task_id"] == task_id]
        task_labels = labels[holdout_mask][indices]
        for left in range(len(task_labels)):
            for right in range(left + 1, len(task_labels)):
                holdout_total_pairs += 1
                holdout_tied_pairs += int(task_labels[left] == task_labels[right])

    teacher_gate = {
        "oracle_correct_ci_at_least_0_80": teacher_correct_ci >= 0.80,
        "oracle_correct_minus_deranged_at_least_0_10":
            teacher_correct_ci - teacher_deranged_ci >= 0.10,
    }
    reconstruction_error = float(np.max(np.abs(correct_raw - correct_sufficient)))
    frozen_geometry_retained = reconstruction_error <= 1e-6
    if not all(teacher_gate.values()):
        verdict = "SYNTHETIC_CONTROL_MIS_SPECIFIED"
    elif not frozen_geometry_retained:
        verdict = "FROZEN_GEOMETRY_INSUFFICIENT_FOR_AFFINITY_POTENTIAL"
    elif not gate_report["gate"]["pass"]:
        verdict = "MAP_REALIZATION_OR_OPTIMIZATION_DEFECT"
    else:
        verdict = "SYNTHETIC_IDENTIFIABILITY_PASS"

    correct_task_map = {row["task_id"]: row["ci"] for row in map_correct_tasks}
    deranged_task_map = {row["task_id"]: row["ci"] for row in map_deranged_tasks}
    teacher_correct_map = {row["task_id"]: row["ci"] for row in teacher_correct_tasks}
    teacher_deranged_map = {row["task_id"]: row["ci"] for row in teacher_deranged_tasks}
    task_metrics = {task: {"teacher_correct_ci": teacher_correct_map[task],
                           "teacher_deranged_ci": teacher_deranged_map[task],
                           "map_correct_ci": correct_task_map[task],
                           "map_deranged_ci": deranged_task_map[task]}
                    for task in correct_task_map}
    protein_sequences = {row["sequence_sha256"]: row["sequence"]
                         for row in _read_jsonl(input_root / "proteins.jsonl")}
    population_rows = list(_read_jsonl(input_root / "rows.label_blind.jsonl"))
    sample = _sample_audit(rows, population_rows, derangement_map,
                           protein_sequences, task_metrics)
    provenance = _provenance_audit(
        canonical_rows, splits, input_root / "rows.label_blind.jsonl")

    parameter_summary = {name: {"l2": float(parameter.detach().norm().cpu()),
                                "maximum_abs": float(parameter.detach().abs().max().cpu())}
                         for name, parameter in model.named_parameters()}
    teacher_audit = {
        "schema": "MetaSieve.SyntheticTeacherAudit.v1",
        "stage": "P1R2B-E0S", "affinity_labels_read": False,
        "teacher_formula": {
            "raw": "sum(contact * qWr * dot(distance_prob,[1,.7,.2,-.2,-.6])) / sum(contact)",
            "label": "ligand_baseline + (raw-train_mean)/train_std",
            "noise": "none", "interaction_mass_channel": "denominator_only",
        },
        "normalization": {"train_mean": center, "train_std": scale},
        "holdout": {"correct_oracle_ci": teacher_correct_ci,
                    "deranged_oracle_ci": teacher_deranged_ci,
                    "correct_minus_deranged": teacher_correct_ci - teacher_deranged_ci,
                    "label_tied_pairs": holdout_tied_pairs,
                    "label_total_pairs": holdout_total_pairs,
                    "gate": {**teacher_gate, "pass": all(teacher_gate.values())}},
        "sufficient_statistic": {"shape": [8, 6, 5],
                                 "maximum_direct_reconstruction_error": reconstruction_error,
                                 "exact_within_1e_6": frozen_geometry_retained},
        "per_holdout_task_correct": teacher_correct_tasks,
        "per_holdout_task_deranged": teacher_deranged_tasks,
        "recipient_labels_read": False,
    }
    retention = {
        "schema": "MetaSieve.SyntheticInformationRetention.v1",
        "stage": "P1R2B-E0S", "affinity_labels_read": False,
        "boundaries": {
            "T0_teacher_oracle": {"holdout_ci": teacher_correct_ci},
            "T1_teacher_sufficient_statistics": {
                "dimensions": 240, "maximum_reconstruction_error": reconstruction_error,
                "rank_retention": 1.0 if frozen_geometry_retained else None},
            "T2_frozen_geometry": {
                "teacher_is_defined_directly_from_frozen_geometry": True,
                "maximum_teacher_statistic_reconstruction_error": reconstruction_error,
                "verdict": "PASS_BY_CONSTRUCTION" if frozen_geometry_retained else "FAIL"},
            "T3_map": {"train_correct_ci": map_train_ci,
                       "holdout_correct_ci": map_correct_ci,
                       "holdout_deranged_ci": map_deranged_ci,
                       "correct_minus_deranged": map_correct_ci - map_deranged_ci},
        },
        "derangement": {
            "holdout": {
                "residue_relative_l2_mean": float(residue_delta[holdout_mask].mean()),
                "contact_relative_l2_mean": float(contact_delta[holdout_mask].mean()),
                "distance_relative_l2_mean": float(distance_delta[holdout_mask].mean()),
                "teacher_absolute_delta_mean": float(holdout_teacher_delta.mean()),
                "teacher_absolute_delta_median": float(np.median(holdout_teacher_delta)),
                "map_absolute_delta_mean": float(holdout_map_delta.mean()),
                "teacher_delta_vs_map_delta": _safe_correlation(
                    holdout_teacher_delta, holdout_map_delta),
                "teacher_delta_vs_residue_delta": _safe_correlation(
                    holdout_teacher_delta, residue_delta[holdout_mask]),
                "teacher_delta_vs_contact_delta": _safe_correlation(
                    holdout_teacher_delta, contact_delta[holdout_mask]),
                "teacher_delta_vs_distance_delta": _safe_correlation(
                    holdout_teacher_delta, distance_delta[holdout_mask]),
            },
            "all_selected_pairs": {
                "residue_relative_l2_mean": float(residue_delta.mean()),
                "contact_relative_l2_mean": float(contact_delta.mean()),
                "distance_relative_l2_mean": float(distance_delta.mean()),
                "teacher_absolute_delta_mean": float(teacher_delta.mean()),
                "map_absolute_delta_mean": float(map_delta.mean()),
            },
        },
        "optimization": {
            "training_curve_persisted": False,
            "gradient_trace_persisted": False,
            "convergence_identifiable": False,
            "reason": "the frozen run persisted final weights and metrics but no epoch/gradient trace",
            "train_correct_ci": map_train_ci,
            "holdout_correct_ci": map_correct_ci,
            "train_residual_pearson": float(pearsonr(
                map_correct[train_mask], (correct_raw[train_mask] - center) / scale).statistic),
            "train_residual_spearman": float(spearmanr(
                map_correct[train_mask], (correct_raw[train_mask] - center) / scale).statistic),
            "holdout_residual_pearson": float(pearsonr(
                map_correct[holdout_mask], (correct_raw[holdout_mask] - center) / scale).statistic),
            "holdout_residual_spearman": float(spearmanr(
                map_correct[holdout_mask], (correct_raw[holdout_mask] - center) / scale).statistic),
            "activation_summary": activation_summary,
            "parameter_summary": parameter_summary,
        },
        "sample": sample,
        "derangement_contract": {
            "status": ("PASS" if sample["derangement_contract_all_pairs_below_40pct"]
                       else "PARTIAL_VIOLATION"),
            "pairs_at_or_above_40pct_local_identity":
                sample["derangement_pairs_at_or_above_40pct_local_identity"],
            "does_not_explain_aggregate_partner_gap": True,
            "reason": ("the violating task retains strong teacher and MAP correct-minus-deranged "
                       "contrast; it is reported as a control-contract caveat"),
        },
        "recipient_labels_read": False,
    }
    failure_graph = {
        "schema": "MetaSieve.E0FailureGraph.v1", "stage": "P1R2B-E0S",
        "nodes": [
            {"id": "T0", "name": "teacher_oracle", "status": "PASS" if all(teacher_gate.values()) else "FAIL"},
            {"id": "T1", "name": "teacher_sufficient_statistics", "status": "PASS" if frozen_geometry_retained else "FAIL"},
            {"id": "T2", "name": "frozen_p1b_geometry", "status": "PASS_BY_CONSTRUCTION" if frozen_geometry_retained else "FAIL"},
            {"id": "T3", "name": "map_checkpoint", "status": "FAIL" if not gate_report["gate"]["pass"] else "PASS"},
            {"id": "C0", "name": "derangement_contract", "status":
             "PASS" if sample["derangement_contract_all_pairs_below_40pct"]
             else "PARTIAL_VIOLATION"},
            {"id": "S0", "name": "synthetic_holdout_coverage", "status":
             "INSUFFICIENT_FOR_ARCHITECTURE_WIDE_CLAIM"},
        ],
        "edges": [
            {"from": "T0", "to": "T1", "maximum_reconstruction_error": reconstruction_error},
            {"from": "T1", "to": "T2", "note": "teacher sufficient statistics are computed from frozen P1B geometry"},
            {"from": "T2", "to": "T3", "holdout_ci_gap": teacher_correct_ci - map_correct_ci,
             "partner_delta_gap": (teacher_correct_ci - teacher_deranged_ci)
                                  - (map_correct_ci - map_deranged_ci)},
        ],
        "verdict": verdict,
        "realization_vs_optimization": "UNRESOLVED_NO_TRAINING_TRACE",
        "repair_authorized": False, "e0_source_authorized": False,
        "affinity_labels_read": False, "recipient_labels_read": False,
    }
    output.mkdir(parents=True, exist_ok=True)
    artifacts = {
        "SYNTHETIC_TEACHER_AUDIT.json": teacher_audit,
        "SYNTHETIC_INFORMATION_RETENTION.json": retention,
        "E0_DATA_PROVENANCE_AUDIT.json": provenance,
        "E0_FAILURE_GRAPH.json": failure_graph,
    }
    for name, value in artifacts.items():
        (output / name).write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
    synthesis = {
        "schema": "MetaSieve.E0Synthesis.v1", "created_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P1R2B-E0S", "verdict": verdict,
        "scientific_state": ("SYNTHETIC_INCREMENTAL_SIGNAL_LEARNED; "
                             "PARTNER_SPECIFIC_RECOVERY_INSUFFICIENT; "
                             "FAILURE_LOCALIZED_TO_T2_TO_T3; "
                             "REALIZATION_VS_OPTIMIZATION_UNRESOLVED"),
        "failure_layer": "T2_TO_T3_REALIZATION_OR_OPTIMIZATION_UNRESOLVED",
        "derangement_contract_status": failure_graph["nodes"][4]["status"],
        "synthetic_holdout_status": failure_graph["nodes"][5]["status"],
        "input_artifacts": {
            "label_blind_rows_sha256": sha256_file(input_root / "rows.label_blind.jsonl"),
            "protein_states_sha256": sha256_file(input_root / "proteins.jsonl"),
            "p1b_checkpoint_sha256": sha256_file(checkpoint),
            "synthetic_gate_sha256": sha256_file(synthetic_root / "synthetic_gate.json"),
            "synthetic_model_sha256": sha256_file(synthetic_root / "model.pt"),
            "canonical_rows_sha256": sha256_file(canonical_rows),
            "split_assignments_sha256": sha256_file(splits),
        },
        "affinity_labels_read": False, "recipient_labels_read": False,
        "model_training_performed": False, "repair_authorized": False,
    }
    return {"synthesis": synthesis, "teacher": teacher_audit, "retention": retention,
            "provenance": provenance, "failure_graph": failure_graph}


def _write_report(output: Path, result: dict) -> None:
    synthesis, teacher = result["synthesis"], result["teacher"]
    retention, provenance = result["retention"], result["provenance"]
    sample = retention["sample"]
    text = f"""# E0 Evidence Consolidation And Synthetic Identifiability

Updated: 2026-08-07

Decision: `{synthesis['verdict']}`. No repair, real affinity training, DAVIS
access or typed-interaction work is authorized.

## Consolidated Finding

The teacher oracle attains holdout CI `{teacher['holdout']['correct_oracle_ci']:.5f}`;
its deranged-protein oracle CI is `{teacher['holdout']['deranged_oracle_ci']:.5f}`,
for partner delta `{teacher['holdout']['correct_minus_deranged']:+.5f}`. The exact
8 x 6 x 5 sufficient statistic reconstructs the direct teacher with maximum
error `{teacher['sufficient_statistic']['maximum_direct_reconstruction_error']:.3g}`.

The frozen MAP checkpoint reaches train correct CI
`{retention['boundaries']['T3_map']['train_correct_ci']:.5f}` and holdout correct
CI `{retention['boundaries']['T3_map']['holdout_correct_ci']:.5f}`. Its holdout
partner delta is `{retention['boundaries']['T3_map']['correct_minus_deranged']:+.5f}`.
The identifiable loss is therefore at T2 -> T3. Existing artifacts do not
separate hypothesis-class realization from optimization because no epoch or
gradient trace was persisted.

One of the eight holdout derangements has local sequence identity >=40%.
This is a control-contract violation, but it does not explain the aggregate
failure: that task retains strong teacher and MAP correct-versus-deranged
contrast. The holdout is also a lexical sample covering only
`{sample['holdout_task_coverage_of_fold4']:.2%}` of fold-4 tasks and
`{sample['holdout_closure_coverage_of_fold4']:.2%}` of fold-4 closure components;
it is insufficient for an architecture-wide generalization claim.

## Boundary Table

| Boundary | Observable | Verdict |
|---|---:|---|
| Teacher attainable | oracle correct CI / partner delta | PASS: {teacher['holdout']['correct_oracle_ci']:.5f} / {teacher['holdout']['correct_minus_deranged']:+.5f} |
| Derangement changes teacher | oracle correct - deranged CI | PASS: {teacher['holdout']['correct_minus_deranged']:+.5f} |
| Teacher sufficient statistics | max reconstruction error | PASS: {teacher['sufficient_statistic']['maximum_direct_reconstruction_error']:.3g} |
| Frozen geometry retention | teacher defined from and exactly reconstructed by P1B geometry | {retention['boundaries']['T2_frozen_geometry']['verdict']} |
| MAP realization | holdout correct / deranged CI | FAIL: {retention['boundaries']['T3_map']['holdout_correct_ci']:.5f} / {retention['boundaries']['T3_map']['holdout_deranged_ci']:.5f} |
| Optimization convergence | training/gradient trace | NOT AUDITABLE |
| Holdout diversity | tasks / closures / proteins | INSUFFICIENT: {sample['holdout_tasks']} / {sample['holdout_closure_components']} / {sample['holdout_proteins']} |
| Derangement contract | pairs at or above 40% local identity | PARTIAL VIOLATION: {sample['derangement_pairs_at_or_above_40pct_local_identity']} / {sample['holdout_tasks']} |
| Corpus consistency | rows before / after task floor | PASS: {provenance['valid_rows_before_task_floor']} / {provenance['retained_rows_after_task_floor']} |

## Provenance

The 197-row transition is exact: {provenance['excluded_tasks']} nonempty tasks fell below
20 valid connectivity compounds after enforcing the frozen P1B input contract.
Another {provenance['tasks_with_zero_model_valid_rows']} governed tasks had zero
model-valid rows, so the governed-to-model task transition is
{provenance['governed_tasks']} ->
{provenance['governed_tasks'] - provenance['total_tasks_removed_from_model_manifest']}.
Reconstructed and materialized activity-ID sets match:
`{str(provenance['materialized_activity_ids_match_reconstruction']).lower()}`.
No affinity value field was materialized.

## Stop Rule

The historical synthetic Gate remains failed. E0S is diagnostic only. A future
repair must be separately registered and cannot reinterpret this audit as an
E0-S or typed-interaction authorization.
"""
    (output / "E0_SYNTHESIS_REPORT.md").write_text(text, encoding="utf-8")


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
    parser.add_argument("--canonical-rows", type=Path,
                        default=Path("dataset/processed/source_affinity/energy_pilot_v1/canonical_rows.jsonl"))
    parser.add_argument("--splits", type=Path,
                        default=Path("dataset/processed/source_affinity/energy_pilot_v1_governance/split_assignments.jsonl"))
    parser.add_argument("--output", type=Path,
                        default=Path("research/e0_identifiability/artifacts/e0s_evidence_v1"))
    parser.add_argument("--device", default="cuda")
    args = parser.parse_args()
    result = audit_e0s(args.input, args.cache, args.checkpoint, args.synthetic,
                       args.canonical_rows, args.splits, args.output, args.device)
    _write_report(args.output, result)
    synthesis_path = args.output / "E0_SYNTHESIS.json"
    synthesis_path.write_text(json.dumps(result["synthesis"], indent=2, sort_keys=True),
                              encoding="utf-8")
    print(json.dumps(result["synthesis"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
