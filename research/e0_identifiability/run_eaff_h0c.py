"""Run E-AFF-H0C support-matched fixed-radial interaction residual."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import torch
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import Ridge

from research.e0_identifiability.eaff_h0c_contract import (
    centered_interaction,
    component_bootstrap,
    component_summary,
)
from research.e0_identifiability.eaff_pilot_contract import fit_pair_ridge
from research.e0_identifiability.metrics import concordance
from research.e0_identifiability.run_eaff_pilot import (
    RIDGE_ALPHA,
    _geometry_basis,
    _hash_key,
    _load_ligand_cache,
    _load_protein_cache,
    _read_jsonl,
    _write_json,
    _write_jsonl,
    build_derangement,
    ligand_oof,
    load_affinity,
    make_observations,
)
from research.e0_identifiability.run_synthetic_pregate import _load_bridge
from scripts.source_affinity.common import sha256_file


STAGE = "P1R2B-E-AFF-H0C_FIXED_RADIAL_INTERACTION_RESIDUAL"
SUPPORT = 20
TEST = 20
LOCAL_FOLDS = 5


def _scaffolds(ligand_path: Path) -> dict[str, str]:
    result = {}
    for row in _read_jsonl(ligand_path):
        molecule = Chem.MolFromSmiles(row["canonical_smiles"])
        if molecule is None:
            raise ValueError(f"invalid governed ligand: {row['ligand_state_key']}")
        core = MurckoScaffold.GetScaffoldForMol(molecule)
        value = Chem.MolToSmiles(core, canonical=True, isomericSmiles=False)
        result[row["ligand_state_key"]] = value or "__EMPTY_MURCKO__"
    return result


def _scaffold_partition(task: str, ligands: dict[str, list[dict]],
                        scaffold: dict[str, str]) -> tuple[list[str], list[str]] | None:
    groups: dict[str, list[str]] = defaultdict(list)
    for ligand in ligands:
        groups[scaffold[ligand]].append(ligand)
    ordered = sorted(groups, key=lambda value: _hash_key("H0C-SCAFFOLD", task, value))
    capacities = [len(groups[value]) for value in ordered]
    total = sum(capacities)
    states: dict[int, tuple[int, ...]] = {0: ()}
    for index, capacity in enumerate(capacities):
        for value, chosen in list(states.items()):
            states.setdefault(value + capacity, chosen + (index,))
    eligible = [value for value in states if SUPPORT <= value <= total - TEST]
    if not eligible:
        return None
    amount = min(eligible, key=lambda value: (abs(value - SUPPORT), states[value]))
    support_groups = {ordered[index] for index in states[amount]}
    support = [ligand for group in support_groups for ligand in groups[group]]
    test = [ligand for group in ordered if group not in support_groups for ligand in groups[group]]
    support = sorted(support, key=lambda value: _hash_key("H0C-SUPPORT", task, value))[:SUPPORT]
    test = sorted(test, key=lambda value: _hash_key("H0C-TEST", task, value))[:TEST]
    if len(support) != SUPPORT or len(test) != TEST:
        return None
    if {scaffold[value] for value in support} & {scaffold[value] for value in test}:
        raise RuntimeError("scaffold partition overlap")
    return support, test


def select_panel(rows: list[dict], scaffold: dict[str, str],
                 h0a_selection: Path) -> tuple[list[dict], dict]:
    used_tasks = {row["task_id"] for row in _read_jsonl(h0a_selection)}
    by_task: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    task_meta = {}
    for row in rows:
        if row["task_id"] in used_tasks:
            continue
        by_task[row["task_id"]][row["ligand_state_key"]].append(row)
        task_meta.setdefault(row["task_id"], row)
    partitions = {}
    by_component: dict[str, list[str]] = defaultdict(list)
    for task, ligands in by_task.items():
        if len(ligands) < SUPPORT + TEST:
            continue
        partition = _scaffold_partition(task, ligands, scaffold)
        if partition is not None:
            partitions[task] = partition
            by_component[task_meta[task]["closure_component_id"]].append(task)
    selected = []
    for component, tasks in sorted(by_component.items()):
        task = min(tasks, key=lambda value: _hash_key("EAFF-H0C", value))
        support, test = partitions[task]
        for partition, ligand_states in (("support", support), ("test", test)):
            for ligand in ligand_states:
                duplicates = by_task[task][ligand]
                row = dict(min(duplicates, key=lambda value: int(value["activity_id"])))
                row["activity_ids"] = sorted(int(value["activity_id"]) for value in duplicates)
                row["h0c_partition"] = partition
                row["murcko_scaffold"] = scaffold[ligand]
                selected.append(row)
    audit = {
        "h0a_tasks_excluded": len(used_tasks),
        "eligible_tasks": len(partitions),
        "tasks": len(by_component),
        "closure_components": len(by_component),
        "rows": len(selected),
        "support_rows": sum(row["h0c_partition"] == "support" for row in selected),
        "test_rows": sum(row["h0c_partition"] == "test" for row in selected),
        "fold_tasks": dict(sorted(Counter(int(row["outer_oof_fold"])
                                           for row in selected[::SUPPORT + TEST]).items())),
        "scaffold_overlap_tasks": sum(bool(
            {row["murcko_scaffold"] for row in selected
             if row["task_id"] == task and row["h0c_partition"] == "support"}
            & {row["murcko_scaffold"] for row in selected
               if row["task_id"] == task and row["h0c_partition"] == "test"})
            for task in {row["task_id"] for row in selected}),
    }
    return selected, audit


def _scaled_ridge(train_x: np.ndarray, train_y: np.ndarray, predict_x: np.ndarray) -> np.ndarray:
    mean = train_x.mean(0)
    scale = train_x.std(0)
    scale = np.where(scale > 1e-8, scale, 1.0)
    model = Ridge(alpha=RIDGE_ALPHA, fit_intercept=True, solver="svd")
    model.fit((train_x - mean) / scale, train_y)
    return model.predict((predict_x - mean) / scale)


def local_ligand_nuisance(task: str, keys: list[str], support_x: np.ndarray,
                          support_y: np.ndarray, test_x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = sorted(range(len(keys)), key=lambda index: _hash_key("H0C-LOCAL-FOLD", task, keys[index]))
    folds = np.empty(len(keys), dtype=int)
    for rank, index in enumerate(order):
        folds[index] = rank % LOCAL_FOLDS
    oof = np.full(len(keys), np.nan, dtype=np.float64)
    for fold in range(LOCAL_FOLDS):
        train, held = folds != fold, folds == fold
        oof[held] = _scaled_ridge(support_x[train], support_y[train], support_x[held])
    if np.any(~np.isfinite(oof)):
        raise RuntimeError("local ligand nuisance left non-finite support predictions")
    return oof, _scaled_ridge(support_x, support_y, test_x)


def _task_metrics(rows: list[dict], labels: np.ndarray, predictions: dict[str, np.ndarray]) -> list[dict]:
    grouped: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        grouped[row["task_id"]].append(index)
    result = []
    for task, indices in sorted(grouped.items()):
        first = rows[indices[0]]
        result.append({
            "task_id": task,
            "closure_component_id": first["closure_component_id"],
            "endpoint_family": first["endpoint_family"],
            **{name: concordance(labels[indices], values[indices])
               for name, values in predictions.items()},
        })
    return result


def run(args) -> dict:
    if not args.device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("E-AFF-H0C requires CUDA")
    output = Path(args.output)
    if output.exists():
        if any(output.iterdir()):
            raise FileExistsError(f"output already exists and is not empty: {output}")
    else:
        output.mkdir(parents=True)
    input_root = Path(args.input)
    all_rows = list(_read_jsonl(input_root / "rows.label_blind.jsonl"))
    scaffold = _scaffolds(input_root / "ligands.jsonl")
    selected, selection_audit = select_panel(all_rows, scaffold, Path(args.h0a_selection))
    if selection_audit["closure_components"] < 32 or selection_audit["scaffold_overlap_tasks"]:
        raise RuntimeError("H0C label-blind selection contract failed")
    _write_jsonl(output / "selection.jsonl", selected)

    labels, label_audit = load_affinity(all_rows, Path(args.canonical_rows))
    observations = make_observations(all_rows, labels)
    lookup = {(row["task_id"], row["ligand_state_key"]): index
              for index, row in enumerate(observations)}
    selected_obs = []
    for row in selected:
        value = dict(observations[lookup[(row["task_id"], row["ligand_state_key"])]])
        value["h0c_partition"] = row["h0c_partition"]
        value["murcko_scaffold"] = row["murcko_scaffold"]
        selected_obs.append(value)

    pooled_keys = {row["ligand_state_key"] for row in observations}
    local_keys = {row["ligand_state_key"] for row in selected_obs}
    pooled, local_ligands, _ = _load_ligand_cache(Path(args.cache), pooled_keys, local_keys)
    global_prediction, ligand_diagnostics = ligand_oof(observations, pooled)
    global_baseline = np.asarray([
        global_prediction[lookup[(row["task_id"], row["ligand_state_key"])]]
        for row in selected_obs], dtype=np.float64)
    pooled_selected = np.stack([pooled[row["ligand_state_key"]] for row in selected_obs]).astype(np.float64)

    wrong_map, wrong_audit = build_derangement(selected_obs, input_root, Path(args.governance))
    _write_jsonl(output / "derangement.jsonl", wrong_audit)
    protein_keys = ({row["protein_sequence_sha256"] for row in selected_obs}
                    | set(wrong_map.values()))
    proteins = _load_protein_cache(Path(args.cache), protein_keys)
    bridge = _load_bridge(Path(args.checkpoint), args.device)
    with np.load(args.tbasis_values, allow_pickle=False) as tbasis:
        correct_phi, deranged_phi = _geometry_basis(
            selected_obs, bridge, proteins, local_ligands, wrong_map,
            tbasis["bin_rbf_expectation"].astype(np.float64),
            tbasis["calibration_coef"].astype(np.float64),
            tbasis["calibration_intercept"].astype(np.float64), args.device)
    correct_psi = centered_interaction(correct_phi).reshape(len(selected_obs), -1)
    deranged_psi = centered_interaction(deranged_phi).reshape(len(selected_obs), -1)
    y = np.asarray([row["p_affinity"] for row in selected_obs], dtype=np.float64)

    local_effect = np.full(len(y), np.nan, dtype=np.float64)
    correct_effect = np.full(len(y), np.nan, dtype=np.float64)
    deranged_effect = np.full(len(y), np.nan, dtype=np.float64)
    direction_tasks, directions, fit_diagnostics = [], [], []
    for task in sorted({row["task_id"] for row in selected_obs}):
        indices = np.asarray([index for index, row in enumerate(selected_obs)
                              if row["task_id"] == task])
        support = indices[[selected_obs[index]["h0c_partition"] == "support" for index in indices]]
        test = indices[[selected_obs[index]["h0c_partition"] == "test" for index in indices]]
        support_target = y[support] - global_baseline[support]
        nuisance_oof, nuisance_test = local_ligand_nuisance(
            task, [selected_obs[index]["ligand_state_key"] for index in support],
            pooled_selected[support], support_target, pooled_selected[test])
        interaction_target = support_target - nuisance_oof
        direction, diagnostic = fit_pair_ridge(
            correct_psi[support], interaction_target,
            np.asarray([task] * len(support)), RIDGE_ALPHA)
        local_effect[test] = nuisance_test
        correct_effect[test] = correct_psi[test] @ direction
        deranged_effect[test] = deranged_psi[test] @ direction
        direction_tasks.append(task)
        directions.append(direction)
        fit_diagnostics.append({"task_id": task, "support": len(support), "test": len(test),
                                "local_oof_rmse": float(np.sqrt(np.mean(
                                    np.square(nuisance_oof - support_target)))), **diagnostic})

    test_indices = np.asarray([index for index, row in enumerate(selected_obs)
                               if row["h0c_partition"] == "test"])
    if any(np.any(~np.isfinite(value[test_indices]))
           for value in (local_effect, correct_effect, deranged_effect)):
        raise RuntimeError("H0C left non-finite test predictions")
    test_rows = [selected_obs[index] for index in test_indices]
    local_score = global_baseline[test_indices] + local_effect[test_indices]
    predictions = {
        "global_ligand": global_baseline[test_indices],
        "local_ligand": local_score,
        "correct": local_score + correct_effect[test_indices],
        "deranged": local_score + deranged_effect[test_indices],
    }
    per_task = _task_metrics(test_rows, y[test_indices], predictions)
    primary = component_summary(per_task)
    confidence = component_bootstrap(per_task)
    interaction = (primary["correct_minus_local_ligand"] >= 0.03
                   and confidence["correct_minus_local_ligand"][0] > 0)
    partner = (primary["correct_minus_deranged"] >= 0.03
               and confidence["correct_minus_deranged"][0] > 0)
    weak_partner = (primary["correct_minus_deranged"] > 0
                    and confidence["correct_minus_deranged"][0] > 0)
    if interaction and partner:
        verdict = "FIXED_RADIAL_INTERACTION_RESIDUAL_AND_PARTNER_SIGNAL_OBSERVED"
    elif interaction:
        verdict = "FIXED_RADIAL_INTERACTION_RESIDUAL_WITHOUT_PARTNER_SPECIFICITY"
    elif weak_partner:
        verdict = "WEAK_PARTNER_CONDITIONED_INCREMENT_BELOW_GATE"
    else:
        verdict = "FIXED_RADIAL_INTERACTION_RESIDUAL_NOT_OBSERVED"

    endpoint = {}
    for name in ("Ki", "Kd"):
        subset = [row for row in per_task if row["endpoint_family"] == name]
        endpoint[name] = component_summary(subset) if subset else {"components": 0}
    chemistry_error = max(float(np.abs(centered_interaction(correct_phi).sum(-1)).max()),
                          float(np.abs(centered_interaction(deranged_phi).sum(-1)).max()))
    radial_error = max(float(np.abs(centered_interaction(correct_phi).sum((-3, -2))).max()),
                       float(np.abs(centered_interaction(deranged_phi).sum((-3, -2))).max()))
    report = {
        "schema": "MetaSieve.EAffH0C.v1", "stage": STAGE,
        "created_utc": datetime.now(timezone.utc).isoformat(), "verdict": verdict,
        "research_only": True, "affinity_labels_read": True,
        "davis_label_reads": 0, "recipient_label_reads": 0,
        "selection": selection_audit,
        "label_audit": {key: value for key, value in label_audit.items()
                        if key != "task_documents"},
        "global_ligand_prior": {"alpha": RIDGE_ALPHA, "folds": ligand_diagnostics},
        "local_nuisance": {"features": "frozen_128D_pooled_ligand_state",
                           "alpha": RIDGE_ALPHA, "support_folds": LOCAL_FOLDS},
        "interaction_head": {"features": "double_centered_288D_radial_tensor",
                             "alpha": RIDGE_ALPHA, "fits": fit_diagnostics},
        "interaction_contract": {"chemistry_marginal_max_abs": chemistry_error,
                                 "radial_marginal_max_abs": radial_error,
                                 "minimum_correct_total": float(correct_phi.sum((1, 2, 3)).min()),
                                 "minimum_deranged_total": float(deranged_phi.sum((1, 2, 3)).min())},
        "primary_component_macro": primary,
        "component_bootstrap_ci95": confidence,
        "conditions": {"interaction_value": interaction, "partner_specificity": partner},
        "endpoint_secondary": endpoint,
        "interpretation_limits": [
            "H0A informed this hypothesis; closure families are development evidence",
            "the centered tensor is algebraic and is not a probability or physical energy",
            "a PASS cannot directly authorize RFSA or production integration",
        ],
    }
    _write_jsonl(output / "task_metrics.jsonl", per_task)
    np.savez_compressed(output / "task_directions.npz", task=np.asarray(direction_tasks),
                        direction=np.stack(directions).astype(np.float64))
    np.savez_compressed(
        output / "features.npz",
        correct_phi=correct_phi.astype(np.float32), deranged_phi=deranged_phi.astype(np.float32),
        global_ligand=global_baseline.astype(np.float32), local_effect=local_effect.astype(np.float32),
        correct_effect=correct_effect.astype(np.float32), deranged_effect=deranged_effect.astype(np.float32),
        label=y.astype(np.float32),
        partition=np.asarray([row["h0c_partition"] for row in selected_obs]))
    _write_json(output / "report.json", report)
    manifest = {
        "stage": STAGE, "research_only": True,
        "inputs": {
            "input_manifest": sha256_file(input_root / "manifest.json"),
            "governance_manifest": sha256_file(Path(args.governance) / "governance_manifest.json"),
            "canonical_rows": sha256_file(Path(args.canonical_rows)),
            "cache_manifest": sha256_file(Path(args.cache) / "manifest.json"),
            "checkpoint": sha256_file(Path(args.checkpoint)),
            "tbasis_values": sha256_file(Path(args.tbasis_values)),
            "h0a_selection": sha256_file(Path(args.h0a_selection)),
            "preregistration": sha256_file(Path(__file__).with_name("EAFF_H0C_PREREGISTRATION.md")),
        },
        "outputs": {name: sha256_file(output / name) for name in (
            "selection.jsonl", "derangement.jsonl", "features.npz",
            "task_directions.npz", "task_metrics.jsonl", "report.json")},
        "label_reads": {"source_affinity": True, "davis": 0, "recipient": 0},
    }
    _write_json(output / "manifest.json", manifest)
    return report


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="dataset/processed/source_affinity/e0_input_v1")
    parser.add_argument("--canonical-rows", default="dataset/processed/source_affinity/energy_pilot_v1/canonical_rows.jsonl")
    parser.add_argument("--governance", default="dataset/processed/source_affinity/energy_pilot_v1_governance")
    parser.add_argument("--cache", default="dataset/processed/source_affinity/e0_local_states_v1")
    parser.add_argument("--checkpoint", default="report/mechanism_refactor/p1b_pilot20k_seed17_v1/best.pt")
    parser.add_argument("--tbasis-values", default="research/e0_identifiability/artifacts/tbasis_r0_v1/basis_values.npz")
    parser.add_argument("--h0a-selection", default="research/e0_identifiability/artifacts/eaff_h0a_v1/selection.jsonl")
    parser.add_argument("--output", default="research/e0_identifiability/artifacts/eaff_h0c_v1")
    parser.add_argument("--device", default="cuda")
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), indent=2, sort_keys=True))
