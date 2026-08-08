"""Run E-AFF-P0 shared fixed-radial source-affinity feasibility."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import torch

from research.e0_identifiability.eaff_pilot_contract import (
    assert_paffinity_direction,
    component_bootstrap,
    component_macro_contrasts,
    coupling_null,
    fit_pair_ridge,
    task_metrics,
)
from research.e0_identifiability.run_synthetic_pregate import _load_bridge
from research.e0_identifiability.run_tbasis_radial import aggregate_basis
from scripts.govern_affinity_homology import local_identity
from scripts.source_affinity.common import sha256_file


STAGE = "P1R2B-E-AFF-P0_FIXED_RADIAL_SOURCE_FEASIBILITY"
SEED = 17
LIGANDS_PER_TASK = 20
RIDGE_ALPHA = 10.0
BATCH_SIZE = 16


def _read_jsonl(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _write_jsonl(path: Path, rows) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")


def _hash_key(prefix: str, *values: str) -> str:
    return hashlib.sha256("|".join((prefix, *values)).encode()).hexdigest()


def _group_observations(rows: list[dict]) -> tuple[dict, dict]:
    by_task: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    task_meta = {}
    for row in rows:
        by_task[row["task_id"]][row["ligand_state_key"]].append(row)
        task_meta.setdefault(row["task_id"], row)
    return by_task, task_meta


def select_panel(rows: list[dict]) -> tuple[list[dict], dict]:
    by_task, task_meta = _group_observations(rows)
    by_component: dict[str, list[str]] = defaultdict(list)
    for task, meta in task_meta.items():
        by_component[meta["closure_component_id"]].append(task)
    selected = []
    for component, tasks in sorted(by_component.items()):
        task = min(tasks, key=lambda value: _hash_key("EAFF-P0", value))
        ligand_states = sorted(
            by_task[task], key=lambda value: _hash_key("EAFF-P0-LIGAND", task, value)
        )[:LIGANDS_PER_TASK]
        if len(ligand_states) != LIGANDS_PER_TASK:
            raise ValueError(f"task {task} lacks 20 ligand states")
        for ligand in ligand_states:
            duplicates = by_task[task][ligand]
            row = dict(min(duplicates, key=lambda value: int(value["activity_id"])))
            row["activity_ids"] = sorted(int(value["activity_id"]) for value in duplicates)
            selected.append(row)
    audit = {
        "selection_policy": "one SHA256-selected task per closure; 20 SHA256-selected states",
        "rows": len(selected),
        "tasks": len({row["task_id"] for row in selected}),
        "closure_components": len({row["closure_component_id"] for row in selected}),
        "fold_tasks": dict(sorted(Counter(int(row["outer_oof_fold"]) for row in selected).items())),
    }
    return selected, audit


def load_affinity(rows: list[dict], canonical_path: Path) -> tuple[dict[int, dict], dict]:
    wanted = {int(row["activity_id"]) for row in rows}
    found = {}
    task_documents: dict[str, set[str]] = defaultdict(set)
    for row in _read_jsonl(canonical_path):
        activity = int(row["activity_id"])
        if activity not in wanted:
            continue
        found[activity] = {
            "p_affinity": float(row["p_affinity"]),
            "standard_value": float(row["standard_value"]),
            "standard_units": row["standard_units"],
            "document_chembl_id": row["document_chembl_id"],
            "assay_chembl_id": row["assay_chembl_id"],
        }
        task_documents[row["task_id"]].add(str(row["document_chembl_id"]))
    if set(found) != wanted:
        raise ValueError(f"canonical affinity mapping missing {len(wanted - set(found))} activities")
    p = np.asarray([value["p_affinity"] for value in found.values()])
    molar = np.asarray([value["standard_value"] * 1e-9 for value in found.values()])
    assert_paffinity_direction(p[np.argsort(molar)], molar[np.argsort(molar)])
    formula_error = max(abs(value["p_affinity"] + np.log10(value["standard_value"] * 1e-9))
                        for value in found.values())
    if formula_error > 1e-8:
        raise ValueError(f"pAffinity formula mismatch: {formula_error}")
    return found, {"mapped_activities": len(found), "paffinity_formula_max_error": formula_error,
                   "stronger_is_larger": True,
                   "task_documents": {key: sorted(value) for key, value in task_documents.items()}}


def make_observations(rows: list[dict], labels: dict[int, dict]) -> list[dict]:
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in rows:
        grouped[(row["task_id"], row["ligand_state_key"])].append(row)
    observations = []
    for (_, _), duplicates in sorted(grouped.items()):
        first = min(duplicates, key=lambda value: int(value["activity_id"]))
        observations.append({
            **first,
            "p_affinity": float(np.median([
                labels[int(row["activity_id"])]["p_affinity"] for row in duplicates
            ])),
            "measurement_count": len(duplicates),
        })
    return observations


def _load_ligand_cache(root: Path, pooled_keys: set[str], local_keys: set[str]):
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    pooled, local = {}, {}
    for item in manifest["ligand_shards"]:
        path = root / item["path"]
        if sha256_file(path) != item["sha256"]:
            raise ValueError(f"ligand cache hash mismatch: {path}")
        with np.load(path, allow_pickle=False) as shard:
            offsets = shard["offsets"]
            for index, raw_key in enumerate(shard["keys"]):
                key = str(raw_key)
                if key in pooled_keys:
                    pooled[key] = shard["pooled"][index].astype(np.float32)
                if key in local_keys:
                    left, right = int(offsets[index]), int(offsets[index + 1])
                    local[key] = {
                        "atoms": shard["atoms"][left:right].astype(np.float32),
                        "chemistry": shard["chemistry"][left:right].astype(np.float32),
                    }
    if set(pooled) != pooled_keys or set(local) != local_keys:
        raise ValueError("ligand cache does not cover E-AFF selection")
    return pooled, local, manifest


def _load_protein_cache(root: Path, keys: set[str]) -> dict:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    result = {}
    for item in manifest["protein_shards"]:
        path = root / item["path"]
        if sha256_file(path) != item["sha256"]:
            raise ValueError(f"protein cache hash mismatch: {path}")
        with np.load(path, allow_pickle=False) as shard:
            for index, raw_key in enumerate(shard["keys"]):
                key = str(raw_key)
                if key in keys:
                    result[key] = {
                        "residues": shard["residues"][index].astype(np.float32),
                        "chemistry": shard["chemistry"][index].astype(np.float32),
                        "mask": shard["mask"][index].astype(bool),
                    }
    if set(result) != keys:
        raise ValueError("protein cache does not cover E-AFF proteins")
    return result


def _weighted_scale(x: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    mean = np.average(x, axis=0, weights=weights)
    variance = np.average(np.square(x - mean), axis=0, weights=weights)
    scale = np.sqrt(variance)
    return mean, np.where(scale > 1e-8, scale, 1.0)


def ligand_oof(observations: list[dict], pooled: dict[str, np.ndarray]) -> tuple[np.ndarray, list[dict]]:
    from sklearn.linear_model import Ridge

    x = np.stack([pooled[row["ligand_state_key"]] for row in observations]).astype(np.float64)
    y = np.asarray([row["p_affinity"] for row in observations], dtype=np.float64)
    folds = np.asarray([int(row["outer_oof_fold"]) for row in observations])
    tasks = np.asarray([row["task_id"] for row in observations])
    counts = Counter(tasks.tolist())
    weights = np.asarray([1.0 / counts[task] for task in tasks])
    prediction = np.full(len(y), np.nan, dtype=np.float64)
    diagnostics = []
    for fold in range(5):
        train, test = folds != fold, folds == fold
        mean, scale = _weighted_scale(x[train], weights[train])
        model = Ridge(alpha=RIDGE_ALPHA, fit_intercept=True, solver="lsqr")
        model.fit((x[train] - mean) / scale, y[train], sample_weight=weights[train])
        prediction[test] = model.predict((x[test] - mean) / scale)
        diagnostics.append({"fold": fold, "train": int(train.sum()), "test": int(test.sum()),
                            "alpha": RIDGE_ALPHA})
    if np.any(~np.isfinite(prediction)):
        raise ValueError("OOF ligand prior left non-finite predictions")
    return prediction, diagnostics


def build_derangement(selected: list[dict], input_root: Path,
                      governance_root: Path) -> tuple[dict[str, str], list[dict]]:
    sequences = {row["sequence_sha256"]: row["sequence"]
                 for row in _read_jsonl(input_root / "proteins.jsonl")}
    homology = {row["protein_sequence_sha256"]: row["homology_component_id"]
                for row in _read_jsonl(governance_root / "homology_assignments.jsonl")
                if not row["excluded_by_davis_protected_homology"]}
    protein_fold: dict[str, int] = {}
    for row in _read_jsonl(input_root / "rows.label_blind.jsonl"):
        protein_fold.setdefault(row["protein_sequence_sha256"], int(row["outer_oof_fold"]))
    selected_proteins = sorted({row["protein_sequence_sha256"] for row in selected})
    used, mapping, audit = set(), {}, []
    candidates_by_fold = defaultdict(list)
    for protein, fold in protein_fold.items():
        if protein in sequences and protein in homology:
            candidates_by_fold[fold].append(protein)
    for correct in selected_proteins:
        candidates = [candidate for candidate in candidates_by_fold[protein_fold[correct]]
                      if candidate != correct and candidate not in used
                      and homology[candidate] != homology[correct]]
        candidates.sort(key=lambda value: _hash_key("EAFF-P0-WRONG", correct, value))
        chosen, identity = None, None
        for candidate in candidates:
            _, value = local_identity((sequences[correct], sequences[candidate]))
            if value < 0.40:
                chosen, identity = candidate, float(value)
                break
        if chosen is None:
            raise ValueError(f"no score-blind derangement for {correct}")
        used.add(chosen)
        mapping[correct] = chosen
        audit.append({"correct_protein": correct, "wrong_protein": chosen,
                      "fold": protein_fold[correct], "local_identity": identity,
                      "correct_homology": homology[correct],
                      "wrong_homology": homology[chosen]})
    return mapping, audit


def _geometry_basis(rows: list[dict], bridge, proteins: dict, ligands: dict,
                    wrong_map: dict[str, str], bin_moments: np.ndarray,
                    calibration_coef: np.ndarray, calibration_intercept: np.ndarray,
                    device: str) -> tuple[np.ndarray, np.ndarray]:
    correct_values, wrong_values = [], []
    bridge.eval()
    with torch.inference_mode():
        for start in range(0, len(rows), BATCH_SIZE):
            subset = rows[start:start + BATCH_SIZE]
            width = max(len(ligands[row["ligand_state_key"]]["atoms"]) for row in subset)
            atom_state = np.zeros((len(subset), width, 128), dtype=np.float32)
            atom_mask = np.zeros((len(subset), width), dtype=np.float32)
            correct_residue, correct_mask, wrong_residue, wrong_mask = [], [], [], []
            for index, row in enumerate(subset):
                ligand = ligands[row["ligand_state_key"]]
                count = len(ligand["atoms"])
                atom_state[index, :count] = ligand["atoms"]
                atom_mask[index, :count] = 1
                correct = proteins[row["protein_sequence_sha256"]]
                wrong = proteins[wrong_map[row["protein_sequence_sha256"]]]
                correct_residue.append(correct["residues"])
                correct_mask.append(correct["mask"])
                wrong_residue.append(wrong["residues"])
                wrong_mask.append(wrong["mask"])
            atoms = torch.from_numpy(atom_state).to(device)
            masks = torch.from_numpy(atom_mask).to(device)
            correct_output = bridge(
                atoms, masks, torch.from_numpy(np.stack(correct_residue)).to(device),
                torch.from_numpy(np.stack(correct_mask).astype(np.float32)).to(device))
            wrong_output = bridge(
                atoms, masks, torch.from_numpy(np.stack(wrong_residue)).to(device),
                torch.from_numpy(np.stack(wrong_mask).astype(np.float32)).to(device))
            correct_distance = correct_output.distance_prob.cpu().numpy()
            wrong_distance = wrong_output.distance_prob.cpu().numpy()
            for index, row in enumerate(subset):
                ligand = ligands[row["ligand_state_key"]]
                count = len(ligand["atoms"])
                atom_channels = ligand["chemistry"][:, 32:40].astype(np.float64)
                correct = proteins[row["protein_sequence_sha256"]]
                wrong = proteins[wrong_map[row["protein_sequence_sha256"]]]
                raw_correct = aggregate_basis(
                    atom_channels, correct["chemistry"].astype(np.float64),
                    np.einsum("isb,bk->isk", correct_distance[index, :count], bin_moments),
                    correct["mask"])
                raw_wrong = aggregate_basis(
                    atom_channels, wrong["chemistry"].astype(np.float64),
                    np.einsum("isb,bk->isk", wrong_distance[index, :count], bin_moments),
                    wrong["mask"])
                correct_values.append(raw_correct.reshape(-1, 6) @ calibration_coef.T
                                      + calibration_intercept)
                wrong_values.append(raw_wrong.reshape(-1, 6) @ calibration_coef.T
                                    + calibration_intercept)
            print(f"basis_rows={min(start + BATCH_SIZE, len(rows))}/{len(rows)}", flush=True)
    correct = np.stack(correct_values).reshape(-1, 8, 6, 6)
    wrong = np.stack(wrong_values).reshape(-1, 8, 6, 6)
    return correct, wrong


def h0_census(observations: list[dict], label_meta: dict) -> dict:
    by_task: dict[str, set[str]] = defaultdict(set)
    task_meta = {}
    for row in observations:
        by_task[row["task_id"]].add(row["ligand_state_key"])
        task_meta.setdefault(row["task_id"], row)
    deep = [task for task, ligands in by_task.items() if len(ligands) >= 40]
    deep_components = {task_meta[task]["closure_component_id"] for task in deep}
    by_target_endpoint: dict[tuple[str, str], set[str]] = defaultdict(set)
    documents: dict[tuple[str, str], set[str]] = defaultdict(set)
    for task, meta in task_meta.items():
        key = (meta["protein_sequence_sha256"], meta["endpoint_family"])
        by_target_endpoint[key].add(task)
        documents[key].update(label_meta["task_documents"].get(task, []))
    transport = [key for key, tasks in by_target_endpoint.items()
                 if len(tasks) >= 2 and len(documents[key]) >= 2]
    return {
        "tasks_ge_40_ligands": len(deep),
        "closure_components_with_task_ge_40": len(deep_components),
        "target_endpoint_groups_with_ge_2_tasks_and_documents": len(transport),
        "h0a_data_supported": len(deep_components) >= 16,
        "h0b_data_supported": len(transport) >= 16,
    }


def finalize_manifest(args, output: Path) -> dict:
    report_path = output / "report.json"
    if not report_path.exists():
        raise FileNotFoundError("cannot finalize E-AFF artifact without report.json")
    report = json.loads(report_path.read_text(encoding="utf-8"))
    selection = list(_read_jsonl(output / "selection.jsonl"))
    report["selection"].pop("fold_tasks", None)
    report["selection"]["fold_rows"] = dict(sorted(Counter(
        int(row["outer_oof_fold"]) for row in selection).items()))
    report["selection"]["fold_tasks"] = dict(sorted(Counter(
        int(next(row["outer_oof_fold"] for row in selection if row["task_id"] == task))
        for task in {row["task_id"] for row in selection}).items()))
    _write_json(report_path, report)
    input_root = Path(args.input)
    manifest = {
        "stage": STAGE, "research_only": True,
        "inputs": {
            "input_manifest": sha256_file(input_root / "manifest.json"),
            "governance_manifest": sha256_file(Path(args.governance) / "governance_manifest.json"),
            "canonical_rows": sha256_file(Path(args.canonical_rows)),
            "cache_manifest": sha256_file(Path(args.cache) / "manifest.json"),
            "checkpoint": sha256_file(Path(args.checkpoint)),
            "tbasis_values": sha256_file(Path(args.tbasis_values)),
            "preregistration": sha256_file(Path(__file__).with_name("EAFF_P0_PREREGISTRATION.md")),
        },
        "outputs": {name: sha256_file(output / name) for name in (
            "selection.jsonl", "derangement.jsonl", "features.npz", "directions.npy",
            "task_metrics.jsonl", "report.json")},
        "label_reads": {"source_affinity": True, "davis": 0, "recipient": 0},
    }
    _write_json(output / "manifest.json", manifest)
    return report


def run(args) -> dict:
    if not args.device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("E-AFF-P0 requires CUDA for frozen P1B inference")
    output = Path(args.output)
    if args.finalize_only:
        return finalize_manifest(args, output)
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    output.mkdir(parents=True)
    input_root = Path(args.input)
    all_rows = list(_read_jsonl(input_root / "rows.label_blind.jsonl"))
    selected, selection_audit = select_panel(all_rows)
    _write_jsonl(output / "selection.jsonl", selected)

    labels, label_audit = load_affinity(all_rows, Path(args.canonical_rows))
    observations = make_observations(all_rows, labels)
    obs_lookup = {(row["task_id"], row["ligand_state_key"]): index
                  for index, row in enumerate(observations)}
    selected_observations = [observations[obs_lookup[(row["task_id"], row["ligand_state_key"])]]
                             for row in selected]

    pooled_keys = {row["ligand_state_key"] for row in observations}
    selected_ligands = {row["ligand_state_key"] for row in selected_observations}
    pooled, local_ligands, cache_manifest = _load_ligand_cache(
        Path(args.cache), pooled_keys, selected_ligands)
    ligand_prediction, ligand_diagnostics = ligand_oof(observations, pooled)
    baseline = np.asarray([ligand_prediction[obs_lookup[(row["task_id"], row["ligand_state_key"])]]
                           for row in selected_observations])

    wrong_map, wrong_audit = build_derangement(
        selected_observations, input_root, Path(args.governance))
    _write_jsonl(output / "derangement.jsonl", wrong_audit)
    protein_keys = ({row["protein_sequence_sha256"] for row in selected_observations}
                    | set(wrong_map.values()))
    proteins = _load_protein_cache(Path(args.cache), protein_keys)
    bridge = _load_bridge(Path(args.checkpoint), args.device)
    with np.load(args.tbasis_values, allow_pickle=False) as tbasis:
        correct, deranged = _geometry_basis(
            selected_observations, bridge, proteins, local_ligands, wrong_map,
            tbasis["bin_rbf_expectation"].astype(np.float64),
            tbasis["calibration_coef"].astype(np.float64),
            tbasis["calibration_intercept"].astype(np.float64), args.device)
    null = coupling_null(correct)
    np.savez_compressed(output / "features.npz", correct=correct.astype(np.float32),
                        deranged=deranged.astype(np.float32), null=null.astype(np.float32),
                        ligand_oof=baseline.astype(np.float32))

    y = np.asarray([row["p_affinity"] for row in selected_observations])
    residual = y - baseline
    folds = np.asarray([int(row["outer_oof_fold"]) for row in selected_observations])
    task_ids = np.asarray([row["task_id"] for row in selected_observations])
    energy_correct = np.full(len(y), np.nan)
    energy_deranged = np.full(len(y), np.nan)
    energy_null = np.full(len(y), np.nan)
    directions, fit_diagnostics = [], []
    flat_correct, flat_deranged, flat_null = [value.reshape(len(value), -1)
                                               for value in (correct, deranged, null)]
    for fold in range(5):
        train, test = folds != fold, folds == fold
        direction, diagnostic = fit_pair_ridge(
            flat_correct[train], residual[train], task_ids[train], RIDGE_ALPHA)
        energy_correct[test] = flat_correct[test] @ direction
        energy_deranged[test] = flat_deranged[test] @ direction
        energy_null[test] = flat_null[test] @ direction
        directions.append(direction)
        fit_diagnostics.append({"fold": fold, "train_rows": int(train.sum()),
                                "test_rows": int(test.sum()), **diagnostic})
    if np.any(~np.isfinite(energy_correct + energy_deranged + energy_null)):
        raise ValueError("shared direction left non-finite OOF predictions")
    np.save(output / "directions.npy", np.stack(directions))
    predictions = {
        "ligand": baseline,
        "correct": baseline + energy_correct,
        "deranged": baseline + energy_deranged,
        "null": baseline + energy_null,
    }
    per_task = task_metrics(selected_observations, y, predictions)
    _write_jsonl(output / "task_metrics.jsonl", per_task)
    primary = component_macro_contrasts(per_task)
    confidence = component_bootstrap(per_task, seed=SEED)
    conditions = {
        "correct_minus_ligand_ge_0_03": primary["correct_minus_ligand"] >= 0.03,
        "correct_minus_ligand_lcb_positive": confidence["correct_minus_ligand"][0] > 0,
        "correct_minus_deranged_ge_0_03": primary["correct_minus_deranged"] >= 0.03,
        "correct_minus_deranged_lcb_positive": confidence["correct_minus_deranged"][0] > 0,
        "correct_minus_null_positive": primary["correct_minus_null"] > 0,
        "correct_minus_null_lcb_positive": confidence["correct_minus_null"][0] > 0,
    }
    h0 = h0_census(observations, label_audit)
    if all(conditions.values()):
        verdict = "SHARED_RADIAL_AFFINITY_FEASIBILITY_OBSERVED"
    elif h0["h0a_data_supported"] or h0["h0b_data_supported"]:
        verdict = "SHARED_DIRECTION_NOT_OBSERVED_H0_DATA_SUPPORTED"
    else:
        verdict = "SHARED_DIRECTION_NOT_OBSERVED_H0_DATA_INSUFFICIENT"
    endpoint = {}
    for name in ("Ki", "Kd"):
        subset = [row for row in per_task if row["endpoint_family"] == name]
        endpoint[name] = component_macro_contrasts(subset) if subset else {"components": 0}
    report = {
        "schema": "MetaSieve.EAffPilot.v1", "stage": STAGE,
        "created_utc": datetime.now(timezone.utc).isoformat(), "verdict": verdict,
        "research_only": True, "affinity_labels_read": True,
        "davis_label_reads": 0, "recipient_label_reads": 0,
        "selection": selection_audit, "label_audit": {key: value for key, value in label_audit.items()
                                                       if key != "task_documents"},
        "ligand_prior": {"type": "closure-OOF pooled-state Ridge", "alpha": RIDGE_ALPHA,
                           "folds": ligand_diagnostics},
        "shared_direction": {"type": "task-balanced residual-difference Ridge",
                             "alpha": RIDGE_ALPHA, "folds": fit_diagnostics},
        "primary_component_macro": primary, "component_bootstrap_ci95": confidence,
        "gate_conditions": conditions, "endpoint_secondary": endpoint,
        "h0_data_census": h0,
        "interpretation_limits": [
            "lightweight source-feasibility panel, not terminal E-AFF",
            "coupling null is an attribution control, not a nonbinder",
            "shared-w failure does not establish radial-basis insufficiency",
            "task-local headroom does not establish target biology",
        ],
    }
    _write_json(output / "report.json", report)
    return finalize_manifest(args, output)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="dataset/processed/source_affinity/e0_input_v1")
    parser.add_argument("--canonical-rows", default="dataset/processed/source_affinity/energy_pilot_v1/canonical_rows.jsonl")
    parser.add_argument("--governance", default="dataset/processed/source_affinity/energy_pilot_v1_governance")
    parser.add_argument("--cache", default="dataset/processed/source_affinity/e0_local_states_v1")
    parser.add_argument("--checkpoint", default="report/mechanism_refactor/p1b_pilot20k_seed17_v1/best.pt")
    parser.add_argument("--tbasis-values", default="research/e0_identifiability/artifacts/tbasis_r0_v1/basis_values.npz")
    parser.add_argument("--output", default="research/e0_identifiability/artifacts/eaff_p0_v1")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--finalize-only", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), indent=2, sort_keys=True))
