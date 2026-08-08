"""Run E-AFF-H0A task-local fixed-radial headroom diagnostic."""
from __future__ import annotations

from collections import Counter, defaultdict
from datetime import datetime, timezone
import argparse
import json
from pathlib import Path

import numpy as np
import torch

from research.e0_identifiability.eaff_pilot_contract import (
    component_bootstrap,
    component_macro_contrasts,
    coupling_null,
    fit_pair_ridge,
    task_metrics,
)
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


STAGE = "P1R2B-E-AFF-H0A_TASK_LOCAL_RADIAL_HEADROOM"
LIGANDS = 40
FIT = 20


def select_h0a(rows: list[dict]) -> tuple[list[dict], dict]:
    by_task: dict[str, dict[str, list[dict]]] = defaultdict(lambda: defaultdict(list))
    task_meta = {}
    for row in rows:
        by_task[row["task_id"]][row["ligand_state_key"]].append(row)
        task_meta.setdefault(row["task_id"], row)
    by_component: dict[str, list[str]] = defaultdict(list)
    for task, ligands in by_task.items():
        if len(ligands) >= LIGANDS:
            by_component[task_meta[task]["closure_component_id"]].append(task)
    selected = []
    for component, tasks in sorted(by_component.items()):
        task = min(tasks, key=lambda value: _hash_key("EAFF-H0A", value))
        ligand_states = sorted(
            by_task[task], key=lambda value: _hash_key("EAFF-H0A-LIGAND", task, value)
        )[:LIGANDS]
        for index, ligand in enumerate(ligand_states):
            duplicates = by_task[task][ligand]
            row = dict(min(duplicates, key=lambda value: int(value["activity_id"])))
            row["activity_ids"] = sorted(int(value["activity_id"]) for value in duplicates)
            row["h0_partition"] = "fit" if index < FIT else "test"
            selected.append(row)
    return selected, {
        "tasks": len(by_component), "closure_components": len(by_component),
        "rows": len(selected), "fit_rows": sum(row["h0_partition"] == "fit" for row in selected),
        "test_rows": sum(row["h0_partition"] == "test" for row in selected),
        "fold_tasks": dict(sorted(Counter(int(rows[0]["outer_oof_fold"])
                                           for task in by_component.values()
                                           for rows in [[task_meta[min(task, key=lambda value: _hash_key("EAFF-H0A", value))]]]).items())),
    }


def run(args) -> dict:
    if not args.device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("E-AFF-H0A requires CUDA")
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(f"output already exists: {output}")
    output.mkdir(parents=True)
    input_root = Path(args.input)
    all_rows = list(_read_jsonl(input_root / "rows.label_blind.jsonl"))
    selected, selection_audit = select_h0a(all_rows)
    if selection_audit["closure_components"] < 16:
        raise ValueError("fewer than 16 independent H0A components")
    _write_jsonl(output / "selection.jsonl", selected)
    labels, label_audit = load_affinity(all_rows, Path(args.canonical_rows))
    observations = make_observations(all_rows, labels)
    lookup = {(row["task_id"], row["ligand_state_key"]): index
              for index, row in enumerate(observations)}
    selected_obs = []
    for row in selected:
        value = dict(observations[lookup[(row["task_id"], row["ligand_state_key"])]] )
        value["h0_partition"] = row["h0_partition"]
        selected_obs.append(value)

    pooled_keys = {row["ligand_state_key"] for row in observations}
    local_keys = {row["ligand_state_key"] for row in selected_obs}
    pooled, local_ligands, _ = _load_ligand_cache(Path(args.cache), pooled_keys, local_keys)
    ligand_prediction, ligand_diagnostics = ligand_oof(observations, pooled)
    baseline = np.asarray([ligand_prediction[lookup[(row["task_id"], row["ligand_state_key"])]]
                           for row in selected_obs])
    wrong_map, wrong_audit = build_derangement(
        selected_obs, input_root, Path(args.governance))
    _write_jsonl(output / "derangement.jsonl", wrong_audit)
    protein_keys = ({row["protein_sequence_sha256"] for row in selected_obs}
                    | set(wrong_map.values()))
    proteins = _load_protein_cache(Path(args.cache), protein_keys)
    bridge = _load_bridge(Path(args.checkpoint), args.device)
    with np.load(args.tbasis_values, allow_pickle=False) as tbasis:
        correct, deranged = _geometry_basis(
            selected_obs, bridge, proteins, local_ligands, wrong_map,
            tbasis["bin_rbf_expectation"].astype(np.float64),
            tbasis["calibration_coef"].astype(np.float64),
            tbasis["calibration_intercept"].astype(np.float64), args.device)
    null = coupling_null(correct)
    y = np.asarray([row["p_affinity"] for row in selected_obs])
    residual = y - baseline
    flat = [value.reshape(len(value), -1) for value in (correct, deranged, null)]
    energy = [np.full(len(y), np.nan) for _ in range(3)]
    fit_diagnostics = []
    direction_tasks, direction_values = [], []
    for task in sorted({row["task_id"] for row in selected_obs}):
        indices = np.asarray([index for index, row in enumerate(selected_obs)
                              if row["task_id"] == task])
        fit = indices[[selected_obs[index]["h0_partition"] == "fit" for index in indices]]
        test = indices[[selected_obs[index]["h0_partition"] == "test" for index in indices]]
        direction, diagnostic = fit_pair_ridge(
            flat[0][fit], residual[fit], np.asarray([task] * len(fit)), RIDGE_ALPHA)
        direction_tasks.append(task)
        direction_values.append(direction)
        for destination, features in zip(energy, flat):
            destination[test] = features[test] @ direction
        fit_diagnostics.append({"task_id": task, "fit": len(fit), "test": len(test), **diagnostic})
    test_indices = np.asarray([index for index, row in enumerate(selected_obs)
                               if row["h0_partition"] == "test"])
    if any(np.any(~np.isfinite(value[test_indices])) for value in energy):
        raise ValueError("H0A left non-finite test predictions")
    test_rows = [selected_obs[index] for index in test_indices]
    predictions = {
        "ligand": baseline[test_indices],
        "correct": baseline[test_indices] + energy[0][test_indices],
        "deranged": baseline[test_indices] + energy[1][test_indices],
        "null": baseline[test_indices] + energy[2][test_indices],
    }
    per_task = task_metrics(test_rows, y[test_indices], predictions)
    primary = component_macro_contrasts(per_task)
    confidence = component_bootstrap(per_task)
    headroom = primary["correct_minus_ligand"] >= 0.03 \
        and confidence["correct_minus_ligand"][0] > 0
    partner = primary["correct_minus_deranged"] >= 0.03 \
        and confidence["correct_minus_deranged"][0] > 0
    coupling = primary["correct_minus_null"] > 0 \
        and confidence["correct_minus_null"][0] > 0
    if headroom and partner:
        verdict = "TASK_LOCAL_RADIAL_HEADROOM_AND_PARTNER_SPECIFICITY_OBSERVED"
    elif headroom:
        verdict = "TASK_LOCAL_RADIAL_HEADROOM_WITHOUT_PARTNER_SPECIFICITY"
    else:
        verdict = "TASK_LOCAL_RADIAL_HEADROOM_NOT_OBSERVED"
    _write_jsonl(output / "task_metrics.jsonl", per_task)
    np.savez_compressed(output / "task_directions.npz",
                        task=np.asarray(direction_tasks),
                        direction=np.stack(direction_values).astype(np.float64))
    np.savez_compressed(output / "features.npz", correct=correct.astype(np.float32),
                        deranged=deranged.astype(np.float32), null=null.astype(np.float32),
                        ligand_oof=baseline.astype(np.float32), label=y.astype(np.float32),
                        partition=np.asarray([row["h0_partition"] for row in selected_obs]))
    endpoint = {}
    for name in ("Ki", "Kd"):
        subset = [row for row in per_task if row["endpoint_family"] == name]
        endpoint[name] = component_macro_contrasts(subset) if subset else {"components": 0}
    report = {
        "schema": "MetaSieve.EAffH0A.v1", "stage": STAGE,
        "created_utc": datetime.now(timezone.utc).isoformat(), "verdict": verdict,
        "research_only": True, "affinity_labels_read": True,
        "davis_label_reads": 0, "recipient_label_reads": 0,
        "selection": selection_audit,
        "label_audit": {key: value for key, value in label_audit.items() if key != "task_documents"},
        "ligand_prior": {"alpha": RIDGE_ALPHA, "folds": ligand_diagnostics},
        "task_heads": {"alpha": RIDGE_ALPHA, "fit_rows": FIT,
                       "test_rows": LIGANDS - FIT, "fits": fit_diagnostics},
        "primary_component_macro": primary, "component_bootstrap_ci95": confidence,
        "conditions": {"headroom": headroom, "partner_specificity": partner,
                       "coupling_attribution": coupling},
        "endpoint_secondary": endpoint,
        "interpretation_limits": [
            "task-local headroom is not target-specific biology",
            "task heads are oracle diagnostics, not deployable few-shot sections",
            "cross-assay target transport was not tested",
        ],
    }
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
            "preregistration": sha256_file(Path(__file__).with_name("EAFF_H0A_PREREGISTRATION.md")),
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
    parser.add_argument("--output", default="research/e0_identifiability/artifacts/eaff_h0a_v1")
    parser.add_argument("--device", default="cuda")
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), indent=2, sort_keys=True))
