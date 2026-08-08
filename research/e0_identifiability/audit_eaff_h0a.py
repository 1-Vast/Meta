"""Independent reconstruction audit for E-AFF-H0A."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from research.e0_identifiability.eaff_pilot_contract import (
    component_bootstrap,
    component_macro_contrasts,
    fit_pair_ridge,
    task_metrics,
)
from scripts.source_affinity.common import sha256_file


def _read_jsonl(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def audit(root: Path) -> dict:
    report = json.loads((root / "report.json").read_text(encoding="utf-8"))
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    preexisting_hashes = all(sha256_file(root / name) == digest
                             for name, digest in manifest["outputs"].items())
    rows = list(_read_jsonl(root / "selection.jsonl"))
    derangement = list(_read_jsonl(root / "derangement.jsonl"))
    stored_metrics = list(_read_jsonl(root / "task_metrics.jsonl"))
    with np.load(root / "features.npz", allow_pickle=False) as values:
        correct = values["correct"].astype(np.float64)
        deranged = values["deranged"].astype(np.float64)
        null = values["null"].astype(np.float64)
        baseline = values["ligand_oof"].astype(np.float64)
        labels = values["label"].astype(np.float64)
        partition = values["partition"].astype(str)
    flat = [value.reshape(len(value), -1) for value in (correct, deranged, null)]
    residual = labels - baseline
    energy = [np.full(len(labels), np.nan) for _ in range(3)]
    direction_tasks, directions = [], []
    for task in sorted({row["task_id"] for row in rows}):
        indices = np.asarray([index for index, row in enumerate(rows) if row["task_id"] == task])
        fit, test = indices[partition[indices] == "fit"], indices[partition[indices] == "test"]
        direction, _ = fit_pair_ridge(flat[0][fit], residual[fit],
                                      np.asarray([task] * len(fit)), alpha=10.0)
        direction_tasks.append(task)
        directions.append(direction)
        for destination, features in zip(energy, flat):
            destination[test] = features[test] @ direction
    np.savez_compressed(root / "task_directions.npz", task=np.asarray(direction_tasks),
                        direction=np.stack(directions).astype(np.float64))
    test = np.flatnonzero(partition == "test")
    test_rows = [rows[index] for index in test]
    predictions = {"ligand": baseline[test], "correct": baseline[test] + energy[0][test],
                   "deranged": baseline[test] + energy[1][test],
                   "null": baseline[test] + energy[2][test]}
    recomputed_metrics = task_metrics(test_rows, labels[test], predictions)
    primary = component_macro_contrasts(recomputed_metrics)
    confidence = component_bootstrap(recomputed_metrics)
    metric_error = max(abs(left[key] - right[key])
                       for left, right in zip(recomputed_metrics, stored_metrics)
                       for key in ("ligand", "correct", "deranged", "null"))
    chemistry_error = float(np.max(np.abs(correct.sum(-1) - null.sum(-1))))
    radial_error = float(np.max(np.abs(correct.sum((-3, -2)) - null.sum((-3, -2)))))
    checks = {
        "preexisting_hashes_match": preexisting_hashes,
        "selection_107_tasks_4280_rows": len(rows) == 4280
        and len({row["task_id"] for row in rows}) == 107,
        "partition_20_20_each_task": all(
            sum(row["task_id"] == task and row["h0_partition"] == part for row in rows) == 20
            for task in {row["task_id"] for row in rows} for part in ("fit", "test")),
        "derangement_one_to_one_below_0_40": len(derangement) == 107
        and len({row["wrong_protein"] for row in derangement}) == 107
        and max(row["local_identity"] for row in derangement) < 0.40,
        "task_metric_max_error_le_1e_7": metric_error <= 1e-7,
        "component_macro_matches": all(abs(primary[key] - report["primary_component_macro"][key]) <= 1e-7
                                       for key in primary),
        "bootstrap_matches": all(np.allclose(confidence[key], report["component_bootstrap_ci95"][key],
                                              atol=1e-7) for key in confidence),
        "coupling_null_marginals_le_1e_6": max(chemistry_error, radial_error) <= 1e-6,
        "no_davis_or_recipient_reads": report["davis_label_reads"] == 0
        and report["recipient_label_reads"] == 0,
    }
    manifest["outputs"]["task_directions.npz"] = sha256_file(root / "task_directions.npz")
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result = {
        "schema": "MetaSieve.EAffH0APostrunAudit.v1", "stage": report["stage"],
        "verdict": "POSTRUN_AUDIT_PASS" if all(checks.values()) else "POSTRUN_AUDIT_FAIL_CLOSED",
        "checks": checks, "recomputed_primary": primary, "recomputed_ci95": confidence,
        "task_metric_max_error": metric_error,
        "derangement_max_identity": max(row["local_identity"] for row in derangement),
        "coupling_null_marginal_max_error": max(chemistry_error, radial_error),
    }
    audit_path = root / "postrun_audit.json"
    audit_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    postrun_manifest = {
        "schema": "MetaSieve.EAffH0APostrunAuditManifest.v1",
        "core_manifest_sha256": sha256_file(manifest_path),
        "postrun_audit_sha256": sha256_file(audit_path),
        "audit_code_sha256": sha256_file(Path(__file__)),
    }
    (root / "POSTRUN_AUDIT_MANIFEST.json").write_text(
        json.dumps(postrun_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    result = audit(Path("research/e0_identifiability/artifacts/eaff_h0a_v1"))
    print(json.dumps(result, indent=2, sort_keys=True))
    raise SystemExit(0 if result["verdict"] == "POSTRUN_AUDIT_PASS" else 2)
