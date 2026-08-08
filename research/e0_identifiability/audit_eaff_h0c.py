"""Independently audit E-AFF-H0C artifacts without refitting models."""
from __future__ import annotations

from collections import Counter, defaultdict
import json
from pathlib import Path

import numpy as np

from research.e0_identifiability.eaff_h0c_contract import (
    centered_interaction,
    component_bootstrap,
    component_summary,
)
from research.e0_identifiability.metrics import concordance
from scripts.source_affinity.common import sha256_file


ROOT = Path("research/e0_identifiability/artifacts/eaff_h0c_v1_run2")


def _read_jsonl(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def run(root: Path = ROOT) -> dict:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    report = json.loads((root / "report.json").read_text(encoding="utf-8"))
    selection = _read_jsonl(root / "selection.jsonl")
    derangement = _read_jsonl(root / "derangement.jsonl")
    recorded_tasks = _read_jsonl(root / "task_metrics.jsonl")
    old_tasks = {row["task_id"] for row in _read_jsonl(
        Path("research/e0_identifiability/artifacts/eaff_h0a_v1/selection.jsonl"))}
    hashes_match = all(sha256_file(root / name) == digest
                       for name, digest in manifest["outputs"].items())

    grouped = defaultdict(lambda: defaultdict(list))
    for row in selection:
        grouped[row["task_id"]][row["h0c_partition"]].append(row)
    partition_ok = all(len(parts["support"]) == 20 and len(parts["test"]) == 20
                       for parts in grouped.values())
    scaffold_ok = all(not ({row["murcko_scaffold"] for row in parts["support"]}
                           & {row["murcko_scaffold"] for row in parts["test"]})
                      for parts in grouped.values())
    wrong_counts = Counter(row["wrong_protein"] for row in derangement)
    wrong_ok = (len(derangement) == len(grouped)
                and max(row["local_identity"] for row in derangement) < 0.40
                and max(wrong_counts.values()) == 1)

    with np.load(root / "features.npz", allow_pickle=False) as data:
        correct_phi = data["correct_phi"].astype(np.float64)
        deranged_phi = data["deranged_phi"].astype(np.float64)
        global_ligand = data["global_ligand"].astype(np.float64)
        local_effect = data["local_effect"].astype(np.float64)
        correct_effect = data["correct_effect"].astype(np.float64)
        deranged_effect = data["deranged_effect"].astype(np.float64)
        labels = data["label"].astype(np.float64)
        partition = data["partition"].astype(str)
    correct_psi = centered_interaction(correct_phi)
    deranged_psi = centered_interaction(deranged_phi)
    marginal_error = max(
        float(np.abs(correct_psi.sum(-1)).max()),
        float(np.abs(correct_psi.sum((-3, -2))).max()),
        float(np.abs(deranged_psi.sum(-1)).max()),
        float(np.abs(deranged_psi.sum((-3, -2))).max()),
    )
    test = np.flatnonzero(partition == "test")
    test_rows = [selection[index] for index in test]
    scores = {
        "global_ligand": global_ligand[test],
        "local_ligand": global_ligand[test] + local_effect[test],
        "correct": global_ligand[test] + local_effect[test] + correct_effect[test],
        "deranged": global_ligand[test] + local_effect[test] + deranged_effect[test],
    }
    by_task = defaultdict(list)
    for index, row in enumerate(test_rows):
        by_task[row["task_id"]].append(index)
    recomputed = []
    for task, indices in sorted(by_task.items()):
        first = test_rows[indices[0]]
        recomputed.append({
            "task_id": task, "closure_component_id": first["closure_component_id"],
            "endpoint_family": first["endpoint_family"],
            **{name: concordance(labels[test][indices], values[indices])
               for name, values in scores.items()},
        })
    task_error = max(abs(row[key] - other[key])
                     for row, other in zip(recomputed, recorded_tasks)
                     for key in scores)
    primary = component_summary(recomputed)
    confidence = component_bootstrap(recomputed)
    primary_error = max(abs(primary[key] - report["primary_component_macro"][key])
                        for key in primary)
    bootstrap_error = max(abs(confidence[key][index]
                              - report["component_bootstrap_ci95"][key][index])
                          for key in confidence for index in (0, 1))
    checks = {
        "preexisting_hashes_match": hashes_match,
        "selection_54_tasks_2160_rows": len(grouped) == 54 and len(selection) == 2160,
        "h0a_task_overlap_zero": not (set(grouped) & old_tasks),
        "partition_20_20_each_task": partition_ok,
        "scaffold_overlap_zero": scaffold_ok,
        "derangement_one_to_one_below_0_40": wrong_ok,
        "centered_marginals_le_1e_6": marginal_error <= 1e-6,
        "task_metrics_match": task_error <= 1e-7,
        "component_macro_matches": primary_error <= 1e-7,
        "bootstrap_matches": bootstrap_error <= 1e-7,
        "no_davis_or_recipient_reads": report["davis_label_reads"] == 0
        and report["recipient_label_reads"] == 0,
    }
    result = {
        "schema": "MetaSieve.EAffH0CPostrunAudit.v1",
        "stage": report["stage"],
        "verdict": "POSTRUN_AUDIT_PASS" if all(checks.values()) else "POSTRUN_AUDIT_FAIL",
        "checks": checks,
        "centered_marginal_max_error": marginal_error,
        "task_metric_max_error": task_error,
        "primary_metric_max_error": primary_error,
        "bootstrap_max_error": bootstrap_error,
        "recomputed_primary": primary,
        "recomputed_ci95": confidence,
    }
    (root / "postrun_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    audit_manifest = {
        "stage": report["stage"],
        "postrun_audit_sha256": sha256_file(root / "postrun_audit.json"),
        "auditor_sha256": sha256_file(Path(__file__)),
    }
    (root / "POSTRUN_AUDIT_MANIFEST.json").write_text(
        json.dumps(audit_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
