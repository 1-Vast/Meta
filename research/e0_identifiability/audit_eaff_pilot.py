"""Independent post-run audit for E-AFF-P0 artifacts."""
from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

import numpy as np

from research.e0_identifiability.eaff_pilot_contract import (
    component_bootstrap,
    component_macro_contrasts,
)
from scripts.source_affinity.common import sha256_file


def _read_jsonl(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def audit(root: Path) -> dict:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    report = json.loads((root / "report.json").read_text(encoding="utf-8"))
    output_hashes = {name: sha256_file(root / name)
                     for name in manifest["outputs"]}
    hash_match = output_hashes == manifest["outputs"]
    selection = list(_read_jsonl(root / "selection.jsonl"))
    derangement = list(_read_jsonl(root / "derangement.jsonl"))
    per_task = list(_read_jsonl(root / "task_metrics.jsonl"))
    primary = component_macro_contrasts(per_task)
    confidence = component_bootstrap(per_task, seed=17)
    with np.load(root / "features.npz", allow_pickle=False) as values:
        correct = values["correct"].astype(np.float64)
        null = values["null"].astype(np.float64)
        feature_shape = list(correct.shape)
        finite = bool(all(np.all(np.isfinite(values[key])) for key in values.files))
        chemistry_error = float(np.max(np.abs(correct.sum(-1) - null.sum(-1))))
        radial_error = float(np.max(np.abs(
            correct.sum((-3, -2)) - null.sum((-3, -2)))))
    directions = np.load(root / "directions.npy", allow_pickle=False)
    checks = {
        "output_hashes_match": hash_match,
        "selection_rows_4900": len(selection) == 4900,
        "selection_tasks_245": len({row["task_id"] for row in selection}) == 245,
        "selection_components_245": len({row["closure_component_id"] for row in selection}) == 245,
        "derangement_rows_equal_unique_proteins": len(derangement) == len({
            row["correct_protein"] for row in derangement}),
        "derangement_wrong_reuse_zero": len(derangement) == len({
            row["wrong_protein"] for row in derangement}),
        "derangement_all_below_0_40": max(row["local_identity"] for row in derangement) < 0.40,
        "features_shape_4900_8_6_6": feature_shape == [4900, 8, 6, 6],
        "features_all_finite": finite,
        "directions_shape_5_288": list(directions.shape) == [5, 288],
        "directions_all_finite": bool(np.all(np.isfinite(directions))),
        "coupling_null_chemistry_marginal_error_le_1e_6": chemistry_error <= 1e-6,
        "coupling_null_radial_marginal_error_le_1e_6": radial_error <= 1e-6,
        "component_macro_exact": all(abs(primary[key] - report["primary_component_macro"][key]) <= 1e-12
                                     for key in primary),
        "bootstrap_exact": all(np.allclose(confidence[key], report["component_bootstrap_ci95"][key])
                               for key in confidence),
        "no_davis_or_recipient_reads": report["davis_label_reads"] == 0
        and report["recipient_label_reads"] == 0,
    }
    result = {
        "schema": "MetaSieve.EAffPilotPostrunAudit.v1",
        "stage": report["stage"],
        "verdict": "POSTRUN_AUDIT_PASS" if all(checks.values()) else "POSTRUN_AUDIT_FAIL_CLOSED",
        "checks": checks,
        "recomputed_primary": primary,
        "recomputed_ci95": confidence,
        "derangement": {
            "pairs": len(derangement),
            "maximum_identity": max(row["local_identity"] for row in derangement),
            "wrong_reuse": len(derangement) - len({row["wrong_protein"] for row in derangement}),
            "fold_counts": dict(sorted(Counter(row["fold"] for row in derangement).items())),
        },
        "coupling_null": {"chemistry_marginal_max_error": chemistry_error,
                          "radial_marginal_max_error": radial_error},
    }
    return result


if __name__ == "__main__":
    root = Path("research/e0_identifiability/artifacts/eaff_p0_v1")
    result = audit(root)
    path = root / "postrun_audit.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    audit_manifest = {
        "schema": "MetaSieve.EAffPilotPostrunAuditManifest.v1",
        "core_manifest_sha256": sha256_file(root / "manifest.json"),
        "postrun_audit_sha256": sha256_file(path),
        "audit_code_sha256": sha256_file(Path(__file__)),
    }
    (root / "POSTRUN_AUDIT_MANIFEST.json").write_text(
        json.dumps(audit_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    raise SystemExit(0 if result["verdict"] == "POSTRUN_AUDIT_PASS" else 2)
