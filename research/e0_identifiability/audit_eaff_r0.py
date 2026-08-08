"""Independently audit the E-AFF-R0 readout-scope artifacts.

Independence: concordance is recomputed here from an explicit pairwise loop that
shares no code with `metrics.concordance`, so an invariance claim cannot be
produced by a defect in the audited implementation.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from research.e0_identifiability.metrics import concordance
from scripts.source_affinity.common import sha256_file


ROOT = Path("research/e0_identifiability/artifacts/eaff_r0_v1")
H0C = Path("research/e0_identifiability/artifacts/eaff_h0c_v1_run2")
FORBIDDEN_KEYS = {
    "p_affinity", "standard_value", "published_value", "pchembl_value_reported",
    "activity_value", "label", "y",
}
SEED = 77010203


def _read_jsonl(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _keys(value) -> set[str]:
    if isinstance(value, dict):
        return set(value) | {key for item in value.values() for key in _keys(item)}
    if isinstance(value, list):
        return {key for item in value for key in _keys(item)}
    return set()


def independent_concordance(labels, predictions) -> float:
    """Explicit pairwise definition, deliberately not the audited implementation."""
    labels = np.asarray(labels, dtype=np.float64)
    predictions = np.asarray(predictions, dtype=np.float64)
    total = credit = 0.0
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            if labels[i] == labels[j]:
                continue
            total += 1.0
            higher = i if labels[i] > labels[j] else j
            lower = j if higher == i else i
            if predictions[higher] > predictions[lower]:
                credit += 1.0
            elif predictions[higher] == predictions[lower]:
                credit += 0.5
    return 0.5 if total == 0.0 else credit / total


def run(root: Path = ROOT) -> dict:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    report = json.loads((root / "report.json").read_text(encoding="utf-8"))
    hashes_match = all(sha256_file(root / name) == digest
                       for name, digest in manifest["outputs"].items())
    inputs_match = (
        sha256_file(H0C / "task_metrics.jsonl") == manifest["inputs"]["h0c_task_metrics"]
        and sha256_file(H0C / "report.json") == manifest["inputs"]["h0c_report"]
        and sha256_file(Path(__file__).with_name("metrics.py"))
        == manifest["inputs"]["metrics_module"])
    forbidden_artifact = bool(_keys(report) & FORBIDDEN_KEYS)

    rng = np.random.default_rng(SEED)
    size = report["test_ligands_per_task"]

    # 1. the audited implementation agrees with the independent definition
    agreement = 0.0
    for _ in range(64):
        labels = rng.normal(size=size)
        predictions = rng.normal(size=size)
        agreement = max(agreement, abs(
            concordance(labels, predictions) - independent_concordance(labels, predictions)))

    # 2. invariance holds under the independent definition too
    invariance = 0.0
    for _ in range(64):
        labels = rng.normal(size=size)
        predictions = rng.normal(size=size)
        base = independent_concordance(labels, predictions)
        shift = float(rng.normal() * 10.0)
        scale = float(abs(rng.normal()) + 0.1)
        for transformed in (
            independent_concordance(labels, predictions + shift),
            independent_concordance(labels, predictions * scale),
            independent_concordance(labels + shift, predictions),
            independent_concordance(labels * scale, predictions),
        ):
            invariance = max(invariance, abs(transformed - base))

    # 3. a constant within-task prediction scores exactly chance
    constant_scores = [independent_concordance(rng.normal(size=size), np.full(size, value))
                       for value in rng.normal(size=32)]
    constant_exactly_half = all(value == 0.5 for value in constant_scores)

    # 4. published H0C contrasts recomputed from the per-task record
    per_task = _read_jsonl(H0C / "task_metrics.jsonl")
    local = np.asarray([row["local_ligand"] for row in per_task])
    correct = np.asarray([row["correct"] for row in per_task])
    deranged = np.asarray([row["deranged"] for row in per_task])
    recomputed = {
        "correct_minus_local": float((correct - local).mean()),
        "deranged_minus_local": float((deranged - local).mean()),
        "correct_minus_deranged": float((correct - deranged).mean()),
    }
    published = report["published_geometry_effect"]["contrasts"]
    contrast_error = max(abs(recomputed[name] - published[name]["mean"]) for name in recomputed)

    level_all_half = all(
        row["within_task_concordance"]["level_oracle"] == 0.5
        for row in report["location_credit"])

    checks = {
        "preexisting_hashes_match": hashes_match,
        "declared_inputs_match": inputs_match,
        "artifact_affinity_fields_absent": not forbidden_artifact,
        "report_declares_no_affinity_label_read": report["affinity_labels_read"] is False,
        "no_davis_or_recipient_reads": report["davis_label_reads"] == 0
        and report["recipient_label_reads"] == 0,
        "no_training_performed": report["training_performed"] is False,
        "audited_metric_matches_independent_definition": agreement == 0.0,
        "invariance_holds_under_independent_definition": invariance == 0.0,
        "constant_prediction_scores_exactly_half": constant_exactly_half,
        "reported_invariance_is_exact": report["invariance"]["exactly_invariant"] is True
        and all(value == 0.0
                for value in report["invariance"]["max_absolute_deviation"].values()),
        "level_oracle_exactly_half_at_every_share": level_all_half,
        "published_contrasts_reproduce": contrast_error <= 5e-6,
    }
    result = {
        "schema": "MetaSieve.EAffR0PostrunAudit.v1",
        "stage": report["stage"],
        "verdict": "POSTRUN_AUDIT_PASS" if all(checks.values()) else "POSTRUN_AUDIT_FAIL",
        "checks": checks,
        "independent_definition_agreement_max_error": agreement,
        "independent_invariance_max_deviation": invariance,
        "published_contrast_max_error": contrast_error,
        "recomputed_h0c_contrasts": {name: round(value, 5)
                                     for name, value in recomputed.items()},
        "scientific_verdict": report["verdict"],
    }
    (root / "postrun_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (root / "POSTRUN_AUDIT_MANIFEST.json").write_text(
        json.dumps({
            "stage": report["stage"],
            "postrun_audit_sha256": sha256_file(root / "postrun_audit.json"),
            "auditor_sha256": sha256_file(Path(__file__)),
        }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
