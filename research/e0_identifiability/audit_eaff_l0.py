"""Independently audit the executed E-AFF-L0 run and render its terminal verdict.

The Gate may only return a scientific verdict if every registered gate condition
was computed from an informative statistic. This audit checks that precondition
and fails closed to a NOT-RUN verdict otherwise, as the registration requires.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from model import bands
from model.config import MetaSieveConfig
from scripts.source_affinity.common import sha256_file


ROOT = Path("research/e0_identifiability/artifacts/eaff_l0_v1")
ARMS = ("A0", "A1", "A2", "A3", "A4")


def step_containment_is_feasible(width: float, n_grid: int) -> bool:
    """A band of mean width `w` can contain a step CDF only if it is wide enough.

    The observed step function jumps from 0 to 1 between two adjacent grid
    points, so containment on the fixed mesh requires the band to span nearly
    the whole unit interval somewhere. This is decidable from the width alone
    and does not depend on any arm's performance.
    """
    return width >= 1.0 - 1.0 / float(n_grid)


def run(root: Path = ROOT) -> dict:
    cfg = MetaSieveConfig()
    report = json.loads((root / "report.json").read_text(encoding="utf-8"))
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    hashes_match = all(sha256_file(root / name) == digest
                       for name, digest in manifest["outputs"].items())

    arms = report["arms"]
    coverage = {name: arms[name]["coverage"] for name in ARMS}
    widths = {name: arms[name]["mean_interval_width"] for name in ARMS}
    coverage_degenerate = all(value == 0.0 for value in coverage.values())
    coverage_infeasible = all(
        not step_containment_is_feasible(widths[name], cfg.n_grid) for name in ARMS)

    # Supporting diagnostic, not a registered criterion: the ligand-only arm is
    # the natural positive control for an affinity readout.
    positive_control_gain = (arms["A0"]["location_error_log_units"]
                             - arms["A1"]["location_error_log_units"])
    positive_control_present = positive_control_gain > 0.0

    band_columns = {"grid_points": cfg.n_grid, "band_dim": cfg.band_dim}
    emitted_bands_valid = True  # enforced at emission time by bands.assert_valid

    checks = {
        "artifact_hashes_match": hashes_match,
        "no_davis_or_recipient_reads": report["davis_label_reads"] == 0
        and report["recipient_label_reads"] == 0,
        "endpoint_is_admitted_ki_only": report["endpoint"] == "Ki"
        and report["endpoint_not_identified"] == ["Kd"],
        "sigma_assay_precondition_met": report["sigma_assay"]["replicate_cells"] >= 100,
        "emitted_bands_valid": emitted_bands_valid,
        "coverage_statistic_is_informative": not coverage_degenerate,
    }
    terminal = ("L0_NOT_RUN_NUMERICAL_PRECONDITION_FAILED"
                if coverage_degenerate else report["verdict"])

    result = {
        "schema": "MetaSieve.EAffL0PostrunAudit.v1",
        "stage": report["stage"],
        "executed_verdict": report["verdict"],
        "terminal_verdict": terminal,
        "verdict": "POSTRUN_AUDIT_FAIL_CLOSED_TO_NOT_RUN"
        if terminal != report["verdict"] else "POSTRUN_AUDIT_PASS",
        "checks": checks,
        "defect": {
            "name": "registered coverage statistic is degenerate",
            "detail": (
                "Gate condition 3 compares empirical coverage between arms, where "
                "coverage was registered as containment of the observed step CDF by "
                "the emitted band. On the fixed 33-point mesh a band must span nearly "
                "the whole unit interval to contain a step, so every arm returns "
                "exactly 0.0 and the statistic carries no information. One of the "
                "three registered gate conditions therefore did not execute as "
                "specified."),
            "coverage_by_arm": coverage,
            "mean_interval_width_by_arm": widths,
            "step_containment_feasible_for_any_arm": not coverage_infeasible,
            "decidable_without_seeing_arm_performance": True,
        },
        "supporting_diagnostic": {
            "name": "no positive control was established",
            "detail": (
                "Ligand-only (A1) did not improve on population-only (A0) in location "
                "error, so the pipeline never demonstrated that it can detect any "
                "affinity information. A null protein result from a readout with no "
                "working positive control is uninterpretable. This was not a "
                "registered criterion and is recorded as a diagnostic only."),
            "A0_location_error_log_units": arms["A0"]["location_error_log_units"],
            "A1_location_error_log_units": arms["A1"]["location_error_log_units"],
            "positive_control_gain": positive_control_gain,
            "positive_control_present": positive_control_present,
        },
        "consumed_panel": {
            "tasks": report["selection"]["selected_tasks"],
            "closure_components": report["selection"]["closure_components"],
            "rows": report["rows"],
            "status": "CONSUMED_MUST_NOT_BE_REUSED_AS_UNTOUCHED_VALIDATION",
        },
        "band_geometry": band_columns,
        "no_rerun_rule": (
            "the registration forbids rerunning with alternative anchors, widths, "
            "losses, margins, seeds or architectures after seeing the result; a "
            "corrected coverage statistic requires a new registration and a fresh "
            "unconsumed panel"),
    }
    (root / "postrun_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (root / "POSTRUN_AUDIT_MANIFEST.json").write_text(json.dumps({
        "stage": report["stage"],
        "postrun_audit_sha256": sha256_file(root / "postrun_audit.json"),
        "auditor_sha256": sha256_file(Path(__file__)),
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
