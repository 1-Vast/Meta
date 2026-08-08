"""Independently audit the E-AFF-X0-B crossed-design artifacts.

Independence: every packed rectangle is checked against the observed panel for
cell existence and cell-disjointness without trusting the construction; cluster,
target and target-pair counts are recounted from the X0 cell record; the design
effect and effective sample size are recomputed from their definitions; and the
breakeven `rho` is found by bisection rather than by the runner's closed form.

Expected values are asserted against the registered constants so a silent drift
in the packing or the aggregation fails the audit closed.
"""
from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path

from research.e0_identifiability.run_eaff_x0b import pack_cell_disjoint
from scripts.source_affinity.common import sha256_file


ROOT = Path("research/e0_identifiability/artifacts/eaff_x0b_v1")
X0 = Path("research/e0_identifiability/artifacts/eaff_x0_v1")
FEAS = Path("research/e0_identifiability/artifacts/eaff_x0_feas_v1")
FORBIDDEN_KEYS = {
    "p_affinity", "standard_value", "published_value", "pchembl_value_reported",
    "activity_value", "label", "y",
}
REQUIRED = 245
EXPECTED = {
    "Ki": {"cell_disjoint_units": 11168, "clusters": 36,
           "distinct_target_pairs": 205, "breakeven_rho_star": 0.0915},
    "Kd": {"cell_disjoint_units": 1041, "clusters": 12,
           "distinct_target_pairs": 49, "breakeven_rho_star": 0.0164},
}
CAP_GRID = (1, 2, 3, 5, 8, 12, 20, 32, 50, 80, 125, 200, 320, 500, 800, 1250, 2000, None)


def _read_jsonl(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _keys(value) -> set[str]:
    if isinstance(value, dict):
        return set(value) | {key for item in value.values() for key in _keys(item)}
    if isinstance(value, list):
        return {key for item in value for key in _keys(item)}
    return set()


def design(sizes: list[int], cap: int | None) -> tuple[int, float]:
    limited = [min(size, cap) if cap else size for size in sizes]
    total = sum(limited)
    if total == 0:
        return 0, 0.0
    return total, sum(size * size for size in limited) / total


def n_eff(total: int, influence: float, rho: float) -> float:
    return total / (1.0 + (influence - 1.0) * rho)


def breakeven_by_bisection(sizes: list[int]) -> tuple[float, int | None]:
    """Largest rho reaching the requirement, found by bisection over each cap."""
    best, best_cap = 0.0, None
    for cap in CAP_GRID:
        total, influence = design(sizes, cap)
        if total < REQUIRED:
            continue
        if influence <= 1.0:
            candidate = 1.0
        else:
            low, high = 0.0, 1.0
            if n_eff(total, influence, high) >= REQUIRED:
                candidate = 1.0
            else:
                for _ in range(200):
                    mid = 0.5 * (low + high)
                    if n_eff(total, influence, mid) >= REQUIRED:
                        low = mid
                    else:
                        high = mid
                candidate = low
        if candidate > best:
            best, best_cap = candidate, cap
    return best, best_cap


def run(root: Path = ROOT) -> dict:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    report = json.loads((root / "report.json").read_text(encoding="utf-8"))
    x0_report = json.loads((X0 / "report.json").read_text(encoding="utf-8"))
    hashes_match = all(sha256_file(root / name) == digest
                       for name, digest in manifest["outputs"].items())
    inputs_match = all(
        sha256_file(path) == manifest["inputs"][key]
        for key, path in (
            ("x0_cells", X0 / "cells.jsonl"),
            ("x0_report", X0 / "report.json"),
            ("x0_manifest", X0 / "manifest.json"),
            ("feasibility_report", FEAS / "report.json"),
        ))
    forbidden_artifact = bool(_keys(report) & FORBIDDEN_KEYS)

    panel_targets: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    panel_endpoint: dict[str, str] = {}
    panel_cluster: dict[str, str] = {}
    for cell in _read_jsonl(X0 / "cells.jsonl"):
        panel = cell["panel_id"]
        panel_endpoint[panel] = cell["endpoint_family"]
        panel_cluster[panel] = cell["closure_component_id"]
        panel_targets[panel][cell["protein_sequence_sha256"]].add(
            cell["ligand_connectivity_key"])

    endpoints = {}
    packing_valid = True
    for endpoint in ("Ki", "Kd"):
        sizes_by_cluster: dict[str, int] = defaultdict(int)
        target_pairs: set[tuple[str, str]] = set()
        targets: set[str] = set()
        ligands: set[str] = set()
        for panel, family in panel_endpoint.items():
            if family != endpoint:
                continue
            observed = panel_targets[panel]
            packed = pack_cell_disjoint(observed)
            used: set[tuple[str, str]] = set()
            for left, right, first, second in packed:
                cells = {(left, first), (left, second), (right, first), (right, second)}
                # every cell must exist in the observed panel and be unused
                if any(ligand not in observed.get(target, set()) for target, ligand in cells):
                    packing_valid = False
                if used & cells:
                    packing_valid = False
                if left == right or first == second:
                    packing_valid = False
                used |= cells
                target_pairs.add((left, right))
                targets.update((left, right))
                ligands.update((first, second))
            if packed:
                sizes_by_cluster[panel_cluster[panel]] += len(packed)

        sizes = sorted(sizes_by_cluster.values(), reverse=True)
        total = sum(sizes)
        rho_star, rho_star_cap = breakeven_by_bisection(sizes)
        endpoints[endpoint] = {
            "cell_disjoint_units": total,
            "clusters": len(sizes),
            "distinct_target_pairs": len(target_pairs),
            "distinct_targets": len(targets),
            "distinct_ligands": len(ligands),
            "breakeven_rho_star": rho_star,
            "breakeven_rho_star_cap": rho_star_cap,
            "cluster_bound_rho_max": len(sizes) / REQUIRED,
            "n_eff_at_rho_one_best_cap": max(
                n_eff(*design(sizes, cap), 1.0) for cap in CAP_GRID if design(sizes, cap)[0] > 0),
        }

    published = report["endpoints"]
    count_error = max(
        abs(endpoints[endpoint][key] - published[endpoint][key])
        for endpoint in ("Ki", "Kd")
        for key in ("cell_disjoint_units", "clusters", "distinct_target_pairs",
                    "distinct_targets", "distinct_ligands"))
    rho_error = max(abs(endpoints[endpoint]["breakeven_rho_star"]
                        - published[endpoint]["breakeven_rho_star"])
                    for endpoint in ("Ki", "Kd"))
    bound_error = max(abs(endpoints[endpoint]["cluster_bound_rho_max"]
                          - published[endpoint]["cluster_bound_rho_max"])
                      for endpoint in ("Ki", "Kd"))
    expected_error = max(
        max(abs(endpoints[endpoint][key] - value)
            for key, value in EXPECTED[endpoint].items() if key != "breakeven_rho_star")
        for endpoint in ("Ki", "Kd"))
    expected_rho_error = max(
        abs(endpoints[endpoint]["breakeven_rho_star"] - EXPECTED[endpoint]["breakeven_rho_star"])
        for endpoint in ("Ki", "Kd"))

    rho_one_matches_x0 = all(
        abs(endpoints[endpoint]["n_eff_at_rho_one_best_cap"]
            - x0_report["endpoints"][endpoint]["dependency_components"]) < 1e-9
        for endpoint in ("Ki", "Kd"))

    frozen = report["frozen_unchanged"]
    checks = {
        "preexisting_hashes_match": hashes_match,
        "declared_inputs_match": inputs_match,
        "artifact_affinity_fields_absent": not forbidden_artifact,
        "report_zero_affinity_value_fields": report["affinity_value_fields_selected"] == 0
        and report["affinity_values_read"] is False,
        "no_davis_or_recipient_reads": report["davis_label_reads"] == 0
        and report["recipient_label_reads"] == 0,
        "no_training_performed": report["training_performed"] is False,
        "source_x0_firewall_clean": x0_report["affinity_value_fields_selected"] == 0,
        "packing_cells_exist_and_are_disjoint": packing_valid,
        "unit_and_cluster_counts_reproduce": count_error == 0,
        "counts_match_registered_expectations": expected_error == 0,
        "breakeven_rho_reproduces_by_bisection": rho_error <= 5e-5,
        "breakeven_rho_matches_registered_expectations": expected_rho_error <= 5e-5,
        "cluster_bound_reproduces": bound_error <= 5e-5,
        "rho_one_reduces_to_x0_components": rho_one_matches_x0,
        "frozen_design_constants_unchanged": (
            frozen["interaction_rms_over_assay_noise"] == 0.5
            and frozen["variance_ratio"] == 1.25
            and frozen["alpha_one_sided"] == 0.05
            and frozen["power"] == 0.80
            and frozen["required_effective_n"] == REQUIRED),
        "verdict_is_conditional_not_evidence_of_interaction": all(
            "CONDITIONAL_DESIGN_SUPPORTED" in published[endpoint]["verdict"]
            for endpoint in ("Ki", "Kd")),
    }
    result = {
        "schema": "MetaSieve.EAffX0BPostrunAudit.v1",
        "stage": report["stage"],
        "verdict": "POSTRUN_AUDIT_PASS" if all(checks.values()) else "POSTRUN_AUDIT_FAIL",
        "checks": checks,
        "recomputed_endpoints": {
            endpoint: {key: (round(value, 6) if isinstance(value, float) else value)
                       for key, value in values.items()}
            for endpoint, values in endpoints.items()},
        "count_max_error": count_error,
        "breakeven_rho_max_error": rho_error,
        "registered_expectation_max_error": expected_error,
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
