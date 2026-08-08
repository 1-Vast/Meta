"""Independently audit the label-blind E-AFF-X0 census artifacts."""
from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path

from research.e0_identifiability.run_eaff_x0 import _dependency_components
from research.e0_identifiability.x0_contract import rectangle_count
from scripts.source_affinity.common import sha256_file


ROOT = Path("research/e0_identifiability/artifacts/eaff_x0_v1")
FORBIDDEN_KEYS = {
    "p_affinity", "standard_value", "published_value", "pchembl_value_reported",
    "activity_value", "label",
}


def _read_jsonl(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _keys(value) -> set[str]:
    if isinstance(value, dict):
        return set(value) | {key for item in value.values() for key in _keys(item)}
    if isinstance(value, list):
        return {key for item in value for key in _keys(item)}
    return set()


def run(root: Path = ROOT) -> dict:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    report = json.loads((root / "report.json").read_text(encoding="utf-8"))
    panels = _read_jsonl(root / "panels.jsonl")
    cells = _read_jsonl(root / "cells.jsonl")
    recorded_dependencies = _read_jsonl(root / "dependency_components.jsonl")
    hashes_match = all(sha256_file(root / name) == digest
                       for name, digest in manifest["outputs"].items())
    sql = Path("research/e0_identifiability/x0_metadata.sql").read_text(encoding="utf-8").lower()
    forbidden_sql = any(value in sql for value in
                        ("act.standard_value", "act.value", "pchembl_value"))
    forbidden_artifact = bool((_keys(panels) | _keys(cells)) & FORBIDDEN_KEYS)

    cells_by_panel = defaultdict(list)
    for cell in cells:
        cells_by_panel[cell["panel_id"]].append(cell)
    panel_count_error = rectangle_error = 0
    for panel in panels:
        target_ligands = defaultdict(set)
        replicate_target_ligands = defaultdict(set)
        for cell in cells_by_panel[panel["panel_id"]]:
            target_ligands[cell["protein_sequence_sha256"]].add(
                cell["ligand_connectivity_key"])
            if cell["exact_assay_replicate_supported"]:
                replicate_target_ligands[cell["protein_sequence_sha256"]].add(
                    cell["ligand_connectivity_key"])
        pairs, rectangles = rectangle_count(target_ligands)
        rep_pairs, rep_rectangles = rectangle_count(replicate_target_ligands)
        panel_count_error = max(panel_count_error,
                                abs(len(cells_by_panel[panel["panel_id"]]) - panel["cell_count"]))
        rectangle_error = max(rectangle_error, abs(rectangles - panel["rectangles"]),
                              abs(pairs - panel["target_pairs"]),
                              abs(rep_rectangles - panel["replicate_rectangles"]),
                              abs(rep_pairs - panel["replicate_target_pairs"]))

    recomputed_dependencies = []
    recomputed_endpoint = {}
    for endpoint in ("Ki", "Kd"):
        total = _dependency_components(panels, endpoint, False)
        replicated = _dependency_components(panels, endpoint, True)
        recomputed_dependencies.extend(total + replicated)
        endpoint_panels = [row for row in panels if row["endpoint_family"] == endpoint]
        endpoint_cells = [row for row in cells if row["endpoint_family"] == endpoint]
        total_rectangles = sum(row["rectangles"] for row in endpoint_panels)
        largest = max((row["rectangle_count"] for row in total), default=0)
        recomputed_endpoint[endpoint] = {
            "panels": len(endpoint_panels),
            "dependency_components": len(total),
            "replicate_supported_dependency_components": len(replicated),
            "targets": len({row["protein_sequence_sha256"] for row in endpoint_cells}),
            "ligands": len({row["ligand_connectivity_key"] for row in endpoint_cells}),
            "cells": len(endpoint_cells),
            "target_pairs": sum(row["target_pairs"] for row in endpoint_panels),
            "rectangles": total_rectangles,
            "replicate_supported_rectangles": sum(
                row["replicate_rectangles"] for row in endpoint_panels),
            "largest_dependency_rectangle_fraction": (
                largest / total_rectangles if total_rectangles else 0.0),
        }
    dependency_match = recomputed_dependencies == recorded_dependencies
    endpoint_error = max(
        abs(recomputed_endpoint[name][key] - report["endpoints"][name][key])
        for name in ("Ki", "Kd") for key in recomputed_endpoint[name])
    checks = {
        "preexisting_hashes_match": hashes_match,
        "sql_affinity_projection_absent": not forbidden_sql,
        "artifact_affinity_fields_absent": not forbidden_artifact,
        "report_zero_affinity_value_reads": report["affinity_values_read"] is False
        and report["affinity_value_fields_selected"] == 0,
        "panel_cells_match": panel_count_error == 0,
        "rectangle_counts_match": rectangle_error == 0,
        "dependency_assignments_match": dependency_match,
        "endpoint_summary_matches": endpoint_error <= 1e-12,
        "no_davis_or_recipient_reads": report["davis_label_reads"] == 0
        and report["recipient_label_reads"] == 0,
    }
    result = {
        "schema": "MetaSieve.EAffX0PostrunAudit.v1",
        "stage": report["stage"],
        "verdict": "POSTRUN_AUDIT_PASS" if all(checks.values()) else "POSTRUN_AUDIT_FAIL",
        "checks": checks,
        "panel_count_error": panel_count_error,
        "rectangle_count_error": rectangle_error,
        "endpoint_summary_max_error": endpoint_error,
        "recomputed_endpoints": recomputed_endpoint,
        "scientific_verdict": report["verdict"],
    }
    (root / "postrun_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _write = {
        "stage": report["stage"],
        "postrun_audit_sha256": sha256_file(root / "postrun_audit.json"),
        "auditor_sha256": sha256_file(Path(__file__)),
    }
    (root / "POSTRUN_AUDIT_MANIFEST.json").write_text(
        json.dumps(_write, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))
