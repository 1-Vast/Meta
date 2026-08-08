"""Independently audit the E-AFF-X0-FEAS unit-feasibility artifacts.

Independence: the closure is rebuilt with a separate union-find, the ceilings
are recounted from the governed inputs, the structural claim is recounted from
the X0 cell record, and the frozen `245` is re-derived by bisection on
`chi2.cdf` rather than by the forward scan on `chi2.sf` used in the runner.
"""
from __future__ import annotations

from collections import defaultdict
import json
from pathlib import Path

from scipy.stats import chi2

from scripts.source_affinity.common import sha256_file


ROOT = Path("research/e0_identifiability/artifacts/eaff_x0_feas_v1")
X0 = Path("research/e0_identifiability/artifacts/eaff_x0_v1")
INPUT = Path("dataset/processed/source_affinity/e0_input_v1")
CORPUS = Path("dataset/processed/source_affinity/energy_pilot_v1")
GOVERNANCE = Path("dataset/processed/source_affinity/energy_pilot_v1_governance")
FORBIDDEN_KEYS = {
    "p_affinity", "standard_value", "published_value", "pchembl_value_reported",
    "activity_value", "label", "y",
}
REQUIRED = 245


def _read_jsonl(path: Path) -> list[dict]:
    with path.open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]


def _keys(value) -> set[str]:
    if isinstance(value, dict):
        return set(value) | {key for item in value.values() for key in _keys(item)}
    if isinstance(value, list):
        return {key for item in value for key in _keys(item)}
    return set()


def _roots(edges: list[tuple[str, str]], nodes: set[str]) -> dict[str, str]:
    """Independent connected-components pass: repeated relabelling to a fixed point."""
    label = {node: node for node in nodes}
    changed = True
    while changed:
        changed = False
        for left, right in edges:
            low = min(label[left], label[right])
            for node in (left, right):
                if label[node] != low:
                    label[node] = low
                    changed = True
        for node in nodes:
            while label[node] != label[label[node]]:
                label[node] = label[label[node]]
                changed = True
    return label


def required_n_by_bisection(ratio: float, alpha: float, power: float) -> int:
    """Smallest n with power >= target, found by bisection on chi2.cdf."""
    def attained(n: int) -> float:
        df = n - 1
        return 1.0 - chi2.cdf(chi2.isf(alpha, df) / ratio, df)
    low, high = 2, 4096
    while low < high:
        mid = (low + high) // 2
        if attained(mid) >= power:
            high = mid
        else:
            low = mid + 1
    return low


def run(root: Path = ROOT) -> dict:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    report = json.loads((root / "report.json").read_text(encoding="utf-8"))
    hashes_match = all(sha256_file(root / name) == digest
                       for name, digest in manifest["outputs"].items())
    inputs_match = all(
        sha256_file(path) == manifest["inputs"][key]
        for key, path in (
            ("label_blind_rows", INPUT / "rows.label_blind.jsonl"),
            ("input_manifest", INPUT / "manifest.json"),
            ("x0_cells", X0 / "cells.jsonl"),
            ("x0_report", X0 / "report.json"),
            ("d0_task_manifest", CORPUS / "task_manifest.json"),
            ("d1_homology_assignments", GOVERNANCE / "homology_assignments.jsonl"),
        ))

    rows = _read_jsonl(INPUT / "rows.label_blind.jsonl")
    forbidden_rows = bool({key for row in rows for key in row} & FORBIDDEN_KEYS)
    forbidden_artifact = bool(_keys(report) & FORBIDDEN_KEYS)

    # --- structural claim, recounted from the X0 cell record -----------------
    closures_by_panel: dict[str, set[str]] = defaultdict(set)
    panel_endpoint: dict[str, str] = {}
    for cell in _read_jsonl(X0 / "cells.jsonl"):
        closures_by_panel[cell["panel_id"]].add(cell["closure_component_id"])
        panel_endpoint[cell["panel_id"]] = cell["endpoint_family"]
    multi_closure = sum(1 for value in closures_by_panel.values() if len(value) > 1)

    # --- E0-input ceilings ---------------------------------------------------
    universe = {row["closure_component_id"] for row in rows}
    with_data: dict[str, set[str]] = defaultdict(set)
    ligands: dict[tuple[str, str], dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for row in rows:
        endpoint, component = row["endpoint_family"], row["closure_component_id"]
        with_data[endpoint].add(component)
        ligands[(endpoint, component)][row["protein_sequence_sha256"]].add(
            row["ligand_connectivity_key"])

    def capable(pool: dict[str, set[str]]) -> bool:
        proteins = sorted(pool)
        return any(len(pool[a] & pool[b]) >= 2
                   for i, a in enumerate(proteins) for b in proteins[i + 1:])

    ceilings = {}
    for endpoint in ("Ki", "Kd"):
        capable_components = sum(
            1 for component in with_data[endpoint] if capable(ligands[(endpoint, component)]))
        ceilings[endpoint] = {
            "with_endpoint_data": len(with_data[endpoint]),
            "panel_free_rectangle_capable": capable_components,
        }

    # --- governed-corpus sweep, rebuilt independently ------------------------
    tasks = json.loads((CORPUS / "task_manifest.json").read_text(encoding="utf-8"))
    homology = _read_jsonl(GOVERNANCE / "homology_assignments.jsonl")
    component_of = {row["protein_sequence_sha256"]: row["homology_component_id"]
                    for row in homology}
    excluded = {row["protein_sequence_sha256"] for row in homology
                if row.get("excluded_by_davis_protected_homology")}
    filters = {
        "e0_core_eligible": lambda task: bool(task["eligible_e0_core"]),
        "full_d0_governed": lambda task: True,
        "min_5_compounds_per_task": lambda task: task["exact_compound_count"] >= 5,
        "min_2_compounds_per_task": lambda task: task["exact_compound_count"] >= 2,
    }
    sweep_error = 0
    recomputed_sweep = {}
    published_sweep = {row["population"]: row
                       for row in report["governed_corpus_population_sweep"]["populations"]}
    for name, keep in filters.items():
        subset = [task for task in tasks
                  if task["protein_sequence_sha256"] not in excluded and keep(task)]
        nodes = {task["protein_sequence_sha256"] for task in subset}
        nodes |= {component_of[node] for node in nodes}
        edges = [(task["protein_sequence_sha256"],
                  component_of[task["protein_sequence_sha256"]]) for task in subset]
        by_document: dict[str, list[str]] = defaultdict(list)
        for task in subset:
            for document in task["document_ids"]:
                by_document[document].append(task["protein_sequence_sha256"])
        for group in by_document.values():
            ordered = sorted(set(group))
            edges.extend((ordered[0], node) for node in ordered[1:])
        label = _roots(edges, nodes)
        entry = {"tasks": len(subset),
                 "proteins": len({task["protein_sequence_sha256"] for task in subset}),
                 "closure_components": len({label[task["protein_sequence_sha256"]]
                                            for task in subset}),
                 "endpoints": {}}
        for endpoint in ("Ki", "Kd"):
            grouped: dict[str, set[str]] = defaultdict(set)
            for task in subset:
                if task["endpoint_family"] == endpoint:
                    protein = task["protein_sequence_sha256"]
                    grouped[label[protein]].add(protein)
            entry["endpoints"][endpoint] = {
                "components_with_endpoint_data": len(grouped),
                "components_with_two_or_more_proteins": sum(
                    1 for value in grouped.values() if len(value) >= 2)}
        recomputed_sweep[name] = entry
        reference = published_sweep[name]
        sweep_error = max(
            sweep_error,
            abs(entry["tasks"] - reference["tasks"]),
            abs(entry["proteins"] - reference["proteins"]),
            abs(entry["closure_components"] - reference["closure_components"]),
            max(abs(entry["endpoints"][endpoint][key]
                    - reference["endpoints"][endpoint][key])
                for endpoint in ("Ki", "Kd")
                for key in ("components_with_endpoint_data",
                            "components_with_two_or_more_proteins")))

    best_ceiling = {
        endpoint: max(entry["endpoints"][endpoint]["components_with_two_or_more_proteins"]
                      for entry in recomputed_sweep.values())
        for endpoint in ("Ki", "Kd")}

    derived = required_n_by_bisection(1.25, 0.05, 0.80)

    checks = {
        "preexisting_hashes_match": hashes_match,
        "declared_inputs_match": inputs_match,
        "input_rows_carry_no_value_field": not forbidden_rows,
        "artifact_affinity_fields_absent": not forbidden_artifact,
        "report_zero_affinity_value_fields": report["affinity_value_fields_selected"] == 0
        and report["affinity_values_read"] is False,
        "no_davis_or_recipient_reads": report["davis_label_reads"] == 0
        and report["recipient_label_reads"] == 0,
        "no_training_performed": report["training_performed"] is False,
        "structural_claim_reproduces": multi_closure == 0
        and report["structural_claim"]["panels_touching_more_than_one_closure_component"] == 0,
        "closure_universe_reproduces": len(universe) == report["closure_component_universe"],
        "e0_input_ceilings_reproduce": all(
            ceilings[endpoint]["with_endpoint_data"]
            == report["endpoints"][endpoint]["closure_components_with_endpoint_data"]
            and ceilings[endpoint]["panel_free_rectangle_capable"]
            == report["endpoints"][endpoint]["closure_components_rectangle_capable_panel_free"]
            for endpoint in ("Ki", "Kd")),
        "population_sweep_reproduces": sweep_error == 0,
        "best_ceiling_reproduces": best_ceiling
        == report["governed_corpus_population_sweep"]["best_ceiling_over_all_populations"],
        "frozen_245_reproduces_by_bisection": derived == REQUIRED
        and report["power_design"]["frozen_required_n"] == REQUIRED,
        "requirement_unattainable_at_every_population": all(
            best_ceiling[endpoint] < REQUIRED for endpoint in ("Ki", "Kd")),
    }
    result = {
        "schema": "MetaSieve.EAffX0FeasPostrunAudit.v1",
        "stage": report["stage"],
        "verdict": "POSTRUN_AUDIT_PASS" if all(checks.values()) else "POSTRUN_AUDIT_FAIL",
        "checks": checks,
        "recomputed_closure_universe": len(universe),
        "recomputed_e0_input_ceilings": ceilings,
        "recomputed_best_ceiling": best_ceiling,
        "recomputed_required_n": derived,
        "population_sweep_max_error": sweep_error,
        "panels_touching_more_than_one_closure_component": multi_closure,
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
