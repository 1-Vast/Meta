"""Label-blind feasibility audit of the frozen E-AFF-X0 independence unit.

X0-FEAS reads only label-blind design metadata that X0 already produced and the
governed label-blind row index. It selects no affinity field, trains nothing and
issues no affinity verdict. Its single question is whether the frozen
`245 effective components per endpoint` requirement is attainable at all under
the unit definition X0 used, on this corpus or on any corpus governed the same
way.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from datetime import datetime, timezone
import json
from pathlib import Path

from scipy.stats import chi2

from scripts.source_affinity.common import sha256_file


STAGE = "P1R2B-E-AFF-X0-FEAS_UNIT_FEASIBILITY_AUDIT"
REQUIRED_COMPONENTS = 245
INTERACTION_RMS_OVER_NOISE = 0.5
ALPHA_ONE_SIDED = 0.05
TARGET_POWER = 0.80
PROHIBITED_ROW_FIELDS = (
    "standard_value", "value", "pchembl_value", "pchembl", "affinity",
    "pki", "pkd", "y", "label",
)
# The D0 task manifest is admitted field by field. Anything outside this set is
# treated as a possible label channel and fails the audit closed.
ALLOWED_TASK_FIELDS = frozenset({
    "task_id", "task_keys", "protein_sequence_sha256", "endpoint_family",
    "assay_chembl_id", "assay_context_sha256", "document_ids",
    "eligible_e0_core", "exact_compound_count", "measurement_count",
    "pair_comparisons", "non_tied_pair_comparisons",
})
# Population filters swept over the full governed D0 corpus. Each is a rule on
# task depth only; none reads a measured value.
POPULATION_FILTERS = (
    ("e0_core_eligible", lambda task: bool(task["eligible_e0_core"])),
    ("full_d0_governed", lambda task: True),
    ("min_5_compounds_per_task", lambda task: task["exact_compound_count"] >= 5),
    ("min_2_compounds_per_task", lambda task: task["exact_compound_count"] >= 2),
)


def _read_jsonl(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


def _write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def chi_square_required_n(ratio: float, alpha: float, power: float) -> int:
    """Smallest n whose one-sided chi-square variance test reaches `power`."""
    n = 2
    while n <= 100000:
        df = n - 1
        if chi2.sf(chi2.ppf(1.0 - alpha, df) / ratio, df) >= power:
            return n
        n += 1
    raise RuntimeError("required n did not converge")


class UnionFind:
    def __init__(self):
        self.parent: dict[str, str] = {}

    def add(self, value: str) -> None:
        self.parent.setdefault(value, value)

    def find(self, value: str) -> str:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: str, right: str) -> None:
        self.add(left)
        self.add(right)
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def closure_over(tasks: list[dict], homology_of: dict[str, str]) -> UnionFind:
    """Rebuild the D1 closure: 40pct homology groups unioned with shared documents."""
    union = UnionFind()
    by_document: dict[str, set[str]] = defaultdict(set)
    for task in tasks:
        protein = task["protein_sequence_sha256"]
        union.add(protein)
        union.union(protein, homology_of[protein])
        for document in task["document_ids"]:
            by_document[document].add(protein)
    for group in by_document.values():
        ordered = sorted(group)
        for protein in ordered[1:]:
            union.union(ordered[0], protein)
    return union


def population_sweep(tasks: list[dict], homology_of: dict[str, str],
                     excluded: set[str]) -> list[dict]:
    """Closure-component ceilings under each task-depth population filter.

    A closure component holding fewer than two distinct proteins for an endpoint
    can never host a rectangle for that endpoint, so the two-protein count is a
    hard upper bound on effective components at that population.
    """
    sweep = []
    for name, keep in POPULATION_FILTERS:
        rows = [task for task in tasks
                if task["protein_sequence_sha256"] not in excluded and keep(task)]
        union = closure_over(rows, homology_of)
        entry = {
            "population": name,
            "tasks": len(rows),
            "proteins": len({task["protein_sequence_sha256"] for task in rows}),
            "closure_components": len({union.find(task["protein_sequence_sha256"])
                                       for task in rows}),
            "endpoints": {},
        }
        for endpoint in ("Ki", "Kd"):
            proteins_by_component: dict[str, set[str]] = defaultdict(set)
            for task in rows:
                if task["endpoint_family"] == endpoint:
                    protein = task["protein_sequence_sha256"]
                    proteins_by_component[union.find(protein)].add(protein)
            multi = sum(1 for value in proteins_by_component.values() if len(value) >= 2)
            entry["endpoints"][endpoint] = {
                "components_with_endpoint_data": len(proteins_by_component),
                "components_with_two_or_more_proteins": multi,
                "requirement": REQUIRED_COMPONENTS,
                "requirement_exceeds_ceiling": REQUIRED_COMPONENTS > multi,
            }
        sweep.append(entry)
    return sweep


def greedy_disjoint_rectangles(target_ligands: dict[str, set[str]]) -> int:
    """Greedy lower bound on target- and ligand-disjoint rectangles in one panel.

    A packed rectangle consumes both of its targets and both of its ligands, so
    no packed rectangle shares a protein or a ligand with any other. This removes
    the target-repeat and ligand-repeat pseudoreplication that inflated the X0
    nominal rectangle counts. It does not remove shared-assay-run dependence.
    """
    remaining = {target: set(ligands) for target, ligands in target_ligands.items()}
    packed = 0
    while True:
        targets = sorted(target for target, ligands in remaining.items() if len(ligands) >= 2)
        chosen = None
        for index, left in enumerate(targets):
            for right in targets[index + 1:]:
                common = remaining[left] & remaining[right]
                if len(common) >= 2:
                    chosen = (left, right, sorted(common)[:2])
                    break
            if chosen:
                break
        if not chosen:
            return packed
        left, right, ligands = chosen
        del remaining[left]
        del remaining[right]
        for ligand in ligands:
            for available in remaining.values():
                available.discard(ligand)
        packed += 1


def audit(args) -> dict:
    output = Path(args.output)
    if output.exists():
        raise FileExistsError(f"output exists: {output}")

    input_root = Path(args.input)
    x0_root = Path(args.x0)
    corpus_root = Path(args.corpus)
    governance_root = Path(args.governance)

    rows = list(_read_jsonl(input_root / "rows.label_blind.jsonl"))
    leaked = sorted({field for row in rows for field in row if field.lower() in PROHIBITED_ROW_FIELDS})
    if leaked:
        raise RuntimeError(f"label firewall violation in X0-FEAS inputs: {leaked}")

    # ---- corpus universe -------------------------------------------------
    universe: set[str] = set()
    by_component: dict[str, dict[str, dict[str, set[str]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(set)))
    for row in rows:
        component = row["closure_component_id"]
        universe.add(component)
        by_component[component][row["endpoint_family"]][
            row["protein_sequence_sha256"]].add(row["ligand_connectivity_key"])

    def rectangle_capable(target_ligands: dict[str, set[str]]) -> bool:
        targets = sorted(target_ligands)
        for index, left in enumerate(targets):
            for right in targets[index + 1:]:
                if len(target_ligands[left] & target_ligands[right]) >= 2:
                    return True
        return False

    # ---- X0 panel/cell replay -------------------------------------------
    cells = list(_read_jsonl(x0_root / "cells.jsonl"))
    panel_targets: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    panel_endpoint: dict[str, str] = {}
    panel_closures: dict[str, set[str]] = defaultdict(set)
    for cell in cells:
        panel = cell["panel_id"]
        panel_endpoint[panel] = cell["endpoint_family"]
        panel_targets[panel][cell["protein_sequence_sha256"]].add(cell["ligand_connectivity_key"])
        panel_closures[panel].add(cell["closure_component_id"])

    # Structural claim: a rectangle needs two proteins inside one document-keyed
    # panel, and the D1 closure unions every pair of proteins that share a
    # document. Both proteins of every rectangle must therefore already sit in
    # the same closure component.
    multi_closure_panels = sorted(
        panel for panel, closures in panel_closures.items() if len(closures) > 1)

    endpoints = {}
    for endpoint in ("Ki", "Kd"):
        components_with_data = {
            component for component, families in by_component.items() if endpoint in families}
        components_capable = {
            component for component in components_with_data
            if rectangle_capable(by_component[component][endpoint])}

        endpoint_panels = [panel for panel, family in panel_endpoint.items() if family == endpoint]
        disjoint_by_component: dict[str, int] = defaultdict(int)
        disjoint_total = 0
        for panel in endpoint_panels:
            packed = greedy_disjoint_rectangles(panel_targets[panel])
            if not packed:
                continue
            disjoint_total += packed
            for closure in panel_closures[panel]:
                disjoint_by_component[closure] += packed

        x0_components = len({
            row["dependency_component_id"]
            for row in _read_jsonl(x0_root / "dependency_components.jsonl")
            if row["endpoint_family"] == endpoint and row["stratum"] == "all_rectangles"})

        endpoints[endpoint] = {
            "x0_effective_components": x0_components,
            "closure_components_with_endpoint_data": len(components_with_data),
            "closure_components_rectangle_capable_panel_free": len(components_capable),
            "attainable_ceiling_under_frozen_unit": len(components_with_data),
            "requirement": REQUIRED_COMPONENTS,
            "requirement_exceeds_ceiling": REQUIRED_COMPONENTS > len(components_with_data),
            "shortfall_against_ceiling": REQUIRED_COMPONENTS - len(components_with_data),
            "design_disjoint_rectangles_greedy": disjoint_total,
            "design_disjoint_components": len(disjoint_by_component),
            "design_disjoint_largest_component_fraction": (
                max(disjoint_by_component.values()) / disjoint_total if disjoint_total else 0.0),
        }

    # ---- governed-corpus population sweep --------------------------------
    tasks = json.loads((corpus_root / "task_manifest.json").read_text(encoding="utf-8"))
    unexpected = sorted({field for task in tasks for field in task} - ALLOWED_TASK_FIELDS)
    if unexpected:
        raise RuntimeError(f"label firewall violation in D0 task manifest: {unexpected}")
    homology = list(_read_jsonl(governance_root / "homology_assignments.jsonl"))
    homology_of = {row["protein_sequence_sha256"]: row["homology_component_id"]
                   for row in homology}
    davis_excluded = {row["protein_sequence_sha256"] for row in homology
                      if row.get("excluded_by_davis_protected_homology")}
    sweep = population_sweep(tasks, homology_of, davis_excluded)

    best_ceiling = {
        endpoint: max(entry["endpoints"][endpoint]["components_with_two_or_more_proteins"]
                      for entry in sweep)
        for endpoint in ("Ki", "Kd")
    }
    by_population = {entry["population"]: entry for entry in sweep}
    anti_monotone = {
        endpoint: (by_population["min_2_compounds_per_task"]["endpoints"][endpoint][
                       "components_with_two_or_more_proteins"]
                   > by_population["full_d0_governed"]["endpoints"][endpoint][
                       "components_with_two_or_more_proteins"])
        for endpoint in ("Ki", "Kd")
    }

    derived_required = chi_square_required_n(
        1.0 + INTERACTION_RMS_OVER_NOISE ** 2, ALPHA_ONE_SIDED, TARGET_POWER)

    attainable = [name for name, ceiling in best_ceiling.items()
                  if ceiling >= REQUIRED_COMPONENTS]
    if attainable:
        verdict = "X0_UNIT_REQUIREMENT_ATTAINABLE_" + "_".join(n.upper() for n in attainable)
    else:
        verdict = "X0_UNIT_REQUIREMENT_UNATTAINABLE_BY_CONSTRUCTION"

    output.mkdir(parents=True)
    report = {
        "schema": "MetaSieve.EAffX0Feas.v1",
        "stage": STAGE,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "verdict": verdict,
        "label_blind": True,
        "affinity_value_fields_selected": 0,
        "affinity_values_read": False,
        "davis_label_reads": 0,
        "recipient_label_reads": 0,
        "training_performed": False,
        "input_rows": len(rows),
        "closure_component_universe": len(universe),
        "closure_component_definition": (
            "union of 40pct-identity homology groups and shared-document target groups"),
        "power_design": {
            "interaction_rms_over_assay_noise": INTERACTION_RMS_OVER_NOISE,
            "alpha_one_sided": ALPHA_ONE_SIDED,
            "power": TARGET_POWER,
            "variance_ratio": 1.0 + INTERACTION_RMS_OVER_NOISE ** 2,
            "frozen_required_n": REQUIRED_COMPONENTS,
            "independently_derived_required_n": derived_required,
            "derivation_reproduced": derived_required == REQUIRED_COMPONENTS,
            "assumed_sampling_model": "independent_identically_distributed_units",
        },
        "structural_claim": {
            "statement": (
                "every rectangle joins two proteins inside one document-keyed panel, "
                "and D1 closure unions all proteins sharing a document, so both "
                "proteins of every rectangle always lie in one closure component"),
            "panels_touching_more_than_one_closure_component": len(multi_closure_panels),
            "claim_holds_on_x0_artifacts": not multi_closure_panels,
            "consequence": (
                "crossing can never create an X0 unit; effective components are "
                "bounded above by the closure-component universe"),
        },
        "endpoints": endpoints,
        "governed_corpus_population_sweep": {
            "note": (
                "components_with_two_or_more_proteins is a hard upper bound on "
                "effective components: a component holding one protein for an "
                "endpoint cannot host a rectangle for it"),
            "best_ceiling_over_all_populations": best_ceiling,
            "requirement_exceeds_best_ceiling": {
                endpoint: REQUIRED_COMPONENTS > ceiling
                for endpoint, ceiling in best_ceiling.items()},
            "ceiling_falls_when_corpus_grows": anti_monotone,
            "populations": sweep,
        },
        "interpretation_limits": [
            "X0-FEAS audits the estimand and its unit, not the biology",
            "it does not weaken, replace or re-register the frozen X0 requirement",
            "design-disjoint rectangle counts remove target/ligand repetition only; "
            "they still share an assay run and are not certified IID units",
            "an unattainable requirement is not evidence that interaction exists",
        ],
    }
    _write_json(output / "report.json", report)
    manifest = {
        "stage": STAGE,
        "label_blind": True,
        "inputs": {
            "label_blind_rows": sha256_file(input_root / "rows.label_blind.jsonl"),
            "input_manifest": sha256_file(input_root / "manifest.json"),
            "x0_cells": sha256_file(x0_root / "cells.jsonl"),
            "x0_dependency_components": sha256_file(x0_root / "dependency_components.jsonl"),
            "x0_report": sha256_file(x0_root / "report.json"),
            "x0_manifest": sha256_file(x0_root / "manifest.json"),
            "d0_task_manifest": sha256_file(corpus_root / "task_manifest.json"),
            "d1_homology_assignments": sha256_file(
                governance_root / "homology_assignments.jsonl"),
            "registration": sha256_file(Path(__file__).with_name("EAFF_X0_FEAS_REGISTRATION.md")),
        },
        "outputs": {"report.json": sha256_file(output / "report.json")},
        "label_reads": {"affinity_values": 0, "davis": 0, "recipient": 0},
    }
    _write_json(output / "manifest.json", manifest)
    return report


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", default="dataset/processed/source_affinity/e0_input_v1")
    parser.add_argument("--x0", default="research/e0_identifiability/artifacts/eaff_x0_v1")
    parser.add_argument("--corpus", default="dataset/processed/source_affinity/energy_pilot_v1")
    parser.add_argument("--governance",
                        default="dataset/processed/source_affinity/energy_pilot_v1_governance")
    parser.add_argument("--output",
                        default="research/e0_identifiability/artifacts/eaff_x0_feas_v1")
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(audit(parse_args()), indent=2, sort_keys=True))
