"""Label-blind D0/D1 census for a release-pinned affinity corpus."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
from statistics import median
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.source_affinity.common import sha256_file, write_canonical_json


def read_jsonl(path: Path):
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                yield json.loads(line)


class Components:
    def __init__(self, nodes):
        self.parent = {node: node for node in nodes}

    def find(self, node):
        while self.parent[node] != node:
            self.parent[node] = self.parent[self.parent[node]]
            node = self.parent[node]
        return node

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[max(left, right)] = min(left, right)


def audit(tasks: list[dict], homology_assignments: list[dict] | None = None) -> dict:
    eligible = [task for task in tasks if task["eligible_e0_core"]]
    targets = sorted({task["protein_sequence_sha256"] for task in eligible})
    documents: dict[str, set[str]] = defaultdict(set)
    for task in eligible:
        for document in task["document_ids"]:
            documents[document].add(task["protein_sequence_sha256"])

    document_components = Components(targets)
    for linked_targets in documents.values():
        ordered = sorted(linked_targets)
        for target in ordered[1:]:
            document_components.union(ordered[0], target)

    homology_by_target = {}
    excluded_targets = set()
    if homology_assignments is not None:
        for row in homology_assignments:
            target = row["protein_sequence_sha256"]
            homology_by_target[target] = row["homology_component_id"]
            if row.get("excluded_by_davis_protected_homology"):
                excluded_targets.add(target)
        missing = set(targets) - set(homology_by_target)
        if missing:
            raise ValueError(f"homology assignments missing {len(missing)} eligible targets")

    union_components = Components(targets)
    for linked_targets in documents.values():
        ordered = sorted(linked_targets)
        for target in ordered[1:]:
            union_components.union(ordered[0], target)
    if homology_by_target:
        grouped: dict[str, list[str]] = defaultdict(list)
        for target in targets:
            grouped[homology_by_target[target]].append(target)
        for linked_targets in grouped.values():
            for target in linked_targets[1:]:
                union_components.union(linked_targets[0], target)

    governed = [task for task in eligible if task["protein_sequence_sha256"] not in excluded_targets]
    governed_targets = {
        task["protein_sequence_sha256"] for task in governed
    }
    component_task_depth = Counter(
        union_components.find(task["protein_sequence_sha256"]) for task in governed
    )
    compounds = [task["exact_compound_count"] for task in governed]
    comparisons = [task["non_tied_pair_comparisons"] for task in governed]
    return {
        "schema": "MetaSieve.AffinityCorpusAudit.v1",
        "stage": "P1R2B-D0/D1",
        "task_definition": "target x endpoint x assay x context",
        "all_tasks": len(tasks),
        "e0_core_tasks_before_governance": len(eligible),
        "e0_core_tasks_after_protected_exclusion": len(governed),
        "targets": len(governed_targets),
        "endpoint_distribution": dict(sorted(Counter(task["endpoint_family"] for task in governed).items())),
        "median_compounds_per_task": median(compounds) if compounds else 0,
        "median_non_tied_comparisons_per_task": median(comparisons) if comparisons else 0,
        "document_count": len(documents),
        "document_components": len({document_components.find(target) for target in governed_targets}),
        "homology_components": (
            len({homology_by_target[target] for target in governed_targets}) if homology_by_target else None
        ),
        "target_document_closure_components": (
            len(component_task_depth) if homology_by_target else None
        ),
        "independent_component_depth": sorted(component_task_depth.values(), reverse=True),
        "largest_component_task_fraction": (
            max(component_task_depth.values()) / len(governed) if component_task_depth and governed else None
        ),
        "davis_protected_excluded_targets_all_candidate_rows": len(excluded_targets),
        "davis_protected_excluded_targets_e0_core": len(set(targets) & excluded_targets),
        "homology_governance_complete": homology_assignments is not None,
        "minimum_task_contract": {
            "exact_compounds": 20,
            "non_tied_pair_comparisons": 12,
            "endpoint_pooling": False,
        },
        "recipient_labels_read": False,
        "training_authorized": False,
    }


def build_split(tasks: list[dict], homology_assignments: list[dict], folds: int = 5) -> list[dict]:
    if folds < 2:
        raise ValueError("at least two OOF folds are required")
    eligible = [task for task in tasks if task["eligible_e0_core"]]
    targets = sorted({task["protein_sequence_sha256"] for task in eligible})
    assignment_by_target = {row["protein_sequence_sha256"]: row for row in homology_assignments}
    if set(targets) - set(assignment_by_target):
        raise ValueError("homology assignments do not cover all eligible targets")
    retained = [
        task for task in eligible
        if not assignment_by_target[task["protein_sequence_sha256"]].get(
            "excluded_by_davis_protected_homology"
        )
    ]
    components = Components(targets)
    by_homology: dict[str, list[str]] = defaultdict(list)
    by_document: dict[str, set[str]] = defaultdict(set)
    for target in targets:
        by_homology[assignment_by_target[target]["homology_component_id"]].append(target)
    for task in eligible:
        for document in task["document_ids"]:
            by_document[document].add(task["protein_sequence_sha256"])
    for group in [*by_homology.values(), *by_document.values()]:
        ordered = sorted(group)
        for target in ordered[1:]:
            components.union(ordered[0], target)
    tasks_by_component: dict[str, list[dict]] = defaultdict(list)
    for task in retained:
        tasks_by_component[components.find(task["protein_sequence_sha256"])].append(task)
    ordered_components = sorted(
        tasks_by_component,
        key=lambda component: (
            -len(tasks_by_component[component]),
            hashlib.sha256(component.encode()).hexdigest(),
        ),
    )
    loads = [0] * folds
    component_fold = {}
    for component in ordered_components:
        fold = min(range(folds), key=lambda index: (loads[index], index))
        component_fold[component] = fold
        loads[fold] += len(tasks_by_component[component])
    return [
        {
            "task_id": task["task_id"],
            "protein_sequence_sha256": task["protein_sequence_sha256"],
            "closure_component_id": components.find(task["protein_sequence_sha256"]),
            "outer_oof_fold": component_fold[components.find(task["protein_sequence_sha256"])],
        }
        for task in sorted(retained, key=lambda row: row["task_id"])
    ]


def verify_split(tasks: list[dict], homology_assignments: list[dict], split: list[dict]) -> dict:
    task_by_id = {task["task_id"]: task for task in tasks if task["eligible_e0_core"]}
    assignment = {row["protein_sequence_sha256"]: row for row in homology_assignments}
    fold_by_task = {row["task_id"]: row["outer_oof_fold"] for row in split}
    retained_ids = {
        task_id for task_id, task in task_by_id.items()
        if not assignment[task["protein_sequence_sha256"]].get(
            "excluded_by_davis_protected_homology"
        )
    }
    if set(fold_by_task) != retained_ids:
        raise ValueError("split does not cover exactly the retained E0-Core tasks")
    homology_folds: dict[str, set[int]] = defaultdict(set)
    document_folds: dict[str, set[int]] = defaultdict(set)
    closure_folds: dict[str, set[int]] = defaultdict(set)
    closure_by_task = {row["task_id"]: row["closure_component_id"] for row in split}
    for task_id in retained_ids:
        task = task_by_id[task_id]
        fold = fold_by_task[task_id]
        homology_folds[assignment[task["protein_sequence_sha256"]]["homology_component_id"]].add(fold)
        closure_folds[closure_by_task[task_id]].add(fold)
        for document in task["document_ids"]:
            document_folds[document].add(fold)
    fold_counts = Counter(fold_by_task.values())
    checks = {
        "retained_tasks_exactly_once": len(fold_by_task) == len(retained_ids),
        "homology_components_straddling": sum(len(folds) > 1 for folds in homology_folds.values()),
        "documents_straddling": sum(len(folds) > 1 for folds in document_folds.values()),
        "closure_components_straddling": sum(len(folds) > 1 for folds in closure_folds.values()),
    }
    if not checks["retained_tasks_exactly_once"] or any(
        checks[key] for key in (
            "homology_components_straddling", "documents_straddling",
            "closure_components_straddling",
        )
    ):
        raise RuntimeError(f"closure split verification failed: {checks}")
    return {**checks, "fold_task_counts": dict(sorted(fold_counts.items()))}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-manifest", type=Path, required=True)
    parser.add_argument("--homology-assignments", type=Path)
    parser.add_argument("--release-manifest", type=Path)
    parser.add_argument("--corpus-manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--split-output", type=Path)
    args = parser.parse_args()
    tasks = json.loads(args.task_manifest.read_text(encoding="utf-8"))
    homology = list(read_jsonl(args.homology_assignments)) if args.homology_assignments else None
    result = audit(tasks, homology)
    if args.split_output:
        if homology is None:
            parser.error("--split-output requires --homology-assignments")
        split = build_split(tasks, homology)
        args.split_output.parent.mkdir(parents=True, exist_ok=True)
        with args.split_output.open("w", encoding="utf-8", newline="\n") as handle:
            for row in split:
                handle.write(json.dumps(row, sort_keys=True, separators=(",", ":")) + "\n")
        result["split"] = {
            "schema": "MetaSieve.AffinityClosureOOFSplit.v1",
            "folds": 5,
            "tasks": len(split),
            "sha256": sha256_file(args.split_output),
            "verification": verify_split(tasks, homology, split),
        }
    result["inputs"] = {
        "task_manifest_sha256": sha256_file(args.task_manifest),
        "homology_assignments_sha256": (
            sha256_file(args.homology_assignments) if args.homology_assignments else None
        ),
        "release_manifest_sha256": (
            sha256_file(args.release_manifest) if args.release_manifest else None
        ),
        "corpus_manifest_sha256": (
            sha256_file(args.corpus_manifest) if args.corpus_manifest else None
        ),
    }
    write_canonical_json(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
