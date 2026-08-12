"""V2-B2: label-free census of measured protein x ligand rectangles.

The crossed-contrast objective consumes 2x2 blocks ``{P1,P2} x {L1,L2}`` in
which **all four cells are measured** under the same document/panel.  This
script counts that supply without reading a single affinity value, so it can be
run before any preregistration is fixed.

Reported per source:

  * eligible rectangles;
  * how many cross a CD-HIT-40 family boundary (the only ones that can teach
    partner identity rather than within-family calibration);
  * how many have scaffold-distinct ligand pairs;
  * how many admit a legitimate k-shot support episode for **both** proteins
    that excludes the rectangle ligands and their scaffolds;
  * the dependency-component structure of the eligible set, since rectangles
    sharing a protein, ligand or document are not independent.

A large rectangle count with one giant dependency component is a development
supply, not a confirmation supply.  That distinction is enforced here.
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
REPORT = ROOT / "report/meta_fewshot/v2_rectangle_census.json"
SCHEMA = "MetaSieve.V2RectangleCensus.v1"
FORBIDDEN_SPLITS = ("meta_test",)
MAX_RECTANGLES_PER_BLOCK = 20000


class DisjointSet:
    def __init__(self) -> None:
        self.parent: dict = {}

    def find(self, value):
        self.parent.setdefault(value, value)
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left, right) -> None:
        left, right = self.find(left), self.find(right)
        if left != right:
            self.parent[right] = left


def load_design(corpus: Path, block_field: str) -> list[dict]:
    """Read identity columns only.  ``pK`` is never touched."""
    with gzip.open(corpus / "cells.jsonl.gz", "rt", encoding="utf-8") as handle:
        rows = [json.loads(line) for line in handle if line.strip()]
    ligands = {
        row["drug_key"]: row.get("scaffold", "")
        for row in map(json.loads,
                       (corpus / "ligands.jsonl").read_text(encoding="utf-8").splitlines())
    }
    design = []
    for row in rows:
        if row["split"] in FORBIDDEN_SPLITS:
            continue
        design.append({
            "target": row["target_id"],
            "ligand": row["ligand_id"],
            "cluster": row["protein_group_40"],
            "scaffold": ligands.get(row["ligand_id"], ""),
            "block": str(row.get(block_field, "ALL")),
        })
    if not design:
        raise RuntimeError("no source-side design rows survived the split filter")
    return design


def census(design: list[dict], k: int) -> dict:
    by_block = defaultdict(list)
    for row in design:
        by_block[row["block"]].append(row)

    ligands_by_target: dict[str, set] = defaultdict(set)
    scaffolds_by_target: dict[str, set] = defaultdict(set)
    for row in design:
        ligands_by_target[row["target"]].add(row["ligand"])
        scaffolds_by_target[row["target"]].add(row["scaffold"])

    cluster_of = {row["target"]: row["cluster"] for row in design}
    scaffold_of = {row["ligand"]: row["scaffold"] for row in design}

    total = cross_family = scaffold_distinct = episode_feasible = 0
    truncated_blocks = []
    dsu = DisjointSet()
    kept = 0

    for block, rows in sorted(by_block.items()):
        measured = defaultdict(set)
        for row in rows:
            measured[row["target"]].add(row["ligand"])
        targets = sorted(measured)
        block_count = 0
        stop = False
        for left, right in itertools.combinations(targets, 2):
            shared = sorted(measured[left] & measured[right])
            if len(shared) < 2:
                continue
            for first, second in itertools.combinations(shared, 2):
                total += 1
                block_count += 1
                if block_count > MAX_RECTANGLES_PER_BLOCK:
                    truncated_blocks.append(block)
                    stop = True
                    break
                if cluster_of[left] == cluster_of[right]:
                    continue
                cross_family += 1
                if scaffold_of[first] == scaffold_of[second]:
                    continue
                scaffold_distinct += 1
                blocked = {scaffold_of[first], scaffold_of[second]}
                feasible = all(
                    len({ligand for ligand in ligands_by_target[target]
                         if ligand not in {first, second}
                         and scaffold_of.get(ligand, "") not in blocked}) >= k
                    for target in (left, right)
                )
                if not feasible:
                    continue
                episode_feasible += 1
                kept += 1
                anchor = ("rect", kept)
                for key in (("protein", left), ("protein", right),
                            ("ligand", first), ("ligand", second),
                            ("block", block)):
                    dsu.union(anchor, key)
            if stop:
                break

    sizes: dict = defaultdict(int)
    for index in range(1, kept + 1):
        sizes[dsu.find(("rect", index))] += 1
    ordered = sorted(sizes.values(), reverse=True)
    largest_share = float(ordered[0] / kept) if kept else 0.0

    return {
        "design_rows": len(design),
        "targets": len(ligands_by_target),
        "blocks": len(by_block),
        "rectangles_total": total,
        "rectangles_cross_family": cross_family,
        "rectangles_cross_family_scaffold_distinct": scaffold_distinct,
        "rectangles_episode_feasible": episode_feasible,
        "dependency_components": len(ordered),
        "largest_component_share": largest_share,
        "truncated_blocks": sorted(set(truncated_blocks)),
        "support_budget_k": k,
    }


def run(corpus: Path, block_field: str, k: int) -> dict:
    result = {
        "schema": SCHEMA,
        "declared_role": "LABEL_FREE_DESIGN_CENSUS_NOT_CONFIRMATORY",
        "affinity_values_read": 0,
        "block_field": block_field,
        "forbidden_splits_dropped_before_any_read": list(FORBIDDEN_SPLITS),
        **census(load_design(corpus, block_field), k),
    }
    feasible = result["rectangles_episode_feasible"]
    components = result["dependency_components"]
    if feasible == 0:
        verdict = "NO_MEASURED_CROSSED_SUPERVISION_AVAILABLE"
    elif components >= 30 and result["largest_component_share"] <= 0.5:
        verdict = "CROSSED_SUPERVISION_SUPPLY_AND_INDEPENDENT_UNITS_AVAILABLE"
    else:
        verdict = "CROSSED_SUPERVISION_DEVELOPMENT_SUPPLY_ONLY"
    result["verdict"] = verdict
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corpus", type=Path, default=CORPUS)
    parser.add_argument("--block-field", default="document_id",
                        help="identity column defining an assay/document block")
    parser.add_argument("--k", type=int, default=5)
    parser.add_argument("--output", type=Path, default=REPORT)
    args = parser.parse_args()
    result = run(args.corpus, args.block_field, args.k)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
