"""Label-blind dependency closure for measured cross-protein selectivity groups."""
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from research.meta_fewshot.fs_corpus_rebuild import homology_candidates, local_identity
from research.meta_fewshot.v1_source_supervision_audit import read_gzip_jsonl, sha256

ROOT = Path(__file__).resolve().parents[2]
MAIN = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
DEV = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_v1_development"
OUT = ROOT / "report/meta_fewshot/a1_selectivity_dependency_audit.json"
Z_SUM = 1.645 + 0.842


class UnionFind:
    def __init__(self, keys):
        self.parent = {key: key for key in keys}

    def find(self, key):
        root = key
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[key] != key:
            key, self.parent[key] = self.parent[key], root
        return root

    def union(self, left, right):
        left, right = self.find(left), self.find(right)
        if left != right:
            keep, drop = sorted((left, right))
            self.parent[drop] = keep


def union_shared(components, groups):
    for values in groups.values():
        ordered = sorted(values)
        for value in ordered[1:]:
            components.union(ordered[0], value)


def audit(main=MAIN, dev=DEV):
    groups = [row for row in read_gzip_jsonl(dev / "source_contrast_groups.jsonl.gz")
              if row["kind"] == "measured_partner"]
    ligand_rows = [json.loads(line) for line in (main / "ligands.jsonl").read_text().splitlines()]
    scaffold = {row["drug_key"]: row["scaffold"] for row in ligand_rows}
    protein_rows = [json.loads(line) for line in (main / "proteins.jsonl").read_text().splitlines()]
    sequence = {row["sequence_sha256"]: row["sequence"] for row in protein_rows}
    group_ids = [f"{index:06d}" for index in range(len(groups))]
    components = UnionFind(group_ids)
    by_document, by_ligand, by_scaffold, by_family = (defaultdict(set) for _ in range(4))
    families_in_groups = set()
    targets_in_groups = set()
    for group_id, group in zip(group_ids, groups):
        document = group["panel_id"].split("|Ki|", 1)[0]
        by_document[document].add(group_id)
        by_ligand[group["ligand_id"]].add(group_id)
        by_scaffold[scaffold[group["ligand_id"]]].add(group_id)
        for member in group["members"]:
            by_family[member["protein_group_40"]].add(group_id)
            families_in_groups.add(member["protein_group_40"])
            targets_in_groups.add(member["target_id"])
    for mapping in (by_document, by_ligand, by_scaffold, by_family):
        union_shared(components, mapping)

    relevant_sequences = {target: sequence[target] for target in targets_in_groups}
    candidates = homology_candidates(relevant_sequences)
    homology_edges = []
    target_groups = {member["target_id"]: member["protein_group_40"]
                     for group in groups for member in group["members"]}
    family_edges = set()
    for left, right in sorted(candidates):
        identity = local_identity(relevant_sequences[left], relevant_sequences[right])
        if identity >= 0.40 and target_groups[left] != target_groups[right]:
            family_edges.add(tuple(sorted((target_groups[left], target_groups[right]))))
            homology_edges.append({"left": left, "right": right, "identity": identity})
    family_components = UnionFind(families_in_groups)
    for left, right in family_edges:
        family_components.union(left, right)
    by_verified_family = defaultdict(set)
    for family, group_set in by_family.items():
        by_verified_family[family_components.find(family)].update(group_set)
    union_shared(components, by_verified_family)

    summary = defaultdict(lambda: {"groups": 0, "families": set(), "ligands": set(),
                                   "scaffolds": set(), "documents": set(), "targets": set()})
    for group_id, group in zip(group_ids, groups):
        value = summary[components.find(group_id)]
        value["groups"] += 1
        value["ligands"].add(group["ligand_id"])
        value["scaffolds"].add(scaffold[group["ligand_id"]])
        value["documents"].add(group["panel_id"].split("|Ki|", 1)[0])
        value["families"].update(member["protein_group_40"] for member in group["members"])
        value["targets"].update(member["target_id"] for member in group["members"])
    component_rows = sorted(({
        "component": root, "groups": value["groups"],
        **{key: len(value[key]) for key in ("families", "ligands", "scaffolds", "documents", "targets")},
    } for root, value in summary.items()), key=lambda row: (-row["groups"], row["component"]))
    count = len(component_rows)
    mde = float(Z_SUM / np.sqrt(count)) if count else None
    identifiable = count >= 18 and mde <= 0.600 and all(row["families"] >= 2 for row in component_rows)
    return {
        "schema": "MetaSieve.A1SelectivityDependencyAudit.v1",
        "label_values_used": 0,
        "eligible_panel_ligand_groups": len(groups),
        "targets": len(targets_in_groups), "frozen_cdhit40_families": len(families_in_groups),
        "homology_candidates": len(candidates),
        "cross_cdhit40_local_identity_edges_ge_0_40": len(homology_edges),
        "components": count, "component_unit_mde": mde,
        "largest_component_group_share": component_rows[0]["groups"] / len(groups) if groups else None,
        "probe_identifiable": identifiable,
        "TERMINAL_VERDICT": ("SELECTIVITY_PROBE_AUTHORIZED" if identifiable
                             else "SELECTIVITY_PROBE_NOT_IDENTIFIABLE"),
        "components_detail": component_rows,
        "group_component_assignment": {
            group_id: components.find(group_id) for group_id in group_ids
        },
        "inputs": {"development_manifest_sha256": sha256(dev / "manifest.json"),
                   "contrast_groups_sha256": sha256(dev / "source_contrast_groups.jsonl.gz")},
    }


if __name__ == "__main__":
    result = audit()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
