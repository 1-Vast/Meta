"""Pure contracts for the label-blind crossed target-ligand census."""
from __future__ import annotations

import hashlib
import json


TARGET_INDEPENDENT_CONTEXT_KEYS = (
    "assay_organism",
    "bao_format",
    "cell_id",
    "tissue_id",
    "subcellular_fraction",
    "relationship_type",
    "parameters",
)


def stable_hash(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)
    return hashlib.sha256(payload.encode()).hexdigest()


def crossed_panel_payload(row: dict, parameters: list[dict]) -> dict:
    if row.get("variant_id") is not None:
        raise ValueError("primary crossed panel excludes variant assays")
    context = {
        "assay_organism": row.get("assay_organism"),
        "bao_format": row.get("bao_format"),
        "cell_id": row.get("cell_id"),
        "tissue_id": row.get("tissue_id"),
        "subcellular_fraction": row.get("subcellular_fraction"),
        "relationship_type": row.get("relationship_type"),
        "parameters": parameters,
    }
    return {
        "document_chembl_id": row["document_chembl_id"],
        "endpoint_family": row["endpoint_family"],
        "context": context,
    }


def rectangle_count(target_ligands: dict[str, set[str]]) -> tuple[int, int]:
    targets = sorted(target_ligands)
    pairs = rectangles = 0
    for index, left in enumerate(targets):
        for right in targets[index + 1:]:
            common = len(target_ligands[left] & target_ligands[right])
            if common >= 2:
                pairs += 1
                rectangles += common * (common - 1) // 2
    return pairs, rectangles


class UnionFind:
    def __init__(self):
        self.parent = {}

    def add(self, value: str) -> None:
        self.parent.setdefault(value, value)

    def find(self, value: str) -> str:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: str, right: str) -> None:
        self.add(left)
        self.add(right)
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[max(a, b)] = min(a, b)
