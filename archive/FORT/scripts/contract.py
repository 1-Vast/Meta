"""Canonical SIMA-DTA data and episode contracts."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import json
from typing import Sequence


@dataclass(frozen=True)
class AffinityRow:
    """One canonical affinity row. Values may remain unavailable during FSA-D0."""

    target_key: str
    ligand_parent_key: str
    scaffold_key: str
    endpoint: str
    assay_key: str
    document_or_provenance_key: str
    affinity_value: float | None
    split_role: str

    def metakey(self) -> tuple[str, str, str, str, str, str, str]:
        return (
            self.target_key,
            self.ligand_parent_key,
            self.scaffold_key,
            self.endpoint,
            self.assay_key,
            self.document_or_provenance_key,
            self.split_role,
        )


@dataclass(frozen=True)
class Episode:
    """A role-closed support/query assignment shared by every model arm."""

    target_key: str
    support_indices: tuple[int, ...]
    query_indices: tuple[int, ...]
    support_size: int
    homology_component: str
    support_scaffolds: tuple[str, ...]
    query_scaffolds: tuple[str, ...]
    provenance_components: tuple[str, ...]
    episode_hash: str

    @classmethod
    def create(
        cls,
        *,
        target_key: str,
        support_indices: Sequence[int],
        query_indices: Sequence[int],
        homology_component: str,
        support_scaffolds: Sequence[str],
        query_scaffolds: Sequence[str],
        provenance_components: Sequence[str],
    ) -> "Episode":
        support = tuple(sorted(support_indices))
        query = tuple(sorted(query_indices))
        if not support or not query:
            raise ValueError("episodes require nonempty support and query sets")
        if set(support).intersection(query):
            raise ValueError("support and query indices must be disjoint")
        if set(support_scaffolds).intersection(query_scaffolds):
            raise ValueError("support/query scaffold closure is violated")
        payload = {
            "target_key": target_key,
            "support_indices": support,
            "query_indices": query,
            "support_size": len(support),
            "homology_component": homology_component,
            "support_scaffolds": tuple(sorted(set(support_scaffolds))),
            "query_scaffolds": tuple(sorted(set(query_scaffolds))),
            "provenance_components": tuple(sorted(set(provenance_components))),
        }
        digest = sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        return cls(episode_hash=digest, **payload)

    def assertrows(self, rows: Sequence[AffinityRow]) -> None:
        indices = self.support_indices + self.query_indices
        if any(index < 0 or index >= len(rows) for index in indices):
            raise ValueError("episode index is outside the row registry")
        if any(rows[index].target_key != self.target_key for index in indices):
            raise ValueError("episode contains more than one target")
        endpoints = {rows[index].endpoint for index in indices}
        if len(endpoints) != 1:
            raise ValueError("episode mixes endpoint strata")
        support_scaffolds = {rows[index].scaffold_key for index in self.support_indices}
        query_scaffolds = {rows[index].scaffold_key for index in self.query_indices}
        if support_scaffolds.intersection(query_scaffolds):
            raise ValueError("episode scaffold closure is violated")
