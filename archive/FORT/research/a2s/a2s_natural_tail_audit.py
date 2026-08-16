"""Metadata-only Gate D0 audit for the A2S-DTA natural pKi tail."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


ROOT = Path("dataset/public/chembl_37/processed/dualcold")
REGISTRY = ROOT / "registry.parquet"
DOCUMENT_METADATA = ROOT / "pcic_o0_document_metadata.json"
DEFAULT_OUTPUT = Path("dataset/processed/a2s_natural_tail_d0.v1.json")
ENDPOINT = "pKi"
SPLIT = "train"
SOURCE_MIN = 100
RECIPIENT_MAX = 30
SUPPORT_KS = (1, 3, 5)
SUPPORT_DRAWS = 5
MIN_QUERY = 10
SEED = 1729
PLANNING_PAIRED_SD = 0.10
MATERIAL_FLOOR = 0.05
NORMAL_95_PLUS_80 = 1.959963984540054 + 0.8416212335729143
METADATA_COLUMNS = (
    "target",
    "conn",
    "endpoint",
    "scaffold",
    "assays",
    "docs",
    "hcluster",
    "dual_cold_split",
)


@dataclass(frozen=True)
class Unit:
    row: int
    target: str
    conn: str
    scaffold: str
    docs: frozenset[str]
    assays: frozenset[str]
    sources: frozenset[str]
    release_min: int
    release_max: int
    homology: str


@dataclass(frozen=True)
class Episode:
    target: str
    homology: str
    release_cut: int
    query_source: str
    draw: int
    support: tuple[int, ...]
    query: tuple[int, ...]


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tokens(value: object) -> frozenset[str]:
    result = frozenset(part.strip() for part in str(value).split("|") if part.strip())
    if not result:
        raise ValueError("document, assay, and source token sets must not be empty")
    return result


def release_number(value: object) -> int:
    prefix, separator, suffix = str(value).rpartition("_")
    if prefix != "CHEMBL" or not separator or not suffix.isdigit():
        raise ValueError(f"invalid ChEMBL release identifier: {value}")
    return int(suffix)


def load_document_metadata(path: Path) -> dict[str, dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    records = payload.get("documents")
    if not isinstance(records, list) or not records:
        raise ValueError("document metadata has no documents")
    result: dict[str, dict[str, Any]] = {}
    for record in records:
        document = str(record["document_chembl_id"])
        if document in result:
            raise ValueError(f"duplicate document metadata: {document}")
        if record.get("src_id") is None:
            raise ValueError(f"document has no source lineage: {document}")
        release_number(record["chembl_release"])
        result[document] = record
    return result


def effective_units(frame: pd.DataFrame) -> pd.DataFrame:
    missing = sorted(set(METADATA_COLUMNS).difference(frame.columns))
    if missing:
        raise ValueError(f"registry is missing D0 metadata columns: {missing}")
    if "source_row" not in frame.columns:
        frame = frame.reset_index(names="source_row")
    else:
        frame = frame.copy()
    frame = frame[
        (frame.endpoint.astype(str) == ENDPOINT)
        & (frame.dual_cold_split.astype(str) == SPLIT)
    ].copy()
    if frame.empty:
        raise ValueError("D0 found no TRAIN pKi metadata rows")
    frame["source_row"] = pd.to_numeric(frame.source_row, errors="raise").astype("int64")
    frame["target"] = frame.target.astype(str)
    keys = ["target", "endpoint", "conn", "docs", "assays"]
    return (
        frame.sort_values(keys + ["source_row"], kind="stable")
        .drop_duplicates(keys, keep="first")
        .reset_index(drop=True)
    )


def decorate_units(
    group: pd.DataFrame, documents: dict[str, dict[str, Any]]
) -> list[Unit]:
    result: list[Unit] = []
    for row in group.itertuples(index=False):
        document_tokens = tokens(row.docs)
        missing = sorted(document_tokens.difference(documents))
        if missing:
            raise ValueError(f"registry documents lack metadata: {missing[:3]}")
        records = [documents[document] for document in document_tokens]
        releases = [release_number(record["chembl_release"]) for record in records]
        sources = frozenset(str(record["src_id"]) for record in records)
        result.append(
            Unit(
                row=int(row.source_row),
                target=str(row.target),
                conn=str(row.conn),
                scaffold=str(row.scaffold),
                docs=document_tokens,
                assays=tokens(row.assays),
                sources=sources,
                release_min=min(releases),
                release_max=max(releases),
                homology=str(row.hcluster),
            )
        )
    return result


def _hash_order(target: str, draw: int, row: int) -> str:
    return sha256(f"{SEED}|{target}|{draw}|{row}".encode("ascii")).hexdigest()


def select_support(
    pool: Iterable[Unit], target: str, draw: int, *, distinct_scaffold: bool = True
) -> tuple[Unit, ...]:
    ordered = sorted(pool, key=lambda unit: _hash_order(target, draw, unit.row))
    chosen: list[Unit] = []
    parents: set[str] = set()
    scaffolds: set[str] = set()
    for unit in ordered:
        if unit.conn in parents or (distinct_scaffold and unit.scaffold in scaffolds):
            continue
        chosen.append(unit)
        parents.add(unit.conn)
        scaffolds.add(unit.scaffold)
        if len(chosen) == max(SUPPORT_KS):
            return tuple(chosen)
    return ()


def _axis_values(unit: Unit, axis: str) -> set[str]:
    if axis == "parent":
        return {unit.conn}
    if axis == "scaffold":
        return {unit.scaffold}
    if axis == "document":
        return set(unit.docs)
    if axis == "assay":
        return set(unit.assays)
    if axis == "source_family":
        return set(unit.sources)
    raise ValueError(f"unknown closure axis: {axis}")


def closed_query(
    pool: Iterable[Unit],
    support: tuple[Unit, ...],
    *,
    axes: tuple[str, ...] = (
        "parent",
        "scaffold",
        "document",
        "assay",
        "source_family",
    ),
) -> tuple[Unit, ...]:
    used = {
        axis: set().union(*(_axis_values(unit, axis) for unit in support))
        for axis in axes
    }
    return tuple(
        unit
        for unit in pool
        if all(_axis_values(unit, axis).isdisjoint(used[axis]) for axis in axes)
    )


def _evaluate_split(
    support_pool: list[Unit],
    query_pool: list[Unit],
    *,
    target: str,
    homology: str,
    release_cut: int,
    query_source: str,
    closure_axes: tuple[str, ...] = (
        "parent",
        "scaffold",
        "document",
        "assay",
        "source_family",
    ),
    distinct_scaffold: bool = True,
) -> list[Episode]:
    episodes: list[Episode] = []
    for draw in range(SUPPORT_DRAWS):
        support = select_support(
            support_pool, target, draw, distinct_scaffold=distinct_scaffold
        )
        if len(support) != max(SUPPORT_KS):
            return []
        query = closed_query(query_pool, support, axes=closure_axes)
        episodes.append(
            Episode(
                target=target,
                homology=homology,
                release_cut=release_cut,
                query_source=query_source,
                draw=draw,
                support=tuple(unit.row for unit in support),
                query=tuple(unit.row for unit in query),
            )
        )
    return episodes


def _candidate_score(episodes: list[Episode]) -> tuple[int, int]:
    depths = [len(episode.query) for episode in episodes]
    return min(depths), sum(depths)


def _strict_candidates(units: list[Unit]) -> list[tuple[list[Episode], int, str, int, int]]:
    target = units[0].target
    homology = units[0].homology
    releases = sorted(
        {unit.release_min for unit in units}.union(unit.release_max for unit in units)
    )
    sources = sorted(set().union(*(unit.sources for unit in units)))
    candidates: list[tuple[list[Episode], int, str, int, int]] = []
    for query_source in sources:
        singleton_source = frozenset((query_source,))
        for release_cut in releases:
            support_pool = [
                unit
                for unit in units
                if unit.release_max <= release_cut and query_source not in unit.sources
            ]
            query_pool = [
                unit
                for unit in units
                if unit.release_min > release_cut and unit.sources == singleton_source
            ]
            if len(support_pool) < max(SUPPORT_KS) or len(query_pool) < MIN_QUERY:
                continue
            episodes = _evaluate_split(
                support_pool,
                query_pool,
                target=target,
                homology=homology,
                release_cut=release_cut,
                query_source=query_source,
            )
            if episodes:
                candidates.append(
                    (episodes, release_cut, query_source, len(support_pool), len(query_pool))
                )
    return candidates


def _best_strict_candidate(
    units: list[Unit],
) -> tuple[list[Episode], int, str, int, int] | None:
    candidates = _strict_candidates(units)
    if not candidates:
        return None
    return max(
        candidates,
        key=lambda item: (
            *_candidate_score(item[0]),
            -item[1],
            item[2],
        ),
    )


def _envelope_flags(units: list[Unit]) -> dict[str, bool]:
    releases = sorted(
        {unit.release_min for unit in units}.union(unit.release_max for unit in units)
    )
    sources = sorted(set().union(*(unit.sources for unit in units)))
    raw = len(units) >= max(SUPPORT_KS) + MIN_QUERY
    temporal = any(
        sum(unit.release_max <= cut for unit in units) >= max(SUPPORT_KS)
        and sum(unit.release_min > cut for unit in units) >= MIN_QUERY
        for cut in releases
    )
    source = any(
        sum(query_source not in unit.sources for unit in units) >= max(SUPPORT_KS)
        and sum(unit.sources == frozenset((query_source,)) for unit in units) >= MIN_QUERY
        for query_source in sources
    )
    temporal_source = any(
        sum(
            unit.release_max <= cut and query_source not in unit.sources
            for unit in units
        )
        >= max(SUPPORT_KS)
        and sum(
            unit.release_min > cut and unit.sources == frozenset((query_source,))
            for unit in units
        )
        >= MIN_QUERY
        for query_source in sources
        for cut in releases
    )
    return {
        "raw_k5_plus_query_envelope": raw,
        "temporal_envelope": temporal,
        "source_family_envelope": source,
        "temporal_source_envelope": temporal_source,
    }


def _optimistic_provenance_upper_bound(
    units: list[Unit], *, distinct_support_scaffold: bool
) -> bool:
    """Return whether five support combinations survive minimal provenance closure."""

    valid_draws = 0
    for support_indices in combinations(range(len(units)), max(SUPPORT_KS)):
        support = tuple(units[index] for index in support_indices)
        if distinct_support_scaffold and len({unit.scaffold for unit in support}) < 5:
            continue
        query = closed_query(
            units,
            support,
            axes=("parent", "document", "assay"),
        )
        if len(query) >= MIN_QUERY:
            valid_draws += 1
            if valid_draws >= SUPPORT_DRAWS:
                return True
    return False


def _full_temporal(
    units: list[Unit],
    *,
    closure_axes: tuple[str, ...],
    distinct_scaffold: bool,
) -> bool:
    target = units[0].target
    homology = units[0].homology
    releases = sorted(
        {unit.release_min for unit in units}.union(unit.release_max for unit in units)
    )
    for release_cut in releases:
        support_pool = [unit for unit in units if unit.release_max <= release_cut]
        query_pool = [unit for unit in units if unit.release_min > release_cut]
        episodes = _evaluate_split(
            support_pool,
            query_pool,
            target=target,
            homology=homology,
            release_cut=release_cut,
            query_source="TEMPORAL_ONLY_DIAGNOSTIC",
            closure_axes=closure_axes,
            distinct_scaffold=distinct_scaffold,
        )
        if episodes and min(len(episode.query) for episode in episodes) >= MIN_QUERY:
            return True
    return False


def _full_source_only(units: list[Unit]) -> bool:
    target = units[0].target
    homology = units[0].homology
    sources = sorted(set().union(*(unit.sources for unit in units)))
    for query_source in sources:
        support_pool = [unit for unit in units if query_source not in unit.sources]
        query_pool = [
            unit for unit in units if unit.sources == frozenset((query_source,))
        ]
        episodes = _evaluate_split(
            support_pool,
            query_pool,
            target=target,
            homology=homology,
            release_cut=0,
            query_source=query_source,
        )
        if episodes and min(len(episode.query) for episode in episodes) >= MIN_QUERY:
            return True
    return False


def _dependency_components(records: list[dict[str, Any]]) -> list[list[str]]:
    parent = {record["target"]: record["target"] for record in records}

    def find(value: str) -> str:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: str, right: str) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    owners: dict[tuple[str, str], str] = {}
    for record in records:
        target = record["target"]
        keys = (("homology", record["homology"]), ("source", record["query_source"]))
        for key in keys:
            if key in owners:
                union(target, owners[key])
            else:
                owners[key] = target
    components: dict[str, list[str]] = {}
    for target in parent:
        components.setdefault(find(target), []).append(target)
    return [sorted(component) for component in components.values()]


def _candidate_dependency_components(
    groups: dict[str, list[Unit]],
) -> list[list[str]]:
    parent = {target: target for target in groups}

    def find(value: str) -> str:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: str, right: str) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    owners: dict[tuple[str, str], str] = {}
    for target, units in groups.items():
        homologies = {unit.homology for unit in units}
        sources = set().union(*(unit.sources for unit in units))
        keys = [("homology", value) for value in homologies]
        keys.extend(("source", value) for value in sources)
        for key in keys:
            if key in owners:
                union(target, owners[key])
            else:
                owners[key] = target
    components: dict[str, list[str]] = {}
    for target in parent:
        components.setdefault(find(target), []).append(target)
    return [sorted(component) for component in components.values()]


def _quantiles(values: Iterable[int]) -> dict[str, float | None]:
    series = pd.Series(list(values), dtype="float64")
    if series.empty:
        return {key: None for key in ("min", "q25", "median", "q75", "q90", "max")}
    return {
        "min": float(series.min()),
        "q25": float(series.quantile(0.25)),
        "median": float(series.median()),
        "q75": float(series.quantile(0.75)),
        "q90": float(series.quantile(0.90)),
        "max": float(series.max()),
    }


def _target_overlap_count(
    groups: dict[str, list[Unit]], source_values: set[str], attribute: str
) -> int:
    count = 0
    for units in groups.values():
        values: set[str] = set()
        for unit in units:
            value = getattr(unit, attribute)
            values.update(value if isinstance(value, frozenset) else (value,))
        count += bool(values.intersection(source_values))
    return count


def build_audit(
    frame: pd.DataFrame,
    documents: dict[str, dict[str, Any]],
    *,
    registry_sha256: str,
    document_metadata_sha256: str,
) -> dict[str, Any]:
    units_frame = effective_units(frame)
    counts = units_frame.groupby("target", sort=True).size().astype(int)
    source_targets = set(counts[counts >= SOURCE_MIN].index.astype(str))
    recipient_targets = set(counts[counts < RECIPIENT_MAX].index.astype(str))
    if source_targets.intersection(recipient_targets):
        raise AssertionError("source and recipient target IDs overlap")

    decorated = {
        str(target): decorate_units(group, documents)
        for target, group in units_frame.groupby("target", sort=True)
    }
    recipient_groups = {
        target: decorated[target] for target in sorted(recipient_targets)
    }
    cutflow = Counter()
    target_diagnostics: list[dict[str, Any]] = []
    roster: list[dict[str, Any]] = []
    for target, target_units in recipient_groups.items():
        flags = _envelope_flags(target_units)
        cutflow.update(key for key, passed in flags.items() if passed)
        provenance_upper_bound = _optimistic_provenance_upper_bound(
            target_units, distinct_support_scaffold=False
        )
        support_scaffold_upper_bound = _optimistic_provenance_upper_bound(
            target_units, distinct_support_scaffold=True
        )
        temporal_provenance_full = _full_temporal(
            target_units,
            closure_axes=("parent", "document", "assay"),
            distinct_scaffold=False,
        )
        temporal_scaffold_full = _full_temporal(
            target_units,
            closure_axes=("parent", "scaffold", "document", "assay"),
            distinct_scaffold=True,
        )
        source_full = _full_source_only(target_units)
        cutflow.update(
            ["optimistic_provenance_upper_bound"] if provenance_upper_bound else []
        )
        cutflow.update(
            ["optimistic_support_scaffold_upper_bound"]
            if support_scaffold_upper_bound
            else []
        )
        cutflow.update(
            ["temporal_provenance_full_closure"] if temporal_provenance_full else []
        )
        cutflow.update(
            ["temporal_scaffold_full_closure"] if temporal_scaffold_full else []
        )
        cutflow.update(["source_full_closure"] if source_full else [])
        best = _best_strict_candidate(target_units)
        admitted = bool(best and min(len(episode.query) for episode in best[0]) >= MIN_QUERY)
        cutflow.update(["strict_admitted"] if admitted else [])
        diagnostic: dict[str, Any] = {
            "target": target,
            "n_eff": len(target_units),
            "homology": target_units[0].homology,
            "source_families": sorted(set().union(*(unit.sources for unit in target_units))),
            **flags,
            "optimistic_provenance_upper_bound": provenance_upper_bound,
            "optimistic_support_scaffold_upper_bound": support_scaffold_upper_bound,
            "temporal_provenance_full_closure": temporal_provenance_full,
            "temporal_scaffold_full_closure": temporal_scaffold_full,
            "source_full_closure": source_full,
            "strict_admitted": admitted,
        }
        if best:
            episodes, release_cut, query_source, support_pool_size, query_pool_size = best
            depths = [len(episode.query) for episode in episodes]
            diagnostic["strict_best"] = {
                "release_cut": release_cut,
                "query_source": query_source,
                "support_pool_units": support_pool_size,
                "late_query_pool_units": query_pool_size,
                "closed_query_depths": depths,
            }
            if admitted:
                episode_records = []
                for episode in episodes:
                    support_by_k = {
                        str(k): list(episode.support[:k]) for k in SUPPORT_KS
                    }
                    episode_payload = {
                        "draw": episode.draw,
                        "support_by_k": support_by_k,
                        "query_rows": list(episode.query),
                    }
                    episode_payload["sha256"] = sha256(
                        json.dumps(
                            episode_payload, sort_keys=True, separators=(",", ":")
                        ).encode("utf-8")
                    ).hexdigest()
                    episode_records.append(episode_payload)
                roster.append(
                    {
                        "target": target,
                        "homology": episodes[0].homology,
                        "release_cut": release_cut,
                        "query_source": query_source,
                        "support_pool_units": support_pool_size,
                        "late_query_pool_units": query_pool_size,
                        "min_closed_query_depth": min(depths),
                        "episodes": episode_records,
                    }
                )
        target_diagnostics.append(diagnostic)

    roster_payload = json.dumps(roster, sort_keys=True, separators=(",", ":")).encode("utf-8")
    components = _dependency_components(roster)
    candidate_components = _candidate_dependency_components(recipient_groups)
    component_count = len(components)
    mde80 = (
        NORMAL_95_PLUS_80 * PLANNING_PAIRED_SD / math.sqrt(component_count)
        if component_count
        else None
    )
    criteria = {
        "at_least_50_recipients": len(roster) >= 50,
        "at_least_10_query_units_all_draws": bool(roster)
        and all(record["min_closed_query_depth"] >= MIN_QUERY for record in roster),
        "component_mde80_at_most_material_floor": mde80 is not None
        and mde80 <= MATERIAL_FLOOR,
        "true_measurement_time_available": False,
    }
    passed = all(criteria.values())

    source_units = [unit for target in source_targets for unit in decorated[target]]
    source_overlap_values = {
        "homology": {unit.homology for unit in source_units},
        "conn": {unit.conn for unit in source_units},
        "scaffold": {unit.scaffold for unit in source_units},
        "docs": set().union(*(unit.docs for unit in source_units)),
        "assays": set().union(*(unit.assays for unit in source_units)),
        "sources": set().union(*(unit.sources for unit in source_units)),
    }
    source_family_target_counts = Counter(
        source
        for units in recipient_groups.values()
        for source in set().union(*(unit.sources for unit in units))
    )
    top_source = source_family_target_counts.most_common(1)
    dependency = {
        "candidate_components": len(candidate_components),
        "candidate_component_sizes": sorted(
            (len(component) for component in candidate_components), reverse=True
        ),
        "candidate_largest_component_fraction": (
            max(map(len, candidate_components)) / len(recipient_groups)
            if candidate_components
            else None
        ),
        "independent_components": component_count,
        "component_sizes": sorted((len(component) for component in components), reverse=True),
        "largest_component_fraction": (
            max(map(len, components)) / len(roster) if roster else None
        ),
        "candidate_top_source_family": top_source[0][0] if top_source else None,
        "candidate_top_source_family_target_fraction": (
            top_source[0][1] / len(recipient_groups) if top_source else None
        ),
    }
    overlap = {
        "source_recipient_target_id_overlap": 0,
        "candidate_recipient_targets_sharing_source_homology": _target_overlap_count(
            recipient_groups, source_overlap_values["homology"], "homology"
        ),
        "candidate_recipient_targets_sharing_source_parent": _target_overlap_count(
            recipient_groups, source_overlap_values["conn"], "conn"
        ),
        "candidate_recipient_targets_sharing_source_scaffold": _target_overlap_count(
            recipient_groups, source_overlap_values["scaffold"], "scaffold"
        ),
        "candidate_recipient_targets_sharing_source_document": _target_overlap_count(
            recipient_groups, source_overlap_values["docs"], "docs"
        ),
        "candidate_recipient_targets_sharing_source_assay": _target_overlap_count(
            recipient_groups, source_overlap_values["assays"], "assays"
        ),
        "candidate_recipient_targets_sharing_source_family": _target_overlap_count(
            recipient_groups, source_overlap_values["sources"], "sources"
        ),
        "admitted_support_query_overlap": {
            "episodes": sum(len(record["episodes"]) for record in roster),
            "row": 0,
            "parent": 0,
            "scaffold": 0,
            "document": 0,
            "assay": 0,
            "source_family": 0,
        },
    }

    return {
        "schema_version": "a2s-natural-tail-d0-v1",
        "program": "A2S-DTA",
        "stage": "GATE_D0_NATURAL_TAIL_CLOSURE_AND_POWER",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "label_policy": "metadata-only; affinity and replicate_sd are not read",
        "input": {
            "registry": str(REGISTRY),
            "registry_sha256": registry_sha256,
            "document_metadata": str(DOCUMENT_METADATA),
            "document_metadata_sha256": document_metadata_sha256,
            "columns_read": ["source_row", *METADATA_COLUMNS],
            "endpoint": ENDPOINT,
            "split": SPLIT,
            "effective_units": len(units_frame),
        },
        "frozen_contract": {
            "source_rule": f"n_eff >= {SOURCE_MIN}",
            "recipient_rule": f"n_eff < {RECIPIENT_MAX}",
            "support_budgets": list(SUPPORT_KS),
            "support_draws": SUPPORT_DRAWS,
            "minimum_query_units": MIN_QUERY,
            "support_selection": "SHA256(seed,target,draw,row), greedy distinct parent and scaffold",
            "time_rule": "all support documents entered ChEMBL by release cut; all query documents entered later",
            "time_metadata_status": "chembl_release is an ingestion proxy; publication and measurement dates are unavailable",
            "source_rule_detail": "query units have one held-out document src_id; support units exclude that src_id",
            "query_closure": ["parent", "scaffold", "document", "assay", "source_family"],
            "common_query_rule": "query is closed against the k=5 support and reused for k=1,3,5",
        },
        "topology": {
            "source_targets": len(source_targets),
            "recipient_candidates": len(recipient_targets),
            "cutflow": {
                "raw_k5_plus_query_envelope": cutflow["raw_k5_plus_query_envelope"],
                "temporal_envelope": cutflow["temporal_envelope"],
                "source_family_envelope": cutflow["source_family_envelope"],
                "temporal_source_envelope": cutflow["temporal_source_envelope"],
                "optimistic_provenance_upper_bound": cutflow[
                    "optimistic_provenance_upper_bound"
                ],
                "optimistic_support_scaffold_upper_bound": cutflow[
                    "optimistic_support_scaffold_upper_bound"
                ],
                "temporal_provenance_full_closure": cutflow[
                    "temporal_provenance_full_closure"
                ],
                "temporal_scaffold_full_closure": cutflow[
                    "temporal_scaffold_full_closure"
                ],
                "source_full_closure": cutflow["source_full_closure"],
                "strict_admitted": cutflow["strict_admitted"],
            },
            "query_depth": _quantiles(
                record["min_closed_query_depth"] for record in roster
            ),
        },
        "overlap": overlap,
        "dependency_concentration": dependency,
        "power": {
            "statistical_unit": "component joining shared homology or held-out source family",
            "planning_paired_sd": PLANNING_PAIRED_SD,
            "two_sided_alpha": 0.05,
            "power": 0.80,
            "mde80": mde80,
            "material_floor": MATERIAL_FLOOR,
            "formula": "(z_0.975 + z_0.80) * planning_paired_sd / sqrt(component_count)",
        },
        "roster": {
            "recipients": len(roster),
            "sha256": sha256(roster_payload).hexdigest(),
            "records": roster,
        },
        "target_diagnostics": target_diagnostics,
        "decision": {
            "status": "READY_FOR_STATIC_PROBE" if passed else "DATA_NOT_READY",
            "pass": passed,
            "criteria": criteria,
            "training_authorized": False,
            "next_action": (
                "run only the preregistered S0 frozen-feature atlas headroom probe"
                if passed
                else "acquire and freeze a provenance-rich pKi natural-tail source with at least 50 closed recipients"
            ),
        },
    }


def run(*, registry: Path, document_metadata: Path, output: Path) -> dict[str, Any]:
    # Read only metadata columns from the full table so source_row remains the
    # canonical global registry row used by the aligned feature cache.
    frame = pd.read_parquet(registry, columns=list(METADATA_COLUMNS)).reset_index(
        names="source_row"
    )
    documents = load_document_metadata(document_metadata)
    report = build_audit(
        frame,
        documents,
        registry_sha256=file_sha256(registry),
        document_metadata_sha256=file_sha256(document_metadata),
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, default=REGISTRY)
    parser.add_argument("--document-metadata", type=Path, default=DOCUMENT_METADATA)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    report = run(
        registry=args.registry,
        document_metadata=args.document_metadata,
        output=args.out,
    )
    print(
        json.dumps(
            {
                "output": str(args.out),
                "status": report["decision"]["status"],
                "cutflow": report["topology"]["cutflow"],
                "components": report["dependency_concentration"]["independent_components"],
                "mde80": report["power"]["mde80"],
            },
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
