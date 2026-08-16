"""Freeze a new metadata-only source lock for the A2S information gate.

This command deliberately does not request the affinity column.  It creates a
model-family-specific, provenance-closed component assignment that is distinct
from the previously inspected CMAL validation/meta-test split.  The emitted
JSON is a preflight certificate only; it is not an experiment result.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "dataset" / "public" / "chembl_37" / "processed" / "dualcold"
REGISTRY = DATA / "registry.parquet"
FEATURES = DATA / "ligand_features.npz"
PROTEINS = DATA / "target_esm2.npz"
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_source_information_gate_lock_2026-08-01.json"

SEED = 20260801
MODEL_FAMILY = "a2s_information_gate_v1"
ROLES = ("fit", "probe", "locked")
ROLE_FRACTIONS = {"fit": 0.60, "probe": 0.20, "locked": 0.20}
PROVENANCE_FIELDS = ("docs", "assays")
REQUIRED_COLUMNS = (
    "target",
    "conn",
    "scaffold",
    "hcluster",
    "docs",
    "assays",
    "endpoint",
    "dual_cold_split",
)


class UnionFind:
    def __init__(self, items: Iterable[str]) -> None:
        self.parent = {item: item for item in items}
        self.size = {item: 1 for item in self.parent}

    def find(self, item: str) -> str:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            parent = self.parent[item]
            self.parent[item] = root
            item = parent
        return root

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.size[left_root] < self.size[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        self.size[left_root] += self.size[right_root]


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_key(seed: int, value: str) -> str:
    return sha256(f"{seed}:{value}".encode("utf-8")).hexdigest()


def tokens(value: object) -> set[str]:
    return {part.strip() for part in str(value).split("|") if part.strip()}


def load_source_metadata() -> pd.DataFrame:
    # Keep this list label-free.  Adding affinity here invalidates the lock.
    frame = pd.read_parquet(REGISTRY, columns=list(REQUIRED_COLUMNS)).reset_index(
        names="source_row"
    )
    if "affinity" in frame.columns:
        raise AssertionError("metadata lock unexpectedly loaded affinity")
    frame = frame[
        (frame.endpoint.astype(str) == "pKi")
        & (frame.dual_cold_split.astype(str) == "train")
    ].copy()
    if frame.empty:
        raise RuntimeError("the TRAIN pKi source role is empty")
    for column in REQUIRED_COLUMNS:
        if frame[column].isna().any():
            raise ValueError(f"source metadata contains null {column}")
    frame["target"] = frame.target.astype(str)
    frame["hcluster"] = frame.hcluster.astype(str)
    for column in PROVENANCE_FIELDS:
        frame[column] = frame[column].astype(str)
    if frame.target.duplicated().any() and not frame.target.is_unique:
        # Duplicated rows are expected; target IDs themselves are the graph nodes.
        pass
    return frame


def build_components(frame: pd.DataFrame) -> tuple[dict[str, str], dict[str, object]]:
    targets = sorted(frame.target.unique())
    graph = UnionFind(targets)

    # Homology closure is mandatory.  Provenance closure uses every token in
    # the pipe-delimited document and assay fields, so no token can cross roles.
    for _, group in frame.groupby("hcluster", sort=True):
        members = sorted(group.target.unique())
        for member in members[1:]:
            graph.union(members[0], member)
    provenance_edges = {field: 0 for field in PROVENANCE_FIELDS}
    for field in PROVENANCE_FIELDS:
        token_targets: dict[str, set[str]] = defaultdict(set)
        for target, value in frame[["target", field]].drop_duplicates().itertuples(
            index=False
        ):
            for token in tokens(value):
                token_targets[token].add(str(target))
        for members in token_targets.values():
            ordered = sorted(members)
            if len(ordered) > 1:
                provenance_edges[field] += len(ordered) - 1
                for member in ordered[1:]:
                    graph.union(ordered[0], member)

    groups: dict[str, list[str]] = defaultdict(list)
    for target in targets:
        groups[graph.find(target)].append(target)
    normalized = {
        f"C{position:04d}": sorted(members)
        for position, members in enumerate(
            sorted(groups.values(), key=lambda values: (min(values), len(values)))
        )
    }
    assignment = {
        target: component
        for component, members in normalized.items()
        for target in members
    }
    stats = {
        "targets": len(targets),
        "components": len(normalized),
        "largest_component_targets": max(map(len, normalized.values())),
        "largest_component_share": max(map(len, normalized.values())) / len(targets),
        "singleton_components": sum(len(members) == 1 for members in normalized.values()),
        "provenance_union_edges": provenance_edges,
        "component_assignment_sha256": sha256(
            canonical(sorted(assignment.items())).encode("utf-8")
        ).hexdigest(),
    }
    return assignment, {"members": normalized, "stats": stats}


def assign_roles(components: dict[str, list[str]]) -> dict[str, str]:
    ordered = sorted(components, key=lambda value: stable_key(SEED, value))
    count = len(ordered)
    fit_count = max(1, int(np.floor(ROLE_FRACTIONS["fit"] * count)))
    probe_count = max(1, int(np.floor(ROLE_FRACTIONS["probe"] * count)))
    if fit_count + probe_count >= count:
        fit_count = max(1, count - 2)
        probe_count = 1
    role_by_component: dict[str, str] = {}
    for position, component in enumerate(ordered):
        role_by_component[component] = (
            "fit" if position < fit_count else "probe" if position < fit_count + probe_count else "locked"
        )
    return role_by_component


def feasibility(frame: pd.DataFrame, target_role: dict[str, str], target_component: dict[str, str]) -> dict[str, object]:
    rows: list[dict[str, object]] = []
    for target, group in frame.groupby("target", sort=True):
        role = target_role[target]
        rows.append(
            {
                "target": str(target),
                "component": target_component[target],
                "role": role,
                "rows": int(len(group)),
                "unique_conn": int(group.conn.nunique()),
                "unique_scaffold": int(group.scaffold.nunique()),
                "unique_docs": int(group.docs.nunique()),
                "unique_assays": int(group.assays.nunique()),
                "eligible_k1": bool(len(group) >= 6),
                "eligible_k3": bool(len(group) >= 8),
                "eligible_k5": bool(len(group) >= 10),
            }
        )
    table = pd.DataFrame(rows)
    by_role: dict[str, object] = {}
    for role in ROLES:
        active = table[table.role == role]
        by_role[role] = {
            "targets": int(len(active)),
            "components": int(active.component.nunique()),
            "rows": int(active.rows.sum()),
            "eligible_targets": {
                f"k{k}": int(active[f"eligible_k{k}"].sum()) for k in (1, 3, 5)
            },
            "median_rows_per_target": float(active.rows.median()) if len(active) else None,
        }
    return {"by_role": by_role, "target_records": rows}


def verify(frame: pd.DataFrame, target_component: dict[str, str], role_by_component: dict[str, str]) -> dict[str, object]:
    frame_component = frame.target.map(target_component)
    frame_role = frame_component.map(role_by_component)
    if frame_component.isna().any() or frame_role.isna().any():
        raise RuntimeError("some source row lacks a component or role")
    overlap: dict[str, dict[str, int]] = {}
    for left in ROLES:
        for right in ROLES:
            if left >= right:
                continue
            left_rows = frame.loc[frame_role == left]
            right_rows = frame.loc[frame_role == right]
            overlap[f"{left}__{right}"] = {
                "target": len(set(left_rows.target).intersection(right_rows.target)),
                "hcluster": len(set(left_rows.hcluster).intersection(right_rows.hcluster)),
                "docs": len(
                    set().union(*(tokens(value) for value in left_rows.docs))
                    & set().union(*(tokens(value) for value in right_rows.docs))
                ),
                "assays": len(
                    set().union(*(tokens(value) for value in left_rows.assays))
                    & set().union(*(tokens(value) for value in right_rows.assays))
                ),
            }
    if any(any(value for value in axes.values()) for axes in overlap.values()):
        raise RuntimeError(f"component/provenance roles overlap: {overlap}")
    return {"overlap": overlap, "all_zero": True}


def build(output: Path) -> dict[str, object]:
    frame = load_source_metadata()
    target_component, component_data = build_components(frame)
    role_by_component = assign_roles(component_data["members"])
    target_role = {target: role_by_component[component] for target, component in target_component.items()}
    checks = verify(frame, target_component, role_by_component)
    role_counts = {
        role: {
            "components": sum(value == role for value in role_by_component.values()),
            "targets": sum(target_role[target] == role for target in target_role),
            "rows": int(sum(target_role[target] == role for target in frame.target)),
        }
        for role in ROLES
    }
    result: dict[str, object] = {
        "schema": "a2s-source-information-gate-lock-v1",
        "status": "METADATA_PREFLIGHT_ONLY",
        "model_family": MODEL_FAMILY,
        "seed": SEED,
        "source_policy": {
            "endpoint": "pKi",
            "dual_cold_split": "train",
            "labels_read": False,
            "requested_columns": ["source_row", *REQUIRED_COLUMNS],
            "component_closure": ["target_hcluster", "pipe_token_docs", "pipe_token_assays"],
            "role_fractions": ROLE_FRACTIONS,
            "roles": ROLES,
        },
        "inputs": {
            "registry": {"path": str(REGISTRY), "sha256": sha256_file(REGISTRY)},
            "features": {"path": str(FEATURES), "sha256": sha256_file(FEATURES)},
            "proteins": {"path": str(PROTEINS), "sha256": sha256_file(PROTEINS)},
        },
        "counts": {
            "source_rows": int(len(frame)),
            "source_targets": int(frame.target.nunique()),
            "source_hclusters": int(frame.hcluster.nunique()),
            "source_doc_tokens": int(len(set().union(*(tokens(value) for value in frame.docs)))),
            "source_assay_tokens": int(len(set().union(*(tokens(value) for value in frame.assays)))),
            "roles": role_counts,
        },
        "components": {
            "assignment": target_component,
            "members": component_data["members"],
            "role_by_component": role_by_component,
            "stats": component_data["stats"],
        },
        "feasibility": feasibility(frame, target_role, target_component),
        "checks": checks,
        "interpretation": {
            "fact": "This certificate uses TRAIN metadata only and does not establish support-label information.",
            "inference": "Document-token closure is a conservative provenance proxy; its large components may reduce independent power.",
            "hypothesis": "A complete OOF base and G0/G1 probes may show exploitable support-label information on the locked roles.",
        },
    }
    result["content_sha256"] = sha256(canonical(result).encode("utf-8")).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, ensure_ascii=True) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = build(args.out.resolve())
    print(json.dumps({
        "status": result["status"],
        "output": str(args.out.resolve()),
        "content_sha256": result["content_sha256"],
        "counts": result["counts"],
        "component_stats": result["components"]["stats"],
        "checks": result["checks"],
    }, indent=2, sort_keys=True, ensure_ascii=True))


if __name__ == "__main__":
    main()
