"""Build a balanced, label-free source lock for the A2S information gate.

The v1 lock closed target roles by unioning every target connected through a
document token. One resulting component dominated the OOF diagnostic. This
version assigns homology components first, then quarantines rows carrying a
document or assay token that crosses roles. The retained rows therefore keep
the provenance firewall without collapsing unrelated targets into one fold.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
from typing import Iterable

import pandas as pd

from research.a2s.a2s_source_lock import canonical, sha256_file, stable_key, tokens


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "dataset" / "public" / "chembl_37" / "processed" / "dualcold"
REGISTRY = DATA / "registry.parquet"
FEATURES = DATA / "ligand_features.npz"
PROTEINS = DATA / "target_esm2.npz"
DEFAULT_OUTPUT = (
    ROOT / "reports" / "active" / "a2s_source_information_gate_lock_v2_2026-08-01.json"
)

SEED = 20260801
MODEL_FAMILY = "a2s_information_gate_v2"
ROLES = ("fit", "probe", "locked")
ROLE_FRACTIONS = {"fit": 0.50, "probe": 0.25, "locked": 0.25}
PROVENANCE_FIELDS = ("docs", "assays")
N_OOF_FOLDS = 5
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


def load_source_metadata() -> pd.DataFrame:
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
    if frame[list(REQUIRED_COLUMNS)].isna().any().any():
        raise ValueError("source metadata contains null required values")
    for column in ("target", "hcluster", *PROVENANCE_FIELDS):
        frame[column] = frame[column].astype(str)
    if int(frame.groupby("target").hcluster.nunique().max()) != 1:
        raise RuntimeError("a target belongs to more than one homology component")
    return frame


def component_records(frame: pd.DataFrame) -> list[dict[str, object]]:
    return [
        {
            "component": str(hcluster),
            "rows": int(len(group)),
            "targets": int(group.target.nunique()),
        }
        for hcluster, group in frame.groupby("hcluster", sort=True)
    ]


def assign_roles(
    records: Iterable[dict[str, object]], *, seed: int = SEED
) -> dict[str, str]:
    records = list(records)
    total_rows = sum(int(record["rows"]) for record in records)
    total_targets = sum(int(record["targets"]) for record in records)
    if len(records) < len(ROLES):
        raise ValueError("too few homology components for all source roles")
    counts = {
        role: {"rows": 0, "targets": 0, "components": 0} for role in ROLES
    }
    ordered = sorted(
        records,
        key=lambda record: (
            -int(record["rows"]),
            stable_key(seed, str(record["component"])),
        ),
    )
    assignment: dict[str, str] = {}
    for record in ordered:
        def deficit(role: str) -> tuple[float, float, int]:
            row_goal = ROLE_FRACTIONS[role] * total_rows
            target_goal = ROLE_FRACTIONS[role] * total_targets
            return (
                (row_goal - counts[role]["rows"]) / row_goal,
                (target_goal - counts[role]["targets"]) / target_goal,
                -counts[role]["components"],
            )

        role = max(ROLES, key=lambda value: (deficit(value), value))
        component = str(record["component"])
        assignment[component] = role
        counts[role]["rows"] += int(record["rows"])
        counts[role]["targets"] += int(record["targets"])
        counts[role]["components"] += 1
    if set(assignment.values()) != set(ROLES):
        raise RuntimeError("role assignment left a source role empty")
    return assignment


def quarantine_cross_role_rows(
    frame: pd.DataFrame, *, fields: tuple[str, ...] = PROVENANCE_FIELDS
) -> tuple[pd.DataFrame, dict[str, set[str]], dict[str, int]]:
    crossing: dict[str, set[str]] = {}
    for field in fields:
        token_roles: dict[str, set[str]] = defaultdict(set)
        for role, value in frame[["role", field]].drop_duplicates().itertuples(
            index=False
        ):
            for token in tokens(value):
                token_roles[token].add(str(role))
        crossing[field] = {
            token for token, roles in token_roles.items() if len(roles) > 1
        }

    reason_counts = {field: 0 for field in fields}
    keep: list[bool] = []
    for row in frame[list(fields)].itertuples(index=False, name=None):
        row_crossing = []
        for field, value in zip(fields, row):
            hit = bool(tokens(value) & crossing[field])
            row_crossing.append(hit)
            reason_counts[field] += int(hit)
        keep.append(not any(row_crossing))
    retained = frame.loc[pd.Series(keep, index=frame.index)].copy()
    return retained, crossing, reason_counts


def assign_oof_folds(
    retained: pd.DataFrame,
    *,
    n_folds: int = N_OOF_FOLDS,
    seed: int = SEED,
) -> tuple[dict[str, int], list[int]]:
    sizes = retained.loc[retained.role == "fit"].groupby("component").size().to_dict()
    if len(sizes) < n_folds:
        raise ValueError("too few retained fit components for OOF folds")
    loads = [0 for _ in range(n_folds)]
    assignment: dict[str, int] = {}
    ordered = sorted(
        sizes.items(),
        key=lambda item: (-int(item[1]), stable_key(seed, str(item[0]))),
    )
    for component, rows in ordered:
        fold = min(range(n_folds), key=lambda value: (loads[value], value))
        assignment[str(component)] = fold
        loads[fold] += int(rows)
    return assignment, loads


def overlap_audit(frame: pd.DataFrame) -> dict[str, object]:
    overlap: dict[str, dict[str, int]] = {}
    for left_index, left in enumerate(ROLES):
        for right in ROLES[left_index + 1 :]:
            left_rows = frame.loc[frame.role == left]
            right_rows = frame.loc[frame.role == right]
            overlap[f"{left}__{right}"] = {
                "target": len(set(left_rows.target) & set(right_rows.target)),
                "hcluster": len(set(left_rows.hcluster) & set(right_rows.hcluster)),
                "docs": len(
                    set().union(*(tokens(value) for value in left_rows.docs))
                    & set().union(*(tokens(value) for value in right_rows.docs))
                ),
                "assays": len(
                    set().union(*(tokens(value) for value in left_rows.assays))
                    & set().union(*(tokens(value) for value in right_rows.assays))
                ),
                "conn": len(set(left_rows.conn) & set(right_rows.conn)),
                "scaffold": len(set(left_rows.scaffold) & set(right_rows.scaffold)),
            }
    mandatory = ("target", "hcluster", "docs", "assays")
    all_mandatory_zero = all(
        axes[axis] == 0 for axes in overlap.values() for axis in mandatory
    )
    return {
        "overlap": overlap,
        "mandatory_zero_axes": list(mandatory),
        "all_mandatory_zero": all_mandatory_zero,
    }


def role_summary(frame: pd.DataFrame) -> dict[str, object]:
    summary: dict[str, object] = {}
    for role in ROLES:
        active = frame.loc[frame.role == role]
        target_sizes = active.groupby("target").size()
        summary[role] = {
            "rows": int(len(active)),
            "targets": int(active.target.nunique()),
            "components": int(active.component.nunique()),
            "eligible_targets": {
                f"k{k}": int((target_sizes >= k + 5).sum()) for k in (1, 3, 5)
            },
        }
    return summary


def build(output: Path) -> dict[str, object]:
    frame = load_source_metadata()
    records = component_records(frame)
    role_by_component = assign_roles(records)
    frame["component"] = frame.hcluster.astype(str)
    frame["role"] = frame.component.map(role_by_component)
    retained, crossing, reason_counts = quarantine_cross_role_rows(frame)
    checks = overlap_audit(retained)
    if not checks["all_mandatory_zero"]:
        raise RuntimeError("retained source roles violate the provenance firewall")
    if set(retained.role) != set(ROLES):
        raise RuntimeError("quarantine emptied a source role")
    oof_assignment, oof_loads = assign_oof_folds(retained)

    target_component = (
        frame[["target", "component"]]
        .drop_duplicates()
        .set_index("target")
        .component.astype(str)
        .to_dict()
    )
    retained_rows = sorted(retained.source_row.astype(int).tolist())
    result: dict[str, object] = {
        "schema": "a2s-source-information-gate-lock-v2",
        "status": "METADATA_PREFLIGHT_ONLY",
        "model_family": MODEL_FAMILY,
        "seed": SEED,
        "source_policy": {
            "endpoint": "pKi",
            "dual_cold_split": "train",
            "labels_read": False,
            "requested_columns": ["source_row", *REQUIRED_COLUMNS],
            "role_assignment": "homology_component_first_weighted_greedy",
            "role_fractions": ROLE_FRACTIONS,
            "cross_role_row_quarantine": list(PROVENANCE_FIELDS),
            "roles": list(ROLES),
        },
        "inputs": {
            "registry": {"path": str(REGISTRY), "sha256": sha256_file(REGISTRY)},
            "features": {"path": str(FEATURES), "sha256": sha256_file(FEATURES)},
            "proteins": {"path": str(PROTEINS), "sha256": sha256_file(PROTEINS)},
        },
        "counts": {
            "source_rows_before_quarantine": int(len(frame)),
            "source_rows_retained": int(len(retained)),
            "source_rows_quarantined": int(len(frame) - len(retained)),
            "retained_fraction": float(len(retained) / len(frame)),
            "roles_before_quarantine": role_summary(frame),
            "roles_retained": role_summary(retained),
        },
        "components": {
            "assignment": target_component,
            "role_by_component": role_by_component,
            "count": len(records),
            "largest_targets": max(int(record["targets"]) for record in records),
            "largest_rows_before_quarantine": max(int(record["rows"]) for record in records),
        },
        "row_policy": {
            "retained_source_rows": retained_rows,
            "retained_source_rows_sha256": sha256(canonical(retained_rows).encode("utf-8")).hexdigest(),
            "cross_role_token_counts": {
                field: len(values) for field, values in crossing.items()
            },
            "quarantine_reason_row_counts_nonexclusive": reason_counts,
        },
        "oof": {
            "n_folds": N_OOF_FOLDS,
            "role_by_component": oof_assignment,
            "fold_rows": oof_loads,
            "max_holdout_share": max(oof_loads) / sum(oof_loads),
        },
        "checks": checks,
        "interpretation": {
            "fact": "The lock is label-free and quarantines every retained row carrying a provenance token shared across roles.",
            "inference": "The balanced homology OOF design repairs the v1 power collapse but changes the source estimand by excluding cross-role provenance rows.",
            "hypothesis": "A repeated G0/G1 diagnostic can now test support-label information without the v1 97-percent OOF fold imbalance.",
        },
    }
    result["content_sha256"] = sha256(canonical(result).encode("utf-8")).hexdigest()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(result, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = build(args.out.resolve())
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.out.resolve()),
                "content_sha256": result["content_sha256"],
                "counts": result["counts"],
                "oof": result["oof"],
                "checks": result["checks"],
            },
            indent=2,
            sort_keys=True,
            ensure_ascii=True,
        )
    )


if __name__ == "__main__":
    main()
