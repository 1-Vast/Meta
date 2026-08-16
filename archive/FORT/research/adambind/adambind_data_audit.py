"""Audit AdaMBind data integrity and the target-level fallback split.

This script is intentionally independent of the FORT affinity registry.  It
only reads the externally downloaded AdaMBind CSVs and reports the split
properties that can be established from their three columns.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _split_targets(targets: list[str], seed: int) -> dict[str, str]:
    ordered = np.asarray(sorted(set(targets)), dtype=object)
    rng = np.random.default_rng(seed)
    ordered = ordered[rng.permutation(len(ordered))]
    n = len(ordered)
    n_train = int(n * 0.8)
    n_val = int(n * 0.1)
    labels = {}
    for name, values in (
        ("train", ordered[:n_train]),
        ("val", ordered[n_train : n_train + n_val]),
        ("test", ordered[n_train + n_val :]),
    ):
        labels.update({str(target): name for target in values})
    return labels


def _audit_file(path: Path, seed: int) -> dict:
    frame = pd.read_csv(path)
    required = {"compound_iso_smiles", "target_sequence", "affinity"}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"{path.name}: missing required columns {missing}")

    finite = pd.to_numeric(frame["affinity"], errors="coerce").replace(
        [np.inf, -np.inf], np.nan
    )
    target_labels = _split_targets(frame["target_sequence"].astype(str).tolist(), seed)
    split = frame["target_sequence"].astype(str).map(target_labels)
    compounds_by_split = {
        name: set(frame.loc[split == name, "compound_iso_smiles"].astype(str))
        for name in ("train", "val", "test")
    }
    cross_split_compounds = {
        "train_val": len(compounds_by_split["train"] & compounds_by_split["val"]),
        "train_test": len(compounds_by_split["train"] & compounds_by_split["test"]),
        "val_test": len(compounds_by_split["val"] & compounds_by_split["test"]),
    }

    # This mirrors train.py's per-target support/query construction: shuffle
    # rows inside each target, then reserve the first `nums` rows as support.
    # We report the smallest useful value (5) without using it for training.
    rng = np.random.default_rng(seed)
    support_pairs = set()
    query_pairs = set()
    support_rows = query_rows = 0
    for _, group in frame.groupby("target_sequence", sort=True):
        idx = np.asarray(group.index, dtype=int)
        idx = idx[rng.permutation(len(idx))]
        support_idx = idx[:5]
        query_idx = idx[5:]
        support_rows += len(support_idx)
        query_rows += len(query_idx)
        support_pairs.update(
            (str(frame.at[i, "compound_iso_smiles"]), str(frame.at[i, "target_sequence"]))
            for i in support_idx
        )
        query_pairs.update(
            (str(frame.at[i, "compound_iso_smiles"]), str(frame.at[i, "target_sequence"]))
            for i in query_idx
        )

    return {
        "file": path.name,
        "bytes": path.stat().st_size,
        "sha256": _sha256(path),
        "rows": int(len(frame)),
        "columns": [str(c) for c in frame.columns],
        "targets": int(frame["target_sequence"].nunique()),
        "compounds": int(frame["compound_iso_smiles"].nunique()),
        "compound_target_pairs": int(
            frame[["compound_iso_smiles", "target_sequence"]].drop_duplicates().shape[0]
        ),
        "duplicate_rows": int(frame.duplicated().sum()),
        "duplicate_compound_target_rows": int(
            frame[["compound_iso_smiles", "target_sequence"]].duplicated().sum()
        ),
        "nonfinite_affinity": int(finite.isna().sum()),
        "fallback_split": {
            "seed": seed,
            "target_counts": {
                name: int(sum(value == name for value in target_labels.values()))
                for name in ("train", "val", "test")
            },
            "row_counts": {
                name: int((split == name).sum()) for name in ("train", "val", "test")
            },
            "cross_split_compounds": cross_split_compounds,
            "target_overlap": 0,
        },
        "per_target_support_query_nums5": {
            "support_rows": int(support_rows),
            "query_rows": int(query_rows),
            "support_query_exact_pair_overlap": int(len(support_pairs & query_pairs)),
            "note": "A duplicated compound-target pair can cross support/query; rows are not independent.",
        },
    }


def audit(data_root: Path, seed: int) -> dict:
    files = sorted(data_root.glob("*-full-data.csv"))
    if not files:
        raise FileNotFoundError(f"no *-full-data.csv files under {data_root}")
    split_files = {
        path.stem.replace("-full-data", ""): [
            (data_root / f"{path.stem.replace('-full-data', '')}_{n}.txt").exists()
            for n in (1, 2, 3)
        ]
        for path in files
    }
    return {
        "data_root": str(data_root.resolve()),
        "seed": seed,
        "cd_hit_available": bool(shutil.which("cd-hit")),
        "official_split_files_present": split_files,
        "split_interpretation": "target-level random 80/10/10 semantics, independently deterministic with NumPy seed; train.py uses the same fallback proportions when *_1.txt..*_3.txt are absent",
        "datasets": [_audit_file(path, seed) for path in files],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=168)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.data_root, args.seed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
