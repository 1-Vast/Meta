"""Create physically isolated runtime label views from canonical DTA rows.

The compiler is a trusted preparation step and may contain train/val/test
labels in one staging file. Model processes must never receive that file. This
module emits independent source and metaval files and no recipient label file.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from scripts.data_contract import read_jsonl, write_jsonl
from contracts.biological_context import DEFAULT_ENDPOINT_VOCABULARY, UNKNOWN_ENDPOINT


SEAL_SCHEMA = "MetaSieve.CompiledDTASeal.v2"
PERMITTED_SPLITS = frozenset({"source", "metaval"})
SPLIT_MAP = {"train": "source", "val": "metaval", "test": "recipient"}
REQUIRED_COLUMNS = frozenset({
    "row_id", "drug_key", "target_key", "pair_key", "task_key", "smiles",
    "sequence", "endpoint_key", "context_key", "context_id", "context_cont",
    "context_mask", "split", "y", "provenance",
})


def _upgrade_legacy_context(row: dict) -> dict:
    """Upgrade trusted v1 staging rows; emitted runtime rows are v2-only."""
    if "context_id" in row:
        return row
    if "gamma" not in row and "gamma_key" not in row:
        return row
    value = dict(row)
    lookup = {key.casefold(): key for key in DEFAULT_ENDPOINT_VOCABULARY}
    endpoint = lookup.get(str(value.pop("gamma_key", "")).casefold(), UNKNOWN_ENDPOINT)
    value.pop("gamma", None)
    value.update({
        "endpoint_key": endpoint,
        "context_key": hashlib.sha256(endpoint.encode("utf-8")).hexdigest(),
        "context_id": DEFAULT_ENDPOINT_VOCABULARY.index(endpoint),
        "context_cont": [],
        "context_mask": [],
    })
    return value


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _validate_rows(rows: list[dict]) -> None:
    if not rows:
        raise ValueError("canonical input contains no rows")
    row_ids, measurement_keys, target_splits = set(), set(), {}
    for row in rows:
        missing = REQUIRED_COLUMNS - set(row)
        if missing:
            raise ValueError(f"canonical row is missing fields: {sorted(missing)}")
        if row["split"] not in SPLIT_MAP:
            raise ValueError(f"unexpected canonical split {row['split']!r}")
        if not 0.0 <= float(row["y"]) <= 1.0:
            raise ValueError(f"row {row['row_id']!r} has y outside [0,1]")
        if row["row_id"] in row_ids:
            raise ValueError(f"duplicate row_id {row['row_id']!r}")
        measurement_key = (row["task_key"], row["pair_key"])
        if measurement_key in measurement_keys:
            raise ValueError(f"unresolved duplicate measurement {measurement_key!r}")
        row_ids.add(row["row_id"])
        measurement_keys.add(measurement_key)
        target_splits.setdefault(row["target_key"], set()).add(row["split"])
    crossed = [target for target, splits in target_splits.items() if len(splits) != 1]
    if crossed:
        raise ValueError(f"{len(crossed)} target(s) cross canonical split boundaries")


def seal_compiled_dataset(compiled_dir: str | Path, output_dir: str | Path) -> dict:
    """Seal a canonical staging directory into single-label runtime views."""
    compiled = Path(compiled_dir).resolve()
    output = Path(output_dir).resolve()
    rows_path = compiled / "rows.jsonl"
    if not rows_path.is_file():
        raise FileNotFoundError(rows_path)
    if output.exists():
        raise FileExistsError(f"sealed output already exists: {output}")
    rows = [_upgrade_legacy_context(row) for row in read_jsonl(rows_path)]
    _validate_rows(rows)

    source_rows, metaval_rows = [], []
    governance: dict[str, dict] = {}
    for row in rows:
        mapped = SPLIT_MAP[row["split"]]
        governance.setdefault(row["target_key"], {
            "target_key": row["target_key"], "split": mapped, "rows": 0,
        })["rows"] += 1
        if mapped in PERMITTED_SPLITS:
            isolated = dict(row)
            isolated["split"] = mapped
            (source_rows if mapped == "source" else metaval_rows).append(isolated)

    output.mkdir(parents=True)
    source_path = output / "source" / "rows.jsonl"
    metaval_path = output / "metaval" / "rows.jsonl"
    governance_path = output / "governance.jsonl"
    write_jsonl(source_path, source_rows)
    write_jsonl(metaval_path, metaval_rows)
    write_jsonl(governance_path, [governance[key] for key in sorted(governance)])

    upstream_manifest = compiled / "manifest.json"
    upstream = json.loads(upstream_manifest.read_text(encoding="utf-8")) \
        if upstream_manifest.is_file() else {}
    artifacts = {
        "source/rows.jsonl": sha256_file(source_path),
        "metaval/rows.jsonl": sha256_file(metaval_path),
        "governance.jsonl": sha256_file(governance_path),
    }
    counts = {
        split: {
            "rows": sum(1 for row in rows if SPLIT_MAP[row["split"]] == split),
            "targets": sum(1 for row in governance.values() if row["split"] == split),
        }
        for split in ("source", "metaval", "recipient")
    }
    manifest = {
        "schema": SEAL_SCHEMA,
        "dataset": upstream.get("dataset"),
        "compiled_rows_sha256": sha256_file(rows_path),
        "compiled_manifest_sha256": sha256_file(upstream_manifest)
        if upstream_manifest.is_file() else None,
        "normalization": upstream.get("normalization"),
        "biological_context": upstream.get("biological_context"),
        "split_mapping": SPLIT_MAP,
        "counts": counts,
        "artifacts": artifacts,
        "recipient_label_artifact_emitted": False,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


class SealedCompiledDataset:
    """Load exactly one permitted label file and verify its seal bindings."""

    def __init__(self, directory: str | Path, split: str, expected_dataset: str | None = None):
        if split not in PERMITTED_SPLITS:
            raise ValueError(f"runtime split must be one of {sorted(PERMITTED_SPLITS)}")
        self.directory = Path(directory).resolve()
        manifest_path = self.directory / "manifest.json"
        self.manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.manifest_sha256 = sha256_file(manifest_path)
        if self.manifest.get("schema") != SEAL_SCHEMA:
            raise ValueError("unsupported compiled dataset seal")
        if expected_dataset and self.manifest.get("dataset") not in {None, expected_dataset}:
            raise ValueError("sealed dataset name does not match the requested dataset")
        if self.manifest.get("recipient_label_artifact_emitted") is not False:
            raise ValueError("seal does not guarantee recipient-label absence")
        self.label_split = split
        self.label_path = self.directory / split / "rows.jsonl"
        expected = self.manifest["artifacts"].get(f"{split}/rows.jsonl")
        if expected != sha256_file(self.label_path):
            raise ValueError(f"sealed {split} label artifact hash mismatch")
        self.rows = read_jsonl(self.label_path)
        if not self.rows:
            raise ValueError(f"sealed {split} view contains no labels")
        if any(row.get("split") != split for row in self.rows):
            raise ValueError(f"sealed {split} view contains another split")

    def audit_snapshot(self) -> dict:
        return {
            "visible_label_split": self.label_split,
            "label_rows_loaded": len(self.rows),
            "recipient_label_reads": 0,
            "recipient_label_artifact_emitted": False,
            "mounted_label_artifact": str(self.label_path),
        }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("compiled_dir")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    print(json.dumps(seal_compiled_dataset(args.compiled_dir, args.output_dir), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
