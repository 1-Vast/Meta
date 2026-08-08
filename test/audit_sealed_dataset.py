"""Trusted preflight audit for a compiled source/metaval runtime seal."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import torch

from scripts.seal_compiled_dataset import SealedCompiledDataset


def _cache_keys(path: str | Path) -> set[str]:
    return set(torch.load(path, map_location="cpu", weights_only=False))


def audit_sealed_dataset(directory: str | Path, *, source_cache: str | Path | None = None,
                         metaval_cache: str | Path | None = None) -> dict:
    """Verify artifact separation before mounting either label view for a run."""
    root = Path(directory)
    source = SealedCompiledDataset(root, "source")
    metaval = SealedCompiledDataset(root, "metaval")
    source_targets = {row["target_key"] for row in source.rows}
    metaval_targets = {row["target_key"] for row in metaval.rows}
    source_rows = {row["row_id"] for row in source.rows}
    metaval_rows = {row["row_id"] for row in metaval.rows}
    errors = []
    if source_targets & metaval_targets:
        errors.append("source and metaval share target keys")
    if source_rows & metaval_rows:
        errors.append("source and metaval share label row IDs")
    if (root / "recipient").exists():
        errors.append("recipient label directory exists")
    if source_cache is not None and _cache_keys(source_cache) != source_targets:
        errors.append("source protein cache keys do not exactly match source targets")
    if metaval_cache is not None and _cache_keys(metaval_cache) != metaval_targets:
        errors.append("metaval protein cache keys do not exactly match metaval targets")
    return {
        "valid": not errors,
        "errors": errors,
        "source_targets": len(source_targets), "metaval_targets": len(metaval_targets),
        "source_rows": len(source_rows), "metaval_rows": len(metaval_rows),
        "recipient_label_artifact_emitted": False,
        "source_audit": source.audit_snapshot(), "metaval_audit": metaval.audit_snapshot(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("sealed_dir")
    parser.add_argument("--source-cache")
    parser.add_argument("--metaval-cache")
    parser.add_argument("--out")
    args = parser.parse_args()
    result = audit_sealed_dataset(
        args.sealed_dir, source_cache=args.source_cache, metaval_cache=args.metaval_cache,
    )
    if args.out:
        Path(args.out).write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
