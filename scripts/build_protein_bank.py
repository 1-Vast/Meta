"""Materialize protein embeddings through an explicit caller-provided backend."""
from __future__ import annotations

import argparse
import csv
import hashlib
from importlib import import_module
import json
from pathlib import Path
from typing import Callable

import torch

from scripts.data_contract import read_jsonl


def residue_slot_mapping(sequence_length: int, slot_count: int) -> dict[str, torch.Tensor | int]:
    """Map each reduced slot to a deterministic half-open sequence interval."""
    if sequence_length < 1 or slot_count < 1:
        raise ValueError("sequence_length and slot_count must be positive")
    edges = torch.div(torch.arange(slot_count + 1) * sequence_length,
                      slot_count, rounding_mode="floor")
    starts, ends = edges[:-1], edges[1:]
    return {"residue_slot_start": starts, "residue_slot_end": ends,
            "residue_mask": ends.gt(starts), "sequence_length": sequence_length}


def _bank_entry(value, sequence: str) -> dict:
    if isinstance(value, dict):
        pooled, residues = value["pooled"], value["residues"]
    elif isinstance(value, tuple) and len(value) == 2:
        pooled, residues = value
    else:
        raise TypeError("protein provider/cache entry must contain pooled and residues tensors")
    if pooled.ndim != 1 or residues.ndim != 2:
        raise ValueError("protein provider/cache tensors have invalid ranks")
    return {"pooled": pooled, "residues": residues,
            **residue_slot_mapping(len(sequence), int(residues.shape[0]))}


def resolve_provider(reference: str) -> Callable[[dict[str, str]], dict[str, tuple[torch.Tensor, torch.Tensor]]]:
    """Resolve ``module:function`` without making model downloads a hidden side effect."""
    module_name, separator, name = reference.partition(":")
    if not separator:
        raise ValueError("embedding provider must use module:function syntax")
    provider = getattr(import_module(module_name), name)
    if not callable(provider):
        raise TypeError("embedding provider is not callable")
    return provider


PROVIDER_PROVENANCE_FIELDS = frozenset({
    "provider", "model_id", "model_revision", "tokenizer_revision",
    "pooling", "slot_policy", "dtype",
})


def _tensor_hash(value: torch.Tensor) -> str:
    tensor = value.detach().cpu().contiguous()
    digest = hashlib.sha256()
    digest.update(str(tensor.dtype).encode("ascii"))
    digest.update(json.dumps(list(tensor.shape)).encode("ascii"))
    digest.update(tensor.view(torch.uint8).numpy().tobytes())
    return digest.hexdigest()


def build_protein_bank(rows_path: str | Path, output_path: str | Path, provider, *,
                       provider_metadata: dict) -> dict:
    missing_metadata = PROVIDER_PROVENANCE_FIELDS - set(provider_metadata)
    if missing_metadata:
        raise ValueError(f"protein provider metadata missing fields: {sorted(missing_metadata)}")
    if any(not str(provider_metadata[key]).strip() for key in PROVIDER_PROVENANCE_FIELDS):
        raise ValueError("protein provider metadata values must be nonempty")
    rows = read_jsonl(rows_path)
    sequences = {row["target_key"]: row["sequence"] for row in rows}
    keys = sorted(sequences)
    values = provider({key: sequences[key] for key in keys})
    if set(values) != set(keys):
        raise ValueError("embedding provider must return every target_key")
    values = {key: _bank_entry(value, sequences[key]) for key, value in values.items()}
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    torch.save(values, target)
    first = values[keys[0]]
    tensor_hashes = {key: {"pooled": _tensor_hash(values[key]["pooled"]),
                           "residues": _tensor_hash(values[key]["residues"])}
                     for key in keys}
    result = {"schema": "MetaSieve.ProteinBank.v3", "targets": len(keys),
              "pooled_dim": int(first["pooled"].numel()),
              "residue_shape": list(first["residues"].shape),
              "provider": dict(provider_metadata),
              "rows_sha256": _hash_file(Path(rows_path)),
              "tensor_hashes": tensor_hashes,
              "output": str(target), "sha256": _hash_file(target)}
    target.with_suffix(".manifest.json").write_text(
        json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result


def _hash_file(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _legacy_sequence_index(path: Path, target_column: str, sequence_column: str,
                           delimiter: str) -> dict[str, str]:
    index = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter=delimiter):
            target, sequence = str(row[target_column]), str(row[sequence_column])
            key = hashlib.sha256("".join(sequence.upper().split()).encode("utf-8")).hexdigest()
            previous = index.setdefault(key, target)
            if previous != target:
                raise ValueError("one canonical sequence maps to multiple legacy target IDs")
    return index


def migrate_protein_cache(rows_path: str | Path, legacy_cache_path: str | Path,
                          legacy_index_path: str | Path, output_path: str | Path, *,
                          target_column: str, sequence_column: str,
                          delimiter: str = "\t") -> dict:
    """Migrate a legacy cache into exactly one sealed label view's target set."""
    rows_path, legacy_cache_path = Path(rows_path), Path(legacy_cache_path)
    legacy_index_path, output_path = Path(legacy_index_path), Path(output_path)
    rows = read_jsonl(rows_path)
    sequences = {str(row["target_key"]): str(row["sequence"]) for row in rows}
    expected = sorted(sequences)
    index = _legacy_sequence_index(legacy_index_path, target_column, sequence_column, delimiter)
    cache = torch.load(legacy_cache_path, map_location="cpu", weights_only=False)
    missing = [key for key in expected if key not in index or index[key] not in cache]
    if missing:
        raise ValueError(f"legacy cache does not cover {len(missing)} sealed targets")
    migrated = {key: _bank_entry(cache[index[key]], sequences[key]) for key in expected}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(migrated, output_path)
    result = {
        "schema": "MetaSieve.ViewProteinCache.v2", "targets": len(migrated),
        "view_rows_sha256": _hash_file(rows_path),
        "legacy_cache_sha256": _hash_file(legacy_cache_path),
        "output": str(output_path), "sha256": _hash_file(output_path),
    }
    output_path.with_suffix(".json").write_text(json.dumps(result, indent=2, sort_keys=True),
                                                  encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("rows")
    parser.add_argument("output")
    parser.add_argument("--provider", required=True, help="module:function")
    parser.add_argument("--provider-metadata", required=True,
                        help="JSON file with immutable provider/model provenance")
    args = parser.parse_args()
    metadata = json.loads(Path(args.provider_metadata).read_text(encoding="utf-8"))
    result = build_protein_bank(args.rows, args.output, resolve_provider(args.provider),
                                provider_metadata=metadata)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
