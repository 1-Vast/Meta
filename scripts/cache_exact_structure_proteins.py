"""Cache exact frozen ESM/P1B residue states in an immutable ragged bank.

Unlike the retained 128-slot bank, this research-only cache preserves one row
per sequence residue.  It does not fit ESM or P1B and stores only the frozen
P1B-projected states needed by R0.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import torch

from scripts.cache_structure_proteins import (
    ESM_MAX_RESIDUES,
    _batches,
    _sequence_chunks,
    _slot_pool,
    _snapshot_hash,
)
from scripts.data_contract import read_jsonl
from scripts.structure_sources.rcsb import sha256_file


def sequence_key(sequence: str) -> str:
    if not sequence:
        raise ValueError("protein sequence cannot be empty")
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def pack_projected_rows(rows: list[tuple[str, np.ndarray, np.ndarray, np.ndarray]]) -> dict:
    """Pack variable-length exact rows and fixed slot rows into one NPZ payload."""
    if not rows:
        raise ValueError("cannot pack an empty exact-protein shard")
    keys, exact, slots, masks = zip(*rows)
    if len(set(keys)) != len(keys):
        raise ValueError("exact-protein shard contains duplicate keys")
    hidden = exact[0].shape[1]
    if any(value.ndim != 2 or value.shape[1] != hidden for value in exact):
        raise ValueError("exact projected states have inconsistent dimensions")
    if any(value.shape != (128, hidden) for value in slots):
        raise ValueError("projected slot states violate the frozen 128-slot contract")
    if any(value.shape != (128,) for value in masks):
        raise ValueError("slot masks violate the frozen 128-slot contract")
    lengths = np.asarray([len(value) for value in exact], dtype=np.int32)
    offsets = np.concatenate(([0], np.cumsum(lengths, dtype=np.int64)))
    return {
        "keys": np.asarray(keys),
        "lengths": lengths,
        "offsets": offsets,
        "exact_projected": np.concatenate(exact, axis=0).astype(np.float16),
        "slot_projected": np.stack(slots).astype(np.float16),
        "slot_mask": np.stack(masks).astype(np.uint8),
    }


def unpack_projected_row(payload: dict[str, np.ndarray], index: int) -> dict[str, np.ndarray]:
    """Read one copied row from an in-memory ragged shard payload."""
    if index < 0 or index >= len(payload["keys"]):
        raise IndexError("exact-protein row index is outside the shard")
    left, right = int(payload["offsets"][index]), int(payload["offsets"][index + 1])
    if right - left != int(payload["lengths"][index]):
        raise ValueError("exact-protein offsets and lengths disagree")
    return {
        "exact_projected": payload["exact_projected"][left:right].copy(),
        "slot_projected": payload["slot_projected"][index].copy(),
        "slot_mask": payload["slot_mask"][index].copy(),
    }


def _encode_batch(model, tokenizer, sequences: list[str], device: str) -> list[torch.Tensor]:
    tokens = tokenizer(sequences, return_tensors="pt", padding=True,
                       add_special_tokens=True)
    tokens = {key: value.to(device) for key, value in tokens.items()}
    with torch.autocast(device_type="cuda", dtype=torch.float16):
        hidden = model(**tokens).last_hidden_state
    return [hidden[index, 1:len(sequence) + 1]
            for index, sequence in enumerate(sequences)]


def _encode_sequence(model, tokenizer, sequence: str, device: str) -> torch.Tensor:
    if len(sequence) <= ESM_MAX_RESIDUES:
        return _encode_batch(model, tokenizer, [sequence], device)[0]
    pieces = []
    for chunk in _sequence_chunks(sequence):
        pieces.append(_encode_batch(model, tokenizer, [chunk], device)[0])
    return torch.cat(pieces, dim=0)


def cache_exact_structure_proteins(records_path: str | Path,
                                   checkpoint_path: str | Path,
                                   output_dir: str | Path, *,
                                   model_id: str, revision: str,
                                   snapshot_path: str | Path,
                                   shard_size: int = 64,
                                   token_budget: int = 2048,
                                   max_batch: int = 4,
                                   device: str = "cuda",
                                   limit: int | None = None) -> dict:
    if not device.startswith("cuda") or not torch.cuda.is_available():
        raise RuntimeError("the exact frozen protein cache is registered for CUDA")
    if min(shard_size, token_budget, max_batch) < 1:
        raise ValueError("cache dimensions must be positive")
    records_path = Path(records_path)
    checkpoint_path = Path(checkpoint_path)
    snapshot = Path(snapshot_path).resolve()
    if not snapshot.is_dir() or snapshot.name != revision:
        raise ValueError("snapshot_path must resolve to the declared immutable revision")
    records = read_jsonl(records_path)
    sequences = {sequence_key(record["sequence"]): record["sequence"] for record in records}
    items = sorted(sequences.items())
    if limit is not None:
        items = items[:limit]
    if not items:
        raise ValueError("no unique sequences selected for exact caching")

    from transformers import AutoTokenizer, EsmModel
    from research.e0_identifiability.run_tdir_pilot import _load_frozen_model

    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    protein_dim = int(checkpoint["protein_dim"])
    frozen, _ = _load_frozen_model(checkpoint_path, protein_dim, device)
    tokenizer = AutoTokenizer.from_pretrained(
        model_id, revision=revision, local_files_only=True)
    esm = EsmModel.from_pretrained(
        model_id, revision=revision, local_files_only=True,
        torch_dtype=torch.float16, add_pooling_layer=False).to(device).eval()
    if int(esm.config.hidden_size) != protein_dim:
        raise ValueError("ESM hidden dimension differs from the frozen P1B checkpoint")
    for parameter in esm.parameters():
        parameter.requires_grad_(False)

    output = Path(output_dir)
    if (output / "manifest.json").exists():
        raise FileExistsError(f"exact protein bank is already complete: {output}")
    output.mkdir(parents=True, exist_ok=True)
    shards = []
    with torch.inference_mode():
        for start in range(0, len(items), shard_size):
            shard_items = items[start:start + shard_size]
            filename = f"shard_{start:06d}.npz"
            path = output / filename
            if not path.exists():
                encoded: dict[str, torch.Tensor] = {}
                ordinary = [(key, sequence) for key, sequence in shard_items
                            if len(sequence) <= ESM_MAX_RESIDUES]
                for batch in _batches(ordinary, token_budget, max_batch):
                    keys, values = zip(*batch)
                    hidden_rows = _encode_batch(esm, tokenizer, list(values), device)
                    encoded.update(zip(keys, hidden_rows))
                for key, sequence in shard_items:
                    if len(sequence) > ESM_MAX_RESIDUES:
                        encoded[key] = _encode_sequence(esm, tokenizer, sequence, device)

                packed_rows = []
                for key, sequence in shard_items:
                    hidden = encoded[key]
                    if hidden.shape != (len(sequence), protein_dim):
                        raise ValueError("ESM exact residue output violates sequence contract")
                    slots_raw, slot_mask = _slot_pool(hidden, len(sequence))
                    exact_projected = frozen.protein.bank_proj(hidden.float()).half().cpu().numpy()
                    slot_projected = frozen.protein.bank_proj(
                        slots_raw.to(device).float()).half().cpu().numpy()
                    packed_rows.append((key, exact_projected,
                                        slot_projected, slot_mask.numpy()))
                payload = pack_projected_rows(packed_rows)
                temporary = path.with_name(path.stem + ".tmp.npz")
                np.savez_compressed(temporary, **payload)
                temporary.replace(path)
            with np.load(path, allow_pickle=False) as shard:
                expected = [key for key, _ in shard_items]
                if shard["keys"].tolist() != expected:
                    raise ValueError(f"cached exact-protein shard has wrong keys: {path}")
                if shard["slot_projected"].shape[1:] != (128, frozen.protein.bank_proj.out_features):
                    raise ValueError(f"cached exact-protein shard has wrong shape: {path}")
                if int(shard["lengths"].sum()) != len(shard["exact_projected"]):
                    raise ValueError(f"cached exact-protein shard has invalid ragged offsets: {path}")
            shards.append({"path": filename, "records": len(shard_items),
                           "sha256": sha256_file(path)})

    manifest = {
        "schema": "MetaSieve.ExactStructureProteinBank.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "records": len(items),
        "raw_hidden_dim": protein_dim,
        "projected_dim": int(frozen.protein.bank_proj.out_features),
        "residue_slots": 128,
        "exact_policy": "one frozen projected row per input sequence residue",
        "slot_policy": "same compact linear 128-bin pooling as frozen P1B",
        "long_sequence_policy": "contiguous nonoverlapping ESM chunks then exact concatenate",
        "tensor_dtype": "float16",
        "model_id": model_id,
        "model_revision": revision,
        "snapshot_path": str(snapshot),
        "snapshot_sha256": _snapshot_hash(snapshot),
        "checkpoint_path": str(checkpoint_path.resolve()),
        "checkpoint_sha256": sha256_file(checkpoint_path),
        "records_sha256": sha256_file(records_path),
        "unique_sequence_set_sha256": hashlib.sha256(
            "\n".join(key for key, _ in items).encode("ascii")).hexdigest(),
        "affinity_labels_used": False,
        "trainable_parameters": 0,
        "shards": shards,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("records")
    parser.add_argument("checkpoint")
    parser.add_argument("output")
    parser.add_argument("--model-id", default="facebook/esm2_t30_150M_UR50D")
    parser.add_argument("--revision", required=True)
    parser.add_argument("--snapshot-path", required=True)
    parser.add_argument("--shard-size", type=int, default=64)
    parser.add_argument("--token-budget", type=int, default=2048)
    parser.add_argument("--max-batch", type=int, default=4)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    result = cache_exact_structure_proteins(
        args.records, args.checkpoint, args.output,
        model_id=args.model_id, revision=args.revision,
        snapshot_path=args.snapshot_path, shard_size=args.shard_size,
        token_budget=args.token_budget, max_batch=args.max_batch,
        device=args.device, limit=args.limit)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
