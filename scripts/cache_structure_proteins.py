"""Cache frozen ESM-2 residue slots for structure or DTA mechanism records."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import torch

from contracts.mechanism import MECHANISM_RESIDUE_SLOTS
from scripts.data_contract import read_jsonl
from scripts.structure_sources.rcsb import sha256_file


ESM_MAX_RESIDUES = 1022


def _snapshot_hash(snapshot: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(value for value in snapshot.rglob("*") if value.is_file()):
        digest.update(str(path.relative_to(snapshot)).replace("\\", "/").encode("utf-8"))
        digest.update(bytes.fromhex(sha256_file(path)))
    return digest.hexdigest()


def _slot_pool(hidden: torch.Tensor, sequence_length: int) -> tuple[torch.Tensor, torch.Tensor]:
    hidden = hidden[:sequence_length].float()
    slots = torch.div(torch.arange(sequence_length, device=hidden.device) *
                      MECHANISM_RESIDUE_SLOTS, sequence_length, rounding_mode="floor")
    pooled = torch.zeros(MECHANISM_RESIDUE_SLOTS, hidden.shape[-1],
                         device=hidden.device, dtype=torch.float32)
    pooled.index_add_(0, slots, hidden)
    counts = torch.bincount(slots, minlength=MECHANISM_RESIDUE_SLOTS)
    mask = counts.gt(0)
    pooled[mask] /= counts[mask].unsqueeze(-1)
    return pooled.half().cpu(), mask.to(torch.uint8).cpu()


def _batches(items: list[tuple[str, str]], token_budget: int,
             max_batch: int) -> list[list[tuple[str, str]]]:
    ordered = sorted(items, key=lambda item: (len(item[1]), item[0]))
    batches, batch = [], []
    longest = 0
    for item in ordered:
        candidate_longest = max(longest, len(item[1]) + 2)
        if batch and (len(batch) == max_batch or
                      candidate_longest * (len(batch) + 1) > token_budget):
            batches.append(batch)
            batch, longest = [], 0
        batch.append(item)
        longest = max(longest, len(item[1]) + 2)
    if batch:
        batches.append(batch)
    return batches


def _sequence_chunks(sequence: str, size: int = ESM_MAX_RESIDUES) -> list[str]:
    if size < 1:
        raise ValueError("sequence chunk size must be positive")
    return [sequence[start:start + size] for start in range(0, len(sequence), size)]


def _record_sequence_key(record: dict) -> str:
    if "sequence_sha256" in record:
        return record["sequence_sha256"]
    if "target_key" in record:
        return record["target_key"]
    raise ValueError("protein record lacks sequence_sha256 or target_key")


def cache_structure_proteins(records_path: str | Path, output_dir: str | Path, *,
                             model_id: str, revision: str, snapshot_path: str | Path,
                             shard_size: int = 128, token_budget: int = 2048,
                             max_batch: int = 8, device: str = "cuda",
                             limit: int | None = None) -> dict:
    if shard_size < 1 or token_budget < 1 or max_batch < 1:
        raise ValueError("cache dimensions must be positive")
    records = read_jsonl(records_path)
    sequences = {_record_sequence_key(record): record["sequence"] for record in records}
    items = sorted(sequences.items())
    if limit is not None:
        items = items[:limit]
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    snapshot = Path(snapshot_path).resolve()
    if not snapshot.is_dir() or snapshot.name != revision:
        raise ValueError("snapshot_path must resolve to the declared immutable revision")

    from transformers import AutoTokenizer, EsmModel
    tokenizer = AutoTokenizer.from_pretrained(
        model_id, revision=revision, local_files_only=True)
    model = EsmModel.from_pretrained(
        model_id, revision=revision, local_files_only=True,
        torch_dtype=torch.float16, add_pooling_layer=False).to(device).eval()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    hidden_dim = int(model.config.hidden_size)
    shards = []
    with torch.inference_mode():
        for start in range(0, len(items), shard_size):
            shard_items = items[start:start + shard_size]
            filename = f"shard_{start:06d}.npz"
            path = output / filename
            if not path.is_file():
                embeddings: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
                for batch in _batches(shard_items, token_budget, max_batch):
                    keys, values = zip(*batch)
                    if len(batch) == 1 and len(values[0]) > ESM_MAX_RESIDUES:
                        pieces = []
                        for chunk in _sequence_chunks(values[0]):
                            tokens = tokenizer(chunk, return_tensors="pt",
                                               add_special_tokens=True)
                            tokens = {key: value.to(device) for key, value in tokens.items()}
                            with torch.autocast(device_type="cuda", dtype=torch.float16,
                                                enabled=device.startswith("cuda")):
                                encoded = model(**tokens).last_hidden_state[0, 1:len(chunk) + 1]
                            pieces.append(encoded)
                        residues = torch.cat(pieces, dim=0)
                        slots, mask = _slot_pool(residues, len(values[0]))
                        pooled = residues.float().mean(0).half().cpu()
                        embeddings[keys[0]] = (pooled.numpy(), slots.numpy(), mask.numpy())
                        continue
                    tokens = tokenizer(list(values), return_tensors="pt", padding=True,
                                       add_special_tokens=True)
                    tokens = {key: value.to(device) for key, value in tokens.items()}
                    with torch.autocast(device_type="cuda", dtype=torch.float16,
                                        enabled=device.startswith("cuda")):
                        hidden = model(**tokens).last_hidden_state
                    for index, (key, sequence) in enumerate(batch):
                        residues = hidden[index, 1:len(sequence) + 1]
                        slots, mask = _slot_pool(residues, len(sequence))
                        pooled = residues.float().mean(0).half().cpu()
                        embeddings[key] = (pooled.numpy(), slots.numpy(), mask.numpy())
                keys = [key for key, _ in shard_items]
                temporary = path.with_name(path.stem + ".tmp.npz")
                np.savez(temporary, keys=np.asarray(keys),
                         pooled=np.stack([embeddings[key][0] for key in keys]),
                         residues=np.stack([embeddings[key][1] for key in keys]),
                         mask=np.stack([embeddings[key][2] for key in keys]))
                temporary.replace(path)
            with np.load(path, allow_pickle=False) as shard:
                expected_keys = [key for key, _ in shard_items]
                if (shard["keys"].tolist() != expected_keys or
                        shard["residues"].shape[1:] != (MECHANISM_RESIDUE_SLOTS, hidden_dim)):
                    raise ValueError(f"cached protein shard has wrong contract: {path}")
            shards.append({"path": filename, "records": len(shard_items),
                           "sha256": sha256_file(path)})
    manifest = {
        "schema": "MetaSieve.StructureProteinBank.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "records": len(items), "hidden_dim": hidden_dim,
        "residue_slots": MECHANISM_RESIDUE_SLOTS,
        "slot_policy": "compact_resolved_sequence_index_linear_bins",
        "long_sequence_policy": "contiguous_nonoverlapping_chunks_then_global_slot_pool",
        "esm_max_residues_per_chunk": ESM_MAX_RESIDUES,
        "tensor_dtype": "float16", "mask_dtype": "uint8",
        "provider": "transformers.EsmModel", "provider_version": __import__("transformers").__version__,
        "model_id": model_id, "model_revision": revision,
        "tokenizer_revision": revision, "snapshot_path": str(snapshot),
        "snapshot_sha256": _snapshot_hash(snapshot),
        "sequence_manifest_sha256": sha256_file(records_path),
        "shards": shards,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("records")
    parser.add_argument("output")
    parser.add_argument("--model-id", default="facebook/esm2_t30_150M_UR50D")
    parser.add_argument("--revision", required=True)
    parser.add_argument("--snapshot-path", required=True)
    parser.add_argument("--shard-size", type=int, default=128)
    parser.add_argument("--token-budget", type=int, default=2048)
    parser.add_argument("--max-batch", type=int, default=8)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--limit", type=int)
    args = parser.parse_args()
    result = cache_structure_proteins(
        args.records, args.output, model_id=args.model_id, revision=args.revision,
        snapshot_path=args.snapshot_path, shard_size=args.shard_size,
        token_budget=args.token_budget, max_batch=args.max_batch,
        device=args.device, limit=args.limit)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
