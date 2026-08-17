"""Build a StructureProteinBank.v1 sharded bank from the local ESM-650M cache.

Input: tools/runtime/esm2_t33_650M_pooled/embeddings.npz (keys = sequence
sha256, pooled [n,1280] fp16, residues [n,128,1280] fp16, mask [n,128] uint8).
Output: a bank directory loadable by QPSMPData with protein_dim = 1280.
External representation lane; the source cache is recorded with hashes and
the bank is a research cache, not a numerical authority.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "tools/runtime/esm2_t33_650M_pooled/embeddings.npz"
OUT = ROOT / "tools/runtime/esm2_t33_650M_protein_bank"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    if OUT.exists():
        print(f"{OUT} already exists");
        return 0
    OUT.mkdir(parents=True)
    with np.load(SOURCE, allow_pickle=False) as store:
        keys = [str(k) for k in store["keys"]]
        pooled = store["pooled"].astype(np.float16)
        residues = store["residues"].astype(np.float16)
        mask = store["mask"].astype(np.uint8)
    shard_path = OUT / "shard_000000.npz"
    np.savez(shard_path, keys=np.asarray(keys), pooled=pooled,
             residues=residues, mask=mask)
    manifest = {
        "schema": "MetaSieve.StructureProteinBank.v1",
        "created_utc": "2026-08-17T00:00:00+00:00",
        "records": len(keys),
        "hidden_dim": int(residues.shape[-1]),
        "residue_slots": int(residues.shape[1]),
        "slot_policy": "compact_resolved_sequence_index_linear_bins",
        "long_sequence_policy": "contiguous_nonoverlapping_chunks_then_global_slot_pool",
        "esm_max_residues_per_chunk": 1022,
        "tensor_dtype": "float16",
        "mask_dtype": "uint8",
        "provider": "transformers.EsmModel",
        "provider_version": "local",
        "model_id": "local:esm2_t33_650M_UR50D",
        "model_revision": "local-snapshot",
        "tokenizer_revision": "local-snapshot",
        "source_cache": str(SOURCE.resolve()),
        "source_cache_sha256": sha256_file(SOURCE),
        "shards": [{
            "path": "shard_000000.npz",
            "records": len(keys),
            "sha256": sha256_file(shard_path),
        }],
    }
    (OUT / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(manifest, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
