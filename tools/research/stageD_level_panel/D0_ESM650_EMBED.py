"""Compute pooled + 128-slot ESM-2 650M embeddings for the governed protein set.

External protein representation lane (reported as external data). The local
snapshot dataset/raw/esm2_t33_650M_UR50D is used; embeddings are mean-pooled
over residues (pooled) and slot-pooled into 128 linear bins (residues), matching
the MECHANISM_RESIDUE_SLOTS convention of the active 150M bank.

Output is a research cache under tools/runtime (never a numerical authority);
the manifest records model identity and hashes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import sys

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

ESM_MAX_RESIDUES = 1022
SLOTS = 128


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def slot_pool(hidden: torch.Tensor, sequence_length: int):
    hidden = hidden[:sequence_length].float()
    slots = torch.div(torch.arange(sequence_length, device=hidden.device) * SLOTS,
                      sequence_length, rounding_mode="floor")
    pooled = torch.zeros(SLOTS, hidden.shape[-1], device=hidden.device,
                         dtype=torch.float32)
    pooled.index_add_(0, slots, hidden)
    counts = torch.bincount(slots, minlength=SLOTS)
    mask = counts.gt(0)
    pooled[mask] /= counts[mask].unsqueeze(-1)
    return pooled.half().cpu(), mask.to(torch.uint8).cpu()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model-dir", type=Path,
                        default=ROOT / "dataset/raw/esm2_t33_650M_UR50D")
    parser.add_argument("--proteins", type=Path,
                        default=ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0/proteins.jsonl")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "tools/runtime/esm2_t33_650M_pooled")
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--limit", type=int, default=None)
    args = parser.parse_args()

    rows = []
    with args.proteins.open(encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    if args.limit is not None:
        rows = rows[:args.limit]
    items = sorted((r["sequence_sha256"], r["sequence"]) for r in rows)
    print(f"proteins: {len(items)}", flush=True)

    from transformers import AutoTokenizer, EsmModel
    tokenizer = AutoTokenizer.from_pretrained(
        str(args.model_dir), local_files_only=True)
    model = EsmModel.from_pretrained(
        str(args.model_dir), local_files_only=True,
        torch_dtype=torch.float16, add_pooling_layer=False).to(args.device).eval()
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    print(f"model loaded, hidden {model.config.hidden_size}", flush=True)

    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    pooled_list, residue_list, mask_list, keys = [], [], [], []
    with torch.inference_mode():
        for key, sequence in items:
            pieces = []
            for start in range(0, len(sequence), ESM_MAX_RESIDUES):
                chunk = sequence[start:start + ESM_MAX_RESIDUES]
                tokens = tokenizer(chunk, return_tensors="pt",
                                   add_special_tokens=True)
                tokens = {k: v.to(args.device) for k, v in tokens.items()}
                with torch.autocast(device_type="cuda", dtype=torch.float16,
                                    enabled=args.device.startswith("cuda")):
                    encoded = model(**tokens).last_hidden_state[0, 1:len(chunk) + 1]
                pieces.append(encoded.float())
            residues = torch.cat(pieces, dim=0)
            slots, mask = slot_pool(residues, len(sequence))
            pooled = residues.mean(0).half().cpu()
            keys.append(key)
            pooled_list.append(pooled.numpy())
            residue_list.append(slots.numpy())
            mask_list.append(mask.numpy())
            if len(keys) % 25 == 0:
                print(f"{len(keys)}/{len(items)}", flush=True)

    np.savez_compressed(out / "embeddings.npz",
                        keys=np.asarray(keys),
                        pooled=np.stack(pooled_list),
                        residues=np.stack(residue_list),
                        mask=np.stack(mask_list))
    manifest = {
        "schema": "MetaSieve.StageD.ESM650Pooled.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "records": len(keys),
        "hidden_dim": int(model.config.hidden_size),
        "layers": int(model.config.num_hidden_layers),
        "residue_slots": SLOTS,
        "pooling": "residue_mean",
        "long_sequence_policy": "contiguous_nonoverlapping_chunks_then_global_slot_pool",
        "esm_max_residues_per_chunk": ESM_MAX_RESIDUES,
        "model_dir": str(args.model_dir.resolve()),
        "model_config_sha256": sha256_file(args.model_dir / "config.json"),
        "model_weights_sha256": sha256_file(args.model_dir / "pytorch_model.bin"),
        "proteins_sha256": sha256_file(args.proteins),
        "embeddings_sha256": sha256_file(out / "embeddings.npz"),
        "tensor_dtype": "float16",
    }
    (out / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(manifest, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
