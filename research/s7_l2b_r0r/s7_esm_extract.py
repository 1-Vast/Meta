"""S7_L2B — frozen ESM2-650M per-residue embedding extraction.

Windowing contract (frozen in PREREG_S7_L2B_UNIFIED.md lineage):
  window 1000 residues, stride 500, an extra final window ending exactly at the
  C-terminus when needed, BOS/EOS excluded from residue output, overlapping
  positions averaged arithmetically, canonical residue indices preserved,
  no silent truncation.

Outputs a single fp16 memmap plus an index, and a runtime manifest with model
revision, weight SHA-256 and environment versions.
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_dataset import build, make_split, protein_components  # noqa: E402

OUT = ROOT / "report" / "s7_l2b_r0r"
CACHE = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "esm2_650M"
# Load from the locally acquired, SHA-256-verified copy so extraction runs
# entirely offline. Acquisition manifest: P0_ESM2_ACQUISITION_MANIFEST.json,
# weights sha256 c874668852c7275a159e2c7ceb6069671d7b1ba2c7b52f59600b34ce0f721008,
# revision 08e4846e537177426273712802403f7ba8261b6c.
MODEL_ID = str(ROOT / "dataset" / "raw" / "esm2_t33_650M_UR50D")
MODEL_NAME = "facebook/esm2_t33_650M_UR50D"
WEIGHTS_SHA256 = "c874668852c7275a159e2c7ceb6069671d7b1ba2c7b52f59600b34ce0f721008"
WINDOW, STRIDE, DIM = 1000, 500, 1280


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    from transformers import AutoTokenizer, EsmModel
    t0 = time.time()
    tok = AutoTokenizer.from_pretrained(MODEL_ID)
    model = EsmModel.from_pretrained(MODEL_ID, torch_dtype=torch.float16)
    model.eval()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model.to(dev)
    print(f"loaded {MODEL_ID} in {time.time()-t0:.0f}s on {dev}", flush=True)

    kept, _q, _c, _f = build()
    comp = protein_components(kept)
    train, _held_all, hA, hB = make_split(kept, comp)
    need = {}
    for r in train + hA + hB:
        need.setdefault(r["seq_key"], r["uniprot_sequence"])
    keys = sorted(need)
    lens = [len(need[k]) for k in keys]
    total = int(sum(lens))
    print(f"unique sequences: {len(keys)}  total residues: {total}  "
          f"max len: {max(lens)}", flush=True)

    memmap_path = CACHE / "esm2_650M_residues.fp16.dat"
    arr = np.memmap(memmap_path, dtype=np.float16, mode="w+", shape=(total, DIM))
    index, off = {}, 0
    checks = {"length_exact": 0, "length_mismatch": 0, "nonfinite": 0}

    with torch.no_grad():
        for n, k in enumerate(keys):
            seq = need[k]
            L = len(seq)
            acc = np.zeros((L, DIM), dtype=np.float32)
            cnt = np.zeros((L, 1), dtype=np.float32)
            starts = list(range(0, max(L - WINDOW, 0) + 1, STRIDE)) or [0]
            if starts[-1] + WINDOW < L:
                starts.append(L - WINDOW)          # final C-terminal window
            for s in starts:
                e = min(s + WINDOW, L)
                sub = seq[s:e]
                enc = tok(sub, return_tensors="pt", add_special_tokens=True)
                enc = {kk: v.to(dev) for kk, v in enc.items()}
                out = model(**enc).last_hidden_state[0]
                # strip BOS and EOS
                res = out[1:1 + (e - s)].float().cpu().numpy()
                acc[s:e] += res
                cnt[s:e] += 1.0
            if (cnt == 0).any():
                raise RuntimeError(f"uncovered residues for {k}")
            emb = acc / cnt
            if emb.shape[0] == L:
                checks["length_exact"] += 1
            else:
                checks["length_mismatch"] += 1
            if not np.isfinite(emb).all():
                checks["nonfinite"] += 1
            arr[off:off + L] = emb.astype(np.float16)
            index[k] = [off, L]
            off += L
            if (n + 1) % 200 == 0:
                print(f"  {n+1}/{len(keys)} sequences, {off}/{total} residues, "
                      f"{time.time()-t0:.0f}s", flush=True)
    arr.flush()
    del arr

    (CACHE / "esm2_650M_index.json").write_text(json.dumps(index), encoding="utf-8")
    h = hashlib.sha256()
    with memmap_path.open("rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""):
            h.update(c)
    manifest = {
        "schema": "MetaSieve.S7L2B.PLMRuntimeManifest.v1",
        "model_id": MODEL_NAME,
        "loaded_from_local_dir": MODEL_ID,
        "weights_sha256": WEIGHTS_SHA256,
        "model_revision": "08e4846e537177426273712802403f7ba8261b6c",
        "offline": True,
        "hidden_dim": DIM,
        "frozen": True,
        "eval_mode": True,
        "dtype": "float16",
        "window": WINDOW, "stride": STRIDE,
        "bos_eos_excluded": True,
        "overlap_policy": "arithmetic mean",
        "sequences": len(keys), "total_residues": total,
        "max_sequence_length": int(max(lens)),
        "cache_path": str(memmap_path.relative_to(ROOT)).replace("\\", "/"),
        "cache_sha256": h.hexdigest(),
        "cache_bytes": memmap_path.stat().st_size,
        "length_checks": checks,
        "torch": torch.__version__,
        "cuda": torch.version.cuda,
        "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "python": sys.version.split()[0],
        "elapsed_sec": round(time.time() - t0, 1),
    }
    (OUT / "GPU_PREFLIGHT_AUDIT.json").write_text(json.dumps(manifest, indent=2),
                                                  encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
