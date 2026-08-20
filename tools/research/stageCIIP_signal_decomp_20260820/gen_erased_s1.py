"""CIIP-S1 S0 item 9: mutation-erasure ESM cache (plan 9.1 item 9; prereg B.1).

For each of the 49 covered pairs: replace the verified mutation position with
'X' in BOTH the WT and variant sequences, ASSERT the erased strings are
exactly equal, run local ESM-2-150M on CPU, store residue states. Because the
erased strings are identical, the WT-erased and variant-erased embeddings are
bit-identical; the cache stores one array per pair (ewt) plus a byte-hash of
the erased string; the <=1e-5 assert is evaluated as max|we - me| == 0.
Writes ERASED_ESM_S1.npz + ERASED_S1_SHA256.txt in THIS directory only.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
BRIDGE = HERE.parent / "stageCIIP_potential_bridge"
SIG = HERE.parent / "stageX_csc_signal"
X0C = SIG / "stageX0c_measurement_qualification_20260818"
sys.path.insert(0, str(SIG))
sys.path.insert(0, str(X0C))

from x0_common import load_duongly, normalize_construct_name  # noqa: E402
from x0_i2 import load_esm, build_pair_records, ESM_MAX_LEN  # noqa: E402

torch.set_num_threads(8)
d1 = json.loads((BRIDGE / "DATA1A.json").read_text(encoding="utf-8"))
d2 = json.loads((BRIDGE / "DATA2X2.json").read_text(encoding="utf-8"))
pair_table = json.loads((SIG / "X0_PAIR_TABLE.json").read_text(encoding="utf-8"))
_, _, seqs = load_duongly()
records = build_pair_records(pair_table, seqs)
by_construct = {normalize_construct_name(r["construct"]): r for r in records}
cov = list(d2["covered_pair_indices"])

tok, model, device = load_esm("cpu")  # CPU frozen by prereg B.1
print("device:", device, flush=True)


def erase_at(seq: str, pos: int) -> str:
    return seq[: pos - 1] + "X" + seq[pos:]


def hidden(seq: str) -> np.ndarray:
    enc = tok(seq, return_tensors="pt", truncation=True, max_length=1022)
    enc = {k: v for k, v in enc.items()}
    with torch.no_grad():
        out = model(**enc)
    return out.last_hidden_state[0].cpu().float().numpy()


out = {}
meta = {}
max_delta = 0.0
for i in cov:
    p = d1["pairs"][i]
    rec = by_construct[normalize_construct_name(p["var_label"])]
    assert p["pos"] == rec["pos"]
    we_seq = erase_at(rec["wt_seq"], p["pos"])
    me_seq = erase_at(rec["mt_seq"], p["pos"])
    assert we_seq == me_seq, f"erased strings differ for pair {i}"
    we = hidden(we_seq)
    me = hidden(me_seq)
    d = float(np.abs(we - me).max())
    max_delta = max(max_delta, d)
    assert d <= 1e-5, f"erased embedding delta {d} exceeds 1e-5 for pair {i}"
    out[f"ewt_{i}"] = we.astype(np.float32)
    meta[str(i)] = {
        "parent": p["parent"], "mutation": p["mutation"], "pos": p["pos"],
        "erased_len": len(we_seq),
        "erased_string_sha256": hashlib.sha256(we_seq.encode()).hexdigest(),
        "embedding_abs_delta": d,
    }
    print(f"pair {i} {p['parent']} {p['mutation']} pos={p['pos']} {we.shape} delta={d}", flush=True)

assert max_delta == 0.0, "identical strings must give identical embeddings"
np.savez_compressed(HERE / "ERASED_ESM_S1.npz", **out)
(HERE / "ERASED_S1_META.json").write_text(json.dumps(
    {"pairs": meta, "max_embedding_abs_delta": max_delta,
     "n_pairs": len(cov), "asserts": ["strings equal", "delta <= 1e-5", "delta == 0"]},
    indent=1), encoding="utf-8")
sha = hashlib.sha256((HERE / "ERASED_ESM_S1.npz").read_bytes()).hexdigest()
(HERE / "ERASED_S1_SHA256.txt").write_text(sha + "\n")
print("wrote ERASED_ESM_S1.npz:", len(cov), "pairs; sha256", sha, "; max delta", max_delta)
