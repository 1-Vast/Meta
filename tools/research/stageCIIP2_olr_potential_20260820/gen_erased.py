"""Generate ERASED_ESM.npz: X-substituted residue states for all 49 covered
pairs (mutation-erasure control, prereg T9). Same ESM model and extraction
rule as the X0c cache (facebook/esm2_t30_150M_UR50D, max_length=1022)."""
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
from x0_i2 import load_esm, build_pair_records  # noqa: E402

d1 = json.loads((BRIDGE / "DATA1A.json").read_text(encoding="utf-8"))
d2 = json.loads((BRIDGE / "DATA2X2.json").read_text(encoding="utf-8"))
pair_table = json.loads((SIG / "X0_PAIR_TABLE.json").read_text(encoding="utf-8"))
_, _, sequences = load_duongly()
records = build_pair_records(pair_table, sequences)
by_construct = {normalize_construct_name(r["construct"]): r for r in records}
cov = list(d2["covered_pair_indices"])

tok, model, device = load_esm("cuda" if torch.cuda.is_available() else "cpu")
print("device:", device, flush=True)


def erase_at(seq, pos):
    return seq[:pos - 1] + "X" + seq[pos:]


def hidden(seq):
    enc = tok(seq, return_tensors="pt", truncation=True, max_length=1022)
    enc = {k: v.to(device) for k, v in enc.items()}
    with torch.no_grad():
        out = model(**enc)
    return out.last_hidden_state[0].cpu().float().numpy()


out = {}
for i in cov:
    p = d1["pairs"][i]
    rec = by_construct[normalize_construct_name(p["var_label"])]
    assert p["pos"] == rec["pos"]
    we = hidden(erase_at(rec["wt_seq"], p["pos"]))
    me = hidden(erase_at(rec["mt_seq"], p["pos"]))
    out[f"we_{i}"] = we.astype(np.float32)
    out[f"me_{i}"] = me.astype(np.float32)
    print(f"pair {i} {p['parent']} {p['mutation']} erased {we.shape}", flush=True)

np.savez_compressed(HERE / "ERASED_ESM.npz", **out)
sha = hashlib.sha256((HERE / "ERASED_ESM.npz").read_bytes()).hexdigest()
(HERE / "ERASED_SHA256.txt").write_text(sha)
print("wrote ERASED_ESM.npz:", len(cov), "pairs; sha256", sha)
