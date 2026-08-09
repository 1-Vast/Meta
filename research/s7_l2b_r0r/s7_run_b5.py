"""S7_L2B — B5 arm: frozen ESM2-650M residue localizer.

Authorised by the registered escalation rule in PREREG_S7_L2B_UNIFIED.md section 7,
which fired on the measured B4 result (B4-BL below the 0.02 effect size and B4
macro-AP 0.0229 << 0.10).

B5 differs from B4 ONLY in the residue representation. Same head, same rank 32,
same projected dimension 128, same negative sampler, same epochs, optimiser,
learning rate, weight decay, seeds, same evaluation mask and the SAME frozen
wrong-partner control maps.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_dataset import MONN, build, make_split, protein_components  # noqa: E402
from s7_localizer import (FeatureStore, build_controls, component_macro,  # noqa: E402
                          evaluate, paired_bootstrap, train_arm)
from s7_run import load_mols  # noqa: E402

OUT = ROOT / "report" / "s7_l2b_r0r"
PRED = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "preds"
ESM = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "esm2_650M"
GATE = 0.02
DIM = 1280


class EsmStore(FeatureStore):
    """Residue features come from the frozen ESM2 cache; atoms unchanged."""

    def __init__(self, memmap, index):
        super().__init__()
        self.mm = memmap
        self.index = index

    def residues(self, rec):
        off, L = self.index[rec["seq_key"]]
        return np.asarray(self.mm[off:off + L], dtype=np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit-train", type=int, default=0)
    args = ap.parse_args()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    t0 = time.time()

    index = json.loads((ESM / "esm2_650M_index.json").read_text())
    total = max(v[0] + v[1] for v in index.values())
    mm = np.memmap(ESM / "esm2_650M_residues.fp16.dat", dtype=np.float16,
                   mode="r", shape=(total, DIM))
    kept, _q, _c, _f = build()
    comp_of = protein_components(kept)
    train, _ha, held_A, held_B = make_split(kept, comp_of)
    if args.limit_train:
        train = train[:args.limit_train]
    train = [r for r in train if r["seq_key"] in index]
    held_A = [r for r in held_A if r["seq_key"] in index]
    held_B = [r for r in held_B if r["seq_key"] in index]
    print(f"train={len(train)} heldA={len(held_A)} heldB={len(held_B)}", flush=True)

    store, mols = EsmStore(mm, index), load_mols()
    print("training B5 (frozen ESM2-650M residue representation) ...", flush=True)
    b5, trace = train_arm(train, store, mols, True, device, log=True)

    prot_map, lig_map = build_controls(held_A, train, comp_of, store, mols)
    by_key = {r["source_key"]: r for r in train}

    def res_wrong(rec):
        k = prot_map.get(rec["source_key"])
        return store.residues(rec) if k is None else store.residues(by_key[k])[:rec["n_res"]]

    def atom_wrong(rec):
        k = lig_map.get(rec["source_key"])
        return store.atoms(rec, mols) if k is None else store.atoms(by_key[k], mols)

    arms = {}
    for name, kw in (("B0", {"prevalence": True}), ("B5", {}),
                     ("BP", {"res_override": res_wrong}),
                     ("BX", {"atom_override": atom_wrong}),
                     ("BM", {"shuffle_residues": True})):
        _pc, cm, macro, rows = evaluate(b5, held_A, store, mols, device, comp_of, **kw)
        arms[name] = {"comp": cm, "macro": macro}
        (PRED / f"heldoutA_B5_{name}.json").write_text(json.dumps(rows), encoding="utf-8")
        print(f"  {name:3s} macro-AP = {macro:.5f}", flush=True)

    # B4 and BL component means from the sealed B4-run predictions
    def comp_from_rows(path):
        rows = json.loads(Path(path).read_text())
        per = {r["key"]: r["ap"] for r in rows}
        cof = {r["key"]: r["component"] for r in rows}
        cm, macro = component_macro(per, cof)
        return cm, macro

    b4_comp, b4_macro = comp_from_rows(PRED / "heldoutA_B4.json")
    bl_comp, bl_macro = comp_from_rows(PRED / "heldoutA_BL.json")
    print(f"  B4 (sealed) macro-AP = {b4_macro:.5f}", flush=True)
    print(f"  BL (sealed) macro-AP = {bl_macro:.5f}", flush=True)

    gates = {}
    for gname, a_comp, b_comp in (
            ("G1_B5_vs_B0", arms["B5"]["comp"], arms["B0"]["comp"]),
            ("G2_B5_vs_BL", arms["B5"]["comp"], bl_comp),
            ("G3_B5_vs_BP", arms["B5"]["comp"], arms["BP"]["comp"]),
            ("G4_B5_vs_BM", arms["B5"]["comp"], arms["BM"]["comp"]),
            ("G5_B5_vs_BX", arms["B5"]["comp"], arms["BX"]["comp"]),
            ("G6_B5_vs_B4", arms["B5"]["comp"], b4_comp)):
        bs = paired_bootstrap(a_comp, b_comp)
        bs["meets_threshold"] = bool(bs["delta"] >= GATE)
        bs["lcb_above_zero"] = bool(bs["lcb95_one_sided"] > 0)
        bs["pass"] = bool(bs["meets_threshold"] and bs["lcb_above_zero"])
        gates[gname] = bs
        print(f"  {gname:14s} delta={bs['delta']:+.5f} lcb={bs['lcb95_one_sided']:+.5f} "
              f"{'PASS' if bs['pass'] else 'FAIL'}", flush=True)

    _pcB, cmB, macroB, _rB = evaluate(b5, held_B, store, mols, device, comp_of)
    allpass = all(g["pass"] for g in gates.values())
    terminal = ("S7L2B_EXACT_RESIDUE_LOCALIZATION_DEVELOPMENT_PASS" if allpass
                else "S7L2B_PLM_BELOW_GATE")

    ck = PRED / "B5_checkpoint.pt"
    torch.save(b5.state_dict(), ck)
    res = {"schema": "MetaSieve.S7L2B.B5DevelopmentGate.v1",
           "created_utc": "2026-08-09",
           "preregistration_sha256": "2c333f223ae450c566cc62b1a3b276ff59c065c38348005ad9504ac1930b9a92",
           "preregistration_commit": "ce186f4",
           "escalation_rule_fired": True,
           "residue_representation": "frozen esm2_t33_650M_UR50D final layer, 1280-d",
           "heldout_A_macro_ap": {k: arms[k]["macro"] for k in arms}
                                 | {"B4_sealed": b4_macro, "BL_sealed": bl_macro},
           "gates_heldout_A": gates,
           "all_gates_pass": allpass,
           "heldout_B_B5_macro_ap": macroB,
           "B5_checkpoint_sha256": hashlib.sha256(ck.read_bytes()).hexdigest(),
           "training_trace_tail": trace[-3:],
           "terminal": terminal,
           "confirmation_cohort": "SEALED - not opened regardless of this outcome; "
                                  "R0R-3 publication/time closure is still unbuilt",
           "elapsed_sec": round(time.time() - t0, 1)}
    (OUT / "S7L2B_B5_GATE.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(f"\nTERMINAL: {terminal}")
    print(f"wrote {OUT / 'S7L2B_B5_GATE.json'}")


if __name__ == "__main__":
    main()
