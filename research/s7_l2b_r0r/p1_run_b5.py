"""Phase 1 — the registered B5 discriminator: frozen ESM2-650M residue features.

Changes ONLY the residue representation. The atom branch, head, rank 32,
projected dim 128, negative sampler, optimizer, learning rate, weight decay,
epochs, seeds, split, evaluation mask, tie policy and Gates are identical to B4.

No attention, no geometry, no typed interactions, no affinity supervision, no
larger PLM, no additional capacity.

B4 and BL are taken from the Phase-0 SEALED tables so every contrast uses the
same post-quarantine rows and the same tie-aware point estimator.
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
from s7_localizer import (FeatureStore, build_controls, component_macro,  # noqa: E402
                          paired_bootstrap, train_arm)
from s7_run import load_mols  # noqa: E402
from p0_seal_predictions import ap_variants, score_matrix, sha_file  # noqa: E402

OUT = ROOT / "report" / "s7_l2b_r0r"
SEALED = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "sealed_preds"
ESM = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "esm2_650M"
PRED5 = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "sealed_preds_b5"
GATE = 0.02
DIM = 1280


class EsmStore(FeatureStore):
    def __init__(self, mm, index):
        super().__init__()
        self.mm, self.index = mm, index

    def residues(self, rec):
        off, L = self.index[rec["seq_key"]]
        return np.asarray(self.mm[off:off + L], dtype=np.float32)


def main():
    device = "cuda" if torch.cuda.is_available() else "cpu"
    PRED5.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    index = json.loads((ESM / "esm2_650M_index.json").read_text())
    total = max(v[0] + v[1] for v in index.values())
    mm = np.memmap(ESM / "esm2_650M_residues.fp16.dat", dtype=np.float16,
                   mode="r", shape=(total, DIM))

    kept, _q, contract, _f = build()
    comp_of = protein_components(kept)
    train, _ha, held_A, held_B = make_split(kept, comp_of)
    train = [r for r in train if r["seq_key"] in index]
    held_A = [r for r in held_A if r["seq_key"] in index]
    held_B = [r for r in held_B if r["seq_key"] in index]
    print(f"train={len(train)} heldA={len(held_A)} heldB={len(held_B)} dev={device}",
          flush=True)

    store, mols = EsmStore(mm, index), load_mols()
    print("training B5 (frozen ESM2-650M residue features) ...", flush=True)
    b5, trace = train_arm(train, store, mols, True, device, log=True, d_res_in=DIM)

    prot_map, lig_map = build_controls(held_A, train, comp_of, store, mols)
    by_key = {r["source_key"]: r for r in train}

    def res_wrong(rec):
        k = prot_map.get(rec["source_key"])
        return store.residues(rec) if k is None else store.residues(by_key[k])[:rec["n_res"]]

    def atom_wrong(rec):
        k = lig_map.get(rec["source_key"])
        return store.atoms(rec, mols) if k is None else store.atoms(by_key[k], mols)

    ARMS = {"B5": dict(), "BP5": dict(res_override=res_wrong),
            "BX5": dict(atom_override=atom_wrong), "BM5": dict(shuffle_residues=True)}
    total_cells = int(sum(r["n_res"] * r["n_atoms"] for r in held_A))
    idx, off = {}, 0
    for rec in held_A:
        idx[rec["source_key"]] = [off, rec["n_res"], rec["n_atoms"]]
        off += rec["n_res"] * rec["n_atoms"]

    tables, hashes = {}, {}
    for name, cfg in ARMS.items():
        p = PRED5 / f"heldoutA_{name}.f16.dat"
        arr = np.memmap(p, dtype=np.float16, mode="w+", shape=(total_cells,))
        rng_c = np.random.default_rng(20260814)
        rng_t = np.random.default_rng(20260816)
        per = {}
        for rec in held_A:
            s = score_matrix(b5, rec, store, mols, device, rng=rng_c, **cfg)
            L, A = rec["n_res"], rec["n_atoms"]
            y = np.zeros((L, A), dtype=np.int8)
            for i, j in rec["edges"]:
                y[i, j] = 1
            o = idx[rec["source_key"]][0]
            arr[o:o + L * A] = s.ravel().astype(np.float16)
            per[rec["source_key"]] = ap_variants(
                s.ravel(), y.ravel(), np.repeat(np.arange(L), A),
                np.tile(np.arange(A), L), rng_t)
        arr.flush()
        del arr
        hashes[name] = sha_file(p)
        tables[name] = per
        _cm, m = component_macro({k: v.get("expected_mc", v["lex"])
                                  for k, v in per.items()}, comp_of)
        print(f"  {name:4s} macro-AP = {m:.6f}", flush=True)

    sealed = json.loads((SEALED / "ap_tables.json").read_text())

    def comp_of_table(tab):
        return component_macro({k: v.get("expected_mc", v["lex"])
                                for k, v in tab.items()}, comp_of)

    macros = {}
    for n, tab in list(tables.items()) + [(n, sealed[n]) for n in
                                          ("B0", "BL", "B4", "BP", "BX", "BM")]:
        cm, m = comp_of_table(tab)
        macros[n] = (cm, m)

    gates = {}
    for gname, other in (("G1_B5_vs_B0", "B0"), ("G2_B5_vs_BL", "BL"),
                         ("G3_B5_vs_BP", "BP5"), ("G4_B5_vs_BM", "BM5"),
                         ("G5_B5_vs_BX", "BX5"), ("G6_B5_vs_B4", "B4")):
        bs = paired_bootstrap(macros["B5"][0], macros[other][0])
        bs["pass"] = bool(bs["delta"] >= GATE and bs["lcb95_one_sided"] > 0)
        gates[gname] = bs
        print(f"  {gname:14s} delta={bs['delta']:+.6f} lcb={bs['lcb95_one_sided']:+.6f} "
              f"{'PASS' if bs['pass'] else 'FAIL'}", flush=True)

    ck = PRED5 / "B5_checkpoint.pt"
    torch.save(b5.state_dict(), ck)
    allpass = all(g["pass"] for g in gates.values())
    res = {
        "schema": "MetaSieve.S7L2B.P1.B5Gate.v1",
        "created_utc": "2026-08-09",
        "preregistration_sha256": "2c333f223ae450c566cc62b1a3b276ff59c065c38348005ad9504ac1930b9a92",
        "preregistration_commit": "ce186f4",
        "phase0_commit": "139effd",
        "residue_representation": "frozen esm2_t33_650M_UR50D final layer 1280-d, "
                                  "weights sha256 c874668852c7275a159e2c7ceb6069671d7b1ba2c7b52f59600b34ce0f721008",
        "only_change_vs_B4": "residue features",
        "split": {"train": len(train), "heldout_A": len(held_A),
                  "components": len({comp_of[r['source_key']] for r in held_A})},
        "macro_ap_tie_aware": {n: macros[n][1] for n in macros},
        "gates": gates,
        "all_gates_pass": allpass,
        "per_pair_prediction_sha256": hashes,
        "B5_checkpoint_sha256": hashlib.sha256(ck.read_bytes()).hexdigest(),
        "training_trace_tail": trace[-3:],
        "elapsed_sec": round(time.time() - t0, 1),
        "confirmation_cohort": "SEALED — not opened",
    }
    (PRED5 / "ap_tables_b5.json").write_text(json.dumps(tables), encoding="utf-8")
    (OUT / "P1_B5_GATE.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("\nmacro-AP:", json.dumps({n: round(macros[n][1], 6) for n in macros}, indent=1))
    print(f"ALL GATES PASS: {allpass}")
    print(f"wrote {OUT / 'P1_B5_GATE.json'}")


if __name__ == "__main__":
    main()
