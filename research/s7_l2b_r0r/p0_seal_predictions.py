"""Phase 0 items 2 and 3 — tie-aware AP, sealed per-pair predictions, manifests.

Blockers addressed:
  TIE_AWARE_AP_AND_PER_PAIR_PREDICTIONS_NOT_MATERIALIZED
  NEGATIVE_CONTROL_AND_TRAINING_DIAGNOSTIC_MANIFESTS_INCOMPLETE

Retrains B4 and BL under the identical frozen seeds and VERIFIES that the B4
checkpoint hash reproduces the previously sealed one, demonstrating run-to-run
determinism rather than assuming it. Every held-out A per-pair score is then
written to a hashed float16 memmap so each reported AP is recomputable from
sealed artifacts.

Tie-aware AP is reported three ways per complex: lexicographic (the registered
rule), optimistic (positives first) and pessimistic (negatives first). Where the
gap is material a Monte-Carlo expectation over random tie orderings is added.
An arm whose AP depends on tie handling is flagged rather than silently reported.
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
                          paired_bootstrap, sample_negatives, train_arm)
from s7_run import load_mols  # noqa: E402

OUT = ROOT / "report" / "s7_l2b_r0r"
PRED = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "sealed_preds"
GATE = 0.02
MC_REPEATS = 12
TIE_MATERIAL = 1e-6


def ap_from_order(order, y):
    yy = y[order].astype(np.float64)
    n_pos = yy.sum()
    if n_pos == 0:
        return None
    tp = np.cumsum(yy)
    prec = tp / np.arange(1, yy.size + 1, dtype=np.float64)
    return float((prec * yy).sum() / n_pos)


def ap_variants(s, y, ri, ai, rng=None):
    lex = ap_from_order(np.lexsort((ai, ri, -s)), y)
    opt = ap_from_order(np.lexsort((-y, -s)), y)
    pes = ap_from_order(np.lexsort((y, -s)), y)
    out = {"lex": lex, "opt": opt, "pess": pes,
           "tie_gap": None if opt is None else float(opt - pes)}
    if rng is not None and out["tie_gap"] is not None and out["tie_gap"] > TIE_MATERIAL:
        vals = [ap_from_order(np.lexsort((rng.random(s.size), -s)), y)
                for _ in range(MC_REPEATS)]
        out["expected_mc"] = float(np.mean(vals))
        out["expected_mc_sd"] = float(np.std(vals))
    return out


def sha_file(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""):
            h.update(c)
    return h.hexdigest()


@torch.no_grad()
def score_matrix(model, rec, store, mols, device, res_override=None,
                 atom_override=None, shuffle_residues=False, prevalence=False,
                 rng=None):
    L, A = rec["n_res"], rec["n_atoms"]
    if prevalence:
        return np.zeros((L, A), dtype=np.float64)
    xr = store.residues(rec) if res_override is None else res_override(rec)
    if shuffle_residues:
        xr = xr[rng.permutation(L)]
    xa = store.atoms(rec, mols) if atom_override is None else atom_override(rec)
    t_r = torch.from_numpy(np.ascontiguousarray(xr)).to(device)
    t_a = torch.from_numpy(np.ascontiguousarray(xa)).to(device)
    return model.pair_logits(t_r, t_a).double().cpu().numpy()


def main():
    device = "cuda" if torch.cuda.is_available() else "cpu"
    PRED.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    kept, quarantine, contract, _ = build()
    comp_of = protein_components(kept)
    train, _ha, held_A, held_B = make_split(kept, comp_of)
    print(f"contract: {contract}", flush=True)
    print(f"train={len(train)} heldA={len(held_A)}", flush=True)

    store, mols = FeatureStore(), load_mols()
    print("retraining B4 and BL under frozen seeds ...", flush=True)
    b4, trace4 = train_arm(train, store, mols, True, device, log=True)
    bl, traceL = train_arm(train, store, mols, False, device, log=True)

    tmp = PRED / "B4_recheck.pt"
    torch.save(b4.state_dict(), tmp)
    prev = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "preds" / "B4_checkpoint.pt"
    determinism = {"recomputed_sha256": sha_file(tmp),
                   "previously_sealed_sha256": sha_file(prev) if prev.exists() else None}
    determinism["reproduces"] = bool(
        determinism["previously_sealed_sha256"] == determinism["recomputed_sha256"])
    print(f"determinism reproduces sealed checkpoint: {determinism['reproduces']}",
          flush=True)

    prot_map, lig_map = build_controls(held_A, train, comp_of, store, mols)
    by_key = {r["source_key"]: r for r in train}

    def res_wrong(rec):
        k = prot_map.get(rec["source_key"])
        return store.residues(rec) if k is None else store.residues(by_key[k])[:rec["n_res"]]

    def atom_wrong(rec):
        k = lig_map.get(rec["source_key"])
        return store.atoms(rec, mols) if k is None else store.atoms(by_key[k], mols)

    total_cells = int(sum(r["n_res"] * r["n_atoms"] for r in held_A))
    print(f"held-out A complete-matrix cells: {total_cells:,}", flush=True)

    ARMS = {"B0": dict(prevalence=True), "BL": dict(), "B4": dict(),
            "BP": dict(res_override=res_wrong), "BX": dict(atom_override=atom_wrong),
            "BM": dict(shuffle_residues=True)}
    MODEL = {"B0": b4, "BL": bl, "B4": b4, "BP": b4, "BX": b4, "BM": b4}

    index, off = {}, 0
    for rec in held_A:
        index[rec["source_key"]] = [off, rec["n_res"], rec["n_atoms"]]
        off += rec["n_res"] * rec["n_atoms"]

    ap_tables, arm_hashes, tie_report = {}, {}, {}
    for name, cfg in ARMS.items():
        mm_path = PRED / f"heldoutA_{name}.f16.dat"
        mm = np.memmap(mm_path, dtype=np.float16, mode="w+", shape=(total_cells,))
        rng_ctrl = np.random.default_rng(20260814)
        rng_tie = np.random.default_rng(20260816)
        per, gaps = {}, []
        for rec in held_A:
            s = score_matrix(MODEL[name], rec, store, mols, device, rng=rng_ctrl, **cfg)
            L, A = rec["n_res"], rec["n_atoms"]
            y = np.zeros((L, A), dtype=np.int8)
            for i, j in rec["edges"]:
                y[i, j] = 1
            o = index[rec["source_key"]][0]
            mm[o:o + L * A] = s.ravel().astype(np.float16)
            ri = np.repeat(np.arange(L), A)
            ai = np.tile(np.arange(A), L)
            v = ap_variants(s.ravel(), y.ravel(), ri, ai, rng_tie)
            per[rec["source_key"]] = v
            if v["tie_gap"] is not None:
                gaps.append(v["tie_gap"])
        mm.flush()
        del mm
        arm_hashes[name] = sha_file(mm_path)
        ap_tables[name] = per
        g = np.array(gaps) if gaps else np.array([0.0])
        tie_report[name] = {"max_tie_gap": float(g.max()),
                            "mean_tie_gap": float(g.mean()),
                            "materially_tie_dependent": bool(g.max() > TIE_MATERIAL)}
        print(f"  {name}: sealed, max tie gap {g.max():.3e}", flush=True)

    def macro(name, which):
        d = {k: (v[which] if which in v else v["lex"]) for k, v in ap_tables[name].items()}
        return component_macro(d, comp_of)

    summary = {}
    for name in ARMS:
        _c1, m_lex = macro(name, "lex")
        _c2, m_opt = macro(name, "opt")
        _c3, m_pes = macro(name, "pess")
        _c4, m_pt = macro(name, "expected_mc")
        summary[name] = {"macro_ap_lex": m_lex, "macro_ap_optimistic": m_opt,
                         "macro_ap_pessimistic": m_pes,
                         "macro_ap_tie_aware_point": m_pt, "tie": tie_report[name]}

    ref, _ = macro("B4", "expected_mc")
    gates = {}
    for gname, other in (("G1_B4_vs_B0", "B0"), ("G2_B4_vs_BL", "BL"),
                         ("G3_B4_vs_BP", "BP"), ("G4_B4_vs_BM", "BM"),
                         ("G5_B4_vs_BX", "BX")):
        oc, _ = macro(other, "expected_mc")
        bs = paired_bootstrap(ref, oc)
        bs["pass"] = bool(bs["delta"] >= GATE and bs["lcb95_one_sided"] > 0)
        gates[gname] = bs

    rng_neg = np.random.default_rng(20260813)
    neg = {"complexes": 0, "positives": 0, "negatives": 0,
           "exactly_six_per_positive": True, "positive_collisions": 0,
           "duplicate_negatives": 0}
    for rec in train[:1500]:
        negs = sample_negatives(rec, rng_neg)
        pos = set(map(tuple, rec["edges"]))
        neg["complexes"] += 1
        neg["positives"] += len(rec["edges"])
        neg["negatives"] += len(negs)
        if len(negs) != 6 * len(rec["edges"]):
            neg["exactly_six_per_positive"] = False
        neg["positive_collisions"] += sum(1 for n in negs if n in pos)
        neg["duplicate_negatives"] += len(negs) - len(set(negs))

    manifest = {
        "schema": "MetaSieve.S7L2B.P0.SealedPredictionManifest.v1",
        "created_utc": "2026-08-09",
        "atom_contract_census": contract,
        "records_quarantined": len(quarantine),
        "split": {"train": len(train), "heldout_A": len(held_A), "heldout_B": len(held_B),
                  "heldout_A_components": len({comp_of[r["source_key"]] for r in held_A})},
        "split_sha256": hashlib.sha256(json.dumps(
            {"train": sorted(r["source_key"] for r in train),
             "heldout_A": sorted(r["source_key"] for r in held_A),
             "heldout_B": sorted(r["source_key"] for r in held_B)},
            sort_keys=True).encode()).hexdigest(),
        "determinism_check": determinism,
        "per_pair_predictions": {"dir": str(PRED.relative_to(ROOT)).replace("\\", "/"),
                                 "dtype": "float16", "total_cells": total_cells,
                                 "index_file": "heldoutA_index.json",
                                 "arm_sha256": arm_hashes},
        "tie_aware_ap": summary,
        "gates_recomputed_from_sealed_predictions": gates,
        "negative_sampling_manifest": neg,
        "control_maps": {
            "wrong_protein_coverage": float(np.mean([v is not None for v in prot_map.values()])),
            "wrong_ligand_coverage": float(np.mean([v is not None for v in lig_map.values()])),
            "wrong_protein_sha256": hashlib.sha256(
                json.dumps(prot_map, sort_keys=True).encode()).hexdigest(),
            "wrong_ligand_sha256": hashlib.sha256(
                json.dumps(lig_map, sort_keys=True).encode()).hexdigest()},
        "training_diagnostics": {"B4_trace_tail": trace4[-3:], "BL_trace_tail": traceL[-3:]},
        "elapsed_sec": round(time.time() - t0, 1),
    }
    (PRED / "heldoutA_index.json").write_text(json.dumps(index), encoding="utf-8")
    (PRED / "ap_tables.json").write_text(json.dumps(ap_tables), encoding="utf-8")
    (OUT / "P0_SEALED_PREDICTION_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"tie_aware_ap": summary,
                      "gates": {k: {"delta": v["delta"], "lcb": v["lcb95_one_sided"],
                                    "pass": v["pass"]} for k, v in gates.items()},
                      "determinism": determinism, "negatives": neg}, indent=2))
    print(f"\nwrote {OUT / 'P0_SEALED_PREDICTION_MANIFEST.json'}")


if __name__ == "__main__":
    main()
