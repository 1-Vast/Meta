"""S7_L2B runner — trainability control, B0/BL/B4 arms, Gates G1-G5.

Order enforced by the global constraints:
  Data/Provenance -> Identifiability audit -> Matched baseline -> Preregistration
  -> Trainability control -> Development training -> Gates.

The preregistration was frozen and committed (ce186f4) before this file existed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_dataset import MONN, build, make_split, protein_components  # noqa: E402
from s7_localizer import (FeatureStore, average_precision, build_controls,  # noqa: E402
                          component_macro, evaluate, paired_bootstrap,
                          train_arm)

OUT = ROOT / "report" / "s7_l2b_r0r"
PRED = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "preds"
PREREG_SHA = "2c333f223ae450c566cc62b1a3b276ff59c065c38348005ad9504ac1930b9a92"
GATE = 0.02


def load_mols():
    md = [pickle.load((MONN / "mol_dict").open("rb"), encoding="bytes"),
          pickle.load((MONN / "independent_dataset_mol_dict").open("rb"),
                      encoding="bytes")]
    cache = {}

    class C(dict):
        def __missing__(self, ccd):
            for d in md:
                m = d.get(ccd.encode("ascii", "ignore"))
                if m is not None:
                    from rdkit import Chem
                    w = Chem.Mol(m)
                    Chem.SanitizeMol(w)
                    self[ccd] = w
                    return w
            raise KeyError(ccd)
    return C(cache)


def evaluator_selftest():
    """Metric contract test: AP must be order-invariant and tie-deterministic."""
    rng = np.random.default_rng(0)
    L, A = 7, 5
    y = (rng.random((L, A)) < 0.2).astype(np.int8)
    if y.sum() == 0:
        y[0, 0] = 1
    s = np.round(rng.random((L, A)), 1)          # many ties by construction
    ri = np.repeat(np.arange(L), A)
    ai = np.tile(np.arange(A), L)
    ap1 = average_precision(s.ravel(), y.ravel(), ri, ai)
    perm = rng.permutation(L * A)
    ap2 = average_precision(s.ravel()[perm], y.ravel()[perm], ri[perm], ai[perm])
    allpos = average_precision(np.ones(10), np.ones(10, np.int8),
                               np.arange(10), np.arange(10))
    none = average_precision(np.ones(4), np.zeros(4, np.int8),
                             np.arange(4), np.arange(4))
    return {"tie_order_invariance": {"ap_natural": ap1, "ap_permuted": ap2,
                                     "identical": bool(ap1 == ap2)},
            "all_positive_ap": allpos,
            "no_positive_returns_none": none is None,
            "verdict": "PASS" if (ap1 == ap2 and allpos == 1.0 and none is None)
                       else "S7L2B_EVALUATOR_CONTRACT_FAIL_CLOSED"}


def trainability_control(train, held, store, mols, device, comp_of, n=400):
    """Fit a KNOWN function of the frozen inputs under the identical pipeline."""
    torch.manual_seed(7)
    g = torch.Generator().manual_seed(7)
    from s7_dataset import ATOM_DIM, RES_DIM
    wr = torch.randn(RES_DIM, generator=g).numpy()
    wa = torch.randn(ATOM_DIM, generator=g).numpy()
    sub_tr, sub_te = train[:n], held[:max(60, n // 4)]

    def synth(rec):
        xr, xa = store.residues(rec), store.atoms(rec, mols)
        s = np.outer(xr @ wr, xa @ wa)
        thr = np.quantile(s, 0.995)
        return (s >= thr).astype(np.int8)

    orig = {}
    for rec in sub_tr + sub_te:
        orig[rec["source_key"]] = rec["edges"]
        yy = synth(rec)
        rec["edges"] = [(int(i), int(j)) for i, j in zip(*np.nonzero(yy))][:400]
    keep_tr = [r for r in sub_tr if r["edges"]]
    keep_te = [r for r in sub_te if r["edges"]]
    model, _ = train_arm(keep_tr, store, mols, True, device)
    _pc, comp_m, macro, _ = evaluate(model, keep_te, store, mols, device, comp_of)
    _pc0, comp_0, macro0, _ = evaluate(model, keep_te, store, mols, device, comp_of,
                                       prevalence=True)
    for rec in sub_tr + sub_te:
        rec["edges"] = orig[rec["source_key"]]
    return {"synthetic_macro_ap": macro, "prevalence_macro_ap": macro0,
            "recovered": bool(macro - macro0 > 0.05),
            "verdict": "PASS" if macro - macro0 > 0.05 else "S7L2B_TRAINABILITY_FAIL"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit-train", type=int, default=0)
    ap.add_argument("--limit-eval", type=int, default=0)
    ap.add_argument("--skip-synthetic", action="store_true")
    args = ap.parse_args()
    device = "cuda" if torch.cuda.is_available() else "cpu"
    OUT.mkdir(parents=True, exist_ok=True)
    PRED.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    self_test = evaluator_selftest()
    print("evaluator self-test:", json.dumps(self_test["verdict"]), flush=True)
    if self_test["verdict"] != "PASS":
        (OUT / "S7L2B_DEVELOPMENT_GATE.json").write_text(
            json.dumps({"terminal": self_test["verdict"], "self_test": self_test},
                       indent=2), encoding="utf-8")
        return 1

    kept, quarantine, contract, _ = build()
    comp_of = protein_components(kept)
    train, held_all, held_A, held_B = make_split(kept, comp_of)
    if args.limit_train:
        train = train[:args.limit_train]
    if args.limit_eval:
        held_A = held_A[:args.limit_eval]
        held_B = held_B[:args.limit_eval]
    print(f"train={len(train)} heldA={len(held_A)} heldB={len(held_B)} dev={device}",
          flush=True)

    store, mols = FeatureStore(), load_mols()
    result = {"schema": "MetaSieve.S7L2B.DevelopmentGate.v1",
              "created_utc": "2026-08-09",
              "preregistration_sha256": PREREG_SHA,
              "preregistration_commit": "ce186f4",
              "device": device,
              "evaluator_self_test": self_test,
              "atom_contract_census": contract,
              "quarantined": len(quarantine),
              "split": {"train": len(train), "heldout_A": len(held_A),
                        "heldout_B": len(held_B),
                        "train_positives": sum(len(r["edges"]) for r in train),
                        "heldout_A_positives": sum(len(r["edges"]) for r in held_A)}}

    if not args.skip_synthetic:
        print("trainability control ...", flush=True)
        result["trainability_control"] = trainability_control(
            train, held_A, store, mols, device, comp_of)
        print("  ", json.dumps(result["trainability_control"]), flush=True)
        if result["trainability_control"]["verdict"] != "PASS":
            result["terminal"] = "S7L2B_TRAINABILITY_FAIL"
            (OUT / "S7L2B_DEVELOPMENT_GATE.json").write_text(
                json.dumps(result, indent=2), encoding="utf-8")
            return 1

    print("training B4 (non-PLM residue baseline) ...", flush=True)
    b4, trace4 = train_arm(train, store, mols, True, device, log=True)
    print("training BL (ligand-only) ...", flush=True)
    bl, traceL = train_arm(train, store, mols, False, device, log=True)

    prot_map, lig_map = build_controls(held_A, train, comp_of, store, mols)
    by_key = {r["source_key"]: r for r in train}
    result["control_maps"] = {
        "wrong_protein_coverage": float(np.mean([v is not None for v in prot_map.values()])),
        "wrong_ligand_coverage": float(np.mean([v is not None for v in lig_map.values()])),
        "wrong_protein_map_sha256": sha_map(prot_map),
        "wrong_ligand_map_sha256": sha_map(lig_map),
    }

    def res_wrong(rec):
        k = prot_map.get(rec["source_key"])
        if k is None:
            return store.residues(rec)
        return store.residues(by_key[k])[:rec["n_res"]]

    def atom_wrong(rec):
        k = lig_map.get(rec["source_key"])
        if k is None:
            return store.atoms(rec, mols)
        return store.atoms(by_key[k], mols)

    arms = {}
    print("evaluating arms on held-out A (complete matrix) ...", flush=True)
    for name, kw, model in (
            ("B0", {"prevalence": True}, b4),
            ("BL", {}, bl),
            ("B4", {}, b4),
            ("BP", {"res_override": res_wrong}, b4),
            ("BX", {"atom_override": atom_wrong}, b4),
            ("BM", {"shuffle_residues": True}, b4)):
        pc, cm, macro, rows = evaluate(model, held_A, store, mols, device, comp_of, **kw)
        arms[name] = {"comp": cm, "macro": macro, "rows": rows}
        print(f"  {name:3s} macro-AP = {macro:.5f}  components={len(cm)}", flush=True)

    for name in arms:
        p = PRED / f"heldoutA_{name}.json"
        p.write_text(json.dumps(arms[name]["rows"]), encoding="utf-8")

    result["heldout_A"] = {n: {"macro_ap": arms[n]["macro"],
                               "components": len(arms[n]["comp"])} for n in arms}
    gates = {}
    for gname, a, b in (("G1_B4_vs_B0", "B4", "B0"), ("G2_B4_vs_BL", "B4", "BL"),
                        ("G3_B4_vs_BP", "B4", "BP"), ("G4_B4_vs_BM", "B4", "BM"),
                        ("G5_B4_vs_BX", "B4", "BX")):
        bs = paired_bootstrap(arms[a]["comp"], arms[b]["comp"])
        bs["meets_threshold"] = bool(bs["delta"] >= GATE)
        bs["lcb_above_zero"] = bool(bs["lcb95_one_sided"] > 0)
        bs["pass"] = bool(bs["meets_threshold"] and bs["lcb_above_zero"])
        gates[gname] = bs
        print(f"  {gname:14s} delta={bs['delta']:+.5f} lcb={bs['lcb95_one_sided']:+.5f} "
              f"{'PASS' if bs['pass'] else 'FAIL'}", flush=True)
    result["gates_heldout_A"] = gates
    result["all_gates_pass"] = all(g["pass"] for g in gates.values())

    print("evaluating held-out B (scaffold-strict robustness) ...", flush=True)
    bB = {}
    for name, kw, model in (("BL", {}, bl), ("B4", {}, b4)):
        _pc, cm, macro, _rows = evaluate(model, held_B, store, mols, device, comp_of, **kw)
        bB[name] = {"comp": cm, "macro": macro}
    result["heldout_B"] = {n: bB[n]["macro"] for n in bB}
    result["heldout_B_B4_minus_BL"] = paired_bootstrap(bB["B4"]["comp"], bB["BL"]["comp"])

    b4_macro = arms["B4"]["macro"]
    g2 = gates["G2_B4_vs_BL"]
    if not result["all_gates_pass"]:
        if not g2["pass"]:
            terminal = "S7L2B_LIGAND_ONLY_SHORTCUT"
        else:
            terminal = "S7L2B_PARTNER_CONTROL_FAIL"
    elif b4_macro < 0.10:
        terminal = "S7L2B_BASELINE_ESTABLISHED_PLM_INDICATED"
    else:
        terminal = "S7L2B_BASELINE_ESTABLISHED_PLM_NOT_INDICATED"
    result["terminal"] = terminal
    result["b5_escalation_authorised"] = bool(
        terminal == "S7L2B_BASELINE_ESTABLISHED_PLM_INDICATED"
        or (not g2["pass"]))
    result["elapsed_sec"] = round(time.time() - t0, 1)

    ckpt = PRED / "B4_checkpoint.pt"
    torch.save(b4.state_dict(), ckpt)
    result["B4_checkpoint_sha256"] = hashlib.sha256(ckpt.read_bytes()).hexdigest()
    result["training_trace_tail"] = trace4[-3:]
    (OUT / "S7L2B_DEVELOPMENT_GATE.json").write_text(json.dumps(result, indent=2),
                                                     encoding="utf-8")
    print(f"\nTERMINAL: {terminal}")
    print(f"wrote {OUT / 'S7L2B_DEVELOPMENT_GATE.json'}")
    return 0


def sha_map(m):
    return hashlib.sha256(json.dumps(m, sort_keys=True).encode()).hexdigest()


if __name__ == "__main__":
    sys.exit(main())
