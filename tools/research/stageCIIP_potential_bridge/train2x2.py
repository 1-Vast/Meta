"""Stage CIIP-1A 2x2 root-cause diagnostic (prereg ee844b2b...).

Four cells: {KLIFS, oracle_local_esm} x {joint, centered-only} on the
matched 49-pair subset, original split assignment, identical budget.
Single seed (1). Root-cause attribution ONLY via effect contrasts;
no single-cell attribution, no CIIP-1A PASS.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import potential as POT  # noqa: E402
import train1a as T  # noqa: E402

PREREG2X2_SHA = "ee844b2b29f0009cf97c1bd18b8a92f68dcd2dc8ea1268731c740c2224f47a8b"
CELLS = ("klifs_joint", "klifs_centered", "esm_joint", "esm_centered")
EPOCHS = 200
BATCH = 512
LR = 1e-3
WD = 1e-4
LAMBDA_ABS = 1.0
DEAD_ZONE = 10.0
SEED = 1


def rng_key(cell, kind, ep=None):
    parts = ["stageCIIP2x2", kind, cell, str(SEED)]
    if ep is not None:
        parts += ["epoch", str(ep)]
    return T.stable_rng(*parts)


def load():
    d1 = json.loads((HERE / "DATA1A.json").read_text(encoding="utf-8"))
    z1 = np.load(HERE / "DATA1A.npz", allow_pickle=False)
    d2 = json.loads((HERE / "DATA2X2.json").read_text(encoding="utf-8"))
    z2 = np.load(HERE / "DATA2X2.npz", allow_pickle=False)
    return d1, z1, d2, z2


def run_cell(cell, d1, z1, d2, z2, device):
    torch.manual_seed(SEED)
    pairs = d1["pairs"]
    targets = d1["targets"]
    covered = list(d2["covered_pair_indices"])
    sp = np.asarray(d1["split"]["pair_split"])
    tr = [i for i in covered if sp[i] == 0]
    va = [i for i in covered if sp[i] == 1]
    te = [i for i in covered if sp[i] == 2]
    is_esm = cell.startswith("esm")
    joint = cell.endswith("joint")
    model = POT.UnifiedPotential(d_p=640 if is_esm else 1700).to(device)
    # per-row protein vectors
    if is_esm:
        Pwt = torch.from_numpy(z2["esm_wt"]).float().to(device)  # (65, 640)
        Pvar = torch.from_numpy(z2["esm_var"]).float().to(device)
    else:
        Prot = torch.from_numpy(z1["prot"]).float().to(device)  # (97, 1700)
    L = torch.from_numpy(z1["lig"]).float().to(device)
    Y = torch.from_numpy(z1["Y"]).float().to(device)

    def row_vec(i, which):
        if is_esm:
            return (Pwt if which == "wt" else Pvar)[i]
        return Prot[pairs[i]["wt_row" if which == "wt" else "var_row"]]

    def pair_hat(i):
        t = targets[i]
        Lm = L[t["lig_idx"]]
        Pw = row_vec(i, "wt").unsqueeze(0).expand(len(t["lig_idx"]), -1)
        Pv = row_vec(i, "var").unsqueeze(0).expand(len(t["lig_idx"]), -1)
        return model.centered_mutation_effect(Pw, Pv, Lm)

    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],
                            lr=LR, weight_decay=WD)
    if joint:
        va_te_rows = ({pairs[i]["wt_row"] for i in va + te}
                      | {pairs[i]["var_row"] for i in va + te})
        abs_rows = sorted({pairs[i]["wt_row"] for i in tr}
                          | {pairs[i]["var_row"] for i in tr})
        if is_esm:
            abs_cells = [(i, l) for i in tr for l in targets[i]["lig_idx"]]
        else:
            abs_rows = [r for r in abs_rows if r not in va_te_rows]
            abs_cells = [(r, l) for r in abs_rows for l in range(Y.shape[1])]
    cells = [(i, j) for i in tr for j in range(len(targets[i]["c"]))]
    c_of = {i: torch.tensor(targets[i]["c"], dtype=torch.float32, device=device)
            for i in tr}
    best_val = None
    best_state = None
    gradcov = None
    for ep in range(EPOCHS):
        rng = rng_key(cell, "order", ep)
        order = np.asarray(cells)[rng.permutation(len(cells))]
        for b0 in range(0, len(order), BATCH):
            batch = order[b0:b0 + BATCH]
            opt.zero_grad()
            hats = {}
            for i in {int(i) for i, _ in batch}:
                hats[i] = pair_hat(i)
            hi = torch.stack([hats[int(i)][j] for i, j in batch])
            ci = torch.stack([c_of[int(i)][j] for i, j in batch])
            loss = ((hi - ci) ** 2).mean()
            if joint and abs_cells:
                rng_a = rng_key(cell, "abs", ep)
                sel = np.asarray(abs_cells)[rng_a.permutation(len(abs_cells))[:BATCH]]
                if is_esm:
                    ra = [pairs[int(i)]["var_row"] for i, _ in sel]
                    Pa = torch.stack([row_vec(int(i), "var") for i, _ in sel])
                    la = [t for _, t in sel]
                    f = model(Pa, L[la])
                    y = Y[ra, la]
                else:
                    ra = [r for r, _ in sel]
                    la = [l for _, l in sel]
                    fin = np.isfinite(Y[ra, la])
                    if fin.any():
                        ra = np.asarray(ra)[fin]
                        la = np.asarray(la)[fin]
                        f = model(Prot[ra], L[la])
                        y = Y[ra, la]
                    else:
                        f, y = None, None
                if f is not None:
                    loss = loss + LAMBDA_ABS * ((f - y) ** 2).mean()
            loss.backward()
            if gradcov is None:
                gradcov = {n: bool(p.grad is not None and float(p.grad.abs().max()) > 0)
                           for n, p in model.named_parameters()}
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
            opt.step()
        # val contrast mse over covered val pairs
        with torch.no_grad():
            errs = []
            for i in va:
                t = targets[i]
                h = pair_hat(i).cpu().numpy()
                errs.append(float(np.mean((h - np.asarray(t["c"])) ** 2)))
            val_mse = float(np.mean(errs))
        if best_val is None or val_mse < best_val:
            best_val = val_mse
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
    model.load_state_dict(best_state)
    # metrics on covered test pairs
    per_pair = []
    with torch.no_grad():
        for i in te:
            t = targets[i]
            h = pair_hat(i).cpu().numpy()
            c = np.asarray(t["c"])
            var_true = float(c.var())
            var_pred = float(h.var())
            mse = float(np.mean((h - c) ** 2))
            nonconst = var_pred > 1e-9
            per_pair.append({
                "pair": int(i), "parent": pairs[i]["parent"],
                "mutation": pairs[i]["mutation"],
                "n_lig": len(c), "var_true": var_true, "var_pred": var_pred,
                "nonconstant": bool(nonconst),
                "r2": float(1 - mse / var_true) if var_true > 0 else None,
                "mse": mse,
                "spearman": T.spearman(h, c) if nonconst else None,
                "pearson": T.pearson(h, c) if nonconst else None,
                "sign_acc": T.dead_zone_sign_acc(h, c),
                "ols_slope": T.ols_scale(c, h) if nonconst else None,
            })
    # variance decomposition at checkpoint
    with torch.no_grad():
        if is_esm:
            sel = np.asarray(abs_cells if joint else cells)[:64]
            Pa = torch.stack([row_vec(int(i), "var") for i, _ in sel])
            ligs = [l if joint else targets[int(i)]["lig_idx"][l]
                    for i, l in sel]
            Ls = L[ligs]
            ep_ = torch.relu(model.p_enc(Pa))
            el_ = torch.relu(model.l_enc(Ls))
        else:
            rows = abs_rows[:64] if joint else sorted({pairs[i]["var_row"] for i in tr})[:64]
            Pa = Prot[rows]
            Ls = L[torch.arange(min(64, Y.shape[1]), device=device)]
            ep_ = torch.relu(model.p_enc(Pa))
            el_ = torch.relu(model.l_enc(Ls))
        var_dec = {
            "var_b_P": float(model.b_P(ep_).var()),
            "var_b_L": float(model.b_L(el_).var()),
            "var_s": float((model.alpha(ep_) * model.psi(el_)).sum(-1).var()),
            "var_f": float((model.mu + model.b_P(ep_).squeeze(-1)
                            + model.b_L(el_).squeeze(-1)
                            + (model.alpha(ep_) * model.psi(el_)).sum(-1)).var()),
        }
    agg = {
        "r2": float(np.nanmean([r["r2"] for r in per_pair])),
        "spearman": float(np.nanmean([r["spearman"] for r in per_pair
                                      if r["spearman"] is not None])),
        "sign_acc": float(np.nanmean([r["sign_acc"] for r in per_pair
                                      if r["sign_acc"] == r["sign_acc"]])),
        "mse": float(np.mean([r["mse"] for r in per_pair])),
        "n_nonconstant": sum(1 for r in per_pair if r["nonconstant"]),
        "n_rank_evaluable": sum(1 for r in per_pair if r["nonconstant"]),
        "n_parents_covered": len({r["parent"] for r in per_pair
                                  if r["nonconstant"]}),
        "n_test_pairs": len(per_pair),
    }
    collapsed = (agg["n_nonconstant"] < 5 or agg["n_rank_evaluable"] < 5
                 or agg["n_parents_covered"] < 4 or agg["r2"] <= 0.02)
    return {"cell": cell, "seed": SEED, "best_val_mse": best_val,
            "grad_cov": gradcov, "var_dec": var_dec, "agg": agg,
            "collapsed": bool(collapsed), "test_rows": per_pair}


def effect_boot(rows_by_cell, key):
    """Bootstrap pair-mean-R2 differences over test parents (6)."""
    rng = T.stable_rng("stageCIIP2x2", "boot", key, "r2", 20260820)
    def cell_means(cell):
        rows = rows_by_cell[cell]
        byp = {}
        for r in rows:
            byp.setdefault(r["parent"], []).append(r["r2"])
        return byp
    parents = sorted({r["parent"] for rows in rows_by_cell.values() for r in rows})
    def diff(ma, mb, idx):
        sa = [x for pi in idx for x in ma[parents[pi]]]
        sb = [x for pi in idx for x in mb[parents[pi]]]
        return float(np.mean(sa) - np.mean(sb))
    ma = cell_means("esm_joint"); mb = cell_means("klifs_joint")
    mc = cell_means("esm_centered"); md = cell_means("klifs_centered")
    draws = {"rep_main_joint": [], "rep_main_centered": [],
             "obj_main_klifs": [], "obj_main_esm": [], "interaction": []}
    for _ in range(2000):
        idx = rng.integers(len(parents), size=len(parents))
        draws["rep_main_joint"].append(diff(ma, mb, idx))
        draws["rep_main_centered"].append(diff(mc, md, idx))
        draws["obj_main_klifs"].append(diff(md, mb, idx))
        draws["obj_main_esm"].append(diff(mc, ma, idx))
        draws["interaction"].append(diff(mc, ma, idx) - diff(md, mb, idx))
    out = {}
    for k, v in draws.items():
        point = float(np.mean(v))
        lo = float(np.percentile(v, 2.5))
        hi = float(np.percentile(v, 97.5))
        status = ("established" if lo > 0 and abs(point) >= 0.05
                  else "absent" if abs(point) < 0.02 else "ambiguous")
        out[k] = {"point": point, "lo2.5": lo, "hi97.5": hi, "status": status}
    # few-pair safeguard: leave-one-parent-out sign stability
    for k in ("rep_main_joint", "obj_main_klifs"):
        signs = set()
        for pi in parents:
            idx = [p for p in parents if p != pi]
            d = (diff(ma, mb, idx) if k == "rep_main_joint" else diff(md, mb, idx))
            signs.add("+" if d > 0 else "-")
        out[k]["leave_one_parent_out_sign_stable"] = len(signs) == 1
    return out


def main() -> int:
    device = "cuda" if torch.cuda.is_available() else "cpu"
    d1, z1, d2, z2 = load()
    results = {}
    for cell in CELLS:
        print(f"=== {cell} ===", flush=True)
        r = run_cell(cell, d1, z1, d2, z2, device)
        results[cell] = r
        print(json.dumps(r["agg"], indent=1), "collapsed:", r["collapsed"], flush=True)
    rows_by_cell = {c: results[c]["test_rows"] for c in CELLS}
    effects = effect_boot(rows_by_cell, "diag")
    print(json.dumps(effects, indent=1), flush=True)
    out = {
        "schema": "MetaSieve.StageCIIP1A.2x2.Result.v1",
        "preregistration_2x2_sha256": PREREG2X2_SHA,
        "data2x2_sha256": hashlib.sha256((HERE / "DATA2X2.json").read_bytes()).hexdigest(),
        "cells": results, "effects": effects,
        "scope": "oracle-covered subset diagnostic; root-cause attribution only; "
                 "no CIIP-1A PASS verdict possible from this stage",
    }
    path = HERE / "RESULT_2X2_DIAG.json"
    path.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
