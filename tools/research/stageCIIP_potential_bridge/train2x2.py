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


def rng_key(cell, kind, ep=None, mb=None):
    parts = ["stageCIIP2x2", kind, cell, str(SEED)]
    if ep is not None:
        parts += ["epoch", str(ep)]
    if mb is not None:
        parts += ["mb", str(mb)]
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
    step1_grads = None
    abs_epoch_stats = {}
    first_step = True
    for ep in range(EPOCHS):
        rng = rng_key(cell, "order", ep)
        order = np.asarray(cells)[rng.permutation(len(cells))]
        abs_vals = []
        for b0 in range(0, len(order), BATCH):
            batch = order[b0:b0 + BATCH]
            opt.zero_grad()
            hats = {}
            for i in {int(i) for i, _ in batch}:
                hats[i] = pair_hat(i)
            hi = torch.stack([hats[int(i)][j] for i, j in batch])
            ci = torch.stack([c_of[int(i)][j] for i, j in batch])
            loss_ctr = ((hi - ci) ** 2).mean()
            loss_abs = None
            if joint and abs_cells:
                # implementation amendment 2026-08-19: the minibatch index b0
                # enters the SHA-256 keyed stream so every contrast minibatch
                # consumes its OWN absolute cells (previous launches re-sampled
                # the same <=512 cells ~12x per epoch, inflating L_abs weight)
                rng_a = rng_key(cell, "abs", ep, b0)
                sel = np.asarray(abs_cells)[rng_a.permutation(len(abs_cells))[:BATCH]]
                if is_esm:
                    ra = torch.tensor([pairs[int(i)]["var_row"] for i, _ in sel],
                                      device=device)
                    Pa = torch.stack([row_vec(int(i), "var") for i, _ in sel])
                    la = torch.tensor([t for _, t in sel], device=device)
                    f = model(Pa, L[la])
                    y = Y[ra, la]
                else:
                    ra = torch.tensor([r for r, _ in sel], device=device)
                    la = torch.tensor([l for _, l in sel], device=device)
                    fin = torch.isfinite(Y[ra, la])
                    if fin.any():
                        ra = ra[fin]
                        la = la[fin]
                        f = model(Prot[ra], L[la])
                        y = Y[ra, la]
                    else:
                        f, y = None, None
                if f is not None:
                    loss_abs = ((f - y) ** 2).mean()
                    abs_vals.append(float(loss_abs.detach()))
            loss = loss_ctr if loss_abs is None else loss_ctr + LAMBDA_ABS * loss_abs
            if first_step:
                # reporting only: per-loss gradient norms on s-params at step 1
                s_pre = ("alpha.", "psi.")
                loss_ctr.backward(retain_graph=True)
                g_ctr_v = {n: p.grad.detach().clone() for n, p in
                           model.named_parameters() if n.startswith(s_pre)}
                g_ctr_s = float(sum((g ** 2).sum() for g in g_ctr_v.values()).sqrt())
                model.zero_grad()
                g_abs_s = None
                dot = None
                if loss_abs is not None:
                    loss_abs.backward(retain_graph=True)
                    g_abs_s = float(sum((p.grad ** 2).sum() for n, p in
                                        model.named_parameters()
                                        if n.startswith(s_pre)).sqrt())
                    dot = float(sum((p.grad * g_ctr_v[n]).sum() for n, p in
                                    model.named_parameters()
                                    if n.startswith(s_pre)))
                step1_grads = {
                    "g_ctr_s": g_ctr_s, "g_abs_s": g_abs_s,
                    "R_g": (g_abs_s / (g_ctr_s + 1e-12)) if g_abs_s is not None else None,
                    "C_g": (dot / (g_ctr_s * g_abs_s + 1e-12)) if g_abs_s is not None else None,
                }
                model.zero_grad()
                first_step = False
            loss.backward()
            if gradcov is None:
                gradcov = {n: bool(p.grad is not None and float(p.grad.abs().max()) > 0)
                           for n, p in model.named_parameters()}
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
            opt.step()
        if joint and abs_vals and (ep == 0 or ep == EPOCHS - 1):
            abs_epoch_stats[f"epoch{ep}"] = {
                "mean": float(np.mean(abs_vals)), "var": float(np.var(abs_vals)),
                "n_batches": len(abs_vals)}
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
            Ls = L[torch.arange(len(rows), device=device)]
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
    # objective-sampling report (reviewer-required)
    if is_esm:
        n_pool = len(abs_cells) if joint else 0
        n_valid = n_pool
        wt_var = [0, len({pairs[i]["var_row"] for i in tr})]
    else:
        n_pool = len(abs_cells) if joint else 0
        n_valid = 0
        wt_var = [0, 0]
        if joint:
            ypool = z1["Y"][np.ix_(abs_rows, np.arange(z1["Y"].shape[1]))]
            n_valid = int(np.isfinite(ypool).sum())
            wt_rows = {pairs[i]["wt_row"] for i in tr} & set(abs_rows)
            var_rows = {pairs[i]["var_row"] for i in tr} & set(abs_rows)
            wt_var = [len(wt_rows), len(var_rows)]
    sampling_report = {
        "n_abs_cells_pool": n_pool, "n_abs_valid_labels": n_valid,
        "wt_variant_rows": wt_var,
        "note": ("ESM joint: variant-row cells only (frozen rule); KLIFS joint: "
                 "WT+variant rows - objective-sampling confound between the two "
                 "joint cells, documented per review") if joint else
                "centered-only: no L_abs",
    }
    # structural information cap of the representation on covered test pairs
    if is_esm:
        nz = sum(1 for i in te if float(np.linalg.norm(
            z2["esm_var"][i] - z2["esm_wt"][i])) > 1e-9)
    else:
        nz = sum(1 for i in te if float(np.linalg.norm(
            z1["prot"][pairs[i]["var_row"]] - z1["prot"][pairs[i]["wt_row"]])) > 0)
    structural_cap = {"n_test_pairs_nonzero_input": int(nz),
                      "n_test_pairs": len(te)}
    return {"cell": cell, "seed": SEED, "best_val_mse": best_val,
            "grad_cov": gradcov, "var_dec": var_dec, "agg": agg,
            "collapsed": bool(collapsed), "test_rows": per_pair,
            "step1_grads": step1_grads, "abs_epoch_stats": abs_epoch_stats,
            "objective_sampling_report": sampling_report,
            "structural_cap": structural_cap}


def effect_boot(rows_by_cell, key):
    """Effect estimates over test parents with parent-cluster bootstrap.

    Primary (frozen) estimand = OBSERVED pair-mean R2 difference; the
    bootstrap mean is NOT a point estimate (implementation amendment
    2026-08-19). observed_parent_mean_effect and bootstrap CI are
    reported alongside; leave-one-parent-out sign stability is computed
    for ALL five effects.
    """
    rng = T.stable_rng("stageCIIP2x2", "boot", key, "r2", 20260820)
    parents = sorted({r["parent"] for rows in rows_by_cell.values() for r in rows})
    def pair_mean(cell, keep=None):
        rows = rows_by_cell[cell] if keep is None else [
            r for r in rows_by_cell[cell] if r["parent"] in keep]
        return float(np.mean([r["r2"] for r in rows]))
    def parent_mean(cell, keep=None):
        rs = rows_by_cell[cell] if keep is None else [
            r for r in rows_by_cell[cell] if r["parent"] in keep]
        byp = {}
        for r in rs:
            byp.setdefault(r["parent"], []).append(r["r2"])
        return float(np.mean([float(np.mean(v)) for v in byp.values()]))
    def resample_diff(cell_a, cell_b, idx):
        sa = [r["r2"] for pi in idx for r in rows_by_cell[cell_a]
              if r["parent"] == parents[pi]]
        sb = [r["r2"] for pi in idx for r in rows_by_cell[cell_b]
              if r["parent"] == parents[pi]]
        return float(np.mean(sa) - np.mean(sb))
    spec = {
        "rep_main_joint": ("esm_joint", "klifs_joint", None),
        "rep_main_centered": ("esm_centered", "klifs_centered", None),
        "obj_main_klifs": ("klifs_centered", "klifs_joint", None),
        "obj_main_esm": ("esm_centered", "esm_joint", None),
        "interaction": ("esm_centered", "esm_joint",
                        ("klifs_centered", "klifs_joint")),
    }
    draws = {k: [] for k in spec}
    for _ in range(2000):
        idx = rng.integers(len(parents), size=len(parents))
        for k, (a, b, sub) in spec.items():
            d = resample_diff(a, b, idx)
            if sub is not None:
                d -= resample_diff(sub[0], sub[1], idx)
            draws[k].append(d)
    out = {}
    for k, (a, b, sub) in spec.items():
        obs_pair = pair_mean(a) - pair_mean(b)
        obs_parent = parent_mean(a) - parent_mean(b)
        if sub is not None:
            obs_pair -= (pair_mean(sub[0]) - pair_mean(sub[1]))
            obs_parent -= (parent_mean(sub[0]) - parent_mean(sub[1]))
        lo = float(np.percentile(draws[k], 2.5))
        hi = float(np.percentile(draws[k], 97.5))
        status = ("established" if lo > 0 and abs(obs_pair) >= 0.05
                  else "absent" if abs(obs_pair) < 0.02 else "ambiguous")
        signs = set()
        for pi in parents:
            keep = [p for p in parents if p != pi]
            has_a = any(r["parent"] in keep for r in rows_by_cell[a])
            has_b = any(r["parent"] in keep for r in rows_by_cell[b])
            if not (has_a and has_b):
                signs.add("nan")
                continue
            d = pair_mean(a, keep) - pair_mean(b, keep)
            if sub is not None:
                d -= (pair_mean(sub[0], keep) - pair_mean(sub[1], keep))
            signs.add("+" if d > 0 else "-")
        out[k] = {
            "observed_pair_mean_effect": obs_pair,
            "observed_parent_mean_effect": obs_parent,
            "bootstrap_ci": {"lo2.5": lo, "hi97.5": hi, "draws": 2000,
                             "cluster": "parent"},
            "bootstrap_mean_not_a_point_estimate": float(np.mean(draws[k])),
            "status": status,
            "leave_one_parent_out_sign_stable": len(signs) == 1,
        }
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
        "interpretation_bounds": {
            "klifs_structural_cap": "among the 9 covered test pairs only 3 have "
                "nonzero KLIFS input difference, so the 5/9 nonconstant gate is "
                "unreachable for KLIFS cells by construction: KLIFS collapse is "
                "a structural statement, NOT an optimization statement; the "
                "objective main effect is therefore not fairly estimable on "
                "KLIFS and the interaction is confounded by the KLIFS null space",
            "objective_sampling_confound": "KLIFS joint L_abs uses WT+variant "
                "row cells; ESM joint uses variant-row cells only (frozen rule) "
                "- per-cell pool sizes and valid label counts are reported; if "
                "they differ materially the stage is a matched representation "
                "diagnostic, not a pure factorial causal estimate",
            "coverage_bias": "16 uncovered pairs are 4 whole families (ALK, MET, "
                "LRRK2, TEK-Y1108F) with higher target variance; every conclusion "
                "is restricted to the oracle-covered subset",
            "no_universal_objective_claim": "this stage cannot establish a "
                "general main effect of joint vs centered-only beyond the "
                "covered subset",
        },
    }
    path = HERE / "RESULT_2X2_DIAG.json"
    path.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
