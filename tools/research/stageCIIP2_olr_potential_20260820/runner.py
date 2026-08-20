"""OLR-Potential stage runner (prereg a7b17e8a... + ADD-1 aa8d06af...).

Arms (single-variable ladder):
  A0-prior   analytic train-mean profile (no training)
  A1-bilinear mean-pool full-sequence bilinear (CIIP-1A form, deployable)
  A2-router  LCRR (deployable)
  A3-oid     A2 + panel centering enforced in loss (structural; protein-axis
             centering is an exact identity on contrasts, asserted in tests)
  A4-cfoie   A3 + cross-fitted residual target (nuisance m_hat)
  A5-gain    A4 + assay-gain weights
  C-perm     A5 trained on within-pair deranged ligand labels
  C-randprot A5 evaluated with same-parent-sibling / random variant swap
  C-erased   A5 evaluated on X-erased residue states (needs ERASED_ESM.npz)
  C-wrongmut A5 predictions scored against same-parent sibling targets
  C-famprior family-mean profile prior
  C-free     free pairwise head (ceiling, non-deployable)

All training is SGD (AdamW). No closed-form solvers in the deployed path.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import olr as O  # noqa: E402

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
EPOCHS = 900          # AM-5: raised after qualification diagnosed undertraining
PATIENCE = 150        # AM-5: with cosine schedule below
LR = 3e-3             # AM-5
WD = 1e-4
RANK = 8
HEADS = 1
LAMBDA_DISTILL = 0.5


def build_all():
    d1, z1, d2, esm = O.load_stage_data()
    states = O.construct_states(d1, esm)
    recs = O.pair_tensors(d1, d2, states, z1["lig"])
    return d1, z1, d2, states, recs


def pad_states(recs):
    """Pad all construct states to (n_pairs, 2, n_max, d) + masks."""
    n_max = max(max(r["res_w"].shape[0], r["res_v"].shape[0]) for r in recs)
    d = recs[0]["res_w"].shape[1]
    res = torch.zeros(len(recs), 2, n_max, d)
    mask = torch.zeros(len(recs), 2, n_max, dtype=torch.bool)
    for j, r in enumerate(recs):
        res[j, 0, :r["res_w"].shape[0]] = r["res_w"]
        mask[j, 0, :r["res_w"].shape[0]] = True
        res[j, 1, :r["res_v"].shape[0]] = r["res_v"]
        mask[j, 1, :r["res_v"].shape[0]] = True
    return res, mask


def crossfit_nuisance(recs, lig_np, train_js, fold_of, seed, epochs=200):
    """Return per-pair m_hat(L) from parent-excluding fold models."""
    lig = torch.from_numpy(lig_np).float().to(DEVICE)
    fold_ids = sorted(set(fold_of.values()))
    models = {}
    for f in fold_ids:
        torch.manual_seed(seed * 1000 + f)
        m = O.NuisanceMLP().to(DEVICE)
        opt = torch.optim.AdamW(m.parameters(), lr=LR, weight_decay=WD)
        fit_js = [j for j in train_js if fold_of[recs[j]["parent"]] != f]
        y = torch.tensor(np.concatenate([np.asarray(recs[j]["c"], dtype=np.float32) for j in fit_js]), device=DEVICE)
        X = lig[torch.tensor(np.concatenate([np.asarray(recs[j]["lig_idx"]) for j in fit_js]), device=DEVICE)]
        for _ in range(epochs):
            opt.zero_grad()
            loss = ((m(X) - y) ** 2).mean()
            loss.backward()
            opt.step()
        models[f] = m
    mhat = {}
    for j, r in enumerate(recs):
        own = fold_of[r["parent"]]
        use = [f for f in fold_ids if f != own]
        with torch.no_grad():
            L = lig[torch.tensor(r["lig_idx"], device=DEVICE)]
            pred = torch.stack([models[f](L) for f in use]).mean(0)
        mhat[j] = pred.cpu().numpy()
    return mhat


def targets_for(recs, js, arm, mhat):
    """Target arrays per pair: c (A1-A3) or residual r (A4+)."""
    out = {}
    for j in js:
        t = np.asarray(recs[j]["c"], dtype=np.float32)
        if arm in ("A4-cfoie", "A5-gain") and mhat is not None:
            t = t - mhat[j]
        out[j] = t
    return out


def train_arm(arm, recs, lig_np, res_pad, mask_pad, split, seed,
              mhat=None, weights=None, teacher=None, site=None, log=None):
    """Train one arm; return fitted model and per-pair predictions."""
    train_js = [j for j in range(len(recs)) if split[j] == 0]
    val_js = [j for j in range(len(recs)) if split[j] == 1]
    test_js = [j for j in range(len(recs)) if split[j] == 2]
    rng = np.random.default_rng(seed)
    torch.manual_seed(seed)
    router = arm != "A1-bilinear"
    model = O.OLRPotential(rank=RANK, heads=HEADS, router=router).to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=EPOCHS)
    tgt = targets_for(recs, range(len(recs)), "A5-gain" if arm == "C-perm" else arm, mhat)
    ligT = torch.from_numpy(lig_np).float().to(DEVICE)
    perm_c = {}
    if arm == "C-perm":
        for j in train_js:
            perm_c[j] = O.permute_within_pair(tgt[j], rng)  # permute the arm's own target
    w_full = torch.ones(lig_np.shape[0], device=DEVICE) if weights is None else torch.from_numpy(weights).float().to(DEVICE)
    best_val, best_state, patience = float("inf"), None, 0
    for ep in range(EPOCHS):
        opt.zero_grad()
        kv = {}
        for j in range(len(recs)):
            kv[(j, 0)] = model.construct_kv(res_pad[j, 0], mask_pad[j, 0])
            kv[(j, 1)] = model.construct_kv(res_pad[j, 1], mask_pad[j, 1])
        total = 0.0
        for j in train_js:
            r = recs[j]
            L = ligT[torch.tensor(r["lig_idx"], device=DEVICE)]
            sw = model.s_from_kv(kv[(j, 0)], mask_pad[j, 0], L)
            sv = model.s_from_kv(kv[(j, 1)], mask_pad[j, 1], L)
            chat = (sv - sw)
            chat = chat - chat.mean()
            t = torch.from_numpy(perm_c[j] if arm == "C-perm" and j in perm_c else tgt[j]).float().to(DEVICE)
            w = w_full[torch.tensor(r["lig_idx"], device=DEVICE)]
            loss = (w * (t - chat) ** 2).mean()
            if teacher is not None:
                with torch.no_grad():
                    tw = teacher["model"].s_from_kv(teacher["kv"][(j, 0)], mask_pad[j, 0], L)
                    tv = teacher["model"].s_from_kv(teacher["kv"][(j, 1)], mask_pad[j, 1], L)
                    tchat = (tv - tw)
                    tchat = tchat - tchat.mean()
                loss = loss + LAMBDA_DISTILL * ((tchat - chat) ** 2).mean()
            total = total + loss
        (total / len(train_js)).backward()
        opt.step()
        sched.step()
        with torch.no_grad():
            val = 0.0
            for j in val_js:
                r = recs[j]
                L = ligT[torch.tensor(r["lig_idx"], device=DEVICE)]
                sw = model.s_from_kv(kv[(j, 0)], mask_pad[j, 0], L)
                sv = model.s_from_kv(kv[(j, 1)], mask_pad[j, 1], L)
                chat = sv - sw
                chat = chat - chat.mean()
                # ADD-2: selection on the deployable FULL centered target
                t_full = torch.from_numpy(np.asarray(recs[j]["c"], dtype=np.float32)).to(DEVICE)
                if arm in ("A4-cfoie", "A5-gain", "C-perm") and mhat is not None:
                    t_full = t_full - torch.from_numpy(np.asarray(mhat[j], dtype=np.float32)).to(DEVICE)
                w = w_full[torch.tensor(r["lig_idx"], device=DEVICE)]
                val += float((w * (t_full - chat) ** 2).mean())
            val /= max(1, len(val_js))
        if val < best_val - 1e-6:
            best_val, patience = val, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            patience += 1
            if patience >= PATIENCE:
                break
        if log is not None and ep % 50 == 0:
            print(f"  [{arm} seed {seed}] ep {ep} train {float(total)/max(1,len(train_js)):.1f} val {val:.1f}", flush=True)
    if best_state is not None:
        model.load_state_dict(best_state)
    preds = {}
    with torch.no_grad():
        for j in range(len(recs)):
            r = recs[j]
            L = ligT[torch.tensor(r["lig_idx"], device=DEVICE)]
            kw = model.construct_kv(res_pad[j, 0], mask_pad[j, 0])
            km = model.construct_kv(res_pad[j, 1], mask_pad[j, 1])
            sw = model.s_from_kv(kw, mask_pad[j, 0], L)
            sv = model.s_from_kv(km, mask_pad[j, 1], L)
            chat = sv - sw
            chat = chat - chat.mean()
            p = chat.cpu().numpy()
            # ADD-2: deployable prediction composes the ligand-only nuisance
            if arm in ("A4-cfoie", "A5-gain", "C-perm") and mhat is not None:
                p = p + np.asarray(mhat[j], dtype=np.float32)
            preds[j] = p
    return {"model": model, "preds": preds, "best_val": best_val,
            "train_js": train_js, "val_js": val_js, "test_js": test_js}


def eval_predictions(preds_map, recs, js, tgt_map):
    pc, tc, pp, tp = [], [], [], []
    for j in js:
        p = preds_map[j]
        t = tgt_map[j]
        pc.append(p); tc.append(t); pp.append(p); tp.append(t)
    return O.all_metrics(np.concatenate(pc), np.concatenate(tc), pp, tp)


def instrument_qualification(seed=20260821):
    """Planted rank-4 interaction on real features; A5 must recover the
    planted residual above the A0 prior (>= 0.25 Delta R2, 90% parent
    bootstrap excluding 0)."""
    d1, z1, d2, states, recs = build_all()
    res_pad, mask_pad = pad_states(recs)
    res_pad, mask_pad = res_pad.to(DEVICE), mask_pad.to(DEVICE)
    lig = z1["lig"]
    split = np.array([r["split1"] for r in recs], dtype=np.int8)
    train_js = [j for j in range(len(recs)) if split[j] == 0]
    parents = [r["parent"] for r in recs]
    fold_of = O.folds_by_parent(parents)
    rng = np.random.default_rng(seed)
    # planted factors from real features (deterministic)
    U = rng.normal(0, 1, (O.D_RES, 4)).astype(np.float32)
    V = rng.normal(0, 1, (O.D_LIG, 4)).astype(np.float32)
    zL = lig.astype(np.float32) @ V                                   # (183,4)
    hbar = {}
    for j, r in enumerate(recs):
        if r["parent"] not in hbar:
            hbar[r["parent"]] = r["res_w"].mean(0).numpy()
    # parent-DEVIATION field (sequence-predictable, linear in hbar deviations),
    # scaled so the empirical cross-parent cell variance matches the measured
    # between-parent component 134.8 (ADD-2 instrument redesign)
    grand = np.mean([h for h in hbar.values()], axis=0)
    fields = {p: ((h - grand) @ U) @ zL.T for p, h in hbar.items()}
    cross_var = float(np.var(np.stack(list(fields.values())), axis=0).mean())
    fields = {p: f * np.sqrt(134.8 / cross_var) for p, f in fields.items()}
    def scale_to(x, var):
        v = x.var()
        return x * np.sqrt(var / v) if v > 0 else x
    s0 = scale_to((zL[:, 0] + zL[:, 1]) * 0.5, 50.0)                   # shared pattern
    m_spec, eps = {}, {}
    for j, r in enumerate(recs):
        u = rng.normal(0, 1, 4).astype(np.float32)
        m = scale_to(u @ zL.T, 44.85)
        m_spec[j] = m
        eps[j] = rng.normal(0, np.sqrt(44.85), 183).astype(np.float32)
    syn = {}
    for j, r in enumerate(recs):
        prof = (s0 + fields[r["parent"]] + m_spec[j] + eps[j])[r["lig_idx"]]
        syn[j] = (prof - prof.mean()).astype(np.float32)   # centered on pair panel
    for j, r in enumerate(recs):
        r["c"] = syn[j]                                    # inject planted targets
    w, _ = O.gain_weights(d1, z1, recs, train_js)
    mhat = crossfit_nuisance(recs, lig, train_js, fold_of, seed)
    out = {"seed": seed, "planted": {"parent_shared_var": 134.8, "specific_signal_var": 44.85,
                                     "noise_var": 44.85, "shared_pattern_var": 50.0}}
    tgt_res = targets_for(recs, range(len(recs)), "A5-gain", mhat)
    tgt_full = {j: np.asarray(recs[j]["c"], dtype=np.float32) for j in range(len(recs))}
    test_js = [j for j in range(len(recs)) if split[j] == 2]
    # ADD-2: primary evaluation on the FULL centered target (deployable quantity)
    prior = O.ligand_prior(recs, train_js, test_js)
    pc = np.concatenate([prior[j] for j in test_js])
    tc = np.concatenate([tgt_full[j] for j in test_js])
    out["A0_prior"] = O.all_metrics(pc, tc)
    run = train_arm("A5-gain", recs, lig, res_pad, mask_pad, split, seed=11,
                    mhat=mhat, weights=w)
    pc = np.concatenate([run["preds"][j] for j in test_js])
    out["A5"] = O.all_metrics(pc, tc)
    tr_pc = np.concatenate([run["preds"][j] for j in run["train_js"]])
    tr_tc = np.concatenate([tgt_full[j] for j in run["train_js"]])
    out["A5_train"] = {"r2": O.r2_cells(tr_pc, tr_tc), "var_recovery": O.var_recovery(tr_pc, tr_tc)}
    out["A5"]["r2"] = O.r2_cells(pc, tc)
    # paired per-parent delta with parent bootstrap
    deltas = []
    for p in sorted({recs[j]["parent"] for j in test_js}):
        js = [j for j in test_js if recs[j]["parent"] == p]
        r2a = O.r2_cells(np.concatenate([run["preds"][j] for j in js]),
                         np.concatenate([tgt_full[j] for j in js]))
        r2b = O.r2_cells(np.concatenate([prior[j] for j in js]),
                         np.concatenate([tgt_full[j] for j in js]))
        deltas.append(r2a - r2b)
    mean, lo, hi = O.paired_parent_bootstrap(deltas)
    out["delta_r2_bootstrap"] = {"mean": mean, "lo90": lo, "hi90": hi,
                                 "per_parent": [float(x) for x in deltas]}
    # prereg rule: recovery of planted residual Delta R2 >= 0.25 at alpha=0.1
    out["verdict"] = "INSTRUMENT_OK" if (mean >= 0.25 and lo > 0) else "INSTRUMENT_UNDERPOWERED"
    return out


def main_smoke(seed=11):
    """Phase 4 single-seed smoke on split S1 with gates (prereg section 8)."""
    d1, z1, d2, states, recs = build_all()
    res_pad, mask_pad = pad_states(recs)
    res_pad, mask_pad = res_pad.to(DEVICE), mask_pad.to(DEVICE)
    lig = z1["lig"]
    split = np.array([r["split1"] for r in recs], dtype=np.int8)
    train_js = [j for j in range(len(recs)) if split[j] == 0]
    test_js = [j for j in range(len(recs)) if split[j] == 2]
    parents = [r["parent"] for r in recs]
    fold_of = O.folds_by_parent(parents)
    w, _ = O.gain_weights(d1, z1, recs, train_js)
    mhat = crossfit_nuisance(recs, lig, train_js, fold_of, seed)
    results = {}
    # A0 prior
    prior = O.ligand_prior(recs, train_js, test_js)
    tc = np.concatenate([recs[j]["c"] for j in test_js])
    results["A0-prior"] = O.all_metrics(np.concatenate([prior[j] for j in test_js]), tc)
    tgt_c = {j: np.asarray(recs[j]["c"], dtype=np.float32) for j in range(len(recs))}
    tgt_res = targets_for(recs, range(len(recs)), "A5-gain", mhat)
    for arm in ["A1-bilinear", "A2-router", "A3-oid", "A4-cfoie", "A5-gain", "C-perm"]:
        m = None if arm in ("A1-bilinear", "A2-router", "A3-oid") else mhat
        ww = w if arm in ("A5-gain", "C-perm") else None
        run = train_arm(arm, recs, lig, res_pad, mask_pad, split, seed=seed,
                        mhat=m, weights=ww, log=True)
        # primary metric for ladder arms: their own target definition;
        # report both raw-centered and residual R2 for comparability
        raw = eval_predictions(run["preds"], recs, test_js, tgt_c)
        resm = eval_predictions(run["preds"], recs, test_js, tgt_res)
        results[arm] = {"raw_centered": raw, "residual": resm, "best_val": run["best_val"]}
    out = {"seed": seed, "device": DEVICE, "n_train": len(train_js), "n_test": len(test_js),
           "arms": results,
           "gates": {}}
    a5 = results["A5-gain"]["raw_centered"]["r2"]      # ADD-2 primary
    cp = results["C-perm"]["raw_centered"]["r2"]
    a0 = results["A0-prior"]["r2"]
    out["gates"]["perm_destroyed"] = bool(cp <= 0.02 or cp <= 0.10 * max(a5, 1e-9))
    out["gates"]["a0_sanity"] = bool(0.08 <= a0 <= 0.18)
    out["gates"]["nonconstant"] = results["A5-gain"]["raw_centered"].get("n_nonconstant", 0)
    (HERE / "PHASE4_SMOKE.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(json.dumps({"gates": out["gates"], "A0": a0, "A5_res": a5, "C-perm": cp}, indent=1))
    return out



# ---------------------------------------------------------------- phase 5

def erased_recs(recs):
    z = np.load(HERE / "ERASED_ESM.npz", allow_pickle=False)
    out = []
    for r in recs:
        we = torch.from_numpy(z[f"we_{r['pair_idx']}"].astype(np.float32))
        me = torch.from_numpy(z[f"me_{r['pair_idx']}"].astype(np.float32))
        r2 = dict(r)
        r2["res_w"], r2["res_v"] = we, me
        out.append(r2)
    return out


def swap_variants(recs, split, seed):
    """C-randprot: permute variant constructs within same-split pairs
    (parent-preserving where siblings exist, else any same-split pair)."""
    rng = np.random.default_rng(seed + 777)
    varmap = {}
    js = list(range(len(recs)))
    by_parent = {}
    for j in js:
        by_parent.setdefault(recs[j]["parent"], []).append(j)
    for j in js:
        sibs = [k for k in by_parent[recs[j]["parent"]] if k != j and split[k] == split[j]]
        pool = sibs if sibs else [k for k in js if k != j and split[k] == split[j]]
        varmap[j] = pool[rng.integers(0, len(pool))] if pool else j
    out = []
    for j in js:
        r2 = dict(recs[j])
        r2["res_v"] = recs[varmap[j]]["res_v"]
        out.append(r2)
    return out


def wrongmut_eval(run, recs, test_js, tgt):
    """C-wrongmut: score predictions against same-parent sibling targets."""
    by_parent = {}
    for j in range(len(recs)):
        by_parent.setdefault(recs[j]["parent"], []).append(j)
    own, wrong = [], []
    for j in test_js:
        sibs = [k for k in by_parent[recs[j]["parent"]] if k != j]
        if not sibs:
            continue
        own.append(O.r2_cells(run["preds"][j], tgt[j]))
        wrong.append(np.mean([O.r2_cells(run["preds"][j], tgt[k]) for k in sibs]))
    if not own:
        return {"n": 0}
    own, wrong = np.array(own), np.array(wrong)
    return {"n": len(own), "own_mean": float(np.nanmean(own)),
            "sibling_mean": float(np.nanmean(wrong)),
            "advantage": float(np.nanmean(own - wrong)),
            "frac_positive": float(np.nanmean((own - wrong) > 0))}


def free_head_run(recs, lig_np, split, seed, mhat, weights):
    """C-free ceiling: free pairwise head on pooled features."""
    train_js = [j for j in range(len(recs)) if split[j] == 0]
    val_js = [j for j in range(len(recs)) if split[j] == 1]
    test_js = [j for j in range(len(recs)) if split[j] == 2]
    torch.manual_seed(seed)
    model = O.FreePairwise().to(R_DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    ligT = torch.from_numpy(lig_np).float().to(R_DEVICE)
    w_full = torch.ones(lig_np.shape[0], device=R_DEVICE) if weights is None else torch.from_numpy(weights).float().to(R_DEVICE)
    tgt = {j: np.asarray(recs[j]["c"], dtype=np.float32) - (mhat[j] if mhat is not None else 0.0) for j in range(len(recs))}
    pooled = {}
    for j, r in enumerate(recs):
        pooled[(j, 0)] = r["res_w"].mean(0).to(R_DEVICE)
        pooled[(j, 1)] = r["res_v"].mean(0).to(R_DEVICE)

    def pred(j):
        L = ligT[torch.tensor(r_lig(recs, j), device=R_DEVICE)]
        pw = pooled[(j, 0)].unsqueeze(0).expand(L.shape[0], -1)
        pv = pooled[(j, 1)].unsqueeze(0).expand(L.shape[0], -1)
        out = model(pw, pv, L)
        return out - out.mean()

    best_val, best_state, patience = float("inf"), None, 0
    for ep in range(EPOCHS):
        opt.zero_grad()
        loss = 0.0
        for j in train_js:
            t = torch.from_numpy(tgt[j]).float().to(R_DEVICE)
            w = w_full[torch.tensor(r_lig(recs, j), device=R_DEVICE)]
            loss = loss + (w * (t - pred(j)) ** 2).mean()
        (loss / len(train_js)).backward()
        opt.step()
        with torch.no_grad():
            val = sum(float((torch.from_numpy(tgt[j]).float().to(R_DEVICE) - pred(j)).pow(2).mean()) for j in val_js) / max(1, len(val_js))
        if val < best_val - 1e-6:
            best_val, patience, best_state = val, 0, {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            patience += 1
            if patience >= PATIENCE:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    with torch.no_grad():
        preds = {j: pred(j).cpu().numpy() for j in range(len(recs))}
    return {"preds": preds, "best_val": best_val, "test_js": test_js}


def r_lig(recs, j):
    return np.asarray(recs[j]["lig_idx"])


def main_phase5(only=None):
    """Phase 5: seeds x {S1,S2,S3,SPB}, ladder arms + full control matrix.
    only: restrict to one split name (for parallel execution)."""
    d1, z1, d2, states, recs = build_all()
    lig = z1["lig"]
    splits = {"S1": np.array([r["split1"] for r in recs], dtype=np.int8),
              "S2": O.split_s2_s3(d1, d2, 2), "S3": O.split_s2_s3(d1, d2, 3)}
    spb, spb_test_parents = O.split_spb(d1, d2)
    splits["SPB"] = spb
    seeds = [11, 22, 33, 44, 55]
    secondary_seeds = [11, 22, 33]
    out = {"spb_test_parents": spb_test_parents, "runs": {}}
    for sname, split in splits.items():
        if only and sname != only:
            continue
        res_pad, mask_pad = pad_states(recs)
        res_pad, mask_pad = res_pad.to(R_DEVICE), mask_pad.to(R_DEVICE)
        train_js = [j for j in range(len(recs)) if split[j] == 0]
        test_js = [j for j in range(len(recs)) if split[j] == 2]
        parents = [r["parent"] for r in recs]
        fold_of = O.folds_by_parent(parents)
        w, _ = O.gain_weights(d1, z1, recs, train_js)
        mhat = crossfit_nuisance(recs, lig, train_js, fold_of, seed=99)
        tgt_c = {j: np.asarray(recs[j]["c"], dtype=np.float32) for j in range(len(recs))}
        tgt_res = targets_for(recs, range(len(recs)), "A5-gain", mhat)
        # A0 prior and family prior
        prior = O.ligand_prior(recs, train_js, test_js)
        tc = np.concatenate([tgt_c[j] for j in test_js])   # ADD-2 full target
        pc = np.concatenate([prior[j] for j in test_js])
        out["runs"][sname] = {"A0-prior": O.all_metrics(pc, tc)}
        fam = O.family_prior(recs, train_js, test_js, O.parent_family)
        pc = np.concatenate([fam[j] for j in test_js])
        out["runs"][sname]["C-famprior"] = O.all_metrics(pc, tc)
        # erased-trained arm (C-erased: train on erased states)
        er = erased_recs(recs)
        er_pad, er_mask = pad_states(er)
        er_pad, er_mask = er_pad.to(R_DEVICE), er_mask.to(R_DEVICE)
        arm_seeds = seeds if sname in ("SPB", "S1") else secondary_seeds
        if sname in ("S2", "S3"):
            arm_list = ["A1-bilinear", "A5-gain"]
        else:
            arm_list = ["A1-bilinear", "A2-router", "A4-cfoie", "A5-gain"]
        for arm in arm_list:
            agg = {"raw_centered": [], "residual": []}
            for seed in arm_seeds:
                m = None if arm == "A1-bilinear" or arm == "A2-router" else mhat
                ww = w if arm == "A5-gain" else None
                run = train_arm(arm, recs, lig, res_pad, mask_pad, split, seed=seed, mhat=m, weights=ww)
                agg["raw_centered"].append(eval_predictions(run["preds"], recs, test_js, tgt_c))
                agg["residual"].append(eval_predictions(run["preds"], recs, test_js, tgt_res))
            out["runs"][sname][arm] = agg
        # controls at seed 11
        run5 = train_arm("A5-gain", recs, lig, res_pad, mask_pad, split, seed=11, mhat=mhat, weights=w)
        # C-erased-site (informative): variant construct replaced by its
        # X-erased states (site residue masked, context preserved). Full
        # erasure is degenerate for potentials (erased WT == erased MT ->
        # contrast == 0 identically) and is asserted in tests instead.
        z_er = np.load(HERE / "ERASED_ESM.npz", allow_pickle=False)
        preds_site = {}
        with torch.no_grad():
            for j in test_js:
                r = recs[j]
                L = torch.from_numpy(lig[np.asarray(r["lig_idx"])]).float().to(R_DEVICE)
                mw = torch.from_numpy(z_er[f"me_{r['pair_idx']}"].astype(np.float32)).to(R_DEVICE)
                kw = run5["model"].construct_kv(res_pad[j, 0], mask_pad[j, 0])
                ke = run5["model"].construct_kv(mw, torch.ones(mw.shape[0], dtype=torch.bool, device=R_DEVICE))
                cs = run5["model"].s_from_kv(ke, torch.ones(mw.shape[0], dtype=torch.bool, device=R_DEVICE), L) \
                    - run5["model"].s_from_kv(kw, mask_pad[j, 0], L)
                preds_site[j] = (cs - cs.mean()).cpu().numpy()
        out["runs"][sname]["C-erased-site"] = {
            "residual": eval_predictions(preds_site, recs, test_js, tgt_res)}
        recs_swap = swap_variants(recs, split, seed=11)
        sw_pad, sw_mask = pad_states(recs_swap)
        preds_swap = {}
        with torch.no_grad():
            for j in test_js:
                L = torch.from_numpy(lig[np.asarray(recs[j]["lig_idx"])]).float().to(R_DEVICE)
                kw = run5["model"].construct_kv(sw_pad[j, 0], sw_mask[j, 0])
                km = run5["model"].construct_kv(sw_pad[j, 1], sw_mask[j, 1])
                cs = run5["model"].s_from_kv(km, sw_mask[j, 1], L) - run5["model"].s_from_kv(kw, sw_mask[j, 0], L)
                preds_swap[j] = (cs - cs.mean()).cpu().numpy()
        out["runs"][sname]["C-randprot"] = eval_predictions(preds_swap, recs, test_js, tgt_res)
        out["runs"][sname]["C-wrongmut"] = wrongmut_eval(run5, recs, test_js, tgt_res)
        free = free_head_run(recs, lig, split, seed=11, mhat=mhat, weights=w)
        pc = np.concatenate([free["preds"][j] for j in test_js])
        out["runs"][sname]["C-free"] = O.all_metrics(pc, tc)
        # paired per-parent delta A5 - A0 (primary deployment statistic, ADD-2:
        # evaluated on the FULL centered target)
        deltas = []
        for p in sorted({recs[j]["parent"] for j in test_js}):
            js = [j for j in test_js if recs[j]["parent"] == p]
            r2a = O.r2_cells(np.concatenate([run5["preds"][j] for j in js]), np.concatenate([tgt_c[j] for j in js]))
            r2b = O.r2_cells(np.concatenate([prior[j] for j in js]), np.concatenate([tgt_c[j] for j in js]))
            deltas.append(r2a - r2b)
        mean, lo, hi = O.paired_parent_bootstrap(deltas)
        out["runs"][sname]["A5_vs_A0_parent_bootstrap"] = {"mean": mean, "lo90": lo, "hi90": hi,
                                                           "per_parent": [float(x) for x in deltas]}
        print(f"[{sname}] A5 res R2 seed11 {out['runs'][sname]['A5-gain']['residual'][0]['r2']:.4f} "
              f"A0 {out['runs'][sname]['A0-prior']['r2']:.4f} delta_boot ({mean:.3f}, {lo:.3f}, {hi:.3f})", flush=True)
        (HERE / f"PHASE5_{sname}.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
        merged = json.loads((HERE / "PHASE5.json").read_text()) if (HERE / "PHASE5.json").exists() else {"runs": {}}
        merged.setdefault("spb_test_parents", spb_test_parents)
        merged["runs"][sname] = out["runs"][sname]
        (HERE / "PHASE5.json").write_text(json.dumps(merged, indent=1), encoding="utf-8")
    print("phase5 complete")


R_DEVICE = DEVICE

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        main_smoke(int(sys.argv[2]) if len(sys.argv) > 2 else 11)
    elif mode == "qualify":
        out = instrument_qualification()
        (HERE / "QUALIFICATION.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
        print(json.dumps({k: out[k] for k in ("verdict", "delta_r2_bootstrap", "A0_prior", "A5") if k in out}, indent=1)[:1500])
    elif mode == "phase5":
        main_phase5(sys.argv[2] if len(sys.argv) > 2 else None)

