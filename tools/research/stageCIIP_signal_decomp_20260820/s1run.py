"""CIIP-S1 runner: arm matrix, forms, controls, RESULT.json (prereg frozen).

Usage:
  python s1run.py smoke      -> 2 pairs, 5 epochs (structure gate)
  python s1run.py seed 1     -> full arm matrix, seed 1
  python s1run.py seeds 1,2,3-> multi-seed (only after single-seed gate)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

import s1lib as S
from x0_i2 import window_mean_esm

torch.set_num_threads(8)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(S.BRIDGE))
from potential import UnifiedPotential  # noqa: E402  (frozen Form-1)

DEVICE = "cpu"
TORCH_SEED_BASE = 910000


# ------------------------------------------------------------- features

def form1_pair_feats(arm, i, feats, esm, aux):
    """(Pw, Pv) 640-d blocks per arm. aux: famshuf map / random pos."""
    if arm in ("F1", "C-perm"):
        return feats[i]["wt_win"], feats[i]["var_win"]
    if arm == "F5":
        j = aux["famshuf"][i]
        return feats[j]["wt_win"], feats[j]["var_win"]
    if arm == "F6":
        q = aux["randpos"][i]
        hw = aux["esm_wt_state"][i]
        hv = aux["esm_var_state"][i]
        return window_mean_esm(hw, q, S.RADIUS), window_mean_esm(hv, q, S.RADIUS)
    raise ValueError(arm)


def form2_protein_feat(arm, i, feats, aux):
    if arm == "F1f":
        return np.concatenate([feats[i]["wt_win"],
                               feats[i]["var_win"] - feats[i]["wt_win"]])
    if arm == "F2":
        return feats[i]["er_pool"]
    if arm == "F2w":
        return feats[i]["er_win"]
    if arm == "F3":
        return np.concatenate([feats[i]["pool_wt"], feats[i]["pool_var"]])
    if arm == "F4":
        return np.concatenate([feats[i]["kl_wt"], feats[i]["kl_var"]])
    if arm == "F7f":
        return np.zeros(640)
    if arm == "F8f":
        return np.concatenate([feats[i]["pool_wt"], feats[i]["pool_var"]])
    raise ValueError(arm)


def form2_ligand_on(arm):
    return arm != "F8f"


# ------------------------------------------------------------- models

class ProbeMLP(nn.Module):
    def __init__(self, d_in, hid=64):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(d_in, hid), nn.GELU(), nn.Linear(hid, 1))

    def forward(self, x):
        return self.net(x).squeeze(-1)


def make_aux(feats, cov, seed, esm):
    famshuf = S.famshuf_map(feats, cov, seed)
    randpos = {i: S.random_window_pos(feats, i, seed) for i in cov}
    d1, z1, d2, _ = S.load_all()
    import x0_common as xc
    from x0_i2 import build_pair_records, ESM_MAX_LEN
    pt = json.loads((S.SIG / "X0_PAIR_TABLE.json").read_text(encoding="utf-8"))
    _, _, seqs_raw = xc.load_duongly()
    records = build_pair_records(pt, seqs_raw)
    from x0_common import normalize_construct_name
    byc = {normalize_construct_name(r["construct"]): r for r in records}
    wk = {xc.normalize_parent_name(k[3:]): k for k in esm.files if k.startswith("wt:")}
    mk = {xc.normalize_construct_name(k[3:]): k for k in esm.files if k.startswith("mt:")}
    esm_wt_state, esm_var_state = {}, {}
    for i in cov:
        p = d1["pairs"][i]
        esm_wt_state[i] = esm[wk[xc.normalize_parent_name(p["wt_label"])]].astype(np.float64)
        esm_var_state[i] = esm[mk[xc.normalize_construct_name(p["var_label"])]].astype(np.float64)
    return {"famshuf": famshuf, "randpos": randpos,
            "esm_wt_state": esm_wt_state, "esm_var_state": esm_var_state}


# ------------------------------------------------------------- targets

def target_vec(feats, prof, i, estimand, ligperm=None):
    if estimand == "T1":
        return feats[i]["c"]
    if estimand == "T0":
        return feats[i]["d"]
    if estimand == "T2":
        return S.t2_target(feats, prof, i)
    if estimand == "T3":
        return feats[i]["c"]
    raise ValueError(estimand)


# ------------------------------------------------------------- Form-1 train

def _tag_seed(tag):
    """Deterministic integer seed for a tag (no Python hash(); keyed sha256)."""
    import hashlib
    return int.from_bytes(hashlib.sha256(("S1.tag." + str(tag)).encode()).digest()[:4], "big")


def rank_loss(h, t, sel):
    """Pairwise logistic rank loss on the batch cells `sel` of ONE pair.
    Margins between adjacent cells in target-sorted order; indices map
    through sel back into the full per-pair prediction vector h (compliance
    fix 2026-08-20: previously ordv indexed h directly -- wrong ligands)."""
    h_sel = h[sel]
    t_sel = t[sel].detach().cpu().numpy()
    ordv = np.argsort(t_sel)
    a = torch.from_numpy(ordv[:-1].copy()).long()
    b = torch.from_numpy(ordv[1:].copy()).long()
    margin = h_sel[b] - h_sel[a]
    return torch.nn.functional.softplus(-margin).mean()


def train_form1(feats, lig, cov, prof, arm, estimand, seed, aux,
                n_pairs=None, epochs=None, tag=""):
    epochs = epochs or S.EPOCHS
    torch.manual_seed(TORCH_SEED_BASE + seed * 100 + _tag_seed(tag) % 97)
    train = [i for i in cov if feats[i]["split"] == 0][:n_pairs]
    val = [i for i in cov if feats[i]["split"] == 1][:n_pairs]
    model = UnifiedPotential(d_p=640).to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=S.LR, weight_decay=S.WD)
    Lmat = torch.from_numpy(lig).float().to(DEVICE)
    perm_c = {}
    if arm == "C-perm":
        for i in train:
            rng = S.rng_for("ligperm", "C-perm", str(seed), str(i))
            perm_c[i] = feats[i]["c"][S.derangement(len(feats[i]["c"]), rng)]
    cells = []
    for i in train:
        t = target_vec(feats, prof, i, estimand)
        if t is None:
            continue
        if arm == "C-perm":
            t = perm_c[i]
        cells += [(i, k) for k in range(len(t)) if np.isfinite(t[k])]
    best_val, best_state = float("inf"), None

    def pair_hat(i):
        Pw, Pv = form1_pair_feats(arm if arm != "C-perm" else "F1", i, feats, None, aux)
        n = len(feats[i]["lig_idx"])
        Pwt = torch.from_numpy(np.asarray(Pw)).float().unsqueeze(0).expand(n, -1).to(DEVICE)
        Pvt = torch.from_numpy(np.asarray(Pv)).float().unsqueeze(0).expand(n, -1).to(DEVICE)
        return model.centered_mutation_effect(Pwt, Pvt, Lmat[torch.tensor(feats[i]["lig_idx"])])

    for ep in range(epochs):
        rng = S.rng_for("order", tag or arm, str(seed), str(ep))
        order = np.asarray(cells)[rng.permutation(len(cells))] if len(cells) else []
        for b0 in range(0, len(order), S.BATCH):
            batch = order[b0:b0 + S.BATCH]
            opt.zero_grad()
            loss = 0.0
            for i in {int(x) for x, _ in batch}:
                h = pair_hat(i)
                t = torch.from_numpy(
                    perm_c.get(i) if arm == "C-perm" and i in perm_c
                    else target_vec(feats, prof, i, estimand)).float().to(DEVICE)
                sel = torch.tensor([k for (jj, k) in batch if jj == i])
                if estimand == "T3":
                    loss = loss + rank_loss(h, t, sel)
                else:
                    loss = loss + ((h[sel] - t[sel]) ** 2).mean()
            (loss / max(1, len({int(x) for x, _ in batch}))).backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), S.CLIP)
            opt.step()
        if val:
            model.eval()
            with torch.no_grad():
                v = 0.0
                for i in val:
                    t = target_vec(feats, prof, i, estimand)
                    if t is None:
                        continue
                    h = pair_hat(i).cpu().numpy()
                    m = np.isfinite(t) & np.isfinite(h)
                    v += float(np.mean((h[m] - t[m]) ** 2))
            model.train()
            if v < best_val:
                best_val, best_state = v, {k: p.detach().clone() for k, p in model.state_dict().items()}
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    preds = {}
    with torch.no_grad():
        for i in cov:
            preds[i] = pair_hat(i).cpu().numpy()
    return {"model": model, "preds": preds, "best_val": best_val}


# ------------------------------------------------------------- Form-2 train

def train_form2(feats, lig, cov, prof, arm, estimand, seed,
                n_pairs=None, epochs=None, tag=""):
    epochs = epochs or S.EPOCHS
    torch.manual_seed(TORCH_SEED_BASE + 10000 + seed * 100 + _tag_seed(tag) % 97)
    train = [i for i in cov if feats[i]["split"] == 0][:n_pairs]
    val = [i for i in cov if feats[i]["split"] == 1][:n_pairs]
    pfeat = {i: form2_protein_feat(arm, i, feats, None) for i in cov}
    d_in = len(pfeat[cov[0]]) + (2048 if form2_ligand_on(arm) else 0)
    model = ProbeMLP(d_in).to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=S.LR, weight_decay=S.WD)
    ligT = torch.from_numpy(lig).float().to(DEVICE)

    def row(i, k):
        pf = torch.from_numpy(pfeat[i]).float()
        if not form2_ligand_on(arm):
            return pf
        li = feats[i]["lig_idx"][k]
        return torch.cat([pf, ligT[li]])

    def row_pair(i):
        # T0m row: protein features + panel-mean ECFP of the pair.
        pf = torch.from_numpy(pfeat[i]).float()
        if not form2_ligand_on(arm):
            return pf
        return torch.cat([pf, torch.from_numpy(meanlig[i]).float()])

    def target_of(i):
        if estimand == "T0m":
            return None  # per-pair scalar handled by t0m_target
        return target_vec(feats, prof, i, estimand)

    t0m_mode = estimand == "T0m"
    if t0m_mode:
        meanlig = {i: lig[feats[i]["lig_idx"]].mean(axis=0) for i in cov}
        t0m_target = {i: float(np.nanmean(feats[i]["d"])) for i in cov}

    cells = []
    if t0m_mode:
        cells = [(i, 0) for i in train]
    else:
        for i in train:
            t = target_of(i)
            if t is None:
                continue
            cells += [(i, k) for k in range(len(t)) if np.isfinite(t[k])]
    best_val, best_state = float("inf"), None
    for ep in range(epochs):
        rng = S.rng_for("order", tag or arm, str(seed), str(ep))
        order = np.asarray(cells)[rng.permutation(len(cells))] if len(cells) else []
        for b0 in range(0, len(order), S.BATCH):
            batch = order[b0:b0 + S.BATCH]
            if len(batch) == 0:
                continue
            if t0m_mode:
                X = torch.stack([row_pair(int(i)) for i, _ in batch])
                y = torch.tensor([t0m_target[int(i)] for i, _ in batch]).float()
                opt.zero_grad()
                h = model(X)
                loss = ((h - y) ** 2).mean()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), S.CLIP)
                opt.step()
                continue
            X = torch.stack([row(int(i), int(k)) for i, k in batch])
            y = torch.stack([torch.tensor(target_of(int(i))[int(k)]).float() for i, k in batch])
            opt.zero_grad()
            h = model(X)
            if estimand == "T3":
                # compliance fix 2026-08-20: pairwise logistic rank loss on
                # batch cells grouped by pair (was MSE on c == T1 duplicate)
                loss = 0.0
                in_batch = sorted({int(i) for i, _ in batch})
                for pi in in_batch:
                    pos = [q for q, (i, k) in enumerate(batch) if int(i) == pi]
                    tv = np.array([float(target_of(pi)[int(batch[q][1])]) for q in pos])
                    ordv = np.argsort(tv)
                    a = [pos[q] for q in ordv[:-1]]
                    b = [pos[q] for q in ordv[1:]]
                    # b-positions hold the LARGER target: margin must be
                    # h[b] - h[a] (ADD-3 sign fix; was flipped in ADD-2)
                    margin = h[b] - h[a]
                    loss = loss + torch.nn.functional.softplus(-margin).mean()
                loss = loss / max(1, len(in_batch))
            else:
                loss = ((h - y) ** 2).mean()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), S.CLIP)
            opt.step()
        if val:
            model.eval()
            with torch.no_grad():
                v = 0.0
                for i in val:
                    if t0m_mode:
                        h1 = model(row_pair(i).unsqueeze(0)).item()
                        v += (h1 - t0m_target[i]) ** 2
                        continue
                    t = target_of(i)
                    if t is None:
                        continue
                    X = torch.stack([row(i, k) for k in range(len(t))])
                    h = model(X).cpu().numpy()
                    m = np.isfinite(t) & np.isfinite(h)
                    if m.sum():
                        v += float(np.mean((h[m] - t[m]) ** 2))
            model.train()
            if v < best_val:
                best_val, best_state = v, {k: p.detach().clone() for k, p in model.state_dict().items()}
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    preds = {}
    with torch.no_grad():
        if t0m_mode:
            for i in cov:
                preds[i] = float(model(row_pair(i).unsqueeze(0)).item())
        else:
            for i in cov:
                X = torch.stack([row(i, k) for k in range(len(feats[i]["lig_idx"]))])
                preds[i] = model(X).cpu().numpy()
    return {"model": model, "preds": preds, "best_val": best_val}


# ------------------------------------------------------------- eval

def eval_arm(feats, cov, prof, preds, estimand, ligperm=None):
    """Per-pair metrics on the 9 test pairs + aggregates."""
    test = [i for i in cov if feats[i]["split"] == 2]
    out = {"per_pair": {}, "estimand": estimand}
    for i in test:
        if estimand == "T2" and prof.get(i) is None:
            out["per_pair"][str(i)] = {"undefined": True}
            continue
        true = (ligperm or {}).get(i) if ligperm is not None else target_vec(
            feats, prof, i, estimand)
        if true is None:
            true = target_vec(feats, prof, i, estimand)
        out["per_pair"][str(i)] = S.per_pair_metrics(preds[i], true)
    vals = [v for v in out["per_pair"].values()
            if not v.get("undefined") and np.isfinite(v.get("cr2", np.nan))]
    out["agg"] = {
        "n": len(vals),
        "mean_cr2": float(np.mean([v["cr2"] for v in vals])) if vals else None,
        "n_nonconstant": sum(1 for v in vals if v["nonconst"]),
        "n_rank_evaluable": sum(1 for v in vals if v["rank_evaluable"]),
        "median_spearman": (float(np.median([v["spearman"] for v in vals
                                             if v["spearman"] is not None]))
                            if any(v["spearman"] is not None for v in vals) else None),
    }
    return out


def eval_severity(feats, cov, preds):
    """T0m: per-pair predicted severity = mean_l of the T0 prediction;
    agreement with true severity via standardized-product contributions."""
    test = [i for i in cov if feats[i]["split"] == 2]
    ps = np.array([float(np.mean(preds[i])) for i in test])
    ts = np.array([float(np.mean(feats[i]["d"])) for i in test])
    contrib = S.severity_contrib(ps, ts)
    pearson = float(np.corrcoef(ps, ts)[0, 1]) if ps.std() > 0 and ts.std() > 0 else None
    from scipy.stats import spearmanr as _sp
    try:
        sp = float(_sp(ps, ts).statistic)
    except Exception:
        sp = None
    return {"per_pair_contrib": {str(i): float(c) for i, c in zip(test, contrib)},
            "pearson": pearson, "spearman": sp,
            "pred_severity": {str(i): float(x) for i, x in zip(test, ps)},
            "true_severity": {str(i): float(x) for i, x in zip(test, ts)}}


def contrast(feats, cov, evalA, evalB, estimand, name, seed):
    """Paired per-pair cr2 contrast A - B with parent-cluster bootstrap."""
    test = [i for i in cov if feats[i]["split"] == 2]
    delta = {}
    for i in test:
        a = evalA["per_pair"].get(str(i), {})
        b = evalB["per_pair"].get(str(i), {})
        if a.get("undefined") or b.get("undefined"):
            continue
        av, bv = a.get("cr2"), b.get("cr2")
        if av is None or bv is None or not (np.isfinite(av) and np.isfinite(bv)):
            continue
        delta[i] = float(av - bv)
    if len(delta) < 3:
        return {"name": name, "n": len(delta), "note": "insufficient pairs"}
    parents = {i: feats[i]["parent"] for i in test}
    boot = S.parent_boot(list(delta.keys()), delta, parents, name, seed)
    boot.update({"name": name, "estimand": estimand, "n": len(delta)})
    return boot


def severity_contrast(feats, cov, sevA, sevB, name, seed, n_draws=2000):
    """A-proposition contrast (compliance fix 2026-08-20): cross-pair
    SPEARMAN over the 9 test pairs, paired F1f-vs-F2 difference, parent-
    cluster bootstrap of the difference + LOPO sign stability (prereg B.8)."""
    from scipy.stats import spearmanr
    test = [i for i in cov if feats[i]["split"] == 2]
    psA = np.array([sevA["pred_severity"][str(i)] for i in test], dtype=float)
    psB = np.array([sevB["pred_severity"][str(i)] for i in test], dtype=float)
    ts = np.array([sevB["true_severity"][str(i)] for i in test], dtype=float)
    parents = [feats[i]["parent"] for i in test]
    uniq = sorted(set(parents))

    def spear(x, y):
        msk = np.isfinite(x) & np.isfinite(y)
        if msk.sum() < 3 or np.allclose(x[msk], x[msk][0]) or np.allclose(y[msk], y[msk][0]):
            return float("nan")
        return float(spearmanr(x[msk], y[msk]).statistic)

    sA, sB = spear(psA, ts), spear(psB, ts)
    point = sA - sB
    rng = S.rng_for("boot", name, str(seed))
    draws = []
    for _ in range(n_draws):
        sel_u = [uniq[q] for q in rng.integers(0, len(uniq), len(uniq))]
        idx = []
        for u in sel_u:
            idx += [p for p, par in enumerate(parents) if par == u]
        d = spear(psA[idx], ts[idx]) - spear(psB[idx], ts[idx])
        if np.isfinite(d):
            draws.append(d)
    if draws:
        lo, hi = float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))
    else:
        lo = hi = float("nan")
    lopo = []
    for u in uniq:
        idx = [p for p, par in enumerate(parents) if par != u]
        lopo.append(spear(psA[idx], ts[idx]) - spear(psB[idx], ts[idx]))
    fin = [x for x in lopo if np.isfinite(x)]
    lopo_stable = bool(fin and all(np.sign(x) == np.sign(point) for x in fin)) if np.isfinite(point) else False
    return {"name": name, "estimand": "T0m", "n": len(test),
            "point": float(point), "lo2.5": lo, "hi97.5": hi,
            "spearman_A": sA, "spearman_B": sB,
            "lopo_excl": lopo, "lopo_sign_stable": lopo_stable,
            "draws": len(draws)}


def run_seed(seed, smoke=False, erased="ERASED_ESM_S1.npz"):
    t0 = time.time()
    feats, lig, cov = S.build_features(HERE / erased)
    _, _, _, esm = S.load_all()
    prof, defined = S.f9_profiles(feats, cov)
    aux = make_aux(feats, cov, seed, esm)
    np_ = 2 if smoke else None
    ep_ = 5 if smoke else None
    runs, evals, severities = {}, {}, {}
    for arm, est in [("F1", "T1"), ("F1", "T3"), ("F5", "T1"), ("F6", "T1"),
                     ("C-perm", "T1")]:
        tag = f"{arm}-{est}-s{seed}"
        t_start = time.time()
        runs[tag] = train_form1(feats, lig, cov, prof, arm, est, seed, aux,
                                n_pairs=np_, epochs=ep_, tag=tag)
        evals[tag] = eval_arm(feats, cov, prof, runs[tag]["preds"], est)
        print(f"[{tag}] {time.time()-t_start:.0f}s", flush=True)
    for arm in ["F1f", "F2", "F2w", "F3", "F4", "F7f"]:
        for est in ["T0", "T0m", "T1", "T2", "T3"]:
            tag = f"{arm}-{est}-s{seed}"
            t_start = time.time()
            runs[tag] = train_form2(feats, lig, cov, prof, arm, est, seed,
                                    n_pairs=np_, epochs=ep_, tag=tag)
            if est == "T0m":
                severities[tag] = eval_severity(feats, cov, runs[tag]["preds"])
                evals[tag] = {"estimand": "T0m", "per_pair": severities[tag],
                              "agg": {"n": 9,
                                      "spearman": severities[tag]["spearman"],
                                      "pearson": severities[tag]["pearson"]}}
            else:
                evals[tag] = eval_arm(feats, cov, prof, runs[tag]["preds"], est)
                if est == "T0":
                    severities[tag] = eval_severity(feats, cov, runs[tag]["preds"])
            print(f"[{tag}] {time.time()-t_start:.0f}s", flush=True)
    for est in ["T0", "T0m"]:
        tag = f"F8f-{est}-s{seed}"
        t_start = time.time()
        runs[tag] = train_form2(feats, lig, cov, prof, "F8f", est, seed,
                                n_pairs=np_, epochs=ep_, tag=tag)
        if est == "T0m":
            severities[tag] = eval_severity(feats, cov, runs[tag]["preds"])
            evals[tag] = {"estimand": "T0m", "per_pair": severities[tag],
                          "agg": {"n": 9, "spearman": severities[tag]["spearman"],
                                  "pearson": severities[tag]["pearson"]}}
        else:
            evals[tag] = eval_arm(feats, cov, prof, runs[tag]["preds"], "T0")
            severities[tag] = eval_severity(feats, cov, runs[tag]["preds"])
        print(f"[{tag}] {time.time()-t_start:.0f}s", flush=True)
    f9pred = {}
    for i in cov:
        f9pred[i] = (prof[i][feats[i]["lig_idx"]] if prof.get(i) is not None
                     else np.zeros(len(feats[i]["c"])))
    evals[f"F9-T1-s{seed}"] = eval_arm(feats, cov, prof, f9pred, "T1")
    test = [i for i in cov if feats[i]["split"] == 2]
    lperm = S.permuted_targets(feats, test, seed, arm="eval")
    perm_evals = {}
    for tag in [f"F1-T1-s{seed}", f"F1f-T1-s{seed}", f"F3-T1-s{seed}", f"F4-T1-s{seed}"]:
        perm_evals[tag] = eval_arm(feats, cov, prof, runs[tag]["preds"], "T1",
                                   ligperm=lperm)
    wrongmut = {}
    wm_choice = S.wrongmut_choice(feats, test, seed)
    import torch as _t
    for tag, arm in [(f"F1-T1-s{seed}", "F1"), (f"F1f-T1-s{seed}", "F1f"),
                     (f"F3-T1-s{seed}", "F3"), (f"F4-T1-s{seed}", "F4")]:
        preds = {}
        for i in test:
            j = wm_choice[i]
            f2 = dict(feats[i])
            if j is not None:
                if arm in ("F1", "F1f"):
                    f2["var_win"] = feats[j]["var_win"]
                elif arm == "F3":
                    f2["pool_var"] = feats[j]["pool_var"]
                elif arm == "F4":
                    f2["kl_var"] = feats[j]["kl_var"]
            if arm == "F1":
                with _t.no_grad():
                    Pw = _t.from_numpy(np.asarray(f2["wt_win"])).float().unsqueeze(0).expand(len(f2["lig_idx"]), -1)
                    Pv = _t.from_numpy(np.asarray(f2["var_win"])).float().unsqueeze(0).expand(len(f2["lig_idx"]), -1)
                    Lm = _t.from_numpy(lig).float()[_t.tensor(f2["lig_idx"])]
                    preds[i] = runs[tag]["model"].centered_mutation_effect(Pw, Pv, Lm).cpu().numpy()
            else:
                pf = (np.concatenate([f2["wt_win"], f2["var_win"] - f2["wt_win"]])
                      if arm == "F1f"
                      else np.concatenate([f2["pool_wt"], f2["pool_var"]]) if arm == "F3"
                      else np.concatenate([f2["kl_wt"], f2["kl_var"]]))
                X = np.concatenate([np.tile(pf, (len(f2["lig_idx"]), 1)),
                                    lig[f2["lig_idx"]]], axis=1)
                with _t.no_grad():
                    preds[i] = runs[tag]["model"](_t.from_numpy(X).float()).cpu().numpy()
        wrongmut[tag] = eval_arm(feats, cov, prof, preds, "T1")
    contrasts = []
    contrasts.append(contrast(feats, cov, evals[f"F2-T1-s{seed}"],
                              evals[f"F7f-T1-s{seed}"], "T1", "B_F2_vs_F7f", seed))
    contrasts.append(contrast(feats, cov, evals[f"F1-T1-s{seed}"],
                              evals[f"F9-T1-s{seed}"], "T1", "Ctotal_F1_vs_F9", seed))
    contrasts.append(contrast(feats, cov, evals[f"F1f-T2-s{seed}"],
                              evals[f"F2-T2-s{seed}"], "T2", "Csharp_F1f_vs_F2", seed))
    contrasts.append(contrast(feats, cov, evals[f"F3-T1-s{seed}"],
                              evals[f"F9-T1-s{seed}"], "T1", "Deploy_F3_vs_F9", seed))
    contrasts.append(contrast(feats, cov, evals[f"F4-T1-s{seed}"],
                              evals[f"F9-T1-s{seed}"], "T1", "Deploy_F4_vs_F9", seed))
    contrasts.append(contrast(feats, cov, evals[f"F2w-T1-s{seed}"],
                              evals[f"F7f-T1-s{seed}"], "T1", "B_secondary_F2w_vs_F7f", seed))
    sev = severity_contrast(feats, cov, severities[f"F1f-T0m-s{seed}"],
                            severities[f"F2-T0m-s{seed}"], "A_F1f_vs_F2_T0m", seed)
    sev_secondary = severity_contrast(
        feats, cov, severities[f"F1f-T0-s{seed}"],
        severities[f"F2-T0-s{seed}"], "A_secondary_T0derived", seed)
    return {"seed": seed,
            "runs": {k: {"best_val": v["best_val"]} for k, v in runs.items()},
            "evals": evals, "severities": severities, "f9": evals[f"F9-T1-s{seed}"],
            "perm_evals": perm_evals, "wrongmut": wrongmut, "contrasts": contrasts,
            "severity_contrast": sev, "severity_contrast_secondary": sev_secondary,
            "elapsed_s": time.time() - t0}


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    if mode == "smoke":
        out = run_seed(1, smoke=True)
        (HERE / "SMOKE_RESULT.json").write_text(json.dumps(out, indent=1, default=str))
        print("smoke done", round(out["elapsed_s"], 1))
    elif mode == "seed":
        seed = int(sys.argv[2])
        out = run_seed(seed)
        (HERE / f"SEED{seed}_RESULT.json").write_text(json.dumps(out, indent=1, default=str))
        print("seed", seed, "done", round(out["elapsed_s"], 1))
    elif mode == "seeds":
        seeds = [int(x) for x in sys.argv[2].split(",")]
        outs = {str(s): run_seed(s) for s in seeds}
        (HERE / "RESULT.json").write_text(json.dumps(outs, indent=1, default=str))
        print("seeds", seeds, "done")


if __name__ == "__main__":
    main()
