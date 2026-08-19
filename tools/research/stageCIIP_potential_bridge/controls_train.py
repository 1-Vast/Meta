"""Stage CIIP-1A control-arm trainer (prereg 39d02166...).

Seven arms on the matched 49-pair covered subset, centered-only
objective, identical budget. PRIMARY question: does the oracle local
ESM nonconstant response depend on the real mutation-centered window,
or on generic local context / annotation shortcuts?
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import torch

HERE = Path(__file__).resolve().parent
SIG = HERE.parent / "stageX_csc_signal"
X0C = SIG / "stageX0c_measurement_qualification_20260818"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(SIG))
sys.path.insert(0, str(X0C))

import potential as POT  # noqa: E402
import train1a as T  # noqa: E402
from x0_i2 import window_mean_esm, ESM_WINDOW_RADIUS  # noqa: E402
from x0_common import stable_rng, normalize_parent_name, normalize_construct_name  # noqa: E402

PREREG_CONTROLS_SHA = "39d02166f69acf235a34d351b649a4cdbf3b828491a0994901bf2378777463f7"
ARMS = ("oracle_local_esm_correct", "family_preserving_shuffle",
        "random_local_window", "ligand_only", "ligand_invariant_shift",
        "random_protein", "free_pairwise")
EPOCHS = 200
BATCH = 512
LR = 1e-3
WD = 1e-4
DEAD_ZONE = 10.0
SEED = 1
BOOT_DRAWS = 2000
BOOT_SEED = 20260820
R = ESM_WINDOW_RADIUS


def rng(arm, kind, ep=None):
    parts = ["stageCIIPcontrols", kind, arm, str(SEED)]
    if ep is not None:
        parts += ["epoch", str(ep)]
    return stable_rng(*parts)


def load():
    d1 = json.loads((HERE / "DATA1A.json").read_text(encoding="utf-8"))
    z1 = np.load(HERE / "DATA1A.npz", allow_pickle=False)
    d2 = json.loads((HERE / "DATA2X2.json").read_text(encoding="utf-8"))
    z2 = np.load(HERE / "DATA2X2.npz", allow_pickle=False)
    esm = np.load(X0C / "q1_esm_cache.npz", allow_pickle=True)
    return d1, z1, d2, z2, esm


def build_windows(arm, d1, d2, z2, esm):
    """(wt, var) windows (n_pairs, 640) for the covered pairs + winpos meta."""
    covered = list(d2["covered_pair_indices"])
    cache = {k: esm[k] for k in esm.files}
    wt_keys = {normalize_parent_name(k[3:]): k for k in cache if k.startswith("wt:")}
    mt_keys = {normalize_construct_name(k[3:]): k for k in cache if k.startswith("mt:")}
    wt = z2["esm_wt"].copy()
    var = z2["esm_var"].copy()
    meta = {"winpos": {}}
    if arm == "oracle_local_esm_correct":
        meta["winpos"] = {int(i): d1["pairs"][i]["pos"] for i in covered}
    elif arm == "family_preserving_shuffle":
        by_parent = {}
        for i in covered:
            by_parent.setdefault(d1["pairs"][i]["parent"], []).append(i)
        rng_rows = rng(arm, "rows")
        for parent in sorted(by_parent):
            idx = sorted(by_parent[parent])
            sh = idx.copy()
            rng_rows.shuffle(sh)
            for a, b in zip(idx, sh):
                var[a] = z2["esm_var"][b]
        meta["winpos"] = {int(i): d1["pairs"][i]["pos"] for i in covered}
    elif arm == "random_protein":
        rng_rows = rng(arm, "rows")
        sh = covered.copy()
        rng_rows.shuffle(sh)
        for a, b in zip(covered, sh):
            var[a] = z2["esm_var"][b]
        meta["winpos"] = {int(i): d1["pairs"][i]["pos"] for i in covered}
    elif arm == "random_local_window":
        rng_pos = rng(arm, "winpos")
        for i in covered:
            p = d1["pairs"][i]
            wk = wt_keys[normalize_parent_name(p["parent"])]
            mk = mt_keys[normalize_construct_name(p["var_label"])]
            Lseq = cache[wk].shape[0] - 1
            lo, hi = 1 + R, Lseq - R
            true_pos = p["pos"]
            cand = [q for q in range(lo, hi + 1) if abs(q - true_pos) > R]
            q = int(rng_pos.choice(cand))
            wt[i] = window_mean_esm(cache[wk], q, R)
            var[i] = window_mean_esm(cache[mk], q, R)
            meta["winpos"][int(i)] = q
    elif arm in ("ligand_only", "ligand_invariant_shift"):
        wt[:] = 0.0
        var[:] = 0.0
    elif arm == "free_pairwise":
        meta["winpos"] = {int(i): d1["pairs"][i]["pos"] for i in covered}
    else:
        raise ValueError(arm)
    return wt, var, meta


def train_arm(arm, d1, z1, d2, z2, esm, device):
    torch.manual_seed(SEED)
    targets = d1["targets"]
    covered = list(d2["covered_pair_indices"])
    sp = np.asarray(d1["split"]["pair_split"])
    tr = [i for i in covered if sp[i] == 0]
    va = [i for i in covered if sp[i] == 1]
    wt, var, meta = build_windows(arm, d1, d2, z2, esm)
    Wt = torch.from_numpy(wt).float().to(device)
    Vt = torch.from_numpy(var).float().to(device)
    L = torch.from_numpy(z1["lig"]).float().to(device)
    c_of = {i: torch.tensor(t["c"], dtype=torch.float32, device=device)
            for i, t in enumerate(targets)}
    if arm == "ligand_invariant_shift":
        return {"arm": arm, "trained": False, "grad_cov": {},
                "best_val_mse": None}, None, wt, var, meta
    model = (POT.FreePairwise(d_p=640).to(device) if arm == "free_pairwise"
             else POT.UnifiedPotential(d_p=640).to(device))
    def pair_hat(i):
        t = targets[i]
        Lm = L[t["lig_idx"]]
        Pw = Wt[i].unsqueeze(0).expand(len(t["lig_idx"]), -1)
        Pv = Vt[i].unsqueeze(0).expand(len(t["lig_idx"]), -1)
        return (model(Pw, Pv, Lm) if arm == "free_pairwise"
                else model.centered_mutation_effect(Pw, Pv, Lm))
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],
                            lr=LR, weight_decay=WD)
    cells = [(i, j) for i in tr for j in range(len(targets[i]["c"]))]
    best_val = None
    best_state = None
    gradcov = None
    for ep in range(EPOCHS):
        rng_ep = rng(arm, "order", ep)
        order = np.asarray(cells)[rng_ep.permutation(len(cells))]
        for b0 in range(0, len(order), BATCH):
            batch = order[b0:b0 + BATCH]
            opt.zero_grad()
            hats = {}
            for i in {int(i) for i, _ in batch}:
                hats[i] = pair_hat(int(i))
            hi = torch.stack([hats[int(i)][j] for i, j in batch])
            ci = torch.stack([c_of[int(i)][j] for i, j in batch])
            loss = ((hi - ci) ** 2).mean()
            loss.backward()
            if gradcov is None:
                gradcov = {n: bool(p.grad is not None and float(p.grad.abs().max()) > 0)
                           for n, p in model.named_parameters()}
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
            opt.step()
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
    return {"arm": arm, "trained": True, "best_val_mse": best_val,
            "grad_cov": gradcov}, model, wt, var, meta


def per_pair_metrics(model, arm, d1, z1, wt, var, device):
    pairs = d1["pairs"]
    targets = d1["targets"]
    d2 = json.loads((HERE / "DATA2X2.json").read_text(encoding="utf-8"))
    covered = list(d2["covered_pair_indices"])
    sp = np.asarray(d1["split"]["pair_split"])
    te = [i for i in covered if sp[i] == 2]
    Wt = torch.from_numpy(wt).float().to(device)
    Vt = torch.from_numpy(var).float().to(device)
    L = torch.from_numpy(z1["lig"]).float().to(device)
    rows = []
    for i in te:
        t = targets[i]
        Lm = L[t["lig_idx"]]
        Pw = Wt[i].unsqueeze(0).expand(len(t["lig_idx"]), -1)
        Pv = Vt[i].unsqueeze(0).expand(len(t["lig_idx"]), -1)
        with torch.no_grad():
            if arm == "ligand_invariant_shift":
                h = np.zeros(len(t["c"]), dtype=np.float32)
            elif model is None:
                h = np.zeros(len(t["c"]), dtype=np.float32)
            elif arm == "free_pairwise":
                h = model(Pw, Pv, Lm).cpu().numpy()
            else:
                h = model.centered_mutation_effect(Pw, Pv, Lm).cpu().numpy()
        c = np.asarray(t["c"])
        var_true = float(c.var())
        var_pred = float(h.var())
        mse = float(np.mean((h - c) ** 2))
        nonconst = var_pred > 1e-9
        rows.append({
            "pair": int(i), "parent": pairs[i]["parent"],
            "mutation": pairs[i]["mutation"], "n_lig": len(c),
            "var_true": var_true, "var_pred": var_pred,
            "scale_ratio": float(np.sqrt(var_pred / var_true)) if var_true > 0 else None,
            "mse": mse, "r2": float(1 - mse / var_true) if var_true > 0 else None,
            "ols_slope": T.ols_scale(c, h) if nonconst else None,
            "sign_acc": T.dead_zone_sign_acc(h, c),
            "spearman": T.spearman(h, c) if nonconst else None,
            "pearson": T.pearson(h, c) if nonconst else None,
            "nonconstant": bool(nonconst),
        })
    agg = {
        "n_nonconstant": sum(1 for r in rows if r["nonconstant"]),
        "n_total": len(rows),
        "n_rank_evaluable": sum(1 for r in rows if r["nonconstant"]),
        "n_parents_covered": len({r["parent"] for r in rows if r["nonconstant"]}),
        "r2": float(np.mean([r["r2"] for r in rows])),
        "spearman": float(np.nanmean([r["spearman"] for r in rows
                                      if r["spearman"] is not None])),
        "sign_acc": float(np.nanmean([r["sign_acc"] for r in rows
                                      if r["sign_acc"] == r["sign_acc"]])),
        "mse": float(np.mean([r["mse"] for r in rows])),
        "slope_median": float(np.nanmedian([r["ols_slope"] for r in rows
                                            if r["ols_slope"] is not None])),
        "scale_ratio_median": float(np.nanmedian([r["scale_ratio"] for r in rows
                                                  if r["scale_ratio"] is not None])),
    }
    return rows, agg


def effect(rows_a, rows_b, key):
    """Observed pair/parent-mean differences + parent-cluster bootstrap CI
    + leave-one-parent-out sign stability (bootstrap mean never a point)."""
    rng_b = stable_rng("stageCIIPcontrols", "boot", key, 20260820)
    parents = sorted({r["parent"] for r in rows_a} | {r["parent"] for r in rows_b})
    def stats(rows, keep=None):
        rs = rows if keep is None else [r for r in rows if r["parent"] in keep]
        pair_mean = float(np.mean([r["r2"] for r in rs]))
        sa = [r.get("sign_acc") for r in rs]
        sign_mean = float(np.mean([s for s in sa if s is not None
                                   and s == s])) if any(s is not None and s == s
                                                        for s in sa) else None
        byp = {}
        for r in rs:
            byp.setdefault(r["parent"], []).append(r["r2"])
        parent_mean = float(np.mean([float(np.mean(v)) for v in byp.values()]))
        return pair_mean, parent_mean, sign_mean
    def resample(rows, idx):
        return [r["r2"] for pi in idx for r in rows if r["parent"] == parents[pi]]
    pair_a, par_a, _ = stats(rows_a)
    pair_b, par_b, _ = stats(rows_b)
    draws = []
    for _ in range(BOOT_DRAWS):
        idx = rng_b.integers(len(parents), size=len(parents))
        sa = resample(rows_a, idx)
        sb = resample(rows_b, idx)
        draws.append(float(np.mean(sa) - np.mean(sb)))
    lo = float(np.percentile(draws, 2.5))
    hi = float(np.percentile(draws, 97.5))
    signs = set()
    for pi in parents:
        keep = [p for p in parents if p != pi]
        if not any(r["parent"] in keep for r in rows_a) or not any(
                r["parent"] in keep for r in rows_b):
            signs.add("nan")
            continue
        pa, _, _ = stats(rows_a, keep)
        pb, _, _ = stats(rows_b, keep)
        signs.add("+" if pa > pb else "-")
    return {
        "observed_pair_mean_effect": pair_a - pair_b,
        "observed_parent_mean_effect": par_a - par_b,
        "bootstrap_ci": {"lo2.5": lo, "hi97.5": hi, "draws": BOOT_DRAWS,
                         "cluster": "parent"},
        "bootstrap_mean_not_a_point_estimate": float(np.mean(draws)),
        "leave_one_parent_out_sign_stable": len(signs) == 1,
    }


def annotation_shortcut_audit(d1, d2, z2, esm):
    """Feature-level: correct-site window delta norm vs random matched
    window delta norm, paired over covered pairs."""
    covered = list(d2["covered_pair_indices"])
    cache = {k: esm[k] for k in esm.files}
    wt_keys = {normalize_parent_name(k[3:]): k for k in cache if k.startswith("wt:")}
    mt_keys = {normalize_construct_name(k[3:]): k for k in cache if k.startswith("mt:")}
    rng_pos = rng("random_local_window", "winpos")
    recs = []
    for i in covered:
        p = d1["pairs"][i]
        wk = wt_keys[normalize_parent_name(p["parent"])]
        mk = mt_keys[normalize_construct_name(p["var_label"])]
        Lseq = cache[wk].shape[0] - 1
        lo, hi = 1 + R, Lseq - R
        cand = [q for q in range(lo, hi + 1) if abs(q - p["pos"]) > R]
        q = int(rng_pos.choice(cand))
        dc = window_mean_esm(cache[mk], p["pos"], R) - window_mean_esm(cache[wk], p["pos"], R)
        dr = window_mean_esm(cache[mk], q, R) - window_mean_esm(cache[wk], q, R)
        recs.append({"pair": int(i), "parent": p["parent"],
                     "correct_delta_norm": float(np.linalg.norm(dc)),
                     "random_delta_norm": float(np.linalg.norm(dr)),
                     "correct_minus_random": float(np.linalg.norm(dc) - np.linalg.norm(dr))})
    diff = [r["correct_minus_random"] for r in recs]
    return {
        "n_pairs": len(recs),
        "per_pair": recs,
        "correct_norm_mean": float(np.mean([r["correct_delta_norm"] for r in recs])),
        "random_norm_mean": float(np.mean([r["random_delta_norm"] for r in recs])),
        "n_correct_exceeds_random": int(sum(1 for d in diff if d > 0)),
        "paired_bootstrap_lo2.5": float(np.percentile(
            [float(np.mean(rng_pos.choice(diff, size=len(diff), replace=True)))
             for _ in range(BOOT_DRAWS)], 2.5)),
        "note": "if correct-site delta norms do not exceed random matched "
                "window deltas, the mutation site carries no measurable "
                "feature-level information",
    }


def main() -> int:
    device = "cuda" if torch.cuda.is_available() else "cpu"
    d1, z1, d2, z2, esm = load()
    results = {}
    for arm in ARMS:
        print(f"=== {arm} ===", flush=True)
        info, model, wt, var, meta = train_arm(arm, d1, z1, d2, z2, esm, device)
        rows, agg = per_pair_metrics(model, arm, d1, z1, wt, var, device)
        results[arm] = {"info": info, "agg": agg, "test_rows": rows,
                        "winpos": meta.get("winpos")}
        print(json.dumps(agg, indent=1), flush=True)
    correct = results["oracle_local_esm_correct"]["test_rows"]
    effects = {}
    for key, arm in (("v_family", "family_preserving_shuffle"),
                     ("v_random_window", "random_local_window"),
                     ("v_random_protein", "random_protein"),
                     ("v_ligand_only", "ligand_only"),
                     ("v_ligand_invariant", "ligand_invariant_shift")):
        effects[key] = effect(correct, results[arm]["test_rows"], key)
    effects["free_v_correct"] = effect(results["free_pairwise"]["test_rows"],
                                       correct, "free_v_correct")
    audit = annotation_shortcut_audit(d1, d2, z2, esm)
    print(json.dumps(effects, indent=1), flush=True)
    print(json.dumps(audit, indent=1)[:400], flush=True)
    out = {
        "schema": "MetaSieve.StageCIIP1A.Controls.Result.v1",
        "preregistration_controls_sha256": PREREG_CONTROLS_SHA,
        "data2x2_sha256": hashlib.sha256(
            (HERE / "DATA2X2.json").read_bytes()).hexdigest(),
        "scope": "oracle-covered subset (49 pairs); annotation-shortcut "
                 "audit; centered-only objective for all arms",
        "arms": results, "effects": effects,
        "annotation_shortcut_audit": audit,
    }
    path = HERE / "CONTROL_RESULT.json"
    path.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
