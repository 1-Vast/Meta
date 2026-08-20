"""CIIP-S1 library: data, features, targets, F9 profiles, metrics, bootstrap.

Frozen by PREREGISTRATION.md (sha 1fb7133b...) + S1_ADDENDUM_THRESHOLDS_20260820.md
(sha beefb620...). Writes nothing; all outputs go through s1run.py.
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
BRIDGE = HERE.parent / "stageCIIP_potential_bridge"
SIG = HERE.parent / "stageX_csc_signal"
X0C = SIG / "stageX0c_measurement_qualification_20260818"
sys.path.insert(0, str(SIG))
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(BRIDGE))

from x0_common import stable_rng, normalize_parent_name, normalize_construct_name  # noqa
from x0_i2 import window_mean_esm, ESM_WINDOW_RADIUS  # noqa

LR, WD, EPOCHS, BATCH, CLIP = 1e-3, 1e-4, 200, 512, 10.0
DEAD = 10.0
RADIUS = ESM_WINDOW_RADIUS
FORM1_ARMS = ("F1", "F5", "F6", "C-perm")
FORM2_ARMS = ("F1f", "F2", "F2w", "F3", "F4", "F7f", "F8f")
ESTIMANDS = ("T0", "T0m", "T1", "T2", "T3")


def rng_for(*parts):
    return stable_rng("S1", *parts)


def load_all():
    d1 = json.loads((BRIDGE / "DATA1A.json").read_text(encoding="utf-8"))
    z1 = np.load(BRIDGE / "DATA1A.npz", allow_pickle=False)
    d2 = json.loads((BRIDGE / "DATA2X2.json").read_text(encoding="utf-8"))
    esm = np.load(X0C / "q1_esm_cache.npz", allow_pickle=True)
    return d1, z1, d2, esm


def _pool(h, seq_len):
    """Mean residue state, residues 1..min(seq_len,1020) (row 0 is <cls>)."""
    n = min(seq_len, 1020)
    return h[1:1 + n].mean(axis=0)


def build_features(erased_path=None):
    """Per covered pair: all frozen feature blocks + legal labels."""
    d1, z1, d2, esm = load_all()
    pairs, targets, sp = d1["pairs"], d1["targets"], np.asarray(d1["split"]["pair_split"])
    Y, prot, lig = z1["Y"], z1["prot"], z1["lig"]
    pt = json.loads((SIG / "X0_PAIR_TABLE.json").read_text(encoding="utf-8"))
    recs = {}
    from x0_i2 import build_pair_records  # lazy import
    _, _, seqs_raw = None, None, None
    import x0_common as xc
    _, _, seqs_raw = xc.load_duongly()
    records = build_pair_records(pt, seqs_raw)
    by_construct = {normalize_construct_name(r["construct"]): r for r in records}
    wk = {normalize_parent_name(k[3:]): k for k in esm.files if k.startswith("wt:")}
    mk = {normalize_construct_name(k[3:]): k for k in esm.files if k.startswith("mt:")}
    erased = np.load(erased_path, allow_pickle=False) if erased_path else None
    feats = {}
    for i in d2["covered_pair_indices"]:
        p = pairs[i]
        rec = by_construct[normalize_construct_name(p["var_label"])]
        L = len(rec["wt_seq"])
        hw = esm[wk[normalize_parent_name(p["wt_label"])]]
        hv = esm[mk[normalize_construct_name(p["var_label"])]]
        f = {
            "parent": p["parent"], "mutation": p["mutation"], "pos": p["pos"],
            "split": int(sp[i]), "seq_len": L,
            "wt_win": np.asarray(window_mean_esm(hw, p["pos"], RADIUS), dtype=np.float32),
            "var_win": np.asarray(window_mean_esm(hv, p["pos"], RADIUS), dtype=np.float32),
            "pool_wt": np.asarray(_pool(hw, L), dtype=np.float32),
            "pool_var": np.asarray(_pool(hv, L), dtype=np.float32),
            "kl_wt": prot[p["wt_row"]].astype(np.float32),
            "kl_var": prot[p["var_row"]].astype(np.float32),
            "lig_idx": np.asarray(targets[i]["lig_idx"], dtype=int),
            "d": np.asarray(targets[i]["d"], dtype=np.float64),
            "c": np.asarray(targets[i]["c"], dtype=np.float64),
        }
        if erased is not None:
            he = erased[f"ewt_{i}"]
            f["er_pool"] = np.asarray(_pool(he, L), dtype=np.float32)
            f["er_win"] = np.asarray(window_mean_esm(he, p["pos"], RADIUS), dtype=np.float32)
        feats[i] = f
    return feats, lig.astype(np.float32), list(d2["covered_pair_indices"])


# ------------------------------------------------------------------ F9 / T2

def f9_profiles(feats, cov):
    """Cross-fitted parent-profile: train pairs leave-self-out; val/test from
    train parents only. Returns (prof_by_pair, defined_by_pair)."""
    train = [i for i in cov if feats[i]["split"] == 0]
    by_par = defaultdict(list)
    for i in train:
        by_par[feats[i]["parent"]].append(i)
    n_lig = 183
    prof, defined = {}, {}
    for i in cov:
        sibs = by_par.get(feats[i]["parent"], [])
        if feats[i]["split"] == 0:
            sibs = [j for j in sibs if j != i]
        if not sibs:
            prof[i], defined[i] = None, False
            continue
        wide = np.full((len(sibs), n_lig), np.nan)
        for k, j in enumerate(sibs):
            wide[k, feats[j]["lig_idx"]] = feats[j]["c"]
        with np.errstate(invalid="ignore"):
            pr = np.nanmean(wide, axis=0)
        prof[i], defined[i] = pr, True
    return prof, defined


def t2_target(feats, prof, i):
    """c minus the pair's cross-fitted profile (NaN where profile undefined)."""
    if prof.get(i) is None:
        return None
    t2 = feats[i]["c"] - prof[i][feats[i]["lig_idx"]]
    return t2


# ------------------------------------------------------------------ metrics

def per_pair_metrics(pred, true, dead=DEAD):
    pred = np.asarray(pred, dtype=float)
    true = np.asarray(true, dtype=float)
    m = np.isfinite(pred) & np.isfinite(true)
    pred, true = pred[m], true[m]
    var_t = float(np.var(true))
    var_p = float(np.var(pred))
    nonconst = var_p > 1e-4
    pc, tc = pred - pred.mean(), true - true.mean()
    cmse = float(np.mean((pc - tc) ** 2)) if len(pred) else float("nan")
    cr2 = float(1 - cmse / var_t) if var_t > 0 else float("nan")
    slope = float((pc * tc).sum() / (tc * tc).sum()) if (tc * tc).sum() > 0 else float("nan")
    dz = np.abs(tc) > dead
    sign = float((np.sign(pc[dz]) == np.sign(tc[dz])).mean()) if dz.sum() else float("nan")
    if len(pred) >= 3 and not np.allclose(pred, pred[0]) and not np.allclose(true, true[0]):
        rho = float(spearmanr(pred, true).statistic)
    else:
        rho = None  # undefined, never 0
    return {"nonconst": bool(nonconst), "var_true": var_t, "var_pred": var_p,
            "scale_ratio": float(np.sqrt(var_p / var_t)) if var_t > 0 else None,
            "cmse": cmse, "cr2": cr2, "ols_slope": slope,
            "sign_acc": sign, "spearman": rho,
            "rank_evaluable": bool(rho is not None), "n": int(len(pred))}


def severity_contrib(pred_scalar, true_scalar):
    """Per-pair standardized-product contribution to the cross-pair Pearson r
    (T0m estimand); deterministic; used for paired contrasts."""
    p = np.asarray(pred_scalar, dtype=float)
    t = np.asarray(true_scalar, dtype=float)
    zp = (p - p.mean()) / (p.std() * np.sqrt(max(len(p), 1)))
    zt = (t - t.mean()) / (t.std() * np.sqrt(max(len(t), 1)))
    return zp * zt


def parent_boot(pairs_idx, delta_by_pair, parents_of, contrast, seed,
                n_draws=2000):
    """Parent-cluster bootstrap of the mean paired delta + LOPO stability."""
    uniq = sorted({parents_of[i] for i in pairs_idx})
    per_parent = {u: float(np.mean([delta_by_pair[i] for i in pairs_idx
                                    if parents_of[i] == u])) for u in uniq}
    arr = np.array([per_parent[u] for u in uniq])
    rng = rng_for("boot", contrast, str(seed))
    draws = np.array([arr[rng.integers(0, len(arr), len(arr))].mean()
                      for _ in range(n_draws)])
    point = float(arr.mean())
    lo, hi = float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))
    lopo = []
    for u in uniq:
        rest = arr[np.array([v != u for v in uniq])]
        lopo.append(float(rest.mean()) if len(rest) else float("nan"))
    lopo_stable = bool(all(np.sign(x) == np.sign(point) for x in lopo if np.isfinite(x)))
    return {"point": point, "lo2.5": lo, "hi97.5": hi,
            "per_parent": per_parent, "lopo_excl_means": lopo,
            "lopo_sign_stable": lopo_stable, "draws": n_draws}


# ------------------------------------------------------------------ controls

def derangement(n, rng):
    idx = np.arange(n)
    while True:
        p = rng.permutation(n)
        if not np.any(p == idx):
            return p


def permuted_targets(feats, pairs_idx, seed, arm="eval"):
    """Within-pair ligand-label permutation (keyed derangement)."""
    out = {}
    for i in pairs_idx:
        rng = rng_for("ligperm", arm, str(seed), str(i))
        out[i] = np.asarray(feats[i]["c"])[derangement(len(feats[i]["c"]), rng)]
    return out


def wrongmut_choice(feats, pairs_idx, seed):
    """Keyed choice of a different same-parent TRAIN pair per eval pair."""
    out = {}
    for i in pairs_idx:
        sibs = [j for j in feats if feats[j]["parent"] == feats[i]["parent"]
                and feats[j]["split"] == 0 and j != i]
        if not sibs:
            out[i] = None
            continue
        rng = rng_for("wrongmut", str(seed), str(i))
        out[i] = sibs[int(rng.integers(len(sibs)))]
    return out


def random_window_pos(feats, i, seed):
    """Keyed non-site position for F6 (excludes [pos-12, pos+12])."""
    L = min(feats[i]["seq_len"], 1020)
    pos = feats[i]["pos"]
    rng = rng_for("winperm", "F6", str(seed), str(i))
    choices = [q for q in range(RADIUS + 1, L - RADIUS) if abs(q - pos) > 2 * RADIUS]
    return choices[int(rng.integers(len(choices)))]


def famshuf_map(feats, cov, seed):
    """Family-preserving shuffle: keyed within-parent permutation of pair
    feature identities (pair i takes the features of another same-parent pair)."""
    by_par = defaultdict(list)
    for i in cov:
        by_par[feats[i]["parent"]].append(i)
    out = {}
    for p, js in sorted(by_par.items()):
        rng = rng_for("famshuf", "F5", str(seed), p)
        if len(js) == 1:
            out[js[0]] = js[0]
            continue
        sh = derangement(len(js), rng)
        out.update({js[k]: js[sh[k]] for k in range(len(js))})
    return out
