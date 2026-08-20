"""OLR-Potential (CIIP-2) implementation. Prereg a7b17e8a... + ADD-1 aa8d06af....

Deployable object: s(P,L) = alpha(P,L)^T beta(L), ligand-conditioned residue
router over frozen ESM-2 residue states. All contrasts are finite differences
of s. No mutation coordinates, no target IDs, no closed-form solvers, no
test-time adaptation anywhere in the deployed path (controls/teachers may use
privileged inputs and are never deployed).
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

HERE = Path(__file__).resolve().parent
BRIDGE = HERE.parent / "stageCIIP_potential_bridge"
SIG = HERE.parent / "stageX_csc_signal"
X0C = SIG / "stageX0c_measurement_qualification_20260818"
sys.path.insert(0, str(SIG))
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(BRIDGE))

from x0_common import stable_rng, normalize_parent_name, normalize_construct_name  # noqa: E402

PREREG_SHA = "a7b17e8a3a6300d1e02bad44233a07d70a91b0b46a7c2a5d4ccd5cdb89489912"
ADD1_SHA = "aa8d06aff536fc107a571d7b722190fffa079910aaf58c85c5c96c7d03d98015"
D_RES = 640
D_LIG = 2048
DEAD_ZONE = 10.0
FAMILIES = {  # Manning groups, frozen
    "TK": {"ABL1", "ALK", "KIT", "MET", "EGFR", "FGFR1", "FGFR2", "FGFR3",
           "FGFR4", "FLT3", "PDGFRA", "RET", "TEK", "SRC", "BTK", "JAK2"},
    "TKL": {"BRAF", "LRRK2"}, "STE": {"MAP2K1"}, "CMGC": {"MAPK14"}, "CAMK": {"CHEK2"},
}
GATEKEEPER = {("ABL1", "T315I"), ("KIT", "T670I"), ("EGFR", "T790M")}
HOTSPOT = {("ABL1", "Y253H"), ("ABL1", "E255K"), ("ABL1", "Q252H"),
           ("KIT", "V560D"), ("EGFR", "L858R"), ("RET", "M918T"),
           ("RET", "C634R"), ("RET", "C634W"), ("FGFR3", "R248C"),
           ("MET", "D1228N"), ("MET", "D1228H"), ("MET", "Y1230C"),
           ("ALK", "F1174L"), ("LRRK2", "G2019S")}


def sha256_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def parent_family(parent: str) -> str:
    for fam, members in FAMILIES.items():
        if parent in members:
            return fam
    return parent


def mutation_class(parent: str, mutation: str) -> str:
    if (parent, mutation) in GATEKEEPER:
        return "gatekeeper"
    if (parent, mutation) in HOTSPOT:
        return "annotated-hotspot"
    return "other"


# ---------------------------------------------------------------- model


class OLRPotential(nn.Module):
    """s(P,L) = alpha(P,L)^T beta(L).

    router=True  : A2+ ligand-conditioned residue routing (LCRR)
    router=False : A1 mean-pooled bilinear (CIIP-1A form on full sequence)
    site_channel : privileged teacher ONLY (never in the deployed path).
    """

    def __init__(self, d_res=D_RES, d_lig=D_LIG, hid=64, rank=8, heads=1,
                 router=True, site_channel=False):
        super().__init__()
        self.router = router
        self.rank = rank
        self.heads = heads
        self.site_channel = site_channel
        self.d_in = d_res + (1 if site_channel else 0)
        if router:
            self.key = nn.Linear(self.d_in, hid * heads, bias=False)
            self.val = nn.Linear(self.d_in, rank * heads, bias=False)
            self.query = nn.Linear(d_lig, hid * heads)
            self.pool = nn.Linear(self.d_in, rank, bias=False)  # AM-4 skip
        else:
            self.alpha = nn.Linear(self.d_in, rank)
            self.psi = nn.Linear(d_lig, rank)
        self.beta = nn.Linear(d_lig, rank)
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.zeros_(m.bias)

    def construct_kv(self, res, mask, site=None):
        """Precompute per-construct keys/values (ligand-independent part).

        res (n, d_in) float, mask (n,) bool, site (n,) float or None.
        Returns dict with keys (n, H*hid) and values (n, H*rank).
        """
        x = res
        if self.site_channel:
            if site is None:
                raise ValueError("site_channel model requires site input")
            x = torch.cat([res, site.unsqueeze(-1)], dim=-1)
        x = x * mask.unsqueeze(-1)
        if self.router:
            pooled = x.sum(0) / mask.sum().clamp(min=1.0)
            return {"k": self.key(x), "v": self.val(x), "pool": self.pool(pooled)}
        return {"h": self.alpha(x.sum(0) / mask.sum().clamp(min=1.0))}

    def s_from_kv(self, kv, mask, lig):
        """Scalar scores for lig batch (B, d_lig) against one construct.

        Mean over heads; softmax over unmasked residues.
        """
        if not self.router:
            a = kv["h"]                                    # (r,)
            b = self.psi(lig)                              # (B, r)
            return (b * a.unsqueeze(0)).sum(-1)
        k = kv["k"]                                        # (n, H*hid)
        v = kv["v"]                                        # (n, H*r)
        q = self.query(lig)                                # (B, H*hid)
        H = self.heads
        hid = k.shape[1] // H
        r = self.rank
        scores = (q.view(-1, H, 1, hid) * k.view(1, -1, H, hid)).sum(-1)  # (B, n, H)
        neg = torch.finfo(scores.dtype).min
        scores = scores.masked_fill(~mask.unsqueeze(0).unsqueeze(-1), neg)
        a = torch.softmax(scores / np.sqrt(hid), dim=1)    # (B, n, H)
        vals = v.view(-1, H, r)                            # (n, H, r)
        routed = torch.einsum("bnh,nhr->bhr", a, vals).mean(1)  # (B, r)
        # AM-4: mean-pool skip + ligand-conditioned routed deviation
        alpha = kv["pool"].unsqueeze(0) + routed           # (B, r)
        beta = self.beta(lig)                              # (B, r)
        return (alpha * beta).sum(-1)

    def s(self, res, mask, lig, site=None):
        return self.s_from_kv(self.construct_kv(res, mask, site), mask, lig)


class FreePairwise(nn.Module):
    """C-free diagnostic head: f(P_w, P_v, L) -> centered-contrast logits.
    NOT deployable (violates potential integrability); ceiling control."""

    def __init__(self, d_res=D_RES, d_lig=D_LIG, hid=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(d_res * 2 + d_lig, hid), nn.ReLU(),
            nn.Linear(hid, hid), nn.ReLU(), nn.Linear(hid, 1))

    def forward(self, pooled_w, pooled_v, lig):
        return self.net(torch.cat([pooled_w, pooled_v, lig], -1)).squeeze(-1)


class NuisanceMLP(nn.Module):
    """Ligand-only nuisance m_hat(L) for cross-fitted residual targets."""

    def __init__(self, d_lig=D_LIG, hid=64):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(d_lig, hid), nn.ReLU(),
                                 nn.Linear(hid, 1))

    def forward(self, lig):
        return self.net(lig).squeeze(-1)


# ---------------------------------------------------------------- data


def load_stage_data():
    d1 = json.loads((BRIDGE / "DATA1A.json").read_text(encoding="utf-8"))
    z1 = np.load(BRIDGE / "DATA1A.npz", allow_pickle=False)
    d2 = json.loads((BRIDGE / "DATA2X2.json").read_text(encoding="utf-8"))
    esm = np.load(X0C / "q1_esm_cache.npz", allow_pickle=True)
    return d1, z1, d2, esm


def construct_states(d1, esm):
    """Residue-state tensors keyed by row index for all covered constructs."""
    rows = d1["rows"]
    covered = d1["covered"] if "covered" in d1 else None
    cache = {k: esm[k] for k in esm.files}
    wt_keys = {normalize_parent_name(k[3:]): k for k in cache if k.startswith("wt:")}
    mt_keys = {normalize_construct_name(k[3:]): k for k in cache if k.startswith("mt:")}
    out = {}
    for i, row in enumerate(rows):
        lab = str(row).strip()
        if "(" in lab:
            key = mt_keys.get(normalize_construct_name(lab))
        else:
            key = wt_keys.get(normalize_parent_name(lab))
        if key is not None:
            out[i] = torch.from_numpy(np.asarray(cache[key], dtype=np.float32))
    return out


def pair_tensors(d1, d2, states, lig):
    """Per covered pair: residue states, target, ligand indices, meta."""
    cov = list(d2["covered_pair_indices"])
    pairs, targets = d1["pairs"], d1["targets"]
    out = []
    for j, i in enumerate(cov):
        p, t = pairs[i], targets[i]
        sw = states[p["wt_row"]]
        sv = states[p["var_row"]]
        out.append({
            "pair_idx": int(i), "j": j,
            "parent": p["parent"], "mutation": p["mutation"],
            "pos": p["pos"], "split1": int(d1["split"]["pair_split"][i]),
            "res_w": sw, "res_v": sv,
            "lig_idx": np.asarray(t["lig_idx"], dtype=np.int64),
            "c": np.asarray(t["c"], dtype=np.float32),
        })
    return out


def covered_rows_raw(d1, z1):
    Y = z1["Y"]
    rows = d1["rows"]
    wt_rows = [i for i, r in enumerate(rows) if "(" not in str(r)]
    return Y, wt_rows, rows


def gain_weights(d1, z1, pair_records, train_js):
    """Assay-gain w_L from TRAIN-parent WT rows only (frozen rule)."""
    Y, wt_rows, rows = covered_rows_raw(d1, z1)
    train_parents = {pair_records[j]["parent"] for j in train_js}
    sel = [i for i in wt_rows if str(rows[i]).strip() in train_parents]
    with np.errstate(invalid="ignore"):
        mu = np.nanmean(Y[sel], axis=0)          # (183,)
    prod = mu * (100.0 - mu)
    med = np.nanmedian(prod)
    prod = np.where(np.isfinite(prod), prod, med)
    w = np.clip(prod, 1e-6, None)
    w = w / w.mean()
    # choose scale s so that clip(s*w, 0.25, 4) has mean exactly 1
    # (monotone in s; 40 bisection steps; preprocessing only, not a model solve)
    lo_s, hi_s = 1e-3, 1e3
    for _ in range(40):
        mid = 0.5 * (lo_s + hi_s)
        m = np.clip(mid * w, 0.25, 4.0).mean()
        if m > 1.0:
            hi_s = mid
        else:
            lo_s = mid
    w = np.clip(0.5 * (lo_s + hi_s) * w, 0.25, 4.0)
    w = (w / w.mean()).astype(np.float32)
    return w, sel


# ---------------------------------------------------------------- splits


def split_s2_s3(d1, d2, which):
    cov = list(d2["covered_pair_indices"])
    parents = [d1["pairs"][i]["parent"] for i in cov]
    rng = stable_rng("stageCIIP2", "split", str(which))
    order = sorted(range(len(cov)), key=lambda j: (parents[j], cov[j]))
    assign = {}
    by_parent = {}
    for j in order:
        by_parent.setdefault(parents[j], []).append(j)
    for parent, js in sorted(by_parent.items()):
        sh = js.copy()
        rng.shuffle(sh)
        n = len(sh)
        n_test = max(1, int(round(n * 0.2)))
        n_val = max(1, int(round(n * 0.2))) if n >= 3 else 0
        assign.update(dict(zip(sh[:n_test], [2] * n_test)))
        if n_val:
            assign.update(dict(zip(sh[n_test:n_test + n_val], [1] * n_val)))
    for j in range(len(cov)):
        assign.setdefault(j, 0)
    return np.array([assign[j] for j in range(len(cov))], dtype=np.int8)


def split_spb(d1, d2):
    """ADD-1 rule: ascending-pair-count greedy over covered multi-pair
    parents until >= 5 test parents and >= 10 test pairs; val = every 7th
    train-side pair in stable (parent, mutation) order."""
    cov = list(d2["covered_pair_indices"])
    parents = [d1["pairs"][i]["parent"] for i in cov]
    counts = {}
    for p in parents:
        counts[p] = counts.get(p, 0) + 1
    multi = sorted(((c, p) for p, c in counts.items() if c >= 2))
    test_parents, held = set(), 0
    for c, p in multi:
        if len(test_parents) >= 5 and held >= 10:
            break
        test_parents.add(p)
        held += c
    assign = {}
    train_side = []
    for j in range(len(cov)):
        if parents[j] in test_parents:
            assign[j] = 2
        else:
            train_side.append(j)
    order = sorted(train_side, key=lambda j: (parents[j], d1["pairs"][cov[j]]["mutation"]))
    for k, j in enumerate(order):
        assign[j] = 1 if (k % 7 == 3) else 0
    return np.array([assign[j] for j in range(len(cov))], dtype=np.int8), sorted(test_parents)


def folds_by_parent(parents, k=3):
    """Parent-grouped folds (T6: no parent spans folds)."""
    uniq = sorted(set(parents))
    rng = stable_rng("stageCIIP2", "folds", *uniq)
    order = uniq.copy()
    rng.shuffle(order)
    fold_of = {p: i % k for i, p in enumerate(order)}
    return fold_of


# ---------------------------------------------------------------- metrics


def _finite(x, y):
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    m = np.isfinite(x) & np.isfinite(y)
    return x[m], y[m]


def r2_cells(pred, true):
    x, y = _finite(pred, true)
    if len(x) < 2:
        return float("nan")
    ss_res = float(((y - x) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def spearman_cells(pred, true):
    from scipy.stats import spearmanr
    x, y = _finite(pred, true)
    if len(x) < 3 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return float("nan")
    r = spearmanr(x, y).statistic
    return float(r)


def pearson_cells(pred, true):
    x, y = _finite(pred, true)
    if len(x) < 3 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def sign_acc_deadzone(pred, true, dz=DEAD_ZONE):
    x, y = _finite(pred, true)
    m = np.abs(y) > dz
    if m.sum() == 0:
        return float("nan")
    return float((np.sign(x[m]) == np.sign(y[m])).mean())


def slope_median(pred, true):
    """Per-pair truth-on-prediction regression slope (diagnostic only;
    uses polyfit in the METRICS module, never in the deployed model)."""
    if len(pred) != len(true):
        raise ValueError("per-pair lists required")
    slopes = []
    for p, t in zip(pred, true):
        x, y = _finite(p, t)
        if len(x) < 3 or np.allclose(x, x[0]):
            continue
        var = float(((x - x.mean()) ** 2).sum())
        slopes.append(float(((x - x.mean()) * (y - y.mean())).sum() / var) if var > 0 else float("nan"))
    return float(np.nanmedian(slopes)) if slopes else float("nan")


def var_recovery(pred, true):
    x, y = _finite(pred, true)
    if len(x) < 2 or np.var(y) == 0:
        return float("nan")
    return float(np.var(x) / np.var(y))


def nonconstant_count(preds_per_pair, thresh=1e-4):
    return int(sum(1 for p in preds_per_pair if np.nanvar(p) > thresh))


def all_metrics(pred_cells, true_cells, preds_per_pair=None, trues_per_pair=None):
    out = {
        "r2": r2_cells(pred_cells, true_cells),
        "spearman": spearman_cells(pred_cells, true_cells),
        "pearson": pearson_cells(pred_cells, true_cells),
        "sign_acc": sign_acc_deadzone(pred_cells, true_cells),
        "var_recovery": var_recovery(pred_cells, true_cells),
        "mse": float(np.nanmean((np.asarray(pred_cells) - np.asarray(true_cells)) ** 2)),
        "n_cells": int(_finite(pred_cells, true_cells)[0].size),
    }
    if preds_per_pair is not None:
        out["slope_median"] = slope_median(preds_per_pair, trues_per_pair)
        out["n_nonconstant"] = nonconstant_count(preds_per_pair)
        out["n_total"] = len(preds_per_pair)
        per = [r2_cells(p, t) for p, t in zip(preds_per_pair, trues_per_pair)]
        per = [v for v in per if v == v]
        out["r2_pair_median"] = float(np.median(per)) if per else float("nan")
    return out


def paired_parent_bootstrap(delta_by_parent, n=2000, seed=20260821, stat=np.mean):
    """Bootstrap over parents of paired per-parent deltas -> (mean, lo10, hi90)."""
    vals = [v for v in delta_by_parent if np.isfinite(v)]
    if len(vals) < 2:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    arr = np.asarray(vals)
    draws = np.array([stat(arr[rng.integers(0, len(arr), len(arr))]) for _ in range(n)])
    return float(arr.mean()), float(np.percentile(draws, 5)), float(np.percentile(draws, 95))


# ---------------------------------------------------------------- controls


def permute_within_pair(c, rng):
    """Derangement of ligand labels within a pair (C-perm)."""
    c = np.asarray(c)
    n = len(c)
    idx = np.arange(n)
    while True:
        p = rng.permutation(n)
        if not np.any(p == idx):
            break
    return c[p]


def ligand_prior(pred_records, train_js, eval_js, n_lig=183):
    """A0-prior: train-mean centered profile m(L) -> per-eval-pair profile."""
    wide = np.full((len(pred_records), n_lig), np.nan)
    items = pred_records.items() if isinstance(pred_records, dict) else enumerate(pred_records)
    for j, rec in items:
        wide[j, np.asarray(rec["lig_idx"], dtype=int)] = np.asarray(rec["c"])
    with np.errstate(invalid="ignore"):
        prof = np.nanmean(wide[list(train_js)], axis=0)
    get = (lambda j: pred_records[j]) if isinstance(pred_records, dict) else (lambda j: pred_records[j])
    return {j: prof[np.asarray(get(j)["lig_idx"], dtype=int)] for j in eval_js}


def family_prior(pred_records, train_js, eval_js, fam_of, n_lig=183):
    """C-famprior: family-mean profile (family-level LOSO by parent)."""
    wide = np.full((len(pred_records), n_lig), np.nan)
    items = pred_records.items() if isinstance(pred_records, dict) else enumerate(pred_records)
    for j, rec in items:
        wide[j, np.asarray(rec["lig_idx"], dtype=int)] = np.asarray(rec["c"])
    fams = {}
    for j in train_js:
        f = fam_of(pred_records[j]["parent"])
        fams.setdefault(f, []).append(j)
    with np.errstate(invalid="ignore"):
        fam_prof = {f: np.nanmean(wide[js], axis=0) for f, js in fams.items()}
    out = {}
    for j in eval_js:
        f = fam_of(pred_records[j]["parent"])
        prof = fam_prof.get(f, np.zeros(n_lig))
        out[j] = prof[np.asarray(pred_records[j]["lig_idx"], dtype=int)]
    return out
