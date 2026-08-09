"""Phase 2B — ligand-conditioned residue residual over a frozen pocket prior.

Registered by research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md
(SHA-256 5e6688f68bf214b6f44c96ef0a2b909eba99da31f1adf17cdde003addb242c96,
committed as b9753db BEFORE this file was written).

Model:
    b_r^P(P) = b + alpha * w_pi(GELU(W_h(LN(h_r))))        frozen B5 branch only
    delta_raw_r(P,L) = sum_k (U h_r)_k (V g(L))_k          K = 8, 10,568 params
    delta = (I - Q_P Q_P^T) delta_raw                      Q_P = onb{1, b^P}
    s_r(P,L) = b_r^P(P) + delta_r(P,L)

b^P is ligand-independent by construction, so it cancels exactly in the
same-protein difference; Q_P is protein-fixed, so it commutes with that
difference and cannot carry ligand information.

No affinity source is opened. ESM2 and B5 are never trained.
"""
from __future__ import annotations

import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

ROOT = Path(r"D:\MetaSieve")
sys.path.insert(0, str(ROOT / "research" / "s7_l2b_r0r"))
from s7_localizer import Localizer  # noqa: E402

PROC = ROOT / "dataset" / "processed" / "s7_l2b_r0r"
ESM = PROC / "esm2_650M"
S5 = PROC / "sealed_preds_b5"
S4 = PROC / "sealed_preds"
OUT = ROOT / "report" / "s7_l2b_r0r"
P2B = PROC / "phase2b"

PREREG_SHA = "5e6688f68bf214b6f44c96ef0a2b909eba99da31f1adf17cdde003addb242c96"
PREREG_COMMIT = "b9753db"

# ---- frozen contract (preregistration sections 3, 8, 13) -------------------
K = 8
D_ESM = 1280
D_ATOM = 41
SEED_PARAM = 20260901
SEED_SAMPLER = 20260902
SEED_BOOT = 20260903
SEED_CTRL = 20260904
SEED_SYNTH = 20260905
EPOCHS = 6
LR = 1e-3
WD = 1e-4
CLIP = 5.0
C_MAX = 2                 # constructs sampled per component per epoch
P_MAX = 8                 # pairs sampled per construct per epoch
BATCH_COMPONENTS = 16
N_BOOT = 10000
ORTHO_TOL = 1e-8
GS_TOL = 1e-10
CANCEL_TOL = 1e-12
N_PARAMS_EXPECTED = K * D_ESM + K * D_ATOM      # 10,568
SYNTH_M = 8
SYNTH_MIN_AP = 0.50
MIN_PARAM_MOVEMENT = 0.05
MIN_ACT_VAR = 1e-8


def sha_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""):
            h.update(c)
    return h.hexdigest()


# ------------------------------------------------------------- frozen B5
def load_b5(device):
    ck = torch.load(S5 / "B5_checkpoint.pt", map_location=device)
    m = Localizer(D_ESM, use_protein=True).to(device)
    m.load_state_dict(ck)
    m.eval()
    for p in m.parameters():
        p.requires_grad_(False)
    return m


@torch.no_grad()
def protein_prior(b5, h_t):
    """b_r^P(P) = b + alpha * w_pi(GELU(W_h(LN(h)))). Ligand-independent by
    construction: the expression contains no atom or ligand term."""
    r = b5.residue_state(h_t)
    pi = b5.w_pi(r).squeeze(-1)
    return (b5.b + b5.alpha * pi).double().cpu().numpy()


# ------------------------------------------------------------- gauge
def nuisance_basis(bP: np.ndarray, tol: float = GS_TOL) -> np.ndarray:
    """Orthonormal basis of {1, b^P} by modified Gram-Schmidt in float64.
    Columns whose residual norm falls below `tol` relative to their original
    norm are dropped, which is also the constant/degenerate-b fallback."""
    L = bP.size
    cols = [np.ones(L, dtype=np.float64), np.asarray(bP, dtype=np.float64)]
    basis = []
    for c in cols:
        n0 = np.linalg.norm(c)
        if n0 <= 0:
            continue
        v = c.copy()
        for q in basis:
            v -= q * float(q @ v)
        n = np.linalg.norm(v)
        if n <= tol * n0:
            continue
        basis.append(v / n)
    return (np.stack(basis, axis=1) if basis
            else np.zeros((L, 0), dtype=np.float64))


def project_np(Q: np.ndarray, d: np.ndarray) -> np.ndarray:
    d = np.asarray(d, dtype=np.float64)
    if Q.shape[1] == 0:
        return d
    return d - Q @ (Q.T @ d)


def ortho_ratio(Q: np.ndarray, d: np.ndarray) -> float:
    if Q.shape[1] == 0:
        return 0.0
    return float(np.linalg.norm(Q.T @ d) / (1e-30 + np.linalg.norm(d)))


# ------------------------------------------------------------- head
class Head(nn.Module):
    """The only trainable object in Phase 2B: U (8x1280) and V (8x41), no bias."""

    def __init__(self, seed: int = SEED_PARAM):
        super().__init__()
        g = torch.Generator().manual_seed(seed)
        self.U = nn.Parameter(torch.empty(K, D_ESM))
        self.V = nn.Parameter(torch.empty(K, D_ATOM))
        for t in (self.U, self.V):
            nn.init.xavier_uniform_(t, gain=1.0)
        # deterministic re-init through an explicit generator
        with torch.no_grad():
            self.U.copy_(_xavier(K, D_ESM, g))
            self.V.copy_(_xavier(K, D_ATOM, g))

    def uh(self, h):                       # (L, D) -> (L, K)
        return h @ self.U.T

    def vg(self, g):                       # (D_ATOM,) -> (K,)
        return self.V @ g


def _xavier(fan_out, fan_in, gen):
    a = float(np.sqrt(6.0 / (fan_in + fan_out)))
    return (torch.rand(fan_out, fan_in, generator=gen) * 2.0 - 1.0) * a


def n_params(head) -> int:
    return sum(p.numel() for p in head.parameters())


def g_of(atom_feats: np.ndarray) -> np.ndarray:
    """g(L) = arithmetic mean over the frozen 41-D atom features. Invariant to
    atom permutation by construction."""
    return np.asarray(atom_feats, dtype=np.float64).mean(axis=0)


# ------------------------------------------------------------- exact tied AP
_HARM = np.concatenate(([0.0], np.cumsum(1.0 / np.arange(1, 8193, dtype=np.float64))))


def _harm(n):
    global _HARM
    if _HARM.size < n + 1:
        _HARM = np.concatenate(([0.0], np.cumsum(1.0 / np.arange(1, n + 1,
                                                                 dtype=np.float64))))
    return _HARM


def ap_exact(scores: np.ndarray, labels: np.ndarray):
    """Exact E[AP] under a uniformly random ordering inside each tied block.
    Identical estimator to Phase 2A amendments 01/02, evaluated in closed form."""
    labels = np.asarray(labels)
    P = float(labels.sum())
    if P == 0 or P == labels.size:
        return None
    N = scores.size
    harm = _harm(N)
    order = np.argsort(-scores, kind="stable")
    s = scores[order]
    y = labels[order].astype(np.float64)
    bounds = np.flatnonzero(np.r_[True, s[1:] != s[:-1], True])
    csum = np.concatenate(([0.0], np.cumsum(y)))
    lo = bounds[:-1].astype(np.int64)
    hi = bounds[1:].astype(np.int64)
    k = csum[hi] - csum[lo]
    m = k > 0
    if not m.any():
        return 0.0
    n = (hi - lo).astype(np.float64)[m]
    a = lo.astype(np.float64)[m]
    b = csum[lo][m]
    k = k[m]
    s1 = harm[hi[m]] - harm[lo[m]]
    grow = np.where(n > 1.0, (k - 1.0) / np.maximum(n - 1.0, 1.0), 0.0)
    return float((((k / n) * ((b + 1.0) * s1 + grow * (n - (a + 1.0) * s1))).sum()) / P)


def chance_ap(labels: np.ndarray):
    """Exact E[AP] of a constant score: one tied block over all residues."""
    return ap_exact(np.zeros(labels.size, dtype=np.float64), labels)


def pair_metrics(dscore: np.ndarray, gain_set, loss_set, L: int):
    """All-residue gain/loss/change AP with per-pair chance levels."""
    gain = np.zeros(L, dtype=np.int8)
    loss = np.zeros(L, dtype=np.int8)
    for r in gain_set:
        gain[r] = 1
    for r in loss_set:
        loss[r] = 1
    change = ((gain + loss) > 0).astype(np.int8)
    if change.sum() == 0:
        return None
    out = {}
    out["ap_gain"] = ap_exact(dscore, gain)
    out["ap_loss"] = ap_exact(-dscore, loss)
    out["ap_change"] = ap_exact(np.abs(dscore), change)
    out["chance_gain"] = chance_ap(gain)
    out["chance_loss"] = chance_ap(loss)
    out["chance_change"] = chance_ap(change)
    g, l = out["ap_gain"], out["ap_loss"]
    cg, cl = out["chance_gain"], out["chance_loss"]
    vals = [v for v in (g, l) if v is not None]
    chs = [c for c, v in ((cg, g), (cl, l)) if v is not None]
    out["ap_bidir"] = float(np.mean(vals)) if vals else None
    out["chance_bidir"] = float(np.mean(chs)) if chs else None
    # secondary, conditional-sign diagnostic on the symmetric-difference set only
    idx = np.flatnonzero(change)
    if idx.size >= 2 and gain[idx].sum() not in (0, idx.size):
        out["ap_symdiff_conditional"] = ap_exact(dscore[idx], gain[idx])
    else:
        out["ap_symdiff_conditional"] = None
    return out


# ------------------------------------------------------- aggregation
def aggregate(pair_values, pair_construct, construct_component):
    """residue -> unordered pair -> construct -> closure component -> macro.
    Invariant to pair duplication because every level is a mean."""
    by_c = defaultdict(list)
    for pid, v in pair_values.items():
        if v is not None:
            by_c[pair_construct[pid]].append(v)
    con = {c: float(np.mean(v)) for c, v in by_c.items() if v}
    by_comp = defaultdict(list)
    for c, v in con.items():
        by_comp[construct_component[c]].append(v)
    comp = {k: float(np.mean(v)) for k, v in by_comp.items()}
    macro = float(np.mean(list(comp.values()))) if comp else float("nan")
    return comp, macro


def component_bootstrap(a_comp, b_comp, n_boot=N_BOOT, seed=SEED_BOOT):
    keys = sorted(set(a_comp) & set(b_comp))
    if not keys:
        return {"delta": float("nan"), "lcb95_one_sided": float("nan"), "units": 0}
    a = np.array([a_comp[k] for k in keys])
    b = np.array([b_comp[k] for k in keys])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(keys), (n_boot, len(keys)))
    d = (a[idx] - b[idx]).mean(1)
    return {"delta": float((a - b).mean()),
            "lcb95_one_sided": float(np.percentile(d, 5)),
            "ci95": [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))],
            "units": len(keys)}


# ------------------------------------------------------- eligible pairs
def build_pairs(records, masks):
    """Unordered, within one exact construct, scaffold-distinct, different
    ligand graph, non-empty symmetric difference. Zero-difference pairs are
    excluded by this frozen rule and counted."""
    by_sk = defaultdict(list)
    for r in records:
        by_sk[r["seq_key"]].append(r)
    pairs, excl = [], {"zero_symmetric_difference": 0, "same_graph": 0,
                       "scaffold_not_distinct": 0}
    for sk, recs in by_sk.items():
        recs = sorted(recs, key=lambda r: r["source_key"])
        for i in range(len(recs)):
            for j in range(i + 1, len(recs)):
                a, b = recs[i], recs[j]
                if a["graph_key"] == b["graph_key"]:
                    excl["same_graph"] += 1
                    continue
                sa, sb = a["scaffold"], b["scaffold"]
                if not (sa and sb and sa != sb):
                    excl["scaffold_not_distinct"] += 1
                    continue
                Ra, Rb = masks[a["source_key"]], masks[b["source_key"]]
                if not (Ra ^ Rb):
                    excl["zero_symmetric_difference"] += 1
                    continue
                pairs.append((sk, a["source_key"], b["source_key"]))
    return pairs, excl


def hierarchical_sample(pairs, construct_component, epoch,
                        c_max=C_MAX, p_max=P_MAX, seed=SEED_SAMPLER):
    """Deterministic hierarchical sampler: every component is visited, at most
    c_max constructs per component, at most p_max pairs per construct."""
    rng = np.random.default_rng(seed + epoch)
    by_sk = defaultdict(list)
    for p in pairs:
        by_sk[p[0]].append(p)
    by_comp = defaultdict(list)
    for sk in by_sk:
        by_comp[construct_component[sk]].append(sk)
    chosen = []
    for comp in sorted(by_comp):
        sks = sorted(by_comp[comp])
        if len(sks) > c_max:
            sks = [sks[i] for i in sorted(rng.choice(len(sks), c_max, replace=False))]
        for sk in sks:
            pl = by_sk[sk]
            if len(pl) > p_max:
                pl = [pl[i] for i in sorted(rng.choice(len(pl), p_max, replace=False))]
            chosen.append((comp, sk, pl))
    return chosen


# ------------------------------------------------------- loss
def pair_loss(ds: torch.Tensor, gain_set, loss_set, L: int, device):
    """Mean over the non-empty groups of BCE against targets 1 / 0 / 0.5.
    Group means make the objective invariant to pocket size and to the
    overwhelming majority of unchanged residues."""
    gain = sorted(gain_set)
    loss = sorted(loss_set)
    if not gain and not loss:
        return None
    ch = np.zeros(L, dtype=bool)
    ch[gain] = True
    ch[loss] = True
    unch = np.flatnonzero(~ch)
    terms = []
    bce = nn.functional.binary_cross_entropy_with_logits
    if gain:
        t = ds[torch.tensor(gain, device=device)]
        terms.append(bce(t, torch.ones_like(t)))
    if loss:
        t = ds[torch.tensor(loss, device=device)]
        terms.append(bce(t, torch.zeros_like(t)))
    if unch.size:
        t = ds[torch.tensor(unch, device=device)]
        terms.append(bce(t, torch.full_like(t, 0.5)))
    return torch.stack(terms).mean()


# ------------------------------------------------------- residue-context control
def context_shuffle(h: np.ndarray, seq: str, seed: int) -> np.ndarray:
    """Shuffle frozen contextual residue states WITHIN the same amino-acid type
    within the same sequence. Composition and per-type state distribution are
    preserved; contextual assignment is destroyed."""
    rng = np.random.default_rng(seed)
    out = h.copy()
    by_aa = defaultdict(list)
    for i, ch in enumerate(seq[:h.shape[0]]):
        by_aa[ch].append(i)
    for aa, idx in sorted(by_aa.items()):
        if len(idx) < 2:
            continue
        perm = rng.permutation(len(idx))
        out[np.asarray(idx)] = h[np.asarray(idx)[perm]]
    return out


def derange(n: int, rng) -> list[int]:
    """A permutation with no fixed point where n >= 2; identity when n < 2."""
    if n < 2:
        return list(range(n))
    for _ in range(200):
        p = rng.permutation(n)
        if not np.any(p == np.arange(n)):
            return [int(x) for x in p]
    p = list(range(1, n)) + [0]        # deterministic cyclic fallback
    return p
