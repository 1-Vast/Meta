"""S7_L2B — exact-residue localizer: arms, training, complete-matrix evaluation.

Registered by PREREG_S7_L2B_UNIFIED.md (sha256 2c333f22..., commit ce186f4,
frozen BEFORE any model existed).

Arms: B0 prevalence, BL ligand-only, B4 non-PLM residue baseline,
      BP wrong protein, BX wrong ligand, BM motif shuffle.
B5 (frozen ESM2) is added only under the registered escalation rule.

Evaluation scores the COMPLETE residue x heavy-atom matrix. Training subsamples
negatives; evaluation never does. Inference units are protein closure
components. No affinity value is read.
"""
from __future__ import annotations

import argparse
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
from s7_dataset import (ATOM_DIM, RES_DIM, atom_features, build,  # noqa: E402
                        make_split, protein_components, residue_features)

OUT = ROOT / "report" / "s7_l2b_r0r"
CACHE = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "preds"
SEED_MODEL = 20260812
SEED_NEG = 20260813
SEED_CTRL = 20260814
SEED_BOOT = 20260811
N_NEG = 6
RANK = 32
D_RES = 128
EPOCHS = 6
LR = 1e-3
WD = 1e-4
N_BOOT = 2000


def sha(b):
    return hashlib.sha256(b).hexdigest()


# --------------------------------------------------------------- metric
def average_precision(scores, labels, res_idx, atom_idx):
    """Deterministic AP. Ties broken by (residue_index, atom_slot), never by
    array order, so AP cannot depend on how the matrix was flattened."""
    order = np.lexsort((atom_idx, res_idx, -scores))
    y = labels[order].astype(np.float64)
    n_pos = y.sum()
    if n_pos == 0:
        return None
    tp = np.cumsum(y)
    prec = tp / np.arange(1, y.size + 1, dtype=np.float64)
    return float((prec * y).sum() / n_pos)


def component_macro(per_complex, comp_of):
    by = defaultdict(list)
    for key, ap in per_complex.items():
        if ap is not None:
            by[comp_of[key]].append(ap)
    comp_mean = {c: float(np.mean(v)) for c, v in by.items()}
    return comp_mean, float(np.mean(list(comp_mean.values()))) if comp_mean else float("nan")


def paired_bootstrap(a_comp, b_comp, n_boot=N_BOOT, seed=SEED_BOOT):
    keys = sorted(set(a_comp) & set(b_comp))
    a = np.array([a_comp[k] for k in keys])
    b = np.array([b_comp[k] for k in keys])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(keys), (n_boot, len(keys)))
    d = (a[idx] - b[idx]).mean(1)
    return {"delta": float((a - b).mean()),
            "lcb95_one_sided": float(np.percentile(d, 5)),
            "ci95": [float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))],
            "units": len(keys)}


# --------------------------------------------------------------- model
class Localizer(nn.Module):
    """r_i = GELU(W_h LN(x_i)); pi=sigmoid(w_pi.r); rho=sigmoid(w_rho.a);
       s = b + alpha*logit(pi) + beta*logit(rho) + (P r).(Q a)"""

    def __init__(self, d_res_in, use_protein=True):
        super().__init__()
        self.use_protein = use_protein
        self.ln = nn.LayerNorm(d_res_in)
        self.W_h = nn.Linear(d_res_in, D_RES)
        self.w_pi = nn.Linear(D_RES, 1)
        self.w_rho = nn.Linear(ATOM_DIM, 1)
        self.P = nn.Linear(D_RES, RANK, bias=False)
        self.Q = nn.Linear(ATOM_DIM, RANK, bias=False)
        self.b = nn.Parameter(torch.zeros(1))
        self.alpha = nn.Parameter(torch.ones(1))
        self.beta = nn.Parameter(torch.ones(1))

    def residue_state(self, x):
        return torch.nn.functional.gelu(self.W_h(self.ln(x)))

    def pair_logits(self, x_res, x_atom, block=()):
        a = x_atom
        rho_logit = self.w_rho(a).squeeze(-1)
        if "atom" in block:
            rho_logit = rho_logit.detach()
        s = self.b + self.beta * rho_logit.unsqueeze(0)
        if self.use_protein:
            r = self.residue_state(x_res)
            if "protein" in block:
                r = r.detach()
            pi_logit = self.w_pi(r).squeeze(-1)
            s = s + self.alpha * pi_logit.unsqueeze(1)
            inter = self.P(r) @ self.Q(a).T
            if "interaction" in block:
                inter = inter.detach()
            s = s + inter
        else:
            s = s + torch.zeros(x_res.shape[0], 1, device=s.device)
        return s


# --------------------------------------------------------------- sampling
def sample_negatives(rec, rng):
    """Exactly N_NEG unique true negatives per positive: residue-local,
    atom-local, uniform, with uniform backfill. Never a positive, never a dup."""
    L, A = rec["n_res"], rec["n_atoms"]
    pos = set(map(tuple, rec["edges"]))
    out = []
    for (ri, aj) in rec["edges"]:
        got = set()
        cats = [("res", 2), ("atom", 2), ("uni", 2)]
        for kind, want in cats:
            tries = 0
            n_added = 0
            while n_added < want and tries < 60:
                tries += 1
                if kind == "res":
                    cand = (ri, int(rng.integers(0, A)))
                elif kind == "atom":
                    cand = (int(rng.integers(0, L)), aj)
                else:
                    cand = (int(rng.integers(0, L)), int(rng.integers(0, A)))
                if cand in pos or cand in got:
                    continue
                got.add(cand)
                n_added += 1
        tries = 0
        while len(got) < N_NEG and tries < 500:
            tries += 1
            cand = (int(rng.integers(0, L)), int(rng.integers(0, A)))
            if cand in pos or cand in got:
                continue
            got.add(cand)
        if len(got) != N_NEG:
            continue
        out.extend(got)
    return out


# --------------------------------------------------------------- features
class FeatureStore:
    def __init__(self):
        self.res = {}
        self.atom = {}

    def residues(self, rec):
        k = rec["seq_key"]
        if k not in self.res:
            self.res[k] = residue_features(rec["uniprot_sequence"])
        return self.res[k]

    def atoms(self, rec, mol_cache):
        k = rec["graph_key"]
        if k not in self.atom:
            self.atom[k] = atom_features(mol_cache[rec["ligand_ccd"]])
        return self.atom[k]


def train_arm(recs, store, mol_cache, use_protein, device, log=None):
    torch.manual_seed(SEED_MODEL)
    model = Localizer(RES_DIM, use_protein=use_protein).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    lossfn = nn.BCEWithLogitsLoss()
    trace = []
    for ep in range(EPOCHS):
        rng = np.random.default_rng(SEED_NEG + ep)
        order = np.random.default_rng(SEED_MODEL + ep).permutation(len(recs))
        tot, nb = 0.0, 0
        for bi in range(0, len(order), 16):
            batch = [recs[i] for i in order[bi:bi + 16]]
            losses = []
            for rec in batch:
                negs = sample_negatives(rec, rng)
                if not negs:
                    continue
                rows = rec["edges"] + negs
                y = torch.zeros(len(rows), device=device)
                y[:len(rec["edges"])] = 1.0
                xr = torch.from_numpy(store.residues(rec)).to(device)
                xa = torch.from_numpy(store.atoms(rec, mol_cache)).to(device)
                ri = torch.tensor([p[0] for p in rows], device=device)
                ai = torch.tensor([p[1] for p in rows], device=device)
                s = model.pair_logits(xr, xa)[ri, ai]
                losses.append(lossfn(s, y))
            if not losses:
                continue
            loss = torch.stack(losses).mean()
            opt.zero_grad(set_to_none=True)
            loss.backward()
            gn = {n: float(p.grad.norm()) for n, p in model.named_parameters()
                  if p.grad is not None}
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            tot += float(loss)
            nb += 1
            if nb % 200 == 0:
                trace.append({"epoch": ep, "batch": nb, "loss": tot / nb,
                              "grad_norms": gn})
                if log:
                    print(f"    ep{ep} b{nb} loss={tot/nb:.5f}", flush=True)
        if log:
            print(f"  epoch {ep} mean loss {tot/max(nb,1):.5f}", flush=True)
    return model, trace


@torch.no_grad()
def evaluate(model, recs, store, mol_cache, device, comp_of,
             res_override=None, atom_override=None, shuffle_residues=False,
             prevalence=False):
    per_complex, rows = {}, []
    rng = np.random.default_rng(SEED_CTRL)
    for rec in recs:
        L, A = rec["n_res"], rec["n_atoms"]
        y = np.zeros((L, A), dtype=np.int8)
        for ri, aj in rec["edges"]:
            y[ri, aj] = 1
        if prevalence:
            s = np.zeros((L, A), dtype=np.float64)
        else:
            xr_np = store.residues(rec)
            if res_override is not None:
                xr_np = res_override(rec)
            if shuffle_residues:
                xr_np = xr_np[rng.permutation(L)]
            xa_np = store.atoms(rec, mol_cache)
            if atom_override is not None:
                xa_np = atom_override(rec)
            xr = torch.from_numpy(np.ascontiguousarray(xr_np)).to(device)
            xa = torch.from_numpy(np.ascontiguousarray(xa_np)).to(device)
            s = model.pair_logits(xr, xa).double().cpu().numpy()
        ri_idx = np.repeat(np.arange(L), A)
        ai_idx = np.tile(np.arange(A), L)
        ap = average_precision(s.ravel(), y.ravel(), ri_idx, ai_idx)
        per_complex[rec["source_key"]] = ap
        rows.append({"key": rec["source_key"], "component": comp_of[rec["source_key"]],
                     "n_res": L, "n_atoms": A, "positives": int(y.sum()), "ap": ap})
    comp_mean, macro = component_macro(per_complex, comp_of)
    return per_complex, comp_mean, macro, rows


def build_controls(held, train, comp_of, store, mol_cache):
    """Frozen, score-blind wrong-partner maps."""
    rng = np.random.default_rng(SEED_CTRL)
    pool = train
    by_len = sorted(pool, key=lambda r: r["n_res"])
    lens = np.array([r["n_res"] for r in by_len])
    prot_map, lig_map = {}, {}
    by_atoms = defaultdict(list)
    for r in pool:
        by_atoms[r["n_atoms"]].append(r)
    for rec in held:
        L = rec["n_res"]
        j = int(np.searchsorted(lens, L))
        cand = [p for p in by_len[j:j + 40]
                if comp_of[p["source_key"]] != comp_of[rec["source_key"]]
                and p["n_res"] >= L]
        prot_map[rec["source_key"]] = cand[0]["source_key"] if cand else None
        opts = [p for p in by_atoms.get(rec["n_atoms"], [])
                if p["graph_key"] != rec["graph_key"]
                and (not p["scaffold"] or p["scaffold"] != rec["scaffold"])]
        lig_map[rec["source_key"]] = (opts[int(rng.integers(0, len(opts)))]["source_key"]
                                      if opts else None)
    return prot_map, lig_map
