"""Canonical loaders for the XP1 crossed affinity panels.

Frozen by `PREREG_XP1.md`.  Nothing here fits a model or computes a contrast.
"""
from __future__ import annotations

import hashlib
import json
import os

import numpy as np
import pandas as pd

RAW = r"D:\MetaSieve\dataset\raw\crossed_panels"
KIN = os.path.join(RAW, "kinase_panels")
PDSP = os.path.join(RAW, "pdsp_kidb")
ANN = os.path.join(RAW, "protein_annotation")

METZ_FLOOR = 4.0
KLAEGER_FLOOR = 5.0

FROZEN_SHA = {
    "metz_matrix.csv": "abe1e3c580478775a352ec5ee78ca565d4c863f0e3e642fdb21d956d8f9d4375",
    "klaeger_matrix.csv": "cdf66c7d4e7c1e3a35aeb6995abbfdaf15be80f3e07715524b2bb4449d871010",
    "KiDatabase.csv": "45c9a18ac30f1fad350d1dde186bc1f226c5a75d474ca50f50713852a5637ac6",
}


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_releases() -> dict:
    out = {}
    for name, want in FROZEN_SHA.items():
        p = os.path.join(KIN if name.endswith("matrix.csv") else PDSP, name)
        got = _sha256(p)
        out[name] = {"expected": want, "observed": got, "match": got == want}
        if got != want:
            raise RuntimeError(f"release drift for {name}: {got} != {want}")
    return out


# --------------------------------------------------------------------------
# Metz
# --------------------------------------------------------------------------
def _peel_to_density(obs: np.ndarray, target: float):
    """Greedy peel; ties broken toward columns (frozen in PREREG section 4)."""
    ri = np.arange(obs.shape[0])
    ci = np.arange(obs.shape[1])
    O = obs.copy()
    while O.mean() < target and O.shape[0] > 5 and O.shape[1] > 5:
        rm = O.mean(axis=1)
        cm = O.mean(axis=0)
        if (1 - rm.min()) * O.shape[1] >= (1 - cm.min()) * O.shape[0]:
            j = int(np.argmin(rm))
            O = np.delete(O, j, axis=0)
            ri = np.delete(ri, j)
        else:
            j = int(np.argmin(cm))
            O = np.delete(O, j, axis=1)
            ci = np.delete(ci, j)
    return ri, ci


def load_metz(density: float = 0.60):
    """Return (Y, mask, compound_ids, kinase_symbols) for BLK-METZ-<density>."""
    df = pd.read_csv(os.path.join(KIN, "metz_matrix.csv"), low_memory=False)
    cid = df.iloc[:, 0].to_numpy()
    mat = df.iloc[:, 1:]
    kin = np.array([c.strip().upper() for c in mat.columns])
    V = mat.to_numpy(dtype=float)
    obs = V > METZ_FLOOR + 1e-9
    ri, ci = _peel_to_density(obs, density)
    return V[np.ix_(ri, ci)], obs[np.ix_(ri, ci)], cid[ri], kin[ci]


def load_klaeger():
    df = pd.read_csv(os.path.join(KIN, "klaeger_matrix.csv"), low_memory=False)
    drug = df.iloc[:, 0].astype(str).to_numpy()
    mat = df.iloc[:, 1:]
    kin = np.array([c.strip().upper() for c in mat.columns])
    W = mat.to_numpy(dtype=float)
    hit = W > KLAEGER_FLOOR + 1e-9
    return W, hit, drug, kin


# --------------------------------------------------------------------------
# PDSP
# --------------------------------------------------------------------------
def load_pdsp():
    """Uncensored human rows of the PDSP Ki database as long-format pKi."""
    df = pd.read_csv(os.path.join(PDSP, "KiDatabase.csv"), low_memory=False,
                     encoding="latin-1")
    df.columns = [c.strip() for c in df.columns]
    v = pd.to_numeric(df["ki Val"], errors="coerce")
    keep = (
        np.isfinite(v)
        & (v > 0)
        & df["ki Note"].isna()                     # drop '>' and '<' censored rows
        & df["species"].astype(str).str.upper().str.strip().eq("HUMAN")
        & df["Unigene"].notna()
    )
    out = df.loc[keep, ["Unigene", "Ligand ID", "Ligand Name", "SMILES",
                        "Hotligand", "source", "Reference"]].copy()
    out["pKi"] = 9.0 - np.log10(v[keep].to_numpy())
    out = out.rename(columns={"Unigene": "target", "Ligand ID": "ligand"})
    return out.reset_index(drop=True)


# --------------------------------------------------------------------------
# KLIFS protein annotation
# --------------------------------------------------------------------------
def load_klifs():
    rec = json.load(open(os.path.join(ANN, "klifs_kinase_information_human.json"),
                         encoding="utf-8"))
    rows = []
    for r in rec:
        pocket = r.get("pocket") or ""
        if len(pocket) != 85:
            continue
        rows.append({
            "name": (r.get("name") or "").strip().upper(),
            "hgnc": (r.get("HGNC") or "").strip().upper(),
            "family": r.get("family") or "",
            "group": r.get("group") or "",
            "uniprot": r.get("uniprot") or "",
            "pocket": pocket,
        })
    return pd.DataFrame(rows)


# legacy gene symbols used by the 2011 Metz panel, resolved against UniProt
SYMBOL_ALIASES = {
    "PRKCN": "PRKD3",
    "STK12": "AURKB",
    "STK6": "AURKA",
    "SGK": "SGK1",
    "KIAA1811": "BRSK1",
}


def map_kinases(symbols, klifs: pd.DataFrame):
    """Map panel kinase symbols onto KLIFS records (HGNC first, then name)."""
    by_hgnc = {}
    by_name = {}
    for _, r in klifs.iterrows():
        by_hgnc.setdefault(r["hgnc"], r.to_dict())
        by_name.setdefault(r["name"], r.to_dict())
    hit, miss = {}, []
    for raw in symbols:
        key = raw.strip().upper()
        s = SYMBOL_ALIASES.get(key, key)
        r = by_hgnc.get(s)
        if r is None:
            r = by_name.get(s)
        if r is None:
            miss.append(key)
        else:
            hit[key] = r          # keyed by the panel's own symbol
    return hit, miss


# --------------------------------------------------------------------------
# additive (two-way) fit with missing cells
# --------------------------------------------------------------------------
def additive_fit(Y, mask, iters: int = 200, tol: float = 1e-10):
    """Least-squares mu + alpha_i + beta_j on observed cells (centred)."""
    Y = np.asarray(Y, float)
    M = np.asarray(mask, bool)
    n, p = Y.shape
    mu = Y[M].mean()
    a = np.zeros(n)
    b = np.zeros(p)
    rc = M.sum(axis=1)
    cc = M.sum(axis=0)
    for _ in range(iters):
        prev = (a.copy(), b.copy(), mu)
        R = np.where(M, Y - mu - b[None, :], 0.0)
        a = np.divide(R.sum(axis=1), np.maximum(rc, 1))
        a[rc == 0] = 0.0
        R = np.where(M, Y - mu - a[:, None], 0.0)
        b = np.divide(R.sum(axis=0), np.maximum(cc, 1))
        b[cc == 0] = 0.0
        a -= a.mean()
        b -= b.mean()
        mu = (Y - a[:, None] - b[None, :])[M].mean()
        if (abs(a - prev[0]).max() < tol and abs(b - prev[1]).max() < tol
                and abs(mu - prev[2]) < tol):
            break
    fit = mu + a[:, None] + b[None, :]
    return mu, a, b, fit
