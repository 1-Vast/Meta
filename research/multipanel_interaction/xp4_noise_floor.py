"""XP4-N — is the within-panel interaction reproducible above measurement noise?

Two independent estimates, neither of which fits a model:

  N1  replicate-based.  BindingDB records multiple rows for many (PMID, target,
      ligand) cells.  Splitting those rows in half gives a direct estimate of the
      per-report measurement sd, and hence of the noise transmitted into gamma by
      within-panel double centring.

  N2  chemistry-based ceiling.  For a held-out ligand in a panel, predict its
      gamma profile as the ECFP-similarity-weighted average of the OTHER ligands'
      gamma profiles in the same panel, restricted to training scaffolds.  This
      is the most favourable non-parametric use of the panel's own data and needs
      no cross-panel generalisation at all.  If it fails, no model can succeed.
"""
from __future__ import annotations

import io
import json
import os
import sys
import zipfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from xp4_build import CACHE, CONSUMED_PMID, KI, NCH, ORG, PMID, SEQ, SMI, UNI, build  # noqa: E402

RAW = r"D:\MetaSieve\dataset\raw\crossed_panels\bindingdb"
ZIP = os.path.join(RAW, "BindingDB_BindingDB_Articles_202608_tsv.zip")
OUT = r"D:\MetaSieve\report\multipanel_interaction"
os.makedirs(OUT, exist_ok=True)
rep = {"stage": "XP4-N"}

# ------------------------------------------------------------------ N1
with zipfile.ZipFile(ZIP) as z:
    name = [n for n in z.namelist() if n.endswith(".tsv")][0]
    with z.open(name) as f:
        df = pd.read_csv(io.TextIOWrapper(f, "utf-8", errors="ignore"), sep="\t",
                         usecols=[SMI, UNI, SEQ, PMID, NCH, ORG, KI],
                         low_memory=False, on_bad_lines="skip")
d = df[df[KI].notna()].copy()
d[PMID] = d[PMID].astype("string").str.strip()
d = d.dropna(subset=[SMI, UNI, SEQ, PMID])
d = d[d[NCH].fillna(1).astype(float) <= 1]
d = d[d[ORG].astype(str).str.contains("Homo sapiens", na=False)]
d = d[~d[PMID].isin(CONSUMED_PMID)]
v = pd.to_numeric(d[KI], errors="coerce")
d = d[np.isfinite(v) & (v > 0)].copy()
d["pKi"] = 9.0 - np.log10(pd.to_numeric(d[KI], errors="coerce")[d.index].to_numpy())

grp = d.groupby([PMID, UNI, SMI])["pKi"]
sizes = grp.size()
multi = sizes[sizes >= 2]
rng = np.random.default_rng(0)
diffs = []
for key in multi.index:
    vals = grp.get_group(key).to_numpy()
    rng.shuffle(vals)
    h = len(vals) // 2
    diffs.append(vals[:h].mean() - vals[h:].mean())
diffs = np.array(diffs)
sd_pair = float(diffs.std(ddof=1)) if len(diffs) > 2 else float("nan")
sigma_rep = sd_pair / np.sqrt(2.0)
rep["N1_replicates"] = {
    "replicated_cells": int(len(multi)),
    "sd_half_difference": sd_pair,
    "per_report_sigma_log_units": float(sigma_rep),
}
print(f"N1  replicated cells {len(multi)}, sd(h1-h2) {sd_pair:.4f} "
      f"-> per-report sigma {sigma_rep:.4f} log units")

# ------------------------------------------------------------------ N2
data, man = build()
pmid, uni, smi = data["pmid"], data["uni"], data["smiles"]
y = data["pki"].astype(float)
scomp = data["scaffold_component"]
lig_ids = {s: i for i, s in enumerate(sorted(set(smi)))}
tgt_ids = {t: i for i, t in enumerate(sorted(set(uni)))}
lid = np.array([lig_ids[s] for s in smi])
tid = np.array([tgt_ids[t] for t in uni])
sys.path.insert(0, HERE)
from xp4_run import ecfp, within_panel_center  # noqa: E402
gamma = within_panel_center(y, pmid, lid, tid)
gsd = float(gamma.std(ddof=1))

# how much of gamma's variance can noise alone account for?
Tbar = float(np.mean([len(set(uni[pmid == p])) for p in set(pmid)]))
noise_into_gamma = sigma_rep * np.sqrt(max(1.0 - 1.0 / Tbar, 0.0))
rep["N2_variance_budget"] = {
    "gamma_sd": gsd, "raw_pki_sd": float(y.std(ddof=1)),
    "mean_targets_per_panel": Tbar,
    "noise_sd_transmitted_into_gamma": float(noise_into_gamma),
    "implied_noise_share_of_gamma_variance": float(min((noise_into_gamma / gsd) ** 2, 1.0)),
}
print(f"N2  gamma sd {gsd:.4f}; noise alone would give {noise_into_gamma:.4f} "
      f"({rep['N2_variance_budget']['implied_noise_share_of_gamma_variance']*100:.0f}% "
      f"of gamma variance)")

# chemistry ceiling: similarity-weighted neighbour transfer inside each panel
F = ecfp(list(smi))
Fn = F / (np.linalg.norm(F, axis=1, keepdims=True) + 1e-9)
scomps = sorted(set(scomp))
rs = np.random.default_rng(1)
order = list(scomps)
rs.shuffle(order)
sfold = {c: i % 5 for i, c in enumerate(order)}
sse_knn = sse_zero = 0.0
n_cell = 0
for p in sorted(set(pmid)):
    m = np.where(pmid == p)[0]
    for f in range(5):
        te = m[np.array([sfold[scomp[i]] == f for i in m])]
        tr = m[np.array([sfold[scomp[i]] != f for i in m])]
        if len(te) < 3 or len(tr) < 10:
            continue
        for i in te:
            same_t = tr[tid[tr] == tid[i]]
            if len(same_t) < 3:
                continue
            sim = Fn[same_t] @ Fn[i]
            w = np.clip(sim, 0, None) ** 3
            pred = float((w * gamma[same_t]).sum() / max(w.sum(), 1e-9))
            sse_knn += (gamma[i] - pred) ** 2
            sse_zero += gamma[i] ** 2
            n_cell += 1
r2_knn = 1.0 - sse_knn / max(sse_zero, 1e-12)
rep["N2_chemistry_ceiling"] = {
    "cells_scored": int(n_cell),
    "r2_gamma_similarity_neighbour_within_panel": float(r2_knn),
    "interpretation": "most favourable non-parametric use of the panel's own data; "
                      "requires no cross-panel generalisation",
}
print(f"N2  within-panel chemistry-neighbour ceiling R2_gamma = {r2_knn:+.4f} "
      f"on {n_cell} cells")

json.dump(rep, open(os.path.join(OUT, "XP4_NOISE_FLOOR.json"), "w"), indent=2,
          default=float)
print("wrote XP4_NOISE_FLOOR.json")
