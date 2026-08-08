"""XP3 STAGE 1 — data and observability census.

No model is fitted here and no arm is scored.  The question is whether ANY
release-pinned public panel can support the estimand that XP2 could not:

    a deployment-observable basis x(P,L) whose protein-specific interaction
    component survives simultaneous protein-group and ligand-scaffold closure
    at k <= 5.

For each candidate panel we report the design facts that determine whether the
question is answerable at all, and then the minimum detectable interaction
effect implied by the design and the measurement noise.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "crossed_panel_identification"))
sys.path.insert(0, os.path.join(ROOT, "crossed_panel_deployability"))
RAW = r"D:\MetaSieve\dataset\raw\crossed_panels"
OUT = r"D:\MetaSieve\report\observability_census"
os.makedirs(OUT, exist_ok=True)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def scaffold_components(smiles, thr=0.5):
    from rdkit import Chem, DataStructs, RDLogger, rdBase  # noqa: F401
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
    mols = [Chem.MolFromSmiles(s) if isinstance(s, str) else None for s in smiles]
    scaf = [MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m is not None else None
            for m in mols]
    uq = sorted({s for s in scaf if s is not None})
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = []
    for s in uq:
        m = Chem.MolFromSmiles(s)
        fps.append(gen.GetFingerprint(m) if m is not None and m.GetNumAtoms() else None)
    lab = list(range(len(uq)))

    def find(x):
        while lab[x] != x:
            lab[x] = lab[lab[x]]
            x = lab[x]
        return x

    for i in range(len(uq)):
        if fps[i] is None:
            continue
        for j in range(i + 1, len(uq)):
            if fps[j] is None:
                continue
            if DataStructs.TanimotoSimilarity(fps[i], fps[j]) >= thr:
                a, b = find(i), find(j)
                if a != b:
                    lab[b] = a
    comp = {s: f"sc{find(i)}" for i, s in enumerate(uq)}
    return [comp.get(s) for s in scaf], len(uq), len({v for v in comp.values()})


def mdi(sigma, n_cells, n_units):
    """Minimum detectable interaction sd at 80% power, alpha 0.05, two-sided,
    using the closure component as the independence unit (design-effect free
    optimistic bound)."""
    if n_units < 2 or n_cells < 2:
        return None
    # detectable mean-square reduction with n_units independent clusters
    return float(2.8 * sigma / np.sqrt(max(n_units, 1)))


report = {"stage": "XP3-1", "panels": {}, "generated": "2026-08-08"}

# ------------------------------------------------------------------ METZ
print("=" * 78)
print("PANEL: METZ 2011 (release-pinned journal supplement)")
from xp2_panel import build as build_metz  # noqa: E402
d, man = build_metz()
Y, M = d["Y"], d["M"]
sc = d["scaffold_component"]
grp = d["group"]
res = np.nan
# residual sd after two-way fit, and the reproducible share from XP1-A
from xp2_core import additive_fit  # noqa: E402
mu, a, b = additive_fit(Y, M)
resid = (Y - mu - a[:, None] - b[None, :])[M]
sigma_tot = float(resid.std(ddof=1))
rect = 0
for j1 in range(Y.shape[1]):
    pass
report["panels"]["METZ"] = {
    "release_sha256": man["source_sha256"],
    "license": "publisher supplementary data (Nature Chemical Biology), public mirror pinned",
    "endpoint": "pKi, single laboratory, single assay family",
    "compounds": man["n_compounds"], "proteins": man["n_kinases"],
    "measured_cells": man["measured_cells"], "density": man["density"],
    "protein_closure_components_group": int(len(set(grp))),
    "protein_closure_components_family": int(len(set(d["family"]))),
    "ligand_scaffold_components": man["n_scaffold_components"],
    "scaffolds": man["n_scaffolds"],
    "exact_sequences": True, "canonical_structures": True,
    "censoring": f"{man['censored_cells_in_block']} left-censored cells at "
                 f"{man['distinct_censoring_thresholds']} distinct thresholds; excluded",
    "replicate_structure": "none within panel (single measurement per cell)",
    "residual_sd_after_two_way": sigma_tot,
    "reproducible_interaction_sd_XP1": 0.442,
    "independent_units_double_closure": int(min(len(set(grp)), 5)) * 1,
    "mdi_group_units": mdi(sigma_tot, int(M.sum()), len(set(grp))),
    "mdi_scaffold_units": mdi(sigma_tot, int(M.sum()), man["n_scaffold_components"]),
    "verdict": "supports the estimand structurally; CONSUMED as XP1/XP2 development panel",
}
for k, v in report["panels"]["METZ"].items():
    print(f"   {k:38s} {v}")

# --------------------------------------------------------------- KLAEGER
print("\n" + "=" * 78)
print("PANEL: KLAEGER 2017 kinobeads")
kp = os.path.join(RAW, "kinase_panels", "klaeger_matrix.csv")
kdf = pd.read_csv(kp, low_memory=False)
W = kdf.iloc[:, 1:].to_numpy(float)
hit = W > 5.0 + 1e-9
smi_map = json.load(open(r"D:\MetaSieve\dataset\processed\crossed_panels_xp2\klaeger_smiles.json"))
drugs = kdf.iloc[:, 0].astype(str).str.strip().tolist()
smis = [(smi_map.get(x) or {}).get("smiles") for x in drugs]
comp_k, n_scaf_k, n_comp_k = scaffold_components(smis)
report["panels"]["KLAEGER"] = {
    "release_sha256": sha(kp),
    "license": "publisher supplementary data (Science), public mirror pinned",
    "endpoint": "apparent pKd (kinobeads chemical proteomics)",
    "compounds": int(W.shape[0]), "proteins": int(W.shape[1]),
    "measured_cells_above_floor": int(hit.sum()),
    "floor_fraction": float(1 - hit.mean()),
    "ligand_scaffold_components": int(n_comp_k), "scaffolds": int(n_scaf_k),
    "structures_resolved": int(sum(1 for s in smis if s)),
    "exact_sequences": True, "canonical_structures": True,
    "replicate_structure": "none published per cell",
    "verdict": "independent platform; 93.3% floor makes it a hit matrix, not a "
               "continuous crossed panel; already CONSUMED as XP1 cross-platform "
               "check and XP2-F external test",
}
for k, v in report["panels"]["KLAEGER"].items():
    print(f"   {k:38s} {v}")

# ------------------------------------------------------------------ PDSP
print("\n" + "=" * 78)
print("PANEL: NIMH PDSP Ki database")
from panels import load_pdsp  # noqa: E402
p = load_pdsp()
cells = p.groupby(["target", "ligand"]).size()
rep = int((cells >= 2).sum())
sm = p.dropna(subset=["SMILES"]).groupby("ligand")["SMILES"].first()
comp_p, n_scaf_p, n_comp_p = scaffold_components(list(sm.values))
report["panels"]["PDSP"] = {
    "release_sha256": sha(os.path.join(RAW, "pdsp_kidb", "KiDatabase.csv")),
    "license": "free public NIMH resource",
    "endpoint": "Ki (nM) from radioligand displacement; heterogeneous radioligands",
    "rows_uncensored_human": int(len(p)),
    "targets": int(p["target"].nunique()), "ligands": int(p["ligand"].nunique()),
    "cells": int(len(cells)), "replicated_cells": rep,
    "ligands_with_structure": int(len(sm)),
    "ligand_scaffold_components": int(n_comp_p), "scaffolds": int(n_scaf_p),
    "distinct_radioligands": int(p["Hotligand"].nunique()),
    "per_report_sigma_log_units": 0.7144,
    "replicate_r": 0.7270,
    "assay_comparability": "radioligand is nearly a function of the target, so an "
                           "assay-by-compound interaction is not separable from gamma",
    "verdict": "continuous and replicated, but the assay confound is structural and "
               "the per-report noise is 0.714 log units; CONSUMED as XP1-C",
}
for k, v in report["panels"]["PDSP"].items():
    print(f"   {k:38s} {v}")

json.dump(report, open(os.path.join(OUT, "XP3_DATA_CENSUS.json"), "w"),
          indent=2, default=float)
print("\nwrote", os.path.join(OUT, "XP3_DATA_CENSUS.json"))
