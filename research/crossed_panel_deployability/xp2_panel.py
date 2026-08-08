"""BLK-METZ-XP2: the structure-bearing crossed panel, built from the journal
supplement itself.  Frozen by PREREG_XP2.md sections 3-6.

Nothing here fits a model or computes a contrast.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
XP1 = os.path.join(os.path.dirname(HERE), "crossed_panel_identification")
sys.path.insert(0, XP1)
RAW = r"D:\MetaSieve\dataset\raw\crossed_panels"
CACHE = r"D:\MetaSieve\dataset\processed\crossed_panels_xp2"
os.makedirs(CACHE, exist_ok=True)

METZ_XLS = os.path.join(RAW, "kinase_panels", "metz.xls")
METZ_SHA = "81731c4004823bd45fa3898e25d6491d799dfd0e0486fcc8c9c821f9419dd591"
META_COLS = ["Cmpd_ID", "PUBCHEM_SID", "Canonical_Smiles", "External_Cmpd_ID",
             "External_Source", "Cluster", "ClusterSize", "Cluster_MCSS",
             "Molecular_Weight", "ALogP", "Num_H_Acceptors", "Num_H_Donors",
             "tPSA", "Promiscuity_1uM"]
MIN_KIN_PER_CMPD = 10
MIN_CMPD_PER_KIN = 50
SCAFFOLD_MERGE_TANIMOTO = 0.5


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def _parse_supplement():
    assert sha256_file(METZ_XLS) == METZ_SHA, "metz.xls release drift"
    d = pd.read_excel(METZ_XLS, sheet_name=0)
    kin = [c for c in d.columns if c not in META_COLS]
    n, p = len(d), len(kin)
    meas = np.zeros((n, p), bool)
    val = np.full((n, p), np.nan)
    thr = np.full((n, p), np.nan)          # left-censoring threshold where known
    for j, c in enumerate(kin):
        col = d[c].values
        for i, v in enumerate(col):
            if isinstance(v, str):
                s = v.strip()
                if s.startswith("<"):
                    thr[i, j] = float(s[1:])
            elif v is not None and not (isinstance(v, float) and np.isnan(v)):
                meas[i, j] = True
                val[i, j] = float(v)
    return d, np.array(kin), meas, val, thr


def build(force=False):
    out = os.path.join(CACHE, "blk_metz_xp2.npz")
    man = os.path.join(CACHE, "blk_metz_xp2_manifest.json")
    if os.path.exists(out) and not force:
        return np.load(out, allow_pickle=True), json.load(open(man))

    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
    from panels import load_klifs, map_kinases

    d, kin, meas, val, thr = _parse_supplement()
    smi = d["Canonical_Smiles"].to_numpy()
    parsable = np.array([isinstance(s, str) and Chem.MolFromSmiles(s) is not None
                         for s in smi])

    keep_c = parsable.copy()
    keep_k = np.ones(len(kin), bool)
    for _ in range(50):
        nk = meas[np.ix_(keep_c, keep_k)].sum(axis=1)
        idx_c = np.where(keep_c)[0]
        keep_c2 = keep_c.copy()
        keep_c2[idx_c[nk < MIN_KIN_PER_CMPD]] = False
        nc = meas[np.ix_(keep_c2, keep_k)].sum(axis=0)
        idx_k = np.where(keep_k)[0]
        keep_k2 = keep_k.copy()
        keep_k2[idx_k[nc < MIN_CMPD_PER_KIN]] = False
        if (keep_c2 == keep_c).all() and (keep_k2 == keep_k).all():
            break
        keep_c, keep_k = keep_c2, keep_k2

    klifs = load_klifs()
    hit, miss = map_kinases(kin[keep_k], klifs)
    ok_k = np.array([k in hit for k in kin[keep_k]])
    idx_k = np.where(keep_k)[0][ok_k]
    idx_c = np.where(keep_c)[0]

    Y = val[np.ix_(idx_c, idx_k)]
    M = meas[np.ix_(idx_c, idx_k)]
    TH = thr[np.ix_(idx_c, idx_k)]
    cmpd = d["Cmpd_ID"].to_numpy()[idx_c]
    smiles = smi[idx_c]
    kinase = kin[idx_k]

    # --- ligand scaffold components -------------------------------------
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    scaf = [MurckoScaffold.MurckoScaffoldSmiles(mol=m) for m in mols]
    uniq = sorted(set(scaf))
    smol = [Chem.MolFromSmiles(s) if s else None for s in uniq]
    from rdkit.Chem import rdFingerprintGenerator
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(m) if m is not None and m.GetNumAtoms() else None
           for m in smol]
    from rdkit import DataStructs
    lab = list(range(len(uniq)))

    def find(x):
        while lab[x] != x:
            lab[x] = lab[lab[x]]
            x = lab[x]
        return x

    for i in range(len(uniq)):
        if fps[i] is None:
            continue
        for j in range(i + 1, len(uniq)):
            if fps[j] is None:
                continue
            if DataStructs.TanimotoSimilarity(fps[i], fps[j]) >= SCAFFOLD_MERGE_TANIMOTO:
                a, b = find(i), find(j)
                if a != b:
                    lab[b] = a
    comp_of_scaffold = {s: f"sc{find(i)}" for i, s in enumerate(uniq)}
    scaf_comp = np.array([comp_of_scaffold[s] for s in scaf])

    prot_group = np.array([hit[k]["group"] for k in kinase])
    prot_family = np.array([hit[k]["family"] for k in kinase])
    pocket = np.array([hit[k]["pocket"] for k in kinase])
    uniprot = np.array([hit[k]["uniprot"] for k in kinase])

    np.savez(out, Y=Y, M=M, THR=TH, cmpd=cmpd, smiles=smiles, scaffold=np.array(scaf),
             scaffold_component=scaf_comp, kinase=kinase, group=prot_group,
             family=prot_family, pocket=pocket, uniprot=uniprot)
    manifest = {
        "panel": "BLK-METZ-XP2",
        "source": "metz.xls Table S1 (journal supplement)",
        "source_sha256": METZ_SHA,
        "filters": {"min_kinases_per_compound": MIN_KIN_PER_CMPD,
                    "min_compounds_per_kinase": MIN_CMPD_PER_KIN,
                    "requires_parsable_smiles": True,
                    "requires_klifs_mapping": True},
        "n_compounds": int(len(cmpd)), "n_kinases": int(len(kinase)),
        "measured_cells": int(M.sum()), "density": float(M.mean()),
        "n_scaffolds": int(len(set(scaf))),
        "n_scaffold_components": int(len(set(scaf_comp))),
        "scaffold_merge_tanimoto": SCAFFOLD_MERGE_TANIMOTO,
        "n_groups": int(len(set(prot_group))),
        "kinases_dropped_no_klifs": [str(k) for k in kin[keep_k][~ok_k]],
        "censored_cells_in_block": int(np.isfinite(TH).sum()),
        "distinct_censoring_thresholds": int(len(np.unique(TH[np.isfinite(TH)]))),
        "index_sha256": hashlib.sha256(
            (",".join(map(str, cmpd)) + "|" + ",".join(kinase)).encode()).hexdigest(),
        "panel_sha256": sha256_file(out),
    }
    json.dump(manifest, open(man, "w"), indent=2)
    return np.load(out, allow_pickle=True), manifest


# -------------------------------------------------------------- ligand features
def ligand_features(smiles, kinds=("desc", "ecfp", "chemberta", "random")):
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    F = {}
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    if "desc" in kinds:
        from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors as rd
        rows = []
        for m in mols:
            rows.append([
                Descriptors.MolWt(m), Crippen.MolLogP(m), rd.CalcTPSA(m),
                rd.CalcNumHBD(m), rd.CalcNumHBA(m), rd.CalcNumRotatableBonds(m),
                rd.CalcNumRings(m), rd.CalcNumAromaticRings(m),
                rd.CalcFractionCSP3(m), m.GetNumHeavyAtoms(),
            ])
        X = np.array(rows, float)
        F["L-DESC"] = (X - X.mean(0)) / (X.std(0) + 1e-9)
    if "ecfp" in kinds:
        from rdkit.Chem import rdFingerprintGenerator
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
        X = np.array([list(gen.GetFingerprint(m)) for m in mols], float)
        F["L-ECFP"] = X
    if "chemberta" in kinds:
        F["L-CHEMBERTA"] = _chemberta(smiles)
    if "random" in kinds:
        F["L-RANDOM"] = np.random.default_rng(0).normal(size=(len(smiles), 64))
    return F


def _chemberta(smiles, batch=32):
    path = os.path.join(CACHE, "chemberta_ligand.npz")
    if os.path.exists(path):
        d = np.load(path, allow_pickle=True)
        m = {k: v for k, v in zip(d["keys"], d["emb"])}
        if all(s in m for s in smiles):
            return np.stack([m[s] for s in smiles])
    import torch
    from transformers import AutoModel, AutoTokenizer
    name = "DeepChem/ChemBERTa-77M-MLM"
    tok = AutoTokenizer.from_pretrained(name)
    mdl = AutoModel.from_pretrained(name).eval()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    mdl = mdl.to(dev)
    embs = []
    with torch.no_grad():
        for i in range(0, len(smiles), batch):
            ch = list(smiles[i:i + batch])
            enc = tok(ch, return_tensors="pt", padding=True, truncation=True,
                      max_length=256).to(dev)
            out = mdl(**enc).last_hidden_state
            mk = enc["attention_mask"].unsqueeze(-1).float()
            embs.append(((out * mk).sum(1) / mk.sum(1)).cpu().numpy())
    E = np.concatenate(embs)
    np.savez(path, keys=np.array(list(smiles)), emb=E)
    return E


if __name__ == "__main__":
    d, man = build(force=True)
    print(json.dumps(man, indent=2))
