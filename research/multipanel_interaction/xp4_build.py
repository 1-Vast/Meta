"""XP4 — build BLK-BDB-PANELS and its closures. Frozen by PREREG_XP4.md."""
from __future__ import annotations

import hashlib
import io
import json
import os
import subprocess
import zipfile

import numpy as np
import pandas as pd

RAW = r"D:\MetaSieve\dataset\raw\crossed_panels\bindingdb"
CACHE = r"D:\MetaSieve\dataset\processed\multipanel"
MMSEQS = r"D:\MetaSieve\tools\mmseqs2\mmseqs\bin\mmseqs.exe"
os.makedirs(CACHE, exist_ok=True)
ZIP = os.path.join(RAW, "BindingDB_BindingDB_Articles_202608_tsv.zip")

CONSUMED_PMID = {"21572424", "29191878", "22037378", "18183025", "21949673", "28767711"}
MIN_TARGETS, MIN_LIGANDS = 2, 20

SMI = "Ligand SMILES"
UNI = "UniProt (SwissProt) Primary ID of Target Chain 1"
SEQ = "BindingDB Target Chain Sequence 1"
PMID = "PMID"
NCH = "Number of Protein Chains in Target (>1 implies a multichain complex)"
ORG = "Target Source Organism According to Curator or DataSource"
KI = "Ki (nM)"


def build(force=False):
    out = os.path.join(CACHE, "blk_bdb_panels.npz")
    man_p = os.path.join(CACHE, "blk_bdb_panels_manifest.json")
    if os.path.exists(out) and not force:
        return np.load(out, allow_pickle=True), json.load(open(man_p))

    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")

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

    # affinity values are read from here on; the panel index is already fixed
    v = pd.to_numeric(d[KI], errors="coerce")
    d = d[np.isfinite(v) & (v > 0)]
    d["pKi"] = 9.0 - np.log10(pd.to_numeric(d[KI], errors="coerce")[d.index].to_numpy())

    ok = d[SMI].map(lambda s: Chem.MolFromSmiles(s) is not None)
    d = d[ok]
    cell = d.groupby([PMID, UNI, SMI], as_index=False)["pKi"].mean()
    keep = []
    for pm, g in cell.groupby(PMID):
        if g[UNI].nunique() >= MIN_TARGETS and g[SMI].nunique() >= MIN_LIGANDS:
            keep.append(pm)
    cell = cell[cell[PMID].isin(keep)].reset_index(drop=True)

    # scaffold components over the retained ligands
    ligs = sorted(cell[SMI].unique())
    scaf = {s: MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(s)) for s in ligs}
    uq = sorted(set(scaf.values()))
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
            if fps[j] is not None and DataStructs.TanimotoSimilarity(fps[i], fps[j]) >= 0.5:
                a, b = find(i), find(j)
                if a != b:
                    lab[b] = a
    sc_of = {s: f"sc{find(uq.index(scaf[s]))}" for s in ligs}
    cell["scaffold_component"] = cell[SMI].map(sc_of)

    # protein clusters by MMseqs2 at 40% identity
    seqs = d.drop_duplicates(UNI).set_index(UNI)[SEQ].to_dict()
    seqs = {k: str(v) for k, v in seqs.items() if k in set(cell[UNI])}
    wd = os.path.join(CACHE, "mmseqs40")
    os.makedirs(wd, exist_ok=True)
    fa = os.path.join(wd, "in.fasta")
    with open(fa, "w") as f:
        for k, s in seqs.items():
            f.write(f">{k}\n{s}\n")
    db, clu, tmp, tsv = (os.path.join(wd, x) for x in ("db", "clu", "tmp", "cluster.tsv"))
    for cmd in ([MMSEQS, "createdb", fa, db, "-v", "1"],
                [MMSEQS, "cluster", db, clu, tmp, "--min-seq-id", "0.4", "-c", "0.5",
                 "--cov-mode", "0", "-v", "1"],
                [MMSEQS, "createtsv", db, db, clu, tsv, "-v", "1"]):
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"{cmd[1]}: {r.stderr[-800:]}")
    pcl = {}
    with open(tsv) as f:
        for line in f:
            rep, mem = line.strip().split("\t")
            pcl[mem] = rep
    cell["protein_cluster"] = cell[UNI].map(pcl)

    man = {
        "panel": "BLK-BDB-PANELS",
        "release": "BindingDB_BindingDB_Articles_202608", "license": "CC BY 3.0",
        "excluded_pmids": sorted(CONSUMED_PMID),
        "filters": {"endpoint": "Ki (nM) -> pKi", "single_chain": True,
                    "organism": "Homo sapiens",
                    "min_targets_per_panel": MIN_TARGETS,
                    "min_ligands_per_panel": MIN_LIGANDS},
        "panels": int(cell[PMID].nunique()), "cells": int(len(cell)),
        "targets": int(cell[UNI].nunique()), "ligands": int(cell[SMI].nunique()),
        "scaffold_components": int(cell["scaffold_component"].nunique()),
        "protein_clusters_mmseqs40": int(cell["protein_cluster"].nunique()),
        "median_targets_per_panel": float(cell.groupby(PMID)[UNI].nunique().median()),
        "median_ligands_per_panel": float(cell.groupby(PMID)[SMI].nunique().median()),
        "index_sha256": hashlib.sha256(
            ("|".join(sorted(cell[PMID].astype(str) + "::" + cell[UNI] + "::" + cell[SMI]))
             ).encode()).hexdigest(),
    }
    np.savez(out, pmid=cell[PMID].to_numpy().astype(str), uni=cell[UNI].to_numpy().astype(str),
             smiles=cell[SMI].to_numpy().astype(str), pki=cell["pKi"].to_numpy(),
             scaffold_component=cell["scaffold_component"].to_numpy().astype(str),
             protein_cluster=cell["protein_cluster"].to_numpy().astype(str))
    json.dump(man, open(man_p, "w"), indent=2)
    json.dump({k: seqs[k] for k in sorted(seqs)},
              open(os.path.join(CACHE, "bdb_sequences.json"), "w"))
    return np.load(out, allow_pickle=True), man


if __name__ == "__main__":
    _, m = build(force=True)
    print(json.dumps(m, indent=2))
