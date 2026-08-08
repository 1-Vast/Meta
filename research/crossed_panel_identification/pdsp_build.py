"""Build the BLK-PDSP-H crossed core and its protein annotation."""
from __future__ import annotations

import json
import os
import ssl
import subprocess
import sys
import urllib.parse
import urllib.request

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from panels import load_pdsp  # noqa: E402

CACHE = r"D:\MetaSieve\dataset\processed\crossed_panels"
MMSEQS = r"D:\MetaSieve\tools\mmseqs2\mmseqs\bin\mmseqs.exe"
os.makedirs(CACHE, exist_ok=True)
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE


def crossed_core(df, min_lig_per_target=40, min_target_per_lig=3):
    d = df.copy()
    for _ in range(30):
        tl = d.groupby("target")["ligand"].nunique()
        keep_t = set(tl[tl >= min_lig_per_target].index)
        d = d[d["target"].isin(keep_t)]
        lt = d.groupby("ligand")["target"].nunique()
        keep_l = set(lt[lt >= min_target_per_lig].index)
        d2 = d[d["ligand"].isin(keep_l)]
        if len(d2) == len(d):
            break
        d = d2
    return d


def to_matrix(d):
    cell = d.groupby(["ligand", "target"])["pKi"].mean().reset_index()
    ligs = sorted(cell["ligand"].unique())
    tgts = sorted(cell["target"].unique())
    li = {v: i for i, v in enumerate(ligs)}
    ti = {v: i for i, v in enumerate(tgts)}
    Y = np.full((len(ligs), len(tgts)), np.nan)
    for l, t, v in cell.itertuples(index=False):
        Y[li[l], ti[t]] = v
    return Y, np.isfinite(Y), np.array(ligs), np.array(tgts)


def fetch_uniprot_by_gene(genes):
    path = os.path.join(CACHE, "pdsp_uniprot.json")
    got = json.load(open(path)) if os.path.exists(path) else {}
    todo = [g for g in genes if g not in got]
    for i, g in enumerate(todo):
        q = f"gene_exact:{g} AND organism_id:9606 AND reviewed:true"
        url = ("https://rest.uniprot.org/uniprotkb/search?query="
               + urllib.parse.quote(q) + "&fields=accession,sequence&format=tsv&size=1")
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            txt = urllib.request.urlopen(req, timeout=60, context=CTX).read().decode()
            lines = [l for l in txt.splitlines() if l.strip()][1:]
            got[g] = {"acc": lines[0].split("\t")[0], "seq": lines[0].split("\t")[1]} if lines else None
        except Exception as e:
            print(f"  !! {g}: {type(e).__name__}")
            got[g] = None
        if (i + 1) % 25 == 0:
            print(f"  uniprot {i+1}/{len(todo)}")
            json.dump(got, open(path, "w"))
    json.dump(got, open(path, "w"))
    return got


def mmseqs_clusters(seqs: dict, min_seq_id=0.4, cov=0.5):
    """Single-linkage-ish homology components via mmseqs easy-cluster."""
    wd = os.path.join(CACHE, f"mmseqs_{int(min_seq_id*100)}")
    os.makedirs(wd, exist_ok=True)
    fa = os.path.join(wd, "in.fasta")
    with open(fa, "w") as f:
        for k, s in seqs.items():
            f.write(f">{k}\n{s}\n")
    # explicit pipeline: easy-cluster segfaults in result2flat on this build
    db = os.path.join(wd, "db")
    clu = os.path.join(wd, "clu")
    tmp = os.path.join(wd, "tmp")
    tsv = os.path.join(wd, "cluster.tsv")
    for cmd in (
        [MMSEQS, "createdb", fa, db, "-v", "1"],
        [MMSEQS, "cluster", db, clu, tmp, "--min-seq-id", str(min_seq_id),
         "-c", str(cov), "--cov-mode", "0", "-v", "1"],
        [MMSEQS, "createtsv", db, db, clu, tsv, "-v", "1"],
    ):
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"{cmd[1]} failed:\n{r.stdout[-1500:]}\n{r.stderr[-1500:]}")
    lab = {}
    with open(tsv) as f:
        for line in f:
            rep, mem = line.strip().split("\t")
            lab[mem] = rep
    return lab


if __name__ == "__main__":
    p = load_pdsp()
    print("uncensored human rows:", len(p))
    core = crossed_core(p)
    Y, M, ligs, tgts = to_matrix(core)
    print(f"crossed core: {Y.shape[0]} ligands x {Y.shape[1]} targets, "
          f"{int(M.sum())} cells, density {M.mean():.4f}")
    print("ligands per target: median", int(np.median(M.sum(0))), "min", int(M.sum(0).min()))
    print("targets per ligand: median", int(np.median(M.sum(1))), "min", int(M.sum(1).min()))

    ud = fetch_uniprot_by_gene([str(t) for t in tgts])
    have = {t: ud[t] for t in map(str, tgts) if ud.get(t)}
    print(f"uniprot resolved: {len(have)}/{len(tgts)}")
    lab = mmseqs_clusters({k: v["seq"] for k, v in have.items()}, 0.40)
    ncl = len(set(lab.values()))
    print(f"mmseqs 40% clusters: {ncl} over {len(lab)} targets")

    keep = np.array([str(t) in have and str(t) in lab for t in tgts])
    np.savez(os.path.join(CACHE, "pdsp_core.npz"),
             Y=Y[:, keep], M=M[:, keep], ligands=ligs, targets=tgts[keep],
             cluster=np.array([lab[str(t)] for t in tgts[keep]]),
             accession=np.array([have[str(t)]["acc"] for t in tgts[keep]]))
    with open(os.path.join(CACHE, "pdsp_sequences.json"), "w") as f:
        json.dump({str(t): have[str(t)]["seq"] for t in tgts[keep]}, f)
    print("wrote pdsp_core.npz with", int(keep.sum()), "targets")
