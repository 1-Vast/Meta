"""S3 + S4 — label-free splits, power design, and the baseline observability audit.

S4 asks one question, and it is deliberately an UPPER BOUND on the whole
sequence+2D input class rather than a test of one checkpoint:

    can the six named 3D mechanism channels, computed from experimental holo
    coordinates, be predicted from deployment-observable inputs (protein
    sequence embedding + 2D ligand structure) under simultaneous protein
    homology and ligand scaffold closure?

P1B's predicted contact/distance geometry is itself a function of exactly these
inputs, so a negative here rules out the whole class, not one model.  A positive
here would then justify running the P1B checkpoint specifically and, after that,
GPU distillation.

No affinity label is read.  No PDB, ligand, target or document identity is a
feature.
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "crossed_panel_deployability"))
from s2_teacher import CHANNELS, extract, teacher  # noqa: E402
from xp2_core import paired_contrast, r2_vs_base, ridge_predict, widest  # noqa: E402

ROOT = r"D:\MetaSieve"
IND = os.path.join(ROOT, "dataset", "raw", "ssl_b2_independent", "mmcif")
CACHE = os.path.join(ROOT, "dataset", "processed", "ssl_b2")
OUT = os.path.join(ROOT, "report", "ssl_b2_structural_observability")
MMSEQS = os.path.join(ROOT, "tools", "mmseqs2", "mmseqs", "bin", "mmseqs.exe")
os.makedirs(CACHE, exist_ok=True)
os.makedirs(OUT, exist_ok=True)
N_FOLDS, LAMS = 5, (1.0, 10.0, 100.0, 1e3, 1e4, 1e5)


def build_dataset(force=False):
    cache_path = os.path.join(CACHE, "teacher_dataset.npz")
    if os.path.exists(cache_path) and not force:
        return np.load(cache_path, allow_pickle=True)
    import gemmi
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

    files = sorted(glob.glob(os.path.join(IND, "*.cif.gz")))
    rows_p = os.path.join(CACHE, f"parsed_rows_{len(files)}.json")
    if os.path.exists(rows_p):
        rows = [{**r, "agg": np.array(r["agg"])} for r in json.load(open(rows_p))]
        print(f"loaded {len(rows)} parsed complexes from cache")
        files = []
    print(f"parsing {len(files)} independent structures")
    rows = rows if not files else []
    for i, f in enumerate(files):
        pid = os.path.basename(f).split(".")[0]
        try:
            prot, lig, comp = extract(f)
            if prot is None:
                continue
            t = teacher(prot, lig)
            if t is None:
                continue
            # protein sequence from the structure itself
            st = gemmi.make_structure_from_block(
                gemmi.cif.read_string(__import__("gzip").open(
                    f, "rt", encoding="utf-8", errors="ignore").read()).sole_block())
            st.setup_entities()
            seq = gemmi.one_letter_code(
                [r.name for ch in st[0] for r in ch
                 if (gemmi.find_tabulated_residue(r.name) or
                     type("x", (), {"is_amino_acid": lambda s: False})()).is_amino_acid()])
            if len(seq) < 40:
                continue
            rows.append({"pdb": pid, "comp": comp, "seq": seq,
                         "agg": t["aggregate"],
                         "n_lig": t["n_ligand_atoms"], "n_prot": t["n_protein_atoms"]})
        except Exception:
            continue
        if (i + 1) % 200 == 0:
            print(f"  parsed {len(rows)}/{i+1}")
    print(f"usable complexes: {len(rows)}")
    if files:
        json.dump([{**r, "agg": list(map(float, r["agg"]))} for r in rows],
                  open(rows_p, "w"))

    # ligand SMILES from the CCD component id via RCSB, cached
    smi_p = os.path.join(CACHE, "ccd_smiles.json")
    smi = json.load(open(smi_p)) if os.path.exists(smi_p) else {}
    import ssl as _ssl
    import urllib.request
    ctx = _ssl.create_default_context()
    ctx.check_hostname = False
    ctx.verify_mode = _ssl.CERT_NONE
    need = sorted({r["comp"] for r in rows} - set(smi))
    print(f"resolving {len(need)} CCD components (parallel)")

    def _one(c):
        for attempt in range(2):
            try:
                u = f"https://files.rcsb.org/ligands/download/{c}.cif"
                txt = urllib.request.urlopen(urllib.request.Request(
                    u, headers={"User-Agent": "Mozilla/5.0"}), timeout=15,
                    context=ctx).read().decode("utf-8", "ignore")
                for line in txt.splitlines():
                    if "SMILES_CANONICAL" in line and "OpenEye" in line:
                        parts = line.split('"')
                        if len(parts) > 1:
                            return c, parts[-2]
                return c, None
            except Exception:
                continue
        return c, None

    if need:
        from concurrent.futures import ThreadPoolExecutor
        done = 0
        with ThreadPoolExecutor(max_workers=8) as ex:
            for c, v in ex.map(_one, need):
                smi[c] = v
                done += 1
                if done % 100 == 0:
                    print(f"  ccd {done}/{len(need)}")
                    json.dump(smi, open(smi_p, "w"))
    json.dump(smi, open(smi_p, "w"))

    keep = [r for r in rows if smi.get(r["comp"])
            and Chem.MolFromSmiles(smi[r["comp"]]) is not None]
    print(f"complexes with a parsable ligand: {len(keep)}")
    L = np.stack([np.array(list(gen.GetFingerprint(
        Chem.MolFromSmiles(smi[r["comp"]]))), float) for r in keep])
    scaf = [MurckoScaffold.MurckoScaffoldSmiles(
        mol=Chem.MolFromSmiles(smi[r["comp"]])) for r in keep]
    Y = np.stack([r["agg"] for r in keep])
    seqs = [r["seq"] for r in keep]
    pdbs = [r["pdb"] for r in keep]
    comps = [r["comp"] for r in keep]

    # protein homology clusters over the independent set, MMseqs 40%
    wd = os.path.join(CACHE, "mmseqs40")
    if os.path.isdir(wd):
        import shutil
        shutil.rmtree(wd, ignore_errors=True)   # stale DBs make `cluster` refuse
    os.makedirs(wd, exist_ok=True)
    fa = os.path.join(wd, "in.fasta")
    with open(fa, "w") as f:
        for p, s in zip(pdbs, seqs):
            f.write(f">{p}\n{s}\n")
    db, clu, tmp, tsv = (os.path.join(wd, x) for x in ("db", "clu", "tmp", "c.tsv"))
    for cmd in ([MMSEQS, "createdb", fa, db, "-v", "1"],
                [MMSEQS, "cluster", db, clu, tmp, "--min-seq-id", "0.4", "-c", "0.8",
                 "--cov-mode", "0", "-v", "1"],
                [MMSEQS, "createtsv", db, db, clu, tsv, "-v", "1"]):
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError(
                f"{cmd[1]} rc={r.returncode} "
                f"STDOUT:{r.stdout[-800:]} STDERR:{r.stderr[-800:]}")
    pc = {}
    with open(tsv) as f:
        for line in f:
            rep, mem = line.strip().split("\t")
            pc[mem] = rep
    pclus = np.array([pc.get(p, p) for p in pdbs])

    np.savez(cache_path, Y=Y, L=L, seqs=np.array(seqs), pdbs=np.array(pdbs),
             comps=np.array(comps), scaffold=np.array(scaf), pclus=pclus)
    return np.load(cache_path, allow_pickle=True)


def esm(seqs, tag="ind"):
    p = os.path.join(CACHE, f"esm2_t30_{tag}.npz")
    if os.path.exists(p):
        z = np.load(p, allow_pickle=True)
        if len(z["emb"]) == len(seqs):
            return z["emb"]
    import torch
    from transformers import AutoTokenizer, EsmModel
    name = "facebook/esm2_t30_150M_UR50D"
    tok = AutoTokenizer.from_pretrained(name)
    mdl = EsmModel.from_pretrained(name).eval()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    mdl = mdl.to(dev)
    out = []
    with torch.no_grad():
        for i in range(0, len(seqs), 4):
            ch = [str(s)[:1022] for s in seqs[i:i + 4]]
            enc = tok(ch, return_tensors="pt", padding=True, truncation=True,
                      max_length=1022).to(dev)
            o = mdl(**enc).last_hidden_state
            mk = enc["attention_mask"].unsqueeze(-1).float()
            out.append(((o * mk).sum(1) / mk.sum(1)).cpu().numpy())
            print(f"  esm {min(i+4, len(seqs))}/{len(seqs)}", end="\r")
    print()
    E = np.concatenate(out)
    np.savez(p, emb=E)
    return E


def main(n_boot=2000):
    d = build_dataset()
    Y, L, pclus, scaf = d["Y"], d["L"], d["pclus"], d["scaffold"]
    P = esm(d["seqs"])
    n = len(Y)
    print(f"\ndataset: {n} complexes, {len(set(pclus))} protein clusters, "
          f"{len(set(scaf))} scaffolds")

    split = {"schema": "MetaSieve.StructuralSplitManifest.v1",
             "independent_source": "RCSB CC0 entries released >= 2024-01-01, "
                                   "disjoint from every pilot20k pdb id",
             "exposed_overlap": 0, "complexes": int(n),
             "protein_clusters_mmseqs40_c80": int(len(set(pclus))),
             "scaffolds": int(len(set(scaf))),
             "closures": ["protein homology (MMseqs 40%/80cov)", "exact sequence",
                          "ligand CCD identity", "Bemis-Murcko scaffold",
                          "PDB entry", "P1B exposure (zero by construction)"],
             "inference_unit": "protein homology cluster and scaffold, never an "
                               "atom-residue pair"}
    json.dump(split, open(os.path.join(OUT, "STRUCTURAL_SPLIT_MANIFEST.json"), "w"),
              indent=2)

    rng = np.random.default_rng(0)
    pcs, scs = sorted(set(pclus)), sorted(set(scaf))
    o1, o2 = list(pcs), list(scs)
    rng.shuffle(o1)
    rng.shuffle(o2)
    pf = {c: i % N_FOLDS for i, c in enumerate(o1)}
    sf = {c: i % N_FOLDS for i, c in enumerate(o2)}

    ARMS = ["B0", "B2-bio", "B3-rand", "B4-deranged"]
    err = {a: {c: {} for c in CHANNELS} for a in ARMS}
    errs = {a: {c: {} for c in CHANNELS} for a in ARMS}
    Ys = (Y - Y.mean(0)) / (Y.std(0) + 1e-9)

    for f in range(N_FOLDS):
        te = np.array([pf[a] == f and sf[b] == f for a, b in zip(pclus, scaf)])
        tr = np.array([pf[a] != f and sf[b] != f for a, b in zip(pclus, scaf)])
        if te.sum() < 5 or tr.sum() < 50:
            continue
        Pz = (P - P[tr].mean(0)) / (P[tr].std(0) + 1e-9)
        Lz = (L - L[tr].mean(0)) / (L[tr].std(0) + 1e-9)
        X = np.hstack([Pz, Lz])
        R = rng.normal(size=X.shape)
        der = rng.permutation(np.where(te)[0])
        Xd = X.copy()
        Xd[np.where(te)[0], :Pz.shape[1]] = Pz[der]
        inner = pclus[tr]
        best, be = LAMS[0], np.inf
        for lam in LAMS:
            e = 0.0
            for c in sorted(set(inner))[:40]:
                m = inner != c
                if m.sum() < 30 or (~m).sum() == 0:
                    continue
                e += float(((ridge_predict(X[tr][m], Ys[tr][m], X[tr][~m], lam)
                             - Ys[tr][~m]) ** 2).sum())
            if e < be:
                be, best = e, lam
        pred = {
            "B0": np.repeat(Ys[tr].mean(0)[None, :], te.sum(), axis=0),
            "B2-bio": ridge_predict(X[tr], Ys[tr], X[te], best),
            "B3-rand": ridge_predict(R[tr], Ys[tr], R[te], best),
            "B4-deranged": ridge_predict(X[tr], Ys[tr], Xd[te], best),
        }
        for a, pr in pred.items():
            for ci, c in enumerate(CHANNELS):
                e = Ys[te][:, ci] - pr[:, ci]
                for u in set(pclus[te]):
                    m = pclus[te] == u
                    s, k = err[a][c].get(u, (0.0, 0))
                    err[a][c][u] = (s + float((e[m] ** 2).sum()), k + int(m.sum()))
                for u in set(scaf[te]):
                    m = scaf[te] == u
                    s, k = errs[a][c].get(u, (0.0, 0))
                    errs[a][c][u] = (s + float((e[m] ** 2).sum()), k + int(m.sum()))

    res = {"stage": "S4", "note": "upper bound on the sequence+2D input class",
           "complexes_scored": int(sum(v[1] for v in err["B0"][CHANNELS[0]].values())),
           "channels": {}}
    print(f"\n{'channel':24s} {'B2-bio R2':>26s} {'vs random':>24s} {'vs deranged':>24s}")
    for c in CHANNELS:
        if not err["B0"][c]:
            continue
        r2 = widest(r2_vs_base(err["B2-bio"][c], err["B0"][c], n_boot, 11),
                    r2_vs_base(errs["B2-bio"][c], errs["B0"][c], n_boot, 12))
        vr = widest(paired_contrast(err["B2-bio"][c], err["B3-rand"][c], n_boot, 21),
                    paired_contrast(errs["B2-bio"][c], errs["B3-rand"][c], n_boot, 22))
        vd = widest(paired_contrast(err["B2-bio"][c], err["B4-deranged"][c], n_boot, 31),
                    paired_contrast(errs["B2-bio"][c], errs["B4-deranged"][c], n_boot, 32))
        res["channels"][c] = {"r2_vs_mean": r2, "vs_random": vr, "vs_deranged": vd}
        print(f"{c:24s} {r2['point']:+.3f}[{r2['ci95'][0]:+.3f},{r2['ci95'][1]:+.3f}]"
              f"   {vr['point']:+.4f}[{vr['ci95'][0]:+.4f},{vr['ci95'][1]:+.4f}]"
              f"   {vd['point']:+.4f}[{vd['ci95'][0]:+.4f},{vd['ci95'][1]:+.4f}]")
    json.dump(res, open(os.path.join(OUT, "S4_OBSERVABILITY_AUDIT.json"), "w"),
              indent=2, default=float)
    print("\nwrote S4_OBSERVABILITY_AUDIT.json")
    return res


if __name__ == "__main__":
    main()
