"""XP3 STAGE 1b — label-blind census of the release-pinned BindingDB articles subset.

Counts design structure only.  Affinity VALUES are never selected: a Ki/Kd column
is used solely through `notna()` to count usable cells, exactly as E-AFF-X0 did.
"""
from __future__ import annotations

import io
import json
import os
import zipfile

import numpy as np
import pandas as pd

RAW = r"D:\MetaSieve\dataset\raw\crossed_panels\bindingdb"
OUT = r"D:\MetaSieve\report\observability_census"
os.makedirs(OUT, exist_ok=True)
ZIP = os.path.join(RAW, "BindingDB_BindingDB_Articles_202608_tsv.zip")

with zipfile.ZipFile(ZIP) as z:
    name = [n for n in z.namelist() if n.endswith(".tsv")][0]
    print("member:", name, z.getinfo(name).file_size, "bytes uncompressed")
    with z.open(name) as f:
        head = pd.read_csv(io.TextIOWrapper(f, "utf-8", errors="ignore"),
                           sep="\t", nrows=5, low_memory=False, on_bad_lines="skip")
cols = list(head.columns)
print("columns:", len(cols))
KEEP = [c for c in cols if any(k in c for k in (
    "BindingDB Reagent Set ID", "Ligand SMILES", "Ki (nM)", "Kd (nM)", "IC50 (nM)",
    "PMID", "Article DOI", "Target Name", "Target Source Organism",
    "BindingDB Target Chain Sequence", "UniProt (SwissProt) Primary ID",
    "Curation/DataSource", "Assay Description", "Temperature", "pH",
    "Number of Protein Chains in Target"))]
print("kept columns:", KEEP)

with zipfile.ZipFile(ZIP) as z:
    with z.open(name) as f:
        df = pd.read_csv(io.TextIOWrapper(f, "utf-8", errors="ignore"), sep="\t",
                         usecols=[c for c in KEEP], low_memory=False,
                         on_bad_lines="skip")
print("rows:", len(df))

SMI = next(c for c in df.columns if "Ligand SMILES" in c)
SEQ = next((c for c in df.columns if "Target Chain Sequence" in c), None)
UNI = next((c for c in df.columns if "UniProt" in c), None)
PMID = next((c for c in df.columns if "PMID" in c), None)
KI = next((c for c in df.columns if c.startswith("Ki (nM)")), None)
KD = next((c for c in df.columns if c.startswith("Kd (nM)")), None)

report = {"stage": "XP3-1b", "release": "BindingDB_BindingDB_Articles_202608",
          "license": "CC BY 3.0", "affinity_values_read": 0, "rows": int(len(df)),
          "endpoints": {}}

for label, col in (("Ki", KI), ("Kd", KD)):
    if col is None:
        continue
    has = df[col].notna()          # PRESENCE only, never the value
    sub = df.loc[has, [c for c in (SMI, UNI, SEQ, PMID) if c]].copy()
    sub = sub.dropna(subset=[SMI] + ([UNI] if UNI else []))
    n_lig = sub[SMI].nunique()
    n_tgt = sub[UNI].nunique() if UNI else None
    cells = sub.groupby([UNI, SMI]).size() if UNI else None
    # crossed structure: ligands seen on >=2 targets, targets with >=20 ligands
    lig_t = sub.groupby(SMI)[UNI].nunique()
    tgt_l = sub.groupby(UNI)[SMI].nunique()
    # document closure: union targets sharing a PMID (the X0/D1 construction)
    comp = {}
    if PMID:
        lab = {t: t for t in sub[UNI].unique()}

        def find(x):
            while lab[x] != x:
                lab[x] = lab[lab[x]]
                x = lab[x]
            return x

        for _, g in sub.dropna(subset=[PMID]).groupby(PMID)[UNI]:
            ts = list(dict.fromkeys(g))
            for t in ts[1:]:
                a, b = find(ts[0]), find(t)
                if a != b:
                    lab[b] = a
        comp = {t: find(t) for t in lab}
    n_comp = len({v for v in comp.values()}) if comp else None
    # panels: documents reporting >=2 targets AND >=5 shared ligands
    panels = 0
    crossed_rect = 0
    if PMID:
        for pm, g in sub.dropna(subset=[PMID]).groupby(PMID):
            nt, nl = g[UNI].nunique(), g[SMI].nunique()
            if nt >= 2 and nl >= 5:
                panels += 1
                crossed_rect += (nt * (nt - 1) // 2) * (nl * (nl - 1) // 2)
    report["endpoints"][label] = {
        "rows_with_value_present": int(has.sum()),
        "distinct_ligands": int(n_lig), "distinct_targets": int(n_tgt or 0),
        "distinct_cells": int(len(cells)) if cells is not None else None,
        "ligands_on_ge2_targets": int((lig_t >= 2).sum()),
        "ligands_on_ge5_targets": int((lig_t >= 5).sum()),
        "targets_with_ge20_ligands": int((tgt_l >= 20).sum()),
        "documents": int(sub[PMID].nunique()) if PMID else None,
        "document_closure_components": n_comp,
        "multi_target_documents_with_ge5_ligands": panels,
        "nominal_crossed_rectangles": int(crossed_rect),
        "sequences_available": bool(SEQ and sub[SEQ].notna().mean() > 0.5)
        if SEQ else False,
    }
    print(f"\n--- {label} ---")
    for k, v in report["endpoints"][label].items():
        print(f"   {k:42s} {v}")

json.dump(report, open(os.path.join(OUT, "XP3_BINDINGDB_CENSUS.json"), "w"),
          indent=2, default=float)
print("\nwrote XP3_BINDINGDB_CENSUS.json")
