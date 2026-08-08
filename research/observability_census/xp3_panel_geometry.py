"""XP3 STAGE 1c — panel geometry of the BindingDB curated-articles Ki subset.

Still label-blind: Ki is used only through notna().  Determines whether the
many-small-panels design can answer the k<=5 double-closure estimand, and checks
for contamination by panels already consumed (Metz, Klaeger, PDSP, DAVIS,
PKIS2, Anastassiadis).
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
ZIP = os.path.join(RAW, "BindingDB_BindingDB_Articles_202608_tsv.zip")

# PMIDs of panels already consumed by this project, plus prohibited sources.
CONSUMED_PMID = {
    "21572424": "Metz 2011 (XP1/XP2 primary)",
    "29191878": "Klaeger 2017 (XP1 cross-platform, XP2-F)",
    "22037378": "Davis 2011 (PROHIBITED)",
    "18183025": "Karaman 2008",
    "21949673": "Anastassiadis 2011 (consumed)",
    "28767711": "Drewry PKIS2 (consumed)",
}

USE = ["Ligand SMILES", "Target Name", "Ki (nM)", "PMID", "Article DOI",
       "Curation/DataSource", "Target Source Organism According to Curator or DataSource",
       "UniProt (SwissProt) Primary ID of Target Chain 1",
       "BindingDB Target Chain Sequence 1",
       "Number of Protein Chains in Target (>1 implies a multichain complex)"]

with zipfile.ZipFile(ZIP) as z:
    name = [n for n in z.namelist() if n.endswith(".tsv")][0]
    with z.open(name) as f:
        df = pd.read_csv(io.TextIOWrapper(f, "utf-8", errors="ignore"), sep="\t",
                         usecols=USE, low_memory=False, on_bad_lines="skip")

SMI, UNI, PMID = "Ligand SMILES", "UniProt (SwissProt) Primary ID of Target Chain 1", "PMID"
SEQ = "BindingDB Target Chain Sequence 1"
NCH = "Number of Protein Chains in Target (>1 implies a multichain complex)"
ORG = "Target Source Organism According to Curator or DataSource"

d = df[df["Ki (nM)"].notna()].copy()                    # PRESENCE only
d[PMID] = d[PMID].astype("string").str.strip()
d = d.dropna(subset=[SMI, UNI, PMID, SEQ])
d = d[d[NCH].fillna(1).astype(float) <= 1]              # single-chain targets only
print("single-chain Ki rows with sequence, structure and PMID:", len(d))

contam = {p: (n, int((d[PMID] == p).sum())) for p, n in CONSUMED_PMID.items()}
print("\ncontamination check (rows in this subset):")
for p, (n, c) in contam.items():
    print(f"   PMID {p:10s} {n:44s} rows={c}")
d = d[~d[PMID].isin(CONSUMED_PMID)]
print("rows after excluding consumed/prohibited PMIDs:", len(d))

d = d[d[ORG].astype(str).str.contains("Homo sapiens", na=False)]
print("rows restricted to Homo sapiens:", len(d))

cell = d.groupby([PMID, UNI, SMI]).size().reset_index(name="n")
geo = []
for pm, g in cell.groupby(PMID):
    nt, nl = g[UNI].nunique(), g[SMI].nunique()
    piv = g.pivot_table(index=SMI, columns=UNI, values="n", aggfunc="size")
    dens = float(piv.notna().to_numpy().mean()) if piv.size else 0.0
    geo.append({"pmid": pm, "targets": int(nt), "ligands": int(nl),
                "cells": int(len(g)), "density": dens})
G = pd.DataFrame(geo)
print(f"\npanels (documents): {len(G)}")

report = {"stage": "XP3-1c", "release": "BindingDB_BindingDB_Articles_202608",
          "affinity_values_read": 0,
          "rows_used": int(len(d)), "panels": int(len(G)),
          "contamination_excluded": {p: contam[p][1] for p in CONSUMED_PMID},
          "distinct_targets": int(d[UNI].nunique()),
          "distinct_ligands": int(d[SMI].nunique())}

print(f"{'requirement':44s} {'panels':>8s} {'targets':>8s} {'ligands':>8s} {'cells':>8s}")
for tmin, lmin in ((2, 5), (2, 10), (3, 10), (2, 20), (3, 20), (4, 20), (3, 30), (5, 30)):
    sel = G[(G.targets >= tmin) & (G.ligands >= lmin)]
    key = f"targets>={tmin}, ligands>={lmin}"
    row = {"panels": int(len(sel)), "targets": int(sel.targets.sum()),
           "ligands": int(sel.ligands.sum()), "cells": int(sel.cells.sum()),
           "median_density": float(sel.density.median()) if len(sel) else None}
    report.setdefault("design_requirements", {})[key] = row
    print(f"{key:44s} {row['panels']:8d} {row['targets']:8d} {row['ligands']:8d} "
          f"{row['cells']:8d}")

# the working design: >=2 targets and enough ligands for k=5 support + test
work = G[(G.targets >= 2) & (G.ligands >= 20)]
report["working_design"] = {
    "definition": "documents with >=2 single-chain human targets and >=20 distinct "
                  "ligands measured with Ki; consumed/prohibited PMIDs removed",
    "panels": int(len(work)), "total_cells": int(work.cells.sum()),
    "median_targets_per_panel": float(work.targets.median()),
    "median_ligands_per_panel": float(work.ligands.median()),
    "median_density": float(work.density.median()),
    "max_targets": int(work.targets.max()) if len(work) else 0,
    "max_ligands": int(work.ligands.max()) if len(work) else 0,
}
print("\nworking design:", json.dumps(report["working_design"], indent=2))
G.sort_values("cells", ascending=False).head(15).to_csv(
    os.path.join(OUT, "xp3_top_panels.csv"), index=False)
json.dump(report, open(os.path.join(OUT, "XP3_PANEL_GEOMETRY.json"), "w"),
          indent=2, default=float)
print("wrote XP3_PANEL_GEOMETRY.json")
