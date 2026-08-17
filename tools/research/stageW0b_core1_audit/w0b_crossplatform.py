"""Stage W0b cross-platform residual reproducibility (read-only audit).

Interpretation per the current instruction: this is ONLY a cross-platform
transfer gate. Low cross-platform reproducibility closes direct cross-platform
residual sharing; it cannot kill single-platform protein-conditioned signal.
"""
from __future__ import annotations
import argparse, csv, json
from collections import defaultdict
from pathlib import Path
import sys
import numpy as np
import pandas as pd
if __package__ in {None, ''}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
from tools.research.stageW0b_core1_audit.w0b_audit import PREREG_SHA
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]

def canon(smiles):
    if not smiles:
        return None
    m=Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(m) if m is not None else None

def load_davis():
    cells=defaultdict(dict)  # kinase -> {smiles: value}
    with (ROOT/'dataset/raw/dta/davis.tab').open(newline='',encoding='utf-8') as f:
        r=csv.reader(f,delimiter='\t'); next(r)
        for line in r:
            if len(line)<5: continue
            kinase=line[2].upper(); smi=canon(line[1]); val=float(line[4])
            if kinase and smi and val<10000.0-1e-9:
                cells[kinase][smi]=val
    return cells

def load_metz():
    xls=pd.read_excel(ROOT/'dataset/raw/crossed_panels/kinase_panels/metz.xls',sheet_name=0)
    smi={}
    for _,row in xls.iterrows():
        if pd.notna(row.get('Cmpd_ID')) and pd.notna(row.get('Canonical_Smiles')):
            smi[int(row['Cmpd_ID'])]=canon(str(row['Canonical_Smiles']))
    matrix=pd.read_csv(ROOT/'dataset/raw/crossed_panels/kinase_panels/metz_matrix.csv')
    matrix=matrix.rename(columns={matrix.columns[0]:'Cmpd_ID'})
    cells=defaultdict(dict)
    for _,line in matrix.iterrows():
        cid=int(line['Cmpd_ID']); s=smi.get(cid)
        if not s: continue
        for col in matrix.columns[1:]:
            v=line[col]
            if pd.isna(v): continue
            v=float(v)
            if v<=4.0+1e-9: continue
            cells[str(col).upper()][s]=v
    return cells

def load_klaeger():
    smi=json.loads((ROOT/'dataset/processed/crossed_panels_xp2/klaeger_smiles.json').read_text())
    matrix=pd.read_csv(ROOT/'dataset/raw/crossed_panels/kinase_panels/klaeger_matrix.csv')
    cells=defaultdict(dict)
    for _,line in matrix.iterrows():
        entry=smi.get(str(line['Drug'])) or {}
        s=canon(entry.get('smiles'))
        if not s: continue
        for col in matrix.columns[1:]:
            v=line[col]
            if pd.isna(v): continue
            v=float(v)
            if v<=5.0+1e-9: continue
            cells[str(col).upper()][s]=v
    return cells

def pair_deltas(cells):
    out=defaultdict(dict)  # kinase -> {(smi_a,smi_b): delta}
    for kinase,ligands in cells.items():
        if len(ligands)<2: continue
        ordered=sorted(ligands)
        for i,a in enumerate(ordered):
            for b in ordered[i+1:]:
                out[kinase][(a,b)]=ligands[b]-ligands[a]
    return out

def compare(left,right,left_name,right_name):
    keys=[]
    for kinase in set(left)&set(right):
        for pair in set(left[kinase])&set(right[kinase]):
            keys.append((kinase,pair,left[kinase][pair],right[kinase][pair]))
    if len(keys)<20:
        return {'n':len(keys),'identifiable':False}
    lx=np.array([k[2] for k in keys]); rx=np.array([k[3] for k in keys])
    def spearman(a,b):
        return float(np.corrcoef(np.argsort(np.argsort(a)),np.argsort(np.argsort(b)))[0,1])
    point=spearman(lx,rx)
    kinases=sorted({k[0] for k in keys})
    rng=np.random.default_rng(20260820)
    vals=[]
    by=defaultdict(list)
    for k in keys: by[k[0]].append(k)
    for _ in range(1000):
        idx=rng.integers(0,len(kinases),size=len(kinases))
        subset=[]
        for i in idx: subset.extend(by[kinases[i]])
        if len(subset)<20: continue
        a=np.array([s[2] for s in subset]); b=np.array([s[3] for s in subset])
        vals.append(spearman(a,b))
    vals=np.array(vals)
    return {'n':len(keys),'kinases':len(kinases),'identifiable':True,
            'spearman':point,'ci_lo':float(np.quantile(vals,.025)),
            'ci_hi':float(np.quantile(vals,.975)),
            'interpretation':'cross-platform transfer gate only'}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,default=HERE/'W0B_CROSSPLATFORM.json'); args=ap.parse_args()
    davis=load_davis(); metz=load_metz(); klaeger=load_klaeger()
    d_delta=pair_deltas(davis); m_delta=pair_deltas(metz); k_delta=pair_deltas(klaeger)
    out={'schema':'MetaSieve.StageW0b.CrossPlatform.v1','stage':'stageW0b_core1_audit',
         'preregistration_sha256':PREREG_SHA,
         'censoring_rule':'detection-floor rows excluded before pair formation',
         'platforms':{
          'metz_vs_klaeger':compare(m_delta,k_delta,'metz','klaeger'),
          'metz_vs_davis':compare(m_delta,d_delta,'metz','davis'),
          'klaeger_vs_davis':compare(k_delta,d_delta,'klaeger','davis'),
         }}
    args.output.write_text(json.dumps(out,indent=1,sort_keys=True)+'\n')
    print(json.dumps(out,indent=1))
    return 0
if __name__=='__main__':
    raise SystemExit(main())
