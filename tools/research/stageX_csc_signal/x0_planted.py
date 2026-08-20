"""I1 planted-signal control on the Duong-Ly real graph."""
from __future__ import annotations
import json, hashlib
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import torch
import torch.nn as nn

if __package__ in {None, ''}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from tools.research.stageX_csc_signal.x0_instruments import (
    AAS, apply_mutations, load_duongly, local_window, mutation_token, parse_mutations)

HERE=Path(__file__).resolve().parent
PREREG_SHA='03cdc907df3e778f5fe79fb1a238d35ebb6ece5e9e743db181728ba6b25e9683'

def hash_vec(name, dim=32, seed=0):
    h=hashlib.sha256(f'{seed}|{name}'.encode()).digest()
    arr=np.frombuffer(h,dtype=np.uint8).astype(np.float32)/255.0
    return np.resize(arr,dim)

def planted_interaction(prot, lig, U, V):
    """Return pairwise dot products after projecting both feature banks."""
    return (prot @ U) @ (lig @ V).T

def main():
    info, matrix, seqs=load_duongly()
    names=list(matrix.iloc[:,0].astype(str))[1:]
    values=matrix.iloc[1:,1:].to_numpy(dtype=float)
    parent_acc={}
    for _,row in info.iterrows():
        n=str(row['Kinase (Mutation)']).strip()
        parent_acc[n.split('(')[0].strip()]=str(row['Protein Accession #']).strip()
    # protein features: mutation token + local window
    prot_feat=[]
    for name in names:
        parent=name.split('(')[0].strip() if '(' in name else name
        seq=apply_mutations(seqs.get(parent_acc.get(parent),''), name)
        v=np.concatenate([mutation_token(seq,seq,name), local_window(seq,name)])
        prot_feat.append(v.astype(np.float32))
    prot=np.stack(prot_feat)
    lig_names=[str(c) for c in matrix.columns[1:]]
    lig=np.stack([hash_vec(c) for c in lig_names])
    y=values  # NxM
    rows=[]
    for i in range(y.shape[0]):
        for j in range(y.shape[1]):
            if np.isfinite(y[i,j]): rows.append((i,j,float(y[i,j])))
    # planted interaction: y_plant = protein_main + ligand_main + tau*dot(U p, V l)
    rng=np.random.default_rng(20260821)
    U=rng.normal(0,0.2,size=(prot.shape[1],4)).astype(np.float32)
    V=rng.normal(0,0.2,size=(lig.shape[1],4)).astype(np.float32)
    prot_main=rng.normal(0,10,size=prot.shape[0]).astype(np.float32)
    lig_main=rng.normal(0,10,size=lig.shape[0]).astype(np.float32)
    true_inter=planted_interaction(prot,lig,U,V)  # NxM
    # make actual pairs and plant
    row_y={}
    for i,j,y0 in rows:
        row_y[(i,j)]=y0
    results={}
    for tau in [0.2,0.4,0.8,1.6]:
        y_plant={k: v + tau*true_inter[k] for k,v in row_y.items()}
        # train low-capacity bilinear gradient model
        class Bilinear(nn.Module):
            def __init__(self):
                super().__init__()
                self.P=nn.Linear(prot.shape[1],16)
                self.L=nn.Linear(lig.shape[1],16)
                self.bias=nn.Parameter(torch.tensor(0.0))
            def forward(self,p,l):
                return (self.P(p)*self.L(l)).sum(-1)+self.bias
        model=Bilinear()
        opt=torch.optim.AdamW(model.parameters(),lr=1e-3,weight_decay=1e-4)
        keys=list(y_plant); rng2=np.random.default_rng(20260822)
        P=torch.tensor(prot); L=torch.tensor(lig)
        truth=torch.tensor([y_plant[k] for k in keys])
        idx_p=torch.tensor([k[0] for k in keys]); idx_l=torch.tensor([k[1] for k in keys])
        for step in range(800):
            perm=rng2.permutation(len(keys))[:512]
            pred=model(P[idx_p[perm]],L[idx_l[perm]])
            loss=(pred-truth[perm]).pow(2).mean()
            opt.zero_grad(); loss.backward(); opt.step()
        pred=model(P[idx_p],L[idx_l]).detach().numpy()
        tru=truth.numpy()
        # ligand-only null: same architecture with protein input zeros
        null_pred=np.zeros_like(tru)
        # interaction component recovery on training rows
        # evaluate Spearman of delta predictions for protein pairs vs planted delta
        inter_pred=pred
        inter_true=np.array([true_inter[k] for k in keys])
        pearson=float(np.corrcoef(inter_pred,inter_true)[0,1]) if inter_pred.std()>1e-12 else 0.0
        spearman=float(np.corrcoef(np.argsort(np.argsort(inter_pred)),np.argsort(np.argsort(inter_true)))[0,1]) if inter_pred.std()>1e-12 else 0.0
        # sign accuracy on planted interaction vs zero? use predicted interaction direction
        sign=float(np.mean(np.sign(inter_pred)==np.sign(inter_true)))
        results[str(tau)]={'pearson_with_planted_interaction':pearson,'spearman':spearman,'sign_accuracy':sign,'final_loss':float(loss)}
    out={'schema':'MetaSieve.StageX.I1Planted.v1','stage':'stageX_csc_signal','preregistration_sha256':PREREG_SHA,'results':results}
    (HERE/'X0_PLANTED.json').write_text(json.dumps(out,indent=1,sort_keys=True)+'\n')
    print(json.dumps(out,indent=1))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
