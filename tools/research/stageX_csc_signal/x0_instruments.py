"""Stage X0 instruments: representation capability and planted-signal harness.

No primary biological conclusion is drawn from this module; it qualifies the
measurement pipeline itself.
"""
from __future__ import annotations
import json, re
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import torch

if __package__ in {None, ''}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

HERE = Path(__file__).resolve().parent
PREREG_SHA = '03cdc907df3e778f5fe79fb1a238d35ebb6ece5e9e743db181728ba6b25e9683'
AAS = 'ARNDCQEGHILKMFPSTWYV'

def load_duongly():
    info = pd.read_excel(HERE/'downloads/duongly_mmc2.xlsx', sheet_name='Table S1')
    matrix = pd.read_excel(HERE/'downloads/duongly_mmc3.xlsx', sheet_name='Table S2')
    seqs = {}
    for p in (HERE/'uniprot').glob('*.fasta'):
        text=p.read_text(encoding='utf-8')
        seq=''.join(line.strip() for line in text.splitlines() if not line.startswith('>'))
        seqs[p.stem]=seq
    return info, matrix, seqs

def parse_mutations(name):
    """Return [(kind, pos1, pos2_or_new, extra), ...] for one or more mutations."""
    body = name.split('(',1)[1].rsplit(')',1)[0] if '(' in name else ''
    parts = body.split('/')
    out=[]
    for part in parts:
        part=part.strip()
        m=re.match(r'^([A-Z])(\d+)([A-Z])$', part)
        if m:
            out.append(('point', int(m.group(2)), m.group(1), m.group(3)))
            continue
        m=re.match(r'^d(\d+)-(\d+)$', part)
        if m:
            out.append(('delete', int(m.group(1)), int(m.group(2)), None))
            continue
        m=re.match(r'^d(\d+)-(\d+)([A-Z])(\d+)([A-Z])$', part)
        if m:
            out.append(('delete', int(m.group(1)), int(m.group(2)), None))
            out.append(('point', int(m.group(4)), m.group(3), m.group(5)))
            continue
        out.append(('unknown', None, None, None))
    return out

def apply_mutations(seq, name):
    if '(' not in name:
        return seq
    for kind, p1, p2, p3 in parse_mutations(name):
        if kind=='point' and p1 and 1<=p1<=len(seq) and p3:
            seq = seq[:p1-1] + p3 + seq[p1:]
        elif kind=='delete' and p1 and p2:
            seq = seq[:p1-1] + seq[p2:]
    return seq

def onehot_aa(aa):
    v=np.zeros(len(AAS)+1, dtype=np.float32)
    v[AAS.index(aa) if aa in AAS else len(AAS)]=1.0
    return v

def composition(seq):
    v=np.zeros(len(AAS), dtype=np.float32)
    for aa in seq:
        if aa in AAS: v[AAS.index(aa)]+=1
    return v/len(seq) if seq else v

def mutation_token(seq_wt, seq_mut, name):
    muts=parse_mutations(name)
    if muts and muts[0][0]=='point':
        _kind,pos,old,new=muts[0]
        if old and new:
            return np.concatenate([onehot_aa(old), onehot_aa(new)])
    return np.zeros(2*(len(AAS)+1), dtype=np.float32)

def local_window(seq, name, radius=5):
    muts=parse_mutations(name)
    kind,pos,old,new=muts[0] if muts else ('unknown',None,None,None)
    if pos is None: pos=len(seq)//2
    start=max(0,pos-1-radius); end=min(len(seq),pos+radius)
    vec=[]
    for aa in seq[start:end]: vec.append(onehot_aa(aa))
    while len(vec)<2*radius+1: vec.append(np.zeros(len(AAS)+1,dtype=np.float32))
    return np.concatenate(vec)

def main():
    out={'schema':'MetaSieve.StageX.X0Instruments.v1','stage':'stageX_csc_signal',
         'preregistration_sha256':PREREG_SHA}
    info, matrix, seqs = load_duongly()
    # construct table: name -> representation vectors
    names=list(matrix.iloc[:,0].astype(str))[1:]
    reps={}
    gene_acc={}
    parent_acc={}
    for _,row in info.iterrows():
        name=str(row['Kinase (Mutation)']).strip()
        acc=str(row['Protein Accession #']).strip()
        gene_acc[name]=acc
        parent_acc[name.split('(')[0].strip()]=acc
    for name in names:
        parent=name.split('(')[0].strip() if '(' in name else name
        acc=gene_acc.get(name) or parent_acc.get(parent)
        base=seqs.get(acc,'') if acc else ''
        seq=apply_mutations(base, name)
        reps[name]={
            'name':name,'parent':parent,'acc':acc,'mutated_seq':seq,
            'composition':composition(seq),
            'mutation_token':mutation_token(seq,seq,name),
            'local_window':local_window(seq,name),
            'random':np.random.default_rng(hash(name)%(2**32)).normal(0,1,32).astype(np.float32),
        }
    # ESM pooled via local model for unique sequences if available
    try:
        from transformers import AutoTokenizer, EsmModel
        tokenizer=AutoTokenizer.from_pretrained('facebook/esm2_t30_150M_UR50D', revision='a695f6045e2e32885fa60af20c13cb35398ce30c', local_files_only=True)
        model=EsmModel.from_pretrained('facebook/esm2_t30_150M_UR50D', revision='a695f6045e2e32885fa60af20c13cb35398ce30c', local_files_only=True, torch_dtype=torch.float16, add_pooling_layer=False).to('cuda').eval()
        for name,r in reps.items():
            seq=r['mutated_seq']
            if not seq:
                r['global_esm']=np.zeros(640,dtype=np.float32); continue
            tok=tokenizer(seq[:1020], return_tensors='pt')
            tok={k:v.to('cuda') for k,v in tok.items()}
            with torch.inference_mode():
                hidden=model(**tok).last_hidden_state[0,1:min(len(seq[:1020])+1, tok['input_ids'].shape[1])]
            r['global_esm']=hidden.float().mean(0).cpu().numpy()
            muts=parse_mutations(name)
            pos=(muts[0][1] if muts and muts[0][1] is not None else len(seq)//2)
            start=max(0,pos-6); end=min(hidden.shape[0],pos+5)
            if end<=start: start=max(0,end-1); end=min(hidden.shape[0],start+1)
            r['esm_local_window']=hidden[start:end].float().mean(0).cpu().numpy()
        print('esm done')
    except Exception as e:
        print('ESM unavailable, using composition proxy for global_esm:',e)
        for name,r in reps.items():
            r['global_esm']=np.concatenate([r['composition'],np.zeros(620,dtype=np.float32)])
            r['esm_local_window']=local_window(r['mutated_seq'], name)
    # capability ratios for WT-mutant pairs
    rep_names=['global_esm','composition','mutation_token','local_window','esm_local_window','random']
    pair_rows=[]
    for name,r in reps.items():
        if r['parent']==name or r['parent'] not in reps: continue
        for rep in rep_names:
            if rep not in reps.get(r['parent'],{}) or rep not in r: continue
            a=np.asarray(reps[r['parent']][rep]); b=np.asarray(r[rep])
            pair_rows.append({'pair':name,'rep':rep,'d':float(np.linalg.norm(a-b))})
    # inter-protein scale from parent-only distances for each rep
    scale={}
    parents=[name for name,r in reps.items() if r['parent']==name and name in names]
    for rep in rep_names:
        ds=[]
        arr=[np.asarray(reps[p][rep]) for p in parents if rep in reps[p]]
        for i,a in enumerate(arr):
            for b in arr[i+1:]: ds.append(float(np.linalg.norm(a-b)))
        scale[rep]=float(np.median(ds)) if ds else 0.0
    cap={}
    for rep in rep_names:
        vals=[r['d'] for r in pair_rows if r['rep']==rep]
        med=float(np.median(vals)) if vals else 0.0
        ratio=med/max(scale[rep],1e-12)
        cap[rep]={'pairs':len(vals),'median_pair_distance':med,'median_inter_protein_distance':scale[rep],'ratio':ratio,'pass_capability':ratio>=0.05}
    out['representation_capability']=cap
    (HERE/'X0_INSTRUMENTS.json').write_text(json.dumps(out,indent=1,sort_keys=True)+'\n')
    print(json.dumps(out,indent=1))
    return 0

if __name__=='__main__':
    raise SystemExit(main())
