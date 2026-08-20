import sys
from pathlib import Path
import numpy as np, torch
HERE = Path('.').resolve(); sys.path.insert(0, str(HERE))
import olr as O, runner as R

d1, z1, d2, states, recs = R.build_all()
lig = z1['lig']
split = np.array([r['split1'] for r in recs], dtype=np.int8)
train_js = [j for j in range(len(recs)) if split[j]==0]
test_js = [j for j in range(len(recs)) if split[j]==2]
rng = np.random.default_rng(20260821)
U = rng.normal(0,1,(O.D_RES,4)).astype(np.float32); V = rng.normal(0,1,(O.D_LIG,4)).astype(np.float32)
zL = lig.astype(np.float32) @ V
hbar = {}
for j,r in enumerate(recs): hbar.setdefault(r['parent'], r['res_w'].mean(0).numpy())
grand = np.mean([h for h in hbar.values()], axis=0)
fields = {p: ((h-grand)@U)@zL.T for p,h in hbar.items()}
cv = float(np.var(np.stack(list(fields.values())), axis=0).mean())
fields = {p: f*np.sqrt(134.8/cv) for p,f in fields.items()}
def scale(x, var):
    v = x.var(); return x*np.sqrt(var/v)
s0 = scale((zL[:,0]+zL[:,1])*0.5, 50.0)
for j,r in enumerate(recs):
    u = rng.normal(0,1,4).astype(np.float32)
    prof = (s0 + fields[r['parent']] + scale(u@zL.T,44.85) + rng.normal(0,np.sqrt(44.85),183).astype(np.float32))[r['lig_idx']]
    r['c'] = (prof - prof.mean()).astype(np.float32)
res_pad, mask_pad = R.pad_states(recs)
res_pad, mask_pad = res_pad.to(R.DEVICE), mask_pad.to(R.DEVICE)
prior = O.ligand_prior(recs, train_js, test_js)
tc = np.concatenate([recs[j]['c'] for j in test_js])
pc = np.concatenate([prior[j] for j in test_js])
print('A0 prior R2:', round(O.r2_cells(pc, tc),4))
for arm in ['A1-bilinear','A2-router','A5-gain']:
    parents = [r['parent'] for r in recs]
    fold_of = O.folds_by_parent(parents)
    w,_ = O.gain_weights(d1, z1, recs, train_js)
    mhat = R.crossfit_nuisance(recs, lig, train_js, fold_of, 99) if arm=='A5-gain' else None
    m = mhat if arm=='A5-gain' else None
    ww = w if arm=='A5-gain' else None
    run = R.train_arm(arm, recs, lig, res_pad, mask_pad, split, seed=11, mhat=m, weights=ww)
    pc = np.concatenate([run['preds'][j] for j in test_js])
    pt = np.concatenate([run['preds'][j] for j in train_js])
    tt = np.concatenate([recs[j]['c'] for j in train_js])
    print(arm, 'test R2', round(O.r2_cells(pc, tc),4), 'train R2', round(O.r2_cells(pt, tt),4), 'var_rec', round(O.var_recovery(pc,tc),3), flush=True)
