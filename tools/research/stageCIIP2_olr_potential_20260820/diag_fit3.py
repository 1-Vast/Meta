import sys
from pathlib import Path
import numpy as np, torch, json
HERE = Path('.').resolve(); sys.path.insert(0, str(HERE))
import olr as O, runner as R

d1, z1, d2, states, recs = R.build_all()
lig = z1['lig']
split = np.array([r['split1'] for r in recs], dtype=np.int8)
train_js = [j for j in range(len(recs)) if split[j]==0]
rng = np.random.default_rng(20260821)
U = rng.normal(0,1,(O.D_RES,4)).astype(np.float32); V = rng.normal(0,1,(O.D_LIG,4)).astype(np.float32)
zL = lig.astype(np.float32) @ V
hbar = {}
for j,r in enumerate(recs): hbar.setdefault(r['parent'], r['res_w'].mean(0).numpy())
def scale(x, var):
    v = x.var(); return x*np.sqrt(var/v)
g_par = {p: scale((h@U)@zL.T, 134.8) for p,h in hbar.items()}
s0 = scale((zL[:,0]+zL[:,1])*0.5, 50.0)
for j,r in enumerate(recs):
    u = rng.normal(0,1,4).astype(np.float32)
    prof = (s0 + g_par[r['parent']] + scale(u@zL.T,44.85) + rng.normal(0,np.sqrt(44.85),183).astype(np.float32))[r['lig_idx']]
    r['c'] = (prof - prof.mean()).astype(np.float32)
res_pad, mask_pad = R.pad_states(recs)
w, _ = O.gain_weights(d1, z1, recs, train_js)
parents = [r['parent'] for r in recs]
fold_of = O.folds_by_parent(parents)
mhat = R.crossfit_nuisance(recs, lig, train_js, fold_of, 20260821)
# residual variance check
rr = np.concatenate([recs[j]['c'] - mhat[j] for j in range(len(recs))])
print('residual total var:', round(float(rr.var()),1))
run = R.train_arm('A5-gain', recs, lig, res_pad.to(R.DEVICE), mask_pad.to(R.DEVICE), split, seed=11, mhat=mhat, weights=w, log=True)
tgt_res = R.targets_for(recs, range(len(recs)), 'A5-gain', mhat)
test_js = run['test_js']
pc = np.concatenate([run['preds'][j] for j in test_js]); tc = np.concatenate([tgt_res[j] for j in test_js])
print('test residual R2:', round(O.r2_cells(pc, tc),4))
pt = np.concatenate([run['preds'][j] for j in train_js]); tt = np.concatenate([tgt_res[j] for j in train_js])
print('TRAIN residual R2:', round(O.r2_cells(pt, tt),4), 'var_rec', round(float(np.var(pt)/np.var(tt)),4))
