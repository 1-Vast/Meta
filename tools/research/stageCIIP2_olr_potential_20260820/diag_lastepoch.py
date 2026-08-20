import sys
from pathlib import Path
import numpy as np, torch
HERE = Path('.').resolve(); sys.path.insert(0, str(HERE))
import olr as O, runner as R
import numpy as np

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
device = R.DEVICE
res_pad, mask_pad = res_pad.to(device), mask_pad.to(device)
prior = O.ligand_prior(recs, train_js, test_js)
tc = np.concatenate([recs[j]['c'] for j in test_js])
pc0 = np.concatenate([prior[j] for j in test_js])
parents = [r['parent'] for r in recs]
fold_of = O.folds_by_parent(parents)
w,_ = O.gain_weights(d1, z1, recs, train_js)
mhat = R.crossfit_nuisance(recs, lig, train_js, fold_of, 99)
tgt = {j: recs[j]['c'] - mhat[j] for j in range(len(recs))}

torch.manual_seed(11)
model = O.OLRPotential(router=True).to(device)
opt = torch.optim.AdamW(model.parameters(), lr=3e-3, weight_decay=1e-4)
sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=900)
ligT = torch.from_numpy(lig).float().to(device)
w_full = torch.from_numpy(w).float().to(device)
for ep in range(900):
    opt.zero_grad()
    kv = {}
    for j in range(len(recs)):
        kv[(j,0)] = model.construct_kv(res_pad[j,0], mask_pad[j,0])
        kv[(j,1)] = model.construct_kv(res_pad[j,1], mask_pad[j,1])
    total = 0.
    for j in train_js:
        r = recs[j]
        L = ligT[torch.tensor(r['lig_idx'], device=device)]
        sw = model.s_from_kv(kv[(j,0)], mask_pad[j,0], L)
        sv = model.s_from_kv(kv[(j,1)], mask_pad[j,1], L)
        chat = sv - sw; chat = chat - chat.mean()
        t = torch.from_numpy(tgt[j].astype(np.float32)).to(device)
        ww = w_full[torch.tensor(r['lig_idx'], device=device)]
        total = total + (ww*(t-chat)**2).mean()
    (total/len(train_js)).backward(); opt.step(); sched.step()
with torch.no_grad():
    preds = {}
    for j in test_js:
        r = recs[j]
        L = ligT[torch.tensor(r['lig_idx'], device=device)]
        kw = model.construct_kv(res_pad[j,0], mask_pad[j,0])
        km = model.construct_kv(res_pad[j,1], mask_pad[j,1])
        cs = model.s_from_kv(km, mask_pad[j,1], L) - model.s_from_kv(kw, mask_pad[j,0], L)
        cs = cs - cs.mean()
        preds[j] = cs.cpu().numpy() + mhat[j]
pc = np.concatenate([preds[j] for j in test_js])
print('last-epoch test R2:', round(O.r2_cells(pc, tc),4), 'delta vs A0:', round(O.r2_cells(pc, tc)-O.r2_cells(pc0, tc),4))
