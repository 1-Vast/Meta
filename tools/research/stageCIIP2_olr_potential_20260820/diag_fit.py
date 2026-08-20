import sys, json
from pathlib import Path
import numpy as np, torch
HERE = Path('.').resolve()
sys.path.insert(0, str(HERE))
import olr as O, runner as R

d1, z1, d2, states, recs = R.build_all()
res_pad, mask_pad = R.pad_states(recs)
lig = z1['lig']
split = np.array([r['split1'] for r in recs], dtype=np.int8)
train_js = [j for j in range(len(recs)) if split[j]==0]
rng = np.random.default_rng(20260821)
U = rng.normal(0,1,(O.D_RES,4)).astype(np.float32)
V = rng.normal(0,1,(O.D_LIG,4)).astype(np.float32)
zL = lig.astype(np.float32) @ V
hbar = {}
for j,r in enumerate(recs):
    hbar.setdefault(r['parent'], r['res_w'].mean(0).numpy())
def scale(x, var):
    v = x.var(); return x*np.sqrt(var/v)
g_par = {p: scale((h@U)@zL.T, 134.8) for p,h in hbar.items()}
s0 = scale((zL[:,0]+zL[:,1])*0.5, 50.0)
for j,r in enumerate(recs):
    u = rng.normal(0,1,4).astype(np.float32)
    prof = (s0 + g_par[r['parent']] + scale(u@zL.T,44.85) + rng.normal(0,np.sqrt(44.85),183).astype(np.float32))[r['lig_idx']]
    r['c'] = (prof - prof.mean()).astype(np.float32)

device = R.DEVICE
res_pad, mask_pad = res_pad.to(device), mask_pad.to(device)
ligT = torch.from_numpy(lig).float().to(device)
tgt = {j: torch.from_numpy(np.asarray(recs[j]['c'])).float().to(device) for j in range(len(recs))}

for name in ['A1-bilinear','A2-router']:
    torch.manual_seed(11)
    router = name=='A2-router'
    model = O.OLRPotential(router=router).to(device)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    for ep in range(301):
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
            total = total + ((tgt[j]-chat)**2).mean()
        (total/len(train_js)).backward()
        opt.step()
        if ep % 60 == 0:
            print(name, 'ep', ep, 'train MSE', round(float(total/len(train_js)),2), flush=True)
