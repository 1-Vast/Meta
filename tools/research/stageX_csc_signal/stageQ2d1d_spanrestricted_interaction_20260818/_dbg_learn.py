import sys, json
from pathlib import Path
import numpy as np
import torch
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "stageX0c_measurement_qualification_20260818"))
sys.path.insert(0, str(HERE))
import truth_d as truth
import runner_d as R
import q2

fz = np.load(HERE / 'q2d1d_features.npz', allow_pickle=False)
P_t = fz['P_t'].astype(np.float32)
L_t = fz['L_t'].astype(np.float32)
splits = json.loads(open(HERE / 'Q2D1D_SPLITS.json', encoding='utf-8').read())
for k in ('train_cells', 'pc', 'lc', 'dc'):
    splits[k] = np.asarray(splits[k], dtype=np.int64)
for k in ('cold_row', 'cold_lig', 'train_row', 'train_lig'):
    splits[k] = np.asarray(splits[k], dtype=bool)
device = "cuda" if torch.cuda.is_available() else "cpu"
print("device:", device, flush=True)
t = truth.generate_truth('M1', 0, P_t, L_t, splits)
P_or = (P_t.astype(np.float64) @ t['A']).astype(np.float32)
Lt_dev = torch.from_numpy(L_t).float().to(device)

model, val, n_cens = R.train_level(P_or, 'oracle_diagnostic', t, 'A', 0, splits, device, 0, Lt_dev)
print('final train loss:', val, flush=True)
res = R.eval_arm(model, P_or, 'oracle_diagnostic', t, 'A', splits, device, Lt_dev)
print('eval:', {s_: {k: round(v, 3) for k, v in res[s_].items()} for s_ in res}, flush=True)
# inspect scale
with torch.no_grad():
    r_ = torch.from_numpy(splits['dc'][:, 0]).to(device)
    l_ = torch.from_numpy(splits['dc'][:, 1]).to(device)
    Pt = torch.from_numpy(P_or).float().to(device)
    o = model(Pt[r_], Lt_dev[l_])
print('dc inter output: mean', round(float(o['inter'].mean()), 3), 'sd', round(float(o['inter'].std()), 3), flush=True)
print('dc truth I: mean', round(float(t['I'][splits['dc'][:, 0], splits['dc'][:, 1]].mean()), 3), 'sd', round(float(t['I'][splits['dc'][:, 0], splits['dc'][:, 1]].std()), 3), flush=True)
print('corr(out, truth):', round(float(np.corrcoef(o['inter'].cpu().numpy(), t['I'][splits['dc'][:, 0], splits['dc'][:, 1]])[0, 1]), 3), flush=True)
