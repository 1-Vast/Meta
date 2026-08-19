"""Attribution diagnostic for the correct arm (M1, level A, seed 0):
trained model's in-fit quality vs cold surfaces, and the null-space energy
of the learned protein map. CPU-only, runs beside the ladder.
"""
import sys, json
from pathlib import Path
import numpy as np
import torch

HERE = Path(__file__).resolve().parent
X0C = HERE.parent / "stageX0c_measurement_qualification_20260818"
sys.path.insert(0, str(X0C))
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
device = 'cpu'
t = truth.generate_truth('M1', 0, P_t, L_t, splits)
Vsp, r = truth._span_projection(P_t, splits)
Lt_dev = torch.from_numpy(L_t).float()
tr = splits['train_cells']
for restart in (0, 1, 2, 3):
    model, val, _ = R.train_level(P_t, 'correct', t, 'A', 0, splits, device, restart, Lt_dev)
    A = model.A.weight.detach().numpy().T  # (32, 4)
    null_frac = 1 - float(np.linalg.norm(Vsp @ (Vsp.T @ A)) ** 2 / np.linalg.norm(A) ** 2)
    with torch.no_grad():
        Pt = torch.from_numpy(P_t).float()
        o_tr = model(Pt[tr[:, 0]], Lt_dev[tr[:, 1]])
    m_tr = q2.eval_metrics(o_tr['inter'].numpy(), t['I'][tr[:, 0], tr[:, 1]])
    res = R.eval_arm(model, P_t, 'correct', t, 'A', splits, device, Lt_dev)
    print('restart', restart, 'train dz/sp', round(m_tr['dead_zone_sign_accuracy'], 3),
          round(m_tr['spearman'], 3), 'null_frac', round(null_frac, 3),
          'pc', round(res['pc']['dz'], 3), 'lc', round(res['lc']['dz'], 3),
          'dc', round(res['dc']['dz'], 3), 'mon', round(val, 3), flush=True)
