"""Debug Q2d-1c oracle: in-fit ALS quality, sd_tr, A-subspace recovery."""
import sys, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import truth_c
import q2

fz = np.load(HERE / 'q2d1c_features.npz', allow_pickle=False)
P_t = fz['P_t'].astype(np.float32)
L_t = fz['L_t'].astype(np.float32)
splits = json.loads((HERE / 'Q2D1C_SPLITS.json').read_text(encoding='utf-8'))
for k in ('train_cells', 'pc', 'lc', 'dc'):
    splits[k] = np.asarray(splits[k], dtype=np.int64)
for k in ('cold_row', 'cold_lig', 'train_row', 'train_lig'):
    splits[k] = np.asarray(splits[k], dtype=bool)

t = truth_c.generate_truth('M1', 0, P_t, L_t, splits)
print('sd_train:', round(t['sd_train'], 3))
tr = splits['train_cells']
r, l = tr[:, 0], tr[:, 1]
A, B = truth_c.als_fit(t['I'][r, l], r, l, P_t, L_t)
# in-fit metrics on train cells
hat = ((P_t[r] @ A) * (L_t[l] @ B)).sum(-1)
m = q2.eval_metrics(hat, t['I'][r, l])
print('ALS in-fit train: dz', round(m['dead_zone_sign_accuracy'], 3), 'sp', round(m['spearman'], 3))
# subspace recovery vs truth A
At = t['A']
from numpy.linalg import svd
_, s1, _ = svd(At, full_matrices=False)
# project: canonical correlation-ish: angle between column spaces
Qa, _ = np.linalg.qr(A)
Qt, _ = np.linalg.qr(At)
_, s_ang, _ = svd(Qa.T @ Qt)
print('subspace angles cos:', np.round(s_ang, 3))
# check raw bilinear truth vs raw bilinear of recovered weights WITHOUT centring
Iraw_hat = (P_t @ A) @ (L_t @ B).T
Iraw = t['I_raw']
m2 = q2.eval_metrics(Iraw_hat[r, l], Iraw[r, l])
print('raw in-fit: dz', round(m2['dead_zone_sign_accuracy'], 3), 'sp', round(m2['spearman'], 3))
# check: does the truth I_train even have variance?
print('I_train sd:', round(t['I'][r, l].std(), 3), 'mean', round(t['I'][r, l].mean(), 3))
