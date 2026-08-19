import sys
sys.path.insert(0, '.')
import numpy as np, json
import truth_c, q2
fz = np.load('q2d1c_features.npz', allow_pickle=False)
P_t = fz['P_t'].astype(np.float32)
L_t = fz['L_t'].astype(np.float32)
splits = json.loads(open('Q2D1C_SPLITS.json', encoding='utf-8').read())
for k in ('train_cells', 'pc', 'lc', 'dc'):
    splits[k] = np.asarray(splits[k], dtype=np.int64)
for k in ('cold_row', 'cold_lig', 'train_row', 'train_lig'):
    splits[k] = np.asarray(splits[k], dtype=bool)
t = truth_c.generate_truth('M1', 1, P_t, L_t, splits)
tr = splits['train_cells']
A, B, w_r, w_c = truth_c.svd_fit(t['I'][tr[:, 0], tr[:, 1]], tr[:, 0], tr[:, 1], P_t, L_t, splits)
cells = splits['pc']
r, l = cells[:, 0], cells[:, 1]
hat = ((P_t[r] @ A) * (L_t[l] @ B)).sum(-1) + P_t[r] @ w_r + L_t[l] @ w_c
m = q2.eval_metrics(hat, t['I'][r, l])
print('oracle pc dz', round(m['dead_zone_sign_accuracy'], 3), 'sp', round(m['spearman'], 3))
hat_t = ((P_t[r] @ t['A']) * (L_t[l] @ t['B'])).sum(-1) + P_t[r] @ t['w_r'] + L_t[l] @ t['w_c']
print('corr(oracle, true) on pc:', round(float(np.corrcoef(hat, hat_t)[0, 1]), 4))
print('scale oracle vs true:', round(float(hat.std() / hat_t.std()), 3))
# check A/B subspace recovery on train
Fr = None
