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
t = truth_c.generate_truth('M1', 0, P_t, L_t, splits)
tr = splits['train_cells']
A, B = truth_c.svd_fit(t['I'][tr[:, 0], tr[:, 1]], tr[:, 0], tr[:, 1], P_t, L_t)
m = truth_c.eval_oracle(A, B, P_t, L_t, tr, t['I'], 'M1')
print('SVD in-fit dz', round(m['dead_zone_sign_accuracy'], 3), 'sp', round(m['spearman'], 3))
for surf in ('pc', 'lc', 'dc'):
    m = truth_c.eval_oracle(A, B, P_t, L_t, splits[surf], t['I'], 'M1')
    print(surf, round(m['dead_zone_sign_accuracy'], 3), round(m['spearman'], 3))
