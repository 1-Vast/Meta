import sys
sys.path.insert(0, '.')
import numpy as np, json
fz = np.load('q2d1c_features.npz', allow_pickle=False)
P_t = fz['P_t'].astype(np.float32)
splits = json.loads(open('Q2D1C_SPLITS.json', encoding='utf-8').read())
trm = np.asarray(splits['train_row'], dtype=bool)
U, S, Vt = np.linalg.svd(P_t[trm].astype(np.float64), full_matrices=False)
print('train-row submatrix singular values:', np.round(S, 3))
print('rank > 1e-6:', int((S > 1e-6).sum()))
# what fraction of TRUE A energy lies in the null space
import truth_c
L_t = fz['L_t'].astype(np.float32)
for k in ('train_cells','pc','lc','dc'):
    splits[k] = np.asarray(splits[k], dtype=np.int64)
t = truth_c.generate_truth('M1', 1, P_t, L_t, splits)
A = t['A'].astype(np.float64)
# null-space projector
Vsp = Vt[:int((S > 1e-6).sum())].T  # (32, r) span basis
Psp = Vsp @ Vsp.T
null_frac = 1 - float(np.linalg.norm(Psp @ A)**2 / np.linalg.norm(A)**2)
print('fraction of true A energy in train-row null space:', round(null_frac, 3))
