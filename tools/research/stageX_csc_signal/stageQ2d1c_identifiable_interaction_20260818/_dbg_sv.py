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
r, l = tr[:, 0], tr[:, 1]
I_tr = t['I'][r, l]
rows_u = sorted(set(r.tolist()))
ligs_u = sorted(set(l.tolist()))
M = np.zeros((len(rows_u), len(ligs_u)))
ri = {v: k for k, v in enumerate(rows_u)}
li = {v: k for k, v in enumerate(ligs_u)}
M[[ri[v] for v in r], [li[v] for v in l]] = I_tr
U, S, Vt = np.linalg.svd(M, full_matrices=False)
print('singular values of raw I_train:', np.round(S[:8], 3))
# residual after offset removal
trm = np.zeros(len(P_t)); tcm = np.zeros(len(L_t))
for i in rows_u:
    trm[i] = I_tr[np.asarray(r) == i].mean()
for j in ligs_u:
    tcm[j] = I_tr[np.asarray(l) == j].mean()
w_r, *_ = np.linalg.lstsq(P_t[rows_u].astype(np.float64), trm[rows_u], rcond=None)
w_c, *_ = np.linalg.lstsq(L_t[ligs_u].astype(np.float64), tcm[ligs_u], rcond=None)
resid = I_tr - P_t[r] @ w_r - L_t[l] @ w_c
M2 = np.zeros_like(M)
M2[[ri[v] for v in r], [li[v] for v in l]] = resid
U2, S2, Vt2 = np.linalg.svd(M2, full_matrices=False)
print('singular values after offset removal:', np.round(S2[:8], 3))
print('energy top-4 / total:', round(float((S2[:4]**2).sum() / (S2**2).sum()), 3))
