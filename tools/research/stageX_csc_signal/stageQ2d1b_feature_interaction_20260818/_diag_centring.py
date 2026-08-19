"""Diagnostic (not frozen): attribute the oracle failure. Variants:
(a) current frozen truth (train-only ID centring), (b) global all-cell
centring, (c) no centring (grand mean only). ALS with 30 iterations.
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import truth
import q2

fz = np.load(HERE / 'q2d1b_features.npz', allow_pickle=False)
P_t = fz['P_t'].astype(np.float64)
L_t = fz['L_t'].astype(np.float64)
splits = __import__('json').loads((HERE / 'Q2D1B_SPLITS.json').read_text(encoding='utf-8'))
for k in ('train_cells', 'pc', 'lc', 'dc'):
    splits[k] = np.asarray(splits[k], dtype=np.int64)


def als(I_tr, r, l, iters=30, seed=0):
    A = np.random.default_rng(seed).normal(0, 0.1, size=(510, 4))
    B = np.random.default_rng(seed + 1).normal(0, 0.1, size=(64, 4))
    y = I_tr.astype(np.float64)
    Xr = P_t[r]
    Xl = L_t[l]
    for _ in range(iters):
        D = Xl @ B
        des = np.hstack([Xr * D[:, k:k + 1] for k in range(4)])
        a, *_ = np.linalg.lstsq(des, y, rcond=None)
        A = a.reshape(510, 4)
        C = Xr @ A
        des2 = np.hstack([Xl * C[:, k:k + 1] for k in range(4)])
        b, *_ = np.linalg.lstsq(des2, y, rcond=None)
        B = b.reshape(64, 4)
    return A, B


def eval_surf(A, B, cells, I_full):
    r, l = cells[:, 0], cells[:, 1]
    hat = ((P_t[r] @ A) * (L_t[l] @ B)).sum(-1)
    return q2.eval_metrics(hat, I_full[r, l])


def make_I(variant):
    t = truth.generate_truth('M1', 0, P_t.astype(np.float32), L_t.astype(np.float32), splits)
    I_raw = t['I_raw']
    tr = splits['train_cells']
    if variant == 'train_only_centred':
        return t['I']
    if variant == 'global_centred':
        rm = I_raw.mean(axis=1)
        cm = I_raw.mean(axis=0)
        Ic = I_raw - rm[:, None] - cm[None, :]
        sd = Ic[tr[:, 0], tr[:, 1]].std()
        return Ic / sd
    if variant == 'grand_mean_only':
        Ic = I_raw - I_raw.mean()
        sd = Ic[tr[:, 0], tr[:, 1]].std()
        return Ic / sd


for variant in ('train_only_centred', 'global_centred', 'grand_mean_only'):
    I_full = make_I(variant)
    tr = splits['train_cells']
    A, B = als(I_full[tr[:, 0], tr[:, 1]], tr[:, 0], tr[:, 1])
    row = {}
    for surf in ('pc', 'lc', 'dc'):
        m = eval_surf(A, B, splits[surf], I_full)
        row[surf] = (round(m['dead_zone_sign_accuracy'], 3), round(m['spearman'], 3))
    print(variant, row, flush=True)
