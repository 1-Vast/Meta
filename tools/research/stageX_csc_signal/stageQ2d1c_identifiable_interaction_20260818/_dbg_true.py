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
for seed in (0, 1, 2):
    t = truth_c.generate_truth('M1', seed, P_t, L_t, splits)
    for surf in ('pc', 'lc', 'dc'):
        cells = splits[surf]
        r, l = cells[:, 0], cells[:, 1]
        hat = ((P_t[r] @ t['A']) * (L_t[l] @ t['B'])).sum(-1) + P_t[r] @ t['w_r'] + L_t[l] @ t['w_c']
        m = q2.eval_metrics(hat, t['I'][r, l])
        print('TRUE-weights', seed, surf, round(m['dead_zone_sign_accuracy'], 3), round(m['spearman'], 3))
    tr = splits['train_cells']
    hat = ((P_t[tr[:, 0]] @ t['A']) * (L_t[tr[:, 1]] @ t['B'])).sum(-1) + P_t[tr[:, 0]] @ t['w_r'] + L_t[tr[:, 1]] @ t['w_c']
    m = q2.eval_metrics(hat, t['I'][tr[:, 0], tr[:, 1]])
    print('TRUE-weights', seed, 'train', round(m['dead_zone_sign_accuracy'], 3), round(m['spearman'], 3))
