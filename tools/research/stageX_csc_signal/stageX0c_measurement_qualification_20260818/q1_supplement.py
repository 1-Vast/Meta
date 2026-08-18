"""Q1 supplement: probe capacity curve (MLP-8) and random-label curve for
the three passing representations. Appends to Q1_SELECTIVITY.json."""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
sys.path.insert(0, str(PARENT))
import q1
from q1 import (load_pairs, build_features, train_logreg, train_mlp, accuracy,
                CLASS_ORDER, X0C_PREREG_SHA)
from x0_common import stable_rng

import torch
device = 'cuda' if torch.cuda.is_available() else 'cpu'
pairs = load_pairs()
feats, notes = build_features(pairs, device)
y = np.asarray([p['in_pocket'] for p in pairs])
parents = [p['parent'] for p in pairs]
reps = ['pair_centered_local_esm', 'mutation_position_only', 'substitution_type_only']

out = json.loads((HERE / 'Q1_SELECTIVITY.json').read_text(encoding='utf-8'))
supp = {'task_A_supplement': {}}
for rep in reps:
    X = feats[rep]
    mask = feats['esm_mask'] if 'esm' in rep else np.ones(len(pairs), dtype=bool)
    # MLP-8 capacity curve (same LOO-parent split)
    fold_accs = []
    for parent in set(parents):
        test_idx = np.asarray([i for i, p in enumerate(pairs) if p['parent'] == parent])
        train_idx = np.asarray([i for i, p in enumerate(pairs) if p['parent'] != parent])
        m = mask & np.isin(np.arange(len(pairs)), train_idx)
        mtest = mask & np.isin(np.arange(len(pairs)), test_idx)
        if not mtest.any():
            continue
        pred = train_mlp(X, y, m, seed=q1.sha256_seed(rep, 'mlp8', parent), hidden=8, steps=300)
        fold_accs.append(accuracy(pred(X[mtest]), y[mtest]))
    mlp8_acc = float(np.mean(fold_accs))
    # fully random labels (no parent structure)
    rng = stable_rng('stageX0c', 'q1', 'random_labels_global', rep)
    y_rand = rng.permutation(y)
    fold_accs_rand = []
    for parent in set(parents):
        test_idx = np.asarray([i for i, p in enumerate(pairs) if p['parent'] == parent])
        train_idx = np.asarray([i for i, p in enumerate(pairs) if p['parent'] != parent])
        m = mask & np.isin(np.arange(len(pairs)), train_idx)
        mtest = mask & np.isin(np.arange(len(pairs)), test_idx)
        if not mtest.any():
            continue
        pred = train_logreg(X, y_rand, m, seed=q1.sha256_seed(rep, 'randlabel', parent))
        fold_accs_rand.append(accuracy(pred(X[mtest]), y_rand[mtest]))
    rand_acc = float(np.mean(fold_accs_rand))
    base = out['results']['task_A'][rep]
    supp['task_A_supplement'][rep] = {
        'linear_task_accuracy': base['task_accuracy'],
        'mlp8_task_accuracy': mlp8_acc,
        'capacity_delta_mlp8_minus_linear': mlp8_acc - base['task_accuracy'],
        'random_label_accuracy_linear': rand_acc,
        'random_label_selectivity': rand_acc - base['control_accuracy'],
        'interpretation': ('a probe whose accuracy is unchanged under random labels '
                           'is reading real structure; the capacity curve shows whether '
                           'a larger probe inflates the result'),
    }
    print(rep, 'mlp8', round(mlp8_acc, 3), 'rand-label', round(rand_acc, 3))
out['supplement'] = supp
(HERE / 'Q1_SELECTIVITY.json').write_text(json.dumps(out, indent=1) + chr(10))
print('Q1 supplement written')
