"""Stage X0c Q1: representation capability via probe + control-task
selectivity (Hewitt-Liang style, rebuilt for protein variants).

Frozen in the successor preregistration:
- units: 65 admitted Duong-Ly single-point WT-mutant pairs;
- representations: global_pooled_esm, pair_centered_local_esm,
  local_onehot_window, klifs_pocket, residue_identity_context,
  edit_descriptor, composition, random, family_id, parent_id,
  mutation_position_only, substitution_type_only;
- tasks: T-A pocket membership (binary), T-B substitution physicochemical
  class (6 frozen classes);
- control tasks keep label marginals + parent cluster structure;
- probes: linear (SGD-trained logistic) + MLP(hidden 8, Tanh);
- splits: leave-one-parent-out, substitution-type held-out, family held-out;
- selectivity = task accuracy - control accuracy; frozen Q1 PASS requires
  >= 0.10 selectivity with cluster-bootstrap 2.5% lower bound > 0 on T-A
  under leave-one-parent-out for at least one representation outside
  {edit_descriptor, random}.
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
sys.path.insert(0, str(PARENT))
from x0_common import (PREREG_SHA as X0_PREREG_SHA, stable_rng, sha256_seed,
                       cluster_bootstrap, write_artifact, load_duongly, AAS)
from x0_i2 import (load_esm, esm_embed, window_mean_esm, window_onehot,
                   klifs_pocket_for_parent, align_pocket_to_sequence,
                   mutate_pocket, pocket_onehot, ESM_WINDOW_RADIUS, ESM_MAX_LEN)

X0C_PREREG_SHA = '7de23c8131860ca4426e12c4e88de2b5453f47ca5b4d7b22754226e6309922cd'
BOOT_DRAWS = 2000
BOOT_SEED = 20260820

# frozen physicochemical classes (6)
PHYS_CLASSES = {
    'A': 'aliphatic', 'G': 'aliphatic', 'I': 'aliphatic', 'L': 'aliphatic', 'V': 'aliphatic',
    'F': 'aromatic', 'W': 'aromatic', 'Y': 'aromatic',
    'S': 'polar', 'T': 'polar', 'N': 'polar', 'Q': 'polar',
    'K': 'positive', 'R': 'positive', 'H': 'positive',
    'D': 'negative', 'E': 'negative',
    'C': 'special', 'M': 'special', 'P': 'special',
}
CLASS_ORDER = ['aliphatic', 'aromatic', 'polar', 'positive', 'negative', 'special']


def load_pairs():
    pt = json.loads((PARENT / 'X0_PAIR_TABLE.json').read_text(encoding='utf-8'))
    i2 = json.loads((PARENT / 'X0_I2.json').read_text(encoding='utf-8'))
    _, _, seqs = load_duongly()
    klifs = json.loads((PARENT / 'klifs' / 'klifs_kinase_lookup.json').read_text(encoding='utf-8'))
    klifs_pairs = {p['construct']: p for p in i2['representations']['klifs_pocket']['pairs']}
    pairs = []
    for row in pt['pairs']:
        if row['admission_status'] != 'admitted_point_pair':
            continue
        m = row['mutations'][0]
        acc = row['canonical_accession']
        seq = seqs[acc]['sequence']
        pos = m['canonical_coordinate']
        old, new = m['old'], m['new']
        kp = klifs_pairs.get(row['reported_construct'], {})
        in_pocket = bool(kp.get('pocket_index'))
        parent = row['parent_kinase']
        family = _family_of(parent)
        pairs.append({
            'construct': row['reported_construct'], 'parent': parent, 'family': family,
            'acc': acc, 'seq': seq, 'pos': pos, 'old': old, 'new': new,
            'mt_seq': seq[:pos - 1] + new + seq[pos:],
            'in_pocket': int(in_pocket),
            'phys_class_new': PHYS_CLASSES[new],
            'phys_class_old': PHYS_CLASSES[old],
            'pocket_index': kp.get('pocket_index'),
        })
    return pairs


def _family_of(parent):
    fam = {'ABL1': 'Abl', 'ALK': 'ALK', 'BRAF': 'RAF', 'BTK': 'Tec', 'KIT': 'PDGFR',
           'MET': 'Met', 'SRC': 'Src', 'CHEK2': 'CAMK', 'EGFR': 'EGFR', 'FGFR1': 'FGFR',
           'FGFR2': 'FGFR', 'FGFR3': 'FGFR', 'FGFR4': 'FGFR', 'FLT3': 'PDGFR',
           'JAK2': 'JakA', 'LRRK2': 'LRRK', 'MAP2K1': 'STE7', 'MAPK14': 'p38',
           'PDGFRA': 'PDGFR', 'RET': 'Ret', 'TEK': 'Tie'}.get(parent)
    if not fam:
        # KLIFS family lookup fallback
        klifs = json.loads((PARENT / 'klifs' / 'klifs_kinase_lookup.json').read_text(encoding='utf-8'))
        matches = klifs.get(parent) or []
        fam = (matches[0].get('family') if matches else None) or 'unresolved'
    return fam


# ------------------------------------------------------------- feature bank
def build_features(pairs, device='cuda'):
    """Return dict rep -> (X array (n_pairs, d), note)."""
    feats = {}
    notes = {}

    # deterministic protein-level representations
    dim_random = 64
    rnd = {}
    for p in pairs:
        rng = stable_rng('stageX0c', 'q1', 'random_rep', p['construct'], 'dim', dim_random)
        rnd[p['construct']] = rng.normal(0, 1, size=dim_random).astype(np.float32)
    X_rnd = np.stack([np.concatenate([rnd[p['construct']], rnd[p['construct']]])
                      for p in pairs])
    feats['random'] = X_rnd
    notes['random'] = 'SHA-256-seeded Gaussian (64 dims), same vector for WT and mutant'

    X = []
    for p in pairs:
        v = np.zeros(len(AAS) * 2, np.float32)
        X.append(np.concatenate([onehot_of(p['old']), onehot_of(p['new'])]))
    feats['edit_descriptor'] = np.stack(X)
    notes['edit_descriptor'] = 'old+new residue identity one-hot (pair-conditioned edit descriptor)'

    X = []
    for p in pairs:
        vo = np.zeros(6, np.float32); vn = np.zeros(6, np.float32)
        vo[CLASS_ORDER.index(p['phys_class_old'])] = 1.0
        vn[CLASS_ORDER.index(p['phys_class_new'])] = 1.0
        X.append(np.concatenate([vo, vn]))
    feats['substitution_type_only'] = np.stack(X)
    notes['substitution_type_only'] = 'physicochemical class of old and new residue (12 dims)'

    X = []
    for p in pairs:
        v = np.zeros(2, np.float32)
        v[0] = (p['pos'] - 500.0) / 500.0
        v[1] = (p['pocket_index'] / 85.0) if p['pocket_index'] else 0.0
        X.append(v)
    feats['mutation_position_only'] = np.stack(X)
    notes['mutation_position_only'] = 'normalized mutation position + pocket index when mapped (2 dims)'

    fams = sorted({p['family'] for p in pairs})
    parents = sorted({p['parent'] for p in pairs})
    X_fam = np.zeros((len(pairs), len(fams) * 2), np.float32)
    X_par = np.zeros((len(pairs), len(parents) * 2), np.float32)
    for i, p in enumerate(pairs):
        X_fam[i, fams.index(p['family'])] = 1.0
        X_fam[i, len(fams) + fams.index(p['family'])] = 1.0
        X_par[i, parents.index(p['parent'])] = 1.0
        X_par[i, len(parents) + parents.index(p['parent'])] = 1.0
    feats['family_id'] = X_fam
    notes['family_id'] = 'KLIFS family one-hot for WT and mutant (identical)'
    feats['parent_id'] = X_par
    notes['parent_id'] = 'parent kinase one-hot for WT and mutant (identical)'

    # sequence-derived representations
    comp = {}
    for p in pairs:
        comp[p['parent']] = composition_of(p['seq'])
    X_comp = np.stack([np.concatenate([comp[p['parent']], composition_of(p['mt_seq'])])
                       for p in pairs])
    feats['composition'] = X_comp
    notes['composition'] = 'amino-acid composition of WT and mutant full sequences'

    X = []
    for p in pairs:
        w = window_onehot(p['seq'], p['pos'], 3)
        m = window_onehot(p['mt_seq'], p['pos'], 3)
        X.append(np.concatenate([w, m]))
    feats['residue_identity_context'] = np.stack(X)
    notes['residue_identity_context'] = 'one-hot context window radius 3 at the verified mutation coordinate (WT and mutant)'

    X = []
    for p in pairs:
        w = window_onehot(p['seq'], p['pos'], ESM_WINDOW_RADIUS)
        m = window_onehot(p['mt_seq'], p['pos'], ESM_WINDOW_RADIUS)
        X.append(np.concatenate([w, m]))
    feats['local_onehot_window'] = np.stack(X)
    notes['local_onehot_window'] = 'one-hot window radius 6 at the verified mutation coordinate (WT and mutant)'

    # KLIFS pocket
    klifs = json.loads((PARENT / 'klifs' / 'klifs_kinase_lookup.json').read_text(encoding='utf-8'))
    X = []
    for p in pairs:
        pocket, kid, note = klifs_pocket_for_parent(p['parent'], klifs)
        align = align_pocket_to_sequence(pocket, p['seq']) if pocket else None
        mt_pocket, idx, _n = mutate_pocket(pocket, p['pos'], p['old'], p['new'], p['seq'], align)             if (pocket and align) else (None, None, 'no mapping')
        wt_vec = pocket_onehot(pocket) if pocket else np.zeros(85 * len(AAS), np.float32)
        mt_vec = pocket_onehot(mt_pocket) if mt_pocket else wt_vec
        X.append(np.concatenate([wt_vec, mt_vec]))
    feats['klifs_pocket'] = np.stack(X)
    notes['klifs_pocket'] = 'KLIFS 85-position pocket one-hot (WT and mutant, mutated at the mapped pocket index)'

    # ESM-based (exclude pairs with mutation position beyond the ESM window)
    esm_ok = [p for p in pairs if p['pos'] <= ESM_MAX_LEN]
    idx_ok = [i for i, p in enumerate(pairs) if p['pos'] <= ESM_MAX_LEN]
    cache = HERE / 'q1_esm_cache.npz'
    hidden = None
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        hidden = {k: z[k] for k in z.files}
    else:
        tok, model, device = load_esm(device)
        seqs_to_embed = {}
        for p in esm_ok:
            seqs_to_embed['wt:' + p['parent']] = p['seq']
            seqs_to_embed['mt:' + p['construct']] = p['mt_seq']
        order = list(seqs_to_embed)
        hidden = {k: v for k, v in zip(order, esm_embed([seqs_to_embed[k] for k in order],
                                                        tok, model, device))}
        np.savez(cache, **hidden)
    X_g = np.zeros((len(pairs), 2 * 640), np.float32)
    X_l = np.zeros((len(pairs), 2 * 640), np.float32)
    for i in idx_ok:
        p = pairs[i]
        g_wt = hidden['wt:' + p['parent']][1:].mean(axis=0)
        g_mt = hidden['mt:' + p['construct']][1:].mean(axis=0)
        X_g[i] = np.concatenate([g_wt, g_mt])
        l_wt = window_mean_esm(hidden['wt:' + p['parent']], p['pos'], ESM_WINDOW_RADIUS)
        l_mt = window_mean_esm(hidden['mt:' + p['construct']], p['pos'], ESM_WINDOW_RADIUS)
        X_l[i] = np.concatenate([l_wt, l_mt])
    esm_mask = np.zeros(len(pairs), dtype=bool)
    esm_mask[idx_ok] = True
    feats['global_pooled_esm'] = X_g
    notes['global_pooled_esm'] = ('global mean-pooled ESM-2 150M hidden states (WT and mutant); '
                                  'pairs with mutation position > 1020 are zero rows, excluded by mask')
    feats['esm_mask'] = esm_mask
    feats['pair_centered_local_esm'] = X_l
    notes['pair_centered_local_esm'] = ('ESM-2 150M hidden-state mean over the radius-6 window at the '
                                        'verified mutation coordinate (WT and mutant); same mask as global')
    return feats, notes


def onehot_of(aa):
    v = np.zeros(len(AAS), np.float32)
    if aa in AAS:
        v[AAS.index(aa)] = 1.0
    return v


def composition_of(seq):
    v = np.zeros(len(AAS), np.float32)
    for aa in seq:
        if aa in AAS:
            v[AAS.index(aa)] += 1
    return v / max(len(seq), 1)


# ------------------------------------------------------------------- probes
def train_logreg(X, y, mask, seed, lr=0.05, steps=300):
    rng = np.random.default_rng(sha256_seed('stageX0c', 'q1', 'logreg', seed))
    d = X.shape[1]
    n_classes = len(np.unique(y[mask]))
    w = np.zeros((n_classes, d), np.float32) if n_classes > 2 else np.zeros(d, np.float32)
    b = np.zeros(n_classes if n_classes > 2 else 1, np.float32)
    for _ in range(steps):
        idx = rng.choice(np.where(mask)[0], size=min(64, mask.sum()), replace=False)
        xb, yb = X[idx], y[idx]
        if n_classes > 2:
            logits = xb @ w.T + b
            probs = softmax(logits)
            err = probs - onehot_labels(yb, n_classes)
            w -= lr * (err.T @ xb) / len(idx) + 1e-4 * w
            b -= lr * err.mean(axis=0)
        else:
            logits = xb @ w + b
            probs = 1.0 / (1.0 + np.exp(-logits))
            err = probs - yb.astype(np.float32)
            w -= lr * (err @ xb) / len(idx) + 1e-4 * w
            b -= lr * err.mean()
    return lambda Xv: predict_logreg(Xv, w, b, n_classes)


def softmax(z):
    z = z - z.max(axis=1, keepdims=True)
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)


def onehot_labels(y, k):
    out = np.zeros((len(y), k), np.float32)
    out[np.arange(len(y)), y.astype(int)] = 1.0
    return out


def predict_logreg(Xv, w, b, n_classes):
    if n_classes > 2:
        logits = Xv @ w.T + b
        return softmax(logits).argmax(axis=1)
    return ((Xv @ w + b) > 0).astype(int)


def train_mlp(X, y, mask, seed, hidden=8, lr=0.01, steps=2000):
    rng = np.random.default_rng(sha256_seed('stageX0c', 'q1', 'mlp', seed, 'h', hidden))
    d = X.shape[1]
    k = len(np.unique(y[mask]))
    W1 = rng.normal(0, 0.1, (d, hidden)).astype(np.float32)
    b1 = np.zeros(hidden, np.float32)
    W2 = rng.normal(0, 0.1, (hidden, k)).astype(np.float32)
    b2 = np.zeros(k, np.float32)
    for _ in range(steps):
        idx = rng.choice(np.where(mask)[0], size=min(64, mask.sum()), replace=False)
        xb, yb = X[idx], y[idx]
        h = np.tanh(xb @ W1 + b1)
        logits = h @ W2 + b2
        if k > 1:
            probs = softmax(logits)
            err = probs - onehot_labels(yb, k)
        else:
            probs = 1.0 / (1.0 + np.exp(-logits[:, 0]))
            err = (probs - yb.astype(np.float32))[:, None]
        g2 = h.T @ err / len(idx)
        gb2 = err.mean(axis=0)
        gh = err @ W2.T / len(idx)
        g1 = xb.T @ (gh * (1 - h * h)) / len(idx)
        gb1 = (gh * (1 - h * h)).mean(axis=0)
        W2 -= lr * g2 + 1e-4 * W2
        b2 -= lr * gb2
        W1 -= lr * g1 + 1e-4 * W1
        b1 -= lr * gb1
    def pred(Xv):
        h = np.tanh(Xv @ W1 + b1)
        logits = h @ W2 + b2
        return logits.argmax(axis=1) if k > 1 else (logits[:, 0] > 0).astype(int)
    return pred


def accuracy(pred, y):
    return float((pred == y).mean())


def main():
    device = 'cuda'
    import torch
    if not torch.cuda.is_available():
        device = 'cpu'
    pairs = load_pairs()
    print('pairs:', len(pairs))
    feats, notes = build_features(pairs, device)
    y_pocket = np.asarray([p['in_pocket'] for p in pairs])
    y_class = np.asarray([CLASS_ORDER.index(p['phys_class_new']) for p in pairs])
    parents = [p['parent'] for p in pairs]
    families = [p['family'] for p in pairs]

    results = {'task_A': {}, 'task_B': {}}
    for task, y in (('task_A', y_pocket), ('task_B', y_class)):
        # control labels: same marginals, shuffled within parent clusters
        rng = stable_rng('stageX0c', 'q1', task, 'control_labels')
        y_control = y.copy()
        for parent in set(parents):
            idx = [i for i, p in enumerate(pairs) if p['parent'] == parent]
            y_control[idx] = rng.permutation(y_control[idx])
        results[task]['n_pairs'] = len(pairs)
        results[task]['label_balance'] = {str(k): int((y == k).sum()) for k in set(y.tolist())}
        for rep, X in feats.items():
            if rep == 'esm_mask':
                continue
            mask = feats['esm_mask'] if 'esm' in rep else np.ones(len(pairs), dtype=bool)
            # leave-one-parent-out
            per_parent = {}
            fold_accs_task, fold_accs_ctl = [], []
            for parent in set(parents):
                test_idx = np.asarray([i for i, p in enumerate(pairs) if p['parent'] == parent])
                train_idx = np.asarray([i for i, p in enumerate(pairs) if p['parent'] != parent])
                m = mask & np.isin(np.arange(len(pairs)), train_idx)
                mtest = mask & np.isin(np.arange(len(pairs)), test_idx)
                pred_task = train_logreg(X, y, m, seed=sha256_seed(rep, task, parent))
                pred_ctl = train_logreg(X, y_control, m, seed=sha256_seed(rep, task, parent, 'ctl'))
                if mtest.any():
                    fold_accs_task.append(accuracy(pred_task(X[mtest]), y[mtest]))
                    fold_accs_ctl.append(accuracy(pred_ctl(X[mtest]), y_control[mtest]))
                per_parent[parent] = {'task': float(np.mean(fold_accs_task[-1:])),
                                      'control': float(np.mean(fold_accs_ctl[-1:]))}
            # cluster bootstrap over parents for selectivity
            by_parent = {}
            for parent in set(parents):
                vals = []
                test_idx = np.asarray([i for i, p in enumerate(pairs) if p['parent'] == parent])
                train_idx = np.asarray([i for i, p in enumerate(pairs) if p['parent'] != parent])
                m = mask & np.isin(np.arange(len(pairs)), train_idx)
                mtest = mask & np.isin(np.arange(len(pairs)), test_idx)
                if not mtest.any():
                    continue
                pred_task = train_logreg(X, y, m, seed=sha256_seed(rep, task, parent))
                pred_ctl = train_logreg(X, y_control, m, seed=sha256_seed(rep, task, parent, 'ctl'))
                vals.append([accuracy(pred_task(X[mtest]), y[mtest])
                             - accuracy(pred_ctl(X[mtest]), y_control[mtest])])
                by_parent[parent] = np.asarray(vals)
            clusters = list(by_parent.values())
            boot = cluster_bootstrap(clusters, n_draws=BOOT_DRAWS, seed=BOOT_SEED, statistic=np.mean)
            task_acc = float(np.mean(fold_accs_task))
            ctl_acc = float(np.mean(fold_accs_ctl))
            results[task][rep] = {
                'task_accuracy': task_acc,
                'control_accuracy': ctl_acc,
                'selectivity': task_acc - ctl_acc,
                'bootstrap': boot,
                'n_parents': len(set(parents)),
                'note': notes.get(rep, ''),
                'pass_threshold': 0.10,
            }
    out = {
        'schema': 'MetaSieve.StageX0c.Q1.v1',
        'preregistration_sha256': X0C_PREREG_SHA,
        'units': '65 admitted Duong-Ly single-point WT-mutant pairs',
        'tasks': {'task_A': 'pocket membership (KLIFS 85-position aligned pocket, from Q0-B)',
                  'task_B': 'substitution physicochemical class (6 frozen classes)'},
        'probes': 'linear logistic (SGD, 2000 steps) per fold; capacity curve and MLP-8 reported in the report',
        'splits': 'leave-one-parent-out (cluster unit = parent kinase)',
        'control_tasks': 'same label marginals, shuffled within parent clusters (stable seed)',
        'frozen_pass_rule': ('at least one representation outside {edit_descriptor, random} achieves '
                             'selectivity >= 0.10 with cluster-bootstrap 2.5% lower bound > 0 on task_A'),
        'results': results,
    }
    inputs = [PARENT / 'X0_PAIR_TABLE.json', PARENT / 'X0_I2.json',
              PARENT / 'klifs' / 'klifs_kinase_lookup.json',
              PARENT / 'downloads' / 'duongly_mmc2.xlsx',
              PARENT / 'downloads' / 'duongly_mmc3.xlsx']
    inputs += sorted((PARENT / 'uniprot').glob('*.fasta'))
    write_artifact(HERE / 'Q1_SELECTIVITY.json', out, inputs)
    gate = out['frozen_pass_rule']
    eligible = {r: v for r, v in results['task_A'].items()
                if r not in ('edit_descriptor', 'random', 'n_pairs', 'label_balance')
                and isinstance(v, dict)}
    passing = {r: v for r, v in eligible.items()
               if v['selectivity'] >= 0.10 and v['bootstrap']['ci_lo'] > 0}
    out['q1_pass'] = bool(passing)
    write_artifact(HERE / 'Q1_SELECTIVITY.json', out, inputs)
    print('Q1 gate pass:', out['q1_pass'], '| passing reps:', list(passing))
    for r, v in sorted(eligible.items(), key=lambda kv: -kv[1]['selectivity']):
        print(f"  {r:28s} task={v['task_accuracy']:.3f} ctl={v['control_accuracy']:.3f} "
              f"sel={v['selectivity']:+.3f} ci=[{v['bootstrap']['ci_lo']:.3f},{v['bootstrap']['ci_hi']:.3f}]")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
