"""Stage Q2c shared feature builder v2: per-ROW mutation-aware pocket/ESM
representation. Each Duong-Ly row (WT parent or mutant construct) maps to
[wt_window_mean; mt_window_mean] (2x640, radius-6 window at the verified
mutation coordinate), i.e. the Q1-passing pair_centered_local_esm
representation generalized to the Q2 row graph. WT rows get the parent
window on both sides; rows whose construct lacks an ESM entry get the
parent window on both sides and are flagged (not masked) in the manifest.
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
X0C = HERE.parent / 'stageX0c_measurement_qualification_20260818'
PARENT = X0C.parent
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(PARENT))
import q1
from x0_i2 import window_mean_esm, ESM_WINDOW_RADIUS
from x0_common import normalize_parent_name

Q2C_PREREG_SHA = '1027ccde8c8946aa8314ebd7642af89a6abbc3366afd965e8ab43f0da5a26a5c'


def build():
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells,         row_of_cell, lig_of_cell = q2.load_duongly_graph()
    pairs = q1.load_pairs()
    cache = np.load(X0C / 'q1_esm_cache.npz', allow_pickle=False)

    wt_of = {}
    mt_of = {}
    pos_of = {}
    parent_seq_of = {}
    for p in pairs:
        wt_of[p['parent']] = p['seq']
        mt_of[p['construct']] = p['mt_seq']
        pos_of[p['construct']] = p['pos']
        parent_seq_of[p['parent']] = p['seq']

    manifest = []
    X = np.zeros((len(rows), 2 * 640), dtype=np.float32)
    for i, row in enumerate(rows):
        parent = normalize_parent_name(row)
        is_wt = (row == parent)
        wt_key = 'wt:' + parent
        mt_key = 'mt:' + row.replace('(', ' (')
        if wt_key not in cache:
            manifest.append({'row': row, 'status': 'missing_wt_cache', 'used': 'zero'})
            continue
        wt_hid = cache[wt_key]
        wt_win = window_mean_esm(wt_hid, pos_of.get(row, None) if is_wt else pos_of.get(row, None), ESM_WINDOW_RADIUS) if False else None
        # WT rows: parent window both sides
        if is_wt:
            seq = parent_seq_of.get(parent)
            pos = None
            # window at mid-pocket is not defined for WT rows; use full-pocket mean
            wt_vec = wt_hid[1:].mean(axis=0)
            X[i] = np.concatenate([wt_vec, wt_vec])
            manifest.append({'row': row, 'status': 'wt_parent_global_mean', 'parent': parent})
            continue
        pos = pos_of.get(row.replace('(', ' ('))
        if mt_key in cache and pos is not None:
            mt_hid = cache[mt_key]
            wt_vec = window_mean_esm(wt_hid, pos, ESM_WINDOW_RADIUS)
            mt_vec = window_mean_esm(mt_hid, pos, ESM_WINDOW_RADIUS)
            X[i] = np.concatenate([wt_vec, mt_vec])
            manifest.append({'row': row, 'status': 'ok_pair_window', 'parent': parent, 'pos': pos})
        else:
            wt_vec = (window_mean_esm(wt_hid, pos, ESM_WINDOW_RADIUS)
                      if pos is not None and pos <= 1020 else None)
            if wt_vec is None:
                wt_vec = wt_hid[1:].mean(axis=0)
            X[i] = np.concatenate([wt_vec, wt_vec])
            manifest.append({'row': row, 'status': 'fallback_parent_window_or_mean',
                             'parent': parent, 'reason': 'no mt cache entry or no verified pos'})
    np.savez(HERE / 'q2c_row_esm.npz', X=X, rows=np.asarray(rows))
    import hashlib
    h = hashlib.sha256((HERE / 'q2c_row_esm.npz').read_bytes()).hexdigest()
    n_ok = sum(1 for m in manifest if m['status'] == 'ok_pair_window')
    out = {'schema': 'MetaSieve.StageQ2c.FEATS.v2',
           'preregistration_sha256': Q2C_PREREG_SHA,
           'representation': ('per-row [wt window mean; mt window mean], radius-6 ESM-2 150M '
                              '(2x640=1280 dim); WT rows and fallback rows use parent mean twice'),
           'n_rows': len(rows), 'dim': int(X.shape[1]),
           'n_ok_pair_window': n_ok, 'manifest': manifest, 'npz_sha256': h}
    (HERE / 'Q2C_FEATS.json').write_text(json.dumps(out, indent=1))
    print(json.dumps({'n_rows': len(rows), 'dim': int(X.shape[1]), 'n_ok_pair_window': n_ok}))
    return 0


if __name__ == '__main__':
    raise SystemExit(build())
