"""Stage Q2c-2: frozen gate rerun with pair_centered_local_esm (per-row
[wt;mt] window, 1280-dim) as the correct-arm protein input. Protocol,
thresholds and negative controls unchanged from the X0c gate (0.30 / 0.70 /
0.05; median over 3 seeds; label permutation, tau*=0, floor-clamp,
no-interaction head, shuffled/family-shuffled/random protein, free-target-id;
cluster bootstrap over eval parents). RUN ONLY AFTER Q2c-1 permits.
"""
import json, sys, time
from pathlib import Path
import numpy as np
import torch

HERE = Path(__file__).resolve().parent
X0C = HERE.parent / 'stageX0c_measurement_qualification_20260818'
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(X0C.parent))
import q2
from q2 import (load_duongly_graph, generate, build_arm_feats, train_arm_with_restarts,
                FreeIDModel, train, predict, eval_metrics, cluster_boot_pairs)
from x0_common import sha256_seed, stable_rng

PREREG_SHA = '1027ccde8c8946aa8314ebd7642af89a6abbc3366afd965e8ab43f0da5a26a5c'
SEEDS = [0, 1, 2]


def main():
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    (rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells,
     row_of_cell, lig_of_cell) = load_duongly_graph()
    n_prot = len(rows)
    fz = np.load(HERE / 'q2c_row_esm.npz', allow_pickle=False)
    esm_prot = fz['X'].astype(np.float32)
    eval_c = splits['eval_cells']

    def train_and_eval(P, lat, seed, arm='correct_protein'):
        model, val = train_arm_with_restarts(
            P, arm, lat, seed, device, n_prot, P, lig_feats,
            row_of_cell, lig_of_cell, splits)
        out = predict(model, P, lig_feats, row_of_cell[eval_c],
                      lig_of_cell[eval_c], device, arm)
        return eval_metrics(out['inter'], lat['I_cells'][eval_c]), out['inter'], model

    # main gate: correct arm = esm_prot; ligand_only baseline
    saved = {}
    results = {}
    for s in SEEDS:
        lat = generate(rows, compounds, prot_feats, lig_feats, cells, 1.0, 4, 'dense', False, s)
        m_c, inter_c, _ = train_and_eval(esm_prot, lat, s, 'correct_protein')
        m_l, inter_l, _ = train_and_eval(np.zeros_like(esm_prot), lat, s, 'ligand_only')
        saved[('correct', s)] = (inter_c, lat['I_cells'][eval_c])
        saved[('ligand_only', s)] = (inter_l, lat['I_cells'][eval_c])
        results[f'seed{s}'] = {'correct': m_c, 'ligand_only': m_l}
        print(f'seed{s}: correct dz={m_c["dead_zone_sign_accuracy"]:.3f} sp={m_c["spearman"]:.3f} '
              f'| ligand_only dz={m_l["dead_zone_sign_accuracy"]:.3f} '
              f'gap={m_c["sign_accuracy"]-m_l["sign_accuracy"]:+.3f}', flush=True)

    dz = float(np.median([results[f'seed{s}']['correct']['dead_zone_sign_accuracy'] for s in SEEDS]))
    sp = float(np.median([results[f'seed{s}']['correct']['spearman'] for s in SEEDS]))
    gap = float(np.median([results[f'seed{s}']['correct']['sign_accuracy'] -
                           results[f'seed{s}']['ligand_only']['sign_accuracy'] for s in SEEDS]))
    gate = {'required': {'spearman': 0.30, 'dead_zone': 0.70, 'gap': 0.05},
            'observed': {'spearman_median': sp, 'dead_zone_median': dz, 'gap_median': gap},
            'pass': bool(sp >= 0.30 and dz >= 0.70 and gap >= 0.05)}

    # negative controls (seed 0, same arm protocols)
    nc = {}
    lat0 = generate(rows, compounds, prot_feats, lig_feats, cells, 1.0, 4, 'dense', False, 0)
    # label permutation
    rng = stable_rng('stageQ2c', 'q2c2', 'label_perm')
    perm = rng.permutation(len(cells))
    lat_perm = dict(lat0)
    lat_perm['I_cells'] = lat0['I_cells'][perm]
    m_p, _, _ = train_and_eval(esm_prot, lat_perm, 0)
    nc['label_permutation'] = {'dz': m_p['dead_zone_sign_accuracy'], 'sp': m_p['spearman']}
    print('label_perm dz=', round(m_p['dead_zone_sign_accuracy'], 3), flush=True)
    # tau*=0
    lat_t0 = generate(rows, compounds, prot_feats, lig_feats, cells, 0.0, 4, 'dense', False, 0)
    m_t0, _, _ = train_and_eval(esm_prot, lat_t0, 0)
    nc['tau0'] = {'dz': m_t0['dead_zone_sign_accuracy'], 'sp': m_t0['spearman']}
    print('tau0 dz=', round(m_t0['dead_zone_sign_accuracy'], 3), flush=True)
    # floor-clamp imputation
    lat_floor = dict(lat0)
    zf = lat0['z_obs'].copy()
    zf[np.isnan(zf)] = np.log(0.5 / 99.5)
    lat_floor['z_obs'] = zf
    lat_floor['determinate'] = np.ones(len(cells), dtype=bool)
    m_f, _, _ = train_and_eval(esm_prot, lat_floor, 0)
    nc['floor_imputation'] = {'dz': m_f['dead_zone_sign_accuracy'], 'sp': m_f['spearman']}
    print('floor dz=', round(m_f['dead_zone_sign_accuracy'], 3), flush=True)
    # shuffled / random protein
    arm_feats = build_arm_feats(rows, prot_feats, lat0)
    for arm in ('shuffled_protein', 'family_preserving_shuffle', 'random_protein', 'no_interaction_head'):
        model, val = train_arm_with_restarts(
            arm_feats[arm], arm, lat0, 0, device, n_prot, arm_feats[arm], lig_feats,
            row_of_cell, lig_of_cell, splits)
        if arm == 'no_interaction_head':
            with torch.no_grad():
                model.inter_scale.fill_(0.0)
        out = predict(model, arm_feats[arm], lig_feats, row_of_cell[eval_c],
                      lig_of_cell[eval_c], device, arm)
        m_nc = eval_metrics(out['inter'], lat0['I_cells'][eval_c])
        nc[arm] = {'dz': m_nc['dead_zone_sign_accuracy'], 'sp': m_nc['spearman']}
        print(arm, 'dz=', round(m_nc['dead_zone_sign_accuracy'], 3), flush=True)

    # cluster bootstrap over eval parents (correct arm, seed 0): dead-zone
    # sign agreement per eval cell, resampled at parent level
    from x0_common import normalize_parent_name
    dz_mask = np.abs(saved[('correct', 0)][1]) > q2.DEAD_ZONE
    sign_ok = (np.sign(saved[('correct', 0)][0]) == np.sign(saved[('correct', 0)][1]))[dz_mask].astype(float)
    par_ok = np.asarray([normalize_parent_name(rows[r]) for r in row_of_cell[eval_c]])[dz_mask]
    by_parent = {}
    for v, p in zip(sign_ok, par_ok):
        by_parent.setdefault(p, []).append(v)
    clusters = list(by_parent.values())
    bs_dz = cluster_boot_pairs(clusters, stat=lambda x: float(np.mean(x)))
    bs = {'dead_zone_sign_agreement': bs_dz, 'n_eval_parents': len(clusters)}

    out = {
        'schema': 'MetaSieve.StageQ2c.Q2C2_GATE.v1',
        'preregistration_sha256': PREREG_SHA,
        'gate': gate,
        'q2c2_pass': gate['pass'],
        'per_seed': results,
        'negative_controls': nc,
        'bootstrap': bs,
        'note': 'run only after Q2c-1 permits; thresholds and protocol unchanged from X0c gate',
    }
    (HERE / 'Q2C2_GATE.json').write_text(json.dumps(out, indent=1))
    print(json.dumps({'gate': gate, 'negative_controls': nc}, indent=1))
    return 0


def _parent_of(row):
    from x0_common import normalize_parent_name
    return normalize_parent_name(row)


if __name__ == '__main__':
    raise SystemExit(main())
