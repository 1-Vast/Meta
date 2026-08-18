"""Q2d-1 diagnostic ONLY (not a deployment): closed-form rank-4
reconstruction upper bound. Fit mu+p_b+l_b by per-parent/per-ligand means on
TRAIN cells (z-scale), take residual, truncate to rank-4 SVD on the TRAIN
observed submatrix, evaluate dead-zone sign accuracy on eval cells.
"""
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
X0C = HERE.parent / 'stageX0c_measurement_qualification_20260818'
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(X0C.parent))
import q2
from q2 import load_duongly_graph, generate, eval_metrics

PREREG_SHA = '4f7e80027b9b82564bf1ea262d360813f23ffe77f8378ca54b490cc6752fa3d1'


def main():
    (rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells,
     row_of_cell, lig_of_cell) = load_duongly_graph()
    n_prot, n_lig = len(rows), len(compounds)
    out = {'per_seed': {}}
    for s in [0, 1, 2]:
        lat = generate(rows, compounds, prot_feats, lig_feats, cells, 1.0, 4, 'dense', False, s)
        z = lat['z']
        tr, ev = splits['train_cells'], splits['eval_cells']
        # main-effect fit on train (z-scale, no censoring)
        mu = float(z[tr].mean())
        pm = np.zeros(n_prot)
        lm = np.zeros(n_lig)
        for i in range(n_prot):
            m = tr[row_of_cell[tr] == i]
            pm[i] = float(z[m].mean() - mu) if len(m) else 0.0
        for j in range(n_lig):
            m = tr[lig_of_cell[tr] == j]
            lm[j] = float(z[m].mean() - mu) if len(m) else 0.0
        resid = z - mu - pm[row_of_cell] - lm[lig_of_cell]
        # rank-4 SVD of the train residual submatrix (diagnostic only)
        R = np.full((n_prot, n_lig), np.nan)
        R[row_of_cell[tr], lig_of_cell[tr]] = resid[tr]
        col_mean = np.nanmean(R, axis=0)
        col_mean = np.where(np.isnan(col_mean), 0.0, col_mean)
        row_mean = np.nanmean(R, axis=1)
        row_mean = np.where(np.isnan(row_mean), 0.0, row_mean)
        Rf = np.where(np.isnan(R), col_mean[None, :], R)
        Rf = Rf - np.nanmean(Rf)
        U, S, Vt = np.linalg.svd(Rf, full_matrices=False)
        R4 = (U[:, :4] * S[:4]) @ Vt[:4, :]
        # evaluate on val cells (train parents x val ligands: row factor known,
        # column factor unseen - a half-cold transfer bound)
        va = splits['val_cells']
        Ihat_val = R4[row_of_cell[va], lig_of_cell[va]]
        Ihat_val = Ihat_val - Ihat_val.mean()
        m_val = eval_metrics(Ihat_val, lat['I_cells'][va])
        # in-sample bound: random train holdout
        rng = np.random.default_rng(777 + s)
        hold = tr[rng.random(len(tr)) < 0.2]
        Ihat_hold = R4[row_of_cell[hold], lig_of_cell[hold]]
        Ihat_hold = Ihat_hold - Ihat_hold.mean()
        m_hold = eval_metrics(Ihat_hold, lat['I_cells'][hold])
        out['per_seed'][str(s)] = {'val_cells_dz': m_val['dead_zone_sign_accuracy'],
                                   'val_cells_sp': m_val['spearman'],
                                   'train_holdout_dz': m_hold['dead_zone_sign_accuracy'],
                                   'train_holdout_sp': m_hold['spearman']}
        print(f'seed{s} closed-form: val dz={m_val["dead_zone_sign_accuracy"]:.3f} '
              f'sp={m_val["spearman"]:.3f} | holdout dz={m_hold["dead_zone_sign_accuracy"]:.3f} '
              f'sp={m_hold["spearman"]:.3f}', flush=True)
    vdz = float(np.median([v['val_cells_dz'] for v in out['per_seed'].values()]))
    vsp = float(np.median([v['val_cells_sp'] for v in out['per_seed'].values()]))
    hdz = float(np.median([v['train_holdout_dz'] for v in out['per_seed'].values()]))
    hsp = float(np.median([v['train_holdout_sp'] for v in out['per_seed'].values()]))
    out['median'] = {'val_cells_dz': vdz, 'val_cells_sp': vsp,
                     'train_holdout_dz': hdz, 'train_holdout_sp': hsp}
    out['role'] = 'diagnostic upper bound only; deployment remains gradient-trained'
    out['preregistration_sha256'] = PREREG_SHA
    import json
    json.dump(out, open(HERE / 'Q2D1_CLOSED_FORM_DIAGNOSTIC.json', 'w'), indent=1)
    print(json.dumps(out['median']))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
