"""Q2c-1b control: oracle_PU arm at tau*=2.0 (rank 4, dense), 3 seeds,
official 8-restart val-selection protocol. Distinguishes '0.70 unreachable
at any tau*' from 'unreachable at tau*=1.0 on this graph'.
"""
import json, sys
from pathlib import Path
import numpy as np
import torch

HERE = Path(__file__).resolve().parent
X0C = HERE.parent / 'stageX0c_measurement_qualification_20260818'
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(X0C.parent))
import q2
from q2 import (load_duongly_graph, generate, train_arm_with_restarts,
                predict, eval_metrics)

PREREG_SHA = '1027ccde8c8946aa8314ebd7642af89a6abbc3366afd965e8ab43f0da5a26a5c'


def main():
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    (rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells,
     row_of_cell, lig_of_cell) = load_duongly_graph()
    n_prot = len(rows)
    eval_c = splits['eval_cells']
    out = {'schema': 'MetaSieve.StageQ2c.Q2C1B_ORACLE_TAU2.v1',
           'preregistration_sha256': PREREG_SHA, 'per_seed': {}}
    for s in [0, 1, 2]:
        lat = generate(rows, compounds, prot_feats, lig_feats, cells, 2.0, 4, 'dense', False, s)
        P = (prot_feats.astype(np.float64) @ lat['U']).astype(np.float32)
        model, val = train_arm_with_restarts(
            P, 'correct_protein', lat, s, device, n_prot, P, lig_feats,
            row_of_cell, lig_of_cell, splits)
        m = eval_metrics(
            predict(model, P, lig_feats, row_of_cell[eval_c], lig_of_cell[eval_c],
                    device, 'correct_protein')['inter'],
            lat['I_cells'][eval_c])
        out['per_seed'][str(s)] = {'dz': m['dead_zone_sign_accuracy'],
                                   'sp': m['spearman'], 'sign_acc': m['sign_accuracy']}
        print(f'seed{s} tau2 oracle: dz={m["dead_zone_sign_accuracy"]:.3f} sp={m["spearman"]:.3f}', flush=True)
    dz = float(np.median([v['dz'] for v in out['per_seed'].values()]))
    sp = float(np.median([v['sp'] for v in out['per_seed'].values()]))
    out['median'] = {'dz': dz, 'sp': sp}
    out['conclusion'] = ('0.70 reachable at tau*=2.0 by the oracle' if dz >= 0.70
                         else '0.70 NOT reachable even at tau*=2.0 by the oracle on this graph')
    json.dump(out, open(HERE / 'Q2C1B_ORACLE_TAU2.json', 'w'), indent=1)
    print(json.dumps(out['median']), out['conclusion'])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
