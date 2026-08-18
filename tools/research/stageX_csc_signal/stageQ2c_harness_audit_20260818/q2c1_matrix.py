"""Stage Q2c-1: representation x learner matrix at the frozen gate point
(tau*=1.0, rank 4, dense, seeds 0-2). Representations: one_hot_pocket,
pocket_esm (Q2C_FEATS), oracle_PU, random, shuffled. Learners: linear
(1 restart, val-selection) and mlp (official 8-restart val-selection
protocol). one_hot/mlp and oracle/mlp cells reuse the frozen X0c artifact
values (identical protocol) and are marked as such.
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
                predict, eval_metrics)
from q2c_lib import LinearInter, train_linear, eval_linear
from x0_common import stable_rng

PREREG_SHA = '1027ccde8c8946aa8314ebd7642af89a6abbc3366afd965e8ab43f0da5a26a5c'
SEEDS = [0, 1, 2]


def main():
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    t0 = time.time()
    (rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells,
     row_of_cell, lig_of_cell) = load_duongly_graph()
    n_prot, n_lig = len(rows), len(compounds)
    fz = np.load(HERE / 'q2c_row_esm.npz', allow_pickle=False)
    pocket_esm = fz['X'].astype(np.float32)
    rng_shuf = stable_rng('stageQ2c', 'q2c1', 'shuffle_rows')
    shuf_perm = rng_shuf.permutation(pocket_esm.shape[0])
    rng_rnd = stable_rng('stageQ2c', 'q2c1', 'random_rep')
    random_rep = rng_rnd.normal(0, 1, size=(n_prot, pocket_esm.shape[1])).astype(np.float32)
    eval_c = splits['eval_cells']

    artifact = json.loads((X0C / 'Q2_PLANTED.json').read_text(encoding='utf-8'))
    art = {}
    for s in SEEDS:
        k = f'1.0,4,dense,seed{s}'
        arms = artifact['results'][k]['arms']
        art[s] = {a: arms[a]['interaction_head'] for a in arms}

    matrix_cells = {}
    for rep in ['one_hot_pocket', 'pocket_esm', 'oracle_PU', 'random', 'shuffled']:
        matrix_cells[rep] = {}
        for learner in ['linear', 'mlp']:
            matrix_cells[rep][learner] = {}
        for s in SEEDS:
            lat = generate(rows, compounds, prot_feats, lig_feats, cells, 1.0, 4, 'dense', False, s)
            I_eval = lat['I_cells'][eval_c]
            if rep == 'one_hot_pocket':
                P = prot_feats
            elif rep == 'pocket_esm':
                P = pocket_esm
            elif rep == 'oracle_PU':
                P = (prot_feats.astype(np.float64) @ lat['U']).astype(np.float32)
            elif rep == 'random':
                P = random_rep
            else:
                P = pocket_esm[shuf_perm]
            for learner in ['linear', 'mlp']:
                if learner == 'mlp':
                    if rep in ('one_hot_pocket', 'oracle_PU') and s in art:
                        # reuse frozen artifact (identical protocol: 8 restarts, val selection)
                        src = 'correct_protein' if rep == 'one_hot_pocket' else 'oracle_protein'
                        m = art[s][src]
                        matrix_cells[rep][learner][s] = {'source': 'Q2_PLANTED.json ' + src,
                                                  'dz': m['dead_zone_sign_accuracy'],
                                                  'sp': m['spearman'],
                                                  'sign_acc': m['sign_accuracy']}
                        continue
                    model, val = train_arm_with_restarts(
                        P, 'correct_protein', lat, s, device, n_prot, P, lig_feats,
                        row_of_cell, lig_of_cell, splits)
                    out = predict(model, P, lig_feats, row_of_cell[eval_c],
                                  lig_of_cell[eval_c], device, 'correct_protein')
                    m = eval_metrics(out['inter'], I_eval)
                    matrix_cells[rep][learner][s] = {'source': 'trained here (8 restarts, val selection)',
                                              'dz': m['dead_zone_sign_accuracy'],
                                              'sp': m['spearman'],
                                              'sign_acc': m['sign_accuracy']}
                else:
                    lin = LinearInter(P.shape[1], n_prot, lig_feats.shape[1], n_lig).to(device)
                    _ep, val = train_linear(lin, P, lig_feats, row_of_cell, lig_of_cell,
                                            splits['train_cells'], lat, device, s,
                                            val_mask=splits['val_cells'])
                    m, _inter = eval_linear(lin, P, lig_feats, row_of_cell[eval_c],
                                            lig_of_cell[eval_c], I_eval, device)
                    matrix_cells[rep][learner][s] = {'source': 'trained here (linear, 1 restart, val selection)',
                                              'dz': m['dead_zone_sign_accuracy'],
                                              'sp': m['spearman'],
                                              'sign_acc': m['sign_accuracy']}
            # identity-link sensitivity branch (z-scale, no sigmoid, no censoring)
            if rep in ('one_hot_pocket', 'oracle_PU', 'pocket_esm'):
                lat_id = dict(lat)
                lat_id['z_obs'] = lat['z']
                lat_id['determinate'] = np.ones(len(cells), dtype=bool)
                lat_id['bounds_lo'] = np.zeros(len(cells))
                lat_id['bounds_hi'] = np.zeros(len(cells))
                model, val = train_arm_with_restarts(
                    P, 'correct_protein', lat_id, s, device, n_prot, P, lig_feats,
                    row_of_cell, lig_of_cell, splits)
                out = predict(model, P, lig_feats, row_of_cell[eval_c],
                              lig_of_cell[eval_c], device, 'correct_protein')
                m = eval_metrics(out['inter'], I_eval)
                matrix_cells[rep].setdefault('mlp_identity', {})[s] = {
                    'source': 'trained here on z-scale identity link (8 restarts, val selection)',
                    'dz': m['dead_zone_sign_accuracy'], 'sp': m['spearman'],
                    'sign_acc': m['sign_accuracy']}
            print(rep, s, 'done', flush=True)
            json.dump({'cells_so_far': {r: {l: {str(k2): v2 for k2, v2 in v.items()}
                                            for l, v in vv.items()} for r, vv in matrix_cells.items()}},
                      open(HERE / 'q2c1_partial.json', 'w'), indent=1)

    # ligand_only baseline for gap (from frozen artifact)
    lonly = {s: art[s]['ligand_only']['sign_accuracy'] for s in SEEDS}
    summary = {}
    for rep in matrix_cells:
        summary[rep] = {}
        for learner in matrix_cells[rep]:
            dzs = [matrix_cells[rep][learner][s]['dz'] for s in SEEDS]
            sps = [matrix_cells[rep][learner][s]['sp'] for s in SEEDS]
            gaps = [matrix_cells[rep][learner][s]['sign_acc'] - lonly[s] for s in SEEDS]
            summary[rep][learner] = {
                'dz_median': float(np.median(dzs)), 'sp_median': float(np.median(sps)),
                'gap_median': float(np.median(gaps)), 'per_seed': matrix_cells[rep][learner],
            }
    out = {
        'schema': 'MetaSieve.StageQ2c.Q2C1_MATRIX.v1',
        'preregistration_sha256': PREREG_SHA,
        'gate_point': 'tau*=1.0, rank 4, dense',
        'ligand_only_sign_acc_by_seed': {str(s): lonly[s] for s in SEEDS},
        'summary': summary,
        'interpretation_rules': ('oracle-linear pass + e2e fail = optimization/routing; '
                                 'oracle pass + pocket_esm fail = representation gap; '
                                 'all probes fail = harness definition/graph power/truth; '
                                 'all pass = proceed to Q2c-2'),
        'budget_note': 'mlp learner = official 8-restart val-selection protocol; linear learner = 1 restart x 6000 steps with val selection',
    }
    json.dump(out, open(HERE / 'Q2C1_MATRIX.json', 'w'), indent=1)
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
