"""Stage Q2c-0: harness self-audit (projection negative controls, oracle
alignment recomputation, endpoint distortion, minimal linear no-censoring Q2,
truth separation). Imports stageX0c/q2.py READ-ONLY; nothing in stageX0c is
modified. Prereg: stageQ2c PREREGISTRATION.md (SHA 1027ccde...).
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
from q2 import (load_duongly_graph, generate, build_arm_feats, Q2Model,
                FreeIDModel, train, train_arm_with_restarts, predict,
                eval_metrics, anova_projection, censored_loss, BATCH, LR, WD,
                sha256_seed)
from x0_common import stable_rng

PREREG_SHA = '1027ccde8c8946aa8314ebd7642af89a6abbc3366afd965e8ab43f0da5a26a5c'
BOUND = 0.15  # near-random bound for projected pearson (frozen in prereg N1)


def perm_p(x, y, n_perm=200, seed=20260818):
    rng = np.random.default_rng(seed)
    def _c(a, b):
        if a.std() == 0 or b.std() == 0:
            return 0.0
        c = float(np.corrcoef(a, b)[0, 1])
        return 0.0 if np.isnan(c) else c
    obs = _c(x, y)
    hits = 0
    for _ in range(n_perm):
        yp = rng.permutation(y)
        if abs(_c(x, yp)) >= abs(obs):
            hits += 1
    return obs, (hits + 1) / (n_perm + 1)


def pearson(a, b):
    if a.std() == 0 or b.std() == 0:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def main():
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    t0 = time.time()
    (rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells,
     row_of_cell, lig_of_cell) = load_duongly_graph()
    n_prot = len(rows)
    print('graph loaded', round(time.time() - t0, 1), 's', flush=True)

    lat = generate(rows, compounds, prot_feats, lig_feats, cells, 1.0, 4, 'dense', False, 0)
    eval_c = splits['eval_cells']
    I_eval = lat['I_cells'][eval_c]
    train_c = splits['train_cells']

    # projection operator: TRAIN cells only (N5)
    proj = anova_projection(row_of_cell, lig_of_cell, train_c, lat['z_obs'], lat['determinate'])
    train_rows_set = set(row_of_cell[train_c].tolist())
    train_ligs_set = set(lig_of_cell[train_c].tolist())
    interior = np.asarray([(row_of_cell[eval_c][i] in train_rows_set and
                            lig_of_cell[eval_c][i] in train_ligs_set)
                           for i in range(len(eval_c))], dtype=bool)
    print('interior eval cells:', int(interior.sum()), '/', len(eval_c), flush=True)

    arms_wanted = ['correct_protein', 'ligand_only', 'no_interaction_head', 'oracle_protein']
    arm_feats = build_arm_feats(rows, prot_feats, lat)
    def save_partial():
        with open(HERE / 'q2c0_partial.json', 'w', encoding='utf-8') as f:
            json.dump({'results_so_far': results, 'stage': 'partial'}, f, indent=1)

    results = {}
    for arm in arms_wanted:
        t1 = time.time()
        model, val = train_arm_with_restarts(
            arm_feats[arm], arm, lat, 0, device, n_prot, prot_feats, lig_feats,
            row_of_cell, lig_of_cell, splits)
        out = predict(model, arm_feats[arm], lig_feats, row_of_cell[eval_c],
                      lig_of_cell[eval_c], device, arm)
        inter = out['inter']
        m = eval_metrics(inter, I_eval)
        if arm == 'no_interaction_head':
            # evidence: the scale parameter was NOT frozen during training
            results[arm + '_scale_drift'] = {'inter_scale_after_training': float(model.inter_scale.item()),
                                             'note': 'no_interaction_head trains with inter_scale=0 at init but the parameter keeps gradients; the artifact zeroed it only at predict time'}
        # N7: pre-projection (raw head) and post-projection (train-fitted operator)
        proj_inter = proj(inter, row_of_cell[eval_c], lig_of_cell[eval_c])
        pp_all, pp_p_all = perm_p(proj_inter, I_eval)
        pp_int = pearson(proj_inter[interior], I_eval[interior])
        results[arm] = {
            'raw_head': {'spearman': m['spearman'], 'pearson': m['pearson'],
                         'dead_zone_sign_accuracy': m['dead_zone_sign_accuracy']},
            'post_projection_all_eval': {'pearson': pp_all, 'perm_p': pp_p_all,
                                         'near_random_bound': BOUND,
                                         'passes_negative_control': abs(pp_all) < BOUND and pp_p_all > 0.05},
            'post_projection_interior_only': {'pearson': pp_int, 'n_cells': int(interior.sum())},
            'best_val_loss': val,
        }
        print(arm, 'raw_dz=', round(m['dead_zone_sign_accuracy'], 3),
              'proj_pearson=', round(pp_all, 3), 'p=', round(pp_p_all, 3), flush=True)
        save_partial()

    # free_target_id arm (uses FreeIDModel)
    fid = FreeIDModel(n_prot).to(device)
    _ep, val_fid = train(fid, prot_feats, lig_feats, row_of_cell, lig_of_cell,
                         train_c, lat, device, 0, 'free_target_id',
                         val_mask=splits['val_cells'])
    out_fid = predict(fid, prot_feats, lig_feats, row_of_cell[eval_c],
                      lig_of_cell[eval_c], device, 'free_target_id')
    inter_fid = out_fid['inter']
    m_fid = eval_metrics(inter_fid, I_eval)
    proj_fid = proj(inter_fid, row_of_cell[eval_c], lig_of_cell[eval_c])
    pp_fid, p_fid = perm_p(proj_fid, I_eval)
    results['free_target_id'] = {
        'raw_head': {'spearman': m_fid['spearman'], 'pearson': m_fid['pearson'],
                     'dead_zone_sign_accuracy': m_fid['dead_zone_sign_accuracy']},
        'post_projection_all_eval': {'pearson': pp_fid, 'perm_p': p_fid,
                                     'passes_negative_control': abs(pp_fid) < BOUND and p_fid > 0.05},
        'best_val_loss': val_fid,
    }
    print('free_target_id proj_pearson=', round(pp_fid, 3), flush=True)

    # N1: tau*=0 - projected correlation must be near-random
    lat0 = generate(rows, compounds, prot_feats, lig_feats, cells, 0.0, 4, 'dense', False, 0)
    I0_eval = lat0['I_cells'][eval_c]
    proj0 = anova_projection(row_of_cell, lig_of_cell, train_c, lat0['z_obs'], lat0['determinate'])
    model0, val0 = train_arm_with_restarts(
        arm_feats['correct_protein'], 'correct_protein', lat0, 0, device, n_prot,
        prot_feats, lig_feats, row_of_cell, lig_of_cell, splits)
    out0 = predict(model0, arm_feats['correct_protein'], lig_feats,
                   row_of_cell[eval_c], lig_of_cell[eval_c], device, 'correct_protein')
    proj0_int = proj0(out0['inter'], row_of_cell[eval_c], lig_of_cell[eval_c])
    pp0, p0 = perm_p(proj0_int, I0_eval)
    results['tau0'] = {'post_projection_all_eval': {'pearson': pp0, 'perm_p': p0,
                       'passes_negative_control': abs(pp0) < BOUND and p0 > 0.05}}
    print('tau0 proj_pearson=', round(pp0, 3), 'p=', round(p0, 3), flush=True)

    # N4: independent random bipartite graph with main effects only
    rng = stable_rng('stageQ2c', 'q2c0', 'rand_graph')
    n_p, n_l, n_obs = len(rows), len(compounds), len(cells)
    r_idx = rng.choice(n_p * n_l, size=n_obs, replace=False)
    r_rows = (r_idx // n_l).astype(np.int64)
    r_ligs = (r_idx % n_l).astype(np.int64)
    r_train = np.zeros(n_obs, dtype=bool); r_train[:len(train_c)] = True
    r_eval = np.zeros(n_obs, dtype=bool); r_eval[-len(eval_c):] = True
    pm_r = rng.normal(0, 1.0, size=n_p)
    lm_r = rng.normal(0, 1.0, size=n_l)
    noise_r = rng.normal(0, 1.0, size=n_obs)
    z_r = 0.5 + pm_r[r_rows] + lm_r[r_ligs] + noise_r
    det_r = np.ones(n_obs, dtype=bool)
    proj_r = anova_projection(r_rows, r_ligs, r_train, z_r, det_r)
    head_r = rng.normal(0, 1.0, size=n_obs)  # independent interaction-head noise
    I_r_truth = rng.normal(0, 1.0, size=n_obs)
    proj_r_head = proj_r(head_r[r_eval], r_rows[r_eval], r_ligs[r_eval])
    pp_r, p_r = perm_p(proj_r_head, I_r_truth[r_eval])
    results['random_bipartite_graph_null'] = {
        'post_projection_pearson': pp_r, 'perm_p': p_r,
        'passes_negative_control': abs(pp_r) < BOUND and p_r > 0.05}
    print('random-graph proj_pearson=', round(pp_r, 3), 'p=', round(p_r, 3), flush=True)

    # endpoint distortion: latent I vs I implied by observable scale
    det_e = lat['determinate'][eval_c]
    z_obs_e = lat['z_obs'][eval_c][det_e]
    I_implied = z_obs_e - (lat['mu'] + lat['pmain'][row_of_cell[eval_c][det_e]]
                           + lat['lmain'][lig_of_cell[eval_c][det_e]])
    epd = {'n_determinate_eval': int(det_e.sum()),
           'pearson_latent_vs_implied': pearson(lat['I_cells'][eval_c][det_e], I_implied),
           'spearman_latent_vs_implied': float(np.corrcoef(
               np.argsort(np.argsort(lat['I_cells'][eval_c][det_e])),
               np.argsort(np.argsort(I_implied)))[0, 1])}
    results['endpoint_distortion'] = epd
    print('endpoint distortion:', {k: round(v, 3) for k, v in epd.items() if isinstance(v, float)}, flush=True)

    # minimal linear no-censoring Q2 (identity link, linear head)
    L_dim = lig_feats.shape[1]
    n_lig = len(compounds)
    class LinearInter2(torch.nn.Module):
        def __init__(self, d_p, n_prot, d_l, n_lig):
            super().__init__()
            self.p_b = torch.nn.Parameter(torch.zeros(n_prot))
            self.l_b = torch.nn.Parameter(torch.zeros(n_lig))
            self.mu = torch.nn.Parameter(torch.zeros(1))
            self.Wp = torch.nn.Linear(d_p, 4)
            self.Wl = torch.nn.Linear(d_l, 4)
            self.inter_scale = torch.nn.Parameter(torch.ones(1))
        def forward(self, pe, le, rows=None, ligs=None):
            return (self.mu + self.p_b[rows] + self.l_b[ligs] +
                    self.inter_scale * ((self.Wp(pe)) * (self.Wl(le))).sum(dim=-1, keepdim=True))
    lin = LinearInter2(prot_feats.shape[1], n_prot, L_dim, n_lig).to(device)
    torch.manual_seed(0)
    opt = torch.optim.AdamW(lin.parameters(), lr=LR, weight_decay=WD)
    Pt = torch.from_numpy(prot_feats).float().to(device)
    Lt = torch.from_numpy(lig_feats).float().to(device)
    z_lin = (lat['mu'] + lat['pmain'][row_of_cell] + lat['lmain'][lig_of_cell]
             + lat['I_cells'] + np.random.default_rng(12345).normal(0, 1.0, size=len(cells)))
    zc = torch.from_numpy(z_lin[train_c]).float().to(device)
    rng_steps = np.random.default_rng(sha256_seed('stageQ2c', 'q2c0', 'lin_steps', 0, 'lin'))
    n_tr = len(train_c)
    for step in range(6000):
        idx = rng_steps.choice(n_tr, size=min(BATCH, n_tr), replace=False)
        r_i, l_i = row_of_cell[train_c][idx], lig_of_cell[train_c][idx]
        out_lin = lin(Pt[r_i], Lt[l_i], r_i, l_i)
        loss = ((out_lin.squeeze(-1) - zc[idx]) ** 2).mean()
        opt.zero_grad(); loss.backward(); opt.step()
    r_e = torch.from_numpy(row_of_cell[eval_c]).to(device)
    l_e = torch.from_numpy(lig_of_cell[eval_c]).to(device)
    with torch.no_grad():
        inter_lin = (lin.inter_scale * ((lin.Wp(Pt[r_e])) *
                     (lin.Wl(Lt[l_e]))).sum(dim=-1)).cpu().numpy()
    m_lin = eval_metrics(inter_lin, I_eval)
    proj_lin = proj(inter_lin, row_of_cell[eval_c], lig_of_cell[eval_c])
    results['minimal_linear_no_censoring'] = {
        'raw_head': {'spearman': m_lin['spearman'], 'pearson': m_lin['pearson'],
                     'dead_zone_sign_accuracy': m_lin['dead_zone_sign_accuracy']},
        'post_projection_pearson': pearson(proj_lin, I_eval)}
    print('minimal-linear raw_dz=', round(m_lin['dead_zone_sign_accuracy'], 3),
          'raw_sp=', round(m_lin['spearman'], 3), flush=True)

    # truth separation (item 5): save + sha
    truth = {'mu': float(lat['mu']), 'pmain': lat['pmain'], 'lmain': lat['lmain'],
             'I_cells': lat['I_cells'], 'z': lat['z'], 'y': lat['y'],
             'determinate': lat['determinate'].astype(np.uint8)}
    np.savez(HERE / 'q2c0_truth_seed0.npz', **truth)
    import hashlib
    h = hashlib.sha256((HERE / 'q2c0_truth_seed0.npz').read_bytes()).hexdigest()

    out = {
        'schema': 'MetaSieve.StageQ2c.Q2C0_PROJECTION_AUDIT.v1',
        'preregistration_sha256': PREREG_SHA,
        'stage': 'stageQ2c_harness_audit_20260818',
        'inputs': {'stageX0c_q2_py_read_only': True, 'gate_point': 'tau*=1.0, rank 4, dense, seed 0'},
        'negative_control_bound': BOUND,
        'arms': results,
        'n5_operator_train_only': 'verified: operator fit exclusively on train mask; eval labels never enter',
        'truth_separated_npz': 'q2c0_truth_seed0.npz',
        'truth_sha256': h,
        'conclusions': {
            'anova_projection_is_safe_diagnostic': all(
                results.get(a, {}).get('post_projection_all_eval', {}).get('passes_negative_control', False)
                for a in ['ligand_only', 'no_interaction_head', 'tau0']) and
                results.get('random_bipartite_graph_null', {}).get('passes_negative_control', False),
            'no_interaction_head_projected_pearson': results.get('no_interaction_head', {}).get('post_projection_all_eval', {}).get('pearson'),
            'oracle_alignment_note': 'see ORACLE_ALIGNMENT_TABLE.md; in-artifact oracle arm dz 0.607-0.674 below gate',
        },
    }
    (HERE / 'Q2C0_PROJECTION_AUDIT.json').write_text(json.dumps(out, indent=1))
    print(json.dumps(out['conclusions'], indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
