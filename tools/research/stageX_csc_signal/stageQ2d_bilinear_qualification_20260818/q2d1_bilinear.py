"""Stage Q2d-1: isomorphic bilinear positive control with stepwise ladder.
A) z-scale identity link; B) + sigmoid % scale loss; C) + 70% panel
missingness; D) + interval censoring; E) + main-effect competition (shared
linear encoders). Arms: exact_bilinear, additive_only, ligand_only,
shuffled_protein, random_protein, oracle_latent, no_interaction_head.
Prereg: 4f7e80027b9b82564bf1ea262d360813f23ffe77f8378ca54b490cc6752fa3d1.
"""
import json, sys, time
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn

HERE = Path(__file__).resolve().parent
X0C = HERE.parent / 'stageX0c_measurement_qualification_20260818'
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(X0C.parent))
import q2
from q2 import (load_duongly_graph, generate, eval_metrics, censored_loss,
                BATCH, LR, WD, sha256_seed)
from x0_common import stable_rng

PREREG_SHA = '4f7e80027b9b82564bf1ea262d360813f23ffe77f8378ca54b490cc6752fa3d1'
SEEDS = [0, 1, 2]
RANK = 4
HID = 64
TOTAL_STEPS = 6000


class BilinearInter(nn.Module):
    """mu + p_b(row) + l_b(lig) + scale*((p A).(l B) + bias)."""
    def __init__(self, d_p, d_l, n_prot, n_lig):
        super().__init__()
        self.mu = nn.Parameter(torch.zeros(1))
        self.p_b = nn.Parameter(torch.zeros(n_prot))
        self.l_b = nn.Parameter(torch.zeros(n_lig))
        self.A = nn.Linear(d_p, RANK)
        self.B = nn.Linear(d_l, RANK)
        self.inter_bias = nn.Parameter(torch.zeros(1))
        self.inter_scale = nn.Parameter(torch.ones(1))

    def inter(self, p, l):
        return self.inter_scale * (((self.A(p)) * (self.B(l))).sum(-1) + self.inter_bias)

    def forward(self, p, l, rows, ligs):
        inter = self.inter(p, l)
        return {'yhat': (self.mu + self.p_b[rows] + self.l_b[ligs] + inter).squeeze(-1),
                'inter': inter}


class SharedEncBilinear(nn.Module):
    """E: pm/lm learned through linear encoders shared with the interaction
    projections (main-effect competition), no nonlinearity."""
    def __init__(self, d_p, d_l, n_prot, n_lig):
        super().__init__()
        self.enc_p = nn.Linear(d_p, HID)
        self.enc_l = nn.Linear(d_l, HID)
        self.p_head = nn.Linear(HID, 1)
        self.l_head = nn.Linear(HID, 1)
        self.p_b = nn.Parameter(torch.zeros(n_prot))
        self.l_b = nn.Parameter(torch.zeros(n_lig))
        self.mu = nn.Parameter(torch.zeros(1))
        self.A = nn.Linear(HID, RANK)
        self.B = nn.Linear(HID, RANK)
        self.inter_bias = nn.Parameter(torch.zeros(1))
        self.inter_scale = nn.Parameter(torch.ones(1))

    def forward(self, p, l, rows, ligs):
        ep, el = self.enc_p(p), self.enc_l(l)
        pm = self.p_head(ep).squeeze(-1) + self.p_b[rows]
        lm = self.l_head(el).squeeze(-1) + self.l_b[ligs]
        inter = self.inter_scale * (((ep @ self.A.weight.T) * (el @ self.B.weight.T)).sum(-1)
                                    + self.inter_bias)
        return {'yhat': (self.mu + pm + lm + inter).squeeze(-1), 'inter': inter}


def train_one(model, P, L, rows_t, ligs_t, mask, lat, device, seed, phase,
              val_mask=None, sigmoid_loss=False):
    torch.manual_seed(seed)
    for p_ in model.parameters():
        if p_.dim() > 1:
            nn.init.xavier_uniform_(p_)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    Pt = torch.from_numpy(P).float().to(device)
    Lt = torch.from_numpy(L).float().to(device)
    m = mask
    n = len(m)
    r_c, l_c = rows_t[m], ligs_t[m]
    zc, dc, lo, hi = (lat['z_obs'][m], lat['determinate'][m],
                      lat['bounds_lo'][m], lat['bounds_hi'][m])
    yc = lat['y'][m]
    v_rows = v_ligs = vz = vd = vlo = vhi = vy = None
    if val_mask is not None and len(val_mask):
        v_rows, v_ligs = rows_t[val_mask], ligs_t[val_mask]
        vz, vd, vlo, vhi = (lat['z_obs'][val_mask], lat['determinate'][val_mask],
                            lat['bounds_lo'][val_mask], lat['bounds_hi'][val_mask])
        vy = lat['y'][val_mask]

    def loss_fn(out, z_t, d_t, lo_t, hi_t, y_t):
        if sigmoid_loss:
            yhat = 100.0 * torch.sigmoid(out)
            yt = torch.from_numpy(np.asarray(y_t, dtype=np.float32)).to(device)
            return ((yhat - yt) ** 2).mean() / 100.0
        return censored_loss({'yhat': out}, z_t, d_t, lo_t, hi_t, device)

    rng = np.random.default_rng(sha256_seed('stageQ2d', 'q2d1', 'steps', seed, phase))
    best, best_ep, best_state = None, 0, None
    for step in range(TOTAL_STEPS):
        idx = rng.choice(n, size=min(BATCH, n), replace=False)
        out = model(Pt[r_c[idx]], Lt[l_c[idx]], r_c[idx], l_c[idx])
        loss = loss_fn(out['yhat'], zc[idx], dc[idx], lo[idx], hi[idx], yc[idx])
        opt.zero_grad(); loss.backward(); opt.step()
        if step % 300 == 0 or step == TOTAL_STEPS - 1:
            model.eval()
            with torch.no_grad():
                if v_rows is not None:
                    outv = model(Pt[v_rows], Lt[v_ligs], v_rows, v_ligs)
                    mon = float(loss_fn(outv['yhat'], vz, vd, vlo, vhi, vy))
                else:
                    out2 = model(Pt[r_c], Lt[l_c], r_c, l_c)
                    mon = float(loss_fn(out2['yhat'], zc, dc, lo, hi, yc))
            if best is None or mon < best - 1e-6:
                best, best_ep = mon, step
                best_state = {k: t.detach().clone() for k, t in model.state_dict().items()}
            model.train()
    model.load_state_dict(best_state)
    model.eval()
    return best


def run_arm(P, arm, lat, phase, device, n_prot, n_lig, prot_feats, lig_feats,
            row_of_cell, lig_of_cell, train_mask, val_mask, eval_mask, seed,
            sigmoid_loss):
    n_restarts = 8 if arm == 'exact_bilinear' else 1
    best_model, best_val = None, None
    for r in range(n_restarts):
        model = (SharedEncBilinear if phase == 'E' else BilinearInter)(
            P.shape[1], lig_feats.shape[1], n_prot, n_lig).to(device)
        if arm == 'additive_only':
            model = _additive_only(P.shape[1], lig_feats.shape[1], n_prot, n_lig).to(device)
        if arm == 'no_interaction_head':
            model.inter_scale.requires_grad_(False)
            with torch.no_grad():
                model.inter_scale.fill_(0.0)
        val = train_one(model, P, lig_feats, row_of_cell, lig_of_cell, train_mask,
                        lat, device, r, phase, val_mask=val_mask, sigmoid_loss=sigmoid_loss)
        if best_val is None or val < best_val:
            best_val, best_model = val, model
    Pt = torch.from_numpy(P).float().to(device)
    Lt = torch.from_numpy(lig_feats).float().to(device)
    rq = torch.from_numpy(row_of_cell[eval_mask]).to(device)
    lq = torch.from_numpy(lig_of_cell[eval_mask]).to(device)
    with torch.no_grad():
        inter = best_model(Pt[rq], Lt[lq], rq, lq)['inter'].cpu().numpy()
    return eval_metrics(inter, lat['I_cells'][eval_mask])


def _additive_only(d_p, d_l, n_prot, n_lig):
    class AdditiveOnly(nn.Module):
        def __init__(self):
            super().__init__()
            self.mu = nn.Parameter(torch.zeros(1))
            self.p_b = nn.Parameter(torch.zeros(n_prot))
            self.l_b = nn.Parameter(torch.zeros(n_lig))
        def forward(self, p, l, rows, ligs):
            y = (self.mu + self.p_b[rows] + self.l_b[ligs]).squeeze(-1)
            z = torch.zeros_like(y)
            return {'yhat': y, 'inter': z}
    return AdditiveOnly()


def main():
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    (rows, compounds, prot_feats, lig_feats, scaffolds, splits, cells,
     row_of_cell, lig_of_cell) = load_duongly_graph()
    n_prot, n_lig = len(rows), len(compounds)
    train_c = splits['train_cells']
    val_c = splits['val_cells']
    eval_c = splits['eval_cells']

    arm_inputs = {}
    rng_missing = stable_rng('stageQ2d', 'q2d1', 'missingness')
    rng_rnd = stable_rng('stageQ2d', 'q2d1', 'random_protein')
    rng_shuf = stable_rng('stageQ2d', 'q2d1', 'shuffle_rows')
    shuf_perm = rng_shuf.permutation(n_prot)
    rand_p = rng_rnd.normal(0, 1, size=prot_feats.shape).astype(np.float32)

    phases = {}
    for phase, tag in [('A', 'z_identity'), ('B', 'sigmoid_pct'), ('C', 'missing70'),
                       ('D', 'censored'), ('E', 'competition')]:
        phases[phase] = {'summary': {}, 'per_seed': {}}

    for s in SEEDS:
        lat = generate(rows, compounds, prot_feats, lig_feats, cells, 1.0, 4, 'dense', False, s)
        latA = dict(lat)
        latA['z_obs'] = lat['z']
        latA['determinate'] = np.ones(len(cells), dtype=bool)
        latA['bounds_lo'] = np.zeros(len(cells))
        latA['bounds_hi'] = np.zeros(len(cells))
        latB = dict(latA)  # sigmoid %-scale loss uses lat['y'] target
        latB['y'] = lat['y']
        # C: 70% observed MCAR (train+val+eval all masked by same pattern)
        obs = rng_missing.random(len(cells)) < 0.70
        trainC = train_c[obs[train_c]]
        valC = val_c[obs[val_c]]
        evalC = eval_c[obs[eval_c]]
        latC = dict(latA)
        # D: censored floor semantics
        latD = dict(lat)  # q2.generate already produced censored z_obs/bounds

        for phase, (lat_p, tr, va, ev, sig) in [
            ('A', (latA, train_c, val_c, eval_c, False)),
            ('B', (latB, train_c, val_c, eval_c, True)),
            ('C', (latC, trainC, valC, evalC, False)),
            ('D', (latD, train_c, val_c, eval_c, False)),
            ('E', (latD, train_c, val_c, eval_c, False)),
        ]:
            if phase == 'C' and ev.sum() == 0:
                phases[phase]['per_seed'][str(s)] = {'note': 'no eval cells'}
                continue
            # oracle latent uses this seed's U
            oracle_p = (prot_feats.astype(np.float64) @ lat['U']).astype(np.float32)
            arm_inputs = {
                'exact_bilinear': prot_feats,
                'additive_only': prot_feats,
                'ligand_only': np.zeros_like(prot_feats),
                'shuffled_protein': prot_feats[shuf_perm],
                'random_protein': rand_p,
                'oracle_latent': oracle_p,
                'no_interaction_head': prot_feats,
            }
            res = {}
            for arm, P in arm_inputs.items():
                m = run_arm(P, arm, lat_p, phase, device, n_prot, n_lig,
                            prot_feats, lig_feats, row_of_cell, lig_of_cell,
                            tr, va, ev, s, sig)
                res[arm] = {'dz': m['dead_zone_sign_accuracy'], 'sp': m['spearman'],
                            'sign_acc': m['sign_accuracy']}
                print(f'{phase} seed{s} {arm}: dz={m["dead_zone_sign_accuracy"]:.3f} '
                      f'sp={m["spearman"]:.3f}', flush=True)
            phases[phase]['per_seed'][str(s)] = res
            with open(HERE / 'q2d1_partial.json', 'w') as f:
                json.dump(phases, f, indent=1)

    for phase in phases:
        ps = phases[phase]['per_seed']
        arms = ['exact_bilinear', 'additive_only', 'ligand_only', 'shuffled_protein',
                'random_protein', 'oracle_latent', 'no_interaction_head']
        for arm in arms:
            vals = [ps[str(s)].get(arm) for s in SEEDS if ps.get(str(s)) and ps[str(s)].get(arm)]
            vals = [v for v in vals if v]
            if not vals:
                continue
            dz = float(np.median([v['dz'] for v in vals]))
            sp = float(np.median([v['sp'] for v in vals]))
            gaps = []
            for s in SEEDS:
                ent = ps.get(str(s), {})
                if ent.get(arm) and ent.get('ligand_only'):
                    gaps.append(ent[arm]['sign_acc'] - ent['ligand_only']['sign_acc'])
            gap = float(np.median(gaps)) if gaps else None
            phases[phase]['summary'][arm] = {'dz_median': dz, 'sp_median': sp, 'gap_median': gap}
        eb = phases[phase]['summary'].get('exact_bilinear', {})
        neg_ok = all(
            not (phases[phase]['summary'].get(a, {}).get('sp_median', 0) >= 0.30 and
                 phases[phase]['summary'].get(a, {}).get('dz_median', 0) >= 0.70)
            for a in ['additive_only', 'ligand_only', 'shuffled_protein',
                      'random_protein', 'no_interaction_head'])
        phases[phase]['gate'] = {
            'exact_bilinear_pass': bool(eb.get('sp_median', 0) >= 0.30 and
                                        eb.get('dz_median', 0) >= 0.70 and
                                        (eb.get('gap_median') or 0) >= 0.05),
            'negative_arms_all_fail': neg_ok,
        }

    out = {
        'schema': 'MetaSieve.StageQ2d.Q2D1_BILINEAR.v1',
        'preregistration_sha256': PREREG_SHA,
        'phases': phases,
        'ladder_decision': 'see Q2D1_REPORT.md',
    }
    json.dump(out, open(HERE / 'Q2D1_BILINEAR.json', 'w'), indent=1)
    print(json.dumps({p: {'summary': phases[p]['summary'], 'gate': phases[p]['gate']}
                      for p in phases}, indent=1))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
