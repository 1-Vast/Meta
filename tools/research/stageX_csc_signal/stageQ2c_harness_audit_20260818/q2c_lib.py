"""Stage Q2c shared library: linear interaction model + training/eval
helpers for the Q2c-1 matrix. Frozen with the Q2c preregistration."""
import numpy as np
import torch

from x0_common import sha256_seed

LR = 5e-3
WD = 1e-4
BATCH = 1024
HEAD_RANK = 4
TOTAL_STEPS = 6000


class LinearInter(torch.nn.Module):
    """Linear interaction learner: mu + p_b(row) + l_b(lig) + scale * sum_k
    ((Wp pe)_k * (Wl le)_k). All encoders linear; no nonlinearity."""
    def __init__(self, d_p, n_prot, d_l, n_lig, rank=HEAD_RANK):
        super().__init__()
        self.mu = torch.nn.Parameter(torch.zeros(1))
        self.p_b = torch.nn.Parameter(torch.zeros(n_prot))
        self.l_b = torch.nn.Parameter(torch.zeros(n_lig))
        self.Wp = torch.nn.Linear(d_p, rank)
        self.Wl = torch.nn.Linear(d_l, rank)
        self.inter_scale = torch.nn.Parameter(torch.ones(1))

    def forward(self, pe, le, rows=None, ligs=None):
        return (self.mu + self.p_b[rows] + self.l_b[ligs] +
                self.inter_scale * ((self.Wp(pe)) * (self.Wl(le))).sum(dim=-1, keepdim=True))


def train_linear(model, P, L, rows_t, ligs_t, mask, lat, device, seed,
                 val_mask=None, val_loss_only=False):
    """Train LinearInter with censored loss on z_obs; returns (best_ep, best_val).
    val_loss_only=True trains on train cells and selects by val loss (same
    selection discipline as the harness); False trains on the given mask only.
    """
    import q2
    torch.manual_seed(seed)
    for p in model.parameters():
        if p.dim() > 1:
            torch.nn.init.xavier_uniform_(p)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WD)
    Pt = torch.from_numpy(P).float().to(device)
    Lt = torch.from_numpy(L).float().to(device)
    m = mask
    n = len(m)
    r_c, l_c = rows_t[m], ligs_t[m]
    zc, dc, lo, hi = (lat['z_obs'][m], lat['determinate'][m],
                      lat['bounds_lo'][m], lat['bounds_hi'][m])
    v_rows, v_ligs = None, None
    if val_mask is not None and len(val_mask):
        v_rows, v_ligs = rows_t[val_mask], ligs_t[val_mask]
        vz, vd, vlo, vhi = (lat['z_obs'][val_mask], lat['determinate'][val_mask],
                            lat['bounds_lo'][val_mask], lat['bounds_hi'][val_mask])
    rng = np.random.default_rng(sha256_seed('stageQ2c', 'q2c1', 'lin_steps', seed, 'lin'))
    best, best_ep, best_state = None, 0, None
    for step in range(TOTAL_STEPS):
        idx = rng.choice(n, size=min(BATCH, n), replace=False)
        out = model(Pt[r_c[idx]], Lt[l_c[idx]], r_c[idx], l_c[idx])
        loss = q2.censored_loss({'yhat': out}, zc[idx], dc[idx], lo[idx], hi[idx], device)
        opt.zero_grad(); loss.backward(); opt.step()
        if step % 300 == 0 or step == TOTAL_STEPS - 1:
            model.eval()
            with torch.no_grad():
                if v_rows is not None:
                    outv = model(Pt[v_rows], Lt[v_ligs], v_rows, v_ligs)
                    mon = float(q2.censored_loss({'yhat': outv}, vz, vd, vlo, vhi, device))
                else:
                    out2 = model(Pt[r_c], Lt[l_c], r_c, l_c)
                    mon = float(q2.censored_loss({'yhat': out2}, zc, dc, lo, hi, device))
            if best is None or mon < best - 1e-6:
                best, best_ep = mon, step
                best_state = {k: t.detach().clone() for k, t in model.state_dict().items()}
    model.load_state_dict(best_state)
    model.eval()
    return best_ep, (best if best is not None else float('nan'))


def eval_linear(model, P, L, rows_q, ligs_q, I_truth, device):
    import q2
    Pt = torch.from_numpy(P).float().to(device)
    Lt = torch.from_numpy(L).float().to(device)
    rq = torch.from_numpy(np.asarray(rows_q, dtype=np.int64)).to(device)
    lq = torch.from_numpy(np.asarray(ligs_q, dtype=np.int64)).to(device)
    with torch.no_grad():
        inter = (model.inter_scale *
                 ((model.Wp(Pt[rq])) * (model.Wl(Lt[lq]))).sum(dim=-1)).cpu().numpy()
    return q2.eval_metrics(inter, I_truth), inter
