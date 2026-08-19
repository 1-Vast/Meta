"""Stage Q2d-1b ladder runner: interaction-only Phase A, then B sigmoid,
C missingness, D censoring (asserted >0), E main effects + competition.
Arms share identical minibatch order, init policy and checkpoint rule.
Prereg SHA: 872bc4402f228d940776e7efe2fee6b91e8310badb4e8830f653ca5e5d2e998e.
"""
import json, sys
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn

HERE = Path(__file__).resolve().parent
X0C = HERE.parent / "stageX0c_measurement_qualification_20260818"
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(HERE))
import q2
from q2 import eval_metrics, censored_loss
from x0_common import stable_rng
import truth

PREREG_SHA = "872bc4402f228d940776e7efe2fee6b91e8310badb4e8830f653ca5e5d2e998e"
SEEDS = [0, 1, 2]
RANK = 4
HID = 32
BATCH = 1024
LR = 5e-3
WD = 1e-4
TOTAL_STEPS = 6000


class InterOnly(nn.Module):
    """scale * ((p A) . (l B) + inter_bias). No main effects, no ID bias."""
    def __init__(self, d_p, d_l):
        super().__init__()
        self.A = nn.Linear(d_p, RANK)
        self.B = nn.Linear(d_l, RANK)
        self.inter_bias = nn.Parameter(torch.zeros(1))
        self.inter_scale = nn.Parameter(torch.ones(1))

    def forward(self, p, l):
        inter = self.inter_scale * ((self.A(p) * self.B(l)).sum(-1) + self.inter_bias)
        return {"yhat": inter, "inter": inter}


class AdditiveOnly(nn.Module):
    def __init__(self, n_rows, n_lig):
        super().__init__()
        self.mu = nn.Parameter(torch.zeros(1))
        self.p_b = nn.Parameter(torch.zeros(n_rows))
        self.l_b = nn.Parameter(torch.zeros(n_lig))

    def forward(self, p, l, rows, ligs):
        y = (self.mu + self.p_b[rows] + self.l_b[ligs]).squeeze(-1)
        z = torch.zeros_like(y)
        return {"yhat": y, "inter": z}


class WithMain(nn.Module):
    """Level E: shared linear encoders for pm/lm (competition) + ID biases."""
    def __init__(self, d_p, d_l, n_rows, n_lig):
        super().__init__()
        self.enc_p = nn.Linear(d_p, HID)
        self.enc_l = nn.Linear(d_l, HID)
        self.p_head = nn.Linear(HID, 1)
        self.l_head = nn.Linear(HID, 1)
        self.p_b = nn.Parameter(torch.zeros(n_rows))
        self.l_b = nn.Parameter(torch.zeros(n_lig))
        self.mu = nn.Parameter(torch.zeros(1))
        self.A = nn.Linear(HID, RANK)
        self.B = nn.Linear(HID, RANK)
        self.inter_bias = nn.Parameter(torch.zeros(1))
        self.inter_scale = nn.Parameter(torch.ones(1))

    def forward(self, p, l, rows, ligs):
        ep = self.enc_p(p)
        el = self.enc_l(l)
        pm = self.p_head(ep).squeeze(-1) + self.p_b[rows]
        lm = self.l_head(el).squeeze(-1) + self.l_b[ligs]
        inter = self.inter_scale * (((ep @ self.A.weight.T) * (el @ self.B.weight.T)).sum(-1)
                                    + self.inter_bias)
        return {"yhat": (self.mu + pm + lm + inter).squeeze(-1), "inter": inter}


def make_level_targets(t, level, cells):
    r = cells[:, 0]
    l = cells[:, 1]
    z = t["I"][r, l] + t["noise"][r, l]
    if level == "E":
        z = z + t["mu"] + t["pm"][r] + t["lm"][l]
    y = 100.0 / (1.0 + np.exp(-z))
    if level in ("A", "C"):
        n = len(cells)
        return z, y, z, np.ones(n, dtype=bool), np.zeros(n), np.zeros(n), 0
    if level == "B":
        n = len(cells)
        return z, y, z, np.ones(n, dtype=bool), np.zeros(n), np.zeros(n), 0
    y_int = np.round(np.clip(y, 0, 100))
    lo_b = np.log(0.5 / 99.5)
    hi_b = np.log(99.5 / 0.5)
    det = (y_int > 0) & (y_int < 100)
    z_obs = np.where(det, np.log(y_int / (100 - y_int)), np.nan)
    blo = np.where(det, 0.0, np.where(y_int <= 0, -np.inf, hi_b))
    bhi = np.where(det, 0.0, np.where(y_int >= 100, np.inf, lo_b))
    return z, y, z_obs, det, blo, bhi, int((~det).sum())


def train_level(P, arm, t, level, seed, splits, device, restart, Lt_dev):
    rng_steps = stable_rng("stageQ2d1b", "steps", "seed", seed, "phase", level,
                            "restart", restart)  # IDENTICAL across arms
    torch.manual_seed(restart)
    n_rows, n_lig = t["I"].shape
    if arm == "additive_only":
        model = AdditiveOnly(n_rows, n_lig).to(device)
    elif level == "E":
        model = WithMain(P.shape[1], 64, n_rows, n_lig).to(device)
    else:
        model = InterOnly(P.shape[1], 64).to(device)
    for p_ in model.parameters():
        if p_.dim() > 1:
            nn.init.xavier_uniform_(p_)
    if hasattr(model, "A"):
        with torch.no_grad():
            model.A.weight.mul_(0.5)
    if hasattr(model, "B") and hasattr(model, "inter_scale"):
        with torch.no_grad():
            model.B.weight.mul_(0.5)
    if arm == "no_interaction_head":
        model.inter_scale.requires_grad_(False)
        with torch.no_grad():
            model.inter_scale.fill_(0.0)
    opt = torch.optim.AdamW([p_ for p_ in model.parameters() if p_.requires_grad],
                            lr=LR, weight_decay=WD)
    Pt = torch.from_numpy(P).float().to(device)
    tr = splits["train_cells"]
    if level == "C":
        obs = stable_rng("stageQ2d1b", "missing").random(len(tr)) < 0.70
        tr = tr[obs]
    zt, yt, z_obs, det, blo, bhi, n_cens = make_level_targets(t, level, tr)
    if level in ("D", "E"):
        assert n_cens > 0, "frozen assertion: censored_count must be > 0"
    zc = torch.from_numpy(z_obs.astype(np.float32)).to(device)
    dc = torch.from_numpy(det).to(device)
    lo = torch.from_numpy(blo.astype(np.float32)).to(device)
    hi = torch.from_numpy(bhi.astype(np.float32)).to(device)
    yc = torch.from_numpy(yt.astype(np.float32)).to(device)
    zfull = torch.from_numpy(zt.astype(np.float32)).to(device)
    r_t = torch.from_numpy(tr[:, 0]).to(device)
    l_t = torch.from_numpy(tr[:, 1]).to(device)
    n = len(tr)

    def loss_fn(out_yhat, idx):
        if level in ("A", "C"):
            return ((out_yhat[idx] - zfull[idx]) ** 2).mean()
        if level == "B":
            yh = 100.0 * torch.sigmoid(out_yhat[idx])
            return ((yh - yc[idx]) ** 2).mean() / 100.0
        return censored_loss({"yhat": out_yhat[idx]}, zc[idx], dc[idx],
                              lo[idx], hi[idx], device)
    best = None
    best_ep = 0
    best_state = None
    for step in range(TOTAL_STEPS):
        idx = rng_steps.choice(n, size=min(BATCH, n), replace=False)
        idx = torch.from_numpy(idx).to(device)
        if arm == "additive_only":
            out = model(None, None, r_t[idx], l_t[idx])
        elif level == "E":
            out = model(Pt[r_t[idx]], Lt_dev[l_t[idx]], r_t[idx], l_t[idx])
        else:
            out = model(Pt[r_t[idx]], Lt_dev[l_t[idx]])
        loss = loss_fn(out["yhat"], torch.arange(len(idx), device=device))
        opt.zero_grad()
        loss.backward()
        opt.step()
        if step % 300 == 0 or step == TOTAL_STEPS - 1:
            model.eval()
            with torch.no_grad():
                if arm == "additive_only":
                    out_all = model(None, None, r_t, l_t)
                elif level == "E":
                    out_all = model(Pt[r_t], Lt_dev[l_t], r_t, l_t)
                else:
                    out_all = model(Pt[r_t], Lt_dev[l_t])
                mon = float(loss_fn(out_all["yhat"], torch.arange(n, device=device)))
            if best is None or mon < best - 1e-9:
                best = mon
                best_ep = step
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            model.train()
    model.load_state_dict(best_state)
    model.eval()
    return model, (best if best is not None else float("nan"))


def eval_arm(model, P, arm, t, level, splits, device, Lt_dev):
    out = {}
    for surf in ("pc", "lc", "dc"):
        cells = splits[surf]
        r = torch.from_numpy(cells[:, 0]).to(device)
        l = torch.from_numpy(cells[:, 1]).to(device)
        Pt = torch.from_numpy(P).float().to(device)
        with torch.no_grad():
            if arm == "additive_only":
                o = model(None, None, r, l)
            elif level == "E":
                o = model(Pt[r], Lt_dev[l], r, l)
            else:
                o = model(Pt[r], Lt_dev[l])
        m = eval_metrics(o["inter"].cpu().numpy(), t["I"][cells[:, 0], cells[:, 1]])
        out[surf] = {"dz": m["dead_zone_sign_accuracy"], "sp": m["spearman"],
                     "sign_acc": m["sign_accuracy"]}
    return out


def main():
    device = "cuda" if torch.cuda.is_available() else "cpu"
    import x0_i1
    rows, compounds, prot_feats, lig_feats, scaffolds, _meta = x0_i1.load_features()
    fz = np.load(HERE / "q2d1b_features.npz", allow_pickle=False)
    P_t = fz["P_t"].astype(np.float32)
    L_t = fz["L_t"].astype(np.float32)
    splits = json.loads((HERE / "Q2D1B_SPLITS.json").read_text(encoding="utf-8"))
    for k in ("train_cells", "pc", "lc", "dc"):
        splits[k] = np.asarray(splits[k], dtype=np.int64)
    for k in ("cold_row", "cold_lig", "train_row", "train_lig"):
        splits[k] = np.asarray(splits[k], dtype=bool)
    n_rows = P_t.shape[0]
    Lt_dev = torch.from_numpy(L_t).float().to(device)

    rng_arm = stable_rng("stageQ2d1b", "arms")
    shuf = rng_arm.permutation(n_rows)
    fams = np.asarray([q2.family_of_parent(x0_i1._parent_of(r)) for r in rows])
    fam_perm = np.arange(n_rows)
    for f in set(fams.tolist()):
        idx = np.where(fams == f)[0]
        fam_perm[idx] = idx[rng_arm.permutation(len(idx))]
    rand_p = rng_arm.normal(0, 1, size=(n_rows, truth.PROT_DIM)).astype(np.float32)
    arm_inputs = {}
    arm_inputs["correct"] = P_t
    arm_inputs["ligand_only"] = np.zeros_like(P_t)
    arm_inputs["additive_only"] = P_t
    arm_inputs["shuffled_protein"] = P_t[shuf]
    arm_inputs["family_preserving_shuffle"] = P_t[fam_perm]
    arm_inputs["random_protein"] = rand_p
    arm_inputs["no_interaction_head"] = P_t
    results = {}
    stored_A = None
    for level in ("A", "B", "C", "D", "E"):
        results[level] = {}
        for seed in SEEDS:
            t = truth.generate_truth("M1", seed, P_t, L_t, splits)
            P_oracle = (P_t.astype(np.float64) @ t["A"]).astype(np.float32)
            arm_inputs["oracle_diagnostic"] = P_oracle
            for arm in arm_inputs:
                P = arm_inputs[arm]
                n_rest = 8 if arm == "correct" else 1
                best_model = None
                best_val = None
                for r_ in range(n_rest):
                    model, val = train_level(P, arm, t, level, seed, splits, device, r_, Lt_dev)
                    if best_val is None or val < best_val:
                        best_val = val
                        best_model = model
                m = eval_arm(best_model, P, arm, t, level, splits, device, Lt_dev)
                results[level].setdefault(str(seed), {})[arm] = m
                print(level, seed, arm, {s_: round(m[s_]["dz"], 3) for s_ in m}, flush=True)
                if level == "A" and seed == 0 and arm == "correct":
                    with torch.no_grad():
                        r_ = torch.from_numpy(splits["dc"][:, 0]).to(device)
                        l_ = torch.from_numpy(splits["dc"][:, 1]).to(device)
                        stored_A = best_model(torch.from_numpy(P).float().to(device)[r_],
                                              Lt_dev[l_])["inter"].cpu().numpy().copy()
    json.dump({"schema": "MetaSieve.StageQ2d1b.LADDER.v1",
               "preregistration_sha256": PREREG_SHA, "results": results},
              open(HERE / "Q2D1B_LADDER.json", "w"), indent=1)
    print("ladder done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
