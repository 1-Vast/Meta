"""Stage Q2d-1e ladder runner: span-initialized protein map + mild L2 on
factor maps (lambda 1e-3). Same truth as Q2d-1d (truth_d imported; truth
streams unchanged). Ladder A-E with value-level reproduction checks, 8
arms, gate on double-cold, NC1/NC2 negative controls, M2/M3 sweep.
Prereg SHA: 61bc0cc50edcd40d581d16e67fd9fb8cef2729d3bbf58c4e9b831b727d0584f7.
Truth inputs: q2d1d_features.npz + Q2D1D_SPLITS.json (identical to Q2d-1d).
"""
import json, sys
from pathlib import Path
import numpy as np
import torch
import torch.nn as nn

HERE = Path(__file__).resolve().parent
X0C = HERE.parent / "stageX0c_measurement_qualification_20260818"
STAGE_D = HERE.parent / "stageQ2d1d_spanrestricted_interaction_20260818"
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(STAGE_D))
import q2
from q2 import eval_metrics, censored_loss
from x0_common import stable_rng
import truth_d as truth

PREREG_SHA = "61bc0cc50edcd40d581d16e67fd9fb8cef2729d3bbf58c4e9b831b727d0584f7"
L2_PEN = 1e-3
SEEDS = [0, 1, 2]
RANK = 4
HID = 32
BATCH = 1024
LR = 5e-3
WD = 1e-4
TOTAL_STEPS = 6000


class InterOnly(nn.Module):
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
    n = len(cells)
    if level in ("A", "C"):
        return z, y, z, np.ones(n, dtype=bool), np.zeros(n), np.zeros(n), 0
    if level == "B":
        return z, y, z, np.ones(n, dtype=bool), np.zeros(n), np.zeros(n), 0
    y_int = np.round(np.clip(y, 0, 100))
    lo_b = np.log(0.5 / 99.5)
    hi_b = np.log(99.5 / 0.5)
    det = (y_int > 0) & (y_int < 100)
    z_obs = np.where(det, np.log(y_int / (100 - y_int)), np.nan)
    blo = np.where(det, 0.0, np.where(y_int <= 0, -np.inf, hi_b))
    bhi = np.where(det, 0.0, np.where(y_int >= 100, np.inf, lo_b))
    return z, y, z_obs, det, blo, bhi, int((~det).sum())


def train_level(P, arm, t, level, seed, splits, device, restart, Lt_dev,
                Vsp=None, max_steps=TOTAL_STEPS):
    rng_steps = stable_rng("stageQ2d1e", "steps", "seed", seed, "phase", level,
                            "restart", restart)  # IDENTICAL across arms
    torch.manual_seed(restart)
    n_rows, n_lig = t["I"].shape
    if arm == "additive_only":
        model = AdditiveOnly(n_rows, n_lig).to(device)
    elif level == "E":
        model = WithMain(P.shape[1], 48, n_rows, n_lig).to(device)
    else:
        model = InterOnly(P.shape[1], 48).to(device)
    for p_ in model.parameters():
        if p_.dim() > 1:
            nn.init.xavier_uniform_(p_)
    if hasattr(model, "A"):
        with torch.no_grad():
            model.A.weight.mul_(0.5)
            if P.shape[1] == Vsp.shape[0] and Vsp is not None and arm != "oracle_diagnostic":
                # frozen span-initialized protein map: A = G @ Vsp^T, G xavier
                g = torch.empty(RANK, Vsp.shape[1])
                nn.init.xavier_uniform_(g)
                model.A.weight.copy_(g.to(model.A.weight.device)
                                     @ torch.from_numpy(Vsp.T).float().to(model.A.weight.device))
    if hasattr(model, "enc_p") and Vsp is not None and P.shape[1] == Vsp.shape[0]:
        with torch.no_grad():
            proj = torch.from_numpy((Vsp @ Vsp.T).astype(np.float32)).to(model.enc_p.weight.device)
            model.enc_p.weight.copy_(proj)
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
        obs = stable_rng("stageQ2d1e", "missing").random(len(tr)) < 0.70
        tr = tr[obs]
    zt, yt, z_obs, det, blo, bhi, n_cens = make_level_targets(t, level, tr)
    zc = torch.from_numpy(z_obs.astype(np.float32)).to(device)
    dc_t = torch.from_numpy(det).to(device)
    lo = torch.from_numpy(blo.astype(np.float32)).to(device)
    hi = torch.from_numpy(bhi.astype(np.float32)).to(device)
    yc = torch.from_numpy(yt.astype(np.float32)).to(device)
    zfull = torch.from_numpy(zt.astype(np.float32)).to(device)
    r_t = torch.from_numpy(tr[:, 0]).to(device)
    l_t = torch.from_numpy(tr[:, 1]).to(device)
    n = len(tr)

    def loss_fn(out_yhat, targets):
        pen = 0.0
        if hasattr(model, "A"):
            pen = L2_PEN * (model.A.weight.square().sum()
                            + model.B.weight.square().sum())
        if level in ("A", "C"):
            return ((out_yhat - targets[0]) ** 2).mean() + pen
        if level == "B":
            yh = 100.0 * torch.sigmoid(out_yhat)
            return ((yh - targets[1]) ** 2).mean() / 100.0 + pen
        return censored_loss({"yhat": out_yhat},
                              targets[2].cpu().numpy(), targets[3].cpu().numpy(),
                              targets[4].cpu().numpy(), targets[5].cpu().numpy(),
                              device) + pen
    best = None
    best_state = None
    for step in range(max_steps):
        idx = rng_steps.choice(n, size=min(BATCH, n), replace=False)
        idx = torch.from_numpy(idx).to(device)
        if arm == "additive_only":
            out = model(None, None, r_t[idx], l_t[idx])
        elif level == "E":
            out = model(Pt[r_t[idx]], Lt_dev[l_t[idx]], r_t[idx], l_t[idx])
        else:
            out = model(Pt[r_t[idx]], Lt_dev[l_t[idx]])
        loss = loss_fn(out["yhat"], (zfull[idx], yc[idx], zc[idx], dc_t[idx],
                                          lo[idx], hi[idx]))
        opt.zero_grad()
        loss.backward()
        opt.step()
        if step % 300 == 0 or step == max_steps - 1:
            model.eval()
            with torch.no_grad():
                if arm == "additive_only":
                    out_all = model(None, None, r_t, l_t)
                elif level == "E":
                    out_all = model(Pt[r_t], Lt_dev[l_t], r_t, l_t)
                else:
                    out_all = model(Pt[r_t], Lt_dev[l_t])
                mon = float(loss_fn(out_all["yhat"], (zfull, yc, zc, dc_t, lo, hi)))
            if best is None or mon < best - 1e-9:
                best = mon
                best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
            model.train()
    model.load_state_dict(best_state)
    model.eval()
    return model, (best if best is not None else float("nan")), n_cens


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


def build_arm_inputs(P_t, t0, shuf, fam_perm, rand_p):
    ai = {}
    ai["correct"] = P_t
    ai["ligand_only"] = np.zeros_like(P_t)
    ai["additive_only"] = P_t
    ai["shuffled_protein"] = P_t[shuf]
    ai["family_preserving_shuffle"] = P_t[fam_perm]
    ai["random_protein"] = rand_p
    ai["no_interaction_head"] = P_t
    ai["oracle_diagnostic"] = (P_t.astype(np.float64) @ t0["A"]).astype(np.float32)
    return ai


def run_seed_arms(P_t, L_t, splits, device, Lt_dev, arm_inputs, mech, level, seed,
                Vsp=None):
    t = truth.generate_truth(mech, seed, P_t, L_t, splits)
    res = {}
    n_cens_total = 0
    for arm in arm_inputs:
        P = arm_inputs[arm]
        n_rest = 8 if arm == "correct" else 1
        best_model = None
        best_val = None
        for r_ in range(n_rest):
            model, val, n_cens = train_level(P, arm, t, level, seed, splits, device, r_,
                                             Lt_dev, Vsp=Vsp)
            n_cens_total += n_cens
            if best_val is None or val < best_val:
                best_val = val
                best_model = model
        res[arm] = eval_arm(best_model, P, arm, t, level, splits, device, Lt_dev)
        print(mech, level, seed, arm, {s_: round(res[arm][s_]["dz"], 3) for s_ in res[arm]},
              flush=True)
    return res, n_cens_total


def main():
    device = "cuda" if torch.cuda.is_available() else "cpu"
    import x0_i1
    rows, compounds, prot_feats, lig_feats, scaffolds, _meta = x0_i1.load_features()
    fz = np.load(STAGE_D / "q2d1d_features.npz", allow_pickle=False)  # identical truth inputs
    P_t = fz["P_t"].astype(np.float32)
    L_t = fz["L_t"].astype(np.float32)
    splits = json.loads(open(STAGE_D / "Q2D1D_SPLITS.json", encoding="utf-8").read())  # identical truth inputs
    for k in ("train_cells", "pc", "lc", "dc"):
        splits[k] = np.asarray(splits[k], dtype=np.int64)
    for k in ("cold_row", "cold_lig", "train_row", "train_lig"):
        splits[k] = np.asarray(splits[k], dtype=bool)
    n_rows = P_t.shape[0]
    Lt_dev = torch.from_numpy(L_t).float().to(device)
    # frozen precondition: the Q2d-1d oracle precheck (same truth) must have
    # passed before any Q2d-1e training
    pre = json.loads(open(STAGE_D / "Q2D1D_ORACLE_PRECHECK.json", encoding="utf-8").read())
    assert pre["M1_identifiable_on_all_surfaces"] is True, "oracle precheck precondition unmet"
    Vsp, _ = truth._span_projection(P_t, splits)
    print("span basis rank:", Vsp.shape[1], flush=True)
    rng_arm = stable_rng("stageQ2d1e", "arms")
    shuf = rng_arm.permutation(n_rows)
    fams = np.asarray([q2.family_of_parent(x0_i1._parent_of(r)) for r in rows])
    fam_perm = np.arange(n_rows)
    for f in set(fams.tolist()):
        idx = np.where(fams == f)[0]
        fam_perm[idx] = idx[rng_arm.permutation(len(idx))]
    rand_p = rng_arm.normal(0, 1, size=(n_rows, P_t.shape[1])).astype(np.float32)
    results = {}
    cens = {}
    stored_A = {}
    for mech, levels in (("M1", ("A", "B", "C", "D", "E")),
                         ("M2", ("A",)), ("M3", ("A",)),
                         ("NC1", ("A",)), ("NC2", ("A",))):
        results[mech] = {}
        cens[mech] = {}
        for level in levels:
            results[mech][level] = {}
            cens[mech][level] = 0
            for seed in SEEDS:
                t0 = truth.generate_truth(mech, seed, P_t, L_t, splits)
                arm_inputs = build_arm_inputs(P_t, t0, shuf, fam_perm, rand_p)
                res, n_cens_total = run_seed_arms(P_t, L_t, splits, device, Lt_dev,
                                                  arm_inputs, mech, level, seed, Vsp=Vsp)
                results[mech][level][str(seed)] = res
                cens[mech][level] += n_cens_total
                if level in ("D", "E"):
                    assert cens[mech][level] > 0, "frozen assertion: censored_count > 0"
                if mech == "M1" and level == "A" and seed == 0:
                    stored_A = {arm: res[arm] for arm in arm_inputs}
    # value-level reproduction: rerun level-A config through the same code path
    # for each level (ingredients disabled) and compare to stored A (seed 0)
    repro = {}
    t0 = truth.generate_truth("M1", 0, P_t, L_t, splits)
    arm_inputs = build_arm_inputs(P_t, t0, shuf, fam_perm, rand_p)
    for level in ("B", "C", "D", "E"):
        resA, _ = run_seed_arms(P_t, L_t, splits, device, Lt_dev, arm_inputs,
                                "M1", "A", 0, Vsp=Vsp)
        ok = all(
            abs(resA[a][s_]["dz"] - stored_A[a][s_]["dz"]) < 1e-9
            for a in arm_inputs for s_ in ("pc", "lc", "dc"))
        repro[level] = bool(ok)
        print("repro", level, ok, flush=True)
    out = {"schema": "MetaSieve.StageQ2d1e.LADDER.v1",
           "preregistration_sha256": PREREG_SHA,
           "results": results, "censored_counts": cens,
           "repro_A_value_level": repro}
    json.dump(out, open(HERE / "Q2D1E_LADDER.json", "w"), indent=1)
    print("ladder done; cens:", cens, "repro:", repro)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
