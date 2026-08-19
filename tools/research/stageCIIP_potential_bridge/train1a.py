"""Stage CIIP-1A trainer (prereg 31d3eeaf...).

All arms share the frozen DATA1A contract (matched rows, split, seeds,
budget). Unified potential f = b_P + b_L + s (s = alpha^T psi) is the
only promotable arm. Free pairwise arms are diagnostic only. Endpoint
stays percent inhibition. End-to-end gradient training only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import potential as POT  # noqa: E402

SIG = HERE.parent / "stageX_csc_signal"
sys.path.insert(0, str(SIG))
from x0_common import stable_rng  # noqa: E402

PREREG_SHA = "31d3eeaf6a0d77c46b3bbbee0fe9d2ff667aadeeb7d9dcabd26ca59ec48d5196"
EPOCHS = 200
BATCH = 512
LR = 1e-3
WD = 1e-4
LAMBDA_ABS = 1.0
DEAD_ZONE = 10.0
BOOT_DRAWS = 2000
BOOT_SEED = 20260820
SCREENING_ARMS = ("unified_local", "ligand_only", "family_shuffle",
                  "free_pairwise")
ALL_ARMS = ("unified_local", "unified_global", "ligand_only",
            "no_interaction", "family_shuffle", "random_protein",
            "ligand_invariant_shift", "free_pairwise", "free_ligand_pair")


class LigandInvariantShift(nn.Module):
    """Per-row scalar mutation shift; the centered contrast is exactly 0."""

    def __init__(self, n_rows, d_p=POT.D_P, d_l=POT.D_L, hid=POT.HID):
        super().__init__()
        self.p_enc = nn.Linear(d_p, hid)
        self.l_enc = nn.Linear(d_l, hid)
        self.b_P = nn.Linear(hid, 1)
        self.b_L = nn.Linear(hid, 1)
        self.a = nn.Parameter(torch.zeros(n_rows))
        self.mu = nn.Parameter(torch.zeros(1))

    def forward(self, P, L, row_idx):
        ep = torch.relu(self.p_enc(P))
        el = torch.relu(self.l_enc(L))
        return (self.mu + self.b_P(ep).squeeze(-1)
                + self.b_L(el).squeeze(-1) + self.a[row_idx])

    def centered_mutation_effect(self, Pw, Pv, Lm):
        return torch.zeros(Lm.shape[0], device=Lm.device)


def load_data():
    art = json.loads((HERE / "DATA1A.json").read_text(encoding="utf-8"))
    z = np.load(HERE / "DATA1A.npz", allow_pickle=False)
    return art, z["Y"], z["prot"], z["lig"], z["pair_split"]


def spearman(a, b):
    from scipy.stats import spearmanr
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    return float(spearmanr(a, b).correlation)


def pearson(a, b):
    from scipy.stats import pearsonr
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        return float("nan")
    return float(pearsonr(a, b)[0])


def dead_zone_sign_acc(a, b, dz=DEAD_ZONE):
    m = np.abs(b) >= dz
    if m.sum() == 0:
        return float("nan")
    return float(np.mean((a[m] > 0) == (b[m] > 0)))


def ols_scale(a, b):
    denom = float((b * b).sum())
    if denom == 0:
        return float("nan")
    return float((a * b).sum() / denom)


def make_arm(name, device, n_rows):
    torch.manual_seed(1)
    if name == "unified_local":
        m = POT.UnifiedPotential().to(device)
    elif name == "unified_global":
        m = POT.GlobalPotential().to(device)
    elif name == "ligand_only":
        m = POT.UnifiedPotential().to(device)
        with torch.no_grad():
            for mm in (m.alpha, m.psi, m.b_P):
                for p_ in mm.parameters():
                    p_.zero_()
                    p_.requires_grad_(False)
    elif name == "no_interaction":
        m = POT.UnifiedPotential().to(device)
        with torch.no_grad():
            for mm in (m.alpha, m.psi):
                for p_ in mm.parameters():
                    p_.zero_()
                    p_.requires_grad_(False)
    elif name in ("family_shuffle", "random_protein"):
        m = POT.UnifiedPotential().to(device)
    elif name == "ligand_invariant_shift":
        m = LigandInvariantShift(n_rows).to(device)
        for p_ in m.parameters():
            if p_.dim() > 1:
                nn.init.xavier_uniform_(p_)
    elif name == "free_pairwise":
        m = POT.FreePairwise().to(device)
    elif name == "free_ligand_pair":
        m = POT.FreeLigandPair().to(device)
    else:
        raise ValueError(name)
    return m


def shuffle_rows(prot, art, rng, family_only):
    perm = np.arange(prot.shape[0])
    if family_only:
        by_parent = {}
        for p in art["pairs"]:
            by_parent.setdefault(p["parent"], []).append(p["var_row"])
        for parent in sorted(by_parent):
            idx = sorted(set(by_parent[parent]))
            sh = idx.copy()
            rng.shuffle(sh)
            for a, b in zip(idx, sh):
                perm[a] = b
    else:
        rng.shuffle(perm)
    return prot[perm].copy()


def pair_hats(model, art, targets, idx, pairs, Pt, Lt):
    """Per-pair centered mutation effect from the SAME s_theta."""
    hats = {}
    with torch.no_grad():
        for i in idx:
            t = targets[i]
            Lm = Lt[t["lig_idx"]]
            Pw = Pt[pairs[i]["wt_row"]:pairs[i]["wt_row"] + 1].expand(
                len(t["lig_idx"]), -1)
            Pv = Pt[pairs[i]["var_row"]:pairs[i]["var_row"] + 1].expand(
                len(t["lig_idx"]), -1)
            hats[i] = model.centered_mutation_effect(Pw, Pv, Lm)
    return hats


def train_arm(name, art, Y, prot, lig, split, device, seed, dry=False):
    torch.manual_seed(seed)
    pairs = art["pairs"]
    targets = art["targets"]
    tr_idx = [i for i, s in enumerate(split) if s == 0]
    va_idx = [i for i, s in enumerate(split) if s == 1]
    te_idx = [i for i, s in enumerate(split) if s == 2]
    rng_rows = stable_rng("stageCIIP1A", "rows", name, seed)
    if name == "family_shuffle":
        P = shuffle_rows(prot, art, rng_rows, True)
    elif name == "random_protein":
        P = shuffle_rows(prot, art, rng_rows, False)
    else:
        P = prot
    Pt = torch.from_numpy(P).float().to(device)
    Lt = torch.from_numpy(lig).float().to(device)
    Yt = torch.from_numpy(Y).float().to(device)
    model = make_arm(name, device, P.shape[0])
    opt = torch.optim.AdamW([p for p in model.parameters() if p.requires_grad],
                            lr=LR, weight_decay=WD)
    va_te_rows = ({pairs[i]["wt_row"] for i in va_idx + te_idx}
                  | {pairs[i]["var_row"] for i in va_idx + te_idx})
    tr_rows = sorted({pairs[i]["wt_row"] for i in tr_idx}
                     | {pairs[i]["var_row"] for i in tr_idx})
    abs_rows = [r for r in tr_rows if r not in va_te_rows]
    # precompute per-pair truth tensors once
    c_of = {i: torch.tensor(t["c"], dtype=torch.float32, device=device)
            for i, t in enumerate(targets)}
    li_of = {i: t["lig_idx"] for i, t in enumerate(targets)}
    # free-pairwise cells: (i, ligand-index-in-target)
    cells = [(i, j) for i in tr_idx for j in range(len(targets[i]["c"]))]
    # free-ligand-pair cells: (i, j1, j2) sampled per epoch
    epochs = 3 if dry else EPOCHS
    best_val = None
    best_state = None
    gradcov = None
    for ep in range(epochs):
        rng = stable_rng("stageCIIP1A", "order", name, seed, "epoch", ep)
        order = np.asarray(cells)[rng.permutation(len(cells))]
        if name == "free_ligand_pair":
            order = np.asarray([(i, int(rng.integers(len(targets[i]["c"]))),
                                 int(rng.integers(len(targets[i]["c"]))))
                                for i, _ in order])
        tot = 0.0
        for b0 in range(0, len(order), BATCH):
            batch = order[b0:b0 + BATCH]
            opt.zero_grad()
            if name == "free_pairwise":
                hats = {}
                for i in {i for i, _ in batch}:
                    t = targets[i]
                    xp_wt = Pt[pairs[i]["wt_row"]:pairs[i]["wt_row"] + 1].expand(len(t["lig_idx"]), -1)
                    xp_v = Pt[pairs[i]["var_row"]:pairs[i]["var_row"] + 1].expand(len(t["lig_idx"]), -1)
                    hats[i] = model(xp_wt, xp_v, Lt[t["lig_idx"]])
                hi = torch.stack([hats[i][j] for i, j in batch])
                ci = torch.stack([c_of[i][j] for i, j in batch])
                loss = ((hi - ci) ** 2).mean()
            elif name == "free_ligand_pair":
                out = []
                for i, j1, j2 in batch:
                    xp = Pt[pairs[i]["var_row"]:pairs[i]["var_row"] + 1]
                    out.append(model(xp, Lt[li_of[i][j1]:li_of[i][j1] + 1],
                                     Lt[li_of[i][j2]:li_of[i][j2] + 1]))
                hi = torch.stack(out)
                ci = torch.stack([c_of[i][j1] - c_of[i][j2] for i, j1, j2 in batch])
                loss = ((hi - ci) ** 2).mean()
            else:
                # per-pair hat WITH gradient (the SAME s_theta), computed
                # once per unique pair in the batch
                hats = {}
                for i in {i for i, _ in batch}:
                    t = targets[i]
                    Lm = Lt[t["lig_idx"]]
                    Pw = Pt[pairs[i]["wt_row"]:pairs[i]["wt_row"] + 1].expand(
                        len(t["lig_idx"]), -1)
                    Pv = Pt[pairs[i]["var_row"]:pairs[i]["var_row"] + 1].expand(
                        len(t["lig_idx"]), -1)
                    hats[i] = model.centered_mutation_effect(Pw, Pv, Lm)
                hi = torch.stack([hats[i][j] for i, j in batch])
                ci = torch.stack([c_of[i][j] for i, j in batch])
                loss = ((hi - ci) ** 2).mean()
            if name not in ("free_pairwise", "free_ligand_pair") and abs_rows:
                rng_a = stable_rng("stageCIIP1A", "abs", name, seed, "epoch", ep)
                ra = np.asarray(abs_rows)[rng_a.permutation(len(abs_rows))[:BATCH]]
                la = rng_a.integers(0, Y.shape[1], size=len(ra))
                fin = np.isfinite(Y[ra, la])
                ra, la = ra[fin], la[fin]
                if len(ra):
                    if name == "ligand_invariant_shift":
                        f = model(Pt[ra], Lt[la], torch.tensor(ra, device=device))
                    else:
                        f = model(Pt[ra], Lt[la])
                    loss = loss + LAMBDA_ABS * ((f - Yt[ra, la]) ** 2).mean()
            loss.backward()
            if gradcov is None:
                gradcov = {n: bool(p.grad is not None and float(p.grad.abs().max()) > 0)
                           for n, p in model.named_parameters()}
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10.0)
            opt.step()
            tot += float(loss.detach())
        val_mse = eval_contrast_mse(model, art, targets, va_idx, pairs, Pt, Lt,
                                    device, name)
        if best_val is None or val_mse < best_val:
            best_val = val_mse
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
    if best_state is not None:
        model.load_state_dict(best_state)
    # variance decomposition over abs cells (sample, checkpoint model)
    var_dec = {}
    with torch.no_grad():
        ra = torch.from_numpy(np.asarray(abs_rows[:min(64, len(abs_rows))])).to(device)
        la = torch.randint(0, Y.shape[1], (len(ra),), device=device)
        if len(ra) and name not in ("free_pairwise", "free_ligand_pair"):
            if name == "ligand_invariant_shift":
                f = model(Pt[ra], Lt[la], ra)
            else:
                f = model(Pt[ra], Lt[la])
            var_dec = {"b_P_var": None, "b_L_var": None, "s_var": None,
                       "total_var": float(f.var().detach())}
    return model, best_val, gradcov, var_dec


def eval_contrast_mse(model, art, targets, idx, pairs, Pt, Lt, device, name):
    with torch.no_grad():
        errs = []
        for i in idx:
            t = targets[i]
            Lm = Lt[t["lig_idx"]]
            Pw = Pt[pairs[i]["wt_row"]:pairs[i]["wt_row"] + 1].expand(
                len(t["lig_idx"]), -1)
            Pv = Pt[pairs[i]["var_row"]:pairs[i]["var_row"] + 1].expand(
                len(t["lig_idx"]), -1)
            if name == "free_pairwise":
                hat = model(Pw, Pv, Lm)
            else:
                hat = model.centered_mutation_effect(Pw, Pv, Lm)
            errs.append(float(((hat.cpu().numpy() - np.asarray(t["c"])) ** 2).mean()))
    return float(np.mean(errs))


def metrics(name, model, art, targets, te_idx, pairs, Pt, Lt, device):
    per_pair = []
    for i in te_idx:
        t = targets[i]
        Lm = Lt[t["lig_idx"]]
        Pw = Pt[pairs[i]["wt_row"]:pairs[i]["wt_row"] + 1].expand(
            len(t["lig_idx"]), -1)
        Pv = Pt[pairs[i]["var_row"]:pairs[i]["var_row"] + 1].expand(
            len(t["lig_idx"]), -1)
        with torch.no_grad():
            if name == "free_pairwise":
                hat = model(Pw, Pv, Lm).cpu().numpy()
            else:
                hat = model.centered_mutation_effect(Pw, Pv, Lm).cpu().numpy()
        c = np.asarray(t["c"])
        per_pair.append({
            "pair": i, "parent": pairs[i]["parent"],
            "mutation": pairs[i]["mutation"],
            "spearman": spearman(hat, c),
            "pearson": pearson(hat, c),
            "sign_acc": dead_zone_sign_acc(hat, c),
            "mse": float(np.mean((hat - c) ** 2)),
            "scale": ols_scale(c, hat),
            "n": len(c),
        })
    return per_pair


def bootstrap_gap(rows_a, rows_b, stat, key):
    """Bootstrap 2.5% lower bound of the pair-mean STAT difference;
    clusters = parents; 2000 draws; SHA-256 keyed."""
    rng = stable_rng("stageCIIP1A", "boot", key, stat, BOOT_SEED)
    parents = sorted({r["parent"] for r in rows_a} | {r["parent"] for r in rows_b})
    vals_a = {p: [] for p in parents}
    vals_b = {p: [] for p in parents}
    for r in rows_a:
        vals_a[r["parent"]].append(r[stat] if r[stat] == r[stat] else 0.0)
    for r in rows_b:
        vals_b[r["parent"]].append(r[stat] if r[stat] == r[stat] else 0.0)
    lo_vals = []
    for _ in range(BOOT_DRAWS):
        idx = rng.integers(len(parents), size=len(parents))
        sa = [x for pi in idx for x in vals_a[parents[pi]]]
        sb = [x for pi in idx for x in vals_b[parents[pi]]]
        if sa and sb:
            lo_vals.append(float(np.mean(sa) - np.mean(sb)))
    return float(np.percentile(lo_vals, 2.5))


def run(name, art, Y, prot, lig, split, device, seed, dry=False):
    te_idx = [i for i, s in enumerate(split) if s == 2]
    pairs = art["pairs"]
    targets = art["targets"]
    Pt = torch.from_numpy(prot).float().to(device)
    Lt = torch.from_numpy(lig).float().to(device)
    if name in ("family_shuffle", "random_protein"):
        rng_rows = stable_rng("stageCIIP1A", "rows", name, seed)
        P = shuffle_rows(prot, art, rng_rows, name == "family_shuffle")
        Pt = torch.from_numpy(P).float().to(device)
    model, best_val, gradcov, var_dec = train_arm(name, art, Y, prot, lig,
                                                  split, device, seed, dry=dry)
    rows = metrics(name, model, art, targets, te_idx, pairs, Pt, Lt, device)
    return {"arm": name, "seed": seed, "best_val_mse": best_val,
            "grad_cov": gradcov, "var_dec": var_dec, "test_rows": rows,
            "agg": {
                "spearman": float(np.nanmean([r["spearman"] for r in rows])),
                "pearson": float(np.nanmean([r["pearson"] for r in rows])),
                "sign_acc": float(np.nanmean([r["sign_acc"] for r in rows])),
                "mse": float(np.mean([r["mse"] for r in rows])),
                "scale_median": float(np.nanmedian([r["scale"] for r in rows])),
            }}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["screening", "full"], default="screening")
    ap.add_argument("--seeds", type=int, nargs="+", default=[1])
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    ap.add_argument("--dry", action="store_true")
    args = ap.parse_args()
    art, Y, prot, lig, split = load_data()
    arms = SCREENING_ARMS if args.mode == "screening" else ALL_ARMS
    out = {"schema": "MetaSieve.StageCIIP1A.Result.v1",
           "preregistration_sha256": PREREG_SHA,
           "mode": args.mode, "dry": bool(args.dry),
           "data1a_sha256": hashlib.sha256(
               (HERE / "DATA1A.json").read_bytes()).hexdigest(),
           "arms": {}}
    for seed in args.seeds:
        for arm in arms:
            print(f"=== {arm} seed {seed} ===", flush=True)
            res = run(arm, art, Y, prot, lig, split, args.device, seed,
                      dry=args.dry)
            out["arms"].setdefault(arm, {})[str(seed)] = res
            print(json.dumps(res["agg"], indent=1), flush=True)
    path = HERE / ("RESULT_DRY.json" if args.dry
                   else f"RESULT_{args.mode.upper()}.json")
    path.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("wrote", path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
