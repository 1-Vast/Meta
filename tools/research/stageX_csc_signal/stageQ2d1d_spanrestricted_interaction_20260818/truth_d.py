"""Stage Q2d-1c truth generator + feature-space oracle precheck.
Prereg SHA: 25b8b9129120d0a770ba353cd56a8a388dd847778e4e5cb2b488ee0cbfee7106.
Fixes vs Q2d-1b: PCA-32 protein features, resolved-only ligand pool (157),
feature-smoothed double centring. ALS is DIAGNOSTIC ONLY.
"""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
X0C = HERE.parent / "stageX0c_measurement_qualification_20260818"
PARENT = HERE.parent
sys.path.insert(0, str(X0C))
sys.path.insert(0, str(PARENT))
from x0_common import stable_rng
import x0_i1
from x0_i1 import (_parent_of, _parent_seq, _mutation_for_row)
from x0_i2 import (klifs_pocket_for_parent, align_pocket_to_sequence,
                    mutate_pocket)

PREREG_SHA = "baf4bb72df02e9411d6b8d4815302ec91c7526cc15447b6e80cd06383d546991"
TAU = 1.0
RANK = 4
COLD_FAMILIES = ["Tec", "FGFR", "LRRK", "STE7", "Src"]
LIG_PROJ_DIM = 48
PROT_PCA_DIM = 32
PROT_DIM = 510

AA_PHYS = {
    "A": [1.8, 88.6, 0, 0, 0, 1], "R": [-4.5, 173.4, 1, 1, 3, 1],
    "N": [-3.5, 114.1, 0, 1, 1, 2], "D": [-3.5, 111.1, -1, 1, 0, 3],
    "C": [2.5, 108.5, 0, 0, 0, 1], "Q": [-3.5, 143.8, 0, 1, 1, 2],
    "E": [-3.5, 138.4, -1, 1, 0, 3], "G": [-0.4, 60.1, 0, 0, 0, 1],
    "H": [-3.2, 153.2, 1, 1, 1, 2], "I": [4.5, 166.7, 0, 0, 0, 1],
    "L": [3.8, 166.7, 0, 0, 0, 1], "K": [-3.9, 168.6, 1, 1, 2, 1],
    "M": [1.9, 162.9, 0, 0, 0, 1], "F": [2.8, 189.9, 0, 0, 0, 1],
    "P": [-1.6, 112.7, 0, 0, 0, 1], "S": [-0.8, 89.0, 0, 1, 1, 2],
    "T": [-0.7, 116.1, 0, 1, 1, 2], "W": [-0.9, 227.8, 0, 0, 1, 1],
    "Y": [-1.3, 193.6, 0, 1, 1, 2], "V": [4.2, 140.0, 0, 0, 0, 1],
}
_AA_MEAN = np.mean([AA_PHYS[a] for a in sorted(AA_PHYS)], axis=0)
_AA_STD = np.std([AA_PHYS[a] for a in sorted(AA_PHYS)], axis=0)
AA_Z = {a: ((np.asarray(v) - _AA_MEAN) / _AA_STD).tolist() for a, v in AA_PHYS.items()}


def build_protein_510(rows):
    pair_table = json.loads((PARENT / "X0_PAIR_TABLE.json").read_text(encoding="utf-8"))
    klifs = json.loads((PARENT / "klifs" / "klifs_kinase_lookup.json").read_text(encoding="utf-8"))
    from x0_common import load_duongly
    _info, _matrix, seqs = load_duongly()
    pocket_cache = {}
    P = np.zeros((len(rows), PROT_DIM), dtype=np.float32)
    for i, row in enumerate(rows):
        parent = _parent_of(row)
        if parent not in pocket_cache:
            seq = _parent_seq(parent, seqs)
            pocket, kid, note = klifs_pocket_for_parent(parent, klifs)
            align = align_pocket_to_sequence(pocket, seq) if pocket else None
            pocket_cache[parent] = (pocket, align, kid, note)
        pocket, align, kid, note = pocket_cache[parent]
        if pocket is None:
            continue
        vec = pocket
        mut = _mutation_for_row(row, pair_table)
        if mut and align is not None:
            mp, pidx, _ = mutate_pocket(pocket, mut["pos"], mut["old"], mut["new"],
                                       _parent_seq(parent, seqs), align)
            if pidx is not None:
                vec = mp
        for k, aa in enumerate(vec):
            if aa in AA_Z:
                P[i, k * 6:(k + 1) * 6] = AA_Z[aa]
    return P


def build_ligand_48(lig_feats_resolved, train_lig_mask):
    rng = stable_rng("stageQ2d1d", "truth", "wl")
    density = 0.1
    W = (rng.random((lig_feats_resolved.shape[1], LIG_PROJ_DIM)) < density).astype(np.float64)
    W = W * rng.normal(0, (1.0 / (density * lig_feats_resolved.shape[1])) ** 0.5, size=W.shape)
    L = lig_feats_resolved.astype(np.float64) @ W
    mu = L[train_lig_mask].mean(axis=0)
    sd = L[train_lig_mask].std(axis=0)
    sd[sd < 1e-9] = 1.0
    L = (L - mu) / sd
    return L.astype(np.float32), W


def make_cold_splits(rows, resolved_idx, scaffolds, family_of, rng):
    n_rows = len(rows)
    n_lig = len(resolved_idx)
    cold_row = np.zeros(n_rows, dtype=bool)
    for i, r in enumerate(rows):
        if family_of(_parent_of(r)) in COLD_FAMILIES:
            cold_row[i] = True
    clus = {}
    for jj, j in enumerate(resolved_idx):
        name = compounds_resolved_names[jj]
        clus[jj] = str(scaffolds.get(name))
    from collections import Counter
    sizes = Counter(clus.values())
    order = sorted(sizes.items(), key=lambda kv: (-kv[1], kv[0]))
    cold_scaffolds = set()
    held = 0
    for sc, n in order:
        if held >= 0.25 * n_lig:
            break
        cold_scaffolds.add(sc)
        held += n
    cold_lig = np.asarray([clus[jj] in cold_scaffolds for jj in range(n_lig)])
    train_row = ~cold_row
    train_lig = ~cold_lig
    tr_r = np.where(train_row)[0]
    tr_l = np.where(train_lig)[0]
    cd_r = np.where(cold_row)[0]
    cd_l = np.where(cold_lig)[0]

    def sample(rs, ls, cap):
        cells = [(int(i), int(j)) for i in rs for j in ls]
        if len(cells) > cap:
            idx = rng.choice(len(cells), size=cap, replace=False)
            cells = [cells[k] for k in sorted(idx)]
        return np.asarray(cells, dtype=np.int64).reshape(-1, 2)
    train_cells = np.asarray([(int(i), int(j)) for i in tr_r for j in tr_l],
                              dtype=np.int64)
    pc = sample(cd_r, tr_l, 800)
    lc = sample(tr_r, cd_l, 800)
    dc = sample(cd_r, cd_l, 800)
    return {"cold_row": cold_row, "cold_lig": cold_lig, "train_row": train_row,
            "train_lig": train_lig, "cold_families": COLD_FAMILIES,
            "cold_scaffolds": sorted(cold_scaffolds),
            "n_cold_rows": int(cold_row.sum()), "n_cold_ligs": int(cold_lig.sum()),
            "train_cells": train_cells, "pc": pc, "lc": lc, "dc": dc,
            "n_lig": n_lig}


def _span_projection(P_t, splits):
    Xtr = P_t[splits["train_row"]].astype(np.float64)
    _U, S, Vt = np.linalg.svd(Xtr, full_matrices=False)
    r = int((S > 1e-6).sum())
    return Vt[:r].T, r


def _qr_scale(X, scales):
    Q, _ = np.linalg.qr(X)
    return Q * np.asarray(scales)


def generate_truth(mechanism, seed, P_t, L_t, splits):
    rng = stable_rng("stageQ2d1d", "truth", mechanism, "seed", seed)
    n_rows, n_lig = len(P_t), len(L_t)
    d_p = P_t.shape[1]
    d_l = L_t.shape[1]
    global VSPAN
    VSPAN, _ = _span_projection(P_t, splits)
    if mechanism in ("M1", "M3"):
        # span-restricted protein map: A_t = V_train @ C, C Gaussian
        # QR-orthonormal with frozen scales (prereg Q2d-1d)
        C = rng.normal(0, 1, size=(VSPAN.shape[1], RANK))
        C = _qr_scale(C, [1.0, 0.8, 0.6, 0.4])
        A = VSPAN @ C
        B = rng.normal(0, 1, size=(d_l, RANK))
        B = _qr_scale(B, [1.0, 1.0, 1.0, 1.0])
    elif mechanism == "M2":
        # block-sparse in the PRE-COMPRESSION 510 space, then compressed
        A510 = np.zeros((PROT_DIM, RANK))
        B = rng.normal(0, 1, size=(d_l, RANK))
        for k in range(RANK):
            lo = k * 21 * 6
            hi = min((k + 1) * 21 * 6, PROT_DIM)
            A510[lo:hi, k] = rng.normal(0, 1, size=hi - lo)
        A510 = _qr_scale(A510, [1.0, 0.8, 0.6, 0.4])
        A = VSPAN @ (VSPAN.T @ (PCA_VT @ A510))
        A = _qr_scale(A, [1.0, 0.8, 0.6, 0.4])
        B = _qr_scale(B, [1.0, 1.0, 1.0, 1.0])
    elif mechanism == "NC2":
        F_r = rng.normal(0, 1, size=(n_rows, RANK))
        F_l = rng.normal(0, 1, size=(n_lig, RANK))
    if mechanism == "NC2":
        I_raw_all = F_r @ F_l.T
    else:
        I_raw_all = (P_t @ A) @ (L_t @ B).T
        if mechanism == "M3":
            I_raw_all = np.tanh(I_raw_all / np.sqrt(RANK))
    tr = splits["train_cells"]
    tr_i = tr[:, 0]
    tr_j = tr[:, 1]
    I_tr = I_raw_all[tr_i, tr_j]
    row_mean = np.zeros(n_rows)
    col_mean = np.zeros(n_lig)
    for i in range(n_rows):
        m = tr_i == i
        row_mean[i] = I_tr[m].mean() if m.any() else 0.0
    for j in range(n_lig):
        m = tr_j == j
        col_mean[j] = I_tr[m].mean() if m.any() else 0.0
    # feature-smoothed double centring: offsets are feature-linear, hence
    # representable by the learner
    train_rows_mask = splits["train_row"]
    train_ligs_mask = splits["train_lig"]
    Xr = P_t[train_rows_mask]
    Xl = L_t[train_ligs_mask]
    w_r, *_ = np.linalg.lstsq(Xr.astype(np.float64), row_mean[train_rows_mask].astype(np.float64),
                              rcond=None)
    w_c, *_ = np.linalg.lstsq(Xl.astype(np.float64), col_mean[train_ligs_mask].astype(np.float64),
                              rcond=None)
    off_r = P_t @ w_r
    off_c = L_t @ w_c
    I_c = I_raw_all - off_r[:, None] - off_c[None, :]
    sd_tr = I_c[tr_i, tr_j].std()
    I = I_c / sd_tr * TAU
    noise_rng = stable_rng("stageQ2d1d", "truth", mechanism, "noise", seed)
    noise = noise_rng.normal(0, 1.0, size=(n_rows, n_lig))
    if mechanism == "NC1":
        I = I * 0.0
        mu = 0.5
        pm = rng.normal(0, 1.0, size=n_rows)
        lm = rng.normal(0, 1.0, size=n_lig)
    else:
        mu, pm, lm = 0.0, np.zeros(n_rows), np.zeros(n_lig)
    return {"I": I, "I_raw": I_raw_all, "noise": noise, "mu": mu,
            "pm": pm, "lm": lm, "A": A, "B": B,
            "w_r": w_r, "w_c": w_c,
            "sd_train": float(sd_tr), "mechanism": mechanism, "seed": seed}


def svd_fit(I_tr, r, l, P_t, L_t, splits):
    """Closed-form rank-4 reconstruction matching the learner's capacity:
    train-estimable feature offsets (row/col means regressed on features)
    are removed first, then truncated SVD on the complete train grid, then
    feature least squares for A, B. Diagnostic only."""
    n_r = P_t.shape[0]
    n_l = L_t.shape[0]
    rows_u = sorted(set(r.tolist()))
    ligs_u = sorted(set(l.tolist()))
    trm = np.zeros(n_r)
    tcm = np.zeros(n_l)
    for i in rows_u:
        trm[i] = I_tr[np.asarray(r) == i].mean()
    for j in ligs_u:
        tcm[j] = I_tr[np.asarray(l) == j].mean()
    w_r, *_ = np.linalg.lstsq(P_t[rows_u].astype(np.float64), trm[rows_u].astype(np.float64),
                              rcond=None)
    w_c, *_ = np.linalg.lstsq(L_t[ligs_u].astype(np.float64), tcm[ligs_u].astype(np.float64),
                              rcond=None)
    off_r = P_t @ w_r
    off_c = L_t @ w_c
    resid = I_tr - off_r[r] - off_c[l]
    M = np.zeros((len(rows_u), len(ligs_u)))
    ri = {v: k for k, v in enumerate(rows_u)}
    li = {v: k for k, v in enumerate(ligs_u)}
    M[[ri[v] for v in r], [li[v] for v in l]] = resid
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    Fr = U[:, :RANK] * S[:RANK]
    Fl = Vt[:RANK, :].T
    A, *_ = np.linalg.lstsq(P_t[rows_u].astype(np.float64), Fr, rcond=None)
    B, *_ = np.linalg.lstsq(L_t[ligs_u].astype(np.float64), Fl, rcond=None)
    return A, B, w_r, w_c


def eval_oracle(A, B, w_r, w_c, P_t, L_t, cells, I_truth, mechanism):
    r, l = cells[:, 0], cells[:, 1]
    hat = ((P_t[r] @ A) * (L_t[l] @ B)).sum(-1) + P_t[r] @ w_r + L_t[l] @ w_c
    if mechanism == "M3":
        hat = np.tanh(hat / np.sqrt(RANK))
    import q2
    return q2.eval_metrics(hat, I_truth[r, l])


def main():
    global compounds_resolved_names, PCA_VT
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, _meta = x0_i1.load_features()
    resolved_idx = np.asarray([j for j in range(len(compounds))
                                if not str(scaffolds.get(compounds[j], "")).startswith("unresolved")],
                               dtype=np.int64)
    compounds_resolved_names = [compounds[j] for j in resolved_idx]
    print("n resolved ligands:", len(resolved_idx), flush=True)
    P510 = build_protein_510(rows)
    Pm = P510 - P510.mean(axis=0)
    _U, _S, Vt = np.linalg.svd(Pm.astype(np.float64), full_matrices=False)
    PCA_VT = Vt[:PROT_PCA_DIM]
    P_t = (Pm @ PCA_VT.T).astype(np.float32)
    rng_split = stable_rng("stageQ2d1d", "splits")
    splits = make_cold_splits(rows, resolved_idx, scaffolds, q2.family_of_parent, rng_split)
    L_t, W = build_ligand_48(lig_feats[resolved_idx], splits["train_lig"])
    np.savez(HERE / "q2d1d_features.npz", P_t=P_t, L_t=L_t, W=W,
             PCA_VT=PCA_VT, resolved_idx=resolved_idx)
    import hashlib
    fsha = hashlib.sha256((HERE / "q2d1d_features.npz").read_bytes()).hexdigest()
    splits_out = {}
    for k, v in splits.items():
        if isinstance(v, np.ndarray):
            splits_out[k] = v.tolist()
        else:
            splits_out[k] = v
    splits_out["feature_npz_sha256"] = fsha
    splits_out["n_train_cells"] = int(len(splits["train_cells"]))
    splits_out["obs_per_param"] = {"M1": len(splits["train_cells"]) / (32 * 4 + 48 * 4.0),
                                   "M2": len(splits["train_cells"]) / (4 * 21 * 6 + 48 * 4.0)}
    splits_out["preregistration_sha256"] = PREREG_SHA
    splits_out["schema"] = "MetaSieve.StageQ2d1d.SPLITS.v1"
    json.dump(splits_out, open(HERE / "Q2D1D_SPLITS.json", "w"), indent=1)
    pre = {"schema": "MetaSieve.StageQ2d1d.ORACLE_PRECHECK.v1",
           "preregistration_sha256": PREREG_SHA, "per_mechanism": {}}
    for mech in ("M1", "M2", "M3"):
        pre["per_mechanism"][mech] = {}
        for seed in (0, 1, 2):
            t = generate_truth(mech, seed, P_t, L_t, splits)
            r = splits["train_cells"][:, 0]
            l = splits["train_cells"][:, 1]
            A, B, w_r, w_c = svd_fit(t["I"][r, l], r, l, P_t, L_t, splits)
            res = {}
            for surf in ("pc", "lc", "dc"):
                m = eval_oracle(A, B, w_r, w_c, P_t, L_t, splits[surf], t["I"], mech)
                res[surf] = {"dz": m["dead_zone_sign_accuracy"],
                             "sp": m["spearman"]}
                print(mech, seed, surf, round(m["dead_zone_sign_accuracy"], 3),
                      round(m["spearman"], 3), flush=True)
            pre["per_mechanism"][mech][str(seed)] = res
    m1_ok = all(
        pre["per_mechanism"]["M1"][str(s)][sf]["dz"] >= 0.70
        for s in (0, 1, 2) for sf in ("pc", "lc", "dc"))
    pre["M1_identifiable_on_all_surfaces"] = bool(m1_ok)
    pre["STOP_before_training"] = not m1_ok
    json.dump(pre, open(HERE / "Q2D1D_ORACLE_PRECHECK.json", "w"), indent=1)
    print("M1 identifiable:", m1_ok)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
