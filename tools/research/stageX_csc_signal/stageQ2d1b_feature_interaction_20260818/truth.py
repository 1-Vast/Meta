"""Stage Q2d-1b truth generator + feature-space oracle precheck.
Prereg SHA: 872bc4402f228d940776e7efe2fee6b91e8310badb4e8830f653ca5e5d2e998e.
Closed-form ALS is a DIAGNOSTIC ONLY; deployment stays gradient-trained.
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

PREREG_SHA = "872bc4402f228d940776e7efe2fee6b91e8310badb4e8830f653ca5e5d2e998e"
TAU = 1.0
RANK = 4
COLD_FAMILIES = ["Tec", "FGFR", "LRRK", "STE7", "Src"]
LIG_PROJ_DIM = 64
PROT_DIM = 510

# frozen amino-acid physicochemical scales: hydrophobicity (Kyte-Doolittle),
# residue volume, charge at pH7, polarity, H-bond donor count, acceptor count
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


def build_protein_features(rows):
    """Per-row 510-dim KLIFS per-position physicochemical descriptor."""
    pair_table = json.loads((PARENT / "X0_PAIR_TABLE.json").read_text(encoding="utf-8"))
    klifs = json.loads((PARENT / "klifs" / "klifs_kinase_lookup.json").read_text(encoding="utf-8"))
    from x0_common import load_duongly
    _info, _matrix, seqs = load_duongly()
    pocket_cache = {}
    P = np.zeros((len(rows), PROT_DIM), dtype=np.float32)
    row_meta = []
    for i, row in enumerate(rows):
        parent = _parent_of(row)
        if parent not in pocket_cache:
            seq = _parent_seq(parent, seqs)
            pocket, kid, note = klifs_pocket_for_parent(parent, klifs)
            align = align_pocket_to_sequence(pocket, seq) if pocket else None
            pocket_cache[parent] = (pocket, align, kid, note)
        pocket, align, kid, note = pocket_cache[parent]
        if pocket is None:
            row_meta.append({"row": row, "status": "no_pocket", "note": note})
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
        row_meta.append({"row": row, "status": "ok"})
    return P, row_meta


def build_ligand_features(lig_feats, train_lig_mask):
    """ECFP4 (2048) through a FROZEN sparse projection to 64-dim, then
    standardized with train-ligand-only per-dim mean/std."""
    rng = stable_rng("stageQ2d1b", "truth", "wl")
    density = 0.1
    W = (rng.random((lig_feats.shape[1], LIG_PROJ_DIM)) < density).astype(np.float64)
    W = W * rng.normal(0, (1.0 / (density * lig_feats.shape[1])) ** 0.5, size=W.shape)
    L = lig_feats.astype(np.float64) @ W
    mu = L[train_lig_mask].mean(axis=0)
    sd = L[train_lig_mask].std(axis=0)
    sd[sd < 1e-9] = 1.0
    L = (L - mu) / sd
    return L.astype(np.float32), W, mu, sd


def make_cold_splits(rows, compounds, scaffolds, family_of, rng):
    n_rows, n_lig = len(rows), len(compounds)
    cold_row = np.zeros(n_rows, dtype=bool)
    for i, r in enumerate(rows):
        if family_of(_parent_of(r)) in COLD_FAMILIES:
            cold_row[i] = True
    # scaffold clusters
    clus = {}
    for j, name in enumerate(compounds):
        sc = scaffolds.get(name)
        if sc is None or str(sc).startswith("unresolved"):
            clus[j] = "unresolved"
        else:
            clus[j] = str(sc)
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
    cold_lig = np.asarray([clus[j] in cold_scaffolds for j in range(n_lig)])
    train_row = ~cold_row
    train_lig = ~cold_lig
    # surfaces
    tr_r = np.where(train_row)[0]
    tr_l = np.where(train_lig)[0]
    cd_r = np.where(cold_row)[0]
    cd_l = np.where(cold_lig)[0]

    def sample(rs, ls, cap, key):
        cells = [(i, j) for i in rs for j in ls]
        if len(cells) > cap:
            idx = rng.choice(len(cells), size=cap, replace=False)
            cells = [cells[k] for k in sorted(idx)]
        arr = np.asarray(cells, dtype=np.int64).reshape(-1, 2)
        return arr
    train_cells = np.asarray([(i, j) for i in tr_r for j in tr_l], dtype=np.int64)
    pc = sample(cd_r, tr_l, 800, "pc")
    lc = sample(tr_r, cd_l, 800, "lc")
    dc = sample(cd_r, cd_l, 800, "dc")
    return {"cold_row": cold_row, "cold_lig": cold_lig, "train_row": train_row,
            "train_lig": train_lig, "cold_families": COLD_FAMILIES,
            "cold_scaffolds": sorted(cold_scaffolds),
            "n_cold_rows": int(cold_row.sum()), "n_cold_ligs": int(cold_lig.sum()),
            "train_cells": train_cells, "pc": pc, "lc": lc, "dc": dc}


def _qr_scale(X, scales):
    Q, _ = np.linalg.qr(X)
    return Q * np.asarray(scales)


def generate_truth(mechanism, seed, P_t, L_t, splits, train_mask_for_truth=True):
    """Returns dict with I (double-centred, sd=tau over train cells),
    per-cell arrays, and mechanism metadata."""
    rng = stable_rng("stageQ2d1b", "truth", mechanism, "seed", seed)
    n_rows, n_lig = len(P_t), len(L_t)
    tr = splits["train_cells"]
    if mechanism in ("M1", "M3"):
        A = rng.normal(0, 1, size=(PROT_DIM, RANK))
        B = rng.normal(0, 1, size=(LIG_PROJ_DIM, RANK))
        A = _qr_scale(A, [1.0, 0.8, 0.6, 0.4])
        B = _qr_scale(B, [1.0, 1.0, 1.0, 1.0])
    elif mechanism == "M2":
        A = np.zeros((PROT_DIM, RANK))
        B = rng.normal(0, 1, size=(LIG_PROJ_DIM, RANK))
        for k in range(RANK):
            lo = k * 21 * 6
            hi = min((k + 1) * 21 * 6, PROT_DIM)
            A[lo:hi, k] = rng.normal(0, 1, size=hi - lo)
        A = _qr_scale(A, [1.0, 0.8, 0.6, 0.4])
        B = _qr_scale(B, [1.0, 1.0, 1.0, 1.0])
    elif mechanism == "NC2":
        F_r = rng.normal(0, 1, size=(n_rows, RANK))
        F_l = rng.normal(0, 1, size=(n_lig, RANK))
    A = None if mechanism == "NC2" else A
    B = None if mechanism == "NC2" else B
    if mechanism == "NC2":
        I_raw_all = F_r @ F_l.T
    else:
        I_raw_all = (P_t @ A) @ (L_t @ B).T  # (n_rows, n_lig)
        if mechanism == "M3":
            I_raw_all = np.tanh(I_raw_all / np.sqrt(RANK))
    # double-centre + scale with TRAIN cells only
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
    I_c = I_raw_all - row_mean[:, None] - col_mean[None, :]
    sd_tr = I_c[tr_i, tr_j].std()
    I = I_c / sd_tr * TAU
    # noise
    noise_rng = stable_rng("stageQ2d1b", "truth", mechanism, "noise", seed)
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
            "sd_train": float(sd_tr), "mechanism": mechanism, "seed": seed}


def _cells_arrays(cells):
    return cells[:, 0], cells[:, 1]


def als_fit(I_tr, r, l, P_t, L_t, rank=4, iters=10):
    """Rank-4 alternating least squares on TRAIN cells (diagnostic only)."""
    n_r = P_t.shape[0]
    n_l = L_t.shape[0]
    A = np.random.default_rng(0).normal(0, 0.1, size=(PROT_DIM, rank))
    B = np.random.default_rng(1).normal(0, 0.1, size=(LIG_PROJ_DIM, rank))
    y = I_tr.astype(np.float64)
    Xr = P_t[r].astype(np.float64)  # (n, 510)
    Xl = L_t[l].astype(np.float64)  # (n, 64)
    for _it in range(iters):
        # solve A given B: vec(y) = sum_k kron(Xl B_k, Xr) a_k
        D = Xl @ B  # (n, rank)
        des = np.hstack([Xr * D[:, k:k + 1] for k in range(rank)])  # (n, 510*rank)
        a, *_ = np.linalg.lstsq(des, y, rcond=None)
        A = a.reshape(PROT_DIM, rank)
        # solve B given A
        C = Xr @ A  # (n, rank)
        des2 = np.hstack([Xl * C[:, k:k + 1] for k in range(rank)])
        b, *_ = np.linalg.lstsq(des2, y, rcond=None)
        B = b.reshape(LIG_PROJ_DIM, rank)
    return A, B


def eval_oracle(A, B, P_t, L_t, cells, I_truth, mechanism):
    r, l = _cells_arrays(cells)
    bil = (P_t[r] @ A) * (L_t[l] @ B)
    hat = bil.sum(axis=-1)
    if mechanism == "M3":
        hat = np.tanh(hat / np.sqrt(RANK))
    import q2
    return q2.eval_metrics(hat, I_truth[r, l])


def main():
    import q2
    rows, compounds, prot_feats, lig_feats, scaffolds, _meta = x0_i1.load_features()
    P_t, row_meta = build_protein_features(rows)
    rng_split = stable_rng("stageQ2d1b", "splits")
    splits = make_cold_splits(rows, compounds, scaffolds, q2.family_of_parent, rng_split)
    train_lig_mask = splits["train_lig"]
    L_t, W, mu, sd = build_ligand_features(lig_feats, train_lig_mask)
    np.savez(HERE / "q2d1b_features.npz", P_t=P_t, L_t=L_t, W=W,
             row_meta=np.asarray(row_meta, dtype=object))
    import hashlib
    fsha = hashlib.sha256((HERE / "q2d1b_features.npz").read_bytes()).hexdigest()
    splits_out = dict(splits)
    for k, v in list(splits_out.items()):
        if isinstance(v, np.ndarray):
            splits_out[k] = v.tolist()
    splits_out["feature_npz_sha256"] = fsha
    splits_out["n_train_cells"] = len(splits["train_cells"])
    splits_out["obs_per_param"] = {"M1": len(splits["train_cells"]) / 2296.0,
                                     "M2": len(splits["train_cells"]) / 760.0}
    splits_out["preregistration_sha256"] = PREREG_SHA
    splits_out["schema"] = "MetaSieve.StageQ2d1b.SPLITS.v1"
    json.dump(splits_out, open(HERE / "Q2D1B_SPLITS.json", "w"), indent=1)
    # oracle precheck (before any training)
    pre = {"schema": "MetaSieve.StageQ2d1b.ORACLE_PRECHECK.v1",
           "preregistration_sha256": PREREG_SHA, "per_mechanism": {}}
    for mech in ("M1", "M2", "M3"):
        pre["per_mechanism"][mech] = {}
        for seed in (0, 1, 2):
            t = generate_truth(mech, seed, P_t, L_t, splits)
            r, l = splits["train_cells"][:, 0], splits["train_cells"][:, 1]
            A, B = als_fit(t["I"][r, l], r, l, P_t, L_t)
            res = {}
            for surf in ("pc", "lc", "dc"):
                m = eval_oracle(A, B, P_t, L_t, splits[surf], t["I"], mech)
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
    json.dump(pre, open(HERE / "Q2D1B_ORACLE_PRECHECK.json", "w"), indent=1)
    print("M1 identifiable:", m1_ok)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
