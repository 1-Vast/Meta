"""Two-stage additive + low-rank interaction fit on a partially observed panel.

Stage 1  mu + alpha_i + beta_j            (least squares on observed cells)
Stage 2  <u_i, v_j> on the stage-1 residual, main effects held FIXED.

Holding the main effects fixed is deliberate: every XP1-B arm must share
identical mu and alpha, so that arm differences are attributable only to the
interaction coordinate and not to a rank-dependent main-effect estimate.
"""
from __future__ import annotations

import numpy as np

from panels import additive_fit


def fit_interaction_basis(Y, mask, rank: int, lam: float = 1.0,
                          iters: int = 80, seed: int = 0, tol: float = 1e-8):
    Y = np.asarray(Y, float)
    M = np.asarray(mask, bool)
    n, p = Y.shape
    mu, a, b, add = additive_fit(Y, M)
    R = np.where(M, Y - add, 0.0)
    if rank == 0:
        return dict(mu=mu, alpha=a, beta=b, U=np.zeros((n, 0)), V=np.zeros((p, 0)),
                    fit=add)
    rng = np.random.default_rng(seed)
    U = rng.normal(0, 0.01, size=(n, rank))
    V = rng.normal(0, 0.01, size=(p, rank))
    Ir = lam * np.eye(rank)
    rows = [np.where(M[i])[0] for i in range(n)]
    cols = [np.where(M[:, j])[0] for j in range(p)]
    prev = np.inf
    for _ in range(iters):
        for i in range(n):
            j = rows[i]
            if j.size == 0:
                U[i] = 0.0
                continue
            Vj = V[j]
            U[i] = np.linalg.solve(Vj.T @ Vj + Ir, Vj.T @ R[i, j])
        for jx in range(p):
            i = cols[jx]
            if i.size == 0:
                V[jx] = 0.0
                continue
            Ui = U[i]
            V[jx] = np.linalg.solve(Ui.T @ Ui + Ir, Ui.T @ R[i, jx])
        sse = float(((R - U @ V.T)[M] ** 2).sum())
        if abs(prev - sse) < tol * max(sse, 1.0):
            break
        prev = sse
    return dict(mu=mu, alpha=a, beta=b, U=U, V=V, fit=add + U @ V.T)


# backwards-compatible alias used by XP1-A
fit_additive_lowrank = fit_interaction_basis


def cv_rank_curve(Y, mask, ranks, lam=1.0, folds=5, seed=0):
    """Random-cell CV: how much held-out interaction variance is reproducible."""
    Y = np.asarray(Y, float)
    M = np.asarray(mask, bool)
    idx = np.argwhere(M)
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(idx))
    fold_id = np.zeros(len(idx), int)
    fold_id[perm] = np.arange(len(idx)) % folds
    out = {r: {"sse": 0.0, "n": 0} for r in ranks}
    base = {"sse": 0.0, "n": 0}
    per_fold = {r: [] for r in ranks}
    per_fold_base = []
    for f in range(folds):
        te = idx[fold_id == f]
        Mtr = M.copy()
        Mtr[te[:, 0], te[:, 1]] = False
        yte = Y[te[:, 0], te[:, 1]]
        mu, a, b, add = additive_fit(Y, Mtr)
        pa = add[te[:, 0], te[:, 1]]
        base["sse"] += float(((yte - pa) ** 2).sum())
        base["n"] += len(te)
        per_fold_base.append(float(((yte - pa) ** 2).mean()))
        for r in ranks:
            m = fit_interaction_basis(Y, Mtr, r, lam=lam, seed=seed + 100 * f)
            pr = m["fit"][te[:, 0], te[:, 1]]
            out[r]["sse"] += float(((yte - pr) ** 2).sum())
            out[r]["n"] += len(te)
            per_fold[r].append(float(((yte - pr) ** 2).mean()))
    res = {}
    for r in ranks:
        res[r] = {"rmse": (out[r]["sse"] / out[r]["n"]) ** 0.5,
                  "r2_gamma": 1.0 - out[r]["sse"] / base["sse"],
                  "per_fold_mse": per_fold[r]}
    res["additive"] = {"rmse": (base["sse"] / base["n"]) ** 0.5,
                       "per_fold_mse": per_fold_base}
    return res
