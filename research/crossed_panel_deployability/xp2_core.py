"""Deterministic estimators shared by the XP2 stages."""
from __future__ import annotations

import numpy as np


def additive_fit(Y, mask, iters=300, tol=1e-11):
    Y = np.asarray(Y, float)
    M = np.asarray(mask, bool)
    n, p = Y.shape
    mu = Y[M].mean()
    a = np.zeros(n)
    b = np.zeros(p)
    rc, cc = M.sum(1), M.sum(0)
    for _ in range(iters):
        pa, pb, pm = a.copy(), b.copy(), mu
        R = np.where(M, Y - mu - b[None, :], 0.0)
        a = np.divide(R.sum(1), np.maximum(rc, 1))
        a[rc == 0] = 0.0
        R = np.where(M, Y - mu - a[:, None], 0.0)
        b = np.divide(R.sum(0), np.maximum(cc, 1))
        b[cc == 0] = 0.0
        a -= a.mean()
        b -= b.mean()
        mu = (Y - a[:, None] - b[None, :])[M].mean()
        if (abs(a - pa).max() < tol and abs(b - pb).max() < tol and abs(mu - pm) < tol):
            break
    return mu, a, b


def interaction_basis(Y, mask, rank, lam=1.0, iters=60, seed=0, tol=1e-8):
    """Two-stage: additive main effects frozen, then ALS on the residual."""
    Y = np.asarray(Y, float)
    M = np.asarray(mask, bool)
    n, p = Y.shape
    mu, a, b = additive_fit(Y, M)
    R = np.where(M, Y - mu - a[:, None] - b[None, :], 0.0)
    if rank == 0:
        return dict(mu=mu, alpha=a, beta=b, U=np.zeros((n, 0)), V=np.zeros((p, 0)))
    rng = np.random.default_rng(seed)
    U = rng.normal(0, 0.01, (n, rank))
    V = rng.normal(0, 0.01, (p, rank))
    I = lam * np.eye(rank)
    rows = [np.where(M[i])[0] for i in range(n)]
    cols = [np.where(M[:, j])[0] for j in range(p)]
    prev = np.inf
    for _ in range(iters):
        for i in range(n):
            j = rows[i]
            if j.size:
                Vj = V[j]
                U[i] = np.linalg.solve(Vj.T @ Vj + I, Vj.T @ R[i, j])
            else:
                U[i] = 0.0
        for k in range(p):
            i = cols[k]
            if i.size:
                Ui = U[i]
                V[k] = np.linalg.solve(Ui.T @ Ui + I, Ui.T @ R[i, k])
            else:
                V[k] = 0.0
        sse = float(((R - U @ V.T)[M] ** 2).sum())
        if abs(prev - sse) < tol * max(sse, 1.0):
            break
        prev = sse
    return dict(mu=mu, alpha=a, beta=b, U=U, V=V)


def ridge_predict(Xtr, Ttr, Xte, lam):
    """Ridge with unpenalised intercept; dual form when d > n."""
    xm, tm = Xtr.mean(0), Ttr.mean(0)
    A, B, C = Xtr - xm, Ttr - tm, Xte - xm
    n, d = A.shape
    if d > n:
        K = A @ A.T
        return tm + C @ (A.T @ np.linalg.solve(K + lam * np.eye(n), B))
    return tm + C @ np.linalg.solve(A.T @ A + lam * np.eye(d), A.T @ B)


def ridge_cv(Xtr, Ttr, Xte, groups, lams=(0.03, 0.1, 1.0, 10.0, 100.0, 1e3, 1e4),
             n_inner=10, seed=0):
    """Group-safe selection of lambda inside the training block.

    DOCUMENTED DEVIATION from PREREG_XP2 section 6: the registration said
    leave-one-closure-component-out.  With 206 training scaffold components that
    is ~1,400 dual ridge solves per arm per fold and is computationally
    infeasible here, so components are binned into `n_inner` grouped folds.
    Whole components still never straddle an inner split, so the procedure
    remains group-safe.  This affects hyperparameter selection inside the
    training block only; it touches no test cell, estimand or threshold.
    """
    gs = sorted(set(groups))
    rng = np.random.default_rng(seed)
    order = list(gs)
    rng.shuffle(order)
    bin_of = {g: i % n_inner for i, g in enumerate(order)}
    binid = np.array([bin_of[g] for g in groups])
    best, best_e = lams[0], np.inf
    for lam in lams:
        e = 0.0
        for b in range(n_inner):
            m = binid != b
            if m.sum() < 5 or (~m).sum() == 0:
                continue
            e += float(((ridge_predict(Xtr[m], Ttr[m], Xtr[~m], lam) - Ttr[~m]) ** 2).sum())
        if e < best_e:
            best_e, best = e, lam
    return ridge_predict(Xtr, Ttr, Xte, best), best


def solve_ligand_params(y_row, mask_row, mu, beta, V, lam):
    """Least-squares (alpha_i, u_i) for one ligand against the frozen basis."""
    j = np.where(mask_row)[0]
    r = y_row[j] - mu - beta[j]
    D = np.column_stack([np.ones(len(j)), V[j]])
    pen = np.diag([0.0] + [lam] * V.shape[1])
    sol = np.linalg.solve(D.T @ D + pen, D.T @ r)
    return float(sol[0]), sol[1:]


def cluster_bootstrap(per_unit, n_boot=2000, seed=0):
    """per_unit: dict unit -> (sse, n). Returns point + CI of MSE."""
    u = sorted(per_unit)
    s = np.array([per_unit[k][0] for k in u])
    n = np.array([per_unit[k][1] for k in u])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(u), (n_boot, len(u)))
    bs = s[idx].sum(1) / np.maximum(n[idx].sum(1), 1)
    return {"mse": float(s.sum() / n.sum()),
            "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
            "n": int(n.sum()), "units": len(u)}


def paired_contrast(a_units, b_units, n_boot=2000, seed=0):
    """MSE(b) - MSE(a): positive means arm a is better."""
    u = sorted(set(a_units) & set(b_units))
    sa = np.array([a_units[k][0] for k in u])
    sb = np.array([b_units[k][0] for k in u])
    nn = np.array([a_units[k][1] for k in u])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(u), (n_boot, len(u)))
    bs = (sb[idx].sum(1) - sa[idx].sum(1)) / np.maximum(nn[idx].sum(1), 1)
    return {"point": float((sb.sum() - sa.sum()) / nn.sum()),
            "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
            "rmse_a": float((sa.sum() / nn.sum()) ** 0.5),
            "rmse_b": float((sb.sum() / nn.sum()) ** 0.5), "units": len(u)}


def r2_vs_base(arm_units, base_units, n_boot=2000, seed=0):
    u = sorted(set(arm_units) & set(base_units))
    sa = np.array([arm_units[k][0] for k in u])
    sb = np.array([base_units[k][0] for k in u])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(u), (n_boot, len(u)))
    bs = 1.0 - sa[idx].sum(1) / np.maximum(sb[idx].sum(1), 1e-12)
    return {"point": float(1.0 - sa.sum() / sb.sum()),
            "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]}


def widest(*intervals):
    """Take the most conservative of several bootstrap intervals."""
    lo = min(i["ci95"][0] for i in intervals)
    hi = max(i["ci95"][1] for i in intervals)
    out = dict(intervals[0])
    out["ci95"] = [lo, hi]
    out["ci_source"] = "widest of %d clusterings" % len(intervals)
    return out
