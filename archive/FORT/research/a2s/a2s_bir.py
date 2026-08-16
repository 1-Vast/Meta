"""A2S-BIR: Budget-constrained Identifiable meta-Residual for A2S-DTA.

Root cause addressed (task.md, 2026-08-01): prediction decomposes into a target
level A(t), a target-INDEPENDENT ligand ordering R(d), and a target-conditioned
interaction residual E(t,d).  A(t) and R(d) are identified; E(t,d) has failed in
five framings.  The diagnosis is not that adapters lack capacity but that k<=5
recipient labels carry too little identifiable interaction information, so any
adapter degenerates to R(d) + A(t).

A2S-BIR attacks identifiability directly, with three coupled parts:

  1. BUDGET.     The residual code dimension is capped by the support budget,
                 m(k) = {1:0, 3:1, 5:2}.  At k=1 nothing beyond the anchor is
                 estimable and the model says so by construction.
  2. BASIS.      The residual basis U is meta-learned for *k-shot recoverability*
                 by backpropagating through the inner ridge solve, instead of the
                 usual variance-maximising SVD/PCA basis.  This is the core
                 innovation: the basis is optimised for what k labels can locate,
                 not for what abundant data can explain.
  3. CERTIFICATE. Per episode, a support-only posterior-contraction statistic
                 decides whether the code is identified.  When it is not, the
                 model abstains to the anchor arm exactly.

  yhat_r(q) = f0(q) + a_r + 1[certified_r] * g_m(q)^T what_r

Protein sequence never enters the mechanism; target specificity is carried only
by the k observed support labels.  Protein-feature injection is a destructive
control that must NOT help.

Endpoint pKi only.  Recipient query labels are read once, by the evaluator.
"""

from __future__ import annotations

import argparse
import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch

ROOT = Path(__file__).resolve().parents[2]
CORPUS = ROOT / "dataset" / "formal_training" / "chembl37_pki_formal.v4"
ROSTER = ROOT / "dataset" / "formal_training" / "a2s_d0r_roster.v2"
DEFAULT_OUT = ROOT / "reports" / "active"

# ---- frozen protocol constants ---------------------------------------------
SUPPORT_K = (1, 3, 5)
CODE_BUDGET = {1: 0, 3: 1, 5: 2}          # m(k): preregistered, not tuned
RIDGE_F0 = 10.0                           # f0 pooled ridge penalty
CERT_CONTRACTION = 0.50                   # posterior-contraction certificate
CERT_MIN_SPREAD = 0.10                    # support residual spread floor (pKi)
META_EPOCHS = 60
META_EPISODES_PER_TARGET = 8
META_LR = 3e-3
N_FOLDS = 5
N_BOOTSTRAP = 5000
SEED = 1729
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


# ============================================================== data loading
@dataclass
class Corpus:
    feat: np.ndarray                       # (n_parents, D) frozen ligand features
    parent_index: dict[str, int]
    source_targets: list[str]
    recipient_targets: list[str]
    rows: pd.DataFrame                     # target_uid, parent, year, doc, pKi
    recipients: pd.DataFrame
    draws: pd.DataFrame
    query: pd.DataFrame
    components: dict[str, int]
    feature_dim: int = field(init=False)

    def __post_init__(self) -> None:
        self.feature_dim = self.feat.shape[1]


def load_corpus(roster: Path | None = None) -> Corpus:
    global ROSTER
    if roster is not None:
        ROSTER = roster
    z = np.load(CORPUS / "features" / "ligand_features.npz", allow_pickle=True)
    parents = [str(p) for p in z["parent_uids"]]
    ecfp = z["ecfp4"].astype(np.float32)
    desc = z["descriptors"].astype(np.float32)
    desc = np.nan_to_num(desc, nan=0.0, posinf=0.0, neginf=0.0)
    feat = np.hstack([ecfp, desc])
    parent_index = {p: i for i, p in enumerate(parents)}

    rows = pd.read_parquet(
        CORPUS / "canonical" / "pki_measurements_exact.parquet",
        columns=["target_uid", "compound_parent_uid", "document_uid",
                 "document_year", "pKi"])
    rows = rows[rows.document_year.notna()].copy()
    rows["document_year"] = rows.document_year.astype(int)
    # one label per (target, parent): median across assay contexts / replicates
    rows = (rows.groupby(["target_uid", "compound_parent_uid"], as_index=False)
            .agg(pKi=("pKi", "median"),
                 document_year=("document_year", "min"),
                 document_uid=("document_uid", "first")))

    recipients = pd.read_parquet(ROSTER / "recipients.parquet")
    draws = pd.read_parquet(ROSTER / "support_draws.parquet")
    query = pd.read_parquet(ROSTER / "query.parquet")
    sources = pd.read_parquet(ROSTER / "sources.parquet")
    comps = dict(zip(sources.target_uid, sources.component_id))
    comps.update(dict(zip(recipients.target_uid, recipients.component_id)))

    return Corpus(feat=feat, parent_index=parent_index,
                  source_targets=sorted(sources.target_uid),
                  recipient_targets=sorted(recipients.target_uid),
                  rows=rows, recipients=recipients, draws=draws,
                  query=query, components=comps)


def feature_matrix(corpus: Corpus, parents: list[str]) -> np.ndarray:
    idx = [corpus.parent_index[p] for p in parents]
    return corpus.feat[idx]


# ============================================================== f0 (pooled prior)
class PooledPrior:
    """Target-balanced pooled ridge on frozen ligand features: the R(d) term."""

    def __init__(self, ridge: float = RIDGE_F0) -> None:
        self.ridge = ridge
        self.mean: np.ndarray | None = None
        self.scale: np.ndarray | None = None
        self.beta: torch.Tensor | None = None
        self.intercept: float = 0.0

    def fit(self, x: np.ndarray, y: np.ndarray, target_ids: np.ndarray) -> "PooledPrior":
        self.mean = x.mean(0)
        self.scale = x.std(0) + 1e-6
        xs = (x - self.mean) / self.scale
        # each source target contributes equal total weight
        _, inv, counts = np.unique(target_ids, return_inverse=True, return_counts=True)
        w = (1.0 / counts[inv]).astype(np.float32)
        w = w / w.mean()
        xt = torch.tensor(xs, dtype=torch.float32, device=DEVICE)
        yt = torch.tensor(y, dtype=torch.float32, device=DEVICE)
        wt = torch.tensor(w, dtype=torch.float32, device=DEVICE)
        self.intercept = float((wt * yt).sum() / wt.sum())
        yc = yt - self.intercept
        xw = xt * wt.unsqueeze(1)
        gram = xt.T @ xw + self.ridge * torch.eye(xt.shape[1], device=DEVICE)
        self.beta = torch.linalg.solve(gram, xw.T @ yc)
        return self

    def predict(self, x: np.ndarray) -> np.ndarray:
        xs = (x - self.mean) / self.scale
        xt = torch.tensor(xs, dtype=torch.float32, device=DEVICE)
        return (xt @ self.beta + self.intercept).cpu().numpy()

    def project(self, x: np.ndarray) -> torch.Tensor:
        xs = (x - self.mean) / self.scale
        return torch.tensor(xs, dtype=torch.float32, device=DEVICE)


# ============================================================== episodes
def document_ordered_episode(frame: pd.DataFrame, k: int, rng: np.random.Generator,
                             qmin: int = 5) -> tuple[list[str], list[str]] | None:
    """Support from earlier documents, query from strictly later ones."""
    years = sorted(frame.document_year.unique())
    if len(years) < 2:
        return None
    order = rng.permutation(len(years) - 1)
    for pos in order:
        split = years[pos]
        pre = frame[frame.document_year <= split]
        post = frame[frame.document_year > split]
        pre_par = list(dict.fromkeys(pre.compound_parent_uid))
        post = post[~post.document_uid.isin(set(pre.document_uid))]
        post_par = [p for p in dict.fromkeys(post.compound_parent_uid) if p not in set(pre_par)]
        if len(pre_par) >= k and len(post_par) >= qmin:
            sup = list(rng.choice(pre_par, size=k, replace=False))
            return sup, post_par
    return None


def random_episode(frame: pd.DataFrame, k: int, rng: np.random.Generator,
                   qmin: int = 5) -> tuple[list[str], list[str]] | None:
    """Exchangeable control protocol used by MAML/ANIL/MetaDTA/AdaMBind."""
    par = list(dict.fromkeys(frame.compound_parent_uid))
    if len(par) < k + qmin:
        return None
    perm = rng.permutation(len(par))
    sup = [par[i] for i in perm[:k]]
    qry = [par[i] for i in perm[k:]]
    return sup, qry


# ============================================================== hierarchy
@dataclass
class Hyper:
    """Hierarchical scales estimated on source targets only.

    tau_b : across-target SD of the target anchor b_t = mean(y - f0)
    tau_z : across-target SD of each task-code coordinate
    sigma : within-target residual noise SD after anchor and code
    drift : mean(later-document residual) - mean(earlier-document residual)
    """
    tau_b: float = 1.0
    tau_z: float = 1.0
    sigma: float = 0.7
    drift: float = 0.0


def hierarchical_posterior(g_sup: torch.Tensor, r_sup: torch.Tensor, hyper: Hyper
                           ) -> tuple[torch.Tensor, torch.Tensor]:
    """Joint exact Gaussian posterior over [b, z]. Differentiable.

    Design is [1, g(d)]; the prior shrinks the anchor toward the global mean with
    scale tau_b and each code coordinate toward zero with scale tau_z.  Shrinking
    the anchor is what makes k=1 well-posed: with one label the anchor is a
    James-Stein estimate, not a raw residual.
    """
    n, m = g_sup.shape[0], g_sup.shape[1]
    design = torch.cat([torch.ones(n, 1, device=g_sup.device, dtype=g_sup.dtype), g_sup], dim=1)
    prior_prec = torch.diag(torch.tensor(
        [1.0 / max(hyper.tau_b ** 2, 1e-8)] + [1.0 / max(hyper.tau_z ** 2, 1e-8)] * m,
        device=g_sup.device, dtype=g_sup.dtype))
    precision = prior_prec + design.T @ design / (hyper.sigma ** 2)
    cov = torch.linalg.inv(precision)
    mean = cov @ (design.T @ r_sup) / (hyper.sigma ** 2)
    return mean, cov


def certificate(g_sup: torch.Tensor, r_sup: torch.Tensor, hyper: Hyper) -> tuple[bool, float]:
    """Support-only interaction-identifiability certificate.

    Posterior contraction of the worst-constrained CODE direction against its
    prior, computed from the joint posterior so that anchor/code collinearity is
    charged to the code.  Reads no query feature and no query label.
    """
    m = g_sup.shape[1]
    if m == 0:
        return False, 0.0
    _, cov = hierarchical_posterior(g_sup, r_sup, hyper)
    code_cov = cov[1:, 1:]
    post_max = float(torch.linalg.eigvalsh(code_cov).max())
    contraction = 1.0 - post_max / max(hyper.tau_z ** 2, 1e-8)
    spread = float(r_sup.std(unbiased=False)) if r_sup.numel() > 1 else 0.0
    ok = contraction >= CERT_CONTRACTION and spread >= CERT_MIN_SPREAD
    return ok, contraction


def estimate_hyper(corpus: Corpus, prior: PooledPrior, target_ids: list[str],
                   u: torch.Tensor | None = None) -> Hyper:
    """Method-of-moments hierarchical scales from source targets only."""
    anchors, within, drifts, codes = [], [], [], []
    for t, f in corpus.rows[corpus.rows.target_uid.isin(target_ids)].groupby("target_uid"):
        pars = list(dict.fromkeys(f.compound_parent_uid))
        if len(pars) < 20:
            continue
        x = feature_matrix(corpus, pars)
        sub = f.groupby("compound_parent_uid").agg(pKi=("pKi", "first"),
                                                   year=("document_year", "min")).reindex(pars)
        r = sub.pKi.to_numpy() - prior.predict(x)
        anchors.append(float(np.mean(r)))
        within.append(float(np.var(r)))
        years = sub.year.to_numpy()
        cut = np.median(years)
        early, late = r[years <= cut], r[years > cut]
        if len(early) >= 3 and len(late) >= 3:
            drifts.append(float(late.mean() - early.mean()))
        if u is not None and u.shape[1] > 0:
            g = prior.project(x) @ u
            rc = torch.tensor(r - r.mean(), dtype=torch.float32, device=DEVICE)
            sol = torch.linalg.lstsq(g, rc.unsqueeze(1)).solution.squeeze(1)
            codes.append(sol.cpu().numpy())
    tau_b = float(np.std(anchors)) if anchors else 1.0
    sigma = float(np.sqrt(np.mean(within))) if within else 0.7
    drift = float(np.mean(drifts)) if drifts else 0.0
    tau_z = float(np.mean(np.std(np.array(codes), axis=0))) if codes else 1.0
    return Hyper(tau_b=max(tau_b, 1e-3), tau_z=max(tau_z, 1e-3),
                 sigma=max(sigma, 1e-2), drift=drift)


# ============================================================== meta-learned basis
def _target_cache(corpus: Corpus, prior: PooledPrior, target_ids: list[str],
                  min_parents: int = 20) -> dict[str, tuple[torch.Tensor, torch.Tensor, dict[str, int], pd.DataFrame]]:
    cache: dict[str, tuple[torch.Tensor, torch.Tensor, dict[str, int], pd.DataFrame]] = {}
    for t, f in corpus.rows[corpus.rows.target_uid.isin(target_ids)].groupby("target_uid"):
        pars = list(dict.fromkeys(f.compound_parent_uid))
        if len(pars) < min_parents:
            continue
        x = feature_matrix(corpus, pars)
        resid = torch.tensor(
            f.groupby("compound_parent_uid").pKi.first().reindex(pars).to_numpy()
            - prior.predict(x), dtype=torch.float32, device=DEVICE)
        cache[t] = (prior.project(x), resid, {p: i for i, p in enumerate(pars)}, f)
    return cache


def meta_learn_basis(corpus: Corpus, prior: PooledPrior, target_ids: list[str],
                     m: int, rng: np.random.Generator, protocol: str, hyper: Hyper,
                     epochs: int = META_EPOCHS, log: list[str] | None = None) -> torch.Tensor:
    """Learn U in R^{D x m} for k-shot recoverability.

    Objective: expected held-out query error of the *deployed* predictor -- the
    joint hierarchical posterior over [b, z] solved on k support points, plus the
    drift term.  Gradients flow through the closed-form solve, so U is optimised
    for what k labels can actually locate rather than for variance explained.
    """
    if m == 0:
        return torch.zeros((corpus.feature_dim, 0), device=DEVICE)

    cache = _target_cache(corpus, prior, target_ids)
    usable = sorted(cache)
    u = torch.nn.Parameter(torch.randn(corpus.feature_dim, m, device=DEVICE) * 0.01)
    opt = torch.optim.Adam([u], lr=META_LR)
    episode_fn = document_ordered_episode if protocol == "ordered" else random_episode
    k_train = max(SUPPORT_K)

    micro = 64
    for epoch in range(epochs):
        running, n_ep, chunk, n_chunk = 0.0, 0, None, 0
        opt.zero_grad()
        for t in usable:
            g_full, resid, pos, frame = cache[t]
            for _ in range(META_EPISODES_PER_TARGET):
                ep = episode_fn(frame, k_train, rng)
                if ep is None:
                    continue
                si = [pos[p] for p in ep[0] if p in pos]
                qi = [pos[p] for p in ep[1] if p in pos]
                if len(si) < k_train or len(qi) < 3:
                    continue
                gs, gq = g_full[si] @ u, g_full[qi] @ u
                post, _ = hierarchical_posterior(gs, resid[si], hyper)
                pred = post[0] + hyper.drift + gq @ post[1:]
                loss = ((resid[qi] - pred) ** 2).mean()
                chunk = loss if chunk is None else chunk + loss
                n_chunk += 1
                n_ep += 1
                if n_chunk >= micro:
                    (chunk / n_chunk).backward()
                    running += float(chunk)
                    chunk, n_chunk = None, 0
        if n_chunk:
            (chunk / n_chunk).backward()
            running += float(chunk)
        if n_ep == 0:
            break
        for prm in [u]:
            if prm.grad is not None:
                prm.grad /= max(n_ep / micro, 1.0)
        opt.step()
        total = running
        with torch.no_grad():                              # orthonormalise for nesting
            q, _ = torch.linalg.qr(u)
            u.copy_(q[:, :m])
        if log is not None and (epoch % 15 == 0 or epoch == epochs - 1):
            log.append(f"    basis m={m} protocol={protocol} epoch {epoch:>3} "
                       f"held-out query MSE {float(total / n_ep):.4f} ({n_ep} episodes)")
    return u.detach()


def svd_basis(corpus: Corpus, prior: PooledPrior, target_ids: list[str], m: int) -> torch.Tensor:
    """Variance-maximising baseline basis (the conventional choice)."""
    if m == 0:
        return torch.zeros((corpus.feature_dim, 0), device=DEVICE)
    mats = []
    for t, f in corpus.rows[corpus.rows.target_uid.isin(target_ids)].groupby("target_uid"):
        pars = list(dict.fromkeys(f.compound_parent_uid))
        if len(pars) < 20:
            continue
        x = feature_matrix(corpus, pars)
        r = f.groupby("compound_parent_uid").pKi.first().reindex(pars).to_numpy() - prior.predict(x)
        r = r - r.mean()
        g = prior.project(x)
        mats.append((g * torch.tensor(r, dtype=torch.float32, device=DEVICE).unsqueeze(1)).mean(0))
    stack = torch.stack(mats)                       # (n_targets, D) residual directions
    _, _, vh = torch.linalg.svd(stack, full_matrices=False)
    return vh[:m].T.contiguous()


# ============================================================== local residual
class LocalResidual(torch.nn.Module):
    """Meta-learned support-local interaction residual with per-query certification.

    Motivated by the measured failure of the global rank-m code: with across-target
    code scale tau_z ~ 0.18 pKi against within-target noise sigma ~ 1.0 pKi, one
    global code dimension needs ~ (sigma/tau_z)^2 ~ 31 labels.  A support-local
    correction instead spends at most k degrees of freedom -- it is a convex
    combination of the k OBSERVED support residuals, so it is identifiable by
    construction.  What must be meta-learned is (a) the metric deciding which
    support compound is informative for a query, and (b) how much local support a
    query actually has.

        yhat(q) = f0(q) + b_hat + rho(q) * sum_i alpha_i(q) * (r_i - b_hat)

    rho is the per-query coverage certificate: a query far from every support
    compound falls back to the anchor exactly.  Nothing here reads the protein.
    """

    def __init__(self, dim: int, rank: int = 16, raw: bool = False,
                 metric: str = "diagonal") -> None:
        super().__init__()
        self.raw = raw
        self.metric = metric
        # diagonal metric initialised at the raw cosine (weights == 1), so the
        # meta-learned metric can only improve on the untrained baseline
        self.log_w = torch.nn.Parameter(torch.zeros(dim))
        self.proj = torch.nn.Parameter(torch.randn(dim, rank) * (1.0 / np.sqrt(dim)))
        self.log_lambda = torch.nn.Parameter(torch.tensor(0.0))   # the identifiability budget
        self.gate_slope = torch.nn.Parameter(torch.tensor(4.0))
        self.gate_bias = torch.nn.Parameter(torch.tensor(0.3))

    def embed(self, g: torch.Tensor) -> torch.Tensor:
        if self.raw:
            z = g
        elif self.metric == "linear":
            return g          # unnormalised linear kernel (no cosine normalisation)
        elif self.metric == "diagonal":
            z = g * torch.nn.functional.softplus(self.log_w + 0.5413)   # softplus(0.5413)~1
        else:
            z = g @ self.proj
        return torch.nn.functional.normalize(z, dim=1)

    def similarity(self, g_q: torch.Tensor, g_s: torch.Tensor) -> torch.Tensor:
        return self.embed(g_q) @ self.embed(g_s).T                # (n_query, k) in [-1,1]

    def forward(self, g_q: torch.Tensor, g_s: torch.Tensor, r_s: torch.Tensor,
                anchor: torch.Tensor, use_gate: bool = True) -> torch.Tensor:
        """Kernel-ridge correction on the k observed support residuals.

        The meta-learned lambda sets the effective degrees of freedom
        tr(K_ss (K_ss + lambda I)^-1) <= k: this IS the identifiability budget,
        learned rather than fixed.
        """
        zs = self.embed(g_s)
        k_ss = zs @ zs.T
        lam = torch.nn.functional.softplus(self.log_lambda) + 1e-3
        eye = torch.eye(k_ss.shape[0], device=k_ss.device, dtype=k_ss.dtype)
        alpha = torch.linalg.solve(k_ss + lam * eye, (r_s - anchor).unsqueeze(1)).squeeze(1)
        local = (self.embed(g_q) @ zs.T) @ alpha
        if use_gate:
            local = self.coverage(g_q, g_s) * local
        return anchor + local

    def coverage(self, g_q: torch.Tensor, g_s: torch.Tensor) -> torch.Tensor:
        sim = self.similarity(g_q, g_s)
        return torch.sigmoid(self.gate_slope * (sim.max(dim=1).values - self.gate_bias))

    def effective_dof(self, g_s: torch.Tensor) -> float:
        zs = self.embed(g_s)
        k_ss = zs @ zs.T
        lam = torch.nn.functional.softplus(self.log_lambda) + 1e-3
        eye = torch.eye(k_ss.shape[0], device=k_ss.device, dtype=k_ss.dtype)
        return float(torch.trace(k_ss @ torch.linalg.inv(k_ss + lam * eye)))


def train_local_residual(corpus: Corpus, prior: PooledPrior, target_ids: list[str],
                         rng: np.random.Generator, protocol: str, hyper: Hyper,
                         epochs: int = META_EPOCHS, rank: int = 16,
                         log: list[str] | None = None, raw: bool = False,
                         metric: str = "diagonal") -> LocalResidual:
    model = LocalResidual(corpus.feature_dim, rank, raw=raw, metric=metric).to(DEVICE)
    opt = torch.optim.Adam(model.parameters(), lr=1e-2)
    cache = _target_cache(corpus, prior, target_ids)
    usable = sorted(cache)
    episode_fn = document_ordered_episode if protocol == "ordered" else random_episode
    k_train = max(SUPPORT_K)
    micro = 32                       # backward every `micro` episodes to bound memory
    for epoch in range(epochs):
        opt.zero_grad()
        running, n_ep, chunk, n_chunk = 0.0, 0, None, 0
        for t in usable:
            g_full, resid, pos, frame = cache[t]
            for _ in range(META_EPISODES_PER_TARGET):
                ep = episode_fn(frame, k_train, rng)
                if ep is None:
                    continue
                si = [pos[p] for p in ep[0] if p in pos]
                qi = [pos[p] for p in ep[1] if p in pos]
                if len(si) < k_train or len(qi) < 3:
                    continue
                r_s = resid[si]
                post, _ = hierarchical_posterior(
                    torch.zeros((len(si), 0), device=DEVICE), r_s, hyper)
                pred = model(g_full[qi], g_full[si], r_s, post[0])
                loss = ((resid[qi] - pred) ** 2).mean()
                chunk = loss if chunk is None else chunk + loss
                n_chunk += 1
                n_ep += 1
                if n_chunk >= micro:
                    (chunk / n_chunk).backward()
                    running += float(chunk)
                    chunk, n_chunk = None, 0
        if n_chunk:
            (chunk / n_chunk).backward()
            running += float(chunk)
        if n_ep == 0:
            break
        for prm in model.parameters():
            if prm.grad is not None:
                prm.grad /= max(n_ep / micro, 1.0)
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        opt.step()
        total = running
        if log is not None and (epoch % 15 == 0 or epoch == epochs - 1):
            lam = float(torch.nn.functional.softplus(model.log_lambda) + 1e-3)
            log.append(f"    local residual (raw={raw}) epoch {epoch:>3} held-out query MSE "
                       f"{float(total / n_ep):.4f} lambda={lam:.4f} "
                       f"gate_bias={float(model.gate_bias):.3f}")
    return model.eval()


# ============================================================== gradient baselines
class MetaMLP(torch.nn.Module):
    """Small head over frozen ligand features, used for the MAML/ANIL baselines.

    Predicts the RESIDUAL y - f0(d), which is centred and O(1); regressing raw
    pKi from a 2058-d input made the inner loop diverge (RMSE ~ 3.8e5).  The
    LayerNorm bounds body activations and the head is zero-initialised so the
    unadapted model starts exactly at the pooled prior.
    """

    def __init__(self, dim: int, hidden: int = 128) -> None:
        super().__init__()
        self.body = torch.nn.Sequential(
            torch.nn.Linear(dim, hidden), torch.nn.ReLU(), torch.nn.LayerNorm(hidden))
        self.head = torch.nn.Linear(hidden, 1)
        torch.nn.init.zeros_(self.head.weight)
        torch.nn.init.zeros_(self.head.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.head(self.body(x)).squeeze(-1)


def _inner_adapt(model: MetaMLP, params: dict, g: torch.Tensor, r: torch.Tensor,
                 adapt: list, steps: int, lr: float, create_graph: bool) -> dict:
    for _ in range(steps):
        loss = ((torch.func.functional_call(model, params, (g,)) - r) ** 2).mean()
        grads = torch.autograd.grad(loss, [params[n] for n in adapt],
                                    create_graph=create_graph)
        clipped = []
        for gr in grads:                      # inner clipping: MAML on k<=5 points
            norm = gr.norm()                  # is otherwise unstable
            clipped.append(gr if norm <= 10.0 else gr * (10.0 / (norm + 1e-12)))
        params = {**params, **{n: params[n] - lr * gr for n, gr in zip(adapt, clipped)}}
    return params


def train_gradient_meta(corpus, prior, target_ids, rng, protocol, mode,
                        inner_steps: int = 5, inner_lr=None, epochs: int = 30,
                        log=None):
    """MAML (adapt all) or ANIL (adapt head only).

    inner_lr is selected on SOURCE episodes only; no recipient label is read.
    """
    cache = _target_cache(corpus, prior, target_ids)
    usable = sorted(cache)
    episode_fn = document_ordered_episode if protocol == "ordered" else random_episode

    def episodes(n, gen):
        out = []
        for t in gen.permutation(usable)[:n]:
            g_full, resid, pos, frame = cache[t]
            ep = episode_fn(frame, max(SUPPORT_K), gen)
            if ep is None:
                continue
            si = [pos[x] for x in ep[0] if x in pos]
            qi = [pos[x] for x in ep[1] if x in pos]
            if len(si) < max(SUPPORT_K) or len(qi) < 3:
                continue
            out.append((g_full[si], resid[si], g_full[qi], resid[qi]))
        return out

    if inner_lr is None:
        best, best_loss = 1e-2, float("inf")
        probe = episodes(40, np.random.default_rng(SEED + 99))
        scratch = MetaMLP(corpus.feature_dim).to(DEVICE)
        base = {n: v.detach() for n, v in scratch.named_parameters()}
        adapt = [n for n in base if mode == "maml" or n.startswith("head")]
        for lr in (1e-3, 1e-2, 5e-2):
            tot = 0.0
            for gs, rs, gq, rq in probe:
                fast = _inner_adapt(
                    scratch, {n: v.clone().requires_grad_(True) for n, v in base.items()},
                    gs, rs, adapt, inner_steps, lr, False)
                with torch.no_grad():
                    tot += float(((torch.func.functional_call(scratch, fast, (gq,))
                                   - rq) ** 2).mean())
            if probe and tot / len(probe) < best_loss:
                best, best_loss = lr, tot / len(probe)
        inner_lr = best
        if log is not None:
            log.append("  %s inner_lr=%g selected on source episodes (held-out MSE %.4f)"
                       % (mode.upper(), inner_lr, best_loss))

    model = MetaMLP(corpus.feature_dim).to(DEVICE)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    for epoch in range(epochs):
        opt.zero_grad()
        batch = episodes(48, rng)
        losses = []
        for gs, rs, gq, rq in batch:
            params = {n: v for n, v in model.named_parameters()}
            adapt = [n for n in params if mode == "maml" or n.startswith("head")]
            fast = _inner_adapt(model, params, gs, rs, adapt, inner_steps, inner_lr, True)
            losses.append(((torch.func.functional_call(model, fast, (gq,)) - rq) ** 2).mean())
        if not losses:
            break
        loss = torch.stack(losses).mean()
        if not torch.isfinite(loss):
            break
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        opt.step()
        if log is not None and (epoch % 10 == 0 or epoch == epochs - 1):
            log.append("    %s epoch %3d outer query MSE %.4f" % (mode.upper(), epoch, float(loss)))
    return model, inner_lr


def adapt_gradient_meta(model, g_sup, r_sup, g_qry, mode, inner_lr, inner_steps: int = 5):
    """Returns the predicted RESIDUAL; the caller adds f0 back."""
    params = {n: v.detach().clone().requires_grad_(True) for n, v in model.named_parameters()}
    adapt = [n for n in params if mode == "maml" or n.startswith("head")]
    fast = _inner_adapt(model, params, g_sup, r_sup, adapt, inner_steps, inner_lr, False)
    with torch.no_grad():
        out = torch.func.functional_call(model, fast, (g_qry,)).cpu().numpy()
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)


# ============================================================== metrics
def ndcg_at_k(y_true: np.ndarray, y_pred: np.ndarray, k: int = 10) -> float:
    if len(y_true) < k:
        return float("nan")
    gain = y_true - y_true.min()
    if gain.sum() <= 0:
        return float("nan")
    disc = 1.0 / np.log2(np.arange(2, k + 2))
    dcg = (gain[np.argsort(-y_pred)][:k] * disc).sum()
    idcg = (np.sort(gain)[::-1][:k] * disc).sum()
    return float(dcg / idcg) if idcg > 0 else float("nan")


def pairwise_accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Constant predictors score at chance (0.5), not 0."""
    dt = y_true[:, None] - y_true[None, :]
    dp = y_pred[:, None] - y_pred[None, :]
    mask = np.triu(np.abs(dt) > 1e-9, k=1)
    if mask.sum() == 0:
        return float("nan")
    correct = ((np.sign(dt) == np.sign(dp)) & mask).sum()
    ties = ((np.abs(dp) <= 1e-12) & mask).sum()
    return float((correct + 0.5 * ties) / mask.sum())


def spearman(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """A constant predictor carries no ordering information: score 0, not NaN."""
    if len(y_true) < 3:
        return float("nan")
    if np.std(y_pred) < 1e-12:
        return 0.0
    from scipy.stats import spearmanr
    r = spearmanr(y_true, y_pred).statistic
    return float(r) if np.isfinite(r) else 0.0


def concordance_index(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """CI, the standard DTA ranking metric (Gonen & Heller; used by DeepDTA et al.).

    CI = (1/Z) sum_{y_i > y_j} h(f_i - f_j), h(x) = 1 if x>0, 0.5 if x==0, else 0.
    """
    dt = y_true[:, None] - y_true[None, :]
    dp = y_pred[:, None] - y_pred[None, :]
    mask = dt > 1e-9                                   # ordered pairs only
    z = mask.sum()
    if z == 0:
        return float("nan")
    conc = ((dp > 1e-12) & mask).sum() + 0.5 * ((np.abs(dp) <= 1e-12) & mask).sum()
    return float(conc / z)


def pearson(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    if len(y_true) < 3 or np.std(y_pred) < 1e-12 or np.std(y_true) < 1e-12:
        return 0.0
    r = float(np.corrcoef(y_true, y_pred)[0, 1])
    return r if np.isfinite(r) else 0.0


def r_squared(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    ss_tot = float(((y_true - y_true.mean()) ** 2).sum())
    if ss_tot <= 0:
        return float("nan")
    return float(1.0 - ((y_true - y_pred) ** 2).sum() / ss_tot)


def rm2(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """Roy's modified squared correlation, standard in the DTA literature.

    rm2 = r2 * (1 - sqrt(|r2 - r02|)), where r02 is the squared correlation of the
    regression through the origin.
    """
    r = pearson(y_true, y_pred)
    r2 = r * r
    denom = float((y_pred ** 2).sum())
    if denom <= 0:
        return float("nan")
    k = float((y_true * y_pred).sum() / denom)         # through-origin slope
    ss_tot = float(((y_true - y_true.mean()) ** 2).sum())
    if ss_tot <= 0:
        return float("nan")
    r02 = 1.0 - float(((y_true - k * y_pred) ** 2).sum()) / ss_tot
    return float(r2 * (1.0 - np.sqrt(abs(r2 - r02))))


def aupr(y_true: np.ndarray, y_pred: np.ndarray, threshold: float = 7.0) -> float:
    """Area under precision-recall for actives at pKi >= threshold."""
    labels = (y_true >= threshold).astype(int)
    if labels.sum() == 0 or labels.sum() == len(labels):
        return float("nan")
    from sklearn.metrics import average_precision_score
    return float(average_precision_score(labels, y_pred))


def episode_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    err = y_pred - y_true
    return {"rmse": float(np.sqrt(np.mean(err ** 2))),
            "mse": float(np.mean(err ** 2)),
            "mae": float(np.mean(np.abs(err))),
            "ci": concordance_index(y_true, y_pred),
            "rm2": rm2(y_true, y_pred),
            "pearson": pearson(y_true, y_pred),
            "r2": r_squared(y_true, y_pred),
            "spearman": spearman(y_true, y_pred),
            "pairwise": pairwise_accuracy(y_true, y_pred),
            "ndcg10": ndcg_at_k(y_true, y_pred, 10),
            "aupr": aupr(y_true, y_pred)}


METRIC_NAMES = ["rmse", "mse", "mae", "ci", "rm2", "pearson", "r2",
                "spearman", "pairwise", "ndcg10", "aupr"]


def component_bootstrap(per_target: dict[str, float], components: dict[str, int],
                        rng: np.random.Generator, n: int = N_BOOTSTRAP) -> dict[str, float]:
    """Paired bootstrap over independent homology components."""
    by_comp: dict[int, list[float]] = {}
    for t, v in per_target.items():
        if np.isfinite(v):
            by_comp.setdefault(components[t], []).append(v)
    comps = sorted(by_comp)
    if not comps:
        return {"mean": float("nan"), "lo": float("nan"), "hi": float("nan"), "n_components": 0}
    means = np.array([np.mean(by_comp[c]) for c in comps])
    idx = rng.integers(0, len(comps), size=(n, len(comps)))
    boots = means[idx].mean(1)
    return {"mean": float(means.mean()),
            "lo": float(np.quantile(boots, 0.025)),
            "hi": float(np.quantile(boots, 0.975)),
            "n_components": len(comps)}


# ============================================================== evaluation
def evaluate(corpus: Corpus, out: Path, protocol: str = "ordered",
             quick: bool = False) -> dict[str, Any]:
    t_start = time.time()
    rng = np.random.default_rng(SEED)
    torch.manual_seed(SEED)
    log: list[str] = []

    src_rows = corpus.rows[corpus.rows.target_uid.isin(corpus.source_targets)]
    x_src = feature_matrix(corpus, list(src_rows.compound_parent_uid))
    prior = PooledPrior().fit(x_src, src_rows.pKi.to_numpy(), src_rows.target_uid.to_numpy())
    log.append(f"  f0 pooled prior fitted on {len(src_rows)} source rows / "
               f"{src_rows.target_uid.nunique()} targets")

    hyper0 = estimate_hyper(corpus, prior, corpus.source_targets)
    log.append(f"  hierarchy (source-only): tau_b={hyper0.tau_b:.4f} "
               f"sigma={hyper0.sigma:.4f} drift={hyper0.drift:+.4f}")

    epochs = 10 if quick else META_EPOCHS
    bases: dict[int, dict[str, torch.Tensor]] = {}
    hypers: dict[int, dict[str, Hyper]] = {}
    for m in sorted({v for v in CODE_BUDGET.values() if v > 0}):
        bases[m] = {
            "meta": meta_learn_basis(corpus, prior, corpus.source_targets, m, rng,
                                     protocol, hyper0, epochs, log),
            "svd": svd_basis(corpus, prior, corpus.source_targets, m),
            "random": torch.linalg.qr(torch.randn(corpus.feature_dim, m, device=DEVICE))[0],
        }
        hypers[m] = {key: estimate_hyper(corpus, prior, corpus.source_targets, u)
                     for key, u in bases[m].items()}
        log.append(f"  m={m} tau_z: " + ", ".join(
            f"{key}={h.tau_z:.4f}" for key, h in hypers[m].items()))

    # identifiability budget: labels needed to resolve one global code dimension
    tau_z_meta = hypers[max(hypers)]["meta"].tau_z if hypers else float("nan")
    labels_needed = float((hyper0.sigma / max(tau_z_meta, 1e-9)) ** 2)
    log.append(f"  global-code identifiability: sigma/tau_z = "
               f"{hyper0.sigma / max(tau_z_meta, 1e-9):.2f}, "
               f"labels needed for 1 dimension ~ {labels_needed:.1f} (available k<=5)")

    local = train_local_residual(corpus, prior, corpus.source_targets, rng,
                                 protocol, hyper0, epochs, 16, log)
    local_raw = train_local_residual(corpus, prior, corpus.source_targets, rng,
                                     protocol, hyper0, epochs, 16, log, raw=True)
    local_rand = LocalResidual(corpus.feature_dim, 16, metric="projection").to(DEVICE).eval()
    local_lin = train_local_residual(corpus, prior, corpus.source_targets, rng,
                                     protocol, hyper0, epochs, 16, log, metric="linear")

    grad_models = {}
    if not quick:
        for mode in ("maml", "anil"):
            grad_models[mode] = train_gradient_meta(corpus, prior, corpus.source_targets,
                                                    rng, protocol, mode, log=log)

    query_by_target = {t: list(f.compound_parent_uid) for t, f in corpus.query.groupby("target_uid")}
    draws_by_target = {t: f for t, f in corpus.draws.groupby("target_uid")}
    labels = corpus.rows.set_index(["target_uid", "compound_parent_uid"]).pKi

    arms = ["recipient_calibration", "f0_only", "f0_anchor", "f0_anchor_shrunk",
            "a2s_bir_global", "a2s_bir_global_nocert", "a2s_bir_nodrift",
            "a2s_bir_svd", "a2s_bir_random_basis",
            "a2s_bir_local", "a2s_bir_local_nogate", "a2s_bir_local_randmetric",
            "a2s_bir_local_rawmetric", "a2s_bir_local_unnorm",
            "knn_dta", "pooled_finetune",
            "wrong_support", "perm_support", "local_wrong_support", "local_perm_support"]
    if not quick:
        arms += ["maml", "anil"]
    FALLBACK = "f0_anchor_shrunk"      # the strongest no-residual arm

    results: dict[int, dict[str, dict[str, dict[str, float]]]] = {}
    coverage: dict[int, dict[str, float]] = {}

    for k in SUPPORT_K:
        m = CODE_BUDGET[k]
        per_target: dict[str, dict[str, list[dict[str, float]]]] = {}
        cert_flags: dict[str, list[float]] = {}
        local_cov: dict[str, list[float]] = {}
        local_dof: dict[str, list[float]] = {}

        for target in corpus.recipient_targets:
            qpars = query_by_target.get(target, [])
            if len(qpars) < 3:
                continue
            x_q = feature_matrix(corpus, qpars)
            y_q = np.array([labels.loc[(target, p)] for p in qpars], dtype=np.float32)
            g_q_full = prior.project(x_q)
            f0_q = prior.predict(x_q)
            dframe = draws_by_target[target]
            dframe_k = dframe[dframe.k == k]
            arm_eps: dict[str, list[dict[str, float]]] = {a: [] for a in arms}
            certs: list[float] = []
            cov_q: list[float] = []
            dof_q: list[float] = []

            for draw_id, sub in dframe_k.groupby("draw_id"):
                spars = list(sub.compound_parent_uid)
                x_s = feature_matrix(corpus, spars)
                y_s = np.array([labels.loc[(target, p)] for p in spars], dtype=np.float32)
                g_s_full = prior.project(x_s)
                f0_s = prior.predict(x_s)
                anchor = float(np.mean(y_s - f0_s))
                r_s = torch.tensor(y_s - f0_s, dtype=torch.float32, device=DEVICE)

                def bir(u: torch.Tensor, resid: torch.Tensor, g_sup: torch.Tensor,
                        hyp: Hyper, certify: bool = True, drift: bool = True
                        ) -> tuple[np.ndarray, bool]:
                    """Deployed predictor: hierarchical [b, z] + certificate + drift."""
                    g_s = g_sup @ u
                    ok = True
                    if certify:
                        ok, _ = certificate(g_s, resid, hyp)
                    post, _ = hierarchical_posterior(g_s if ok else g_s[:, :0], resid, hyp)
                    out = f0_q + float(post[0]) + (hyp.drift if drift else 0.0)
                    if ok and u.shape[1] > 0:
                        out = out + ((g_q_full @ u) @ post[1:]).cpu().numpy()
                    return out, ok

                empty = torch.zeros((corpus.feature_dim, 0), device=DEVICE)
                h0 = hypers[m]["meta"] if m > 0 else hyper0
                pred = {
                    "recipient_calibration": np.full(len(qpars), float(np.mean(y_s))),
                    "f0_only": f0_q,
                    "f0_anchor": f0_q + anchor,
                    "f0_anchor_shrunk": bir(empty, r_s, g_s_full, h0, certify=False)[0],
                }

                # ---- A2S-BIR global-code family (measured negative) -------
                for arm, key in (("a2s_bir_global", "meta"), ("a2s_bir_svd", "svd"),
                                 ("a2s_bir_random_basis", "random")):
                    if m == 0:
                        pred[arm] = pred["f0_anchor_shrunk"]
                        continue
                    pred[arm], ok = bir(bases[m][key], r_s, g_s_full, hypers[m][key])
                    if key == "meta":
                        certs.append(float(ok))

                if m == 0:
                    pred["a2s_bir_global_nocert"] = pred["f0_anchor_shrunk"]
                    pred["a2s_bir_nodrift"] = bir(empty, r_s, g_s_full, h0,
                                                  certify=False, drift=False)[0]
                else:
                    pred["a2s_bir_global_nocert"] = bir(bases[m]["meta"], r_s, g_s_full,
                                                        hypers[m]["meta"], certify=False)[0]
                    pred["a2s_bir_nodrift"] = bir(bases[m]["meta"], r_s, g_s_full,
                                                  hypers[m]["meta"], drift=False)[0]

                # ---- A2S-BIR support-local family -------------------------
                post_b, _ = hierarchical_posterior(
                    torch.zeros((len(spars), 0), device=DEVICE), r_s, h0)
                anchor_t = post_b[0]
                with torch.no_grad():
                    pred["a2s_bir_local"] = (
                        f0_q + local(g_q_full, g_s_full, r_s, anchor_t).cpu().numpy())
                    pred["a2s_bir_local_nogate"] = (
                        f0_q + local(g_q_full, g_s_full, r_s, anchor_t,
                                     use_gate=False).cpu().numpy())
                    pred["a2s_bir_local_randmetric"] = (
                        f0_q + local_rand(g_q_full, g_s_full, r_s, anchor_t).cpu().numpy())
                    pred["a2s_bir_local_rawmetric"] = (
                        f0_q + local_raw(g_q_full, g_s_full, r_s, anchor_t).cpu().numpy())
                    pred["a2s_bir_local_unnorm"] = (
                        f0_q + local_lin(g_q_full, g_s_full, r_s, anchor_t).cpu().numpy())
                    cov_q.append(float(local.coverage(g_q_full, g_s_full).mean()))
                    dof_q.append(local.effective_dof(g_s_full))

                # ---- destructive: permuted support labels ---------------
                yp = y_s[rng.permutation(len(y_s))]
                rp = torch.tensor(yp - f0_s, dtype=torch.float32, device=DEVICE)
                pred["perm_support"] = bir(bases[m]["meta"] if m else empty, rp, g_s_full,
                                           hypers[m]["meta"] if m else h0)[0]
                pb, _ = hierarchical_posterior(
                    torch.zeros((len(spars), 0), device=DEVICE), rp, h0)
                with torch.no_grad():
                    pred["local_perm_support"] = (
                        f0_q + local(g_q_full, g_s_full, rp, pb[0]).cpu().numpy())

                # ---- destructive: another recipient's support -----------
                other = corpus.recipient_targets[
                    (corpus.recipient_targets.index(target) + 7) % len(corpus.recipient_targets)]
                osub = draws_by_target[other]
                osub = osub[(osub.k == k) & (osub.draw_id == draw_id)]
                o_pars = list(osub.compound_parent_uid)
                x_o = feature_matrix(corpus, o_pars)
                y_o = np.array([labels.loc[(other, p)] for p in o_pars], dtype=np.float32)
                g_o_full = prior.project(x_o)
                r_o = torch.tensor(y_o - prior.predict(x_o), dtype=torch.float32, device=DEVICE)
                pred["wrong_support"] = bir(bases[m]["meta"] if m else empty, r_o,
                                            g_o_full, hypers[m]["meta"] if m else h0)[0]
                po, _ = hierarchical_posterior(
                    torch.zeros((len(o_pars), 0), device=DEVICE), r_o, h0)
                with torch.no_grad():
                    pred["local_wrong_support"] = (
                        f0_q + local(g_q_full, g_o_full, r_o, po[0]).cpu().numpy())

                # ---- kNN-DTA over the recipient support ------------------
                sim = (g_q_full @ g_s_full.T /
                       (g_q_full.norm(dim=1, keepdim=True) * g_s_full.norm(dim=1) + 1e-8))
                wgt = torch.softmax(sim * 4.0, dim=1).cpu().numpy()
                pred["knn_dta"] = wgt @ y_s

                # ---- equal-budget pooled fine-tune -----------------------
                lam = 1.0
                gs = g_s_full
                delta = torch.linalg.solve(
                    gs.T @ gs + lam * torch.eye(gs.shape[1], device=DEVICE),
                    gs.T @ torch.tensor(y_s - f0_s, dtype=torch.float32, device=DEVICE))
                pred["pooled_finetune"] = f0_q + (g_q_full @ delta).cpu().numpy()

                # ---- gradient meta-learners -----------------------------
                if not quick:
                    for mode in ("maml", "anil"):
                        mdl, ilr = grad_models[mode]
                        pred[mode] = f0_q + adapt_gradient_meta(
                            mdl, g_s_full, r_s, g_q_full, mode, ilr)

                for arm in arms:
                    arm_eps[arm].append(episode_metrics(y_q, np.asarray(pred[arm], dtype=float)))

            per_target[target] = {a: arm_eps[a] for a in arms}
            cert_flags[target] = certs
            local_cov[target] = cov_q
            local_dof[target] = dof_q

        # ---- aggregate: average draws within recipient, then macro -------
        metric_names = ["rmse", "mae", "spearman", "pairwise", "ndcg10"]
        summary: dict[str, dict[str, float]] = {}
        target_means: dict[str, dict[str, dict[str, float]]] = {}
        for arm in arms:
            target_means[arm] = {}
            for target, arms_eps in per_target.items():
                eps = arms_eps[arm]
                target_means[arm][target] = {
                    mn: float(np.nanmean([e[mn] for e in eps])) for mn in metric_names}
            summary[arm] = {mn: float(np.nanmean(
                [v[mn] for v in target_means[arm].values()])) for mn in metric_names}

        # ---- paired gains vs the no-transfer control ---------------------
        gains: dict[str, dict[str, Any]] = {}
        base = target_means["recipient_calibration"]
        for arm in arms:
            if arm == "recipient_calibration":
                continue
            d_rmse = {t: base[t]["rmse"] - target_means[arm][t]["rmse"] for t in base}
            entry = {"rmse_gain": component_bootstrap(d_rmse, corpus.components, rng)}
            for mn in ("spearman", "pairwise", "ndcg10"):
                d = {t: target_means[arm][t][mn] - base[t][mn] for t in base}
                entry[f"{mn}_gain"] = component_bootstrap(d, corpus.components, rng)
            vals = np.array([v for v in d_rmse.values() if np.isfinite(v)])
            entry["negative_transfer_rate"] = float((vals < 0).mean())
            entry["benefiting_fraction"] = float((vals > 0).mean())
            entry["median_gain"] = float(np.median(vals))
            gains[arm] = entry

        # ---- gain of every arm over the strongest no-residual fallback ----
        fb = target_means[FALLBACK]
        for arm in arms:
            if arm in (FALLBACK, "recipient_calibration"):
                continue
            d = {t: fb[t]["rmse"] - target_means[arm][t]["rmse"] for t in base}
            gains[arm]["rmse_gain_vs_fallback"] = component_bootstrap(d, corpus.components, rng)
            ds = {t: target_means[arm][t]["spearman"] - fb[t]["spearman"] for t in base}
            gains[arm]["spearman_gain_vs_fallback"] = component_bootstrap(
                ds, corpus.components, rng)
            vals = np.array([v for v in d.values() if np.isfinite(v)])
            gains[arm]["negative_transfer_rate_vs_fallback"] = float((vals < 0).mean())

        # ---- gain against the strongest external baseline -----------------
        strongest = min((a for a in arms if a not in
                         ("recipient_calibration", "f0_only", "f0_anchor")
                         and not a.startswith("a2s_bir")
                         and not a.endswith("support")),
                        key=lambda a: summary[a]["rmse"])
        sb = target_means[strongest]
        for arm in arms:
            if arm == strongest:
                continue
            d = {t: sb[t]["rmse"] - target_means[arm][t]["rmse"] for t in base}
            gains.setdefault(arm, {})["rmse_gain_vs_strongest_baseline"] =                 component_bootstrap(d, corpus.components, rng)
            ds = {t: target_means[arm][t]["spearman"] - sb[t]["spearman"] for t in base}
            gains[arm]["spearman_gain_vs_strongest_baseline"] =                 component_bootstrap(ds, corpus.components, rng)
        gains["_strongest_baseline"] = strongest

        results[k] = {"summary": summary, "gains": gains,
                      "fallback_arm": FALLBACK}
        cov = [np.mean(v) for v in cert_flags.values() if v]
        lcov = [np.mean(v) for v in local_cov.values() if v]
        ldof = [np.mean(v) for v in local_dof.values() if v]
        coverage[k] = {"global_code_certified_fraction": float(np.mean(cov)) if cov else 0.0,
                       "local_mean_coverage": float(np.mean(lcov)) if lcov else 0.0,
                       "local_effective_dof": float(np.mean(ldof)) if ldof else 0.0,
                       "code_dimension": m,
                       "recipients_evaluated": len(per_target)}
        print(f"  k={k} (m={m}) evaluated on {len(per_target)} recipients, "
              f"global-code certified {coverage[k]['global_code_certified_fraction']:.3f}, "
              f"local coverage {coverage[k]['local_mean_coverage']:.3f}")

    report = {
        "schema": "a2s-bir-evaluation-v1",
        "model": "A2S-BIR: budget-constrained identifiable meta-residual",
        "endpoint": "pKi",
        "protocol": protocol,
        "seed": SEED,
        "device": DEVICE,
        "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "roster": {"recipients": len(corpus.recipient_targets),
                   "sources": len(corpus.source_targets),
                   "components": len(set(corpus.components[t] for t in corpus.recipient_targets))},
        "code_budget": CODE_BUDGET,
        "certificate": {"contraction_threshold": CERT_CONTRACTION,
                        "spread_threshold": CERT_MIN_SPREAD},
        "coverage": coverage,
        "identifiability": {
            "tau_b_across_target_anchor_sd": round(hyper0.tau_b, 4),
            "sigma_within_target_noise_sd": round(hyper0.sigma, 4),
            "tau_z_across_target_code_sd": round(float(tau_z_meta), 4),
            "document_ordered_drift": round(hyper0.drift, 4),
            "labels_needed_for_one_global_code_dim": round(labels_needed, 1),
            "labels_available": max(SUPPORT_K),
        },
        "results": {str(k): v for k, v in results.items()},
        "statistical_unit": "independent homology component (paired bootstrap, 5000)",
        "wall_time_s": round(time.time() - t_start, 1),
        "peak_torch_mem_mb": (round(torch.cuda.max_memory_allocated() / 2**20, 1)
                              if torch.cuda.is_available() else None),
        "log": log,
    }
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"a2s_bir_{protocol}_seed{SEED}.json"
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"\nwrote {path}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="A2S-BIR training and evaluation")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--protocol", choices=("ordered", "random"), default="ordered")
    parser.add_argument("--quick", action="store_true", help="smoke run")
    parser.add_argument("--roster", type=Path, default=None)
    args = parser.parse_args()
    corpus = load_corpus(args.roster)
    print(f"corpus: {len(corpus.source_targets)} sources, "
          f"{len(corpus.recipient_targets)} recipients, device={DEVICE}")
    report = evaluate(corpus, args.out, args.protocol, args.quick)
    for k, block in report["results"].items():
        print(f"\n=== k={k} (code dim {report['coverage'][int(k)]['code_dimension']}) ===")
        print(f"{'arm':<24}{'RMSE':>8}{'Spearman':>10}{'pairwise':>10}{'NDCG@10':>9}")
        for arm, s in block["summary"].items():
            print(f"{arm:<24}{s['rmse']:>8.4f}{s['spearman']:>10.4f}"
                  f"{s['pairwise']:>10.4f}{s['ndcg10']:>9.4f}")


if __name__ == "__main__":
    main()
