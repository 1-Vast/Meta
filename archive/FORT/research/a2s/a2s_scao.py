"""A2S-SCAO: Support-Conditioned Adaptation Operator for abundant-to-scarce DTA.

The prediction mechanism itself is a trainable adaptation operator meta-learned
from source episodes -- not an episodic wrapper around a closed-form solve.

    T_phi : (support set S_r) -> (adaptation state)  ->  predictor for target r

Given a recipient's k<=5 support pairs, a set encoder produces a task context,
a hypernetwork emits FiLM modulation and an anchor shift for the prediction
head, and a cross-attention path lets each query compound read the observed
support residuals directly.  Everything in phi is trained end to end on
document-ordered episodes drawn from abundant source targets.

Core innovation retained from the measured identifiability analysis
(task.md 2026-08-01): with across-target interaction scale tau_z ~ 0.18 pKi
against within-target noise sigma ~ 1.0 pKi, roughly 30 labels are needed to
resolve one global interaction direction.  So the operator does not emit a free
task vector.  It emits a BUDGETED adaptation whose magnitude is gated by a
learned identifiability head reading support-only statistics, and it falls back
exactly to the shrunk-anchor prior when the support cannot identify a direction.
This is the "interaction-identifiable, budget-constrained meta-residual".

Engineering contract: the ligand feature matrix, all labels and all episode
index tensors live on the GPU for the whole run; episodes are processed in
padded batches so a training step is one fused set of kernels.  No per-episode
host transfer, no pandas lookup inside the loop.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

from research.a2s.a2s_bir import (METRIC_NAMES, Corpus, Hyper, PooledPrior, adapt_gradient_meta,
                              component_bootstrap, episode_metrics, estimate_hyper,
                              hierarchical_posterior, load_corpus, train_gradient_meta)

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
SEED = 1729
SUPPORT_K = (1, 3, 5)
N_BOOTSTRAP = 5000

# ---- frozen training protocol ------------------------------------------------
EMBED_DIM = 256
HIDDEN_DIM = 512
N_HEADS = 4
DROPOUT = 0.1
BATCH_EPISODES = 96          # episodes per optimiser step (GPU saturation)
MAX_QUERY = 32               # padded query slots per episode
STEPS = 3000
LR = 3e-4
WEIGHT_DECAY = 1e-4
RANK_LOSS_WEIGHT = 0.3       # pairwise ranking term; CI/Spearman are endpoints too
SPEC_WEIGHT = 1.0            # counterfactual support-contrast weight
SPEC_MARGIN = 0.05           # required margin of correct support over mismatched
WARMUP = 200


# ============================================================ GPU-resident data
class GpuEpisodeStore:
    """All features, labels and episode indices resident on the GPU.

    Nothing in the training or evaluation loop touches the host except the final
    metric reduction.
    """

    def __init__(self, corpus: Corpus, prior: PooledPrior) -> None:
        x = corpus.feat
        mean = torch.tensor(prior.mean, dtype=torch.float32, device=DEVICE)
        scale = torch.tensor(prior.scale, dtype=torch.float32, device=DEVICE)
        feats = torch.tensor(x, dtype=torch.float32, device=DEVICE)
        self.g = (feats - mean) / scale                       # (P, D) standardised
        self.f0 = (self.g @ prior.beta + prior.intercept)     # (P,) pooled prior
        del feats
        self.dim = self.g.shape[1]
        self.parent_index = corpus.parent_index

        # per (target, parent) label and residual, as flat GPU tensors
        rows = corpus.rows
        self.row_parent = torch.tensor(
            [corpus.parent_index[p] for p in rows.compound_parent_uid],
            dtype=torch.long, device=DEVICE)
        self.row_y = torch.tensor(rows.pKi.to_numpy(), dtype=torch.float32, device=DEVICE)
        self.row_resid = self.row_y - self.f0[self.row_parent]
        self.row_year = torch.tensor(rows.document_year.to_numpy(),
                                     dtype=torch.long, device=DEVICE)

        # target -> row slice
        codes, uniques = pd.factorize(rows.target_uid, sort=True)
        self.target_ids = list(uniques)
        self.target_of_row = torch.tensor(codes, dtype=torch.long, device=DEVICE)
        order = np.argsort(codes, kind="stable")
        self.row_order = torch.tensor(order, dtype=torch.long, device=DEVICE)
        counts = np.bincount(codes, minlength=len(uniques))
        starts = np.concatenate([[0], np.cumsum(counts)[:-1]])
        self.slice_start = torch.tensor(starts, dtype=torch.long, device=DEVICE)
        self.slice_len = torch.tensor(counts, dtype=torch.long, device=DEVICE)
        self._doc = rows.document_uid.to_numpy()
        self._year = rows.document_year.to_numpy()
        self._parent = rows.compound_parent_uid.to_numpy()
        self._codes = codes

    def target_rows(self, target: str) -> np.ndarray:
        idx = self.target_ids.index(target)
        return np.flatnonzero(self._codes == idx)


def build_source_episodes(store: GpuEpisodeStore, source_targets: list[str],
                          rng: np.random.Generator, k: int, protocol: str,
                          n_per_target: int, max_query: int = MAX_QUERY
                          ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Pre-materialise episodes as GPU index tensors.

    Returns (support_rows (E,k), query_rows (E,Q), query_mask (E,Q)) where entries
    index the flat row tables. Document-ordered protocol: support strictly precedes
    query in publication year and shares no document.
    """
    sup_all, qry_all, mask_all = [], [], []
    for target in source_targets:
        rows = store.target_rows(target)
        if len(rows) < k + 5:
            continue
        years = store._year[rows]
        docs = store._doc[rows]
        uniq_years = np.unique(years)
        for _ in range(n_per_target):
            if protocol == "ordered":
                if len(uniq_years) < 2:
                    continue
                split = uniq_years[rng.integers(0, len(uniq_years) - 1)]
                pre = rows[years <= split]
                post = rows[(years > split) & ~np.isin(docs, np.unique(docs[years <= split]))]
                if len(pre) < k or len(post) < 5:
                    continue
                sup = rng.choice(pre, size=k, replace=False)
                pool = post[~np.isin(store._parent[post], store._parent[sup])]
            else:
                perm = rng.permutation(rows)
                sup, pool = perm[:k], perm[k:]
            if len(pool) < 5:
                continue
            q = rng.choice(pool, size=min(max_query, len(pool)), replace=False)
            pad = np.full(max_query, q[0])
            pad[:len(q)] = q
            mask = np.zeros(max_query, dtype=np.float32)
            mask[:len(q)] = 1.0
            sup_all.append(sup)
            qry_all.append(pad)
            mask_all.append(mask)
    if not sup_all:
        raise SystemExit("no source episodes could be built")
    return (torch.tensor(np.array(sup_all), dtype=torch.long, device=DEVICE),
            torch.tensor(np.array(qry_all), dtype=torch.long, device=DEVICE),
            torch.tensor(np.array(mask_all), dtype=torch.float32, device=DEVICE))


# ============================================================ the operator
class LigandEncoder(nn.Module):
    def __init__(self, dim: int, embed: int = EMBED_DIM, hidden: int = HIDDEN_DIM) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(dim, hidden), nn.GELU(), nn.LayerNorm(hidden), nn.Dropout(DROPOUT),
            nn.Linear(hidden, embed), nn.GELU(), nn.LayerNorm(embed))

    def forward(self, g: torch.Tensor) -> torch.Tensor:
        return self.net(g)


class MetaDeepKernelAdapter(nn.Module):
    """A2S-MDK: the adaptation operator IS a meta-learned deep-kernel ridge solve.

    Measured motivation.  A plain ridge head fitted on the k support residuals is
    the only arm that shows target-specific ranking (destroyed by wrong-target
    support and by permuting the support labels).  Its adaptation is
        w_r = argmin_w ||Phi_S w - (r_S - b_r)||^2 + lambda ||w||^2,
    which in dual form needs only the kernel.  So instead of replacing that
    estimator with a hypernetwork, A2S-MDK keeps it as the inner solve and
    meta-learns the space it is solved in:

        local_r(q) = k(q,S) (K(S,S) + lambda I)^-1 (r_S - b_r)
        k(a,b)     = (1-alpha) <g_a,g_b>/D + alpha <phi(g_a),phi(g_b)>/E

    phi, alpha, lambda and the budget gate are trained end to end on source
    episodes by differentiating through the solve.  The raw-feature kernel is
    retained as a residual branch so the operator starts from the identified
    estimator and the meta-learned part can only add to it.

    lambda is the identifiability budget: the realised effective degrees of
    freedom tr(K(K+lambda I)^-1) is bounded by k by construction, which is what
    the tau_z/sigma analysis says the support can afford.
    """

    def __init__(self, dim: int, embed: int = EMBED_DIM) -> None:
        super().__init__()
        self.encoder = LigandEncoder(dim, embed)
        self.dim, self.embed = dim, embed
        self.logit_alpha = nn.Parameter(torch.tensor(-2.0))    # start near the raw kernel
        # lambda is applied to a TRACE-NORMALISED kernel so it is scale-free; the
        # init corresponds to near-interpolation, matching the ridge head that was
        # measured to carry the target-specific signal.
        self.log_lambda = nn.Parameter(torch.tensor(-4.6))
        self.budget = nn.Sequential(nn.Linear(5, 64), nn.GELU(), nn.Linear(64, 1))
        nn.init.zeros_(self.budget[-1].weight)
        nn.init.constant_(self.budget[-1].bias, 2.0)           # start ~ fully open

    def kernel(self, a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
        alpha = torch.sigmoid(self.logit_alpha)
        raw = torch.bmm(a, b.transpose(1, 2)) / self.dim
        pa = self.encoder(a.reshape(-1, self.dim)).reshape(a.shape[0], a.shape[1], -1)
        pb = self.encoder(b.reshape(-1, self.dim)).reshape(b.shape[0], b.shape[1], -1)
        deep = torch.bmm(pa, pb.transpose(1, 2)) / self.embed
        return (1 - alpha) * raw + alpha * deep

    def _budget(self, g_sup: torch.Tensor, r_sup: torch.Tensor) -> torch.Tensor:
        b, k, _ = g_sup.shape
        zs = F.normalize(g_sup, dim=-1)
        sim = torch.bmm(zs, zs.transpose(1, 2))
        off = (sim.sum((1, 2)) - k) / max(k * (k - 1), 1)
        stats = torch.stack([
            torch.full((b,), float(k), device=g_sup.device),
            r_sup.std(dim=1) if k > 1 else torch.zeros(b, device=g_sup.device),
            r_sup.mean(dim=1), off, g_sup.norm(dim=-1).mean(dim=1) / self.dim ** 0.5,
        ], dim=1)
        return torch.sigmoid(self.budget(stats)).squeeze(-1)

    def forward(self, g_sup: torch.Tensor, r_sup: torch.Tensor, g_qry: torch.Tensor,
                anchor: torch.Tensor, use_budget: bool = True,
                use_deep: bool = True) -> torch.Tensor:
        b, k, _ = g_sup.shape
        if use_deep:
            k_ss = self.kernel(g_sup, g_sup)
            k_qs = self.kernel(g_qry, g_sup)
        else:
            k_ss = torch.bmm(g_sup, g_sup.transpose(1, 2)) / self.dim
            k_qs = torch.bmm(g_qry, g_sup.transpose(1, 2)) / self.dim
        # trace-normalise so that lambda has the same meaning at every k and for
        # both kernel branches
        scale = torch.diagonal(k_ss, dim1=1, dim2=2).mean(-1).clamp_min(1e-6)
        k_ss = k_ss / scale[:, None, None]
        k_qs = k_qs / scale[:, None, None]
        lam = F.softplus(self.log_lambda) + 1e-4
        eye = torch.eye(k, device=g_sup.device).unsqueeze(0).expand(b, -1, -1)
        alpha = torch.linalg.solve(k_ss + lam * eye,
                                   (r_sup - anchor.unsqueeze(1)).unsqueeze(-1))
        local = torch.bmm(k_qs, alpha).squeeze(-1)
        if use_budget:
            local = self._budget(g_sup, r_sup).unsqueeze(1) * local
        return anchor.unsqueeze(1) + local

    def effective_dof(self, g_sup: torch.Tensor) -> torch.Tensor:
        b, k, _ = g_sup.shape
        k_ss = self.kernel(g_sup, g_sup)
        scale = torch.diagonal(k_ss, dim1=1, dim2=2).mean(-1).clamp_min(1e-6)
        k_ss = k_ss / scale[:, None, None]
        lam = F.softplus(self.log_lambda) + 1e-4
        eye = torch.eye(k, device=g_sup.device).unsqueeze(0).expand(b, -1, -1)
        return torch.diagonal(torch.bmm(k_ss, torch.linalg.inv(k_ss + lam * eye)),
                              dim1=1, dim2=2).sum(-1)


class SupportConditionedAdaptationOperator(nn.Module):
    """T_phi: support set -> adaptation state -> recipient predictor.

    Three trainable paths, all conditioned on the support set and nothing else:

      1. set encoder      : permutation-invariant attention over support tokens
                            (h(d_i), r_i) -> task context c
      2. hypernetwork     : c -> FiLM (gamma, beta) for the head + anchor shift
      3. cross-attention  : query attends over support tokens, reading the
                            OBSERVED residuals -- the identifiable local path

    The budget head reads support-only statistics and emits a scalar in [0,1]
    that scales the whole adaptation. At zero the model returns the shrunk-anchor
    prior exactly, which is the registered no-transfer fallback.
    """

    def __init__(self, dim: int, embed: int = EMBED_DIM) -> None:
        super().__init__()
        self.encoder = LigandEncoder(dim, embed)
        self.resid_embed = nn.Linear(1, embed)
        self.set_attn = nn.MultiheadAttention(embed, N_HEADS, dropout=DROPOUT,
                                              batch_first=True)
        self.set_norm = nn.LayerNorm(embed)
        self.pool_query = nn.Parameter(torch.randn(1, 1, embed) * 0.02)

        self.hyper = nn.Sequential(
            nn.Linear(embed, embed), nn.GELU(), nn.LayerNorm(embed),
            nn.Linear(embed, 2 * embed + 1))
        self.cross_attn = nn.MultiheadAttention(embed, N_HEADS, dropout=DROPOUT,
                                                batch_first=True)
        self.head = nn.Sequential(
            nn.Linear(embed, embed), nn.GELU(), nn.LayerNorm(embed), nn.Linear(embed, 1))
        self.local_head = nn.Linear(embed, 1)
        self.budget = nn.Sequential(
            nn.Linear(5, 64), nn.GELU(), nn.Linear(64, 1))
        for module in (self.head[-1], self.local_head):
            nn.init.zeros_(module.weight)
            nn.init.zeros_(module.bias)

    def adapt(self, g_sup: torch.Tensor, r_sup: torch.Tensor
              ) -> tuple[torch.Tensor, ...]:
        """Produce the adaptation state from the support set alone."""
        b, k, _ = g_sup.shape
        h_sup = self.encoder(g_sup.reshape(b * k, -1)).reshape(b, k, -1)
        tokens = self.set_norm(h_sup + self.resid_embed(r_sup.unsqueeze(-1)))
        ctx, _ = self.set_attn(self.pool_query.expand(b, -1, -1), tokens, tokens)
        ctx = ctx.squeeze(1)                                        # (B, E)
        params = self.hyper(ctx)
        embed = h_sup.shape[-1]
        gamma = 1.0 + params[:, :embed]
        beta = params[:, embed:2 * embed]
        anchor_shift = params[:, -1]

        # support-only statistics for the identifiability budget
        sim = F.normalize(h_sup, dim=-1) @ F.normalize(h_sup, dim=-1).transpose(1, 2)
        off = (sim.sum((1, 2)) - k) / max(k * (k - 1), 1)
        stats = torch.stack([
            torch.full((b,), float(k), device=g_sup.device),
            r_sup.std(dim=1) if k > 1 else torch.zeros(b, device=g_sup.device),
            r_sup.mean(dim=1),
            off,                                                    # support redundancy
            h_sup.norm(dim=-1).mean(dim=1),
        ], dim=1)
        budget = torch.sigmoid(self.budget(stats)).squeeze(-1)      # (B,) in [0,1]
        return h_sup, tokens, gamma, beta, anchor_shift, budget

    def forward(self, g_sup: torch.Tensor, r_sup: torch.Tensor, g_qry: torch.Tensor,
                anchor: torch.Tensor, use_budget: bool = True,
                use_cross: bool = True, use_film: bool = True) -> torch.Tensor:
        b, q, _ = g_qry.shape
        h_sup, tokens, gamma, beta, shift, budget = self.adapt(g_sup, r_sup)
        h_qry = self.encoder(g_qry.reshape(b * q, -1)).reshape(b, q, -1)

        modulated = h_qry * gamma.unsqueeze(1) + beta.unsqueeze(1) if use_film else h_qry
        out = self.head(modulated).squeeze(-1)                      # (B, Q)

        if use_cross:
            values = tokens * r_sup.unsqueeze(-1)
            ctx, _ = self.cross_attn(h_qry, tokens, values)
            out = out + self.local_head(ctx).squeeze(-1)

        scale = budget.unsqueeze(1) if use_budget else torch.ones_like(budget).unsqueeze(1)
        return anchor.unsqueeze(1) + scale * (out + shift.unsqueeze(1))


# ============================================================ training
def shrunk_anchor(r_sup: torch.Tensor, hyper: Hyper) -> torch.Tensor:
    """Batched James-Stein anchor: tau_b << sigma makes the raw support mean unusable."""
    k = r_sup.shape[1]
    prec = 1.0 / max(hyper.tau_b ** 2, 1e-8) + k / (hyper.sigma ** 2)
    return r_sup.sum(dim=1) / (hyper.sigma ** 2) / prec


def ranking_loss(pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor,
                 per_episode: bool = False) -> torch.Tensor:
    """Pairwise logistic loss; CI and Spearman are the reported endpoints."""
    dp = pred.unsqueeze(2) - pred.unsqueeze(1)
    dt = target.unsqueeze(2) - target.unsqueeze(1)
    pair = mask.unsqueeze(2) * mask.unsqueeze(1) * (dt.abs() > 1e-6).float()
    loss = F.softplus(-torch.sign(dt) * dp) * pair
    if per_episode:
        return loss.sum((1, 2)) / pair.sum((1, 2)).clamp_min(1.0)
    if pair.sum() == 0:
        return pred.sum() * 0.0
    return loss.sum() / pair.sum()


def train_operator(store: GpuEpisodeStore, episodes: tuple[torch.Tensor, ...],
                   hyper: Hyper, rng: np.random.Generator, steps: int = STEPS,
                   log: list[str] | None = None, kind: str = 'mdk',
                   spec_weight: float = SPEC_WEIGHT):
    sup_idx, qry_idx, qry_mask = episodes
    model = (MetaDeepKernelAdapter(store.dim) if kind == 'mdk'
             else SupportConditionedAdaptationOperator(store.dim)).to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    sched = torch.optim.lr_scheduler.OneCycleLR(
        opt, max_lr=LR, total_steps=steps, pct_start=min(0.3, WARMUP / max(steps, 1)))
    scaler = torch.amp.GradScaler("cuda", enabled=DEVICE == "cuda")
    n_ep = sup_idx.shape[0]
    model.train()
    for step in range(steps):
        sel = torch.randint(0, n_ep, (BATCH_EPISODES,), device=DEVICE)
        s_rows, q_rows, mask = sup_idx[sel], qry_idx[sel], qry_mask[sel]
        g_sup = store.g[store.row_parent[s_rows]]
        g_qry = store.g[store.row_parent[q_rows]]
        r_sup = store.row_resid[s_rows]
        r_qry = store.row_resid[q_rows]
        anchor = shrunk_anchor(r_sup, hyper)

        # COUNTERFACTUAL SUPPORT CONTRAST: pair every query set with a support set
        # from a different source target. Query loss alone is dominated by the
        # target-INDEPENDENT ligand term, so an operator trained on it degenerates
        # into support-agnostic generic ranking (measured: wrong-target support
        # ranked as well as correct support). Requiring the correct support to
        # beat its own mismatched counterfactual makes target-specificity the
        # trained objective rather than a hoped-for by-product.
        shift = 1 + int(torch.randint(0, BATCH_EPISODES - 1, (1,)).item())
        roll = torch.roll(torch.arange(BATCH_EPISODES, device=DEVICE), shift)
        g_cf, r_cf = g_sup[roll], r_sup[roll]
        anchor_cf = shrunk_anchor(r_cf, hyper)

        with torch.amp.autocast("cuda", enabled=DEVICE == "cuda"):
            pred = model(g_sup, r_sup, g_qry, anchor)
            mse = (((pred - r_qry) ** 2) * mask).sum() / mask.sum()
            rank = ranking_loss(pred, r_qry, mask)
            fit = mse + RANK_LOSS_WEIGHT * rank

            pred_cf = model(g_cf, r_cf, g_qry, anchor_cf)
            rank_ep = ranking_loss(pred, r_qry, mask, per_episode=True)
            rank_cf_ep = ranking_loss(pred_cf, r_qry, mask, per_episode=True)
            mse_ep = (((pred - r_qry) ** 2) * mask).sum(1) / mask.sum(1).clamp_min(1.0)
            mse_cf_ep = (((pred_cf - r_qry) ** 2) * mask).sum(1) / mask.sum(1).clamp_min(1.0)
            contrast = (F.relu(SPEC_MARGIN + rank_ep - rank_cf_ep).mean()
                        + F.relu(SPEC_MARGIN + mse_ep - mse_cf_ep).mean())
            loss = fit + spec_weight * contrast

        opt.zero_grad(set_to_none=True)
        scaler.scale(loss).backward()
        scaler.unscale_(opt)
        torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
        scaler.step(opt)
        scaler.update()
        sched.step()
        if log is not None and (step % 500 == 0 or step == steps - 1):
            log.append(f"    {kind.upper()} step {step:>5} loss {float(loss):.4f} mse {float(mse):.4f} "
                       f"rank {float(rank):.4f} contrast {float(contrast):.4f} "
                       f"spec-gap(rank) {float((rank_cf_ep - rank_ep).mean()):+.4f}")
    return model.eval()


# ============================================================ evaluation
STRATA = ("all", "near", "mid", "far")


def _prepare_recipients(corpus: Corpus, store: GpuEpisodeStore):
    """Precompute every index/label tensor ONCE, on the GPU.

    Pandas lookups inside the episode loop were the dominant cost and the reason
    GPU utilisation collapsed during evaluation.
    """
    label = corpus.rows.set_index(["target_uid", "compound_parent_uid"]).pKi
    query_by_target = {t: list(f.compound_parent_uid)
                       for t, f in corpus.query.groupby("target_uid")}
    draws_by_target = {t: f for t, f in corpus.draws.groupby("target_uid")}
    prepared = {}
    for t in corpus.recipient_targets:
        qp = query_by_target.get(t, [])
        if len(qp) < 3:
            continue
        qi = torch.tensor([corpus.parent_index[x] for x in qp], dtype=torch.long, device=DEVICE)
        yq = torch.tensor([label.loc[(t, x)] for x in qp], dtype=torch.float32, device=DEVICE)
        per_k = {}
        frame = draws_by_target[t]
        for k in SUPPORT_K:
            fk = frame[frame.k == k]
            sup, ys = [], []
            for _, sub in fk.groupby("draw_id"):
                pars = list(sub.compound_parent_uid)
                sup.append([corpus.parent_index[x] for x in pars])
                ys.append([float(label.loc[(t, x)]) for x in pars])
            if sup:
                per_k[k] = (torch.tensor(sup, dtype=torch.long, device=DEVICE),
                            torch.tensor(ys, dtype=torch.float32, device=DEVICE))
        prepared[t] = {"q_idx": qi, "y_q": yq, "draws": per_k}
    return prepared


def _similarity_strata(store: GpuEpisodeStore, q_idx: torch.Tensor,
                       s_idx: torch.Tensor) -> dict:
    """Max cosine similarity of each query compound to the support set."""
    zq = F.normalize(store.g[q_idx], dim=-1)
    zs = F.normalize(store.g[s_idx], dim=-1)
    sim = (zq @ zs.T).max(dim=1).values.cpu().numpy()
    return sim


@torch.no_grad()
def evaluate(corpus: Corpus, out: Path, protocol: str = "ordered",
             steps: int = STEPS, quick: bool = False) -> dict:
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    torch.manual_seed(SEED)
    log = []
    if DEVICE == "cuda":
        torch.cuda.reset_peak_memory_stats()

    src_rows = corpus.rows[corpus.rows.target_uid.isin(corpus.source_targets)]
    with torch.enable_grad():
        from research.a2s.a2s_bir import feature_matrix
        x_src = feature_matrix(corpus, list(src_rows.compound_parent_uid))
        prior = PooledPrior().fit(x_src, src_rows.pKi.to_numpy(),
                                  src_rows.target_uid.to_numpy())
    del x_src
    store = GpuEpisodeStore(corpus, prior)
    hyper = estimate_hyper(corpus, prior, corpus.source_targets)
    log.append("  hierarchy: tau_b=%.4f sigma=%.4f" % (hyper.tau_b, hyper.sigma))

    n_per = 4 if quick else 24
    eps = build_source_episodes(store, corpus.source_targets, rng, max(SUPPORT_K),
                                protocol, n_per)
    log.append("  %d source episodes materialised on %s" % (eps[0].shape[0], DEVICE))
    n_steps = 200 if quick else steps
    with torch.enable_grad():
        model = train_operator(store, eps, hyper, rng, n_steps, log, kind='mdk')
        model_nospec = train_operator(store, eps, hyper, rng, n_steps, log,
                                      kind='mdk', spec_weight=0.0)
        model_scao = train_operator(store, eps, hyper, rng, n_steps, log, kind='scao')
        grad_models = {}
        if not quick:
            for mode in ("maml", "anil"):
                grad_models[mode] = train_gradient_meta(
                    corpus, prior, corpus.source_targets, rng, protocol, mode, log=log)

    prepared = _prepare_recipients(corpus, store)
    arms = ["recipient_calibration", "f0_only", "f0_anchor_shrunk",
            "pooled_finetune", "pooled_finetune_wrong", "pooled_finetune_perm",
            "a2s_mdk", "a2s_mdk_wrong", "a2s_mdk_perm",
            "a2s_mdk_nospec", "a2s_mdk_nospec_wrong",
            "a2s_mdk_nodeep", "a2s_mdk_nobudget",
            "a2s_scao", "a2s_scao_wrong"]
    if not quick:
        arms += ["maml", "anil"]

    # arm -> (correct arm, wrong-support counterpart) for the target-specificity test
    SPECIFICITY = {"a2s_mdk": "a2s_mdk_wrong", "a2s_mdk_nospec": "a2s_mdk_nospec_wrong",
                   "a2s_scao": "a2s_scao_wrong",
                   "pooled_finetune": "pooled_finetune_wrong"}
    PERMUTED = {"a2s_mdk": "a2s_mdk_perm", "pooled_finetune": "pooled_finetune_perm"}

    results, diagnostics = {}, {}
    targets = list(prepared)
    for k in SUPPORT_K:
        # stratum -> arm -> target -> metric
        acc = {st: {a: {} for a in arms} for st in STRATA}
        budgets, dofs = [], []
        for ti, t in enumerate(targets):
            rec = prepared[t]
            if k not in rec["draws"]:
                continue
            q_idx, y_q = rec["q_idx"], rec["y_q"]
            s_idx, y_s = rec["draws"][k]
            n_draw = s_idx.shape[0]
            g_q = store.g[q_idx].unsqueeze(0).expand(n_draw, -1, -1)
            f0_q = store.f0[q_idx]
            g_s = store.g[s_idx]
            r_s = y_s - store.f0[s_idx]
            anchor = shrunk_anchor(r_s, hyper)

            other = targets[(ti + 7) % len(targets)]
            o_idx, o_y = prepared[other]["draws"].get(k, (s_idx, y_s))
            o_idx, o_y = o_idx[:n_draw], o_y[:n_draw]
            r_o = o_y - store.f0[o_idx]
            g_o = store.g[o_idx]
            anchor_o = shrunk_anchor(r_o, hyper)

            perm = torch.stack([r[torch.randperm(k, device=DEVICE)] for r in r_s])

            pred = {}
            pred["recipient_calibration"] = y_s.mean(1, keepdim=True).expand(-1, len(q_idx))
            pred["f0_only"] = f0_q.unsqueeze(0).expand(n_draw, -1)
            pred["f0_anchor_shrunk"] = f0_q.unsqueeze(0) + anchor.unsqueeze(1)
            pred["a2s_mdk"] = f0_q.unsqueeze(0) + model(g_s, r_s, g_q, anchor)
            pred["a2s_mdk_wrong"] = f0_q.unsqueeze(0) + model(g_o, r_o, g_q, anchor_o)
            pred["a2s_mdk_perm"] = f0_q.unsqueeze(0) + model(g_s, perm, g_q,
                                                             shrunk_anchor(perm, hyper))
            pred["a2s_mdk_nodeep"] = f0_q.unsqueeze(0) + model(g_s, r_s, g_q, anchor,
                                                               use_deep=False)
            pred["a2s_mdk_nobudget"] = f0_q.unsqueeze(0) + model(g_s, r_s, g_q, anchor,
                                                                 use_budget=False)
            pred["a2s_mdk_nospec"] = f0_q.unsqueeze(0) + model_nospec(g_s, r_s, g_q, anchor)
            pred["a2s_mdk_nospec_wrong"] = f0_q.unsqueeze(0) + model_nospec(
                g_o, r_o, g_q, anchor_o)
            pred["a2s_scao"] = f0_q.unsqueeze(0) + model_scao(g_s, r_s, g_q, anchor)
            pred["a2s_scao_wrong"] = f0_q.unsqueeze(0) + model_scao(g_o, r_o, g_q, anchor_o)
            budgets.append(float(model._budget(g_s, r_s).mean()))
            dofs.append(float(model.effective_dof(g_s).mean()))

            eye = torch.eye(store.dim, device=DEVICE)
            for name, gg, rr in (("pooled_finetune", g_s, r_s),
                                 ("pooled_finetune_wrong", g_o, r_o),
                                 ("pooled_finetune_perm", g_s, perm)):
                gram = torch.baddbmm(eye.expand(n_draw, -1, -1),
                                     gg.transpose(1, 2), gg)
                delta = torch.linalg.solve(gram, torch.bmm(gg.transpose(1, 2),
                                                           rr.unsqueeze(-1)))
                pred[name] = f0_q.unsqueeze(0) + torch.bmm(g_q, delta).squeeze(-1)

            if not quick:
                with torch.enable_grad():
                    for mode in ("maml", "anil"):
                        mdl, ilr = grad_models[mode]
                        rows = [adapt_gradient_meta(mdl, g_s[d], r_s[d], g_q[d], mode, ilr)
                                for d in range(n_draw)]
                        pred[mode] = f0_q.unsqueeze(0) + torch.tensor(
                            np.array(rows), dtype=torch.float32, device=DEVICE)

            y_np = y_q.cpu().numpy()
            sim = _similarity_strata(store, q_idx, s_idx.reshape(-1))
            lo, hi = np.quantile(sim, [1 / 3, 2 / 3])
            masks = {"all": np.ones(len(sim), bool), "near": sim >= hi,
                     "mid": (sim > lo) & (sim < hi), "far": sim <= lo}
            for arm in arms:
                block = pred[arm].float().cpu().numpy()
                for st, m in masks.items():
                    if m.sum() < 3:
                        continue
                    per_draw = [episode_metrics(y_np[m], block[d][m]) for d in range(n_draw)]
                    acc[st][arm][t] = {mm: float(np.nanmean([e[mm] for e in per_draw]))
                                       for mm in METRIC_NAMES}

        strata_out = {}
        for st in STRATA:
            summary = {a: {mm: float(np.nanmean([v[mm] for v in acc[st][a].values()]))
                           for mm in METRIC_NAMES} for a in arms if acc[st][a]}
            gains = {}
            for reference in ("f0_anchor_shrunk", "recipient_calibration", "pooled_finetune"):
                ref = acc[st].get(reference, {})
                blk = {}
                for arm in arms:
                    if arm == reference or not acc[st][arm] or not ref:
                        continue
                    shared = [t for t in acc[st][arm] if t in ref]
                    entry = {}
                    for mm in ("ci", "spearman", "ndcg10", "rmse", "rm2", "pearson"):
                        sign = -1.0 if mm == "rmse" else 1.0
                        d = {t: sign * (acc[st][arm][t][mm] - ref[t][mm]) for t in shared}
                        entry[mm + "_gain"] = component_bootstrap(d, corpus.components, rng,
                                                                  N_BOOTSTRAP)
                    blk[arm] = entry
                gains[reference] = blk

            # ---- PRIMARY: target-specific ranking gain -------------------
            specificity = {}
            for arm, wrong in SPECIFICITY.items():
                if not acc[st].get(arm) or not acc[st].get(wrong):
                    continue
                shared = [t for t in acc[st][arm] if t in acc[st][wrong]]
                entry = {}
                for mm in ("ci", "spearman", "ndcg10"):
                    d = {t: acc[st][arm][t][mm] - acc[st][wrong][t][mm] for t in shared}
                    entry["vs_wrong_support_" + mm] = component_bootstrap(
                        d, corpus.components, rng, N_BOOTSTRAP)
                perm_arm = PERMUTED.get(arm)
                if perm_arm and acc[st].get(perm_arm):
                    shared_p = [t for t in acc[st][arm] if t in acc[st][perm_arm]]
                    for mm in ("ci", "spearman", "ndcg10"):
                        d = {t: acc[st][arm][t][mm] - acc[st][perm_arm][t][mm]
                             for t in shared_p}
                        entry["vs_permuted_labels_" + mm] = component_bootstrap(
                            d, corpus.components, rng, N_BOOTSTRAP)
                specificity[arm] = entry
            strata_out[st] = {"summary": summary, "gains": gains,
                              "target_specific_ranking": specificity,
                              "n_targets": len(acc[st]["f0_anchor_shrunk"])}

        results[k] = strata_out
        diagnostics[k] = {"recipients_evaluated": len(acc["all"]["f0_anchor_shrunk"]),
                          "mean_learned_budget": float(np.mean(budgets)) if budgets else 0.0,
                          "mean_effective_dof": float(np.mean(dofs)) if dofs else 0.0,
                          "kernel_alpha": float(torch.sigmoid(model.logit_alpha)),
                          "kernel_lambda": float(F.softplus(model.log_lambda) + 1e-4)}
        ts = strata_out["all"]["target_specific_ranking"].get("a2s_mdk", {})
        ci = ts.get("vs_wrong_support_ci", {})
        print("  k=%d: %d recipients, budget %.3f, target-specific CI gain %s"
              % (k, diagnostics[k]["recipients_evaluated"],
                 diagnostics[k]["mean_learned_budget"],
                 ("%+.4f [%+.4f,%+.4f]" % (ci.get("mean", float("nan")),
                                            ci.get("lo", float("nan")),
                                            ci.get("hi", float("nan"))))
                 if ci else "n/a"))

    report = {
        "schema": "a2s-scao-evaluation-v2",
        "model": "A2S-MDK: meta-learned deep-kernel support-conditioned adaptation operator",
        "primary_estimand": "target-specific within-target ranking gain = "
                            "metric(correct support) - metric(wrong-target support), "
                            "paired per recipient, bootstrapped over homology components",
        "endpoint": "pKi", "protocol": protocol, "seed": SEED, "device": DEVICE,
        "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "roster": {"recipients": len(corpus.recipient_targets),
                   "sources": len(corpus.source_targets),
                   "components": len({corpus.components[t] for t in corpus.recipient_targets})},
        "training": {"steps": 200 if quick else steps, "batch_episodes": BATCH_EPISODES,
                     "source_episodes": int(eps[0].shape[0]),
                     "parameters": sum(p.numel() for p in model.parameters()),
                     "operator": "MetaDeepKernelAdapter",
                     "rank_loss_weight": RANK_LOSS_WEIGHT,
                     "spec_weight": SPEC_WEIGHT, "spec_margin": SPEC_MARGIN},
        "hierarchy": {"tau_b": round(hyper.tau_b, 4), "sigma": round(hyper.sigma, 4)},
        "diagnostics": {str(k): v for k, v in diagnostics.items()},
        "strata": list(STRATA),
        "results": {str(k): v for k, v in results.items()},
        "metrics_reported": METRIC_NAMES,
        "statistical_unit": "independent homology component (paired bootstrap, 5000)",
        "wall_time_s": round(time.time() - t0, 1),
        "peak_torch_mem_mb": (round(torch.cuda.max_memory_allocated() / 2 ** 20, 1)
                              if torch.cuda.is_available() else None),
        "log": log,
    }
    out.mkdir(parents=True, exist_ok=True)
    path = out / ("a2s_mdk_%s_seed%d.json" % (protocol, SEED))
    path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print("\nwrote %s" % path)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="A2S-SCAO training and evaluation")
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).resolve().parents[2] / "reports" / "active")
    parser.add_argument("--protocol", choices=("ordered", "random"), default="ordered")
    parser.add_argument("--steps", type=int, default=STEPS)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--roster", type=Path, default=None)
    args = parser.parse_args()
    corpus = load_corpus(args.roster)
    print(f"corpus: {len(corpus.source_targets)} sources, "
          f"{len(corpus.recipient_targets)} recipients, device={DEVICE}")
    report = evaluate(corpus, args.out, args.protocol, args.steps, args.quick)
    key = ["rmse", "ci", "rm2", "spearman", "pearson", "ndcg10", "aupr"]
    for k, block in report["results"].items():
        allb = block["all"]
        print("\n=== k=%s | all query compounds (n_targets=%d) ===" % (k, allb["n_targets"]))
        print("%-26s" % "arm" + "".join("%9s" % m for m in key))
        for arm, row in allb["summary"].items():
            print("%-26s" % arm + "".join("%9.4f" % row[m] for m in key))
        print("  PRIMARY target-specific ranking gain (correct - wrong-target support):")
        for arm, ent in allb["target_specific_ranking"].items():
            for name, v in ent.items():
                flag = "*" if v["lo"] > 0 else ""
                print("    %-24s %-28s %+.4f [%+.4f, %+.4f]%s"
                      % (arm, name, v["mean"], v["lo"], v["hi"], flag))
        for st in ("near", "far"):
            sb = block.get(st)
            if not sb or not sb["target_specific_ranking"]:
                continue
            print("  stratum=%s (support-query similarity):" % st)
            for arm, ent in sb["target_specific_ranking"].items():
                v = ent.get("vs_wrong_support_ci")
                if v:
                    flag = "*" if v["lo"] > 0 else ""
                    print("    %-24s CI gain %+.4f [%+.4f, %+.4f]%s"
                          % (arm, v["mean"], v["lo"], v["hi"], flag))


if __name__ == "__main__":
    main()
