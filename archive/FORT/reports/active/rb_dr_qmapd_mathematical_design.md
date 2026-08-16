# RB-DR-QMAPD — mathematical design (2026-07-26)

Source of truth for the estimand, the oracle decomposition, the teacher family, and the (gated)
completion-marginalized / design-robust / student mathematics. Written before the O1 result. Downstream
sections (Stages D–one-seed) are DESIGN ONLY and are not implemented unless O1 and O2 pass. Novelty and
prior art: `reports/active/rb_dr_qmapd_novelty_audit.md`.

## 1. Frozen estimand

```
y(t,d) = b(t,d) + c_t + g_t(d) + eps(t,d)
```
`b(t,d)` frozen zero-shot base; `c_t` scalar target calibration; `g_t(d)` the target-specific
within-target SAR ranking function (the object of interest); `eps` noise. `g_t` must exclude target
intercept, global ligand potency, target-only effects, and ligand-only fixed effects — i.e. it lives in
the **interaction quotient**. The exact weighted-Hodge projector `M_X^W` (`research/orrc_v2.py`,
audited) defines and fits the *train-only teacher* interaction estimand `r_int = M_X^W y`. Predictions
are **never** Hodge-projected using the inference query set; the student ranking output is **pointwise**
and independent of which unrelated queries share a batch.

## 2. Base freeze (Section 4)

The matched base question (ligand-only vs protein-conditioned vs shuffled/random/additive) has been
answered NEGATIVE by the existing record and is not re-run (adding compute after a settled negative is
inadmissible): Gate PB (`panel_gate_pb.json`) — `CFRI−B0 = +0.0180 [−0.0145,+0.0534]`, the jointly
trained interaction `I0 = +0.0376 [LCB −0.0033]` (under-powered, not significant), and the same-feature
zero-shot head does not clear its controls; `HQGBMA_STAGE_D_FAIL_STOP` — frozen ESM predicts no
transferable interaction subspace, and `PROTEIN_CONDITIONED_PRIOR_NOT_LOAD_BEARING`. No protein-
conditioned base satisfies the admissibility criteria. Therefore **`b(t,d) = b(d)`** (ligand-only,
heteroscedastic, cross-fitted OOF), and the protein path is removed from RB-DR-QMAPD. ESM/cross-attention
are not retained for cosmetic conventionality.

## 3. Oracle decomposition (Stage O1 — the decisive gate)

Same frozen teacher at both support sizes; identical rows:
```
Delta_info(k,K)  = rho(T_K) − rho(T_k)      # value of MORE same-target evidence
Delta_arch(k)    = rho(T_k) − rho(B_k)      # value of the teacher at k over the baseline
Delta_total(k,K) = rho(T_K) − rho(B_k) = Delta_info + Delta_arch  (exact, identical rows)
```
`B_k = B0` = the frozen zero-shot ligand base (the strongest robustly-established non-teacher baseline;
the historical k-shot methods are the same posterior family, so relative to them `Delta_arch≈0` and
`Delta_total≈Delta_info` — using B0 only makes the total bar *easier*, while the decisive `Delta_info`
is invariant to this choice). `rho` = within-component (tunit) macro Spearman on the fixed dual-cold
query. Primary pair (4,16); secondary (4,32). Pass conditions and thresholds are frozen in
`rb_dr_qmapd_oracle_preregistration.md`; empirical MDE80 for info/total comes from a support-resampling
null (same teacher, two independent support draws).

## 4. Teacher family (Section 6)

Registered before O1: **Teacher A** — the existing validated exact-Cholesky Bayesian residual posterior
(`model/bayes_posterior.py` via `model/descriptor_baseline.py`), design `[1, z(rep_d)]` (intercept =
calibration `c_t`; linear = ranking `g_t`), meta-trained with support sizes drawn up to `K` so it is
competent at both `k` and `K` (fairness to `T_K`), then frozen; the closed-form posterior is identical
architecture at every support size. **Teacher B** — corrected hierarchical measurement-error model
`β̂_t|β_t ~ N(B_t β_t, V_t)`, `β_t ~ N(μ,Γ)`, `μ,Γ` by REML/EM with Helmert-transformed noise
`R̃ = H_k R_s H_k^T` (registered, implemented only if A/C are inconclusive). **Teacher C** — deep-kernel
GP on frozen ligand reps with a PD kernel, no target id / no query label (registered). Selection: one
family by nested inner components on a frozen score (component-macro Spearman + RMSE safety + proper
score + wrong-support specificity); no post-hoc maximum. For the decisive O1 the primary registered
teacher is A (validated, cheapest); B/C are the registered alternatives.

## 5. (GATED) Completion-marginalized teacher law — Stages D/§10

Implemented only after O1 and O2 pass. For legal label-blind completion designs
`V={ν_uniform, ν_scaffold_balanced, ν_distance_stratified, ν_assay_document_balanced}` and samples
`A_j^(m)`, the **Rao–Blackwellized design-specific teacher law** is `P̄_j = Σ_m w_jm P_jm` (weights
fixed by the design, sum 1). For Gaussian components `P_jm = N(μ_jm, Σ_jm)`, the **law of total
covariance** gives
```
μ̄_j = Σ_m w_jm μ_jm ,   Σ̄_j = Σ_m w_jm Σ_jm + Σ_m w_jm (μ_jm−μ̄_j)(μ_jm−μ̄_j)^T ,
```
within-completion + between-completion (missing-evidence) uncertainty; the **full** matrix is used, not
its diagonal. `M` is fixed by a train-only Monte-Carlo error audit. `P̄_j` is a completion-design-
marginalized law, **not** the true posterior `p(g_Q|S)` (grounds: Nearly-Optimal Bayesian Inference for
Structural Missingness, arXiv 2601.18500 — integrate over missingness, never plug-in impute).

## 6. (GATED) Design-robust objective — §11

`D_j = D(q_φ, P̄_j)`; entropic robust objective
`L_DR = τ log[(1/J) Σ_j exp(D_j/τ)]` (τ→0 ⇒ worst design; larger τ ⇒ average risk), τ from an inner-
component grid (grounds: DRO-NPE, arXiv 2605.28516). Controls: average / single-uniform / hard-worst /
entropic. Retained only if it improves worst-design behaviour or calibration without materially lowering
average ranking.

## 7. (GATED) Student — §12/§13

Pointwise, capacity-controlled set-to-function (grounds: Neural Operator Processes, arXiv 2606.22946):
permutation-invariant support encoder, query-to-support cross-attention (no query-to-query attention),
pointwise head producing `μ_i, v_i>0, u_i∈R^{r_cov}` per query; joint `Σ_Q = diag(v) + U_Q U_Q^T` so
adding unrelated queries never changes a query's prediction. No intercept in the ranking head; exact
`k=0` base fallback; exact zero ranking correction at `k=1`. Objective
`L = L_query + λ_DR L_DR + λ_vario L_variogram + λ_rank L_rank + λ_inv L_invariance`, `L_query`
(strictly proper score on genuinely held query labels) at coefficient 1 — the privileged teacher is
never the sole gradient direction (grounds: Self-Distilled RLVR leakage finding, arXiv 2604.03128); an
ablation uses the teacher only as a bounded stop-gradient weighting on `L_query`. Forbidden inputs:
target/cluster id, protein embedding as controller, privileged completion rows/labels, query labels,
split/source/document/assay id in the prediction path.

## 8. Claim discipline

No claim of true-posterior recovery, Bayes optimality, external validation, large improvement, general
DTA superiority, or SOTA is made unless the corresponding gate is satisfied. Mathematical identities
(the LTC decomposition, `Delta_total = Delta_info + Delta_arch`), completion-design teacher laws,
train-only mechanism evidence, internal predictive evidence, and independent confirmation evidence are
kept strictly separate. The Metz development rows are spent → any predictive number is
`ARCHITECTURE_MECHANISM_RESULT_ONLY`; confirmation/Davis/sealed labels are not read.
