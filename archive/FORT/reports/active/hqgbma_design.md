# HQ-GBMA — mathematical design report (2026-07-26)

**HQ-GBMA = Hodge-Quotient interaction Transformer + error-corrected Grassmann Bayesian Meta-Adapter.**
A redesign of the dual-cold few-shot DTA system in which every documented failure is turned into a
structural mathematical constraint, on a mature cross-attention DTA backbone. This report is the source
of truth for the estimand, derivations, identifiability analysis, posterior equations, and failure
cases. It is written before the code and before any Stage result is read.

## 0. Failures reverse-engineered into constraints

| Documented failure (history.md) | Root cause | Structural constraint adopted |
|---|---|---|
| CFRI soft centering / base-orthogonality could bias a real interaction; shuffled protein beat real | penalty ≠ projection; the interaction head can still contain main effects | **exact weighted Hodge quotient** `g = M_X^W h`, `XᵀW g = 0` to registered tolerance (Module 1) |
| BM0/BM0P/BM1-RR/PC: intercept/scale adaptation improves RMSE without reordering; "support-conditioned ligand kernel" | one posterior mixed calibration (mean) and ranking evidence | **two independent posteriors**; ranking uses Helmert contrasts, calibration a separate scalar (§6) |
| Protein-conditioned covariance prior *not* load-bearing (2026-07-26); shuffled/random protein ≥ true | a soft covariance orientation can be matched by any orientation; the correction is not confined | **hard Grassmann confinement**: correction ∈ `range(P_t)`, `P_t=R_t R_tᵀ` predicted from ESM (Module 2) |
| adapter treats `k(k−1)/2` pairs as independent; duplicate support inflates evidence | non-orthonormal contrast basis | **orthonormal Helmert** `H_k`, `H_k 1=0`, `H_k H_kᵀ=I_{k−1}` (§6) |
| noisy target coefficients treated as exact interaction functions | plug-in bias | **error-corrected containment** subtracts `tr((I−P_t)V_t)` (§5) and **error-corrected EB** marginalises `V_t` |
| effective-rank ceiling inherited as `r=8` | unjustified | **rank selected by train-only nested folds**, `r∈{3,4,6}` (§5) |

Prohibited rescues (full MAML, support-label tokens, learned adaptation gates, another
covariance-only adapter, ordinary shuffle-contrast as the core, Student-t, expert banks, AttnRes,
mandatory 3D/docking/PLIF, OSA-shrinker retuning, post-failure capacity/epoch/threshold increases) are
excluded by construction.

## 1. Formal estimand

Targets are homology components (the statistical unit). For an unseen target `t` (UniProt/accession,
homology-cluster, and pocket/family disjoint from train) and unseen ligand `d` (parent-connectivity,
Bemis–Murcko scaffold, and Tanimoto<0.95 disjoint from train and from the target's support), with `k`
few-shot support measurements `S_t={(d_i,y_i)}`, predict `y(t,d)`. The scientific target is the
**within-target ligand ordering** (paired target-component macro Spearman at `k=4`), which changes only
through a protein-specific interaction, not through target-level location/scale. Random-split metrics
are excluded: per Zhang et al., ICLR 2025 (Oral), random test sets are dominated by high-similarity
pairs; the firewalls below implement similarity-aware evaluation.

## 2. Functional-ANOVA decomposition and identifiability

On the observed target–ligand bipartite graph,
```
y(t,d) = a_t + b(d) + g(t,d) + eps .
```
`a_t` (target main effect) and `b(d)` (ligand main effect) cannot reorder ligands within a target;
`g(t,d)` (the interaction) is the only reordering term. Functional-ANOVA components are **not unique**
without an explicit orthogonality/sum-to-zero constraint (Limmer et al., MLSP 2025, Neural-ANOVA; Park
et al., ICML 2025, ANOVA-TPNN). On a **complete** grid the identifying condition is sum-to-zero. On an
**incomplete** panel it must be the **weighted observed-edge** condition:
```
X = [X_T , X_D]  (target/ligand incidence),   W  fixed train-only weights,
M_X^W = I − X (Xᵀ W X)^+ Xᵀ W,     admissible interaction ⇔  Xᵀ W g = 0 .
```
Separate two-sided feature centering is **not** observed-space orthogonality under missingness
(blueprint counter-example: 3×3 panel, one missing cell → cosine −0.293 with the interaction). The
weighted pseudo-inverse projector `M_X^W` is the exact fix and is already audited by Gate PA
(`research/orrc_v2.py::project_block_w`; relative KKT `<1e-8`, idempotence `<1e-7`, LSMR/LSQR
disagreement `<1e-6`).

## 3. Module 1 — Hodge-Quotient Interaction Transformer

Backbone (mature, capacity-controlled; MolXProt / CAPLA lineage): ligand features (audited 1034-d
Morgan+descriptor tensor, mandatory control; optional GIN only if audited graphs exist) attend to a
light trainable Transformer over frozen ESM-2 segment tokens via shallow PreNorm cross-attention
(existing `model/cfri.py::InteractionHead`, no AttnRes). It emits a per-edge ambient interaction
feature `u(t,d) ∈ R^m` and a raw scalar `h_θ(t,d)`.

**Exact quotient in the forward path.** Stack ambient features over the `N` observed training edges,
`U ∈ R^{N×m}`, and project columnwise `U⊥ = M_X^W U`. Any linear readout `w^T u⊥` then satisfies
`Xᵀ W (U⊥ w) = 0` exactly, so the interaction head **cannot** contain target intercepts, ligand
potency, target/ligand scale, or any main-effect combination. The interaction is supervised against the
equivalently projected label `r_int = M_X^W y` (the exact residual `Substrate.residual`), so the
projection is in the forward path *and* in the (equivalent) profiled objective — not a penalty.

**Mixed-difference supervision.** For an observed 2×2 rectangle `(t_a,t_b)×(d_i,d_j)`,
```
Δy = y(t_a,d_i) − y(t_a,d_j) − y(t_b,d_i) + y(t_b,d_j)
```
algebraically annihilates `a_{t_a},a_{t_b},b(d_i),b(d_j)` and equals `Δg` exactly (proof: substitute
the decomposition; all four main-effect terms cancel in the double difference). Losses: a robust
(Huber) mixed-difference regression `Huber(Δĝ − Δy)`, a rank-reversal sign loss
`softplus(−sign(Δy)·Δĝ)`, and an optional margin scaled only by train-set measurement uncertainty.
Rectangles use only observed cells; **no missing affinity is imputed** (incomplete regions rely on the
observed-edge projection, not completion).

The fixed nuclear-norm ORRC estimator (`research/orrc_v2.py::fit_orrc`) is retained as a **train-only
mathematical teacher and control** — not the primary model, and the OSA monotone shrinker
(`OSA_ORRC_ARCHITECTURE_FAIL_REVIEW`) is not reused.

## 4. Ambient interaction coefficient and its covariance (teacher)

For each **training** target `t`, its empirical interaction direction in ligand-feature space is the
ridge solution over that target's own observed edges (cross-fitted across homology folds):
```
β̂_t = argmin_β Σ_{i∈obs(t)}  ( r_int(t,d_i) − βᵀ v_{d_i} )²  + ρ ‖β‖²
     = (Vᵀ_t V_t + ρI)^{-1} Vᵀ_t r_t ,     V_t = stacked ligand features of t's edges,
V̂arβ̂_t = σ̂²_t (Vᵀ_t V_t + ρI)^{-1} (Vᵀ_t V_t) (Vᵀ_t V_t + ρI)^{-1}   (sandwich).
```
`β̂_t` is a **genuinely noisy, target-specific** estimate from the target's own data — deliberately
**not** `Bᵀu_t` (which, being a deterministic function of the ESM feature, would make the protein→
subspace question trivial). The ORRC bilinear `B` is the low-rank denoised aggregate of `{β̂_t}` across
targets and is used only as a control/teacher for effective rank.

## 5. Module 2 — Error-Corrected Grassmann Bayesian Meta-Adapter

**Protein → Grassmann point.** A light Transformer over frozen ESM-2 tokens maps the target to
`R_t ∈ R^{m×r}` with `R_tᵀ R_t = I_r`, realised by a QR/Householder retraction (StelLA, NeurIPS 2025:
USVᵀ with Stiefel factors; here the point is *predicted* from ESM, not learned per task — an
amortised manifold map in the spirit of Forgione et al., 2025). The biological object is the projector
`P_t = R_t R_tᵀ`. The Bayesian feature is `φ(t,d) = R_tᵀ u(t,d) ∈ R^r`. Because the ranking posterior
(below) lives in `R^r`, every correction `R_t m_t ∈ range(P_t)`: **support cannot rotate the correction
outside the protein-selected subspace.** This is a *hard* confinement, unlike the falsified soft
covariance orientation. Rank `r∈{3,4,6}` is chosen by train-only nested folds; `r=8` is not inherited.

**Grassmann invariance.** `R_t` and `R_t O` (`O∈O(r)`) give the same `P_t`, `range(R_t)`, and posterior
predictive; supervision and metrics therefore use **projector-level** quantities (`P_t`, principal
angles), never raw basis columns.

**Error-corrected containment (Stage-D objective).** Write `β̂_t = β_t + ξ`, `ξ∼N(0,V_t)`. Then
`E[β̂_tᵀ(I−P_t)β̂_t] = β_tᵀ(I−P_t)β_t + tr((I−P_t)V_t)`. The **unbiased** estimator of the true
outside-subspace energy is therefore
```
L_ECG(t) = [ β̂_tᵀ(I−P_t)β̂_t − tr((I−P_t)V_t) ]₊ / max( β̂_tᵀβ̂_t − tr(V_t), ε ) .
```
The subtraction (audited correction — the spec's "+" is treated as a typo, because only the minus makes
`E[·]=β_tᵀ(I−P_t)β_t` and only it yields the "negative noise-corrected energy" the spec says to clamp
conservatively; `[·]₊` is that clamp). A direction is worth including in `range(P_t)` only when its
observed energy exceeds its noise variance — this stops the protein from "recovering" pure noise.

**Error-corrected empirical-Bayes prior.** With `β̂_t | β_t ∼ N(β_t,V_t)` and `β_t ∼ N(0,Ω_0)`,
`Ω_0 = L Lᵀ + diag(exp δ)`, fit `Ω_0` on **training targets only** by marginal likelihood of
`β̂_t ∼ N(0, Ω_0 + V_t)` (never treating `β̂_t` as noise-free, never touching development/confirmation).
For a new target restrict the prior to its subspace: `Σ_{0,t} = R_tᵀ Ω_0 R_t ∈ R^{r×r}`.

## 6. Separate calibration and ranking posteriors

**Ranking posterior (subspace, contrast space).** Helmert `H_k∈R^{(k−1)×k}`, `H_k1=0`,
`H_kH_kᵀ=I_{k−1}`. With support residuals `r_s`, subspace features `Φ_s=[φ(t,d_i)]`, reliability `R`,
noise `σ²`:
```
r̃ = H_k r_s ,   Φ̃ = H_k Φ_s ,
Σ_t = ( Σ_{0,t}^{-1} + Φ̃ᵀ R Φ̃ / σ² )^{-1} ,   m_t = Σ_t Φ̃ᵀ R r̃ / σ² ,
Δ_q = φ(t,d_q)ᵀ m_t  (query ranking correction) .
```
All solves FP32 Cholesky; never an explicit inverse (mirrors validated `model/bayes_posterior.py`).
`H_k 1 = 0` gives residual-offset invariance (adding a constant to `r_s` leaves `r̃` unchanged) and,
with the symmetric Gram, support-permutation invariance. `k=0`: `m_t=0`, exact base fallback.
`k=1`: `H_1` has zero rows ⇒ ranking correction **exactly 0** (calibration may still act).

**Calibration posterior (scalar).** A 1-d Gaussian posterior for the target residual intercept `c_t`
from the support residual mean. It shifts absolute affinity only and **cannot** change within-target
order. It is reported and ablated separately and never feeds the ranking Bayes factor.

## 7. Evidence-derived adaptation (no neural gate)

Two separate Bayes factors from Gaussian marginal likelihoods:
- ranking: null (`r̃` is noise) vs alternative (`r̃` from the subspace-restricted posterior) →
  `π_rank = σ(log BF_rank)`;
- calibration: null vs intercept model → `π_cal`.

Reported predictor:
```
ŷ(t,d) = b(d) + π_cal(t)·c_t + π_rank(t)·Δ_t(d) ,   π_cal(0)=π_rank(0)=0 .
```

## 8. Staged training (fail-fast; each stage has a mechanism gate)

* **A base** — ligand-only `b(d)`, heteroscedastic; provenance-audited out-of-component cross-fitted
  predictions; frozen thereafter.
* **B teacher** — fixed nuclear-norm ORRC on train components; per-target `β̂_t`, `V̂arβ̂_t`; effective
  rank, bootstrap principal-angle stability, mask/document sensitivity.
* **C quotient representation** — cross-attention interaction under exact `M_X^W`; projected-label +
  mixed-difference + rank-reversal objectives; **gate: true pairing beats target and ligand
  derangement** before proceeding (PD-M already shows 3/8 feature-explainable directions).
* **D protein→Grassmann** — train `R_t=f(ESM_t)` by error-corrected containment; **gate: held-component
  subspace recovery with true ESM beats shuffled, random and pooled-only protein** (this report's
  decisive train-only experiment; see `hqgbma_preregistration.md`).
* **E Bayesian episodic meta-training** — one target = one task; matched nested support/query; support
  labels enter only the exact posterior; shared neural weights trained through query loss; **no
  test-time neural update**.
* **F optional joint refinement** — low-LR, only if C–E pass, preserving the exact projection and every
  Bayesian invariant.

## 9. Identifiability analysis and failure cases

* **Identified**: `g` in the quotient space (`XᵀWg=0`), `β_t` up to `V_t` noise (error-corrected),
  `range(P_t)` on the Grassmann (invariant to `O(r)`), the two posteriors (exact, permutation- and
  offset-invariant), `k=0/1` fallbacks.
* **Not identified / out of scope**: a signed zero-shot prior mean (locked route; would break exact
  `k=0`); target identity features; anything reading development/confirmation/sealed labels.
* **Failure cases the design will honestly surface**: (i) if the panel's interaction has effective rank
  ≫ 6, a small `r` cannot contain `β_t` and Stage D fails on nested selection; (ii) if `β_t` is not a
  generalisable function of ESM across homology components, true ESM will **not** beat shuffled — the
  central hypothesis is then falsified (a real, expected possibility given the program's null history);
  (iii) if `V_t` is mis-scaled, the error correction can over- or under-subtract — audited via the
  noise-corrected-energy sign and a `V_t`-scale sensitivity check.

## 10. Admissibility of predictive claims

Per the OPEN-S audit (`NO_OPEN_POWERED_INDEPENDENT_PANEL`) and the spent status of the Metz development
rows, **no adequately powered independent open panel exists**. Therefore the decisive experiment in this
session is the **train-only** Stage-C/D mechanism test. Any predictive score on spent development rows
is reported as `ARCHITECTURE_MECHANISM_RESULT_ONLY__NO_INDEPENDENT_PREDICTIVE_CLAIM`, never as external
validation. All firewalls (target id, accession, homology component, parent connectivity, scaffold,
Tanimoto<0.95, document, assay) are retained.

## Literature audited (inspiration, not adoption)

AdaMBind (Nat Commun 2026) — reject full MAML core, keep task-adaptivity motivation; MolXProt
(JCTC 2026) — mature cross-attention backbone; Neural-ANOVA (MLSP 2025) & ANOVA-TPNN (ICML 2025) —
interaction identifiability requires an explicit orthogonality/sum-to-zero constraint (→ weighted Hodge
quotient); StelLA (NeurIPS 2025) — Stiefel/QR subspace factor (→ predicted Grassmann point); Manifold
Meta-Learning (2025) — amortised low-dim manifold map, no second-order/test-time updates; Similarity-
Aware Evaluation (ICLR 2025 Oral) — random splits inflate; similarity-aware firewalls mandatory.
