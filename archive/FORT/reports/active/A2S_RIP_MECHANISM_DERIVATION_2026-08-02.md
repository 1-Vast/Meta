# A2S-RIP — a transferable target adaptation state as a certified ranking-intervention policy

**Literature-grounded mechanism derivation. Pre-implementation. No mechanism code written.**

Date: 2026-08-02
Scope: source-only. `locked` role and the A2S recipient roster sealed and not requested.
Measured inputs: `A2S_TRACE_Q1_STRATUM_DECISION`, `A2S_TRACE_Q2_MECHANISM_DECISION`,
`A2S_MODE_GATES_A0_A4_DECISION`, `A2S_MODE_GENERALIZATION_DECISION` (all 2026-08-01/02).
Labels: every substantive statement is **FACT**, **INFERENCE** or **HYPOTHESIS**.

---

# 1. The constraint set the mechanism must satisfy

These are this programme's own measurements. A proposal that contradicts one of them is rejected.

| # | Measured | Consequence |
|---|---|---|
| M1 | Support transport `Σ_i w(x_q,x_i,p)·r_i` is admitted only at nearest-Tanimoto ≥ 0.55; null below 0.35 in every policy | the whole transport family is distance-limited |
| M2 | Learned per-pair transport reliability adds −0.0001 CI [−0.0006, +0.0005] over globally-scaled KRR, with a positive control at +0.016–0.026 | similarity reweighting is exhausted |
| M3 | One global transport scale is worth **+0.009 CI**; the frozen base carries ~2× the within-episode spread of the correction while ordering at chance | **correction magnitude is a free scalar and a standing artefact risk** |
| M4 | A per-target head on a 26-dim label-free basis survives a within-target scaffold-disjoint split: **+0.052 CI [+0.029, +0.075]** | a query-only target-conditioned object is real |
| M5 | Source-target heads have a nearly flat spectrum; a rank-2 projection retains **−6 %** of the gain | **no low-dimensional shared structure exists**; every "compact target code" route is closed |
| M6 | Protein → head, zero-shot: −0.019 [−0.073, +0.019] | no protein shortcut; adaptation must come from the support labels |
| M7 | Label learning curve, empirical-Bayes head: k=3 +0.001, **k=5 +0.011 [−0.003, +0.025]**, k=10 +0.026 [+0.012, +0.041], k=40 +0.052 | at k ≤ 5 the head is **estimable but too noisy to apply wholesale**; knee at k ≈ 10 |
| M8 | Per-query support-subset hindsight oracle reaches CI 0.93–0.95 | enormous reordering capacity exists; nothing label-free has predicted which edit is right |
| M9 | The episode-constant level channel is worth exactly 0.0000 CI in 45/45 cells | calibration is not adaptation |

**INFERENCE — the shape of the remaining opening.** M4 and M7 together say the adaptation object
exists and is *partially* estimable at k=5 (point estimate +0.011, interval crossing zero). M5 says
you cannot make it easier by compressing it. So the only remaining move is not to reduce the
*dimension of the object* but to reduce the **number of decisions you commit to**: apply the noisy
estimate only where it is certifiably right, and abstain elsewhere.

That converts an estimation problem below its learning-curve knee into a **selection problem with a
meta-learned threshold** — and selection with a threshold learned across tasks is exactly the one
statistical device the literature has that works at `k ≤ 5`.

---

# 2. Literature: what exists, what does not

## 2.1 Few-shot molecular / DTA adaptation

| Work | Adaptation object | Why it does not meet this contract |
|---|---|---|
| [FS-CAP](https://pubs.acs.org/doi/10.1021/acs.jcim.4c00485) (JCIM 2024) — context compounds + activities aggregated into an assay encoding, separate query encoder, continuous activity output | a continuous context embedding | closest prior art on *ranking* by affinity from few context compounds, but it emits a full correction for **every** query with no abstention, no certification and no reported harm rate. Its own stated weaknesses are a simple representation and a limited-expressiveness aggregation — not the decision layer |
| [MHNfs](https://pubs.acs.org/doi/10.1021/acs.jcim.4c02373) (JCIM 2025) — modern Hopfield in-context retrieval over a large context set | retrieved context | retrieval is a similarity operation ⇒ M1/M2 apply |
| [PACIA](https://www.ijcai.org/proceedings/2024/576) (IJCAI 2024) — parameter-efficient GNN adapter, task-level and query-level modulation | a small set of adaptive parameters | a compact continuous task code ⇒ **M5** closes it on this substrate |
| [AdaMBind](https://www.nature.com/articles/s41467-026-70554-5) (Nat Commun 2026) — meta-learning with easy-to-hard adaptive task scheduling | meta-initialisation + task adaptation | gradient adaptation of many parameters from ≤5 labels ⇒ M7; curriculum changes the training distribution, not the identifiability budget |
| [ADKF-IFT](https://arxiv.org/abs/2205.02708) | representation such that a task-specific GP solve is best | target-specific estimand + similarity kernel ⇒ M1 + M7 |
| [APN](https://academic.oup.com/bib/article/25/5/bbae394/7731658), [CFS-HML](https://pmc.ncbi.nlm.nih.gov/articles/PMC12510055/), Meta-Mol hypernetwork | prototypes / property-specific encoders / task-conditioned weights | all are compact continuous task states ⇒ M5 |

**FACT.** Every one of these emits a dense correction for every query. **None** of them abstains per
query, certifies an intervention, or reports the rate at which its adaptation makes the ranking worse.

## 2.2 The machinery that does work at k ≤ 5

- [**Few-shot Conformal Prediction with Auxiliary Tasks**](https://arxiv.org/abs/2102.08898)
  (Fisch, Schuster, Jaakkola, Barzilay, ICML 2021). Casts conformalisation itself as **meta-learning
  over exchangeable auxiliary tasks**, so a new task with very few labels still gets tight, valid
  sets. Validated partly on computational chemistry. **This is the key import:** the transferable
  object is a *calibration rule*, not a task function — and M5/M6 say a task function is exactly what
  this substrate refuses to share.
- [**Conformal Risk Control**](https://arxiv.org/abs/2208.02814) (Angelopoulos et al.) and
  **Learn Then Test** (Angelopoulos et al., 2021): finite-sample control of any bounded monotone
  risk. The rate of *harmful* interventions is bounded and monotone in a threshold, so it is directly
  controllable.
- [**Selective Conformal Risk Control**](https://arxiv.org/html/2512.12844v2): unifies conformal
  prediction with selective classification, gives simultaneous selective-coverage and
  conditional-risk guarantees — and **explicitly names ranking as future work**. That is the door.
- **Conformal selection in drug discovery** — [Optimized Conformal
  Selection](https://arxiv.org/pdf/2411.17983) (Bai & Jin, 2024), SCoRE
  ([2026](https://arxiv.org/pdf/2603.24704)), [ConfBiXtCPI](https://pubs.acs.org/jcisd8/article-abstract/66/6/3013/5080337/Trustworthy-Compound-Protein-Interaction),
  and [ML-guided docking screens with conformal
  error control](https://www.nature.com/articles/s43588-025-00777-x). All select **which compounds to
  test**, controlling false leads. **A2S-RIP selects which *interventions on an existing ranking* to
  commit to — a different estimand, and one nobody has taken.**
- **SelectiveNet / risk–coverage** (Geifman & El-Yaniv, ICML 2019): abstention trained with an
  explicit coverage constraint.
- A calibration caution worth importing honestly: cross-task transfer of calibration
  ([MARGIN, 2026](https://arxiv.org/pdf/2605.22949)) reliably beats an uncalibrated baseline but is
  reported as several times weaker than adaptation fitted on the deployment distribution itself.
  **INFERENCE:** expect a real but modest effect, which is consistent with this programme's MDE of
  0.005 CI rather than with a large win.

**INFERENCE — the gap, stated exactly.** Meta-learned conformal calibration exists (Fisch). Selective
risk control exists (SCRC) and names ranking as open. Conformal selection in chemistry exists but
selects *compounds*. Few-shot DTA exists but applies dense corrections without abstention. **Nobody
has meta-learned a certified, budgeted policy of ranking interventions for an unseen target.**

---

# 3. A2S-RIP

## 3.1 Mechanism statement

> **A2S-RIP meta-learns, across abundant source targets, the calibration of its own adaptation
> uncertainty — a monotone map from an observable per-compound evidence margin to the probability
> that intervening on that compound improves the ranking — and at meta-test uses `k ≤ 5` recipient
> measurements to place each compound on that curve, committing to a bounded, certified subset of
> ranking interventions and abstaining on the rest.**

The **transferable adaptation state** is the pair `z_t = (ĥ_t, Σ_t)`: the closed-form empirical-Bayes
posterior over the target's response head. The **intervention** is a sparse, bounded score edit. The
**transferable object learned from source targets** is neither `ĥ` nor a compressed code (M5 forbids
both) but the **certification rule**.

## 3.2 Definition

Let `g(x) ∈ R^d` be the compact label-free basis (`d = 26`, already built and measured), `f_0` the
frozen base, `r_i = y_i − f_0(p, x_i)` the `k` support residuals, `Λ = diag(τ_j²)` the source head
prior (measured on `fit` targets, not fitted), `σ` the measured within-target residual scale.

```
Adaptation state (closed form, zero free parameters)
    Σ_t   = σ² ( G_Sᵀ G_S + σ² Λ⁻¹ )⁻¹                       posterior covariance, d × d
    ĥ_t   = Σ_t G_Sᵀ r_S / σ²  +  ( I − Σ_t Λ⁻¹ ) h̄          posterior mean, shrunk to the source mean

Per-compound score and its honest uncertainty
    s_q   = ⟨ ĥ_t , g(x_q) ⟩                                  proposed edit
    v_q   = g(x_q)ᵀ Σ_t g(x_q)                                epistemic variance (from the k labels)
    z_q   = s_q / sqrt( v_q + σ_a² )                          evidence margin  (epistemic + aleatoric)

Certified, bounded, sparse intervention  ← the ranking action
    δ_q   = clip( s_q , ±B ) · 1[ ψ_θ( z_q , e_t , g(x_q) ) ≥ τ_α ]
    ŷ_q   = f_0(p, x_q) + δ_q
```

`ψ_θ` is a meta-learned monotone conformity score (a learned refinement of `z_q`, in the spirit of
[optimised conformal selection](https://arxiv.org/pdf/2411.17983)); `τ_α` is fixed by **meta-conformal
calibration across source targets** (Fisch et al.), *not* by the recipient's 5 labels; `B` is a
declared bound, `e_t` is label-free episode evidence (support spread, Gram conditioning, `k`).

**The decision is per compound.** No candidate-set statistic enters, so query-permutation,
query-subset and library-size invariance are **structural**, not empirical — this programme's standing
requirement, and the reason not to use a top-`m` selection rule.

## 3.3 Why each measured constraint is respected

| constraint | how |
|---|---|
| M1 (transport is distance-limited) | `δ_q` contains no support compound; it is a function of `x_q` and the state. Defined for every query |
| M2 (similarity reweighting exhausted) | no similarity weight anywhere in the action |
| M3 (magnitude is a free scalar) | `B` is declared, and the **magnitude-matched wholesale control** in §3.6 is mandatory |
| M5 (no low-rank shared structure) | the head is used at **full rank**; nothing is compressed. The shared object is the calibration curve, which is one monotone function |
| M6 (protein predicts nothing) | protein may enter `ψ_θ` only as a certification feature, never as the correction |
| M7 (k ≤ 5 is below the knee) | **the entire point**: the noisy head is applied only where certified |
| M9 (level is rank-null) | the level is profiled out and reported on its own channel |

## 3.4 Identifiability at k = 1, 3, 5 — the crux

**FACT.** The target-specific quantities are `ĥ_t` and `Σ_t`, both closed-form functions of the `k`
labels with **no free parameters**. `τ_α` and `ψ_θ` are estimated from hundreds of source targets ×
thousands of episodes.

> **INFERENCE — why this escapes the trap that closed A2S-MODE and A2S-IDA.** Those routes needed the
> `k` labels to *identify a target-specific object*, and `ρ_k = τ²/(τ²+σ²/k) ≈ 0.15` says they cannot.
> A2S-RIP asks the `k` labels only to **order compounds by certainty** and to supply a `Σ_t` that is
> honest about how bad the estimate is. The *cut point* is meta-learned. Ranking a handful of
> candidates by an observable statistic is a far weaker demand than estimating the statistic's target.

- **k = 1**: `Σ_t ≈ Λ`, so `v_q` is near-prior and `z_q` is near zero for almost every compound ⇒ the
  policy abstains almost everywhere. The k=1 rank silence the user chose to keep is preserved and
  becomes a *derived* property rather than a design constraint.
- **k = 3**: a small certified set is expected; coverage should be low and precision high.
- **k = 5**: the regime the mechanism targets — M7 says the wholesale head is +0.011 with the interval
  crossing zero; the claim is that its certified subset is materially better than that.

## 3.5 Nested restrictions (the claim is exactly the measured delta)

| setting | recovers |
|---|---|
| `τ_α = ∞` (nothing certified) | the frozen base, **exactly** |
| `τ_α = −∞`, `B = ∞`, `ψ_θ` inert | the wholesale empirical-Bayes head measured in M7, **exactly** |
| `Λ → 0` | the source mean head (a target-independent control) |
| `r_S ≡ 0` | `ĥ_t = h̄` and `s_q` is an episode-independent function ⇒ no *target-specific* action |

**The headline claim is `RIP − wholesale EB head` at k = 3 and 5.**

## 3.6 Controls, and the one this programme uniquely needs

Structural (provable): per-compound decision ⇒ query permutation / subset / library-size invariance;
`τ_α = ∞` ⇒ exact no-op; bounded by `B`.

Empirical, all mandatory:

1. **Magnitude-matched wholesale control — non-negotiable.** Selective intervention lowers the mean
   `|δ|`. M3 measured that correction magnitude alone is worth +0.009 CI. So RIP must be compared with
   a wholesale head **rescaled to the same mean `|δ|`**. Without this, any "selection" gain is a
   shrinkage artefact. *No paper in §2 runs this control, because none of them measured M3.*
2. **Random-selection-at-matched-coverage.** Intervene on the same *number* of compounds, chosen at
   random. This separates "choosing the right compounds" from "intervening less". This is the
   derangement-equivalent for a selection mechanism and it is the primary falsifier.
3. Residual derangement and norm-matched wrong-target support.
4. Protein shuffle / protein zero inside `ψ_θ`.
5. Validity: realised harmful-intervention rate on probe must be ≤ the certified `α`.
6. Convergence: as `k → 20/40`, RIP must approach the wholesale head (nothing left to abstain from).
   A mechanism that keeps winning at k = 40 is doing shrinkage, not selection.

## 3.7 Registered predictions

| # | Prediction | Falsifies |
|---|---|---|
| P1 | `RIP − wholesale EB head` LCB > 0.005 CI at k = 3 and k = 5 | the mechanism |
| P2 | RIP ≈ wholesale head at k = 20–40 | the selection story (if RIP still wins, it is shrinkage) |
| P3 | Realised harmful-intervention rate ≤ certified `α` on unseen components | the meta-conformal transfer (Fisch's premise on this substrate) |
| P4 | Random selection at matched coverage destroys the gain | the whole idea |
| P5 | Magnitude-matched wholesale control destroys the gain ⇒ **retract**; the effect was M3's scalar | the whole idea |
| P6 | Risk–coverage is monotone; conditional harm falls as coverage falls | the certification |
| P7 | Effect size is small — of order 0.005–0.02 CI, not 0.05 | nothing; it is the honest prior from MARGIN's cross-task calibration finding, stated in advance |

## 3.8 Maximum scientific risk

**HYPOTHESIS at risk.** That the margin→correctness relation is *shared across targets* when the head
itself is not (M5) and the protein map is not (M6). It is entirely possible that reliability is as
idiosyncratic as the head, in which case `τ_α` transfers no better than a prototype and P3 fails.
That is a clean, cheap falsification, and it is the first thing Gate R0 measures.

Second risk: the certified subset may be **empty** at k ≤ 5 for a defensible `α`. Then the honest
deliverable is the coverage-at-validity curve — i.e. the measured statement "at k = 5, x % of
compounds can be intervened on with ≤ α harm", which is itself a publishable bound on few-shot
adaptation and directly useful to a screening campaign.

---

# 4. Gate R0 — the decisive measurement, before any implementation

Consistent with the discipline that has now killed two mechanisms cheaply, nothing is trained until
the ceiling is measured. R0 is one pass over the existing probe episodes and reuses code already
written (`a2s_mode_generalization.py` supplies `ĥ_t`; `Σ_t` is one extra line).

**R0a — selection ceiling.** At k = 3/5, compute the wholesale EB head, then apply it only to an
oracle-chosen subset (hindsight). Sweep coverage 0→100 %. If the oracle risk–coverage curve peaks at
or below the wholesale head, there is nothing selection can buy and A2S-RIP dies here.

**R0b — is the margin informative at all?** Measure AUC of `|z_q|` for predicting whether the edit on
`q` is in the right direction, on **unseen** probe components. Chance = 0.5. This is the transferable
object's existence test, and it needs no training.

**R0c — does the threshold transfer?** Fit `τ_α` on `fit` targets; measure the realised harm rate on
`probe` at the same `τ_α`. Validity gap = P3, pre-tested without a model.

**R0d — magnitude confound, measured up front.** Compare the oracle-selected curve against the
magnitude-matched wholesale head at every coverage. If they coincide, the effect is M3 again.

**Stop rule.** R0a peak ≤ wholesale, or R0b AUC ≈ 0.5, ⇒ report the null with its measured ceiling
and stop. Both pass ⇒ implement `ψ_θ` and the meta-conformal calibration, and run §3.6 in full.

---

# 5. Positioning, stated plainly

Architecturally A2S-RIP is a recombination: an empirical-Bayes posterior (classical), a conformity
score (classical), and meta-learned conformalisation ([Fisch et al.
2021](https://arxiv.org/abs/2102.08898)). Its contribution is (i) moving the transferable object from
the task function — which this substrate measurably refuses to share (M5, M6) — to the *calibration
of the estimator's own uncertainty*; (ii) making the adaptation action a **sparse, bounded, certified
ranking intervention** with a per-compound decision rule, which is the open item in
[SCRC](https://arxiv.org/html/2512.12844v2) and is absent from every few-shot DTA method in §2.1; and
(iii) a control set — magnitude-matched wholesale and random-selection-at-matched-coverage — that
this programme's own measurements show is necessary and that the literature does not run.

If R0 fails, the deliverable is the measured selection ceiling at `k ≤ 5` plus the margin-AUC, which
together bound what *any* selective few-shot adaptation can achieve on this substrate. That is a
defensible result and it is preregistered as acceptable.
