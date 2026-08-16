# A2S-TRACE — mechanism analysis and design

**From the programme record and the Q1 measurement to one load-bearing mechanism**

Date: 2026-08-01
Scope: source-only. `locked` role and A2S recipient labels sealed and never requested.
Companions: `A2S_TRACE_MECHANISM_EXPLORATION_PROMPT_2026-08-01.md` (objective + preregistration),
`a2s_trace_q1_stratum_2026-08-01.json` (Q1), `research/a2s_trace.py` (implementation).

---

# Part A — What the record actually constrains

## A.1 The five numbers that survive audit

| # | Measured | Consequence for any mechanism |
|---|---|---|
| 1 | Shrunk anchor beats `f0_only` on RMSE at every k; unshrunk anchor is harmful (C1/C2) | The level channel is real, identifiable and **rank-null**. Keep it, shrink it, report it separately, never call it adaptation |
| 2 | `τ_z ≈ 0.185`, `σ ≈ 0.997` ⇒ `ρ₅ ≈ 0.147`; the identifiability certificate fired on 0.000 of episodes (C3) | Any mechanism that must **estimate a target-specific code** from `k ≤ 5` labels is dead on arrival in this function class |
| 3 | Frozen base: probe CI ≈ 0.51, Spearman ≈ 0.03 (C4) | There is no good ordering to perturb. A correction is not a small perturbation |
| 4 | Label-free, protein-free classifier identifies the correct support arm at 51.6–54.0 % vs 25 % chance (C5) | An adapter whose weights read labels or residuals cannot be cleanly controlled. Make the weights **label-free** and the control becomes structural |
| 5 | Effective dof 0.99/2.89/4.74 at k=1/3/5; at k=5 a support-local smoother ties a pooled-ridge fine-tune (C6) | Extra gain at k=3/5 must come from a **prior over the correction's shape**, learned across targets — not from extracting more from the k labels |

## A.2 The measurement that unblocked the programme (Q1, this session)

**FACT.** `research/a2s_trace_stratum.py`, probe role, 110 homology components, 12,246 episodes,
three declared support policies, five support→query relation strata, one frozen
component-cross-fitted base, one fixed Tanimoto KRR (λ = 0.1), one true residual derangement, one
paired component bootstrap.

Fixed Tanimoto KRR minus frozen base, target-macro CI, paired component 95 % LCB:

| policy | k | `t<0.20` | `0.20–0.35` | `0.35–0.55` | `t≥0.55` | all |
|---|---:|---:|---:|---:|---:|---:|
| `random_within_target` | 3 | −0.003 | −0.003 | −0.003 | **+0.036** | +0.023 |
| `random_within_target` | 5 | −0.007 | −0.006 | +0.006 | **+0.048** | +0.033 |
| `scaffold_disjoint` | 3 | −0.003 | −0.003 | −0.000 | **+0.023** | +0.015 |
| `scaffold_disjoint` | 5 | −0.000 | +0.001 | +0.012 | **+0.031** | +0.026 |
| `provenance_disjoint` (v2 policy) | 3 | −0.012 | −0.006 | +0.004 | −0.020 | +0.001 |
| `provenance_disjoint` (v2 policy) | 5 | −0.006 | −0.013 | −0.006 | −0.004 | +0.009 |

**FACT.** In every admitted cell the magnitude-matched controls agree: the correct-minus-deranged CI
LCB is positive (+0.021 to +0.062), and the correct-minus-**norm-matched-wrong-target** CI LCB is
positive (+0.019 to +0.057). In every null cell both controls are also null.

**FACT.** The episode-constant (level) channel scored a CI delta of **exactly 0.0000** with a
degenerate bootstrap in all 45 policy × k × stratum cells. C1 is now a verified structural property
of the harness, not an assumption.

**FACT.** The Nadaraya–Watson smoother recovers only ~50 % of the KRR gain
(`t≥0.55`, k=5, random policy: +0.024 vs +0.048 LCB). The support-Gram inverse is load-bearing.

**INFERENCE — the resolution of C9.** Support-label information is a property of the
**support→query chemical relation**, not of the corpus. The v2 `provenance_disjoint` policy draws
queries at mean nearest-Tanimoto 0.19–0.30, i.e. almost entirely inside the two measured null bins;
its global null was a correct measurement of a stratum that contains no information. The BindingDB
positive was a correct measurement of a stratum that does. Both stand.

**INFERENCE — k=1 is not silent here, and that is informative.** At k=1 the stratum-conditional gains
are null, but the pooled `all` cell is positive (+0.0066 LCB, random policy). A single support
compound produces a *query-dependent* correction `T(q,s)·r_s/(1+λ)`, so ordering can change through
the between-query variation of similarity alone. This is a genuine rank action, unlike the
target-code constructions where k=1 was structurally rank-null.

## A.3 The inference nobody had drawn

Every previous mechanism — raw anchor, shrunk anchor, global code, SVD/random basis,
FiLM/hypernetwork (SCAO), kernel-ridge posterior (MDK/BIR), attention operator (CMAL), episode-level
kernel router (TAMSK), identifiability-shaped basis (IDA) — either

- estimated a **target-specific object** from `k ≤ 5` labels (anchor, code, IDA), which C3 kills; or
- selected an **episode-level scalar or mixture** (TAMSK, ECMK, static MKL), which is one number per
  episode and therefore cannot express *which pairs to trust*.

> **INFERENCE — the reframing.** Q1 shows that a fixed, isotropic chemical similarity already
> transports residuals usefully at short range. What it cannot do is know **when short range lies.**
> Tanimoto 0.7 sometimes means a bioisosteric swap that preserves potency and sometimes means an
> activity cliff. Which one it is depends on *what changed* between the two molecules and on *which
> protein* is asked — and that is exactly the kind of regularity that is shared across targets and
> therefore learnable in the outer loop with **zero target-specific parameters**.

That is the mechanism this phase tests.

---

# Part B — Candidate screen

Scored against the seven admissibility conditions of the mechanism prompt.

| Candidate | Learned | Identifiable at k≤5 | Query-dep. | Structural abstention | Bounded | Nested-falsifiable | Shortcut-proof | Verdict |
|---|---|---|---|---|---|---|---|---|
| **B1. TRACE — amortised label-free per-pair transport reliability** | ✅ | ✅ **zero** target-specific parameters; C3 does not bind | ✅ per (query, support) pair | ✅ `r_S≡0 ⇒ Δ≡0` from the functional form | ✅ `\|Δ_q\| ≤ max_i\|r_i\|` | ✅ log-Tanimoto scorer + no null slot + no whitening = NW smoother exactly | ✅ weights never read labels ⇒ correct and deranged get *identical* weights | **CORE** |
| B2. Null-slot abstention (part of B1) | ✅ | — | ✅ | ✅ by definition | ✅ | ✅ recovers B1 with `a₀≡0` | ✅ converts C5's chemistry signal into a confidence channel, not a correction | **ADOPT as sub-module** |
| B3. Label-free Gram whitening of the residual channel | ➖ one meta-learned scalar `λ` | ✅ | ✅ | ✅ (linear in `r`) | ✅ | ✅ `λ→∞` recovers the unwhitened channel | ✅ `K_SS` is label-free | **ADOPT as sub-module** |
| B4. IDA (identifiability-shaped target code) | ✅ | ❌ C3 binds; eight unrepaired defects in the audit | ✅ | ⚠ loss-level only | ⚠ | ✅ | ⚠ | **DEFER** — it estimates the thing Q1 says need not be estimated |
| B5. TAMSK-style episode kernel router | ✅ | ✅ | ❌ episode-level mixture, not pair-level | ❌ convex mixture has no exact zero | ✅ | ⚠ | ❌ router reads residual evidence | **BASELINE, not mechanism** |
| B6. ADRO (active measurement selection) | ✅ | ✅ | ✅ | — | — | — | — | **REJECT for this contract** — changes the observation process; separate track |

**INFERENCE.** B1 is the only candidate that attacks what Q1 actually measured (the relation), carries
no `k`-shot identifiability burden, and makes both of the programme's historical failure modes
(C2 unbounded scale, C5 chemistry shortcut) structurally impossible rather than empirically
monitored. B2 and B3 are not independent mechanisms; they are the abstention state and the residual
conditioning that B1 needs. One model, one headline claim.

---

# Part C — TRACE

## C.1 Mechanism statement

> **TRACE learns a single amortised, protein-conditioned, label-free reliability function that says
> how much *more or less* than isotropic chemical similarity a measured residual on one support
> compound should transport to a specific query compound, together with a per-query gate that lets
> the transport abstain.**

The learned object is one function `g_ψ : (x_q, x_i, p) → ℝ` (the per-pair reliability) plus one gate
`h_ψ : (x_q, p, k, ·) → ℝ`. **Number of target-specific quantities estimated at meta-test: zero.**

## C.2 Definition (as implemented in `research/a2s_trace.py`)

Let `f₀` be the frozen support-free base, `S = {(x_i, y_i)}_{i=1..k}` the support,
`r_i = y_i − f₀(p, x_i)` the frozen-base residuals, `K_SS` the support Tanimoto Gram, `k_qS` the
Tanimoto kernel row of query `q`, and `λ = softplus(θ_λ)` one meta-learned scalar.

```
label-free residual conditioning   r̃    = (K_SS + λI)^{-1} r_S        (linear in r; K_SS is label-free)

per-pair reliability               m_qi = 2·σ( g_ψ( φ(x_q, x_i), c_p ) )   ∈ (0, 2), zero-init ⇒ m ≡ 1
per-query gate                     α_q  = σ( h_ψ( x_q, c_p, k, T_q· ) + b ) ∈ (0, 1), b ≫ 0 ⇒ α ≈ 1
global transport scale             c    = softplus(θ_c) > 0               one target-independent scalar

observed bound                     B_q  = max_i |r_i|
bounded transport                  Δ_q  = clip( c · α_q · Σ_i m_qi · k_qi · r̃_i , ±B_q )

adapted prediction                 ŷ_q  = f₀(p, x_q) + Δ_q       ← RANK channel
                                          + ẑ₀                    ← LEVEL channel, rank-null, RMSE only
```

`φ` is the label-free pair-relation descriptor: Tanimoto, log-Tanimoto, Dice, three normalised
bit-overlap ratios, the ten signed and ten absolute physicochemical descriptor deltas, and (in the
full variant) a learned bilinear interaction `e_q ⊙ e_i`, `|e_q − e_i|` over a learned 1024→64
fingerprint projection. `c_p` is a learned 1280→32 projection of the pooled ESM-2 target embedding.

**Why `clip` and not a smooth squash.** With `m ≡ 1`, `α ≡ 1`, `c = 1`, whitening on, the expression
is *exactly* `k_qS (K_SS + λI)^{-1} r_S` — fixed Tanimoto KRR — as long as the bound is inactive, and
the bound is measured to be inactive at `c = 1`. A `tanh` squash would have made the strongest
analytic baseline unreachable as a restriction. With `weights = k_qS/Σk_qS` and whitening off it is
*exactly* the Nadaraya–Watson smoother. Both restrictions are asserted to floating-point equality in
`tests/test_a2s_trace.py`.

**`c` is a baseline parameter, not a claim.** `c` is an episode constant. It rescales the transport
relative to the frozen base and therefore *does* move a ranking metric — measured on probe, `c`
alone lifts the admitted-stratum CI at k=5 from 0.5730 to 0.5831. Because it is target-independent
and query-independent, it is granted to the analytic bar (rung R2c) and to the static-mixture
baseline, and the mechanism's claim is measured strictly on top of it.

## C.3 Identifiability at k = 1, 3, 5

**FACT.** The meta-test estimand is `Δ_q`, a deterministic function of `(x_q, S, p)` and the frozen
base. No parameter is fitted at meta-test. `ρ_k = τ²/(τ² + σ²/k)` therefore has no argument: there is
no target-specific quantity whose reliability could be computed. **C3 does not bind.**

The burden moves, correctly, to *transfer*: does `ψ` learned on `fit` components hold on unseen
`probe` components? That is measured by the paired component bootstrap, not asserted.

- **k = 1.** One support, one weight, one null slot. The mechanism can only decide *whether* and *how
  strongly* to transport a single residual. It cannot compare evidence. Q1 says the stratum-conditional
  k=1 channel is null, so TRACE is expected to add little here; the k=1 arm is retained as a
  falsification arm for P2, not as a claim.
- **k = 3.** Two centred contrasts exist. The mechanism selects among three supports and may abstain.
- **k = 5.** Four centred contrasts. Full abstain-or-select behaviour with the whitened channel.

## C.4 Why it is not the neighbouring things

| Neighbour | Link | The exact increment |
|---|---|---|
| Fixed Tanimoto KRR / NW smoother | classical | Their weight is a *fixed isotropic* function of one similarity scalar. TRACE's weight is a learned function of *what changed* between the two molecules and of the protein, and it has an explicit abstain state. NW is recovered exactly as a restriction |
| Centred-alignment MKL, CKA-NNLS | <https://jmlr.org/papers/v13/cortes12a.html> | Those choose *episode-level* kernel weights from residual evidence. TRACE's weights are **per pair** and **never read residuals**, which is what makes the derangement control structural rather than empirical |
| TAMSK / MetaVRF / MetaKernel | <https://proceedings.mlr.press/v119/zhen20a.html> | Same objection: one mixture vector per episode cannot express which pair to trust. TAMSK is run here as a baseline family, not as an ancestor |
| ADKF-IFT | <https://arxiv.org/abs/2205.02708> | Meta-learns a representation so that a *task-specific GP solve* performs best — a target-specific estimand at meta-test. TRACE removes the target-specific estimand entirely |
| CNP / ANP / MetaDTA / FS-CAP | <https://proceedings.mlr.press/v80/garnelo18a.html>, <https://openreview.net/forum?id=yzlif16IASM> | Set encoders read `(x, y)` jointly, so `r ≡ 0` is not a structural no-op and correct-vs-deranged support cannot be isolated. Q1's own record shows a plain CNP's "gain" survives residual-null. TRACE's weight path is label-free by construction |
| MetaFun | <https://proceedings.mlr.press/v119/xu20i.html> | A free functional update with no statable abstention or bound. TRACE's update is a convex, bounded, abstaining transport |
| A2S-IDA (internal) | — | IDA estimates a target code; TRACE estimates nothing per target. They are complementary answers to Q2, and Q1 favours TRACE |
| A2S-CMAL (internal) | — | CMAL multiplied a residual aggregate by an unbounded learned scale and its attention read residuals. Both failure modes are structurally excluded here |

## C.5 Structural vs empirical controls

**Structural — provable from the functional form, asserted in tests, not measured:**

| Control | Why it holds | Test |
|---|---|---|
| Residual-null exact no-op | `r_S = 0 ⇒ r̃ = 0` and `B_q = 0 ⇒ Δ_q := 0`. Bitwise zero | `test_residual_null_is_a_bitwise_no_op` |
| Derangement uses identical weights | `m` and `α` depend only on `(x_q, x_i, p)`. Permuting `r` cannot change one weight, so correct−deranged isolates assignment exactly | `test_weights_do_not_read_the_residual` |
| Bounded correction | `\|Δ_q\| ≤ max_i \|r_i\|`, an observed quantity. C2's unbounded-scale failure is impossible | `test_transport_is_bounded_by_an_observed_quantity` |
| Linear in the evidence | Doubling `r` doubles `Δ` below the bound — the operator has no hidden magnitude channel | `test_transport_is_linear_in_the_residual` |
| Support permutation equivariance | no ordering input anywhere | `test_support_order_is_irrelevant` |
| Query permutation / subset / library-size invariance | `Δ_q` depends on `q` and `S` only; no candidate-set statistic enters | `test_queries_do_not_influence_each_other` |
| Exact nesting | R2 reproduces NW and R2b reproduces KRR to floating point | `test_*_restriction_reproduces_the_analytic_smoother` |
| Protein-zero really removes protein | output invariant to a protein shuffle when `c_p ≡ 0` | `test_protein_zero_removes_the_protein_channel` |
| Level/rank separation | `ẑ₀` is an episode constant on a separately reported channel | `test_level_channel_is_rank_null` |
| No future metadata | `φ` reads fingerprints and descriptors only | by construction |

**Empirical — must be measured:** norm-matched wrong-target support, label-noise dose–response,
protein shuffle, protein zero, target shuffle, distractor insertion, query-subset stability, and the
null-stratum non-inferiority check.

## C.6 Ablation ladder

| Rung | Configuration | Falls back to | Claim tested |
|---|---|---|---|
| R0 | frozen base | — | reference |
| R1 | level channel only | shrunk anchor | *sanity*: RMSE gain, **exactly zero** ranking change |
| R2 | `weights=nw`, no whitening, no learned part | NW smoother, exactly | *sanity*: reproduces the Q1 NW number |
| R2b | `weights=krr`, whitening, no learned part | fixed Tanimoto KRR, exactly | *sanity*: reproduces the Q1 KRR number |
| **R2c** | **+ global transport scale `c`** | R2b | the **bar to beat** — one target-independent scalar |
| R3 | + per-query gate `α_q` | R2c | learned abstention |
| **R4** | **+ per-pair reliability `m_qi`** | R3 | **HEADLINE — transport reliability is learnable** |
| R5 | `m` without `α` | R2c | is the gate or the modulation load-bearing? |
| R6 | R4 with `c_p ≡ 0` | R4 | the protein claim, honestly sized |
| R7 | R4 with interpretable pair scalars only | R2c | is a null at R4 capacity or signal? |

Rungs are nested in **optimisation** as well as in parameterisation: a three-epoch warm-up trains
only the analytic scalars, so every learned rung starts from the converged restriction below it and
the measured delta is the learned increment, not an optimisation artefact.

If R4 ≤ R2c, TRACE is not a mechanism and the deliverable is the null plus Q1 — but only if the
**synthetic positive control** shows the identical learner recovering an injected pair-reliability
signal on unseen probe components. Without that, the null measures power, not biology.

## C.7 Maximum scientific risk

**Stated plainly.** Q1 shows that at nearest-Tanimoto ≥ 0.55 a zero-parameter isotropic kernel
already captures a +0.03 to +0.05 CI gain. The residual headroom for *learned reliability* on top of
that may be small — plausibly of the same order as the 0.005 MDE at 66–84 components. The honest
failure mode is therefore not "TRACE breaks" but "TRACE ties KRR", in which case the correct output
is: Q1's stratum map, a measured upper bound on learned transport reliability in this construction,
and the component count a confirmatory test would need. That is a defensible deliverable and it is
preregistered as acceptable.

The second risk is that the admitted stratum is *chemically local by construction*, so a positive
result is a claim about **local SAR transport under a declared support policy**, not about global
dual-cold DTA. Every reported number carries its policy and stratum for that reason.
