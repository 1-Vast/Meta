# A2S-TRACE Exploration Prompt

Date: 2026-08-01
Status: exploration prompt — defines the objective, the fixed contract, the admissibility bar, and the
preregistration surface for the next meta-adaptation mechanism in the A2S-DTA programme.
Supersedes: `A2S_EXPLORATION_PROMPT.md` (2026-08-01) for the substrate, the power budget and §2 constraints.
It does **not** relax any admissibility condition in that document.

Use: self-contained. Read together with
`reports/active/A2S_IDA_TAMSK_ADRO_DEEP_RESEARCH_SYNTHESIS_2026-08-01.md`,
`reports/active/A2S_FINAL_PI_META_MECHANISM_REDESIGN_2026-08-01.md`, and
`reports/active/A2S_TRACE_MECHANISM_ANALYSIS_2026-08-01.md`.

---

## 0. What is fixed and what is open

**FIXED — do not renegotiate.**

1. **The task.** Abundant-to-scarce drug–target affinity transfer. Meta-train on abundant source
   targets. At meta-test a strictly unseen recipient target supplies `k ∈ {1,3,5}` measured support
   affinities drawn under a **declared support policy**. The model ranks that target's query compounds.
2. **The paradigm.** Meta-learning. The adaptation behaviour is learned across source episodes. No
   recipient-specific SGD, no recipient-specific analytic refit of a free parameter that was not
   meta-learned, no recipient label read at any point outside the k support labels.
3. **The primary endpoint.** Target-macro **compound ranking**: CI and NDCG@10, with pairwise proper
   log loss. RMSE/MAE are secondary and may never be substituted for the primary endpoint.
4. **The label firewall.** Locked-role and recipient labels are sealed. No sealed label may inform
   diagnosis, architecture, hyperparameters, stopping or mechanism admission.
5. **The statistical unit.** The protein-homology component. Aggregation order is fixed:
   episode draws → seed mean within target → target mean within component → paired component bootstrap.
6. **The support policy is part of the estimand.** Any result must name the policy under which support
   was drawn. A gain under one policy is not evidence for another.

**OPEN — actively encouraged to change.**

- The support-free base: architecture, capacity, schedule, cross-fitting, stopping rule.
- The ligand and protein representation, including *pairwise relation* features between a support
  compound and a query compound.
- The adapter's functional form, its state, and its output parameterisation.
- The objective, the counterfactual construction and the loss weighting.
- Episode construction **within** the sealed data contract, including the support policy, provided the
  policy is declared before the run and reported with every number.
- The metric family reported *alongside* the primary endpoint.

---

## 1. Objective

Produce **one mechanism with one load-bearing claim** that makes correctly assigned support labels
produce an **absolute** target-macro ranking improvement over the identical frozen support-free base on
unseen source components, that **beats the strongest analytic support baseline** at equal information,
and that fails cleanly when the claimed innovation is ablated.

A mechanism is admissible only if all seven hold:

1. **Learned.** The adaptation rule is meta-learned from source episodes.
2. **Identifiable at the deployment budget.** The number of target-specific quantities the mechanism
   must estimate at meta-test is stated explicitly and justified against `k ≤ 5` by a reliability
   argument, not by hope. Stating "zero target-specific free parameters" is an acceptable and
   preferred answer, but then the burden moves to showing the amortised object transfers.
3. **Query-dependent.** Its action is not an episode-constant. An episode-constant is rank-null and
   worth exactly zero on the primary endpoint.
4. **Structurally abstaining.** Null, deranged or magnitude-matched wrong evidence produces an
   *equivalent no-op* by construction, not by a trained tendency. `r_S ≡ 0 ⇒ Δ ≡ 0` must hold to
   floating-point exactness, provable from the functional form, not measured.
5. **Bounded.** The correction magnitude is bounded by an observed quantity. No unbounded learned
   scale multiplying a residual aggregate (this reproduced the CMAL failure, constraint C2).
6. **Nested-falsifiable.** The claimed innovation is a single component whose removal recovers a
   *named existing baseline exactly*, so the claim is exactly the measured delta.
7. **Not shortcut-driven.** It survives magnitude-matched, chemistry-preserving controls, and the
   controls are structural where possible rather than empirical.

---

## 2. Binding measured constraints

Established results in this programme. Treat as facts; do not rediscover them; do not propose anything
that contradicts them without first falsifying them.

**C1 — The well-identified channel is rank-null.** CI, Spearman and NDCG are invariant to an
episode-constant shift. The target level `A(t)` is worth exactly zero on the primary endpoint. The
programme's replicated positive result (shrunk anchor, RMSE 1.3511/1.3243/1.3132 at k=1/3/5 vs
`f0_only` 1.3724) lives entirely in that channel. **Any proposal whose gain can be produced by an
episode-constant is not a ranking mechanism.**

**C2 — The unshrunk anchor is actively harmful.** `f0_anchor` scores RMSE 1.6159/1.4680/1.4394.
Any operator that multiplies a raw support-residual aggregate by an unbounded learned scale reproduces
the A2S-CMAL failure.

**C3 — Reliability at the deployment budget.** For a target-specific code of dimension `m` estimated
from `k` labels, `ρ_k = τ²/(τ² + σ²/k)`. Measured post-hoc: `τ_z ≈ 0.185`, `σ ≈ 0.997`, so
`ρ₅ ≈ 0.147` and the identifiability certificate fired on **0.000** of episodes. `σ` and `τ` are
**estimator properties, not data properties**, but under the frozen-feature rank-`m` code class the
target-specific-code route is measured to be unidentifiable at `k ≤ 5`.
**Corollary that this prompt makes binding: a mechanism that estimates zero target-specific
parameters is not subject to C3 at all. That is a legitimate and preferred design response.**

**C4 — The base has no ordering skill to correct.** Frozen base on source holdout: CI 0.5112,
Spearman 0.0334, within-episode R² −1.432. A residual correction on top of chance-level ordering is not
a small-perturbation problem.

**C5 — There is a strong chemistry shortcut.** A label-free, protein-free classifier identifies the
correct support arm at 51.6 %/54.0 % against 25 % chance. Wrong-target support recovers 61.2 % of the
correct arm's training ranking gain. Consequence made binding here: **the per-pair transport weight
function must not read support labels or residuals.** If weights are label-free, correct and deranged
support receive *identical* weights and any measured difference isolates assignment exactly.

**C6 — At k=5 a support-local smoother saturates the label budget.** A2S-MDK realised effective dof
0.99/2.89/4.74 at k=1/3/5; at k=5 it is indistinguishable from a pooled-ridge fine-tune. Any advantage
at k=3/5 must come from a **prior over the correction's shape**, learned across source targets.

**C7 — Power, restated for the new substrate.** Component-level SD of the paired ranking difference is
≈0.005–0.012 CI. At 80 % power: ≈12 components detect 0.010, ≈47 detect 0.005, ≈120 detect 0.003.
The previous source splits had 12/15 components. **The chembl37 pKi formal v4 source pool has 267
targets with ≥64 compounds each**; after homology clustering and a 50/25/25 role split the development
(probe) role is expected to hold on the order of 60+ independent components and the fit role 130+.
This raises detectable effect sizes into the 0.003–0.005 range for the first time in the programme.
An observed effect materially larger than theory permits at the given n is a shortcut until proven
otherwise.

**C8 — Assignment controls must be magnitude-matched.** Permuting support **labels** does not preserve
the residual multiset: `Var(r^π) − Var(r) ≈ 2·Cov(y, μ) > 0` whenever the base has any skill. Only
permuting **residuals** preserves every moment and isolates the compound↔evidence assignment. At k=1 no
within-support assignment permutation exists; norm-matched sign flip tests residual *sensitivity* only.

**C9 — Information admission is stratum-dependent, and the stratum has never been resolved.**
The balanced ChEMBL v2 gate (222/110/107 fit/probe/locked components) returned no positive lower bound
for `Δ_label` at k=1/3/5 nor `Δ_assign` at k=3/5, while its synthetic label channel was strongly
detectable (`+0.23` to `+0.39`) and a high-data target oracle showed positive ranking headroom
(LCB `+0.026/+0.070/+0.071`). A separate BindingDB branch, sampling support randomly *within* a target
and without assay/document closure, reported fixed Tanimoto KRR at `+0.0485 [0.0325, 0.0651]` CI.
**INFERENCE:** these are not contradictory; they are different support policies over different
support–query relations. Nearest support–query Tanimoto in the ChEMBL passive construction is ≈0.223
with ≈81 % scaffold-cold queries.
**Binding consequence: no mechanism may be trained before the admissible stratum is measured. A
mechanism trained where the analytic information gate is null is uninterpretable.**

---

## 3. Prohibited claims and language

Do not assert, and reject any proposal that asserts:

- that KRR/GP predictions are confined to the convex hull of support labels;
- that `(σ/τ)² ≈ 29` is a universal sample-complexity bound;
- that effective degrees of freedom near `k` prove information exhaustion;
- that matching the number of experts to `k−1` establishes identifiability (expert count and
  observation rank are different quantities);
- that a derangement-penalty *loss* is a structural no-op *guarantee*;
- that a train-minibatch/holdout-aggregate comparison establishes an overfit ratio;
- that forward-vs-reverse temporal asymmetry proves campaign causality;
- that a failure to detect information with a finite probe is an information-theoretic impossibility;
- that a split repeatedly used in development is a holdout;
- that a mechanism is *protein-conditioned* because a protein embedding appears in its input;
- that correct-support beating wrong-support demonstrates beneficial adaptation;
- that sharing parameters across many pairs removes a k-shot identifiability constraint.

Also prohibited: query document-year gaps or any other future metadata as inference-time features.

---

## 4. Baselines the proposal must beat, and nested restrictions it must recover

**Must beat, at equal data, equal representation and equal support information:** the identical frozen
support-free base; ligand-only DTA; intercept/slope calibration; ridge; **fixed Tanimoto KRR with a
support-only tuned bandwidth and ridge**; static learned kernel mixture; CKA-linear and CKA-NNLS
kernel weighting; the closed-form Bayesian/MDK shrunk-anchor-plus-local-kernel posterior; an
episode-level kernel router of the TAMSK family; MAML; ANIL; CNP/ANP/MetaDTA; MetaFun; ADKF-IFT where
feasible; an equal-capacity generic listwise reranker; and A2S-CMAL frozen as a failed but mechanically
active baseline.

**Must recover exactly as nested restrictions.** For each, state the parameter setting that reduces the
proposal to it. If a baseline cannot be recovered as a restriction, explain why the comparison is fair.

---

## 5. Required output

1. **Mechanism statement in one sentence**, naming the single learned object.
2. **Mathematical definition**, including the exact quantity estimated at meta-test and its dimension.
3. **Identifiability argument** against `k = 1, 3, 5` separately. State what the mechanism abstains from
   at each budget, and whether C3 binds at all.
4. **Why it is not** calibration, interpolation, retrieval, ridge/KRR, a fixed GP/Bayesian posterior,
   fine-tuning, a bigger encoder, a renamed CNP/ANP/set encoder, a renamed MKL/CKA weighting, a generic
   learned functional update, or a generic listwise decoder. Name the closest prior art with a link and
   state the exact increment.
5. **Source meta-training procedure**, including how base residuals are made strictly out-of-fold
   (*every label-trained layer of a fold model excludes that fold* — head-only cross-fitting is not OOF)
   and how representation drift across folds is prevented.
6. **Meta-test procedure**, one pass, no recipient fitting.
7. **Ablation ladder**, one rung per claim, each falling back to a named baseline.
8. **Shortcut and leakage controls**: magnitude-matched residual derangement (k≥3), norm-matched sign
   flip (k=1), residual-null, chemistry/norm/assay-matched wrong support, label-noise dose–response,
   protein shuffle **and** protein zero, target shuffle, query permutation, distractor insertion,
   query-subset and library-size stability. Say which are structural and which are empirical.
9. **Registered predictions**, stated before running, including at least one that would falsify the
   whole idea.
10. **One decisive source-only falsification experiment**, with metrics, component-level uncertainty,
    preregistered MDE and a stop rule.
11. **Maximum scientific risk**, including the outcome in which the honest deliverable is a null.

Label every substantive statement **FACT**, **INFERENCE**, or **HYPOTHESIS**.

Do not propose a broad hyperparameter sweep. Do not force a winner. If nothing survives §1's seven
admissibility conditions, say so — a rigorous, positively-controlled null with a measured upper bound
and a stated required sample size is an acceptable and valuable outcome.

---

## 6. The two questions to answer, in order

> **Q1 (measurement, must run first).** Under which declared support policy and support–query relation
> stratum does *correctly assigned* support contain transferable ranking information at all, as
> measured by a fixed analytic smoother against the identical frozen base and against a
> magnitude-matched residual derangement?

The previous programme answered a coarser question — "does support information exist?" — with a global
null, then observed a large positive under a different policy on a different corpus. Q1 resolves that
conflict on one corpus, one base, one estimator and one bootstrap, and it defines which task the next
model would actually solve. Nothing may be trained before Q1 returns.

> **Q2 (mechanism).** Given a stratum in which information is admitted, should the mechanism try harder
> to estimate the target's residual field, or should it change *what must be estimated* — replacing the
> target-specific code with an amortised, query-dependent rule for deciding how much each measured
> residual transports to each query compound, with structural abstention when none does?

Every mechanism tried so far — anchor, global code, FiLM/hypernetwork, kernel-ridge posterior,
attention operator, kernel router — either took the base's residual geometry as given and attacked the
*estimation* problem, or attacked the *kernel-selection* problem at episode level. None asked whether
the target-specific estimand could be removed entirely, leaving only an amortised transport rule whose
identifiability burden at `k ≤ 5` is zero. A proposal that ignores Q2 must say why.
