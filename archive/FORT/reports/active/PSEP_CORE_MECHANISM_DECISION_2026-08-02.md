# PSEP — core-mechanism decision document

Date: 2026-08-02 · Branch `research/a2s-transfer-object-20260802`
Substrate `dataset/processed/psep.v1` · Runners `research/psep_{substrate,d0,m0,m4,m2,m3}.py`
Artifacts `reports/active/psep_{d0,m0,m4,m2,m3}_2026-08-02.json` (+ record parquets)
Roles read: **`discover` only**. `validate` and `confirm` untouched.

**Headline: a core mechanism target is identified, and one of its two premises is
already falsified.** The A2S programme's terminal negative is **partly
overturned** — a target-specific adaptation object does exist under simultaneous
provenance separation — but the object is high-dimensional, and the fraction of
it reachable at k ≤ 5 is small and close to already-saturated. The honest verdict
is stated in §9 and §16.

---

## 1. Fundamental mathematical diagnosis

### 1.1 The estimand, restated

For target `t` write the affinity of ligand `x` as

    y_t(x) = f_0(x) + h_t(x) + b_{c(t,x)} + eps

where `f_0` is target-agnostic chemistry, `h_t` the target-specific ligand
response, `b_c` a per-measurement-context (document/assay) term, and `c(t,x)` the
context a row was measured in. The A2S programme measured (N0) that **`b_c`
explains 68 % of residual variance with SD ≈ 1.7 pKi**, and that its group is
affine, not offset-only (log-scale SD 0.40).

Three facts follow immediately and they organise everything below.

1. **Any estimator consuming absolute residuals is fitting `b_c` two-thirds of
   the time.** The signal `h_t` is roughly half the size of the nuisance.
2. **Conventional splits do not separate `b_c` from `h_t`.** Under
   scaffold-only splitting this substrate leaves **92.0 % document overlap and
   87.4 % assay overlap** between support and query.
3. **Therefore the only admissible estimand is a within-context contrast.**
   Within-document pair concordance is exactly invariant to any additive `b_c`,
   which is why it is the headline metric throughout.

### 1.2 Identifiability of `h_t` from k observations

Take the working form the brief proposes but does not mandate,
`h_t(x) = B(x)^T z_t`, with `B: X -> R^d` a shared feature map and `z_t` the task
state. Given support `S_t = {(x_i,y_i)}_{i=1}^k`, ridge/GP/MAML-inner-loop
adaptation all produce an estimate confined to

    z_hat_t  in  rowspace(B(x_1),...,B(x_k)),     dim <= k.

So the recoverable part of `z_t` is its projection onto a ≤ k-dimensional
subspace **chosen by the support chemistry, not by the task**. Writing the
population design second moment `Sigma = E[B(x)B(x)^T]`, the excess risk of any
such estimator decomposes as

    E ||z_hat - z||^2_Sigma  =  ||(I - P_k) z||^2_Sigma   +   sigma^2 * tr(...)/k
                                \___ bias: unreached ___/    \__ variance __/

The **bias term is the binding one**, and it is governed by how much of `z_t`'s
`Sigma`-mass lies outside a random k-dimensional support span. If the spectrum of
the task-state distribution is flat, that mass is `~(1 - k/d)` — and no amount of
meta-learning of the *estimator* helps, because the obstruction is in the
*geometry of `z_t`*, not the inference rule.

**This yields the programme's sharpest falsifiable question**, and it is the one
the experiments below answer:

> Is the effective dimension of the task-state distribution small enough that
> k ∈ {1,3,5} spans a useful fraction of it — and if not, can a learned `B`
> make it so?

### 1.3 The impossibility conditions, made concrete

- **k = 1 is structurally silent.** After centring, a single support label
  determines an intercept only. Within-document concordance is invariant to
  intercepts, so *any* k=1 method must score exactly the target-agnostic value.
  Measured: exactly `+0.0000` (§6), not assumed.
- **If `Sigma`-effective-rank of the task distribution ≫ k and protein cannot
  predict the missing coordinates, few-shot adaptation cannot recover `h_t`.**
  Prior programme measurement: pooled ESM-2 → head is `-0.019 [-0.073,+0.019]`,
  i.e. protein predicts nothing zero-shot. That leg is already closed and was not
  re-litigated here.
- **If `h_t` is not expressible in the feature space at all, the question is
  moot.** Tested and *rejected* (§6, `same_doc_cv` = +0.100).

---

## 2. Literature-derived mechanism map

Searched: meta-learning/few-shot, Bayesian & amortised inference, neural
processes, low-rank adaptation & task subspaces, hypernetworks, information
bottleneck, optimal design, inverse problems, domain generalisation, and the
bioactivity-benchmark-integrity literature.

**Directly load-bearing prior art**

| Work | What it establishes | Bearing here |
|---|---|---|
| [Assay-Based ML: Rethinking Evaluation in Drug Discovery](https://chemrxiv.org/engage/chemrxiv/article-details/6806b355e561f77ed42b88a3) | predicting a training assay's mean is a surprisingly strong baseline; proposes assay-disjoint splits | **Independently anticipates the document-oracle result.** Our +0.0860 oracle vs +0.0756 chemical head is a sharper, quantified form of the same phenomenon |
| [Clever Hans in Chemistry](https://arxiv.org/html/2512.20924) | chemist/author "style" signals confound activity prediction; author-disjoint splits | Same family of confound; document ≈ author proxy |
| [Systematic Data Leakage in Protein–Ligand Affinity Benchmarks](https://www.biorxiv.org/content/10.64898/2026.06.29.735309v1) | >6 000 leaking assay pairs, persisting to 0.2 sequence identity | Justifies our homology-component firewall at τ = 0.20 |
| [Provable Meta-Learning with Low-Rank Adaptations](https://arxiv.org/pdf/2410.22264) | sample-complexity gains **conditional on** a genuinely low-rank task family | Supplies the theory our measurement tests the premise of |
| [Meta-learning of shared linear representations](https://arxiv.org/pdf/2501.18975) | shared-subspace recovery needs enough tasks *and* samples per task | Our 553 heads × 266 dims is in the regime it describes |
| [MT-net](https://arxiv.org/pdf/1801.05558) | meta-learns a *metric* warping activation space + a subspace for inner-loop descent | **Closest prior art to our M3 whitening idea** |
| [CG-LoRA / K-FAC subspace](https://arxiv.org/html/2602.16456) | whitens with Kronecker curvature before extracting the subspace, on a function-space alignment argument | **Same mathematical move as M3, published 2026** |
| [AdaMBind](https://www.nature.com/articles/s41467-026-70554-5) | current DTA meta-learning SOTA claim, adaptive task scheduling | The baseline any DTA claim must beat — and it is evaluated on splits that do *not* separate document/assay |

**Transfer audit for the one idea we would import (function-space subspace):**
*Original problem* — inner-loop adaptation wastes its budget on parameter
directions that barely change the function. *Structural correspondence* —
parameter → per-target head `w_t`; curvature/Fisher → ligand design second moment
`Sigma`; task → (protein, endpoint) unit. *New information* — none; it is a
change of metric, not of information. *Failure condition* — it buys nothing once
rank exceeds the point where `Sigma`'s conditioning stops mattering. **Measured
to fail exactly there (§6, M3).**

---

## 3–4. Candidate core innovations, with derivations

Five were formulated; two were killed before implementation by existing
programme measurements, three were tested here.

### C1 — Nuisance-Equivariant Adaptation (NEA)
*Hypothesis:* since per-context offsets carry 68 % of residual variance, fitting
the head on within-document **contrasts** rather than absolute residuals recovers
signal the ordinary head loses. Formally, quotient the label by the group
`G0: y -> y + b_c` (or affine `G`) and estimate `h_t` on the quotient.
*Learned:* nothing extra. *Adapted:* the same `z_t`, on invariant statistics.
*Why not trivially MAML/ridge:* those consume absolute residuals.
*Shortcut:* none obvious. *Falsifier:* contrast-fitted head must beat
raw-residual head. **STATUS: FALSIFIED, §6.**

### C2 — Shared low-rank task subspace (the "IDA/mode dictionary" family)
*Hypothesis:* per-target heads lie near a common m-dimensional subspace, so k
labels need only locate a coordinate in `R^m`, m ≲ k.
*Learned:* the basis `U in R^{d x m}`. *Adapted:* `z_t in R^m`.
*Falsifier (registered by the prior programme):* rank-m retention must rise well
above the −6 % previously measured. **STATUS: PREMISE MEASURED AND FOUND WEAK —
retention is positive but only 20.7 % at rank 2 and 29.6 % at rank 16, §6.**

### C3 — Function-space (Σ-whitened) task subspace
*Hypothesis:* C2 failed because the subspace was taken in the wrong metric.
Compress in `<w1,w2>_Sigma = w1^T Sigma w2`, under which head distance equals
expected squared prediction difference.
*Derivation:* minimising `E_t || B(x)^T(w_t - P w_t) ||^2` over rank-m projectors
`P` gives the eigenproblem of `Sigma^{1/2} Cov(w) Sigma^{1/2}`, not `Cov(w)`.
*Falsifier:* whitened rank-16 retention must exceed raw by > 10 points.
**STATUS: FALSIFIED AT THE OPERATIVE RANK — +5.2 points at rank 16 and converging
to zero advantage by rank 32, §6. It *does* win decisively at rank ≤ 8.**

### C4 — Protein-conditioned prior over the task state
Killed pre-implementation: the prior programme measured pooled ESM-2 → head at
`-0.019 [-0.073,+0.019]`. Not re-run.

### C5 — Predictable per-episode transport scale (TAMSK)
Killed pre-implementation: distance-limited by construction, and the Q1/Q2
measurements bound its ceiling. Not re-run.

---

## 5. Substrate and protocol

The prior programme's registered reopening condition was **not a new
architecture** but **≈ 445 homology components** (vs 92), naming Papyrus 05.7.
That corpus was on disk and unused.

`dataset/processed/psep.v1`:

| | |
|---|---:|
| rows | 420 076 |
| units (accession × endpoint) | 1 371 |
| **homology components** | **828** |
| unique structures | 259 144 |
| documents / assays | 17 451 / 40 096 |
| endpoints | pKi 185 625 · pIC50 219 533 · pKd 14 918 |
| role split (components) | discover **508** · validate 147 · confirm 173 |

**Firewall — strictly stronger than the original.** Every accession holding a
`probe` or `locked` role in the A2S source lock (266 accessions) plus the entire
recipient roster is sealed, **and so is every accession in the same homology
component as a sealed one**. Cost: 1 809 → 1 371 units. No sealed labels were
read; only accessions and sequences, to exclude them.

**Homology components** — single-linkage on 4-mer Jaccard at τ = 0.20,
*calibrated against the existing ChEMBL `hcluster` partition* rather than chosen
freely: on 863 shared accessions within-cluster minimum Jaccard is 0.224 and
between-cluster maximum is 0.248. τ = 0.20 therefore never splits an existing
cluster and merges slightly more — **fewer components, less power, the
conservative direction**.

**Base** — ridge on 1024 Morgan bits + 10 descriptors, cross-fitted by homology
component, deliberately **target-agnostic**, so every target- and context-specific
effect lands in the residual.

**Independent unit** — the homology component, throughout. Bootstrap aggregates
draws → unit → component. Never rows, pairs, episodes or seeds.

---

## 6. Quantitative results

### 6.1 D0 — the registered reopening test (`psep_d0`)

Same 553 units, two splits differing only in what they separate. 379 components.

| | scaffold-only | **separated** (scaffold+doc+assay) |
|---|---:|---:|
| document / assay overlap | 92.0 % / 87.4 % | **0.0 % / 0.0 %** |
| **document-mean oracle** (chemistry-free) | **+0.0860** [+0.0760,+0.0958] | **+0.0000** [0,0] |
| own head (d=26), all pairs | +0.0756 [+0.0648,+0.0869] | +0.0102 [−0.0002,+0.0203] |
| own head (d=26), within-document | +0.0415 [+0.0318,+0.0517] | **+0.0029** [−0.0063,+0.0117] |

Two things to read here.

**(a) The harness self-validates.** On the separated split the document oracle is
*structurally* powerless and measures **exactly zero** — a structural check, not a
statistical one.

**(b) The conventional-split pathology replicates on an independent corpus.** A
chemistry-free oracle that knows only a compound's document scores **+0.0860**,
**beating the entire 26-dimensional chemical head (+0.0756)**. This independently
reproduces the A2S T0 finding on Papyrus at 427 components.

**At d = 26, D0 fails again** — and *more* decisively than before: the prior point
estimate of +0.0123 is now excluded by the upper bound (+0.0117).

### 6.2 M0 — the null was a mis-parameterisation (`psep_m0`)

Within-document gain over base, best ridge, 379 components:

| capacity d | 10 | 26 | 74 | **266** |
|---|---:|---:|---:|---:|
| cross-document head | +0.0022 | +0.0046 | +0.0077 [+0.0019] | **+0.0154 [+0.0090,+0.0220]** |

**Monotone in capacity across every ridge and every target variant** — a
systematic effect, not selection noise. **D0's d = 26 was simply too small.**

*Locality.* A head fitted inside the same documents (cross-fitted over compounds)
scores **+0.1001 [+0.0896,+0.1106]**. So the object is richly *expressible*; only
~15 % of it crosses a document boundary.

*C1 falsified.* raw **+0.0154** > centred +0.0115 > studentised +0.0109.
**Removing the 68 % nuisance from the training target does not help — it hurts.**

### 6.3 M4 — is it target-specific? (`psep_m4`)

All arms identical in capacity (d=266), ridge, and evaluation rows; only the rows
the head is *fitted on* differ. 379 components.

| arm | within-document gain vs base |
|---|---:|
| **global** (target-agnostic, same capacity, 200k+ rows) | **−0.0036** [−0.0073,+0.0002] |
| **own** (unit's own training documents) | **+0.0154** [+0.0090,+0.0220] |
| **wrong** (deranged donor, different component) | **−0.0062** [−0.0123,−0.0000] |
| own − global | **+0.0190** [+0.0123,+0.0260] |
| own − wrong | **+0.0215** [+0.0138,+0.0300] |

**This is the programme's first clean positive.** The gain is not base capacity
(the target-agnostic head of identical capacity is *negative*), and it is not a
generic shortcut (a wrong-target head is *actively harmful*).

### 6.4 M2 — the few-shot budget and the rank obstruction (`psep_m2`)

Gain **over the target-agnostic head**, best ridge per k:

| k | 1 | 3 | 5 | 10 | 20 | 100* | 200* |
|---|---:|---:|---:|---:|---:|---:|---:|
| gain | **+0.0000** | +0.0043 [+0.0005] | **+0.0055** [+0.0015] | +0.0092 [+0.0045] | **+0.0132** [+0.0072] | +0.0192 | +0.0334 |

\* fewer components (176, 98) — deep units only, not comparable to k ≤ 20.

k=1 measures exactly zero, confirming the structural prediction of §1.3.
Smallest k clearing the 0.005 MDE on its lower bound: **k = 20**.

**Head spectrum (553 heads in 266 dims): participation ratio 115.4.** Top-1
direction 2.7 % of variance, top-3 6.8 %, top-16 25.8 %. The task family is
genuinely close to isotropic.

Rank-m retention of the +0.0190 full-head gain (basis from other components):

| rank | 1 | 2 | 4 | 8 | 16 | 32 | 64 | 128 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| raw PCA | 18.4 % | 20.7 % | 21.3 % | 25.5 % | 29.6 % | 40.4 % | 50.7 % | 78.6 % |

The prior programme's **−6 % at rank 2 is corrected to +20.7 %** — low rank is not
*harmful*, as previously recorded; it is merely *weak*.

### 6.5 M3 — the metric, and where it stops helping (`psep_m3`)

Retained fraction, subspace estimated with the evaluated component held out:

| rank | 1 | 2 | 4 | 8 | **16** | 32 | 64 |
|---|---:|---:|---:|---:|---:|---:|---:|
| raw | 18.4 % | 20.7 % | 21.3 % | 25.5 % | **29.6 %** | 40.4 % | 50.7 % |
| **whitened (Σ-metric)** | **33.2 %** | **33.1 %** | **34.4 %** | **33.8 %** | **34.8 %** | 40.8 % | 50.1 % |
| precision-weighted | 24.3 % | 24.0 % | 29.0 % | 29.2 % | 31.7 % | 27.3 % | 47.2 % |

**The single most informative number in this document: whitened rank-1 retains
33.2 %, which raw PCA does not reach until rank ~32.** In function-space geometry
the object is *one dominant shared direction plus an essentially isotropic
remainder*. But the advantage **saturates at ~34 % and vanishes by rank 32**, so
C3's registered falsifier (> 10 points at rank 16) **fails at +5.2 points**.

---

## 7. Falsification and control results

| Control | Purpose | Result |
|---|---|---|
| document-mean oracle on separated split | structural harness validity | **exactly 0.0000** — split is what it claims |
| document-mean oracle on scaffold-only split | positive control | **+0.0860 > chemical head** — confound reproduced |
| target-agnostic head, identical capacity | is the gain just base capacity? | **−0.0036** — no |
| wrong-target head (component-deranged) | is the gain a generic shortcut? | **−0.0062** — actively harmful |
| k = 1 | structural no-op prediction | **exactly 0.0000** as predicted |
| contrast-fitted heads (C1) | does nuisance removal recover signal? | **no** — raw beats centred beats studentised |
| leave-component-out subspace | does the rank basis leak? | enforced throughout |
| three endpoints separately (D0) | is it one endpoint's artefact? | all three null at d=26; consistent |

---

## 8. Rejected mechanisms and exact reasons

- **C1 NEA — REJECTED.** Contrast-fitted heads score *below* raw-residual heads
  (+0.0115 / +0.0109 vs +0.0154). The 68 % nuisance is real, but the within-
  document *metric* already neutralises it; removing it from the *training
  target* only discards usable signal. This closes a mechanism the A2S programme
  designed, registered, and never implemented.
- **C2 shared low-rank subspace — REJECTED as a stand-alone mechanism.** Premise
  measured directly: participation ratio 115.4/266, rank-16 retention 29.6 %.
  The task family is not low-rank.
- **C3 function-space subspace — REJECTED against its registered falsifier**
  (+5.2 points at rank 16, threshold 10). Retained as a *component* (§12) because
  its rank ≤ 8 behaviour is strong and it is the correct geometry.
- **C4 protein-conditioned prior — not re-run**, closed by prior measurement.
- **C5 transport scale — not re-run**, closed by prior measurement.

---

## 9. The strongest surviving mechanism — and its honest status

**Surviving object (validated):** a target-specific, query-dependent, document-
transferable chemical adaptation object. `own − global = +0.0190
[+0.0123,+0.0260]`, `own − wrong = +0.0215 [+0.0138,+0.0300]`, at 379 homology
components, under simultaneous scaffold+document+assay separation, with a
structurally self-validating harness. It passes admission-gate conditions
1,2,3,4,6,7,8,9,10 (§F of the brief). Condition 5 (protein conditioning) is not
claimed.

**Surviving mechanism candidate (NOT yet validated):**

> **The k ≤ 5-identifiable part of the target-adaptation object is one
> coordinate along a single shared direction in function-space (Σ-)geometry;
> everything beyond it is isotropic and unreachable at that budget.**

Evidence for: whitened rank-1 retains 33.2 % of the full-head gain — and free
k=5 ridge fitting achieves +0.0055, which is **29 % of the same gain**. Those two
numbers agreeing is the substantive finding: *k = 5 free fitting is already
operating at approximately the one-shared-direction ceiling.*

**What this means, stated plainly.** The bottleneck is **not** a bad inference
rule that a better meta-learner could fix. At k = 5 the estimator is already
extracting essentially all of what a rank-1 shared structure can offer. The
remaining +0.0135 of the object (from +0.0055 to +0.0190) sits in ~115 effectively
isotropic dimensions and is **information-theoretically out of reach at k ≤ 5**
unless a learned representation changes the task-state geometry itself — which
C3 shows the natural metric change does *not* do beyond rank ~8.

---

## 10. Novelty — assessed honestly, at three levels

**Mathematical novelty: LOW for the mechanism, HIGH for the diagnosis.** The
Σ-whitened subspace is the same mathematical move as
[MT-net](https://arxiv.org/pdf/1801.05558) (2018) and
[CG-LoRA/K-FAC subspace](https://arxiv.org/html/2602.16456) (2026). *We should
not claim it as a new mechanism.* What is new is the **identifiability
measurement**: participation ratio, rank-retention curve, and label-budget curve
of a real few-shot adaptation family measured under provenance separation.

**DTA novelty: HIGH.** No DTA work we found separates scaffold, document and
assay simultaneously and then asks whether a target-specific object survives. The
closest works ([Assay-Based ML](https://chemrxiv.org/engage/chemrxiv/article-details/6806b355e561f77ed42b88a3),
[Clever Hans](https://arxiv.org/html/2512.20924)) establish the confound but stop
there; they do not go on to measure what remains.

**Scientific novelty: HIGH, and it is a negative-plus-positive result.** It
resolves a previously unresolved failure mode: few-shot DTA methods are evaluated
on splits where a chemistry-free document oracle beats the chemical model
(+0.0860 vs +0.0756), and the residual real object is high-dimensional, making
k ≤ 5 adaptation structurally near-saturated. **This reframes the adaptation
problem rather than proposing another architecture** — which is what the brief
asked for.

---

## 11. Why this should outperform the existing approach

It should not be claimed that it will, by much. The defensible claim is:

- Against the **current programme baseline**: +0.0190 over a target-agnostic head
  of identical capacity, where the previous best measurement was +0.0029 and
  indistinguishable from zero. The gain comes from **capacity in the right place**
  (d=266 per-target, not global), which is a diagnosis nobody had.
- Against **published DTA meta-learning** (e.g. AdaMBind): the comparison is not
  yet meaningful, because those are evaluated on splits this substrate shows to be
  ~92 % document-contaminated. The first contribution is the evaluation, not a
  win on the old one.

---

## 12–13. Minimal architecture, if it is built

Only if §15's Stage-1 passes. Hierarchy, ≤ 3 modules:

1. **Core** — per-target head at d ≈ 266 with a **Σ-geometry rank-≤8 shared
   component plus a shrunk idiosyncratic remainder** (empirical Bayes: the shared
   direction is the prior mean, the remainder is shrunk by its own reliability).
   Nested baselines: rank-0 (global head) and full-rank (own head), both measured.
2. **Supporting representation** — *only if* Stage-1 shows the learned basis
   moves the retention curve; otherwise omitted. Its specific job: raise rank-≤8
   retention above 34 %.
3. **Supporting uncertainty** — per-unit reliability weighting for the remainder,
   whose specific job is preventing the negative transfer visible in the wrong-
   target arm (−0.0062).

Frozen encoders, GNN/Transformer layers, MLP heads and schedulers are
**infrastructure**, not modules.

---

## 14–15. Ablation and experimental plan

**Stage 1 (decides whether there is a paper about a mechanism at all).**
Meta-learn a representation `B_phi` on `discover` with the explicit objective of
*concentrating the head spectrum*, and re-measure exactly the M2/M3 curves.
- **Falsifier: Σ-rank-8 retention must exceed 50 %** (vs 33.8 % fixed) **and k=5
  gain must exceed +0.0100** (vs +0.0055), on `validate`.
- If it fails, the paper is the measurement paper (§10), not a mechanism paper.

**Stage 2 (only if Stage 1 passes).** Confirm pre-specified on `confirm` (173
components, never touched). Baselines: support-free DTA, intercept/slope
calibration, kNN/Tanimoto retrieval, ridge/KRR with a meta-learned global
transport scale, MAML, ANIL, AdaMBind, deep-kernel GP. Metrics: within-document
CI (primary), target-macro Spearman/Pearson, NDCG, RMSE/MAE, negative-transfer
rate, across k ∈ {1,3,5,10,20}. Multiple seeds and component bootstrap **only
here** — not during discovery.

**Mandatory controls in every table:** document-mean oracle, wrong-target head,
target-agnostic head of matched capacity, k=1 structural zero.

---

## 16. Risks and explicit STOP conditions

| Risk | STOP condition |
|---|---|
| The +0.0190 object is real but too small to matter | If Stage-1 k=5 gain stays < +0.0100, **stop building models**; publish the measurement |
| Optimistic selection inflated M0/M4 | d=266/ridge=100 were selected on `discover`. **If they do not reproduce on `validate` within its CI, the mechanism claim is void** |
| Capacity effect is a kernel-smoother artefact | If a Tanimoto-kNN baseline at matched k reproduces +0.0154, the "head" framing is wrong |
| Σ-whitening is prior art | Already established (§10). **Do not claim it as novel** |
| pIC50 dominates the corpus | pIC50 is 314/379 components. Any headline must be shown to survive on pKi alone |
| Chasing the same-document +0.1001 | It is inflated by congeneric-series analogues within a document. **Not a target** |

**Not permitted to rescue a failure:** larger models, more modules, more epochs,
more seeds, post-hoc thresholds, unregistered losses.

---

## Verdict

**A CORE OBJECT IS IDENTIFIED AND VALIDATED; A CORE *MECHANISM* IS NOT YET.**

The A2S terminal negative is corrected: a target-specific, document-transferable
adaptation object exists (+0.0190 over target-agnostic, +0.0215 over wrong-target,
379 components, fully provenance-separated). It was missed because it is
high-dimensional and was measured at d = 26.

But the object is close to isotropic (participation ratio 115/266), and free k=5
ridge already reaches ~29 % of it — approximately the one-shared-direction
ceiling. Three candidate mechanisms were derived and tested; **all three were
rejected against their own registered falsifiers**. The single remaining
mechanism hypothesis — that a *learned* representation can concentrate the head
spectrum where the natural metric cannot — is now precisely specified with a
numerical falsifier, and is Stage 1.

No mechanism was invented to satisfy the task.
