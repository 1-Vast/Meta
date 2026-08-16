# A2S-DTA Stage 1 — meta-learning mechanism discovery

Date: 2026-08-02 · Branch: `research/a2s-transfer-object-20260802`
Status: **design only. No implementation. No code is authorized by this document.**
Evidence base: `A2S_TRANSFER_OBJECT_GATE_T0_DECISION_2026-08-02.md` (rev 2),
`A2S_CFES_C0B_DECISION_2026-08-02.md`, and the nine falsified predecessors.

> ## Binding precondition
> Gate T0 (rev 2) measured that a chemistry-free document-mean oracle (+0.0610)
> **outscores a full per-target chemical head** (+0.0519), that within-target
> scaffold-disjoint splits leave 91.1 % of query rows in a support-seen document,
> that the same-document chemical remainder (+0.0290 [+0.0046, +0.0557]) **fails**
> its 0.005 admission bar, and that cross-target head transfer is **negative**
> (−0.0183 [−0.0435, +0.0054]).
>
> Every mechanism below is therefore **conditional on Gate D0**: rebuild the split
> with simultaneous target / scaffold / document / assay separation and re-measure
> the surviving headroom. If D0 returns a headroom whose lower bound is under the
> MDE, **no mechanism should be selected** and the terminal deliverable is the
> negative result. Stage 1 proceeds so that a mechanism is ready if D0 passes; it
> does not assume D0 passes.

---

# 1. The scientific bottleneck

## 1.1 Four candidate diagnoses, separated

| Diagnosis | Test that would show it | Verdict on this substrate |
|---|---|---|
| **Optimization failure** — the learner cannot fit | synthetic positive control fails | **Ruled out.** Every gate's planted-structure control passed: TRACE recovered +0.016–0.026 of an injected signal; PIRS recovered +0.31 at k=5; T0 recovered the generating head at 18.2 % vs 0.9 % chance |
| **Representation failure** — the basis cannot express the object | full-support *oracle* over the representation fails | **Real but not sufficient.** PIRS's learned coordinates failed at the oracle (+0.0046, CI crosses zero). But the fixed 26-dim basis *does* express a full-support object, so representation is not the only barrier |
| **Insufficient information** — k ≤ 5 labels are too few | learning curve knee far above 5 | **Real.** G4's dense curve knees at k ≈ 10; T0's discrete curve is null or negative at every k ≤ 10 |
| **Incorrect adaptation object** — the estimand is not what we think | the object is reproduced by a non-chemical oracle | **Confirmed, and dominant.** A document-mean oracle beats the chemical head |

## 1.2 The answer to the question as posed

The question was: *"the model cannot learn the target state"* or *"the target state
is not identifiable from k ≤ 5"*?

**Neither, as stated. The target state was mis-defined.**

The quantity the programme has been estimating — the per-target residual response
`r_t(x) = y_t(x) − μ(x)` fitted on a target's own rows — is a **mixture of two
objects with different transfer semantics**:

```
r_t(x)  =  c_t(x)          target-specific chemical response  (transferable, wanted)
        +  δ_t,d(x)        per-document / per-assay offset    (nuisance, not transferable)
        +  ε
```

`δ` is recoverable from chemistry because a ChEMBL document reports a congeneric
series, so ligand features act as a document classifier. Measured: `δ` alone is
worth **more** than the entire fitted head. Nine mechanisms then competed to
estimate a quantity that is majority nuisance — and the nuisance is precisely the
component that **cannot** transfer to an unseen target, because a new target's
documents are new.

## 1.3 Why "k ≈ 10" is the observed number

Not a property of chemistry or of biology. It is where a dense estimator's
shrinkage begins to recover a signal whose ranking-relevant part is small relative
to residual dispersion. Three consequences follow, and all three were measured:

1. the knee moves with the *estimand*, not the estimator — dense EB (G4) and
   discrete selection (T0) knee in the same region despite unrelated hypothesis
   spaces;
2. the apparent gain concentrates where support and query share context (Q1's
   Tanimoto ≥ 0.55 stratum is plausibly the same-document stratum);
3. controls that break chemistry but preserve context (magnitude-matched,
   pooled-head, random-selection) repeatedly matched the mechanism.

**Therefore the bottleneck is estimand mis-specification first, information
scarcity second, representation third, optimization not at all.**

---

# 2. Why previous meta-learning approaches failed

| Method | Transferable object learned | Why it fails on strict unseen-target k-shot DTA | Remaining gap |
|---|---|---|---|
| **MAML** | an initialization `θ₀` | the inner loop is SGD in a space of dimension ≫ k−1; at the measured SNR the gradient step is noise. Nothing in the objective constrains the adapted space to be identifiable | no control of *how many* effective parameters adaptation spends |
| **ANIL** | frozen features + head initialization | reduces to linear-head adaptation, i.e. exactly the empirical-Bayes head whose curve knees at k ≈ 10 | same estimand as the failed baseline |
| **MetaDTA** | support-conditioned regression head on a frozen encoder | inherits the mixed estimand; published protocols split by target/compound, not by document or assay | no measurement-context control in the task definition |
| **AdaMBind** | task-adaptive representation for binding | same estimand; adaptation is continuous and high-dimensional | identifiability never quantified against the label budget |
| **Metric / ProtoNet / matching** | an embedding metric or class prototypes | is similarity weighting; measured null below Tanimoto 0.35, and the surviving local effect coincides with shared context | cannot act on chemically distant queries |
| **Neural processes (CNP/ANP)** | amortized map: support set → global latent | amortization changes the estimator, not the information; the latent is continuous and high-dimensional | no mechanism to bound the latent's identifiable dimension |
| **Hypernetworks / HyperPCM** | descriptor → target-specific parameters | requires label-free target conditioning to carry the load; measured null **four** times here (G3 pooled ESM, PIRS segments, CFES-C0B pocket, T0C shortlist) | the conditioning signal does not exist on this substrate |
| **Bayesian adaptation (VERSA, ALPaCA)** | posterior over head parameters | this *is* the optimal estimator for a fixed hypothesis space; T0 shows the space is the problem | optimality inside the wrong space |
| **Kernel / meta-kernel (R2D2, MetaOptNet, meta-KRR)** | a learned kernel or closed-form inner solve | TRACE measured the learned component at −0.0001 [−0.0006, +0.0005] over scaled KRR | learned transport adds nothing beyond a global scale |

**The gap, stated once.** Every method above treats the *estimand* as given and
competes on the estimator, the representation, or the prior. **None designs the
adaptation target so that it is (i) free of measurement-context nuisance and
(ii) identifiable within a stated label budget.** That is the opening.

---

# 3. The missing transferable object

> **What sources can teach, and a recipient can use, is not a target's response
> function. It is how a target responds to *chemical contrast within a fixed
> measurement context*.**

Three reasons this is the right abstraction:

1. **It is nuisance-free by construction.** `δ_t,d` is additive within a context,
   so any within-context contrast `r_i − r_j` cancels it *exactly*. The dominant
   confound is removed by the definition of the estimand, not by a post-hoc
   control.
2. **It matches how the science is actually done.** Medicinal chemistry reasons in
   ΔΔG within a series; SAR is inherently comparative. A target's identity shows up
   as *which contrasts it amplifies or reverses*.
3. **It is cheap in bits.** A contrast is a signed, bounded quantity. A mechanism
   that consumes contrasts and emits a bounded ranking intervention needs to
   identify a direction and a sign, not a 26-dimensional vector.

---

# 4. Eight candidate mechanisms

Notation throughout: source task `T_i = (S_i, Q_i)`; support `S_t = {(x_j, y_j, c_j)}`
with `c_j` the measurement context (document/assay); frozen support-free base
`μ(p_t, x)`; residual `r_j = y_j − μ(p_t, x_j)`; ligand features `g(x)`; protein
`p_t`; adaptation state `z_t = A_θ(S_t)`; prediction
`ŷ_q = f_θ(x_q, p_t, z_t)`.

---

## M1 — Contrast Response Operator (CRO)

**A. Transferable object.** A learned map from *within-context chemical contrasts*
to a bounded ranking intervention: "when a target prefers B over A by this much,
and the query differs from A in this direction, move the query this way."

**B. Formulation.** Contrast set `P_t = {(j,l) : c_j = c_l}`, contrasts
`d_jl = g(x_j) − g(x_l)`, `Δ_jl = r_j − r_l` (offset-free, exactly).

```
z_t = A_θ(S_t) = Σ_{(j,l)∈P_t} α_θ(d_jl, Δ_jl, k) · ψ_θ(d_jl) · sign(Δ_jl)
ŷ_q = μ(p_t, x_q) + b · tanh( ⟨ z_t , φ_θ(x_q) ⟩ / b )
```

`ψ_θ` embeds a contrast direction; `α_θ` is a learned, budget-aware weight on how
much a contrast of that magnitude should count; `φ_θ` maps the query into the same
space; `b` bounds the intervention. `z_t = 0` when `P_t = ∅`, so support removal is
an exact no-op.

**C. Identifiability at k = 1, 3, 5.** k = 5 within one context yields up to 10
contrasts (k = 3 → 3; k = 1 → 0, an honest structural no-op). The estimand is a
*direction plus sign* in `ψ`-space, not a 26-dim vector, and contrast noise is the
within-context dispersion rather than the total residual dispersion.

**D. Why genuinely meta-learning.** `α_θ`, `ψ_θ`, `φ_θ` are learned across sources:
which contrast directions are *reliable evidence* and how strongly to act on them.
A fixed KRR/ridge/LASSO on contrasts has no way to learn that a 3-log contrast
along one direction is informative while the same magnitude along another is noise.
Retrieval cannot express it because the query need not resemble any support.

**E. Failure risks.** Degenerates to similarity weighting if `φ ≈ ψ` and `α`
saturates — testable by the Tanimoto < 0.35 stratum. Could collapse to
ligand-only if `z_t` becomes support-independent — testable by label permutation.
**Nearest negative evidence: the assay-coherence gate already failed** when
restricting data to exact assays; CRO differs by changing the *estimand* to
contrasts rather than restricting the data fed to a dense head, but this is a real
warning.

---

## M2 — Transductive Library-Conditioned Adaptation (TLCA)

**A. Object.** A map from the *unlabeled query library* plus k labels to an
intervention. The library carries the high-dimensional part; labels carry ~1 bit.

**B. Formulation.** `U_t = {x_q}` the unlabeled candidate set.
`h_t = E_θ({g(x) : x ∈ U_t ∪ S_t^x})` (a permutation-invariant set encoding, no
labels); `z_t = A_θ(h_t, {r̃_j})` where `r̃` are centred residuals;
`ŷ_q = μ + b·tanh(⟨z_t, φ_θ(x_q)⟩/b)`.

**C. Identifiability.** The label-dependent part of `z_t` is constrained to ≤ 2
scalars; everything high-dimensional comes from label-free data available at
deployment.

**D. Meta-learning.** Learns what a *library composition* implies about a target's
pharmacology — unavailable to any estimator that sees only supports.

**E. Risks.** Target memorization via library fingerprinting; leakage if the query
library encodes the answer. T0C's null on coarse chemotype overlap is direct
negative evidence, though it tested a fixed statistic, not a learned encoder.

---

## M3 — Adaptation Program Selection (APS)

**A. Object.** A small meta-learned library of adaptation operators plus a
support-conditioned routing posterior.

**B. Formulation.** `Ω = {ω_1..ω_M}`, `M ≤ 4`;
`z_t = softmax(log π − Σ_j ℓ_θ(r̃_j, ω̃_m(x_j)) / T_θ(k))`;
`ŷ_q = μ + Σ_m (z_m − π_m) ω_m(x_q)`.

**C. Identifiability.** Routing among `M ≤ 4` costs ≤ 2 bits.

**D. Meta-learning.** The library is shaped for k-shot discriminability.

**E. Risks.** **Novelty is the problem, not feasibility** — this is a modular /
mixture-of-experts task-mode meta-learner (Modular Meta-Learning, MMAML, CNAPs,
online task mixtures), and information-constrained task inference is already
explicit in the literature. Also degenerates if operators collapse together.

---

## M4 — Constraint / Exclusion Adaptation (CEA)

**A. Object.** A learned map from supports to *inequality constraints* on the
target's ranking, rather than to a predicted response. Adaptation as feasibility,
not estimation.

**B. Formulation.** Each within-context contrast asserts a half-space
`⟨w, d_jl⟩ · sign(Δ_jl) > 0` on the unknown response direction `w`.
`z_t = A_θ({half-spaces}, p_t)` returns the analytic centre of the feasible cone,
with a learned prior cone from sources; `ŷ_q = μ + b·tanh(⟨z_t, φ_θ(x_q)⟩/b)`.

**C. Identifiability.** A half-space is 1 bit and is robust to label magnitude —
only the *sign* of a contrast is used, so measurement scale and offset both drop
out. k = 5 gives up to 10 constraints.

**D. Meta-learning.** Sources teach the prior cone and which constraints deserve
weight — a fixed convex solver has neither.

**E. Risks.** Sign information may be too weak at the measured noise level;
degenerates to a linear classifier on contrasts.

---

## M5 — Identifiability-Shaped Coordinates (ISC)

**A. Object.** A ligand representation trained so that the per-target response is
1–2 dimensional *in that representation*, optimized through the k-shot inner solve.

**B. Formulation.** `φ_θ : x → R^m`, `m ∈ {1,2}`;
`z_t = argmin_a Σ_j (r̃_j − ⟨a, φ̃_θ(x_j)⟩)² + λ‖a‖²` (differentiable);
outer objective is the *k-shot* ranking loss, so `φ_θ` is pushed toward
identifiability rather than variance.

**C. Identifiability.** `m ≤ 2` against k−1 ≥ 2 contrasts at k = 3.

**D. Meta-learning.** The representation is the meta-object; PCA/low-rank optimize
variance, this optimizes few-shot adaptation error.

**E. Risks.** **PIRS is the direct precedent and it failed at the oracle** — a 2-D
learned coordinate could not express the headroom. Reviving this requires
explaining what changes.

---

## M6 — Contextual Nuisance Factorization (CNF)

**A. Object.** A learned factorization separating context offset from chemical
response, with only the chemical factor transferred.

**B. Formulation.** `r = δ_c + c_t(x) + ε`; meta-learn an amortized decomposition
on sources; at deployment infer `δ` per context from supports and apply only `ĉ_t`.

**C. Identifiability.** `δ` is a scalar per context; the chemical part is bounded.

**D. Meta-learning.** Learns what context effects look like across many targets.

**E. Risks.** **Likely forbidden** — inferring and removing an offset is close to
calibration, and an offset cannot change within-context ranking at all. Best value
is as a *diagnostic and control*, not a mechanism. Retained for completeness.

---

## M7 — Ordinal Response-Mode Transition (ORMT)

**A. Object.** A learned finite-state program: the recipient starts in a generic
prior state and each support measurement triggers one learned discrete transition.

**B. Formulation.** `s_0 = s_prior`; for each support `j` in a canonical order,
`s_{j} = τ_θ(s_{j-1}, g(x_j), r̃_j)` over a finite state set `|S| ≤ 8`;
`z_t = emb(s_k)`; `ŷ_q = μ + b·tanh(⟨z_t, φ_θ(x_q)⟩/b)`.

**C. Identifiability.** k transitions over ≤ 8 states cost ≤ 3 bits total, and the
budget is structural — the state set cannot expand.

**D. Meta-learning.** The transition function is the transferable object: a learned
*adaptation program*, not a learned estimator.

**E. Risks.** Order-dependence (must be made permutation-invariant or canonically
ordered); state collapse; hard to interpret; close to MMAML's mode inference.

---

## M8 — Cross-Context Consistency Operator (CCCO)

**A. Object.** A learned operator that extracts only the component of the support
signal that is *reproducible across two or more measurement contexts*, and applies
that component.

**B. Formulation.** For supports spanning contexts `c ≠ c'`, form within-context
contrast statistics separately, then
`z_t = A_θ(agreement between contexts)`; disagreement is discarded rather than
averaged.

**C. Identifiability.** Requires k ≥ 4 spanning ≥ 2 contexts; the estimand is the
agreeing component only, which is low-dimensional by construction.

**D. Meta-learning.** Sources teach what cross-context agreement looks like.

**E. Risks.** **Coverage** — most k ≤ 5 draws will not span two contexts with
enough rows; this is a measurable precondition, and the assay census suggests it
will bind hard.

---

# 5. Comparison

Scores are 1–5, assigned against the criteria in the brief. "Evidence" cites
measurements in this programme, not opinion.

| | M1 CRO | M2 TLCA | M3 APS | M4 CEA | M5 ISC | M6 CNF | M7 ORMT | M8 CCCO |
|---|---|---|---|---|---|---|---|---|
| Novelty | **4** | 4 | 2 | **5** | 2 | 1 | 3 | 3 |
| Meta-learning contribution | **5** | 4 | 4 | 4 | 4 | 2 | **5** | 3 |
| Identifiability at k ≤ 5 | **5** | 4 | **5** | **5** | 3 | 4 | 4 | 2 |
| Biological plausibility | **5** | 3 | 3 | 4 | 2 | 4 | 2 | 4 |
| Mathematical clarity | **5** | 3 | 4 | 4 | 4 | 4 | 2 | 3 |
| Falsifiability | **5** | 4 | 4 | **5** | 4 | 4 | 3 | 4 |
| Open-data reproducibility | **5** | 4 | 5 | 5 | 4 | 5 | 4 | 2 |
| **Nuisance-immune by construction** | **yes** | no | no | **yes** | no | n/a | no | **yes** |
| Nearest negative evidence | assay-coherence null | T0C null | prior art | none direct | **PIRS oracle fail** | forbidden class | none direct | coverage |
| **Total** | **34** | 26 | 27 | **32** | 23 | 24 | 23 | 21 |

---

# 6. Selected mechanism — M1, Contrast Response Operator

**One-sentence contribution.** *We meta-learn how a protein target responds to
chemical contrast within a single measurement context, so that a few paired
comparisons — rather than a few absolute affinities — identify a bounded,
target-specific re-ranking of unseen compounds.*

**Definition.**

```
P_t  = {(j,l) : c_j = c_l}                         within-context support pairs
d_jl = g(x_j) − g(x_l)          Δ_jl = r_j − r_l   contrast, offset-free exactly
z_t  = A_θ(S_t) = Σ_{(j,l)∈P_t} α_θ(‖d_jl‖, |Δ_jl|, k) · sign(Δ_jl) · ψ_θ(d_jl)
ŷ_q  = μ(p_t, x_q) + b · tanh( ⟨ z_t , φ_θ(x_q) ⟩ / b )
```

Outer objective: `min_θ E_t E_{k∈{1,3,5}} L_rank(Q_t ; f_θ(·, A_θ(S_t)))`, with
`L_rank` the smoothed-CI surrogate (the convex RankNet logistic was measured to
mismatch the metric), plus a permuted-label term driven to zero and a bound on
`‖z‖` so a wrong adaptation is survivable.

**Why it solves identifiability.** Three reductions, each measured rather than
assumed. (i) The additive context offset — the component worth more than the entire
chemical head — cancels exactly in `Δ_jl`, so the estimand contains no nuisance.
(ii) Only the *sign and direction* of a contrast enter, not its calibrated
magnitude, so measurement scale drops out. (iii) k = 5 in one context supplies up
to 10 contrasts against a bounded low-dimensional `z_t`, instead of 4 level-free
contrasts against a 26-dim head.

**Why it differs from existing methods.** ActFound and pairwise DTA use differences
to *cancel assay bias in prediction*; CRO makes the contrast the **adaptation
state's input** and meta-learns `α_θ` — which contrasts count as evidence. MAML/ANIL
adapt parameters, not estimands. Neural processes and hypernetworks map supports to
a high-dimensional latent. KRR/ridge/LASSO on contrasts are the required baselines
and cannot learn contrast-dependent reliability. Metric learning acts through
support–query similarity; CRO's query need not resemble any support, because it is
matched to a *contrast direction*, not to a compound.

**Expected experiments.** Gate D0 first (split with simultaneous target / scaffold /
document / assay separation; re-measure the surviving headroom; **stop if it does
not clear the MDE**). Then: synthetic positive control with a planted contrast
response; contrast-coverage census (how many k ≤ 5 draws yield ≥ 1 within-context
pair — a hard precondition); CRO vs contrast-KRR, contrast-ridge, contrast-LASSO,
MAML, ANIL, MetaDTA, AdaMBind, TRACE, pooled head, ligand-only, and the
document-mean oracle at matched capacity; controls for wrong-target support, label
permutation, sign permutation, support removal (exact no-op), and protein ablation;
evaluation on same-document pairs from documents **absent from the support**;
replication on Papyrus 05.7 and SPD 2023 (family transfer only — SPD is
query-depth underpowered), with Davis/KIBA for comparability.

**Registered stop conditions.** Coverage census shows most k ≤ 5 draws have no
within-context pair; gains confined to Tanimoto ≥ 0.55; contrast-KRR matches CRO;
gain vanishes under sign permutation; or D0 returns no admissible headroom.

**Honest risk.** The assay-coherence gate already returned a null when data were
restricted to exact assays. CRO's claim is that the failure there was estimating a
*dense head* on restricted data, not the restriction itself. If a contrast-native
estimand also returns null, the correct conclusion is that open ChEMBL affinity
does not support target-specific few-shot ranking adaptation at k ≤ 5 — and that,
with the document-oracle result, is the programme's terminal deliverable.

---

# 7. Implementation status

**Not started, and not authorized by this document.** Stage 2 may design
architecture and code only after (a) Gate D0 returns an admissible headroom on a
document- and assay-separated split, and (b) the CRO coverage census shows the
mechanism has a deployment path at k ≤ 5.
