# A2S-DTA Stage 1 (v2) — search for a new meta-learning adaptation principle

Date: 2026-08-02 · Branch: `research/a2s-transfer-object-20260802`
Status: **mechanism discovery. No code. No training. No implementation authorized.**
Supersedes `A2S_STAGE1_MECHANISM_DISCOVERY_2026-08-02.md`, whose selected mechanism
(Contrast Response Operator) is retained below as an *instance* of the principle
selected here, not as the contribution.

> ## Precondition that binds every mechanism below
> Gate T0 (rev 2) measured: a chemistry-free document-mean oracle scores **+0.0610**,
> beating a full per-target chemical head (**+0.0519**); within-target
> scaffold-disjoint splits leave **91.1 %** of query rows in a support-seen document;
> the same-document chemical remainder **+0.0290 [+0.0046, +0.0557]** *fails* its
> 0.005 bar; cross-target head transfer is **−0.0183 [−0.0435, +0.0054]**.
> **Gate D0 must re-establish an admissible chemical headroom on a
> target/scaffold/document/assay-separated split before any mechanism is built.**
> If it does not, the terminal deliverable is the negative result.

---

# Part 1 — The current scientific bottleneck

## 1.1 The estimand, decomposed

```
r_t(x) = y_t(x) − μ(x)  =  c_t(x)   target-specific chemical response  [transferable]
                        +  δ_{t,ctx} per-context offset/scale          [nuisance, NOT transferable]
                        +  ε
```

Nine mechanisms competed to estimate `r_t`. Measurement says `δ` alone is worth
more than the whole fitted head. **The programme has been spending its entire
k ≤ 5 label budget estimating a quantity that is majority nuisance — and the
nuisance component is exactly the part that cannot transfer, because an unseen
target's assays and documents are new.**

## 1.2 Why "k ≈ 10" is the observed number

It is not a chemical or biological constant. It is where a dense estimator's
shrinkage begins to recover a signal whose ranking-relevant part is small relative
to residual dispersion — and it moves with the **estimand**, not the estimator:
dense empirical Bayes (G4) and discrete selection (T0) knee in the same region
despite unrelated hypothesis spaces. An estimator cannot fix an estimand.

## 1.3 The bottleneck in one sentence

> **Few-shot target adaptation in open affinity data is nuisance-limited before it
> is information-limited.** The labels are few *and* they are spent on a quantity
> that is not invariant to the measurement process.

---

# Part 2 — Why the previous nine mechanisms failed

| # | Mechanism | Adaptation object | Immediate cause | Root cause |
|---|---|---|---|---|
| 1 | **TRACE** | learned per-pair transport reliability | learned part −0.0001 [−0.0006, +0.0005] over scaled KRR | Failure 1: transport is bounded by support–query distance, and the surviving local effect coincides with shared context |
| 2 | **MODE** | discrete mode dictionary | k-shot inference not separable (A2/A4) | Failure 3: mode selection not identifiable from passive support |
| 3 | **IDA / rank-`m`** | low-rank target code | spectrum flat; rank-2 retains −6 % | Failure 2: no shared low-dimensional structure exists |
| 4 | **RIP** | certified subset of ranking edits | ceiling +0.075 real, margin AUC 0.555 | Failure 5: selection is a decision layer, not an adaptation object |
| 5 | **HOTSPOT** | sparse target-specific coordinates | killed by an exact rotation control | Failure 4: truncation performance was a generic property of a noisy dense ridge |
| 6 | **assay-coherence** | context-restricted residuals | σ 1.26→0.83, still not admitted | restricted the **data** but kept estimating a **dense head** — the estimand was never changed |
| 7 | **MMP / TCRS** | explicit transformation grammar | zero coverage below Tanimoto 0.35 | no deployment path to distant queries |
| 8 | **PIRS** | protein-conditioned interaction coordinates | full-support **oracle** failed (+0.0046, crosses 0) | Failure 6: sequence-derived conditioning carries no signal here |
| 9 | **CFES** | conformational state population | ligand-only beat ligand×pocket by 0.092 | Failure 6 again, now on experimental structures |

**Plus a tenth, self-inflicted:** Gate T0 rev 1 reported positive transfer. It was
retracted — the basis was irreproducible (`torch.svd_lowrank` takes no generator)
and the same-document control was applied to the wrong arm. **Every gate built on
that basis (A0–A4, G1–G4, R0, HOTSPOT) has seed-dependent recorded intervals.**

**Unified account.** Failures 1, 2, 3, 4 are all *estimator or hypothesis-space*
variations on one estimand. Failure 6 removes the label-free conditioning that
might have subsidised the label budget. Nothing in the sequence ever changed
**what is being estimated**. That is the untouched degree of freedom.

---

# Part 3 — Ten candidate mechanisms

Notation: source tasks `T_i=(S_i,Q_i)`; support `S_t={(x_j,y_j,c_j)}` with `c_j` the
measurement context; frozen base `μ(p_t,x)`; residual `r_j`; ligand features `g(x)`;
protein `p_t`; state `z_t=A_θ(S_t)`; adaptation `Δ_t(x)=g_θ(x,p_t,z_t)`; objective
`min_θ E_T[L_Q(f_θ(Q,z_t))]` throughout.

---

### M1 — Adaptation-Complexity Minimizing Representation *(Direction A)*

**Hypothesis.** Meta-learning can transform the representation so that few labels
suffice, by optimizing adaptation complexity directly rather than variance.

**Formulation.** `φ_θ(x,p)`; per-target head solved in closed form inside the
episode; outer loss is the *k-shot* ranking loss plus a differentiable penalty on
the effective degrees of freedom `tr(H_t)` of the inner solve. The learned object
is `φ_θ` with minimal `C(f_t|φ_θ)`.

**Identifiability.** Complexity is the optimization target, so `dof ≤ k−1` is
enforced during meta-training rather than hoped for.

**vs prior.** PIRS learned coordinates but optimized prediction, not complexity;
MAML optimizes an initialization; MetaOptNet optimizes through a solver without
constraining dof.

**Falsification.** Planted low-complexity world; fixed-representation baseline;
report measured dof vs k-curve; ligand-only; protein shuffle.

**Risk.** Complexity penalties can be satisfied by discarding signal — must show
the ceiling does not fall with the dof.

---

### M2 — Contrast-Sufficient Representation *(Direction A)*

**Hypothesis.** A representation exists in which *within-context contrasts* are a
sufficient statistic for the target's ranking function.

**Formulation.** Train `φ_θ` with an explicit sufficiency loss: a head predicted
from contrasts alone must match the full-support head. `z_t` is then a function of
contrasts only.

**Identifiability.** Sufficiency means no information is lost by discarding
absolute values — so the nuisance-free estimand is not a sacrifice.

**vs prior.** No DTA method states or tests a sufficiency claim for contrasts.

**Falsification.** Sufficiency gap measured directly; affine-relabel invariance;
wrong support; k-curve.

**Risk.** Sufficiency may simply fail — a real and informative outcome.

---

### M3 — Coordinate-Generating Adaptation *(Direction B)*

**Hypothesis.** Coordinates themselves adapt: `φ_t(x)=g_θ(x,z_t)`, with the
*identity* of the coordinate chosen label-free and only its amplitude/sign read
from labels.

**Formulation.** `z_t = A_θ(S_t^x, U_t)` (support inputs + unlabeled query library,
no labels) selects a coordinate; a scalar `a_t` from labels sets amplitude;
`Δ_t(x) = a_t·⟨φ_θ(x), v(z_t)⟩`.

**Identifiability.** Label burden is **one scalar**, independent of ambient
dimension — the reason a learned coordinate can need k ≤ 5 where a fixed basis
needs k ≈ 10.

**vs prior.** MODE selected among fixed modes using labels; here selection is
label-free and labels only scale.

**Risk.** T0C measured label-free shortlisting as null (protein +0.005, chemotype
−0.010). Direct negative evidence; a learned encoder is not obviously better.

---

### M4 — Budget-Nested Coarse-to-Fine Coordinates *(Direction B)*

**Hypothesis.** The transferable object is a *nested* coordinate hierarchy in which
every prefix is usable, so adaptation resolution is set by `k`.

**Formulation.** Ordered coordinates `v_1,v_2,…`; at budget `k` only the first
`⌊(k−1)/2⌋` are activated by a structural mask. `k=1` → exact no-op.

**Identifiability.** Structural, not regularized: the model cannot spend dof it
does not have.

**vs prior.** IDA/rank-`m` fixed one rank; this makes rank a function of budget and
trains all prefixes jointly.

**Risk.** Reduces to low-rank adaptation, which G2 refuted, unless the nesting is
load-bearing.

---

### M5 — Amortized Adaptation Operator *(Direction C)*

**Hypothesis.** The transferable object is the operator `𝒜: S ↦ Δ`, never a state.

**Formulation.** A set-encoder maps `S_t` directly to a bounded function `Δ_t(·)`
in a small learned function class.

**Identifiability.** Only if the output class is explicitly bounded.

**vs prior.** This is a neural-process/CNAPs-style amortization.

**Risk.** **Degenerates into a hypernetwork generating a large head** — the brief's
named failure. Low novelty.

---

### M6 — Nuisance-Equivariant Adaptation *(Direction C)* — **SELECTED**

**Hypothesis.** The transferable object is not a state or an operator but the
**symmetry class the operator belongs to**: adaptation must be *equivariant to the
nuisance group acting on measurements*, so meta-learning happens in the quotient
space. Full treatment in Part 5.

---

### M7 — Learned Intervention Axes *(Direction D)*

**Hypothesis.** k labels identify which of a few learned chemical intervention axes
is active for this target.

**Formulation.** Learned axes `{i_1..i_m}`, `m ≤ 4`, each a bounded response
direction; `z_t` is a posterior over which axis is active; `Δ_t(x)=Σ z_m·i_m(x)`.

**Identifiability.** ≤ 2 bits.

**vs Free-Wilson.** Axes are learned, continuous and target-conditioned rather than
substituent-position indicators fitted per series.

**Risk.** **Failure 3 applies**: mixture/mode selection without a *new
identifiability principle* is explicitly rejected. Needs one to survive.

---

### M8 — Reversal-Detector Interventions *(Direction D)*

**Hypothesis.** Targets differ mainly in **which SAR trends reverse**; a reversal is
one bit and is the natural unit of selectivity.

**Formulation.** Learned global SAR axes; `z_t ∈ {−1,0,+1}^m`, `m ≤ 3`, inferred
from contrast signs; `Δ_t(x)=Σ_m z_{t,m}·s_m(x)`.

**Identifiability.** Sign-only estimand: robust to offset *and* scale; ≤ 3 bits.

**Biology.** Selectivity between related targets is often a reversed preference at a
few positions.

**Risk.** Sign information may be too weak at measured noise; overlaps M6's
invariance argument without M6's generality.

---

### M9 — Design-Conditioned Adaptation *(Direction E)*

**Hypothesis.** How much to adapt should depend on the *geometry* of the realized
support set, which is label-free and currently ignored.

**Formulation.** `z_t = A_θ(contrasts, Γ_t)` where `Γ_t` summarizes the support
design — Gram spectrum, leverage, context spread, chemical diversity.

**Identifiability.** `Γ_t` is label-free, so it costs no bits.

**Motivation from measurement.** Support draws are chemically coherent (same series,
same document), which is precisely the rank-deficient regime; nothing has ever
conditioned on it.

**vs prior.** Not calibration — it modulates the *main path's* adaptation, not a
confidence score.

**Risk.** A multiplier, not a source: if there is no adaptation signal, conditioning
on geometry gains nothing. **Best used as a module inside another mechanism.**

---

### M10 — Identifiability-Aware Episode Geometry *(Direction E)*

**Hypothesis.** Meta-training should span the deployment distribution of episode
geometries, and learn geometry-indexed operators.

**Formulation.** Episode sampler stratified by design geometry; operators indexed by
geometry stratum; compared against random, diversity and D-optimal support draws.

**Identifiability.** Unchanged per episode; improves the *estimator's* match to the
realized geometry.

**Risk.** In the sealed passive A2S protocol supports cannot be chosen at
deployment, so D-optimal is a **diagnostic upper bound**, not a deployable method.

---

# Part 4 — Novelty comparison with the literature

## Meta-learning

| Method | Transferable object | k ≤ 5 unseen-target DTA | Gap |
|---|---|---|---|
| MAML | initialization `θ₀` | inner SGD in dim ≫ k−1 is noise-dominated | no control of adaptation dof |
| ANIL | features + head init | = linear-head adaptation = our EB head, knee k≈10 | same estimand |
| MetaOptNet / R2D2 | features + convex inner solver | optimal estimator in a fixed space | space is the problem |
| CNP / NP / ANP | amortized support → global latent | latent continuous, high-dimensional | no bound on identifiable dimension |
| CNAPs | FiLM parameters from a support encoder | support-conditioned modulation of a frozen backbone | no nuisance model |
| Hypernetworks / HyperPCM | descriptor → target parameters | needs label-free target conditioning; null ×4 here | conditioning signal absent |
| Modular Meta-Learning / MMAML | reusable modules + task-mode inference | **covers M7 directly** | mode selection unidentifiable (Failure 3) |
| Information-Theoretic Meta-Learning | task representation under an information constraint | constrains *how much* task info, not *which nuisance* to exclude | **constrains quantity, not invariance** |

## Representation and invariance

| Line | Object | Relation to this proposal |
|---|---|---|
| Invariant representation learning (IRM, DANN) | representations invariant to environment | invariance imposed on **features**; here it is imposed on the **adaptation operator** |
| Sufficient statistics / sufficient representation | minimal sufficient encoding | M2 uses it; M6 pairs it with a *maximal invariant* |
| Task-aware / disentangled adaptation | separate task and content factors | disentanglement is learned and soft; M6's separation is a **group action, exact by construction** |

## Biological ML

| Line | Object | Relation |
|---|---|---|
| Free–Wilson | additive substituent contributions per series | M7/M8 differ by learning axes; falsified in its explicit form by the MMP census |
| SAR modelling / matched pairs / ActFound | pairwise differences to cancel assay bias **in prediction** | M6 makes the invariant the **input to adaptation** and formalizes the group |
| Binding-determinant / interaction-fingerprint learning | explicit contact features | measured null here (PIRS, CFES-C0B) |
| Protein–ligand foundation models | large pretrained joint encoders | forbidden as a contribution; ligand-only beat ligand×pocket on contact prediction |

**The gap, stated once.** Invariance is standard for *features* and for
*distribution shift*. **No meta-learning method defines the adaptation operator by
the nuisance group of the measurement process, and no DTA method identifies the
maximal invariant statistic of that group as the information ceiling for few-shot
adaptation.**

---

# Part 5 — Mathematical formulation of the selected mechanism

## 5.1 The nuisance group

Affinity is a context-conditional observable. Within measurement context `c`,
reference conditions shift the scale and the readout compresses the range:

```
G :  y ↦ a_c · y + b_c ,      a_c > 0 ,   independently per context c
```

`b_c` is measured and dominant (the document-mean oracle). `a_c` is expected from
assay dynamic-range compression and is testable.

## 5.2 The principle

> **An adaptation operator must be invariant under `G`.**
> `A_θ(g · S_t) = A_θ(S_t)` for all `g ∈ G`.

Meta-learning then happens in the quotient space `S/G` — the space of what the
measurement process has *not* destroyed.

## 5.3 The maximal invariant statistic

For a context `c` with `k_c ≥ 3` supports, the maximal invariant of the affine
group acting on `(r_1..r_{k_c})` is the **studentized contrast vector**

```
m_c = ( r_j − r̄_c ) / s_c ,   s_c = within-context dispersion
```

equivalently the ordering plus normalized gaps. Any `G`-invariant function of the
supports factors through `m_c`. **This is the information ceiling under nuisance,
and it is a theorem, not a modelling choice.**

Dimension accounting per context: `k_c − 2` invariant degrees of freedom (one lost
to offset, one to scale).

| `k` in one context | raw | offset-free | **`G`-invariant** |
|---:|---:|---:|---:|
| 1 | 1 | 0 | 0 |
| 3 | 3 | 2 | **1** |
| 5 | 5 | 4 | **3** |

Under the offset-only subgroup `G₀` (`a_c ≡ 1`) the budget is `k_c − 1`. Which
group applies is **an empirical question with a preregistered test**, not an
assumption.

## 5.4 The mechanism

```
m_t   = { studentized within-context contrasts of S_t }        (maximal invariant)
Γ_t   = label-free design geometry of S_t                       (Gram spectrum, spread)
z_t   = A_θ( m_t , Γ_t )                    ∈ R^d_z ,  d_z ≤ 3  (bounded by construction)
Δ_t(x)= b · tanh( ⟨ z_t , φ_θ(x) ⟩ / b )
ŷ_q   = μ(p_t, x_q) + Δ_t(x_q)
```

`z_t = 0` whenever no context has `k_c ≥ 3`, so support removal and `k=1` are
**exact no-ops**. Objective:

```
min_θ  E_T E_{k∈{1,3,5}} [ L_rank( Q_T ; f_θ(·, A_θ(S_T)) ) ]
```

with `L_rank` the smoothed-CI surrogate (the convex RankNet logistic was measured
to mismatch the metric), plus a permutation term driven to zero on label-shuffled
batches. Three modules only: contrast encoder, invariant aggregator, bounded head.

---

# Part 6 — Selected mechanism

## M6 — Nuisance-Equivariant Adaptation (NEA)

**One-sentence contribution.** *We define few-shot target adaptation on the quotient
space of the measurement nuisance group — meta-learning an operator that consumes
exactly the maximal invariant statistic of the support set — turning the dominant
confound in open affinity data from a control to be checked into a symmetry the
mechanism is built from.*

**Transferable object.** The `G`-equivariant adaptation operator: *how a target's
invariant contrast pattern maps to a bounded re-ranking*. Not a state, not a
similarity, not a target embedding.

**Why k ≤ 5 is possible.** Three reductions, each measured rather than assumed.
(i) The nuisance component — worth **more than the entire fitted chemical head** —
is removed by the *definition* of the estimand, not by a post-hoc control, so the
whole label budget is spent on transferable signal. (ii) The estimand is the
maximal invariant, of dimension `k_c − 2`, against a state bounded at `d_z ≤ 3` —
so the model cannot spend degrees of freedom it does not have. (iii) `Γ_t` supplies
episode-geometry conditioning at **zero label cost**.

**Why previous methods cannot do this.**
- **TRACE** acts through support–query similarity; NEA's query need not resemble any
  support, because it is matched to an invariant *contrast pattern*.
- **MODE / IDA** select or project in a fixed space using absolute residuals — the
  nuisance enters the selection statistic itself.
- **RIP** is a decision layer over an already-contaminated estimate.
- **PIRS / CFES** change the representation while keeping the estimand.
- **CMAL / MAML / ANIL** adapt parameters, not estimands, and their inner loops
  consume absolute labels.
- **AdaMBind / MetaDTA** are support-conditioned heads on the same contaminated
  quantity.
- **ActFound** uses pairwise differences for *prediction*; NEA makes the invariant
  the input to *adaptation*, formalizes the group, and tests which group applies.

**Minimal architecture (3 modules).** Contrast encoder `ψ_θ` over studentized
within-context contrasts → invariant aggregator `A_θ` (permutation-invariant, also
consuming `Γ_t`) → bounded intervention head over `φ_θ(x)`. Frozen base untouched.

**Falsifiability — the property no competitor has.** `G`-invariance is a
**numerical, not statistical, test**: apply a random per-context affine transform to
the support labels and predictions must be **bit-identical**. A mechanism that
fails this is not the mechanism.

---

# Part 7 — Validation gates, in order (design only)

| Gate | Question | Pass condition | If it fails |
|---|---|---|---|
| **D0** | does *any* chemical headroom survive proper separation? | simultaneous target/scaffold/document/assay-separated split; headroom lower bound > MDE on same-document pairs from documents absent from support | **stop the programme**; write the negative result |
| **N0** | which nuisance group actually acts? | test `G₀` (offset) vs `G` (affine): does per-context scale vary materially? | adopt `G₀` and gain one dof per context |
| **N1** | coverage | fraction of k ≤ 5 draws with some context holding ≥ 3 supports; must support a deployment path | mechanism has no deployment path at k ≤ 5 — stop |
| **N2** | invariance | bit-identical predictions under random per-context affine relabeling | implementation defect, not evidence |
| **N3** | sufficiency (M2's claim) | contrast-only head vs full-support head: measured sufficiency gap | quantifies the price of invariance |
| **N4** | synthetic positive control | planted invariant contrast response recovered at k = 3, 5; k = 1 exactly chance | no power — no negative may be reported |
| **N5** | mechanism vs matched baselines | beats contrast-KRR, contrast-ridge, contrast-LASSO, MAML, ANIL, MetaDTA, AdaMBind, TRACE, pooled head, ligand-only, **document-mean oracle**, at matched capacity and budget | not a meta-learning contribution |
| **N6** | shortcut battery | wrong-target support, label permutation, **sign permutation**, support removal (exact no-op), protein shuffle, random representation | gain is not adaptation |
| **N7** | k-curve | monotone over k ∈ {1,3,5,10,20}, with k = 1 an exact no-op | mis-specified budget scaling |
| **N8** | generalization | Papyrus 05.7; SPD 2023 (family transfer only — query-depth underpowered); Davis/KIBA for comparability; multi-seed | benchmark-specific result |

**Registered stop conditions.** N1 shows most k ≤ 5 draws have no usable context;
contrast-KRR matches NEA (then it is an estimator, not a mechanism); gains confined
to Tanimoto ≥ 0.55; gain survives sign permutation; or D0 returns no admissible
headroom.

**Data.** ChEMBL 37 dual-cold (present, hashed), Papyrus 05.7, SPD 2023 (CC BY-4.0),
Davis/KIBA by public download. Deterministic preprocessing; the basis is now an
exact sign-fixed eigendecomposition after the `svd_lowrank` defect.

**Honest prior negative.** The assay-coherence gate returned null when *data* were
restricted to exact assays. NEA's claim is that the failure was fitting a **dense
head on restricted data**, not the restriction itself — the estimand was never
changed. If an invariant-native estimand also nulls, the correct conclusion is that
open affinity data does not support k ≤ 5 target-specific ranking adaptation, and
that — with the document-oracle result — is the programme's terminal deliverable.

**Implementation remains unauthorized** until D0, N0 and N1 pass.
