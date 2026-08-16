# A2S few-shot meta-adaptation — mechanism design

Date: 2026-08-02 · **Revision 2 (supersedes revision 1 of the same date)**
Branch: `research/a2s-transfer-object-20260802`
Status: **withdrawn as a Stage 1 proposal; retained as a design record under
revision.** Nothing is trained. **Gate F1 does not run as written.**
Evidence base: `A2S_TRANSFER_OBJECT_GATE_T0_DECISION_2026-08-02.md` (revision 2)
and the nine falsified predecessors listed in §1.

> ## What changed in revision 2
>
> Revision 1 selected A2S-FBA on the strength of T0's transfer result. **That
> result is withdrawn** — selected-head transfer does not survive the
> same-document control (−0.018 [−0.044, +0.005]), and T0A's chemical remainder
> does not clear its own admission threshold. The premise the mechanism was built
> on is gone, so the selection is void.
>
> Three further defects in revision 1's own design, independent of T0:
>
> 1. **The entropy term was wrong.** Minimizing `+λ_e·H(z)` rewards *confident*
>    routing, including confidently wrong routing and always choosing one
>    operator. It does not reward support-dependent information. §7 replaces it.
> 2. **The harmlessness term did not certify the deployed predictor.** It scored
>    each `μ + ω_m` while deployment uses `μ + Σ_m (z_m − π_m) ω_m`. Individually
>    harmless operators do not make the combination harmless, and near-identical
>    near-zero operators would make every operator "useful" while adapting
>    nothing. §7 replaces it.
> 3. **k = 1 was structurally excluded.** The exact-no-op property means FBA
>    cannot satisfy an "improvement across k = 1, 3, 5" criterion at all. §11
>    states this as a scope limitation rather than leaving it implicit.
>
> The candidate search in §3 is also insufficient: six entries, of which one is a
> nuisance control, one is out of scope and one is a forbidden larger model. §3b
> adds the required method-specific analysis (MetaDTA, AdaMBind, neural processes,
> hypernetworks, and the modular/mixture prior art), which shows FBA as specified
> is **not yet distinguishable from existing modular meta-learners**.
>
> *Note on scope:* the brief I was working from asked for a "candidate mechanism
> search" without a count; the eight-candidate requirement reached me only through
> the review. §3b supplies the analysis either way.

---

# 1. Problem diagnosis

## 1.1 Nine mechanisms, one failure

| # | Mechanism | Adaptation object | Verdict |
|---|---|---|---|
| 1 | TRACE | learned per-pair transport reliability | learned part −0.0001 CI over scaled KRR |
| 2 | MODE | discrete mode dictionary | k-shot inference not separable (A2/A4) |
| 3 | IDA / rank-`m` code | low-rank target code | source-head spectrum flat; rank-2 retains −6 % (G2) |
| 4 | RIP | certified subset of ranking edits | ceiling +0.075 real, margin AUC 0.555 — unreachable |
| 5 | HOTSPOT | sparse coordinates on a target-specific support | falsified by an exact rotation control |
| 6 | assay-coherence | context-restricted residuals | σ 1.26 → 0.83, still not admitted |
| 7 | MMP grammar / TCRS | transferable chemical transformation rules | zero coverage below Tanimoto 0.35 |
| 8 | PIRS | protein-conditioned interaction coordinates | full-support oracle failed (+0.0046, CI crosses 0) |
| 9 | CFES | conformational state population | C0B: ligand-only beats ligand×pocket by 0.092 |

Every one of them changed **the estimator or the basis**, and left the *information
available to the estimator* untouched. Gate T0 measured that information directly.

## 1.2 What Gate T0 established

**(a) 60 % of the prize was never a chemical object.** A chemistry-free oracle
that knows only a compound's ChEMBL document and that document's mean residual
scores **+0.0610 [+0.0386, +0.0824]**, beating the full 26-dimensional per-target
chemical head (+0.0527). Ligand features recover document identity because a
document reports a congeneric series. The offset-free chemical headroom, measured
on same-document query pairs, is **+0.0313 [+0.0056, +0.0601]**.

This also dissolves the programme's oldest puzzle: Q1's "information only at
Tanimoto ≥ 0.55" is the same stratum as "same measurement context".

**(b) The transfer premise is true — at full support.** Selecting one source
target's whole head using ~64 recipient labels gives **+0.0266 [+0.0054, +0.0484]**
over the frozen base and **+0.0346 [+0.0152, +0.0551]** over a single pooled head.
This is the programme's first positive transfer result.

**(c) The library is mostly poison.** Only 40.7 % of source heads beat the base;
the median head scores **−0.0257**; on average **12.9 of 110** heads reach half the
best head's gain.

**(d) Selection break-even is k ≈ 20, and the deficit is exactly accountable.**

```
available(k) = (k-1) · ½ · log2(1 + τ²/σ²)          τ = 0.642, σ = 1.343
             = (k-1) · 0.1486 bits          →  0.59 bits at k=5
required     = log2(M / M_useful) = log2(110 / 12.9) = 3.09 bits
```

Equating gives k ≈ 21.8; the **measured** break-even is k = 20. No free parameter.

## 1.3 The diagnosis

> The binding constraint is not the estimator, not the ligand basis, not the
> protein representation, and not the optimizer. It is that the **hypothesis
> space of adaptations is 5× too large for the information the support labels
> carry**. Every previous mechanism inherited a hypothesis space it did not
> design — 26 descriptor coordinates, 110 source heads, a k-means dictionary, a
> rank-`m` subspace — and then built a better estimator inside it.

A better estimator cannot create bits. Only a smaller, better-shaped hypothesis
space can reduce the bits required.

---

# 2. Literature positioning

| Family | Representative | What it adapts | Why it does not answer this problem |
|---|---|---|---|
| Gradient meta-learning | MAML, Meta-SGD, ANIL, LEO | continuous parameters via inner SGD | inner step is a noisy gradient in a space of dimension ≫ k−1; nothing constrains the space to be identifiable. ANIL's linear-head adaptation is dominated by our EB head, whose knee is k≈10 (G4) |
| Closed-form inner solves | R2D2, MetaOptNet, iMAML, Bayesian/ALPaCA | ridge / posterior over a fixed feature space | this is the *optimal* estimator given the space. T0 shows the space is the problem; our EB head is exactly this and fails |
| Metric / retrieval few-shot | ProtoNet, matching nets, Matsy, ActFound pairwise | similarity-weighted support labels | ruled out below Tanimoto 0.35 by Q1/TRACE, and T0A shows the surviving local effect is measurement context |
| Target-conditioned DTA | HyperPCM, DrugBAN, PSICHIC, DeepDTA/GraphDTA | protein-conditioned parameters or interaction maps | zero-shot protein conditioning measured null three times here (G3 pooled ESM; PIRS segments; CFES-C0B pocket composition) |
| Structure/physics | ensemble docking, IPBind, PReorg-FEP | state-specific energies | requires state semantics that C0B could not admit on this substrate |
| Sparse / compressed adaptation | LASSO, OMP, LISTA, spike-and-slab | sparse coordinates | falsified here by a rotation control: rotating the basis and truncating does *better* (76.8 % vs 63.7 %) |
| Conformal / selective prediction | Fisch et al. meta-conformal, RIP | which predictions to trust | ceiling real (+0.075) but margin AUC 0.555 — no observable is sharp enough |

**The gap.** Every family above optimizes an estimator, a representation, or a
decision rule **given** a hypothesis space. None of them optimizes *the hypothesis
space itself for identifiability under a fixed measurement budget*. That is the
open position this design occupies.

The closest prior ideas — and the honest prior-art boundary — are modular
meta-learning / mixture-of-experts task routing (a discrete code chosen per task)
and information-bottleneck task embeddings. Neither is claimed as novel. What is
new is the **meta-objective**: making the operator library's *k-shot posterior
entropy* and its *worst-case harm* first-class training terms, with a
measurement-derived target (`M/M_useful ≤ 1.5`) rather than a hyperparameter.

---

# 3. Candidate mechanism search

All candidates are constrained by T0: they must reduce required bits below ~0.6,
must not rely on label-free protein/library conditioning (T0C null), and must be
evaluated on same-document pairs (T0A).

| # | Candidate | Transferable object | Bits required | Assessment |
|---|---|---|---|---|
| **F1** | **Few-Bit Adaptation (FBA)** — meta-learned library of `M ≤ 4` operators, jointly optimized for expressiveness, harmlessness and k-shot discriminability | the library + its selection posterior | `log2(M/M_useful)` → target ≤ 0.6 | **selected**. Directly attacks the measured constraint; the only candidate whose success criterion is a number T0 already produced |
| F2 | Context-normalized adaptation | a learned map from `(k` supports, provenance`)` to a measurement-context correction | ~1 | T0A says this is 60 % of the historical gain — but it is a *nuisance*, is unavailable for a genuinely new compound, and the user's brief forbids post-hoc calibration. **Retained as a mandatory control arm, not a mechanism** |
| F3 | Amplitude-only adaptation | one scalar on a label-free predicted direction | ~0 (1 continuous dof) | needs a label-free direction; T0C/G3/PIRS/C0B measured that as null four times. **Rejected on measured grounds** |
| F4 | Pairwise-comparison adaptation | a learned operator on within-support *orderings* rather than values | ~0.6 at k=5 | orderings are offset-free, so this is immune to T0A's confound. Weaker signal (`k(k−1)/2` binary contrasts at low SNR). **Retained as the F1 fallback and as an ablation of the residual encoder** |
| F5 | Active/diverse support policy | the support-selection rule itself | raises `τ²` | changes the task from passive to active; out of scope for the sealed A2S protocol. **Deferred**, but registered: T0's incoherence caveat predicts diverse supports help |
| F6 | Larger encoder / adapter / Transformer | none | unchanged | explicitly excluded by the brief and by nine prior results |

## 3b. Method-specific novelty analysis (added in revision 2)

| Prior method | What it adapts | Distance from FBA as specified |
|---|---|---|
| [Modular Meta-Learning](https://arxiv.org/abs/1806.10166) (Alet et al.) | selects a composition of reusable modules per task by structured search | **Nearest prior art.** A small library of reusable modules chosen per task *is* the FBA object. FBA differs only in using a soft posterior instead of discrete search |
| [MMAML](https://arxiv.org/abs/1910.13616) | modulates initialization by an inferred discrete-ish task mode, then adapts | covers "infer a task mode from support, then modulate" — FBA's inner loop with the gradient step removed |
| [Online mixtures of tasks](https://arxiv.org/abs/1812.06080) | Dirichlet-process mixture over task clusters | covers the mixture-over-task-modes posterior directly |
| [CNAPs](https://arxiv.org/abs/1906.07697) | FiLM parameters generated from a support-set encoder | a support-conditioned adaptation state modulating a frozen backbone — FBA's architecture with a richer state |
| [VERSA](https://arxiv.org/abs/1805.09921) | amortized posterior over task-specific head parameters | the amortized version of our EB head |
| [Information-Theoretic Meta-Learning](https://arxiv.org/abs/2009.03228) | task representation under an explicit information constraint | **pre-empts FBA's headline claim** that constraining task-inference information is the contribution |
| [Meta-Learning without Memorization](https://arxiv.org/abs/1912.03820) | forces support-dependence via an information penalty | pre-empts "make the state genuinely support-dependent" as a novelty |
| MetaDTA / few-shot DTA | support-conditioned regression heads on frozen DTA encoders | the direct application-domain baseline; must be run, not just cited |
| AdaMBind | adaptive binding-affinity meta-learning | same domain, same k-shot framing; a required baseline |
| Neural processes (CNP/ANP) | amortized conditioning on a support set to a global latent | a general form of "support set → adaptation state"; FBA is a discrete-latent special case |
| Hypernetworks / HyperPCM | generate target-specific parameters from a target descriptor | excluded here on measured grounds (label-free target conditioning null ×4), not on novelty grounds |

**Conclusion.** A small operator library plus a support-conditioned posterior is a
modular mixture-of-experts task-mode meta-learner, and information-constrained task
inference is already explicit in the literature above. **The entropy-plus-harm
combination is an application-specific objective, not a new meta-learning
mechanism.** Any future claim must be narrowed to something the prior art does not
already cover — the most defensible remaining candidate is the *measurement-context
identifiability* result itself, which is a property of the data, not of a model.

---

# 4. Selected innovation

> **A2S-FBA — Few-Bit Adaptation.** Meta-learn a *small library of adaptation
> operators together with the selector that identifies them*, under an objective
> that explicitly minimizes the number of bits the recipient's support labels
> must supply, subject to each operator being non-harmful.

**The transferable object** is the pair `(Ω, A_θ)`:

- `Ω = {ω_1 … ω_M}`, `M ≤ 4` bounded ranking operators over ligand features,
  learned once on source targets and frozen at deployment;
- `A_θ`, a support-conditioned selector that maps `k` label contrasts to a
  posterior over `Ω`.

**What is learned from source targets.** Not "how targets differ" (G2: no shared
low-rank structure) and not "which compounds are similar" (Q1/TRACE: local only).
What is learned is **a partition of target-adaptation behaviour into a few modes
that are simultaneously (i) useful, (ii) mutually distinguishable from four noisy
contrasts, and (iii) individually safe.** T0 proves such a structure is not
automatic: the natural library of 110 source heads has `M/M_useful = 8.5`.

**Why it transfers.** Because the library is chosen so that *most members help
most targets*. Transfer does not require the recipient to resemble any particular
source target — a requirement T0C showed is unpredictable — only that the
recipient's adaptation lies near one of a few broad behavioural modes.

**Why not MAML/ANIL/ridge/KRR/retrieval/Bayes.** All six are estimators over a
hypothesis space they do not shape. Formally: for a fixed space `H` and channel
capacity `C(k)`, the minimax excess risk is bounded below by a function of
`log|H_eff| − C(k)`, and that quantity is **estimator-independent**. T0 measured
`log|H_eff| = 3.09` bits and `C(5) = 0.59` bits. FBA is the only candidate that
changes the *first* term. Each of the six is retained as a required baseline (§9).

---

# 5. Mathematical formulation

Frozen support-free base `μ(p_t, x)` (unchanged; never retrained). Ligand features
`g(x) ∈ R^d`, `d = 26`, label-free and fit-role-only.

**Operators.** Each operator is a bounded scalar field on ligand space:

```
ω_m(x) = b · tanh( u_mᵀ g(x) / b ),      m = 1 … M,   b fixed
```

`u_m ∈ R^d` are meta-parameters. Boundedness makes a wrong selection survivable —
a direct response to the measured fact that the median library member is harmful.

**Support evidence — level-free by construction.** With support `S_t = {(x_i, y_i)}`,
residuals `r_i = y_i − μ(p_t, x_i)`, define centred contrasts

```
r̃_i = r_i − mean_j r_j ,        ω̃_m(x_i) = ω_m(x_i) − mean_j ω_m(x_j)
```

so the unidentifiable target level never enters, and `k = 1` yields exactly zero
contrasts.

**Inner loop — the adaptation state.** The state is a posterior over the library:

```
z_t = A_θ(S_t) = softmax_m [ log π_m − ( 1 / T_θ(k) ) · Σ_i ℓ_θ( r̃_i , ω̃_m(x_i) ) ]
```

- `π` is a learned label-free prior over operators (uniform unless T0C is ever
  overturned);
- `ℓ_θ` is a **learned** discrepancy, not a fixed squared error;
- `T_θ(k)` is a learned budget-dependent temperature, so the mechanism knows how
  much to trust `k` labels. `T_θ(1) = ∞` is enforced structurally.

**Prediction.** The intervention is relative to the prior, so support removal is
an exact no-op:

```
ŷ_q = μ(p_t, x_q) + Σ_m ( z_{t,m} − π_m ) · ω_m(x_q)
```

At `k = 1`, `z_t = π` and `ŷ_q ≡ μ`. **The base is recovered exactly**, not
approximately.

**Outer loop.**

```
min_θ,Ω   E_{t ~ source} E_{k ~ {1,3,5}} E_{S_t, Q_t} [ L_rank( Q_t ; f_θ(·, z_t) ) ]
        + λ_h · E_t Σ_m [ harm_m(t) ]_+          (harmlessness)
        + λ_e · E_t [ H(z_t) ]                    (discriminability)
```

- `L_rank` is the smoothed-CI surrogate already validated in TRACE (the convex
  RankNet logistic was measured to mismatch the metric).

**Both regularizers as written in revision 1 were wrong, and are replaced.**

*Identifiability.* `+λ_e·H(z_t)` rewards confident routing — including confidently
wrong routing, and the degenerate solution of always selecting one operator. The
quantity that actually measures support-dependence is the **mutual information
between the support set and the adaptation state**, estimated across a batch of
targets:

```
I(z ; S) ≈ H( E_t[z_t] ) − E_t[ H(z_t) ]                (maximize)
```

The first term forbids collapse onto one operator; the second rewards decisiveness
*per target*. Only the difference is optimized. A permuted-label batch must drive
this to zero, and that is a training-time assertion, not a post-hoc check.

*Harmlessness.* `harm_m(t)` scored each `μ + ω_m` in isolation, but deployment uses
the mixture. Individually harmless operators do not certify the combination, and
`M` near-identical near-zero operators would score perfectly while adapting
nothing. The penalty is therefore applied to the **deployed predictor** under the
realized posterior, plus an explicit anti-degeneracy term:

```
harm(t) = [ L_rank(Q_t ; μ) − L_rank(Q_t ; μ + Σ_m (z_{t,m} − π_m) ω_m) ]_+
separation = min_{m≠m'} || ω_m − ω_{m'} ||   (lower-bounded, or the library collapsed)
```

`λ_h`, `λ_e` and the separation floor are the only tuning knobs, selected on
**nested fit-only inner validation** and frozen before any held evaluation.

---

# 6. Identifiability argument

**What k ≤ 5 measurements can identify.** After level removal, `k − 1` contrasts at
measured SNR 0.229 carry `(k−1)·0.1486` bits: **0, 0.30, 0.59** bits at k = 1, 3, 5.

**What is impossible.** A dense 26-dof head (needs ≈ 26 dof against 4). Selection
among 110 heads (3.09 bits). Any continuous target latent of dimension > 1.
Any ranking change at k = 1 from support labels — structurally zero contrasts,
confirmed by the T0 synthetic control recovering exactly chance at k=1.

**What FBA needs.** `log2(M / M_useful) ≤ 0.59`, i.e. **`M/M_useful ≤ 1.5`**.
With `M = 3` this means ≥ 2 of the 3 operators must be useful for a typical
target. The measured library value is 8.5, so the training objective must improve
this by ~5.7×; the harmlessness penalty is precisely the term that does it.

**Why meta-learning reduces sample complexity here — the statistical statement.**
Training across many source tasks does not sharpen an estimator. It **selects the
hypothesis class**. Sample complexity for identifying a member of a finite class
scales as `log|H| / SNR`; by shrinking `|H_eff| = M/M_useful` from 8.5 to ≤ 1.5,
the required `k` falls from ≈ 20 to ≈ 5 **at unchanged SNR and unchanged
estimator**. This is the whole mechanism, and it is falsifiable as a number.

**The tension that makes it a real research question.** Shrinking `M` also shrinks
the achievable ceiling. The quantity to maximize is
`ceiling(M) × identifiability(M, k)`, and whether its maximum at k = 5 exceeds the
0.005 MDE is **not known** — it is exactly what Gate F1 measures before any
training.

---

# 7. Training objective and 8. minimal architecture

Frozen: the base `μ`, the label-free basis `g`. Trainable: `{u_m}` (`M·26 ≤ 104`
parameters), `ℓ_θ` (a 2-layer MLP on `(r̃, ω̃)`, ≤ 600 parameters), `T_θ` (one
scalar per `k`), `π` (`M − 1` logits). **Total ≈ 750 parameters.** Capacity is not
the contribution and must not become one; every baseline in §9 is run at matched
or greater capacity.

Episodes: source `fit` targets only, `k ∈ {1,3,5}`, three support policies
(`random_within_target`, `scaffold_disjoint`, `provenance_disjoint`), seeds
{1729, 1730, 1731}. Aggregation order: episode draws → seed/target mean →
component mean → paired component bootstrap (≥ 2000 draws). MDE 0.005 at ~50
components.

---

# 9. Open-data reproducible experimental design

**Primary substrate (present, hashed).** ChEMBL 37 dual-cold registry,
`dataset/public/chembl_37/processed/dualcold/` — 343 211 rows, lock content hash
`6bcf6edc…`; roles `fit`/`probe` open, `locked` and the recipient roster sealed.
Pipeline: `main.py dataset-run` → `a2s-source-lock-v2` → gate runners. Splits are
target-disjoint by homology component; within-target splits are Murcko-scaffold
disjoint; **document disjointness is now added** (see §10 F1).

**External replication (on disk).** Papyrus 05.7 (`dataset/public/papyrus_05_7`)
for cross-dataset robustness; Novartis SPD 2023 (`dataset/public/spd_2023`,
CC BY-4.0) for **multi-family** transfer. SPD is recorded as dual-cold
*underpowered* on query depth (median 14 compounds/assay), so it tests family
generalization, **not** a powered dual-cold effect — that limitation is stated in
advance, not discovered afterwards.

**External replication (download).** Davis and KIBA via the standard public
releases, for comparability with the DTA literature. Both are small and
target-shallow; they are a comparability check, not the confirmation set.

Seeds, aggregation, and compute: RTX 4060 (`D:\anaconda\envs\drug\python.exe`).
Gate T0 ran in 38 s; FBA training at ~750 parameters is minutes, not hours.

**Required baselines**, all at matched parameter count, training budget and data
access: frozen base; ligand-only; protein-only; full ligand–protein; per-target
ridge; fixed LASSO; OMP; PCA low-rank adaptation; **random-basis adaptation**;
MAML; ANIL; TRACE / scaled KRR; kernel adaptation; **pooled single head**;
**document-mean oracle** (new, from T0A); **magnitude-matched wholesale**;
**random selection at matched coverage**.

**Mandatory gates** (the user's Gates 1–6, mapped):

| Gate | Content |
|---|---|
| **G1 synthetic positive control** | planted `M`-mode world at measured σ; must recover the mode and the k-shot gain, as T0's control already does (18.2 % vs 0.9 % chance) |
| **G2 representation improvement** | `M/M_useful` must fall from 8.5 to ≤ 1.5; posterior entropy at k=5 must fall; report effective adaptation dimension and the k-curve |
| **G3 k-shot improvement** | break-even must move from the measured k ≈ 20 to k ≤ 5, with a k=5 lower bound > 0.005 **over the frozen base and over the pooled head** |
| **G4 wrong-support control** | correct support must beat wrong-target support, label permutation and random support |
| **G5 shortcut prevention** | ligand-only, protein shuffle, support permutation, query–support mismatch, random representation, matched capacity — **plus same-document evaluation and the document-mean oracle** |
| **G6 protein contribution** | protein enters only if it adds measurable incremental value; on current evidence (G3/PIRS/C0B/T0C null ×4) the design **excludes protein from the main path** and it must earn re-entry |

---

# 10. Failure modes, and the gate that fires first

1. **The ceiling collapses with small `M`.** Most likely failure. Three operators
   may not span enough of the +0.031 offset-free headroom to clear the MDE.
   *Detected by Gate F1, before any training.*
2. **Harmlessness and expressiveness are irreconcilable.** Operators forced to be
   non-harmful may converge to the pooled head, which T0 measured at −0.008.
   *Detected by G2 (`M/M_useful` falls but the ceiling falls with it).*
3. **The document confound re-enters.** Any gain must be shown on same-document
   pairs. *Detected by G5.*
4. **Selection collapses to retrieval.** If `A_θ` degenerates to nearest-support
   similarity, the stop criterion fires. *Detected by the Tanimoto < 0.35
   stratum and by the TRACE baseline matching it.*
5. **It reduces to a closed-form posterior.** If a fixed squared-error `ℓ` and a
   fixed temperature match the learned ones, the meta-learning claim is void.
   *That comparison is a required arm, not an afterthought.*
6. **Power.** ~50 probe components, MDE80 ≈ 0.005–0.010, on a prize now restated
   at +0.031. The k=5 effect must be ≥ 0.005 to be detectable at all.

## Gate F1 — WITHDRAWN as written

F1 asked whether a small operator library has a useful oracle ceiling. It is
withdrawn for two independent reasons:

1. **It decided go/no-go on `probe`**, which the ledger reserves — probe outcomes
   were consumed once by PIRS and may not drive model selection.
2. **Its premise is gone.** It presupposed that whole-head transfer is real; the
   corrected T0B measures that transfer at −0.018 [−0.044, +0.005] on
   same-document pairs.

## Gate D0 — what must run first instead

Before any mechanism is selected, the substrate must be rebuilt so that a positive
result would mean something:

> **D0a — a genuinely separated split.** Nested, deterministic, `fit`-only
> development tasks with **simultaneous** target, scaffold, **document** and
> **assay** separation. The measured baseline to beat is 91.1 % document overlap
> and 88.8 % assay overlap.
>
> **D0b — the estimand.** Score adaptation on same-document query pairs drawn from
> documents **absent from the support**, so that neither the support nor the
> evaluation can carry the offset.
>
> **D0c — how much is left.** Report the per-target chemical headroom that survives
> D0a+D0b, with the document-mean oracle as a mandatory control. This number
> replaces +0.053 and +0.031 as the programme's prize, and every power calculation
> is restated on it.
>
> **D0d — prospective utility, not heuristic bits.** Replace the retracted bit
> account with (i) an operationally defined utility per candidate adaptation and
> (ii) an *empirical* estimate of task/operator information — a decoder-error
> curve or a measured mutual-information estimate — not a Gaussian-channel proxy.

**If D0c returns a headroom whose lower bound is below the MDE, the programme's
terminal deliverable is the negative result, and no mechanism is selected.** That
is a live possibility: T0A's chemical remainder already fails to clear its
threshold under the confounded split, and D0 can only reduce it.

---

# 11. Go/no-go decision criteria

## Scope limitation: k = 1 (stated, not implied)

FBA enforces exact base recovery at k = 1, because after level removal a single
label yields zero contrasts — confirmed by T0's control recovering exactly chance
at k = 1. **FBA therefore cannot satisfy an "improvement across k = 1, 3, 5"
criterion.** Only two routes exist to a k = 1 ranking gain: label-free target
conditioning (measured null four times here — G3, PIRS, CFES-C0B, T0C), or an
explicitly revised k = 1 objective that admits the level channel into ranking
(measured rank-null by constraint C1). If the criterion is binding, this design
does not meet it and must be replaced rather than tuned.

**Proceed to implementation** only after Gate D0 (§10) returns a headroom whose
lower bound clears the MDE on a target-, scaffold-, document- and
assay-separated split, and only with a mechanism whose novelty survives §3b.

**Admit the mechanism** only if all hold on untouched probe components:
k=5 gain over the frozen base **and** over the pooled head with lower 95 % bound
> 0.005; positive on same-document pairs; correct support beats wrong-target,
permuted-label and random support; the learned `ℓ_θ`/`T_θ` beat their closed-form
counterparts; matched-capacity and random-basis controls do not match it; support
removal is an exact no-op; no degradation of the support-free path.

**Stop immediately** if the mechanism reduces to retrieval or similarity
weighting; is matched by a closed-form posterior; is matched by fixed sparse
regression; gains only in chemically local or cross-document regimes; fails to
move break-even below k = 5; or fails any shortcut control.

**If NO-GO or stop — the terminal deliverable.** The bit-budget half of revision
1's proposed negative is **retracted** (§5 of the T0 decision): it depended on an
arbitrary definition of a useful head and its apparent agreement with the measured
break-even used the same data on both sides. What remains is measured,
deterministic and basis-independent:

> On open ChEMBL affinity data, a chemistry-free oracle that knows only a
> compound's source document and that document's mean residual scores +0.061 CI,
> **exceeding a full per-target chemical response head (+0.052)**. Within-target
> Murcko-scaffold-disjoint splits on this substrate leave 91.1 % of query rows in a
> document already seen in support and 88.8 % in a seen assay, with every one of 52
> targets sharing documents across the split. Once same-document pairs are scored,
> the per-target chemical remainder does not clear a 0.005 admission threshold, and
> cross-target transfer of a fitted response head is negative
> (−0.018 [−0.044, +0.005]).
>
> Consequently, few-shot DTA results on such corpora that are evaluated without
> document- and assay-disjoint controls are not safely interpretable as
> target-specific chemical adaptation.

This is a claim about the evaluation of a field, supported by a deterministic
measurement and a reproducible audit. It does not require a mechanism to be true.

Both branches of this decision produce a publishable, reproducible, open-data
result. That is the point of running F1 first.

---

## Standing constraints

`locked` and the A2S recipient roster stay sealed; only `fit` and `probe` may be
opened. `probe` is a development role — confirmation requires freezing the
protocol and opening `locked` once. No promotion to `model/` before a complete
gate sequence passes. Every ranking result is reported beside the frozen base,
the pooled head, a magnitude-matched control, a random-selection control and the
document-mean oracle.
