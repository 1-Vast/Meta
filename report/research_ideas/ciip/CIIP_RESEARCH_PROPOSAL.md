# Corrected Interaction Identifiability Probe (CIIP)

**Document class:** research design input  
**Status:** proposal only; not preregistered, not executed, and not evidence of model performance  
**Authority:** `task.md` and `report/EVIDENCE_LEDGER.md` override this document  
**Scope:** determine whether a transferable protein-conditioned ligand-effect signal is identifiable before integrating a new interaction module into MetaSieve

## 1. Research question

The project should no longer ask whether a sufficiently expressive network can fit
protein and ligand identifiers. The relevant question is:

> After controlling for ligand transformation, local chemical context, assay
> context, target and ligand main effects, does local protein information predict
> heterogeneous ligand effects on held-out protein components?

The two load-bearing hypotheses are:

- **H1, heterogeneity:** a defensible interaction-effect surface contains variation
  beyond its estimand-specific repeat-measurement noise.
- **H2, transferability:** on unseen protein components, the correct local protein
  representation predicts that variation better than chemical-context-only,
  shuffled-protein, family-preserving-shuffle, and ligand-invariant protein
  controls.

H1 is only an admission check. Non-zero residual variance alone is not evidence of
a protein-conditioned interaction because it can also be produced by missing assay
covariates, censoring, chemical-context mismatch, or correlated measurement error.
H2 is the primary scientific gate.

## 2. Existing evidence that constrains the proposal

CIIP must extend, not repeat, Stages T, U, and V.

- Stage U found 37,945 same-panel fit observations across 243 targets and 30,463
  core-inclusive exact keys. Only 1,001 keys were rich across at least three
  targets/components.
- One target, alone in its component, contributed 29.63% of fit observations,
  exceeding the frozen 25% concentration cap.
- The designated Stage V internal repeated-key surface contained only 32 rows over
  four protein components, below the inherited 100-row evaluability floor, with
  zero internal rich exact keys.
- Meta-train and the development-validation split shared zero core-inclusive exact
  keys under the double-cold closure.
- Exact-MMP deployment coverage was 0.226, 0.362, 0.442, and 0.526 for
  k = 1, 2, 3, and 5.
- Only 88 of 42,534 same-panel MMP pairs had repeat measurements; 40 had zero
  observed range. The available repeat subset did not provide a stable,
  estimand-matched noise floor.
- Stage V therefore stopped before neural training. Its scoped conclusion was
  **not estimable on the current exact-key surface**, not biological absence.

Consequently, a new BindingDB exact-rectangle census is not a new experiment.
Any CIIP work on the current corpus must identify precisely what information is
new relative to Stage U/V.

## 3. Three distinct estimands

The following surfaces must not be pooled or described as equivalent.

### A. Exact chemical rectangle

For proteins `P1`, `P2` and ligands `A`, `B` measured under compatible assay
conditions:

```text
D = [y(P1, B) - y(P1, A)] - [y(P2, B) - y(P2, A)]
```

This is the cleanest interaction estimand because protein and ligand main effects
cancel algebraically. On the current BindingDB double-cold surface, however, Stage
U/V already established that confirmatory support is insufficient.

### B. Same-ligand protein contrast

For a ligand measured against a controlled WT/mutant or matched protein pair:

```text
C(L; P1, P2) = y(P1, L) - y(P2, L)
```

This is not automatically a ligand-dependent interaction. A ligand-invariant
protein or mutation shift can predict it. Interaction evidence requires multiple
ligands per protein pair and rejection of a ligand-invariant shift control.

This surface is best suited to a same-platform dense WT/variant panel and is the
preferred real positive control.

### C. Contextual matched molecular-pair pseudo-rectangle

Ligand transformations are grouped by increasingly permissive context:

- `C0`: exact transformation and exact shared core;
- `C1`: transformation plus attachment environment;
- `C2`: transformation plus bounded local molecular neighborhood;
- `C3`: whole-core or series context.

`C1` may be used as a primary exploratory screen. `C2` and `C3` reintroduce
chemical-context nuisance and must remain exploratory until reproduced on a
confirmatory same-platform surface. `C3` additionally risks scaffold or series
memorization and requires a scaffold-cold evaluation.

The ladder is not assumed to be monotonic in scientific quality. Too little
context can make a transformation non-exchangeable across scaffolds, while too
much context can memorize a chemical series. CIIP-0B must therefore estimate a
Pareto curve over context levels:

```text
context residual variance V_C = Var(D - E[D | C])
effective cross-component support N_eff(C)
```

The preferred context is the smallest level that materially reduces residual
heterogeneity without collapsing cross-component support. `C1` is not privileged
in advance; it is only the first candidate to audit.

## 4. Assay compatibility contract

`assay_block_id` must be deterministic and auditable. A document identifier alone
does not establish exchangeability. The strongest available block should include:

```text
source/document
+ endpoint type
+ assay format
+ measurement orientation and units
+ target construct
+ substrate/cofactor and concentration context
+ protocol signature
```

Missing fields must be counted and reported. Confirmatory analysis may not use a
learned text embedding to infer that two protocols are equivalent. If construct or
protocol equivalence cannot be established, the comparison is downgraded or
excluded rather than silently repaired.

## 5. Transformation representation

The transformation representation should separate even context from odd signed
change:

```text
u(A -> B) = [u_even(shared context), u_odd(signed edit)]
u_even(B -> A) = u_even(A -> B)
u_odd(B -> A)  = -u_odd(A -> B)
g(P, B -> A)   = -g(P, A -> B)
```

Antisymmetry belongs to the predicted transformation effect and its signed edit
subspace, not to the entire hidden representation. Attachment, core, and local
environment information are symmetric context and should not be forced to change
sign.

For unseen transformations, a categorical transformation lookup is invalid. The
model must use a continuous chemical-context encoder learned only from training
components. Seen- and unseen-transformation results are reported separately.

## 6. Protein representation and router

The first admissible protein route is deliberately low capacity:

- frozen residue-level protein features;
- preserved residue order and masks;
- local or pocket-aware aggregation when the mapping is governed;
- bounded interaction head with no free target-ID embedding;
- no test-time gradient, closed-form solver, ridge estimator, or query-label input.

Frozen ESM features and low capacity do not by themselves prevent target-key
memorization. Required controls include:

- matched random-protein input;
- family-preserving protein shuffle;
- local-window versus matched random-window comparison;
- protein-length, document-frequency, construct and mutation-status matching;
- ligand-invariant protein or mutation-effect baseline;
- attention/localization audit against governed pocket or mutation positions.

KLIFS-based claims are kinase-specific. They must not be generalized to the full
BindingDB target population.

## 7. Cross-fitted nuisance model

A nuisance model may estimate transformation and context main effects, but it must
be cross-fitted after the protein-component split:

```text
b_hat = nuisance_model(transformation, chemical_context, assay_context)
D_res = D - b_hat
```

The transformation vocabulary, normalization statistics, hyperparameters and
nuisance fit may use meta-train components only. Both raw `D` and residualized
`D_res` must be reported because residualization can remove real interaction when
transformation and protein component are confounded.

The audit must report:

- variance absorbed by the nuisance model;
- transformation coverage by protein component;
- seen- and unseen-transformation performance;
- raw and residualized effect metrics;
- whether the nuisance model used any evaluation-component labels.

## 8. Required model and control arms

The exact arm set depends on the estimand, but the confirmatory minimum is:

1. zero interaction;
2. chemical context only;
3. free target or family identifier, diagnostic only and never deployable;
4. correct local protein representation;
5. matched shuffled protein;
6. family-preserving shuffled protein;
7. ligand-invariant protein/mutation shift;
8. optional oracle representation, diagnostic only.

The correct arm must improve its own prediction of the true estimand. A contrast
cannot pass merely because training makes a corrupted arm catastrophically worse.
For strict component-held-out evaluation, a shuffle that cannot be constructed
without violating matching strata is recorded as not evaluable rather than relaxed.

For a dense WT/variant panel, the positive-control question is split into two
claims:

- **CIIP-1A, representation positive control:** within a parent protein, does the
  local protein representation explain ligand-dependent WT/variant deviations
  beyond a ligand-invariant mutation shift?
- **CIIP-1B, transfer positive control:** after holding out complete parents or
  pocket groups, does the same representation transfer to an unseen parent?

CIIP-1A passing does not imply CIIP-1B. It demonstrates representational capacity,
not cross-parent transfer.

## 9. Noise, censoring, and effective support

### Noise

The noise reference must match the exact estimand and assay surface. A cell-level
or cross-condition variance cannot be subtracted from a four-measurement rectangle
without validating its covariance structure. If no estimand-matched repeat surface
exists, H1 remains unresolved.

### Censoring

For intervals `y in [l, u]`, interval subtraction and addition must preserve bound
orientation through every contrast. Report:

- interval-likelihood performance on all admissible observations;
- sign-identifiable coverage;
- determinate-only sensitivity;
- floor-imputation as an explicit negative control.

The sign-identifiable subset cannot be the only headline result because selecting
it changes the evaluation population.

### Effective support

No single inverse-Herfindahl number is sufficient. Report unique counts and
concentration across:

- protein components and parent proteins;
- transformations;
- scaffolds;
- documents and assay blocks;
- connected components of the protein-transformation bipartite graph.

Primary uncertainty uses protein-component cluster bootstrap. Transformation-block
permutation and parent/scaffold sensitivity analyses are secondary. Row bootstrap
is not an admissible substitute.

### Nuisance specification sensitivity

The primary biological estimand remains the raw contrast `D`. Residualized contrasts
are attribution diagnostics. At minimum, report three nuisance specifications:

```text
N0: no residualization
N1: transformation + chemical context
N2: transformation + chemical context + assay context
```

Protein incremental value must be stable across these specifications. A gain that
appears only under `N2` is not automatically biological; protein features and assay
residualization may be absorbing the same correlated confound. Report variance
absorbed, raw/residual effect sizes, and the direction of every primary contrast.

## 10. Corrected staged programme

### CIIP-0A: historical equivalence audit

Reuse Stage U/V artifacts and map every proposed exact-rectangle statistic to an
existing result. Do not re-run an equivalent census. Output only the gaps not
already answered.

### CIIP-0B: new contextual census

Build `C1` and, secondarily, `C2` transformation contexts after split creation.
Measure graph connectivity, support concentration, assay compatibility, censoring,
seen/unseen transformation coverage, and power. Do not train a neural model.

### CIIP-0C: estimand-specific information audit

Estimate repeat noise only where the repeat surface matches the intended contrast.
Use planted-signal simulations on the observed graph to determine detectable effect
sizes under the intended split and cluster structure. A synthetic pass qualifies
the instrument, not the biological hypothesis.

### CIIP-1A: same-parent representation positive control

Use a governed WT/variant panel with multiple identical ligands per protein pair.
Estimate the ligand-invariant mutation shift and test whether the centered
ligand-dependent deviation is predictable. This tests local protein information
capacity and is not yet a transfer claim.

### CIIP-1B: held-out-parent transfer positive control

Hold out complete parent proteins or governed pocket groups. Test whether the local
protein representation improves ligand-dependent contrasts beyond chemical-only,
global-protein, family-shuffle, and ligand-invariant controls. Preserve endpoint
semantics; functional inhibition is not renamed as Ki, Kd, pK, or DTA.

### CIIP-2: low-capacity transferable probe

Only after CIIP-1 passes, train the bounded local protein-by-transformation probe on
held-out parent or pocket components. Pre-register primary contrasts and correct
for multiple comparisons.

### CIIP-3: MetaSieve integration

Only a representation that passes the real positive control and held-out-component
gate may enter `model/` and production `scripts/`. Integration is an ablation with
the verified interaction branch on/off; it does not replace the practical few-shot
baseline or the fixed ligand-similarity comparator.

## 11. Admission and stopping rules

### GO

Advance only when all of the following hold:

- the surface has adequate component-level support according to a frozen power
  analysis, not an arbitrary raw row count;
- assay and construct compatibility are auditable;
- the noise model matches the estimand or the effect is identifiable without an
  unsupported noise subtraction;
- the correct protein arm beats chemical-only and the pre-registered protein
  controls on unseen protein components;
- the primary component-bootstrap lower bound is above zero;
- the lower bound also exceeds a pre-specified minimum meaningful effect
  `delta_min`, derived from the estimand-matched noise audit and planted-signal
  power curve rather than selected after seeing results;
- the gain survives scaffold, document, censoring, and transformation-novelty
  sensitivity analyses;
- no query label, evaluation-component statistic, or target identifier enters the
  deployed input path.

### WEAK-GO

A contextual `C1/C2` result may authorize an external replication, but not a claim
of exact protein-conditioned interaction and not production integration.

### NO-GO

Stop the local route when exact support is already exhausted, contextual support is
too concentrated, assay/noise equivalence is unresolved, correct and shuffled
protein are indistinguishable, or the apparent gain is explained by a
ligand-invariant protein shift.

Statistical detectability alone is not sufficient for GO. A positive but negligible
effect below `delta_min` is recorded as a weak result and does not authorize model
integration.

## 12. Relationship to the performance programme

CIIP is an attribution gate, not the sole model-admission gate. It must not block
the practical few-shot programme:

- practical target-cold evaluation may use realistic support sizes and related
  chemical series;
- ligand-only Tanimoto, ordinary fine-tuning, and meta-learning remain valid
  comparators;
- strict protein-component-cold plus scaffold-cold evaluation remains the
  fundamental stress test;
- failure of CIIP forbids attributing performance to transferable
  protein-conditioned interaction, but does not imply that a ligand-driven
  few-shot model has no practical value.

The final DTA performance claim still requires independent training and evaluation
on each dataset under its native label semantics, using MSE/RMSE, CI and rank
metrics with target/component-level uncertainty.

## 13. Current feasibility assessment

| Candidate route | Feasibility | Current status |
|---|---:|---|
| BindingDB exact rectangles | 2/10 | Existing Stage U/V evidence is effectively a local NO-GO |
| BindingDB contextual `C1/C2` | 5/10 | New exploratory census is defensible after governance reopening |
| Same-platform WT/variant panel | 7/10 | Preferred real positive-control surface |
| Immediate neural CIIP training | 1/10 | Not authorized |
| Integration after a positive control | Conditional | Requires all CIIP gates |

## 14. Highest-priority real positive-control candidate

The 2026 KiRHub kinase-inhibitor profiling study is a strong candidate for CIIP-1A
and CIIP-1B. The primary article reports 92 clinical kinase inhibitors profiled
against 758 kinases, including 409 wild-type kinases and 349 mutant kinases or
fusions; 86 compounds were FDA approved. The full panel used a HotSpot radiometric
functional kinase assay at 1 micromolar with K_m ATP and duplicate measurements.
The paper reports approximately 290,000 kinase-drug measurements and within-dataset
replicate agreement of R2 = 0.99. Primary source:

<https://www.nature.com/articles/s41587-026-03090-8>

This is well matched to CIIP-1A because identical inhibitors are measured across
WT and variants in one profiling campaign. The estimand remains functional
inhibition/selectivity, not Kd, Ki, pK, or DTA. Assay optimization varies by kinase,
including substrate, cofactor, and ATP-related conditions, so those fields must be
imported into the assay-block census rather than assumed identical merely because
the platform name is shared.

Before training, audit:

```text
usable WT-variant pairs
identical-ligand count per pair
point mutations versus fusions
same substrate and construct coverage
duplicate and saturation rates
mutation-specific centered-effect variance
parent and pocket-group connectivity
held-out-parent support
```

The dense panel is a representation positive control by default. It becomes a
transfer positive control only after complete parent or pocket-group holdouts are
feasible. Its result must not be used to claim that BindingDB Ki prediction is
solved.

Duong-Ly 2016 remains a historical replication surface, not an automatically
exchangeable training set. Cross-study WT values and protocol differences must be
audited before any paired replication is reported.

## 15. Instrument qualification simulations

CIIP-0C must not use only an iid planted interaction. On the observed graph, inject
at least three fully synthetic regimes with independent train/validation/test
components:

1. **Protein main-effect regime:** `D = u_p + noise`. The interaction probe must
   not pass.
2. **Family-shortcut regime:** `D = u_family(p,tau) + noise`. A
   family-preserving shuffle must expose the shortcut; the probe must not call it
   transferable local interaction.
3. **Local interaction regime:** `D = z_p^T M z_tau + noise`, with `z_p` derived
   from a governed local protein representation. Only this regime may authorize a
   protein-conditioned interaction pass.

Report interaction recovery, dead-zone sign accuracy, scale recovery, false-positive
rate on regimes 1/2, and power over a range of effect sizes. A planted pass
qualifies the instrument only; it does not establish biology.

## 16. Immediate recommendation

Do not increase model complexity. First close the active synthetic qualification
chain with a terminal summary, then explicitly authorize a new governance cycle if
CIIP is to proceed. The first actionable CIIP task is a CPU-only equivalence and
contextual census that reuses Stage U/V, quantifies only genuinely new `C1/C2`
support, and decides whether an external same-platform positive control is required.

No code should be promoted into `model/` or production `scripts/` on the basis of
this proposal alone.

The most informative next action is a CPU-first KiRHub compatibility census followed
by a frozen CIIP-1A/1B preregistration. In parallel, CIIP-0A should map already-
answered Stage U/V questions and CIIP-0B may audit only genuinely new contextual
support. If KiRHub is unavailable or fails its assay/coverage census, retain the
current bounded conclusion rather than substitute a weaker cross-platform panel.

## 17. Evidence and literature notes

- BindingDB is a document-level aggregation of heterogeneous affinity and activity
  records. Its provenance fields do not guarantee arbitrary cross-target
  exchangeability: <https://pmc.ncbi.nlm.nih.gov/articles/PMC11701568/>.
- ActFound motivates within-assay pairwise relative activity because absolute
  measurements across assays are not directly comparable; this does not itself
  prove protein-conditioned transfer:
  <https://www.nature.com/articles/s42256-024-00876-w>.
- Matched-molecular-pair uncertainty depends on pair count, experimental source,
  measurement uncertainty and physical-effect variability:
  <https://pubs.acs.org/doi/10.1021/jm500317a>.

## 18. What CIIP can and cannot complete

CIIP must not be presented as the complete MetaSieve solution. The final task is
the prediction of an absolute endpoint for an unseen target with a variable-size
support set:

```text
zero-shot:  y_hat_q = F(P*, L_q)
few-shot:   y_hat_q = F(P*, L_q, {(L_i, y_i)}[i=1..k])
```

CIIP primarily tests whether a source-learned protein-conditioned ligand effect
exists and transfers. It does not by itself solve three separate problems:

1. **Absolute level:** for any additive target and ligand functions,
   `y_tilde(P,L) = y(P,L) + a(P) + b(L)`, an exact double difference is unchanged.
   Therefore a perfect interaction estimator cannot recover the absolute affinity
   gauge. A DTA model still needs a level/baseline path.
2. **Zero-shot utility:** a statistically real interaction may be too small or too
   sparse to improve unseen-target MSE, concordance, or ranking.
3. **Few-shot adaptation:** CIIP does not specify how one to five support labels
   modify a query-specific interaction state.

The deployable prediction therefore needs a factorized form with an explicit
identifiability convention:

```text
y_hat(P,L,S) = level(P,L,S) + shape(P,L,S)
```

`level` handles absolute calibration. `shape` is centered within the target or
support-defined reference set so it cannot absorb a constant target offset. A
simple initial constraint is a zero-mean shape over the training/query panel; a
shrinkage version is preferable when the panel is small. The constraint is an
engineering identifiability device, not a claim that biological effects have zero
mean.

## 19. Practical completion path

The shortest implementable route is a bridge, not direct CIIP-loss integration.

### Bridge A: frozen interaction utility

After a real positive-control representation is available, freeze it and compare
capacity-matched models on unseen target components:

```text
B0: ligand-only baseline
B1: ligand + global protein baseline
B2: ligand + verified local interaction representation
```

Primary utility metrics are within-target centered MSE/RMSE, concordance index and
Spearman correlation. Absolute MSE remains mandatory but is secondary for judging
the interaction branch. Also compute an oracle interaction ceiling using the
measured interaction effect. This separates three cases:

- oracle and learned both help: integration is justified;
- oracle helps but learned does not: the data contain useful signal, but the
  representation or learner is inadequate;
- oracle does not help: the interaction is not the current DTA bottleneck.

### Bridge B: support adaptation utility

Use nested supports `S1 subset S2 subset S3 subset S5` on the same query panels and
compare:

```text
Z0: zero-shot
Z1: level-only support calibration
Z2: level + fixed Morgan/Tanimoto transport
Z3: level + verified protein-conditioned transport
```

The central quantity is not merely `MSE(k) - MSE(0)`. It is the incremental shape
utility of `Z3` over `Z1` and `Z2`, with label-shuffled, wrong-target and
structure-only support controls. This preserves the project's strongest existing
ligand-side comparator instead of assuming the protein branch replaces it.

For `k=1`, the model should primarily update level and confidence. A single support
point does not identify a new within-target ligand slope. A large query-specific
shape correction at `k=1` requires a source-learned prior and must be compared
against a no-shape level-only control. For `k>=2`, support residual differences
provide direct SAR evidence and can update a centered shape field.

## 20. Realistic interpretation of KiRHub

KiRHub is useful because its common inhibitor panel and WT/variant structure can
test whether local protein information has ligand-dependent functional effects. It
does not make BindingDB Ki prediction easier by itself, and it does not provide a
cross-domain guarantee outside kinases. Its most practical use is:

```text
KiRHub: validate representation and transfer capacity
BindingDB-Ki: validate native DTA level and few-shot utility
```

The two label spaces must remain separate. A successful KiRHub probe authorizes a
small verified protein representation experiment; it does not justify merging the
datasets or pretraining on one and reporting the result as a single Ki benchmark.

## 21. Feasibility without overengineering

The project is feasible if the claims are staged according to what the data can
support:

| Claim | Current feasibility | Required evidence |
|---|---:|---|
| Protein-local information exists in a controlled panel | High | KiRHub CIIP-1A |
| It transfers to unseen kinase parents | Moderate | KiRHub CIIP-1B |
| It improves BindingDB zero-shot DTA | Unresolved | frozen Bridge A |
| It improves few-shot SAR beyond level/Tanimoto | Moderate for k>=2; low for k=1 | Bridge B |
| Strict double-cold MSE <= 1.0 | High risk | independent multi-seed native-DTA evaluation |

This path is practical because CIIP-1A/1B can use a low-capacity probe and CPU/GPU
smoke tests, while Bridge A/B reuse the existing MetaSieve evaluation and support
contracts. No Cartesian module, diffusion model, closed-form solver, or large
meta-adapter is needed before a measurable utility signal appears.

The correct success statement is therefore conditional:

> CIIP can establish a transferable interaction prerequisite. Only Bridge A and
> Bridge B can establish that the prerequisite completes the cold-target DTA task.

If CIIP passes but both bridges fail, retain the interaction result as a valid
scientific finding but stop treating it as the main performance route. If the
bridges pass without CIIP, retain the practical ligand-driven model but do not
claim protein-conditioned mechanism. Only the conjunction of signal transfer,
zero-shot utility and few-shot incremental utility warrants integrating the branch
into production MetaSieve.

## 22. Deployable Potential Bridge

The main correction introduced by this revision is functional, not another
evaluation gate. A free pairwise transformation predictor is not automatically a
single-ligand DTA scorer. It may predict `g(P,A,B)` accurately while violating
cycle consistency and while having no scalar function whose differences equal the
predictions.

The CIIP supervision should therefore be attached to a deployable scalar
protein-conditioned SAR potential:

```text
s_theta(P, L)                  # single-ligand shape score
g_theta(P, A -> B) = s_theta(P, B) - s_theta(P, A)
D_hat = [s(P_a,B)-s(P_a,A)] - [s(P_b,B)-s(P_b,A)]
```

This parameterization gives, by construction:

- antisymmetry under `A <-> B`;
- zero effect for an identity transformation;
- cycle consistency for `A -> B -> C -> A`;
- a score that can be evaluated for one query ligand at deployment.

The potential is a hypothesis to be tested, not a claimed result. The first
diagnostic should compare a low-capacity potential against a free pairwise head on
synthetic and real positive-control surfaces. If the potential cannot recover a
known local interaction while the free pairwise head can, the function-level
restriction is too strong and must not be forced into production.

### Absolute endpoint decomposition

The proposed DTA form is:

```text
y_hat_0(P,L) = b_L(L) + b_P(P) + s_theta(P,L)
y_hat_k(P,L_q) = y_hat_0(P,L_q) + b_S + delta_s(P,L_q,S)
```

Here `b_L` is a ligand prior, `b_P` is a protein-derived level prior, `s_theta`
is the only branch allowed to carry protein-by-ligand interaction, and `b_S` is a
support-level calibration. `b_L` and `b_P` cancel from exact double differences.
The level/shape split must be implemented as an information-path constraint, not
only as a weak regularizer; otherwise the absolute loss can make the shape branch
relearn target level.

The shape score needs a gauge convention. A query-panel mean is transductive and
must not be used to define a single-query score. Prefer a fixed training reference
ligand measure or reference basis, with support-specific centering applied only to
the few-shot residual update.

### Local interaction potential versus global bilinear probe

The compact probe

```text
s(P,L) = alpha(P)^T psi(L)
```

is useful because it is integrable and cheap. It is not the default production
backbone: Stage S showed that early global protein compression can become a target
or family key. The preferred successor, only after positive-control evidence, is a
bounded local field:

```text
H_P = residue-level protein states
H_L = atom or functional-group ligand states
e_ij = phi(H_P[i], H_L[j])
a_ij = sparse, masked ligand-conditioned residue weighting
c_j  = phi_I(H_L[j], sum_i a_ij H_P[i], H_L[j] * sum_i a_ij H_P[i])
s(P,L) = scalar_pool(c_j)
```

Without a common complex frame, `a_ij` means learned ligand-conditioned residue
relevance, not contact probability or an atomic 3D interaction. The low-dimensional
space should be applied after local interaction formation (`z(P,L)`), not by
compressing the entire protein into `alpha(P)` before interaction.

### Few-shot update on the same field

For support examples, compute:

```text
r_i = y_i - y_hat_0(P*, L_i)
b_S = shrink_k(mean_i(r_i))
r_tilde_i = r_i - mean_i(r_i)
```

The shape update must be label-directed:

```text
delta_s(P*, L_q, S)
  = eta_k * sum_i K(q_q, q_i) * r_tilde_i
```

or an equivalent low-rank update derived from the same local field. It must satisfy
`r_tilde_i = 0 for all i => delta_s = 0`. Consequently, `k=1` has no centered
shape evidence and should initially be level-only; `k>=2` is the first setting in
which support residual differences can identify target-specific SAR. Fixed
Morgan/Tanimoto transport remains a comparator and prior, not something to discard.

### Training objectives share one field

The first trainable prototype should use one model and three supervised views:

```text
L_abs   = loss(b_L + b_P + s, y)
L_delta = loss((b_L(B)+s(P,B)) - (b_L(A)+s(P,A)), y_B-y_A)
L_D     = loss(Delta_P Delta_L s, D)
```

Exact rectangles receive higher confidence than contextual pseudo-rectangles.
Contextual pairs use cross-fitted nuisance residuals and lower weights; they must
not force `s` to absorb chemical-context or assay mismatch. After CIIP
initialization, retain a small `L_D` term or equivalent representation retention
term during native DTA training so the absolute objective cannot erase the
interaction property. This remains ordinary joint training, not a closed-form
solver, MAML inner loop, or deployment-time optimization.

## 23. Revised execution order

The proposal is now ordered by implementation risk and information value:

1. Finish and archive the currently running Q2d chain. Do not start a new
   synthetic successor merely to avoid its terminal decision.
2. Run a CPU-only KiRHub compatibility census. Do not train until usable pairs,
   endpoint semantics, saturation, assay fields, and parent/pocket connectivity
   are known.
3. Run CIIP-1A on a low-capacity potential and its free-pairwise diagnostic, then
   CIIP-1B only if 1A is interpretable.
4. Run the frozen Potential Bridge on native BindingDB-Ki with B0/B1/B2 and the
   oracle utility ceiling. This is the first test of DTA usefulness.
5. Run the support bridge with level-only, fixed Tanimoto, and potential-conditioned
   transport for `k=1,2,3,5`.
6. Only if the potential is useful and controls are clean, implement the local
   interaction field in `model/` and `scripts/` as a single-variable successor.

No result from CIIP-1A alone authorizes production integration. No result from a
functional inhibition panel is reported as BindingDB affinity performance.

## 24. Current execution status (2026-08-19)

The current repository evidence changes the data order, but not the potential
bridge itself:

- **Q2d-1e:** terminal `GATE FAIL` has been adjudicated. The authorized
  span-parameterization diagnostic is still running; no further synthetic
  successor is allowed. This is a synthetic learner diagnosis and makes no
  biological claim.
- **KiRHub census:** `DATA BLOCKER`. The raw profiling tables and a verified
  download URL are not locally available, so CIIP-1A/1B training on KiRHub is not
  authorized. Published headline counts are not a substitute for a local usable
  pair census.
- **Davis census:** `INSUFFICIENT ALONE`. It has 67 usable WT/variant pairs and a
  median of 33 common ligands per pair, but only seven parents with at least two
  mutant pairs. Its high missingness and unresolved NA semantics also limit a
  transfer claim. It is useful for a representation probe or replication audit,
  not as the sole CIIP-1B surface.
- **P-line:** practical few-shot training is parallel and must remain separate
  from the mechanism attribution conclusion.

The immediate data action is therefore a combined, read-only census of the local
Davis, Anastassiadis, and Duong-Ly panels. It must reuse the frozen ten-item
compatibility checklist, report intersection and union support, and distinguish
same-parent representation capacity from held-out-parent transfer. Only a
combined surface that actually supplies enough valid parents, identical ligands,
assay-compatible observations, and uncensored/interval-aware labels may authorize
a new CIIP preregistration.

If the combined local surface remains insufficient, the correct outcome is
`UNRESOLVED/DATA-BLOCKED`, not a weaker cross-platform merge and not a return to
BindingDB exact-MMP training. Acquiring KiRHub or another governed same-platform
panel is an external data prerequisite, not a model-code fix.

## 25. Updated status after the combined local census

The next census materially changes the data-admission decision:

- The combined local census is **SUFFICIENT for a new CIIP-1 preregistration**.
- Duong-Ly alone supplies the usable same-platform surface: 70 eligible single-
  mutation WT/variant pairs, a median of 183 identical ligands per pair, and 12
  parent proteins that satisfy the held-out-parent requirement. Six multi-mutant
  rows are excluded from the single-mutation estimand and retained as a separate
  count. The endpoint remains percent inhibition; it must never be relabeled as
  pK, Ki, Kd, or affinity.
- Davis remains useful for same-parent representation evidence but does not meet
  the held-out-parent transfer requirement by itself.
- Anastassiadis is a cross-endpoint historical replication surface only. Kd and
  percent-inhibition observations must not be merged into one CIIP estimand.
- KiRHub remains a blocked external option, but it is no longer a prerequisite for
  the first real CIIP experiment.

This authorizes the following practical order after Q2d terminal adjudication:

1. Freeze a Duong-Ly CIIP-1A/1B preregistration with the exact endpoint, mutation
   eligibility, parent holdouts, saturation/overflow policy, and all controls.
2. Run structural tests and a single-seed screening pass before spending the full
   multi-seed budget.
3. If the screening pass is interpretable, run the pre-registered multi-seed
   CIIP-1A and CIIP-1B probe.
4. Only after a positive control, run the Potential Bridge on native BindingDB-Ki.

The first practical few-shot baseline also changes the interpretation of the next
stage. Ordinary ESM plus bilinear fine-tuning is worse than fixed Morgan/Tanimoto
transport on the current k=5/10/20/40 tests, with significant paired degradation at
k=5 and k=10. This is not evidence that protein information is absent; it is a
strong utility baseline that any protein-conditioned potential must beat. The
potential bridge must therefore include:

```text
ordinary fine-tuning
fixed Morgan/Tanimoto
potential-conditioned transport
```

and must report incremental gains over Tanimoto, not only gains over ligand-only.
The practical line remains separate from CIIP attribution and should continue with
the frozen arm-4--7 screening plan after the diagnostic releases the GPU.

## 26. Unified scientific object for the next cycle

The next cycle must treat CIIP, Potential Bridge, and few-shot adaptation as
different observations of one deployable function, not as three independent
modules. Define a scalar field:

```text
s_theta(P,L)
```

and derive every relative target from its finite differences:

```text
Delta_P s(L)       = s(P_variant,L) - s(P_WT,L)
Delta_L s(P)       = s(P,L2) - s(P,L1)
Delta_P Delta_L s  = [s(Pv,L2)-s(Pv,L1)]
                     - [s(Pw,L2)-s(Pw,L1)]
```

The same `s_theta` must be used for Duong-Ly CIIP-1A/1B, the BindingDB
zero-shot bridge, and few-shot shape correction. This removes the weak bridge in
which one pairwise model is validated and a different DTA head is later trained.

For Duong-Ly, define the observed protein contrast:

```text
d_vl = y(P_variant,L) - y(P_WT,L)
c_vl = d_vl - mean_L(d_vl)
```

The model predicts `s(P_variant,L)-s(P_WT,L)` and applies the same centering. The
centered target rejects a ligand-invariant mutation shift. CIIP-1A tests this
within parent; CIIP-1B changes only the split by holding out complete parents.

The potential must be compared with a free protein-pair predictor and a free
ligand-pair predictor. The purpose is not to assume the potential is always more
accurate, but to measure the information lost when a pairwise signal is constrained
to be deployable as a scalar field. A large free-pairwise advantage means the
potential representation is inadequate or the input omits important context; it
does not prove that biological effects are non-integrable.

## 27. Architecture and training scope

The core model innovation is one object only:

> **Integrable Local Protein-Conditioned SAR Potential**

The training/adaptation innovation is one object only:

> **Centered Evidence Transport on the Potential Field**

The global bilinear form `alpha(P)^T psi(L)` is retained as a low-capacity
diagnostic. It is not the production backbone because Stage S showed that early
global protein compression can become a target or family key. A production
candidate, if authorized, retains residue-level protein states and ligand
atom/functional-group states until after local interaction coordinates are formed.

The native DTA decomposition is:

```text
y_hat_q^(k) = b_L(L_q) + b_P(P) + s_theta(P,L_q)
              + lambda_k * mean(r_i)
              + eta_k * sum_i [K_Tanimoto(L_q,L_i)
                               + gamma*K_int(P,L_q,L_i)] * r_tilde_i
```

with:

```text
r_i       = y_i - [b_L(L_i)+b_P(P)+s_theta(P,L_i)]
r_tilde_i = r_i - mean(r_i)
```

`K_int` is an incremental correction to the fixed Tanimoto prior, never a
replacement for it. The first implementation must enforce
`r_tilde = 0 => delta_s = 0`; thus k=1 is level-only by default and k>=2 is the
first setting with direct centered SAR evidence.

The ligand-difference loss must not force `s_theta` to absorb the complete ligand
main effect. Either include the ligand prior difference explicitly in the target
or use a cross-fitted ligand/context nuisance before applying the potential loss.
Exact rectangles get higher confidence than contextual pseudo-rectangles.

## 28. Stage 1 single-seed screening adjudication (2026-08-19)

The first Duong-Ly CIIP-1A screening completed normally under the frozen
preregistration. It did not authorize the three-seed run or CIIP-1B. The
screening result is a model/estimand diagnostic, not a biological falsification.

The frozen adjudicator returned:

```text
unified_local Spearman:       -0.0409
unified_local sign accuracy:   0.4869
vs family-shuffle lower bound: -0.0623
vs ligand-only lower bound:    -0.0348
free-pairwise Spearman:        0.2605
free-pairwise - unified gap:   0.3014
```

The four screened arms were matched in rows, split, seed, optimizer, and budget.
The data and structure contracts were green (17 tests), and the run preserved
the raw `% inhibition` endpoint and centered mutation target. Thus there is no
evidence that the screening failure was caused by an episode-contract or
identity/antisymmetry bug.

The failure has four important qualifications.

### 28.1 The tested representation is narrower than the proposed local route

`unified_local` in this screening uses the frozen KLIFS pocket one-hot matrix
(`97 x 1700`) and ECFP4 ligand features. It does not use the per-position local
ESM representation that previously passed the representation-capability probe.
Therefore the result closes the hypothesis

```text
KLIFS one-hot + low-capacity integrable potential
```

but does not close the broader hypothesis

```text
validated local protein representation + integrable potential
```

This distinction must be retained in all summaries.

### 28.2 The free-pairwise gap is a representation/constraint diagnostic

The free pairwise arm reached a higher mean Spearman than the unified potential,
while the unified arm remained near chance. This is consistent with either an
overly restrictive scalar potential, inadequate protein/ligand coordinates, or
optimization interference. It does not prove a biological interaction signal,
because only three of thirteen test pairs had finite Spearman for both the
unified and free-pairwise arms. The gap is therefore a hypothesis for diagnosis,
not a promotion criterion.

The next audit must report the number of nonconstant prediction pairs before
aggregating rank metrics. Constant predictions must be represented explicitly,
not silently converted to zero or omitted without a denominator.

The free-pairwise diagnostic has one expected zero-gradient parameter: its final
output bias is subtracted from the same network under reversed inputs, so that
bias cancels identically in the antisymmetric output. A `grad_cov=false` entry
for this bias is therefore an identifiability property, not evidence that the
whole diagnostic branch is dead. Gradient coverage reports must distinguish
structurally cancelling parameters from trainable parameters that fail to
receive useful gradients.

### 28.3 The centered contrast is informative but statistically sparse

The 13 pair-level test surface is clustered by parent and contains many nearly
constant mutation effects. `ligand_only` is exactly the zero centered baseline,
so its Spearman is undefined; the current adjudicator maps undefined values to
zero only inside the bootstrap comparison. This is conservative but is not a
fully informative ranking comparison. Sign accuracy, centered MSE, explained
variance, and the count of informative pairs must be reported alongside
Spearman. A later decision must not be driven by a single aggregate Spearman
number.

### 28.4 Absolute and centered objectives may compete

The frozen training loss combines centered contrast loss with raw absolute
`% inhibition` loss. The endpoint spans approximately `-12.5..191.3`, whereas
the centered mutation target is on a different and much smaller scale. The
absolute term can therefore dominate optimization or encourage the nuisance
heads to explain the signal while leaving `s_theta` nearly constant. The
screening result does not establish that this happened, but the negative/near
zero potential scale and frequent constant predictions make it a required
diagnostic.

The audit must inspect per-loss gradient norms and potential variance before any
hyperparameter change. It must not tune the loss after seeing the screening
result and call the tuned run a continuation of the frozen gate.

The current screening also exposes a scale problem: raw `% inhibition` reaches
approximately 191 while the centered contrast is bounded by the paired panel
variation. Before changing the loss, the audit must report the numerical scale
and gradient norm of `L_contrast` and `L_abs` separately. This is a diagnosis of
objective competition, not permission to silently rescale labels or alter the
frozen gate. Any rescaled or centered-only successor must receive a new
preregistration and be labeled a successor diagnostic.

## 29. Authorized follow-up after the screening failure

The immediate follow-up is a read-only qualification report, not a new model
family. It must answer:

1. Are constant predictions caused by the representation, the potential
   parameterization, the absolute-loss scale, or an evaluation bug?
2. Does the free-pairwise advantage survive when finite-pair denominators and
   parent-cluster bootstrap are reported transparently?
3. Does the previously qualified local ESM representation alter potential
   identifiability, without changing split, labels, seed policy, or gates?
4. Does a centered-only diagnostic (no production promotion) separate objective
   competition from representation failure?

The follow-up may run CPU audits and a preregistered successor diagnostic only
after the current screening artifacts are archived. It may not enter CIIP-1B,
BindingDB Potential Bridge, or production integration until a new, frozen
decision is made. Any successor must preserve the original Stage 1 result and
state explicitly which factor it changes.

The current scientific status is therefore:

```text
Stage CIIP-1A screening: FAIL for the tested one-hot potential arm
Biological protein-conditioned signal: UNRESOLVED
CIIP-1B: not authorized
BindingDB Potential Bridge: not authorized
Production model changes: not authorized
```

## 30. Collapse-first diagnostic plan

The next action is a read-only failure attribution, not another model trial. The
current failure must be separated into four possibilities:

```text
input information -> potential coordinates -> optimization gradients
                         plus independent evaluation audit
```

For the current bilinear potential,

```text
s(P,L) = alpha(P)^T psi(L)
Delta_P s_v(L) = [alpha(Pv)-alpha(Pw)]^T psi(L)
Var_L(Delta_P s_v) = Delta_alpha_v^T Cov_L(psi) Delta_alpha_v
```

This identity gives a direct, training-free decomposition of output collapse:

1. `Delta_alpha` is near zero: the protein branch does not encode WT/variant
   differences or the representation is not reaching the potential;
2. `Cov(psi)` is near zero: the ligand interaction coordinates have collapsed;
3. both are nonzero but the quadratic form is near zero: protein and ligand
   coordinates are misaligned or rank-limited;
4. both coordinates and gradients are healthy but predictions are wrong: the
   representation or interaction basis is biologically mismatched.

The audit must report, per mutation pair, true target variance, predicted
variance, prediction-to-target scale ratio, slope, centered MSE, explained
variance, finite-rank denominator, and parent identity. Rank metrics must state
`N_rank_evaluable/N_total`; undefined Spearman for a constant predictor must not
be silently treated as observed zero correlation. Centered MSE, sign accuracy,
slope, and explained variance are the defined comparisons for the zero baseline.

The audit must also compute the interaction-parameter gradients separately:

```text
g_abs = grad_theta_s(L_abs)
g_ctr = grad_theta_s(L_contrast)
R_g   = ||g_abs|| / (||g_ctr|| + eps)
C_g   = cosine(g_abs, g_ctr)
```

These quantities must be checked at initialization and across the training
trajectory when checkpoints exist. A large `R_g` with negative `C_g` supports an
objective-conflict hypothesis. Healthy gradients with collapsed potential
variance shift suspicion to representation or parameterization. A final loss
scalar comparison alone is insufficient.

The input audit must compare KLIFS WT/variant differences with the previously
qualified per-position local ESM differences, but must not train a successor.
It must also report the effective rank of the admissible KLIFS mutation-difference
matrix. The local ESM route is eligible only because it has independent Q1
capability evidence; it is not admitted merely because it is a larger model.

The decision after the read-only audit is fixed:

```text
representation bottleneck -> ESM-only successor
objective conflict         -> centered-only successor
both supported              -> preregistered 2x2 (KLIFS/ESM x joint/centered)
neither supported            -> potential-rank/capacity successor
```

No successor may change representation and objective simultaneously unless the
2x2 design is explicitly frozen before training. No successor may enter CIIP-1B
until it has nonconstant prediction coverage, centered-MSE improvement,
positive scale/slope, sign accuracy, and consistency across parents. A stable
free-pairwise advantage can identify a limitation of the current observable
potential, but it cannot be interpreted as evidence that biology is
non-integrable.
