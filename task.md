# Current task

Canonical project status, retained modules, terminal decisions and archive
boundaries are summarized in [`PROJECT_SUMMARY.md`](PROJECT_SUMMARY.md). This
file defines requirements; `history.md` remains the chronological record.

## Non-negotiable objective

Build and validate a **trainable few-shot drug-target affinity predictor for
unseen targets**. Targets are meta-learning tasks. The model must use large
open datasets to learn transferable biological knowledge, then adapt from
`k=1/2/3/5` support affinities without reading query labels.

Success requires all three layers:

1. **Biology:** a partner-specific protein-ligand measurement that adds
   affinity information beyond ligand-only and wrong-protein controls.
2. **Meta-learning:** source target episodes learn a shared low-dimensional
   mechanism basis; support labels estimate only the new target's identifiable
   section.
3. **Mathematics:** the admitted bounded statistic enters the unchanged
   `A(F,z)=K(B(z)F(z))` operator with explicit rank, conditioning, query
   coverage and abstention.

## Minimal model contract

```text
phi(P,L) in R^288                 frozen audited biology candidate
U in R^(288 x d), d <= 5          source/meta-learned task subspace
m(P,L) = U^T phi(P,L)
f0 = f_L(L, endpoint) + w0^T m(P,L)
a_t = positive-ridge support section in row(M_support)
y_hat(P_t,L_q) = f0 + a_t^T m(P_t,L_q)
```

`U` and `w0` are learned through target-wise support/query episodes. This is
the trainable meta-learning core. The closed-form section replaces a free
MAML inner network; it is smaller, deterministic and aligned with the support
identifiability requirement.

## Core reference: AdaMBind

Wan et al., [A meta learning and task adaptive approach for drug target
affinity prediction](https://www.nature.com/articles/s41467-026-70554-5)
(Nature Communications, 2026; [official
repository](https://github.com/Moohyun-w/AdaMBind)) is a core methodological
reference for this programme.

Adopted principles:

- define one target and its measured ligands as one task;
- keep support adaptation and query evaluation explicit during source
  meta-training;
- use protein-sequence-disjoint evaluation as a minimum, strengthened here by
  ligand-scaffold and source-document dependency closure;
- treat source-task selection as a learnable training intervention rather than
  test-time freedom;
- test training-label perturbation prospectively as a robustness intervention.

Translation boundaries:

- retain the uncentered positive analytic ridge and Identifiable Meta-Section;
  do not import full-model MAML or test-time neural fine-tuning;
- `ATS-Section` may change only source-task sampling, must learn from dependency-
  isolated meta-validation targets, sample without replacement, and must not be
  loaded at test time;
- label noise is a separate arm from task scheduling and biological-pair work;
  it is applied only during source training, never to sealed evaluation labels,
  and its scale must be endpoint/assay aware rather than silently mixing Davis,
  KIBA, Ki, Kd or Kdapp semantics;
- the public AdaMBind 8:1:1 and CD-HIT-40 protocols are engineering references,
  not sufficient evidence of MetaSieve's stricter leakage control;
- scheduler or noise gains cannot admit a biological representation unless
  ligand-only, correct-versus-wrong-protein and support-specificity Gates also
  pass.

The scheduler experiment uses a matched-capacity null that permutes task
statistics across candidate tasks while preserving scheduler size, candidate
pool, selected-task count and compute. The noise experiment must separate
support-noise from the zero-noise arm. Independent zero-mean query-label noise
is excluded because under squared loss it is an expected-gradient null and
would contaminate the clean evaluation path, whereas support noise changes the
adapted section. Per the current development priority, V1 executes the paired
uniform/ATS x clean/support-noise experiment first while the exact-pair R0
work remains paused. Admission still requires absolute Cold Target quality,
ligand-only, support-specificity and wrong-protein Gates; an engineering PASS
does not override those biological Gates.

### Meta-learning Transfer & Robustness Audit (2026-08-11)

The static audit of the AdaMBind article, complete Supplementary Information,
official commit `01a169a6`, current `TaskScheduler`, support-noise path and
Tier 1/2/3 literature is complete. It did not train a model, open confirmation
labels or change frozen theory.

Current ATS and fixed support noise remain failed mechanisms: ATS did not beat
its equal-capacity permuted-statistics null at every k, support noise did not
beat clean training, and neither repaired support specificity or partner
identity. The official AdaMBind code is reference-only because its published
flags and task-role separation do not faithfully support direct scheduler and
noise ablations. Full-model MAML, TaskNorm and task-conditioned metric
replacement are not authorized.

Three research-only candidates remain, in order: (A) an offline cross-fitted
task reliability/transferability audit with biological corruptions and
familiarity covariates; (B) assay/replicate-aware uncertainty propagated through
the existing positive ridge; and (C) a support-resampling stability objective
for the same low-dimensional section, conditional on A/B finding a residual
instability. Gate 0 is complete. Gates 1-3 and fresh confirmation are not
passed; no candidate may enter production or override F-132. The unified audit
is archived, with its binding decisions consolidated in `PROJECT_SUMMARY.md`.

### Task reliability/transferability Gate 1 (2026-08-11)

Candidate A has now been executed as a formal development-only CUDA audit.
Gate 0 passed (`285` k5 source targets/`207` components; `37` meta-validation
targets/`9` components; largest validation component share `0.243`). The
cross-fitted scorer used only query difficulty and support/query gradient
alignment. Scaffold overlap, task size, label spread and provenance density
were nuisance-audit columns and did not enter fitting; empty Murcko strings
were missing, not a shared scaffold. Empirical replicate disagreement and
continuous protein familiarity remain `NA`.

Gate 1 failed every k. Clean score-to-first-order-utility Spearman was
`-0.060/-0.163/0.010/-0.089` at `k=1/2/3/5`; no k met the frozen absolute
`rho>=0.20, LCB>0` criterion, matched-null superiority failed all k, and the
protein-shuffle necessity LCB was negative all k. Therefore Candidate A is
`REJECT_TASK_SCHEDULER_GATE1_FAIL_CLOSED`: do not modify the main training
scheduler and do not run the short end-to-end Gate. The code remains a
falsification harness. Candidate B remains a separate future axis only if real
assay/replicate uncertainty is available; it may not use fixed synthetic noise
as a proxy. Full evidence is archived; the binding verdict is consolidated in
`PROJECT_SUMMARY.md` and `history.md` F-133.

### Cold Target V1 development result (2026-08-11)

The full three-seed run was executed in the `drug` environment on CUDA. It
completed 24 hash-bound checkpoints and scored 1,251,600 prediction rows on a
shared 33-target cohort with nested `k=1/2/3/5` support and a fixed query set.
Meta-validation selected `uniform_clean`; ATS-clean beat its matched null on
aggregate meta-validation but not under the all-k cluster-bootstrap Gate, and
support-noise was not identified against clean training. Target-macro RMSE was
`1.582/1.429/1.363/1.313` and CI was `0.529/0.534/0.538/0.536` for
`k=1/2/3/5`. Absolute-quality, support-specificity and wrong-protein Gates all
failed. The terminal verdict is
`COLD_TARGET_FEWSHOT_V1_NOT_YET_GOOD`; no production or biological admission
is authorized.

### GPU vectorization result

Candidate episodes, dual Gram solves, per-task ATS gradients, validation and
prediction are now batched on CUDA; source/noise RNG order and the uncentered
positive-ridge mathematics remain fixed. Under the exact formal 1000-step
configuration, a uniform model fell from 57.3 to 4.97 seconds and an ATS model
from 106.1 to 20.86 seconds. The full three-seed/24-model reproduction fell
from 2,066 to 388 seconds (`5.32x`) with peak reserved memory below 150 MiB.
Frozen-checkpoint inference preserved all 20,860 episode identities and differed
by at most `5.25e-6` in standardized prediction.

Float32 training trajectories are not bitwise invariant: the full reproduction
selected `ats_clean` instead of `uniform_clean`. It nevertheless reproduced the
same terminal scientific failure: absolute quality, support specificity,
partner identity and ATS-vs-null all-k Gates remain failed. Vectorized kernels
therefore require fresh end-to-end checkpoints and may not be mixed with the
legacy scalar run. On larger machines, independent seed/arm jobs should be
distributed across devices or concurrent workers; model width must not be
increased solely to inflate utilization.

### Targeted V1 repair result

The post-V1 repair ladder retained only changes that improved the paired Cold
Target development task. A support-only section, a 64-wide nonlinear ligand
population, and meta-validation-selected `d=4, ridge=1` reduced RMSE by
`7.6-8.9%` versus vectorized V1 and increased CI by `0.025-0.037` across
`k=1/2/3/5`. Final RMSE was `1.495/1.357/1.290/1.230`; CI was
`0.549/0.553/0.558/0.567`. The nonlinear pair map was rejected after a severe
regression; wider population networks and population pretraining were rejected
by meta-validation. The final verdict remains
`V1_TARGETED_REPAIRS_MATERIAL_IMPROVEMENT_NOT_YET_GOOD` because absolute,
support-specificity and all-k partner Gates did not pass.

R0-B completed exact geometry and frozen-prior pre-fit scoring without affinity
labels. Its additive headroom Gate passed, but `MDE80=0.00746` exceeded the
registered `delta*=0.00616` with only 53 heldout components. It stopped before
training as `R0B_NOT_RUN_FAIL_CLOSED`. The next biological action is to add
fresh dependency-closed protein components; lowering the effect or power
threshold is forbidden.

## Open-data roles

- BindingDB curated Ki/Kd and Klaeger Kdapp: endpoint-specific quantitative
  constraints; never merge endpoint scales.
- Kinobeads, PKIS and PKIS2: within-panel ordinal/profile pretraining, not
  absolute Ki/Kd calibration.
- PDSP: non-kinase development stratum after panel/provenance census.
- The public AdaMBind Davis/KIBA CSV snapshots were consumed by K1 as
  cross-dataset engineering Gates and can no longer serve as independent
  confirmation for this architecture. Sealed recipient and future time-split
  data remain closed confirmation only.

Training data may be large and dependent. Population claims still require
document, assay, protein-homology, ligand-scaffold and publication-time
controls. Optimization authorization and scientific-admission authorization
are separate Gates.

## Phase 0 outcome (2026-08-10)

The episodic stage was attempted and stopped at its first precondition.

```text
FEWSHOT_EPISODE_DATA_NOT_IDENTIFIABLE
```

Label-blind census, zero affinity label reads, Ki only. Leakage is zero across
target, ligand, scaffold, document and protein-homology-40. Source supply is
ample (442 targets; 220 usable at `k=5`). The evaluation split is the binding
constraint: 68 targets, `24/19/18/16` usable at `k=1/2/3/5` and `24/18/9/8`
scaffold-disjoint, against a declared minimum of 30 and a declared power
ceiling `MDE_d <= 0.600` versus the observed `0.622`.

No model was preregistered, none was trained, and no threshold was moved. The
hypothesis of target-specific coefficient heterogeneity is **untested**, not
refuted.

A label-blind follow-up localized the cause and **rejected the data-absent
explanation**. Recounting the same governed projection under few-shot rules
only — single chain, Ki, `>= k+3` ligands — and dropping the cycle-positive
quotient requirement yields 25,072 rows and 910 targets, of which
584/499/459/**394** support `k=1/2/3/5` and 218 have `>=8` ligands across `>=2`
documents. The quotient-shaped corpus kept 12,457 cells and 236 `k=5` targets.

```text
EVALUATION_PANEL_LIMITED_BY_ESTIMAND_MISMATCHED_CORPUS
```

The cycle-positive filter belongs to the crossed-rectangle estimand, which needs
closed rectangles; few-shot needs only per-target ligand depth. `MDE_d` falls
from `0.622` at 16 targets to `0.249` at 100.

Not yet established: independent **component** depth after protein-40, scaffold
and document closure. Target depth is demonstrated; component depth is not, and
the 85.86% giant-component pathology may persist.

## Immediate next stage (revised)

R0-C has now closed the structural prerequisite with 219 fresh, independent
components and adequate registered power. The CUDA Full residual greatly beat
the frozen P1B prior, but the capacity-matched additive-exact N2 arm was better
than Full (`RPS 0.039302` versus `0.040522`). The exact-pair incremental and NLL
Gates failed with terminal verdict:

```text
MARGINAL_OR_SLOT_RECALIBRATION_ONLY
```

Do not proceed to the measured-affinity R1 field or connect its typed summary
to V1. The next biological proposal must explain what new observable breaks
the additive atom/residue shortcut before another structural training run is
authorized. Scheduler/noise work remains a separate robustness lane and cannot
override this partner-identifiability failure.

MetaSieve-main v0 has now run under its separate, preregistered method-level
protocol: exact Ki, one protein per task, CD-HIT 40% complete-cluster 8:1:1,
k=5 and five seeds. This did not open or replace the strict confirmation lane.

The registered target-level M1/M2/M3 Gates passed. A dependence sensitivity
over the six eligible test CD-HIT clusters retained the meta-effect, all three
support controls, and full versus ligand-only, but failed correct versus wrong
protein (`cluster-macro MSE reduction -0.081`, one-sided 95% LCB `-0.227`).
One cluster contains 21/33 test targets and drives the target-level biological
contrast.

```text
BIOLOGICAL_SPECIFICITY_NOT_IDENTIFIED_CLUSTER_SENSITIVITY
```

Production migration is not authorized. The next experiment must be
preregistered to address partner specificity in the trainable biological
coordinate using source/meta-validation data only, and must use fresh held-out
protein clusters for confirmation. Q-PMA and CSMO integration remain closed.

## V1 development outcome (2026-08-11)

The source/meta-validation repair experiment was run without loading any
main-v0 test value. It compared the cluster-balanced v0, a shared full-rank
pair prior (V1-A), and the same architecture plus measured within-panel ligand
and cross-family partner pKi differences (V1-B). No unmeasured pair was called
a non-binder.

V1-B enlarged the cluster-macro correct-versus-permuted reduction from `0.213`
to `0.390` and moved its one-sided LCB from `-0.034` to `+0.039`. It was not
selected: correct MSE `3.890` is worse than v0 `1.800` and ligand d0 `3.084`;
the full-rank pair d0 prior is itself worse than ligand d0. Correct-versus-wrong
query is positive, but using the same wrong protein for both support and query
recovers correct-arm performance (`v0 1.765` versus `1.800`; V1-B `3.866`
versus `3.890`).

```text
V1_DEVELOPMENT_REPAIR_NOT_SELECTED
```

This localizes the current mechanism as support/query-coordinate consistency,
not correct biological partner identity. End-to-end frontend work, Q-PMA,
CSMO and production migration remain unauthorized. The next admissible step is
a source-only information audit of the frozen interaction statistic, not more
readout capacity or reuse of the consumed test set.

## Biological gauge and selectivity audit outcome (2026-08-11)

The sequential source-only audit is complete. A0 did not identify one exact,
cross-cluster orthogonal gauge: real support Procrustes, held-out query transfer,
Gram and kernel residuals are all far from zero. Local support/query kernel
alignment remains a possible contributor to wrong/wrong recovery, so no causal
"primary explanation" is claimed.

A1 first closed 1,820 measured selectivity groups into 21 dependency
components (`MDE=0.543`) and then ran the registered component-held-out,
capacity-matched ridge probes. Calibrated T-BASIS was worse than zero, ESM
additive and a rewired coupling null. Both primary loss reductions were
negative and a 999-draw fixed-hyperparameter diagnostic gave uncalibrated
`p=1.0`; the planted positive control passed.

Independent review found that the single coupling rewiring leaves 39.6% of
rows fixed and the label permutation does not preserve repeated-family
incidence; ridge hyperparameters were also fixed after selection on observed
labels. Therefore `p=1.0` is an uncalibrated diagnostic, not a publication-ready
coupling test. The A1 result remains fail-closed because T-BASIS is directionally
worse than every control under component-macro, group-weighted and
leave-giant-out sensitivities.

```text
NO_SINGLE_EXACT_ORTHOGONAL_GAUGE_IDENTIFIED
TBASIS_SELECTIVITY_SIGNAL_NOT_IDENTIFIED
A2_NOT_RUN_GATE_CLOSED
```

The next admissible stage is a new, prospectively governed, assay-matched dense
crossed-selectivity cohort with independent document/scaffold/protein-family
components. Do not rescue the current result by changing probe capacity,
running A2, adding a partner anchor, unfreezing the frontend, manufacturing
wrong-partner labels or reopening main-v0 test. Production migration, Q-PMA and
CSMO remain closed.

## R2 Cowork diagnosis resolution (2026-08-11)

The Cowork report correctly identified a predominantly calibration-driven
section, but its stronger identifiability-regime headline is rejected. The
registered E1 `gauge_ratio` is `1.0306` on meta-val and `1.0987` on meta-test;
both exceed the preregistered `>1` H0 falsifier. Positive ridge and the learned
population-coordinate term also prevent extending the unregularized fixed-
residual GL invariance to the complete wrong-protein predictor.

E0 shows that a matched pair-intercept is competitive with or better than full
on average; on meta-test cluster macro the full-specific gain over that
intercept is only about 2.1% of the pair-support gain. This is descriptive on
consumed splits and is not uniform across clusters. E2, now using a physically
label-redacted index and a bipartite 2-core, finds 98.07% additive explained
feature variance and 5.13% fixed-ligand partner dispersion. This describes the
observed sparse design; it does not prove no nonlinear capacity. E3 identifies
historical crossed-panel development supply but no fresh confirmation supply.

```text
R2_H0_REGIME_FALSIFIED_BY_ITS_OWN_PREREGISTRATION
META_SECTION_PREDOMINANTLY_CALIBRATION_DESCRIPTIVELY
RFMS_TRAINING_NOT_AUTHORIZED
```

RFMS remains blocked: nonconstant partner coefficients do not prevent
wrong/wrong recovery unless the reserved quotient exposure is demonstrably
nonzero. The consolidated decision record is in `history.md`.

## Completed main-v0 experiment

The experiment changed only the failed sharing assumption:

```text
shared w FAIL
  -> source-learn U, w0 across target episodes
  -> evaluate unseen-target k=1/2/3/5 sections
```

It compared population d=0, correct, zero, foreign and permuted support,
capacity-matched ligand-only and wrong-protein arms. Full-correct target-macro
MSE was `1.916`, versus `8.711` for d=0 and `3.426` for ligand-only. Absolute
generalization remains weak (`R2=-1.244`, Pearson `0.097`), so this is a
mechanism screen rather than a performance claim. Full evidence is consolidated
in `history.md`.

## Complexity boundary

The user has now authorized replacement of failed modules, but not simultaneous
uncontrolled expansion. Preserve the target-as-task closed-form section and
change one scientific axis per registered arm:

1. first separate the support intercept from a centered `d<=5` section;
2. admit an exact-residue pair field only after structural and measured-crossed
   Gates;
3. expand the admitted linear feature kernel to a rich PSD kernel only after the
   retained linear arm passes biology Gates.

Do not add a new PLM, Q-PMA, MAML/support Transformer, knowledge graph or free
pair prior. A typed local field is a conditional replacement for failed
T-BASIS, not an automatically admitted module. No raw pair tensor enters `z`.

The first migration experiment has now rejected item 1 as a predictor change.
The centered section beat a pure intercept across public BindingDB, Davis and
KIBA, but was worse than the original uncentered ridge in all three datasets.
Keep the intercept as an evaluation control and retain the uncentered solver.
The next replaceable module is the failed biological pair representation, still
subject to R0/R1 admission before any production migration.
