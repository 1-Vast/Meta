# MetaSieve current research status

Updated: 2026-08-11.

## Objective

The sole final task is unseen-target few-shot DTA at `k=1/2/3/5`. Structural
localisation and open-data quotient learning are upstream evidence, not the
product objective.

## Current verdict

```text
OPEN_BINDINGDB_QUOTIENT_TRAINING_PIPELINE_EXECUTABLE
CQ_TBASIS_LINEAR_AFFINITY_WITNESS_NOT_OBSERVED
REAL_META_SECTION_META_EFFECT_IDENTIFIED
REAL_META_SECTION_SUPPORT_SPECIFICITY_IDENTIFIED
REAL_META_SECTION_LIGAND_ONLY_CONTROL_BEATEN
BIOLOGICAL_SPECIFICITY_NOT_IDENTIFIED_CLUSTER_SENSITIVITY
NO_SINGLE_EXACT_ORTHOGONAL_GAUGE_IDENTIFIED
TBASIS_SELECTIVITY_SIGNAL_NOT_IDENTIFIED
R2_H0_REGIME_FALSIFIED_BY_REGISTERED_E1
META_SECTION_PREDOMINANTLY_CALIBRATION_DESCRIPTIVELY
LOCAL_CROSSED_TRAINING_SUPPLY_EXISTS_NO_FRESH_CONFIRMATION
RFMS_TRAINING_NOT_AUTHORIZED
CENTERED_KERNEL_ALGEBRA_EVIDENCE_RETAINED
CENTERED_SECTION_CROSS_DATASET_FAIL
NO_BIOLOGICAL_Z_ADMISSION

Compact index: `report/meta_fewshot/METASIEVE_RESEARCH_SUMMARY.md`.
```

## MetaSieve-main v0 real training (2026-08-11)

The strict 29-target confirmation lane no longer gates the separate main
method-level protocol. The preregistered main corpus applies exact-Ki cleaning,
one protein sequence per task, and CD-HIT 40% complete-cluster 8:1:1 splitting.
It materialized 17,717 observations across 499 targets; k=5 train/validation/
test task counts are `285/37/33`. The complete frozen 288D feature bank was
generated before training.

Five real training repeats selected `d=2`, ridge `1.0`. The target-level
preregistered Gates all passed: full correct beat d=0, zero, foreign, permuted,
ligand-only and wrong-protein controls with positive one-sided 95% lower
bounds. Full-correct target-macro MSE was `1.916`, versus `8.711` for d=0 and
`3.426` for ligand-only.

The scientific admission nevertheless fails. The 33 test targets occupy only
six eligible CD-HIT clusters. Cluster-bootstrap sensitivity retains every
contrast except correct versus wrong protein, whose cluster-macro MSE reduction
is `-0.081` with one-sided 95% LCB `-0.227`. The target-level positive result is
driven by one cluster containing 21 targets; all five other cluster means favor
the wrong-protein arm. Correct-protein ranking metrics also do not beat wrong
protein. Therefore partner-specific biology is not independently identified.

```text
registered target-level verdict  REAL_BIOLOGICAL_META_SECTION_V0_PASS
scientific admission verdict     BIOLOGICAL_SPECIFICITY_NOT_IDENTIFIED_CLUSTER_SENSITIVITY
production migration             NOT AUTHORIZED
```

Absolute generalization remains weak (`R2=-1.244`, Pearson `0.097`, Spearman
0.086). The frozen frontend was not trained end-to-end, and CSMO law metrics
remain `NA_NOT_ADMITTED` because no 28D biological z map exists. The next
amendment must target partner specificity on source/meta-validation only and
confirm on fresh protein clusters. The consumed main-v0 test set cannot support
a confirmatory retest of an amended model.

## V1 source/meta-validation repair (2026-08-11)

A physical development seal excluded all 50 main-v0 test targets and 1,934
test cells. Meta-validation model inputs contain support pKi but no query pKi;
predictions were persisted before the independent truth file was read. New
wrong-protein controls use only length/composition-matched meta-train donors
from different CD-HIT40 groups.

A source census found 1,002 measured within-panel ligand groups and 1,820
same-panel/same-ligand partner groups crossing CD-HIT40 families. These support
scientifically valid measured-difference losses without treating missing
BindingDB edges as non-binders.

Five-seed development compared cluster-balanced v0, V1-A with a shared
full-rank pair prior and residual pair adapter, and V1-B with additional
measured within-panel/partner difference losses. V1-B stabilized the permuted
support contrast, but failed the frozen selection criteria:

```text
target-macro MSE
  v0 correct       1.800
  ligand d0        3.084
  pair d0          3.806
  V1-A correct     4.204
  V1-B correct     3.890
```

The decisive 2x2 protein factorial shows that changing only support or only
query protein is destructive, while changing both to the same wrong protein
restores performance. V0 wrong/wrong MSE is `1.765` versus correct `1.800`;
V1-B wrong/wrong is `3.866` versus correct `3.890`. Thus the section currently
uses a self-consistent coordinate system but does not require correct protein
identity.

```text
V1_DEVELOPMENT_REPAIR_NOT_SELECTED
END_TO_END_FRONTEND_NOT_AUTHORIZED
PRODUCTION_MIGRATION_NOT_AUTHORIZED
```

The full report is `report/meta_fewshot/v1_development/V1_DEVELOPMENT_REPORT.md`.

## Biological gauge and selectivity audit (2026-08-11)

A preregistered source-only audit tested the proposed gauge explanation before
changing the architecture. The exact orthogonal invariance is correct in
synthetic controls, but the real correct/wrong coordinates are not one shared
orthogonal transform: median support Procrustes, held-out query-transfer, Gram
and `H_CC-H_WW` relative errors are `0.445/0.478/0.729/0.383`; the global
leave-one-cluster-out residual is `0.801`. Partial local kernel alignment and
population/section cancellation remain possible.

For measured selectivity, 1,820 source groups close into 21 dependency
components (`MDE=0.543`); one giant component contains 86.43% of groups. Under
component-held-out, capacity-matched ridge probes, T-BASIS MSE is `0.926`
versus `0.585` for zero and `0.635` for the rewired coupling null. The
registered coupling-null and ESM-additive loss reductions are `-0.291` and
`-0.083`, with one-sided LCBs `-0.522/-0.452`; a 999-draw fixed-hyperparameter
diagnostic gives uncalibrated `p=1.0`. The planted positive control passes.

The randomization result is diagnostic rather than confirmatory. A single
rewiring leaves 39.6% of rows fixed; groupwise label permutations do not
preserve repeated-family incidence; and ridge hyperparameters were selected on
the observed labels. Directional failure is nevertheless robust: 15/21
components favor controls, all group-weighted contrasts are negative, and the
rewiring/zero contrasts remain negative after removing the giant component.

```text
NO_SINGLE_EXACT_ORTHOGONAL_GAUGE_IDENTIFIED
TBASIS_SELECTIVITY_SIGNAL_NOT_IDENTIFIED
A2_NOT_RUN_GATE_CLOSED
```

This does not prove that selectivity information is absent or localize a loss
inside the frozen frontend. It closes the current 288D linear-selectivity claim
fail-closed, without claiming a calibrated coupling-specific p-value, and
forbids result-driven A2 probing. The next scientific input must be a fresh,
assay-matched dense crossed-selectivity cohort with independent dependency
components, not another readout or architecture amendment. Full details are in
`report/meta_fewshot/BIOLOGICAL_GAUGE_AUDIT_REPORT.md`.

## R2 identifiability-regime resolution (2026-08-11)

Cowork's narrower calibration diagnosis is supported descriptively, but its
main H0-regime claim is rejected by the preregistration's own E1 falsifier.
`gauge_ratio` is `1.0306` on meta-val and `1.0987` on meta-test, rather than the
registered near-gauge threshold `<=0.5`; calibration share is `0.694/0.779`,
rather than the registered near-uniform threshold `>=0.9`. The complete v0
predictor also contains a learned population-coordinate term, so a fixed-
residual unregularized GL identity is not a proof about wrong-protein outcomes.

A matched pair-intercept yields MSE `1.4408` versus full `1.5780` on meta-val
target macro and `1.8965` versus `1.9162` on meta-test target macro. At
meta-test cluster macro, full improves over the pair-intercept by only `0.0453`,
about 2.1% of the total pair-support gain. The direction is not consistent in
all clusters, so the bounded conclusion is predominantly calibration on
average, not universal absence of ligand-specific work.

E2 now reads a physically label-redacted structural index. On the 11,278-row
bipartite 2-core, 98.07% of frozen-feature variance is additively explained and
fixed-ligand partner dispersion is 5.13%. This is an observed-design linear
description, not a capacity impossibility; the frontend does use ESM residue
states to condition distance logits before its coarse radial aggregation.

E3 finds ample historical crossed-panel development interaction df in Metz and
PDSP, but no fresh confirmation: each closes to one dependency component, and
the BindingDB panel supply has only eight components with 91.2% in the largest.
Substantial protein/ligand overlap with main-v0 further rules out calling these
packages fresh.

RFMS is not authorized. Its wrong/wrong guarantee fails unless
`Xi(c_correct-c_wrong)` is shown nonzero, where
`Xi=b-a(A^T A+lambda I)^-1 A^T C`; nonconstant `c0` alone is insufficient.
The next minimal candidate is an explicit pair-intercept plus centered
residual/coordinate section, preregistered on new source-side development
components. See `report/meta_fewshot/R2_MULTI_AGENT_RESOLUTION.md`.

BindingDB Articles 202608 yields 12,457 governed Ki cells, 320 panels, train
quotient rank 6,608 and development quotient rank 220. Strict closure has only
31 conflict components and its largest component holds 85.86% of cells, so the
corpus supports source optimization but not a population-wide claim.

The first real training run fitted one shared linear response on the frozen
288D T-BASIS. It explained `0.000709` of development quotient variance, and no
correct-versus-zero/foreign/deranged confidence interval excluded zero. The
population-shared radial direction failed. Target-specific coefficient
heterogeneity remains untested.

## Phase 0 episode feasibility — failed closed

A label-blind episode census ran on the governed Ki corpus with zero affinity
label reads. The split is clean: target, ligand, scaffold, document and
protein-homology-40 leakage are all exactly zero. The source side is ample —
442 targets, 220 usable at `k=5`.

The evaluation side is not. Only 68 development targets exist, of which
`24/19/18/16` can carry `k=1/2/3/5`, and only `24/18/9/8` with a
scaffold-disjoint support set. At `k=5` that is 16 held-out targets against a
declared requirement of 30, giving `MDE_d = 0.622` against a declared ceiling of
`0.600`.

```text
FEWSHOT_EPISODE_DATA_NOT_IDENTIFIABLE
```

No model was preregistered and none was trained. The scientific hypothesis —
target-specific coefficient heterogeneity — was **not** tested and is neither
supported nor refuted by this run.

A label-blind follow-up localized the cause and **rejected the data-absent
explanation**. Recounting the same governed projection under few-shot rules
only, dropping the cycle-positive quotient requirement, yields 25,072
single-chain Ki rows and 910 targets, of which 584/499/459/**394** support
`k=1/2/3/5` and 218 have `>=8` ligands across `>=2` documents. The
quotient-shaped corpus retained 12,457 cells and 236 `k=5` targets.

```text
EVALUATION_PANEL_LIMITED_BY_ESTIMAND_MISMATCHED_CORPUS
```

The CQ corpus required cycle-positive panels because the crossed-rectangle
quotient estimand needs closed rectangles; the few-shot estimand needs only
per-target ligand depth. Independent *component* depth after protein-40,
scaffold and document closure is still unmeasured, so the giant-component
pathology may persist. The next stage is a preregistered few-shot-shaped corpus
rebuild under unchanged closure and unchanged Phase 0 thresholds.

## O1 corpus rebuild (2026-08-10)

The preregistered label-blind rebuild has passed FS-C0 and FS-C1. From the same
metadata projection, exact single-chain Ki admission and scaffold-valid
canonicalization yielded 21,473 cells across 880 targets. Strict document,
Murcko-scaffold and protein-40 closure yielded 66 components; 22 contain a
`k=5` eligible target. The frozen rule assigns the largest component to source
and every remaining complete component to evaluation: source/evaluation `k=5`
targets are `336/33`, with target-unit `MDE_d=0.433`.

```text
FEWSHOT_CORPUS_STRUCTURALLY_AND_POWER_FEASIBLE_FROZEN_SPLIT
```

No affinity label was read and no model was trained. The large source component
contains 736/880 targets and 94.54% of cells, so component-level sensitivity
must remain prominent in later inference.

### Numeric-corpus amendment

A reproducible identity-only, hash-matched availability audit showed that this
PASS is structural rather than executable. It rebuilt the same 66 frozen roots
and used no affinity value or aggregation. Only 18,509/21,473 canonical pairs have
an exact numeric Ki row. Under the frozen dependency assignment this reduces
`k=5` source/evaluation targets from `336/33` to `313/29`; requiring the old
288D feature bank reduces them further to `224/12`. The value 29 is diagnostic,
not a new formal Gate, because numeric admission and aggregation of the 1,389
replicated admitted pairs have not yet been preregistered. The machine-readable
result is `report/meta_fewshot/FS_NUMERIC_AVAILABILITY_AUDIT.json`.

The minimal Meta-Section passed its synthetic control and the separate
MetaSieve-main method-level real-training run is complete. The O1 strict lane
remains closed and unchanged; it did not authorize or supply labels to main-v0.

The proposed Theory-Projected Q-PMA failed the architecture Gate and was not
moved into production. AdaMBind's three public CSVs were downloaded from pinned
commit `01a169a6...` for audit, but exact reproduction remains blocked by the
restricted versioned archive and missing CD-HIT split manifests. The literature
and MetaSieve protocols remain separate.

The main benchmark policy is now reference-based rather than copy-based.
AdaMBind supplies a comparable protein-task/novel-target precedent and CARA
supplies assay-aware numeric-cleaning and per-assay evaluation precedents.
MetaSieve retains its own Ki estimand, biological interaction representation,
support-identifiable section and gated law output. The former dependency-closed
split remains a strict confirmation stress test rather than the sole source of
training tasks. CARA's assay-level median and AdaMBind's protein-level task do
not compose automatically; their BindingDB assay-to-target reconciliation must
be frozen before a new main corpus is materialized.

## Active route

Freeze the negative R2 result and do not train RFMS. The user has authorized a
staged replacement of failed modules. The retained core is target-as-task
episodic learning plus support-only positive ridge. The first executable change
is an explicit support intercept plus centered `d<=5` section, with the original
linear arm retained. An exact-residue pair field may replace T-BASIS only after
structural and measured-crossed admission; a rich PSD kernel is locked behind
that linear biological Gate. The research-only centered kernel core and tests
are implemented, but no new biological model has been trained.

The preregistered K1 cross-dataset Gate subsequently rejected migration of the
centered predictor. Across 94 BindingDB, 75 Davis and 42 KIBA evaluation targets,
the centered section beat the pure intercept in every dataset but lost to the
retained uncentered ridge in every dataset. Therefore the uncentered positive
ridge remains the predictor; the support intercept becomes a mandatory
calibration-only control. No code moved to `model/` or `scripts/`.

Historical crossed panels may support additive versus low-rank PCM/IMC
development under protein-family and scaffold hold-out; they are not
confirmation supplies. The complete retain/replace/Gate contract is
`report/meta_fewshot/METASIEVE_V2_INTEGRATED_ARCHITECTURE.md`.

Evaluation must be target-family, scaffold and document held out. Correct
support must beat zero, foreign and permuted support while reporting rank,
conditioning, query coverage and abstention. Only a partner-specific,
affinity-incremental and independently replicated statistic can later enter
the unchanged law operator.

## Boundaries

- The 288D basis is structurally validated but not affinity-admitted.
- Research-only main-v0 target subspaces were trained and evaluated, but none
  is admitted to `model/`; further real training requires a new preregistration.
- Public AdaMBind Davis/KIBA CSVs are consumed for the K1 engineering Gate;
  sealed Davis/KIBA/recipient and external confirmation labels remain closed.
- `model/`, production `z`, CSMO, Band and `A(F,z)=K(B(z)F(z))` are unchanged.
- Terminated S7/SSL/correspondence/X1A details were removed from the active tree;
  their verdicts remain in `history.md`, the evidence ledger and Git commit
  `c05d3f95fe59f1f0b1e1cc34163ba473f16ea008`.

## Current evidence

1. `report/crossed_interaction/OPEN_DATA_TRAINING_AND_FEWSHOT_ROUTE.md`
2. `report/crossed_interaction/cq_r2_tbasis_linear/weights.report.json`
3. `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`
4. `task.md`
5. `history.md`
