# Experimental evidence ledger

Updated: 2026-08-11. This is the compact human-readable ledger. Machine results
remain in their current report directories; removed historical details are
recoverable from Git commit `c05d3f95fe59f1f0b1e1cc34163ba473f16ea008`.

| Stage | Evidence | Verdict | Consequence |
|---|---|---|---|
| Frozen law | Operator and contract tests | PASS | Keep `A(F,z)=K(B(z)F(z))` unchanged. |
| P1B geometry | Contact/distance bridge and partner controls | PASS | Geometry is a valid upstream input, not affinity evidence. |
| T-BASIS | 288D radial chemistry reconstruction and partner dependence | PASS structural | Candidate biological measurement only. |
| S7/SSL residue routes | Generic pocket, marginal and correspondence audits | TERMINATED | Exact transferable affinity mechanism was not identified; do not restore larger heads. |
| ChEMBL X1A/X1A-R | Direct-DD dependence and effective information | FAIL precondition | No population interaction claim or model training from that route. |
| BindingDB CQ-R0/R1 | 12,457 Ki cells, 320 panels, quotient interaction | PASS development | Open source optimization is executable; strict closure is not claim-ready. |
| BindingDB CQ-R2 | Shared 288D panel-balanced linear witness | FAIL | One universal radial affinity direction is rejected. |
| Target coefficient meta-learning | Main-v0, `d=2`, k=5, five seeds | PASS target and cluster sensitivity | Correct support beats d=0, zero, foreign and permuted controls. |
| Biological specificity | Full versus ligand-only and wrong protein | PARTIAL/FAIL admission | Ligand-only contrast survives cluster sensitivity; wrong-protein contrast does not. |
| V1 pair-prior repair | Full-rank pair prior and residual adapter on source/meta-val | FAIL development | Pair d0 and V1-A/V1-B harm absolute generalization. |
| V1 measured contrasts | Within-panel ligand and measured cross-family partner differences | PARTIAL | Permuted gap improves, but wrong/wrong recovery rejects correct-partner identification. |
| A0 gauge audit | Synthetic invariance plus sealed correct/wrong episode coordinates | FAIL exact global gauge | Partial local alignment remains possible; do not claim one orthogonal gauge or causal dominance. |
| A1 selectivity probe | 1,820 groups, 21 closed components, matched ridge/null probes | FAIL | Frozen calibrated T-BASIS selectivity signal not identified; A2 remains closed. |
| Biological `z` admission | Partner, affinity, transfer and support Gates | NOT RUN | Production assembly remains closed. |

## Current quantitative evidence

```text
BindingDB Ki cells                 12,457
quotient-positive panels             320
strict dependency components           31
largest component share              0.8586
train/development quotient rank       6608 / 220
shared-linear explained fraction      0.000709
correct-zero loss reduction           +0.000239 [-0.000981, +0.001496]
correct-foreign loss reduction        +0.000870 [-0.001045, +0.002847]
correct-deranged loss reduction       -0.000817 [-0.003419, +0.001498]
```

The failed shared witness does not test whether target coefficients occupy a
small transferable subspace. The next experiment must isolate that hypothesis
through target-held-out episodes and must not modify the biological basis and
optimizer architecture simultaneously.

## Meta-fewshot Phase 0 (2026-08-10)

```text
verdict                        FEWSHOT_EPISODE_DATA_NOT_IDENTIFIABLE
affinity label reads                                              0
endpoint                                                    Ki only
source targets / usable at k=5                          442 /  220
evaluation targets / usable at k=5                       68 /   16
evaluation usable at k=5, scaffold-disjoint                       8
declared minimum evaluation targets                              30
MDE_d at k=5 (targets) / declared ceiling             0.622 / 0.600
leakage: target, ligand, scaffold, document, homology-40           0
```

The corpus trains but cannot evaluate. Episodic source supply is ample and the
closure governance is clean; the constraint is held-out target depth. No model
was preregistered or trained, so target-coefficient heterogeneity is untested
rather than refuted.

A label-blind follow-up rejected the data-absent explanation:

```text
verdict            EVALUATION_PANEL_LIMITED_BY_ESTIMAND_MISMATCHED_CORPUS
same projection, few-shot rules only:
  single-chain Ki rows                         25,072   (corpus kept 12,457)
  distinct Ki targets                             910   (corpus kept    510)
  targets usable at k=5                           394   (corpus had     236)
  >=8 ligands and >=2 documents                   218
MDE_d: 16 -> 0.622, 30 -> 0.454, 50 -> 0.352, 100 -> 0.249
```

The cycle-positive quotient filter belongs to the crossed-rectangle estimand,
not the few-shot one. Independent component depth after closure is still
unmeasured.

## O1 corpus rebuild (2026-08-10)

`research/meta_fewshot/fs_corpus_rebuild.py` ran in `drug` on the governed
metadata projection with zero affinity-label reads. The preregistered exact-Ki,
single-chain, scaffold-valid corpus has 21,473 canonical cells and 880 targets.
After document, Murcko-scaffold and protein-40 closure it has 66 components,
22 with a `k=5`-eligible target. The frozen largest-component-to-source split
has 336 source and 33 evaluation `k=5` targets (`MDE_d=0.433`), passing the
unchanged `>=30` and `<=0.600` requirements. The giant source component still
contains 94.54% of cells; retain component-level sensitivity analysis.

This is a structural availability result. A reproducible hash-matched,
identity-only diagnostic rebuilt the same 66 roots and found
exact numeric Ki for 18,509/21,473 pairs and only 313/29 source/evaluation
targets at `k=5`; exact numeric Ki plus the old quotient-shaped 288D bank covers
224/12. It used no affinity value or aggregation and counted 1,389 replicated
admitted pairs. Therefore numeric FS-C0/FS-C1 and full feature materialization are open
preconditions. No real Meta-Section training is authorized. The separately
implemented synthetic Meta-Section passed all registered implementation checks;
that result validates the algebra and code path, not biological affinity.

The proposed Theory-Projected Q-PMA failed architecture review: its attention
prior can inject support-residual null-space components and a separate task
state per query. A row-space-only repair is a new support-feature-conditioned
family, not an implementation of the frozen linear family. It remains a future
amendment and was not moved into `model/` or `scripts/`.

The benchmark policy is literature-referenced rather than copied. AdaMBind's
protein-as-task and 40% novel-target precedent and CARA's assay-aware cleaning
and task-level evaluation are candidate conventions for the main benchmark;
MetaSieve retains the Ki estimand, biological representation, identifiable
section and law output. The existing dependency-closed split is retained as a
strict confirmation lane. The unresolved assay-to-protein-task aggregation is
registered as a corpus-contract blocker, not chosen from downstream counts.

## MetaSieve-main v0 (2026-08-11)

The preregistered method-level lane materialized 17,717 exact-Ki observations,
499 protein tasks and `285/37/33` k=5 train/validation/test tasks under CD-HIT
40% complete-cluster 8:1:1 assignment. A full 17,717-row frozen T-BASIS bank
was generated, and five CUDA training repeats selected `d=2`, ridge `1.0`.

Target-level one-sided 95% lower bounds for MSE reduction were positive against
d=0 (`4.179`), zero (`3.750`), foreign (`3.199`), permuted (`0.060`),
ligand-only (`1.046`) and wrong protein (`0.020`). Cluster sensitivity over six
eligible test clusters passed the first five but failed wrong protein (mean
`-0.081`, LCB `-0.227`). One cluster contains 21/33 test targets and is the
only cluster favoring correct protein strongly enough to reverse the target
average. Hence meta-learning and support specificity are identified on this
benchmark, but independently replicated partner biology is not.

The registered target verdict is retained as an audit fact; the governing
scientific verdict is `BIOLOGICAL_SPECIFICITY_NOT_IDENTIFIED_CLUSTER_SENSITIVITY`.
No production migration, Q-PMA, CSMO bridge or strict confirmation is
authorized. See `report/meta_fewshot/main_v0/MAIN_V0_REPORT.md`.

## V1 biological-axis development (2026-08-11)

The physically sealed run used no main-v0 test values. Source supervision
contains 1,820 measured same-panel/same-ligand cross-family partner groups, so
V1-B used observed pKi differences rather than random missing-edge negatives.

V1-B moved the correct-permuted cluster LCB above zero (`0.039`) but had MSE
`3.890`, worse than v0 (`1.800`) and ligand d0 (`3.084`). Pair d0 was also worse
than ligand d0 by cluster-macro `1.083`. Most importantly, a 2x2 biology
factorial found that wrong support plus wrong query restored the correct error:
v0 `1.765` versus correct `1.800`, and V1-B `3.866` versus `3.890`.

```text
V1_DEVELOPMENT_REPAIR_NOT_SELECTED
```

The active failure is not insufficient section capacity. It is that the frozen
statistic supports self-consistent adaptation under an arbitrary partner
coordinate without identifying the correct biological partner. No frontend
unfreezing, Q-PMA, CSMO bridge or production migration is authorized.

## Biological gauge and selectivity audit (2026-08-11)

The exact orthogonal ridge invariance passes synthetic controls, but no single
exact correct-to-wrong orthogonal gauge transfers across real clusters. Median
support/query/Gram/kernel residuals are `0.445/0.478/0.729/0.383`; the global
leave-one-cluster-out residual is `0.801`.

The A1 dependency Gate opened at 21 components (`MDE=0.543`), despite an 86.43%
giant component. Calibrated T-BASIS then failed both fixed-sequence contrasts:
rewired-null minus T-BASIS is `-0.291` (LCB `-0.522`), and ESM-additive minus
T-BASIS is `-0.083` (LCB `-0.452`); a 999-draw fixed-hyperparameter diagnostic
gives uncalibrated `p=1.0`.
The planted positive control passes. A2 was not run.

The reported `p=1.0` is an uncalibrated fixed-hyperparameter diagnostic, not a
confirmatory randomization p-value: the single rewiring leaves 39.6% of rows
unchanged, the label permutation does not preserve repeated-family incidence,
and hyperparameters were selected on observed labels. The fail-closed decision
is robust to this defect because 15/21 components, all group-weighted
contrasts, and leave-giant-out rewiring/zero sensitivities disfavor T-BASIS.

```text
NO_SINGLE_EXACT_ORTHOGONAL_GAUGE_IDENTIFIED
TBASIS_SELECTIVITY_SIGNAL_NOT_IDENTIFIED
A2_NOT_RUN_GATE_CLOSED
```

The bounded conclusion is lack of registered linear/ridge decodability from
the frozen 288D statistic, not absence of all selectivity or local interaction
information. A fresh assay-matched dense crossed-selectivity cohort is required
before any confirmatory A1 rerun or A2 localization.

## R2 identifiability-regime correction (2026-08-11)

| Check | Result | Evidence status |
|---|---:|---|
| E0 meta-val target, pair-intercept vs full MSE | 1.441 vs 1.578 | Descriptive; consumed split |
| E0 meta-test cluster, full share beyond pair-intercept | 2.1% | Predominantly calibration; not uniform across clusters |
| E1 meta-val/meta-test gauge ratio | 1.031 / 1.099 | Both trigger preregistered H0 falsifier `>1` |
| E1 calibration share | 0.694 / 0.779 | Does not reach registered `>=0.9` support threshold |
| E2 2-core additive / interaction fraction | 0.981 / 0.019 | Label-redacted observed-design description, not capacity proof |
| E2 fixed-ligand partner dispersion | 0.051 | Partner variation is weak after current aggregation |
| E3 local fresh confirmation panels | 0 | Historical development supply only |
| RFMS quotient-exposure guarantee | Counterexample found | Training blocked |

```text
R2_H0_REGIME_FALSIFIED_BY_ITS_OWN_PREREGISTRATION
META_SECTION_PREDOMINANTLY_CALIBRATION_DESCRIPTIVELY
LOCAL_CROSSED_TRAINING_SUPPLY_EXISTS_NO_FRESH_CONFIRMATION
RFMS_TRAINING_NOT_AUTHORIZED
```

The accepted correction is narrower than the Cowork headline. For a fixed
support residual, the unregularized full-rank ridge kernel is GL-invariant and
positive ridge is orthogonally invariant. The complete v0 predictor includes a
learned population-coordinate term, and real wrong-protein replacement was not
shown to be one shared transform. `k>=d` also does not make auxiliary biology
useless under noisy, regularized or misspecified episodes.

The only retained model candidate is an explicit pair support-intercept plus a
centered residual/coordinate ridge section. Its post-hoc consumed-split values
select it for a future preregistered source-side comparison; they do not count
as model evidence. RFMS, Q-PMA, CSMO integration and production migration stay
closed. Full audit: `report/meta_fewshot/R2_MULTI_AGENT_RESOLUTION.md`.

## Historical lessons retained

- More attention or a larger pair head did not establish affinity direction.
- Generic pocket localisation can survive wrong ligands and is not sufficient.
- Binary/local structural teachers can be dominated by marginals.
- Crossed rectangles and quotient rank do not by themselves create independent
  scientific units.
- A failed numerical or synthetic control is not a biological failure.
- Optimization authorization on a large open corpus is distinct from admission
  of a population claim.
- A filter inherited from a previous estimand can masquerade as a data limit;
  check that the corpus was built for the question being asked.
- Clean leakage governance does not imply an evaluable panel: a corpus can be
  well-governed, leakage-free and large enough to train on while still being too
  shallow on the held-out side to evaluate the estimand it was built for.

Full chronology and exact historical verdicts are in `history.md`.

K1 cross-dataset closure: centered residual correction beat the intercept arm but
lost to the existing uncentered ridge on BindingDB, Davis, and KIBA. Its
implementation branch was retired after the failed Gate; preregistration and
result artifacts remain. See `history.md` F-124 and
`report/meta_fewshot/METASIEVE_RESEARCH_SUMMARY.md`.
