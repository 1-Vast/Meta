# DCST-R14 label-blind transport-support preregistration

Date: 2026-07-28  
Status: frozen before implementation and execution

## Question

R6 learned a privileged-specific PLINDER structural mechanism, while R12 and
R13 found no usable ChEMBL strict dual-cold effect. R14 tests the missing
premise: whether the high-quality Stage-1 source has support in the Stage-2
feature domain and whether the frozen teacher remains pair-responsive there.

This is a diagnostic and route-selection gate, not an affinity experiment.
The implementation must not request or load the ChEMBL `affinity` column.
Confirmation and sealed features are excluded. ChEMBL development features
may be audited but may not fit an adapter, anchor bank, scaler, or source
weight.

## Frozen inputs

- Stage 1: firewalled, exact-target PLINDER train rows admitted by R6;
- Stage 2 fit domain: ChEMBL train entities and features only;
- audit-only target domain: ChEMBL development entities and features;
- protein view: frozen 32-segment ESM-2 cache and sequence 4-mer containment;
- ligand view: radius-2 1024-bit Morgan cache;
- mechanism view: frozen R6 privileged and NoPriv teachers.

All sampling is deterministic at seed 1729. Domain classifiers use balanced
source/ChEMBL samples and five-fold out-of-fold ROC AUC. Nearest-source
coverage is reported for every eligible target and for a deterministic sample
of at most 20,000 unique ChEMBL ligands.

## Mechanism responsiveness

For each sampled pair, record the R6 `32 x 8` distribution and the
content-weighted 256-dimensional structural moment. For every target, compare
the observed ligand moment with its mean under a fixed, diverse bank of 16
ChEMBL-train ligand anchors:

```text
delta_m(t,d) = m(t,d) - mean_a m(t,a).
```

The anchor bank is selected without affinity labels by farthest-first
Tanimoto traversal on ChEMBL-train Morgan fingerprints. It is frozen once and
is used for source, train, and development diagnostics. Report raw and
centered moment domain AUC, normalized mechanism entropy, and within-target
centered-moment RMS. No head is trained on affinity.

## Frozen route selection

Define:

- target overlap: at least 20% of ChEMBL-train targets have maximum
  source sequence 4-mer containment at least `0.40`;
- ligand overlap: at least 20% of sampled ChEMBL-train ligands have maximum
  source Morgan Tanimoto at least `0.40`;
- pair responsiveness: ChEMBL-train median target-level centered-moment RMS
  is at least 20% of the corresponding PLINDER value;
- strong domain separation: target, ligand, or centered-moment domain AUC is
  at least `0.80`.

Select exactly one route:

1. `ADVANCE_ACMP`: both overlap gates and pair responsiveness pass, with no
   strong separation. Train-domain anchor-centered mechanism prompting is the
   next model.
2. `ADVANCE_WEIGHTED_ACMP`: both overlap gates and pair responsiveness pass,
   but any domain axis is strongly separated. Fit train-only density-ratio
   weights for Stage-1 examples, retain the R6 source mechanism gate, then use
   the same anchor-centered prompt.
3. `STOP_PLINDER_SOURCE_EXPAND_STAGE1`: either entity overlap gate or pair
   responsiveness fails. No further PLINDER-to-ChEMBL architecture variant is
   authorized until the high-quality source is expanded or replaced.

Development metrics are audit-only: disagreement between train and
development is reported as an inductive-risk warning and never changes the
route selected from source plus ChEMBL train.

## Required audit

The JSON artifact must include file hashes, exact row/entity counts, selected
anchors, sampling indices or hashes, all three domain results, the selected
route, `downstream_affinity_columns_loaded: false`,
`confirmation_features_loaded: false`, and `sealed_test_consumed: false`.

