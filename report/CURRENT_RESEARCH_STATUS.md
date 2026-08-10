# MetaSieve current research status

Updated: 2026-08-10.

## Objective

The sole final task is unseen-target few-shot DTA at `k=1/2/3/5`. Structural
localisation and open-data quotient learning are upstream evidence, not the
product objective.

## Current verdict

```text
OPEN_BINDINGDB_QUOTIENT_TRAINING_PIPELINE_EXECUTABLE
CQ_TBASIS_LINEAR_AFFINITY_WITNESS_NOT_OBSERVED
FEWSHOT_EPISODE_DATA_NOT_IDENTIFIABLE
TARGET_COEFFICIENT_META_LEARNING_NOT_YET_TESTED
K_SHOT_SECTION_NOT_IDENTIFIED
NO_BIOLOGICAL_Z_ADMISSION
```

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

## Active route

Use target-wise episodes over governed open datasets to learn a single
`d<=5` mechanism subspace on the frozen biological basis. Adapt a new target
only by a positive-ridge section in the support-observable row space. This is
the one trainable meta-learning hypothesis; no parallel adapter or additional
representation branch is authorized.

Evaluation must be target-family, scaffold and document held out. Correct
support must beat zero, foreign and permuted support while reporting rank,
conditioning, query coverage and abstention. Only a partner-specific,
affinity-incremental and independently replicated statistic can later enter
the unchanged law operator.

## Boundaries

- The 288D basis is structurally validated but not affinity-admitted.
- The target subspace and few-shot section are planned, not implemented.
- Davis, KIBA, recipient and external confirmation labels remain closed.
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
