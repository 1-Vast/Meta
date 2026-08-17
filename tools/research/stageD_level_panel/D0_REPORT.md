# Stage D0 report — re-audit of the Stage C boundary and the level/shape anatomy

Date: 2026-08-17. No training in this stage. Baseline: the leak-free Stage B T
checkpoint (report/meta_fewshot/stageB_complementary_20260817/T/checkpoint.pt).
Artifacts: D0_AUDIT_DECOMPOSITION.json, D0_LEVEL_IDENTIFIABILITY.json,
D0_LEVEL_ANATOMY.json, D0_OCCUPANCY_STRATA.json.

## Q1. Is the level/shape decomposition defined strictly per episode and per target?

**It is per episode (per draw), not per canonical target.** FEASIBILITY.json
decomposes each drawn 16-query panel separately and aggregates with equal
weight per target then per component. Re-running with the identical nested
bank reproduces the structure (level share ~62-65%; small numeric differences
come from the raw-episode versus target-weighted aggregation). Per-target
aggregation across the two draws moves level^2 by only ~0.01-0.04 pK^2 and the
panel-sampling variance of the drawn truth means is 0.013-0.034 pK^2, so the
drawn-panel component of 'level' is small. The level term is a between-target
quantity, but it is the level of the (protein, assay, tested-ligand-panel)
triple, not of the protein alone.

## Q2. Is the calibrated constant fitted on meta_train only?

**No.** The calibrated-constant REFERENCE (1.3471 pK^2) is the meta_val
target-level mean; it reads meta_val labels and is disclosed as a reference
line, not a legitimate predictor. The meta_train-only constant is 2.1547
(per-episode target) / 2.1703 (canonical target), and the honest comparison
is against that: ESM-150M linear probe 1.7478, ESM-650M linear probe 1.6875,
panel-feature MLP 1.8868. So features do beat the train-only constant, but no
legal predictor approaches the 0.1239 pK^2 level budget required for k=0
MSE <= 1.00 at the measured shape term.

## Q3-Q5. Does target level mix assay history / panel composition / lab or measurement-type bias, and is it a protein property at all?

Measured on meta_train (346 targets, 258 components, 466 documents,
between-target level variance 2.3619 pK^2, 5-fold CV by component):

| covariate family | held-out MSE | held-out R^2 | in-fold R^2 |
|---|---|---|---|
| component identity (one-hot) | 2.3884 | -0.011 | +0.463 |
| document identity (one-hot) | 2.2006 | +0.068 | +0.702 |
| panel composition (ligand-set stats) | 1.7980 | +0.239 | +0.273 |
| protein (ESM-150M pooled + length) | 2.0810 | +0.119 | +0.429 |
| joint | 1.7501 | +0.259 | +0.865 |

Interpretation:
- 46% of the level variance is within-component (protein-family/assay
  confounded) and does NOT transfer across components (held-out R^2 -1.1%).
- 70% is explained in-fold by document (assay) identity, but only 6.8%
  transfers to unseen documents; the double-cold split gives meta_val zero
  document overlap, so raw document identity cannot transfer at all.
- **Panel composition is the best transferring signal (23.9%)** — consistent
  with 'a BindingDB target mean affinity depends on which ligands were
  tested against it'.
- The protein sequence embedding alone transfers 11.9%.

**Conclusion: level is a joint property of protein, assay and the tested
ligand population.** The protein-sequence share that generalizes across
homology components is small but nonzero; the assay/document share is large
in-sample but document-identity itself does not generalize; panel composition
is the strongest legal, generalizing covariate.

## Level identifiability on the frozen meta_val banks (legal inputs only)

Target = episode truth mean (the quantity the level term supervises).
Selection on meta_train component folds; meta_val read once:

| feature family | MLP level MSE | linear level MSE |
|---|---|---|
| grand mean (meta_train) baseline | 2.1547 | — |
| protein ESM-150M pooled | 3.2472 | 1.7478 |
| protein ESM-650M pooled (external, local snapshot) | 3.2927 | 1.6875 |
| panel composition + assay covariates | **1.8868** | 1.9400 |
| protein + panel | 2.0115 | 1.7993 |
| panel, shuffled-panel control | 5.0747 | — |

The shuffled-panel control confirms the panel association is real signal,
not probe capacity. Still, the best legitimate level predictor (1.6875) is
13.6x the 0.1239 budget, and 25% worse than the meta_val-calibrated
reference (1.3471). The Stage C boundary — k=0 MSE <= 1.00 unreachable by
level calibration with available legal inputs — is **confirmed and
strengthened**: it now covers ESM-150M/650M pooled embeddings, panel
composition, assay covariates, and their joints, not only the 150M bank.

## Q4 (occupancy signal) — stratification

Frozen linear probe on within-target centered affinity (Stage C selected
weight decay 1.0), meta_train component folds, meta_val read once:

| stratum | n | mean r |
|---|---|---|
| meta_train held-out folds, scaffold-shared | 137 | +0.236 |
| meta_train held-out folds, scaffold-novel | 93 | +0.154 |
| meta_train held-out folds, recall-high | 87 | +0.173 |
| meta_train held-out folds, recall-low | 143 | +0.221 |
| meta_val overall | 41 panels | +0.203 |

The signal **survives scaffold novelty and low ligand recall**, and is
positive across components on average (per-component r heterogeneous:
several +0.3..+0.8, a few near zero or negative singletons). It remains a
small, real, unexploited shape lever (~0.04 pK^2 if fully exploited).

## Q5b (decisive): does assay history transfer level ACROSS targets?

Within meta_train only, leave-one-target-out transfer inside shared documents
(same DOI, different target): the mean level of a target's document-mates
predicts its level with MSE 1.010 vs the grand-mean 1.912 — R^2 = +0.451
(D0b_DOC_TRANSFER.json; 466 documents, targets sharing documents). So nearly
half of the level variance transfers between different proteins **within one
document**, and the double-cold split (document_overlap = 0) makes exactly this
signal unavailable at inference. Zero-shot cold-target level is therefore
underdetermined by construction: the strongest legitimate level predictor
(assay history) does not transfer to new documents, and the legal inputs that
do transfer (panel composition 23.9%, protein sequence 11.9%) cover at most
about a quarter of the level variance.

## What Stage E therefore tests

Preregistered in PREREGISTRATION.md: a panel-set level head (I1, framework)
trained by orthogonal level/shape routing (I2, training), against the
leak-free T2 baseline and two single-variable ablations. The expected honest
outcome band: level MSE ~1.6-1.9 and centered MSE ~0.83 give k=0 MSE
~2.4-2.7 — i.e., the panel head can only shave a fraction of the level gap,
and the <=1.00 k=0 target stays out of reach unless the level head finds
signal beyond what the D0 probes found.
