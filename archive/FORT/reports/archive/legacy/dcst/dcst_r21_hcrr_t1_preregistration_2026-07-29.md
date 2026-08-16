# DCST-R21 HCRR Stage-1 preregistration

Date: 2026-07-29  
Status: frozen before loading TRAIN affinity

## Candidate

`High-Confidence Residual-Rank pretraining` (HCRR) is the minimal
implementation of the user's proposed two-stage design:

1. Stage 1 learns the target-ligand residual `g(t,l)` only from the R19
   high-confidence ChEMBL-TRAIN pairs;
2. Stage 2, if authorized, fine-tunes the same interaction branch on all
   ChEMBL TRAIN pairs before development-only evaluation.

R21-T1 tests Stage 1 only. It does not reuse the stopped exact-rectangle
objective and does not add a new attention/router. It uses the existing CFRI
interaction architecture so any change is attributable to evidence selection
and curriculum.

## TRAIN-only split

Use the R19 high-confidence definition:

- at least two replicate records;
- finite replicate SD at most `0.30` pK;
- at least two assays or two documents.

Assign entire ChEMBL homology components to Stage-1 fit (80%) or audit (20%)
by SHA-256 of seed `1729` and component ID. No component may occur in both.
Targets with fewer than four eligible audit ligands do not contribute a
Spearman unit.

The ligand-only B0 residual target must be out-of-fold by ChEMBL homology
component. B0 may use all TRAIN rows in the other folds but cannot train on
the audited component whose residual it predicts.

## Model and loss

Reuse the existing CFRI `Interaction` module and frozen ESM-2 target features.
Train only on eligible Stage-1 fit rows using the already-defined
cross-fitted residual-rank objective:

- rank reversal against out-of-fold B0;
- residual affinity MSE;
- within-target centering;
- B0 orthogonality.

Run all trainable tensor work on CUDA in the `drug` environment.

## Frozen controls

On the held-component high-confidence audit:

- B0: out-of-fold ligand-only prediction;
- HCRR: correct target and ligand pairing;
- HCRR-Tshuffle: same trained model, deterministically deranged target
  features;
- HCRR-Lshuffle: same trained model, cyclically deranged ligands within each
  target episode.

No development, confirmation, or sealed row may be loaded.

## Frozen gate

R21-T1 passes only if:

1. fit side has at least 10,000 rows, 200 targets, and 180 homology
   components; audit has at least 2,000 rows, 50 targets with four or more
   ligands, and 40 homology components;
2. HCRR minus B0 macro within-target Spearman is at least the existing frozen
   ChEMBL MDE (`0.0586`) and its target-block bootstrap 95% lower bound is
   positive;
3. HCRR beats both Tshuffle and Lshuffle with positive 95% lower bounds;
4. each destruction removes at least 70% of the HCRR-minus-B0 effect;
5. HCRR macro RMSE is at most 1.02 times B0 RMSE;
6. zero development/confirmation/sealed rows are loaded.

Failure returns `STOP_HCRR_STAGE1_NOT_IDENTIFIED`.

Pass returns `REQUEST_HCRR_STAGE2_PREREGISTRATION`; it does not itself
authorize development scoring. A Stage-2 preregistration must match total
optimization steps against a one-stage CFRI control and retain the target and
ligand destructions.
