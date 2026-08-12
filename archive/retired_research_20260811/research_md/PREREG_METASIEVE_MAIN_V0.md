# MetaSieve-main v0 preregistration

Date frozen: 2026-08-11

## Scope

This protocol tests whether a frozen biological pair statistic supports useful
protein-task adaptation through a source-learned, support-identifiable section.
It is a method-level literature-aligned MetaSieve experiment, not an exact
AdaMBind split reproduction and not the strict dependency-closed confirmation.

The strict 29-target result does not gate this experiment. Davis, KIBA,
recipient, the O1 strict evaluation labels and external confirmation remain
closed. Only the open BindingDB Articles 202608 exact-Ki development artifact is
used.

## Numeric corpus

The cleaning contract references CARA but retains the MetaSieve Ki estimand:

1. require one protein chain, one unambiguous sequence, canonicalizable SMILES,
   ligand molecular weight at most 1000 Da, at most 128 graph atoms and an
   implicit-hydrogen canonical graph compatible with the frozen T-BASIS atom
   ordering;
2. retain exact positive uncensored Ki only and use the existing
   `pKi = 9 - log10(Ki[nM])` extractor;
3. treat the existing `document + Ki + normalized protocol` panel ID as the
   BindingDB assay proxy;
4. take the median of replicates within `(panel, protein, ligand)`;
5. to obtain one observation per AdaMBind-style protein task, take the median
   of panel medians for repeated `(protein, ligand)` pairs. Panels receive equal
   weight, so prolific replicated panels do not dominate;
6. exclude identity conflicts rather than selecting a value by model outcome.

The cross-panel median is a declared MetaSieve reconciliation, not attributed
to CARA or AdaMBind. A later sensitivity analysis must compare panel-specific
observations without cross-panel aggregation.

## Task and split

One unique protein sequence is one task. Targets with at least four distinct
ligands enter split construction; legal support sizes are those satisfying
`n >= k + 3` for `k in {1,2,3,5}`. The v0 headline trains and evaluates k=5,
therefore only targets with at least eight observations contribute episodes.

Run the pinned CD-HIT 4.8.1 binary as:

```text
cd-hit -i proteins.fasta -o clusters -c 0.40 -n 2 -d 0 -T 1 -M 0
```

Complete clusters are assigned 8:1:1 to meta-train/meta-validation/meta-test by
a deterministic target-count balancing algorithm with split seed `20260811`.
The cluster file, command, binary hash, assignment and corpus hashes are saved.
Affinity values and model scores do not participate in assignment.

## Frozen biological v0

The first real v0 uses the frozen, geometry-validated 288D T-BASIS as `phi`.
The learned objects are the ligand population predictor, the orthonormal
`U in R^(288 x d)` and population coordinate coefficient `w0`. The support
section is the positive-ridge closed form. Query loss differentiates through
the solve into `U`, but not into the frozen ESM/GINE/mechanism frontend.

This is not claimed to be end-to-end biological-coordinate training. Unfreezing
that frontend requires a differentiable replacement for the current NumPy
T-BASIS aggregation and a separate stability Gate.

Candidate `d=1..5` and ridge values `{0.01, 0.1, 1.0}` are selected on
meta-validation target-macro MSE, ties favoring smaller d then larger ridge.
The d=0 population arm is never used for hyperparameter selection.

## Episodes and inference

- k=5 support and at least three disjoint query ligands;
- five frozen training seeds: 20260811 through 20260815;
- target-uniform episode sampling;
- five fixed support draws per test target, with identical draws across arms;
- query labels are evaluator-only and never enter adaptation;
- target-macro point estimates and paired target bootstrap confidence intervals.

## Arms and primary Gates

Train capacity-matched full and ligand-only sections. Evaluate the same full
checkpoint with correct, zero, foreign-target and cyclically permuted support,
and with the frozen wrong-protein feature arm.

Primary paired squared-loss reductions at k=5 are:

```text
M1: full correct > d=0
M2: full correct > zero, foreign and permuted support
M3: full correct > ligand-only and wrong-protein
```

Each Gate requires a one-sided 95% target-bootstrap lower bound above zero after
averaging draws and seeds within target. Report MSE, RMSE, R2, CI, Spearman and
Pearson target-macro where defined. Undefined correlations remain NA.

## CSMO boundary

The production 28D z interface has no admitted biological coordinate map.
Therefore CSMO/band loss is `NA_NOT_ADMITTED` in v0 rather than populated by an
invented state. A point-level PASS authorizes preregistration of the compact z
bridge; it does not establish law-valued prediction.

## Verdicts

```text
REAL_BIOLOGICAL_META_SECTION_V0_PASS
META_EFFECT_NOT_IDENTIFIED
SUPPORT_SPECIFICITY_NOT_IDENTIFIED
BIOLOGICAL_SPECIFICITY_NOT_IDENTIFIED
```

No v0 outcome authorizes tuning on the strict confirmation set.
