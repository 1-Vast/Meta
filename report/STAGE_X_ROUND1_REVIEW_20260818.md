# Stage X round-1 independent review

Date: 2026-08-18. Scope: Stage X0 data acquisition, live integrity tests,
representation-capability instrumentation, and the uncommitted I1 draft. This
review does not modify the frozen preregistration or any Stage X implementation.

## Executive decision

Stage X is a justified reopening of Core Task 1. X0-D materially improves the
evidence base: the original Duong-Ly, Anastassiadis, Davis and PKIS2 supplements
are now governed local inputs, and the Duong-Ly surface contains 183 compound
columns, 76 mutant-construct metadata rows and 21 parent accessions. The Davis
matrix preserves blanks rather than silently converting them to exact 10 uM
measurements.

X0 has **not passed**. X0-D acquisition is complete enough to continue and the
seven current tests are useful smoke checks, but I2 is not valid as currently
computed, I1 has neither run nor implemented the preregistered estimand, and
I3-I5 remain incomplete. No X0-P or X1 biological inference is authorized.

## What is established

1. The external-data blocker from the preceding W0-P cycle has been removed for
   development: official supplementary assets are present with URLs, hashes and
   label semantics. Licensing restrictions correctly keep raw/derived panel
   data out of Git.
2. Global mean-pooled ESM is weakly sensitive to the constructed point mutants
   under the current calculation: median WT-mutant distance 0.0204 versus an
   inter-protein scale of 1.9278, ratio 0.0106. This supports excluding that
   representation from interpreting a point-mutation null, subject to the
   sequence-numbering audit below.
3. Composition is likewise insensitive (ratio 0.0188).
4. The Duong-Ly matrix and metadata are not a simple `183 x 76` rectangle:
   Table S2 has 97 assay rows after the CAS row, comprising parent and mutant
   constructs, while Table S1 describes 76 mutant constructs. Parsing and
   power calculations must keep these roles separate.

## Load-bearing defects in I2

### Pair coordinates are inconsistent

`local_window` and `esm_local_window` are precomputed independently for each
row. A mutant row uses its mutation coordinate; the parent WT row has no
mutation annotation and therefore uses the sequence midpoint. The reported
WT-mutant distance consequently compares two different sequence locations.
Ratios 1.0 and 1.323 do not demonstrate mutation sensitivity.

The corrected computation must be pair-conditioned: for each WT-mutant pair,
extract both WT and mutant representations at the same verified mutation
coordinate and with the same mask/window contract.

### Mutation-token denominator is degenerate

Parent rows receive a zero mutation token, so every parent-parent distance is
zero. Dividing 1.414 by the `1e-12` floor creates the reported
`1.41e12` ratio. This is an edit descriptor with an undefined inter-protein
scale, not an admitted protein representation. It must not count toward the
three-representation pass requirement.

### Mutation numbering is not fully governed

The mutation application routine does not assert that the reference residue at
the supplied coordinate equals the annotated old residue. A direct audit finds
72 of 76 point mutations match the downloaded sequence and four do not:
`BRAF(V599E)` and three `PDGFRalpha` variants. These may reflect historical
numbering, isoforms, construct offsets or a wrong accession, but must be mapped
explicitly before representation extraction. Silent substitution at a
mismatched coordinate is invalid.

### Reproducibility and coverage gaps

The random control is seeded with Python `hash(name)`, which changes across
processes unless `PYTHONHASHSEED` is fixed. Use a stable cryptographic hash.
The preregistration names a KLIFS aligned-pocket representation, but the current
instrument does not implement it. Multi-mutation/deletion constructs, mutations
outside the 1020-token ESM truncation and unknown parse forms need explicit
admission/exclusion counts.

## I6 scope

The seven tests verify file presence, a few label anchors and the frozen hash,
but do not yet satisfy the complete I6 contract. The CSC test exercises a local
toy function only with `reference=0`; it does not test a production CSC
implementation or reference-sign reversal. There are no current assertions for
dead regularizers, gradient coverage or the information destruction of each
permutation control. The suite should be described as seven initial integrity
smokes, not as complete I6 qualification.

## I1 draft review

The uncommitted `x0_planted.py` currently fails at its `einsum` call because
the equation specifies four operands and supplies three. More importantly,
repairing that line alone would not make the test valid:

- generated protein and ligand main effects are unused;
- the planted interaction is added to the real percent-remaining matrix, so
  unknown biological structure remains a confound;
- training and evaluation use the same rows;
- the model predicts the raw endpoint, but raw predictions are scored directly
  against the planted interaction without removing learned main effects;
- the ligand-only arm is an untrained zero vector rather than a matched model;
- sign accuracy uses raw predictions, whose endpoint offset can dominate sign;
- the same invalid local representation from I2 is reused.

I1 must be rebuilt around an exactly known synthetic label decomposition on the
real observation/censoring graph, with component-held-out evaluation, matched
arms, extraction of the fitted interaction component and recovery measured
against planted interaction truth.

## Required next order

1. Preserve the frozen X0 preregistration and record this as an implementation
   correction, not a threshold change.
2. Correct I2 pair-conditioned extraction, mutation numbering, stable seeds,
   KLIFS coverage and exclusions; rerun I2 before claiming a pass.
3. Replace the I1 draft with a test that implements the frozen estimand and add
   regression tests that fail on the current confounded construction.
4. Expand I6 to the full preregistered contract.
5. Implement I3, I4 and I5 only after the shared identifiers and censoring
   representation are governed.
6. Produce `X0_RESULT.json` and `X0_REPORT.md`. Only an all-instrument X0 pass
   authorizes a separately frozen X0-P preregistration.

## Scientific status

The round is promising because it acquired the right class of matched-variant
data and confirmed why global pooling is a poor mutation instrument. It has not
yet shown that local representations are capable, that the diagnostic pipeline
can recover a planted interaction, or that any real protein-conditioned signal
exists. The correct current label is **X0 active; data acquired; instrumentation
requires correction**.
