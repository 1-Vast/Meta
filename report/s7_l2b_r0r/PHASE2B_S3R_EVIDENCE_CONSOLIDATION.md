# Phase 2B S3R evidence consolidation

## Terminal result

```text
S2R synthetic ordinal trainability ........ PASS
S3R real structural transfer .............. REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED
heldout-B .................................. NOT OPENED
R6 amplitude/B5 integration ............... NOT RUN
affinity value reads ....................... 0
```

S2R showed that the gauge-free direct-`W` ordinal estimator is trainable on
fresh synthetic closure components. The sealed seed reached held-out
component-macro `AP_bidir = 0.6620`, so the old factorized optimizer defect is
not a valid explanation for S3R.

S3R then trained the same estimator once on the real MONN residue-differential
labels. The trainable object was one `1280 x 41` matrix, constrained to unit
Frobenius norm, with the frozen ESM2 residue states and frozen mean-pooled 41-D
ligand atom features. Training used 210 fixed updates and the frozen hierarchical
sampler. Candidate, repeat and permuted-label learners used the same stream.

## Primary evidence

The primary panel contains 46,818 unordered ligand pairs across 112 closure
components. Every arm used the same pair, construct and component map.

| arm | component-macro AP_bidir |
|---|---:|
| candidate | 0.035880 |
| zero-`W` chance | 0.025472 |
| frozen B5 differential | 0.031582 |
| foreign ligand pair | 0.035735 |
| residue-context corruption | 0.032336 |
| trained permuted-label learner | 0.037125 |

| Gate | observed delta | one-sided LCB95 | required | result |
|---|---:|---:|---:|:---:|
| R1 candidate - chance | +0.010408 | +0.006920 | +0.05 | FAIL |
| R2 candidate - B5 | +0.004298 | -0.001630 | +0.03 | FAIL |
| R3 candidate - foreign pair | +0.000145 | -0.002651 | +0.03 | FAIL |
| R4 candidate - context corruption | +0.003544 | -0.003156 | +0.03 | FAIL |
| R5 candidate - permuted learner | -0.001245 | -0.006595 | +0.05 | FAIL |

R1 is statistically above chance but far below the preregistered practical
margin. It does not replicate beyond B5, foreign ligands, contextual corruption
or a capacity-matched permuted-label learner. The earliest-failure verdict is
therefore `REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED`.

## What the failure means

This is not a flat-data verdict. Phase 2A established real same-construct,
scaffold-distinct ligand conditionality in the MONN residue masks. It is also
not a numerical or participation failure: gradients, movement, unit norm,
variance, zero-`W` chance, stream equality and exact repeat prediction replay all
passed.

The result is scoped to the current measurement basis and estimator:

```text
frozen ESM2 residue states
  x mean-pooled 41-D ligand atom features
  x one gauge-free direct-W ordinal head
```

The strongest current interpretation is that mean pooling removes ligand
topology and atom-local correspondence needed to recover the real
ligand-conditioned residue direction under closure shift. It does not close all
sequence-plus-2D models and does not establish that a larger network would fix
the problem.

## Governance notes

The parent S3R preregistration's R6 was superseded before execution. Pairwise
ordinal training identifies neither absolute amplitude, ligand-feature origin,
nor directions outside the ligand-difference span; adding its raw output to B5
would therefore be arbitrary. R6 remains unopened.

`PHASE2B_S3R_FAIL_CLOSED.json` records a duplicate prepare invocation rejected
by no-clobber after the valid prepare invocation had already completed. The raw
event is preserved and adjudicated in
`PHASE2B_S3R_ORCHESTRATION_ADJUDICATION.json`; it did not modify the frozen
manifests or control the scientific verdict.

`unit_norm_pass` was serialized as the string `"True"` in one machine artifact.
The underlying norm is `1.0000000116`, the aggregate participation field is a
Boolean PASS, and the defect does not affect any Gate. The artifact is preserved
rather than rewritten.

## Remaining boundary

No biological statistic is admitted to `z`. Affinity, selectivity, few-shot
sectioning, heldout-B, R6, DAVIS/KIBA/recipient labels and the frozen operator
`A(F,z)=K(B(z)F(z))` remain untouched.
