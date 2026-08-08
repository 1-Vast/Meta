# E-AFF-L0R Corrected Protein Affinity Location Gate Result

## Terminal Verdict

```text
L0R_NOT_RUN_POSITIVE_CONTROL_ABSENT
```

The registered positive-control precondition failed, so **no protein verdict was
computed or is reported**. Claim A remains untested.

## The Two Repairs Worked

L0R was registered to fix exactly the two defects `EAFF_L0_RESULT.md` named.
Both fixes behaved as intended.

**Repair 2, informative coverage.** Mean-interval coverage replaced
step-containment coverage. It now varies across arms — `0.0674` for the
population arm up to `0.1664` for ligand-only — instead of being identically
`0.0`. Gate condition 3 is now computable.

**Repair 3, positive control.** The precondition caught a pipeline that is not
yet able to support a protein claim, and stopped the Gate before any protein
contrast could be misread as evidence. That is the whole point of registering it.

The estimator change from a fixed-bandwidth kernel to `k`-nearest neighbours
(`k = ceil(sqrt(n_train)) = 56`) also worked: ligand-only now moves in the right
direction against the population arm, which it did not in L0.

## Executed Design

`195` fresh tasks across `78` closure components, at most `3` tasks per
component, `20` ligands each, `3,900` observations. All `470` tasks consumed by
P0, H0A, H0C and L0 were excluded. Five-fold cross-fitting by closure component,
identical folds and estimator for every arm, seven-dimensional bounded inputs
throughout.

| Arm | Band loss | Location error (log units) | Coverage | Mean width | `p_0` mass |
|---|---:|---:|---:|---:|---:|
| A0 population | 0.24210 | 1.12972 | 0.0674 | 0.0456 | 1.0000 |
| A1 ligand-only | 0.21298 | 1.09552 | 0.1664 | 0.0892 | 0.0104 |
| A2 sequence-only | 0.22674 | 1.18342 | 0.1579 | 0.0899 | 0.0095 |
| A3 correct + geometry | 0.22121 | 1.15148 | 0.1590 | 0.0901 | 0.0151 |
| A4 deranged + geometry | 0.22119 | 1.15165 | 0.1500 | 0.0916 | 0.0139 |

## Why It Stopped

```text
positive control  =  location_error(A0) - location_error(A1)
                  =  +0.03421 log units,  95% CI [-0.03304, 0.10793]
required          =  0.1 * sigma_assay = 0.04797
```

The point estimate has the right sign, but it is below the required threshold and
its closure-component bootstrap interval contains zero. The pipeline therefore
cannot demonstrate, at this panel size, that it detects **any** affinity-location
information — not even from the ligand, which is the easiest available signal.

Under the registration, protein contrasts computed on such a readout are not
interpretable, and none is reported as a finding.

## The Useful Measurement

L0R did produce one solid, reusable number. With `sigma_assay = 0.47971` log
units, the best available cross-component location signal in governed ChEMBL37
Ki — ligand identity itself — is worth about `+0.034` log units over a population
band, roughly **7% of assay noise**, and is not separable from zero across `78`
closure components.

This bounds what any L0-style location Gate can achieve on this corpus at this
scale. It is not a statement about proteins; it is a statement about the power of
the design. A future Gate must either enlarge the panel substantially, condition
within assay strata to remove between-stratum location variance, or accept that
cross-component absolute-location prediction is not the discriminating estimand
here.

## What Is Not Claimed

The arm table shows `A3` below `A1` and `A3` indistinguishable from `A4`
(`+0.00017` log units, CI `[-0.05078, 0.05102]`). **These numbers are not
evidence about proteins** and are recorded only for completeness, because they
come from a readout that failed its positive control. Reporting them as a
protein finding is exactly what the precondition exists to prevent.

`PROTEIN_SPECIFIC_AFFINITY_LOCATION_NOT_YET_TESTED` stands. No historical verdict
is overturned or confirmed, nothing is admitted to `z`, and DAVIS, recipient
labels, X1, X2, angular work, RFSA, theory changes and P2-P4 remain frozen.

## Consumed Panel

The `195` tasks are consumed. Closure components overlap L0's, since only `78`
unconsumed components remained, so L0R is **development evidence, not untouched
validation**. No task was reused.

## No Rerun

No rerun with alternative anchors, widths, losses, margins, seeds or estimators
is permitted now that the result is seen. A further attempt requires a new
registration, and on this evidence it should change the **design power or the
estimand**, not the model.
