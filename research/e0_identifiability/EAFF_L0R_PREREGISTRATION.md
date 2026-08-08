# E-AFF-L0R Corrected Protein Affinity Location Gate

Status: registered before selection, fitting or scoring. L0R is the corrected
repeat that `EAFF_L0_RESULT.md` prescribed. It inherits every frozen element of
`EAFF_L0_PREREGISTRATION.md` and `EAFF_L0_DATA_CONTRACT.md` unchanged, including
the operator, the anchor ladder, the arms, `margin_L0 = 0.5 * sigma_assay`, the
Ki-only endpoint admission, and every prohibition.

Three things change, each named in the L0 result as a required repair, and each
fixed here **before** any L0R outcome exists.

## Repair 1: fresh panel

L0's `115` tasks are consumed. L0R excludes every task used by P0, H0A, H0C
**and L0**, leaving `3,070` eligible unconsumed Ki tasks in `78` closure
components. L0R takes at most `3` tasks per component, `20` ligands each, for a
target of `195` tasks and `3,900` observations.

Closure components necessarily overlap L0's, because only `78` unconsumed
components remain. L0R is therefore **development evidence, not untouched
validation**, and is recorded as such. No task is reused.

## Repair 2: informative coverage statistic

L0's coverage was containment of the observed step CDF by the emitted band,
which is identically `0.0` for any band narrower than the whole unit interval.
It is replaced by **mean-interval coverage**: the indicator that the observed
`Y` lies inside the band-induced mean interval
`[a_max - integral U, a_max - integral L]`.

This is informative for narrow bands, is monotone in width so it cannot be
gamed by shrinking, and is the coverage notion appropriate to a location claim.
Gate condition 3 keeps its form: `coverage(A3) >= coverage(C) - 0.05`.

## Repair 3: registered positive-control precondition

L0 produced a null from a pipeline that could not detect a ligand effect either.
L0R therefore requires a working positive control **before** any
protein-versus-control contrast is interpreted:

```text
location_error(A0) - location_error(A1) >= 0.1 * sigma_assay,
with a 95% closure-component-bootstrap lower bound above zero.
```

If this fails, L0R stops with

```text
L0R_NOT_RUN_POSITIVE_CONTROL_ABSENT
```

and reports **no** protein verdict. A readout that cannot see ligand information
cannot be asked about protein information. The `0.1` ratio is deliberately far
below the `0.5` Gate margin: it asks only that the pipeline detect *something*.

## Estimator change, declared as a consequence of Repair 3

L0 used a fixed-bandwidth Nadaraya-Watson conditional CDF with bandwidth
`n ** (-1/(d+4))`, which in `d = 7` returns nearly the marginal and is the
likely reason no arm carried signal. L0R uses a **k-nearest-neighbour
conditional CDF** in the same bounded seven-dimensional feature space, with

```text
k = ceil(sqrt(n_train))
```

identical for every arm, declared here and not tuned. This is a
dimension-robust estimator of the same object; it is not an architecture search,
and no alternative will be tried after seeing the result.

## Unchanged

Arms `A0..A4`, capacity matching, identical outer folds, five-fold
cross-fitting by closure component, deranged arm evaluation-only and scored with
`A3`'s fitted estimator, `sigma_assay` estimated from within-assay replicates
before scoring, `band_loss` primary, closure-component macro with bootstrap
intervals, Ki only, and every prohibition in the data contract.

## Gate

After the positive control passes, for each control `C` in `{A1, A2, A4}`, `A3`
passes against `C` only if all three hold:

1. `band_loss(C) - band_loss(A3) > 0`, 95% bootstrap lower bound above zero;
2. `location_error(C) - location_error(A3) >= margin_L0`, lower bound above zero;
3. `coverage(A3) >= coverage(C) - 0.05`.

`A3` passes overall only if it passes against `A1`, `A2` and `A4`
simultaneously.

## Permitted terminal verdicts

- `L0R_NOT_RUN_POSITIVE_CONTROL_ABSENT`
- `L0R_NOT_RUN_NUMERICAL_PRECONDITION_FAILED`
- `PROTEIN_SPECIFIC_AFFINITY_LOCATION_NOT_IDENTIFIED`
- `PROTEIN_SPECIFIC_AFFINITY_LOCATION_IDENTIFIED_IN_SOURCE`

## Scope

A pass establishes only
`PROTEIN_SPECIFIC_AFFINITY_LOCATION_IDENTIFIED_IN_SOURCE`, on development
evidence. It does not authorize biological `z` admission, `model/` promotion,
DAVIS or recipient evaluation, X1, angular or many-body bases, RFSA, theory
modification, P2-P4, an end-to-end DTA claim, or a free-energy interpretation.
No rerun with alternative anchors, widths, losses, margins, seeds or estimators
is permitted after the result is seen.
