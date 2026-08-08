# E-AFF-L0 Protein Affinity Location Gate Result

## Terminal Verdict

```text
L0_NOT_RUN_NUMERICAL_PRECONDITION_FAILED
```

The run executed end to end and produced arm scores, but one of the three
registered gate conditions was computed from a degenerate statistic. Under the
registered fail-closed rule a numerical or contract prerequisite failure returns
a NOT-RUN verdict rather than a scientific failure. **Claim A remains untested.**

The executed report retains its raw computed verdict for the record; the
terminal verdict is the one above, rendered by `audit_eaff_l0.py`.

## Preconditions That Passed

| Stage | Verdict |
|---|---|
| operator and anchor contract | `L0_OPERATOR_AND_ANCHOR_CONTRACT_FROZEN` |
| location estimand | `L0_LOCATION_ESTIMAND_IDENTIFIED_KI` |
| evidence audits (R0, X0-FEAS, X0-B) | all `POSTRUN_AUDIT_PASS` |

Ki was admitted and Kd was excluded (`C1 = 10 < 30`), so L0 ran on Ki alone.

`sigma_assay` was estimated **before** any arm was scored, from `4,261`
replicate cells and `4,840` degrees of freedom:

```text
sigma_assay = 0.47971  [0.47034, 0.48946]   log units
margin_L0   = 0.5 * sigma_assay = 0.23985
```

That value is a useful independent by-product: the governed ChEMBL37 Ki
within-assay replicate noise is close to half a log unit, which is the scale any
future affinity claim must exceed.

## Executed Design

`115` tasks from `115` distinct closure components, `20` ligands each, `2,300`
observations, five-fold cross-fitting by closure component so evaluation
proteins never share a component with their fitting pool. `355` tasks consumed
by P0, H0A and H0C were excluded before selection. All five arms used identical
folds, the same conditional-CDF estimator, the same declared bandwidth rule and
the same seven-dimensional bounded input.

| Arm | Band loss | Location error (log units) | Coverage | Mean width | `p_0` mass |
|---|---:|---:|---:|---:|---:|
| A0 population | 0.25662 | 1.24758 | **0.0000** | 0.0575 | 1.0000 |
| A1 ligand-only | 0.23153 | 1.26684 | **0.0000** | 0.0956 | 0.0000 |
| A2 sequence-only | 0.23099 | 1.26647 | **0.0000** | 0.0958 | 0.0003 |
| A3 correct + geometry | 0.23107 | 1.26346 | **0.0000** | 0.0960 | 0.0000 |
| A4 deranged + geometry | 0.23108 | 1.26351 | **0.0000** | 0.0960 | 0.0000 |

## The Defect

Gate condition 3 compares empirical coverage between arms. Coverage was
registered as containment of the observed step CDF by the emitted band. On the
fixed `33`-point mesh a step function jumps from `0` to `1` between adjacent
grid points, so a band can contain it only by spanning nearly the whole unit
interval. Every arm returns exactly `0.0`, including the population band.

The statistic is therefore identically zero and carries no information, so one
of the three registered gate conditions did not execute as specified. This is
decidable from the mean interval widths alone and does **not** depend on any
arm's performance, which is why it is admissible as a precondition failure
rather than a post-hoc reinterpretation of a null.

## Supporting Diagnostic, Not A Registered Criterion

Ligand-only (A1) did **not** improve on population-only (A0) in location error:
the gain is `-0.01926` log units. The pipeline therefore never demonstrated that
it can detect any affinity information at all. A null protein result from a
readout with no working positive control is uninterpretable independently of the
coverage defect. This was not registered as a criterion and is recorded as a
diagnostic only.

It is consistent with prior evidence: in H0C the global ligand prior, fitted on
`147k` rows with the full `128`-dimensional ligand state, reached a
component-macro CI of only `0.55487`. A seven-component projection fitted on
`~1,840` rows carrying no power is unsurprising.

## What This Does And Does Not Mean

It does **not** mean the correct protein carries no affinity-location
information. That question is untouched. The contrasts

```text
A3 - A1  band loss +0.00046 [-0.00055, 0.00146]   location +0.00338 [-0.00411, 0.01074]
A3 - A2  band loss -0.00008 [-0.00133, 0.00125]   location +0.00301 [-0.00655, 0.01307]
A3 - A4  band loss +0.00001 [-0.00001, 0.00003]   location +0.00005 [-0.00011, 0.00020]
```

are all far inside their intervals and far below `margin_L0 = 0.23985`, but they
were produced by a design that could not detect a ligand effect either, so they
carry no evidential weight about proteins.

It does **not** overturn or confirm any historical verdict, and it does not
change the project state: `PROTEIN_SPECIFIC_AFFINITY_LOCATION_NOT_YET_TESTED`
stands.

## Consumed Panel

The `115` tasks and `115` closure components used here are now **consumed**.
They must not be reused as untouched validation by any future stage, exactly as
the P0, H0A and H0C panels are treated.

## No Rerun

The registration forbids rerunning with alternative anchors, widths, losses,
margins, seeds or architectures after seeing the result. Repairing the coverage
statistic and repeating on this panel would be precisely the post-hoc selection
the contract prohibits.

A corrected Gate requires a **new registration** and a **fresh unconsumed
panel**. Two repairs are indicated by this run and must be registered in
advance, not chosen after seeing outcomes:

1. Replace step-containment coverage with a statistic that is informative for
   narrow bands, for example the mean band violation already inside `band_loss`,
   or interval coverage of a calibrated central region.
2. Register an explicit **positive-control precondition**: the ligand-only arm
   must beat the population arm by a preregistered amount before any
   protein-versus-control contrast is interpretable. The present run shows this
   precondition is not automatic and must be enforced.

Neither repair may be applied to this panel.
