# MetaSieve v1 biological-axis development report

Date: 2026-08-11

## Verdict

```text
V1_DEVELOPMENT_REPAIR_NOT_SELECTED
scientific confirmation          false
production migration             false
main-v0 test values used         0
```

The measured-contrast auxiliary objective improves support-label specificity,
but the proposed shared pair prior and residual pair adapter substantially harm
absolute meta-validation performance. Correct protein identity remains
unidentified because using the same wrong protein for both support and query
recovers the correct-arm error.

## Protocol integrity

- optimization used 285 k=5 source tasks in 207 CD-HIT40 clusters;
- cluster-uniform then target-uniform episode sampling for every arm;
- five fixed seeds, k=5, d=2, ridge 1.0, 1,000 steps;
- meta-validation contains 37 k=5 tasks in nine eligible clusters;
- 50 main-v0 test targets and 1,934 test cells are physically blacklisted;
- meta-validation episodes contain support pKi but no query pKi;
- predictions were written before the independent truth file was opened;
- wrong-protein donors are all meta-train proteins, cross-family and matched by
  sequence length/composition; none comes from meta-validation or test.

The source census found legitimate measured supervision: 1,002 within-panel
ligand groups and 1,820 same-panel/same-ligand cross-family partner groups.
V1-B uses measured pKi differences only. No unmeasured pair is labelled as a
non-binder.

## Main results

Target-macro values on development meta-validation:

| Arm | MSE down | R2 up | Pearson up | Spearman up |
|---|---:|---:|---:|---:|
| Ligand d0 | 3.084 | -16.734 | **0.161** | **0.164** |
| Pair d0 | 3.806 | -40.384 | 0.125 | 0.128 |
| V0 correct | **1.800** | **-4.579** | 0.151 | 0.146 |
| V1-A correct | 4.204 | -14.765 | 0.113 | 0.121 |
| V1-B correct | 3.890 | -9.878 | 0.098 | 0.095 |
| V1-B permuted | 4.482 | -10.159 | 0.028 | 0.030 |
| V1-B wrong query | 4.960 | -13.390 | 0.091 | 0.088 |
| V1-B wrong support/query | 3.866 | -11.859 | 0.099 | 0.094 |
| V1-B foreign support | 19.423 | -144.568 | 0.017 | 0.014 |

The full-rank pair prior alone is harmful: pair d0 is worse than ligand d0 by
cluster-macro MSE `1.083` (one-sided LCB for improvement `-2.015`). The residual
pair adapter therefore does not provide a transferable shared biological prior
on the frozen T-BASIS.

## Mechanism contrasts

Values are cluster-macro `MSE(control)-MSE(correct)` with one-sided 95% LCB.

| Model/control | Reduction | LCB | Result |
|---|---:|---:|---:|
| V0 zero | 1.497 | 0.406 | PASS |
| V0 permuted | 0.145 | 0.056 | PASS |
| V0 wrong query | 2.171 | 0.393 | PASS |
| V0 wrong support | 2.470 | 1.174 | PASS |
| V0 wrong support/query | **-0.033** | **-0.090** | FAIL |
| V0 foreign | 4.494 | 2.783 | PASS |
| V1-A permuted | 0.213 | -0.034 | FAIL |
| V1-A wrong query | 1.107 | 0.645 | PASS |
| V1-B permuted | 0.390 | 0.039 | PASS |
| V1-B wrong query | 0.898 | 0.582 | PASS |
| V1-B wrong support/query | **0.066** | **-0.123** | FAIL |
| V1-B foreign | 17.853 | 8.859 | PASS |
| V1-B versus ligand d0 | **-0.234** | **-1.141** | FAIL |

V1-B successfully enlarges and stabilizes the correct-versus-permuted gap, but
it reduces the wrong-query gap relative to V1-A and remains worse than ligand
d0. It fails the frozen candidate criteria.

The 2x2 factorial localizes the deeper defect. Replacing only support protein
or only query protein is destructive, but replacing both with the same wrong
protein restores performance. The solver is exploiting a self-consistent
support/query coordinate system; it is not requiring the biologically correct
protein identity. This explains why support specificity can coexist with weak
partner biology.

## Decision

Do not select V1-B, unfreeze ESM/GINE/bridge, add Q-PMA, connect CSMO, or move
research code into `model/` or `scripts/`. More readout capacity is not the next
repair: both pair d0 and the residual pair adapter worsen generalization.

The next permissible work is a source-only information audit that asks whether
the frozen interaction statistic contains correct-partner information after
controlling for support/query self-consistency, ligand identity, protein length
and family. A new trainable frontend requires a separate preregistration and
real interaction/geometry supervision with a fresh confirmation supply; this
development result does not authorize it.

## Audit artifacts

- result SHA256: `9eb06aab27c0409802ec33c368f5b8248bc0de81ca37f7bab8c40a12c53f5933`
- prediction SHA256: `b70108d2afa72a1aa045af5b99090add1449fcf1142906e39bc2936f8d44576e`
- development seal SHA256:
  `5ee1fc9e857c20c452e5e7f23732425b678cb6deaf8edd283d1b78769dfb1298`
- wrong-feature manifest SHA256:
  `b370d4cb8634664fb59ba8ba3e0f5e284e453692508df3f873c4bb7ec604786a`

