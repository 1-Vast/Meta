# PD-M2A architecture decision (2026-07-25)

## Scope

This is train-only method development on the 12,574 Metz TRAIN cells in 101 frozen target-homology
components. No Metz development label, Davis label, confirmation label, or sealed asset was read.
The result is an architecture-selection result, not evidence of external dual-cold prediction.

## Architecture tested

OSA-ORRC preserved the exact observed-edge main-effect projection and frozen product readout. A
ridge-stabilised bilinear interaction coefficient was followed by the two-parameter monotone map

```text
eta(s) = s * sigmoid((log(s / s1) - log(threshold)) / width).
```

Both global parameters and the ridge were selected inside each outer homology-component fold. The
controls were the unshrunk ridge coefficient, a single soft threshold, fixed hard-rank truncation,
and the current convex fixed-nuclear-norm ORRC.

The inverse-pilot weighted-nuclear-norm candidate was not trained. With singular values ordered
descending, inverse-pilot weights are ordered ascending; that penalty does not have the convexity
certificate claimed in blueprint v3. Treating its weighted SVT fixed point as a certified convex
solution would be a mathematical error.

## Runs and adjustment

The initial frozen ridge grid is recorded in `panel_pdm2.json`. OSA obtained mean held-component
correlation 0.1884 versus 0.1874 for identity, but its paired 95% lower bounds were below zero for
every control and all selected ridges were on the upper grid boundary.

The boundary identified the only justified adjustment: test stronger global ridge rather than add
spectral capacity. `panel_pdm2_adjusted.json` records that run. OSA fell to 0.1896 while identity rose
to 0.1923 and soft thresholding to 0.1938. This adjusted run is explicitly post-result development,
not a second independent test.

`panel_pdm2_final.json` adds the current fixed-nuclear-norm ORRC with the preregistered expanded
penalty grid. Final held-component results were:

| arm | mean correlation | bootstrap LCB95 | bootstrap UCB95 |
|---|---:|---:|---:|
| monotone OSA | 0.1896 | 0.1342 | 0.2443 |
| ridge identity | 0.1923 | 0.1405 | 0.2442 |
| single soft threshold | 0.1938 | 0.1403 | 0.2478 |
| hard rank | 0.1925 | 0.1414 | 0.2441 |
| fixed nuclear norm ORRC | 0.1933 | 0.1394 | 0.2468 |

OSA minus fixed-nuclear-norm ORRC was -0.0037 with paired component interval
[-0.0124, 0.0045]. OSA also selected boundary values and failed the frozen paired-LCB criterion.
One fixed-nuclear-norm outer fold selected the new 0.01 lower boundary, so the spectral grid issue is
not fully closed by merely widening that grid.

## Decision

```text
OSA_ORRC_ARCHITECTURE_FAIL_REVIEW
RETAIN_FIXED_NUCLEAR_ORRC_AS_TRAIN_ONLY_REFERENCE
```

The extra monotone spectral module is rejected. The soft-threshold ridge arm exceeded fixed nuclear
ORRC by only 0.0005 in mean correlation and was not an independently registered replacement test;
that difference does not justify replacing the existing convex, KKT-audited reference. The retained
architecture is therefore exact observed-edge projection plus fixed-nuclear-norm ORRC. This is a
model-development decision, not a successful predictive model claim.

PD-M2 transferable-subspace closure remains incomplete. A document-block sensitivity analysis is
also structurally unavailable on Metz because all rows come from one document. External predictive
success remains unidentified until an open, powered, independent panel exists.
