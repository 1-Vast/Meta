# Stage X0 governance verdict: INVALID INSTRUMENT

Date: 2026-08-18. Authority: user directive (round 2). This document does not
modify any frozen artifact.

## Verdict

The original Stage X0 (`STAGE_X0_PREREGISTRATION.md`, SHA-256
`03cdc907df3e778f5fe79fb1a238d35ebb6ece5e9e743db181728ba6b25e9683`) is
formally ruled **INVALID INSTRUMENT** — not a biological PASS and not a
biological FAIL.

## Reason

The round-1 review and this round's audit established that the I2 distance
ratio r_pair = median_pair_distance / median_inter_protein_distance cannot
certify representation capability: the ratio is sensitive to the edit
descriptor, the window contract and the denominator construction, but a pass
does not imply the representation contains readable, transferable
variant-relevant information, and a fail does not imply incapability (the
global pooled ESM fail is informative, but the local one-hot pass is vacuous —
a single substitution changes exactly two one-hot positions). Replacing the
distance ratio with probe/control-task selectivity changes the load-bearing
estimand; that is a measurement redefinition, not a bug fix, and cannot be
retroactively applied to the frozen preregistration.

## What is a bug vs a measurement-definition failure

Implementation bugs (repaired in the successor, evidence preserved):

| # | Defect | Class |
|---|---|---|
| 1 | WT local window at sequence midpoint vs mutant window at mutation coordinate | implementation bug (coordinate contract) |
| 2 | same coordinate defect in ESM local windows | implementation bug |
| 3 | mutation_token parent-parent denominator = 0 -> ratio ~1.41e12 | implementation bug (degenerate denominator) |
| 4 | mutation_token counted as an admissible representation | definition error (edit descriptor is pair-conditioned, not a protein representation) |
| 5 | BRAF V599E and PDGFRalpha variants mismatched the reference (historical numbering + wrong-species accession Q9DE49) | data-mapping defect |
| 6 | Python hash() for the random control | implementation bug (process-dependent seed) |
| 7 | KLIFS not implemented | implementation gap |
| 8 | no admission rules for multi/deletion/insertion/truncation/fusion/unknown/long-sequence | definition gap |
| 9 | I1 draft einsum failure; main effects unused; real endpoints mixed into planted truth; train=eval; raw endpoint vs interaction truth; zero-vector ligand-only; split without parent/scaffold isolation | implementation + statistical-design bugs |

Measurement-definition failures (cannot be repaired in place):

| # | Defect | Consequence |
|---|---|---|
| M1 | distance ratio as the load-bearing capability gate | gate invalid as a capability certificate -> X0 INVALID INSTRUMENT |
| M2 | "tau in {0.2,0.4,0.8,1.6} log units" unanchored to noise SD and silent on rank/locality | planted grid under-specified; superseded by the successor's tau* x rank x locality grid |
| M3 | no ANOVA/projection operator shared between truth and prediction | interaction recovery could be conflated with main-effect removal |

## Disposition of artifacts

- `STAGE_X0_PREREGISTRATION.md`, `X0_PREREGISTRATION_SHA256.txt`,
  `X0_INSTRUMENTS.json`, `x0_instruments.py`, `x0_planted.py`,
  `X0_PLANTED.json`, `x0_instruments.log`, `x0_planted.log` are preserved
  unchanged as negative evidence.
- Round-2 correction artifacts (`x0_pair_table.py`, `X0_PAIR_TABLE.json`,
  `x0_i2.py`, `X0_I2.json`, `X0_I2_REPORT.md`, `tests/test_x0_corrections.py`)
  are retained as diagnostic-QC evidence and as the seed of the successor's
  Q0-B/Q1 layers; their PASS/FAIL labels are NOT load-bearing for X0.
- Corrected successor:
  `tools/research/stageX_csc_signal/stageX0c_measurement_qualification_20260818/`
  with its own frozen preregistration (see PREREGISTRATION.md + SHA there).

## Successor gates (order enforced)

Q0 variant-coordinate layer -> Q1 representation capability (probe
selectivity) -> Q2 fully synthetic planted harness -> Q3 biological panel
qualification (Saifudeen 2026) -> B1 same-study WT->mutant positive control ->
B2 localization -> C cold-protein interaction test -> D final DTA. No gate may
be entered before its predecessors pass. Q0 failure => UNRESOLVED(mapping);
Q1 failure => representation incapability reported as such; Q2 failure => the
pipeline may not interpret any real-data negative; Q3 informs pairability only.
