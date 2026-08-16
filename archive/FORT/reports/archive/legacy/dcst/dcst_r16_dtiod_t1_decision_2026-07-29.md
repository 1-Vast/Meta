# DCST-R16 DTIOD T1 decision

Date: 2026-07-29  
Decision: `STOP_DTIOD_T1_NO_PRIVILEGED_TANGENT`

## Result

T1 evaluated 220 PLINDER development pairs across 84 exact targets using four
active-Morgan interventions per pair. All seven operator arms ran on CUDA;
no ChEMBL affinity column was loaded.

The privileged local mixed response was not numerical noise. Its
target-macro median RMS was `2.45251e-4`, compared with:

- NoPriv: `3.78443e-8` (`6480.53x` smaller);
- matched random teacher: `1.72029e-7` (`1425.63x` smaller);
- random segment masks: `1.15093e-4` (`2.1309x` smaller);
- random active-bit masks: `1.90654e-5` (`12.8637x` smaller).

Wrong-target destruction reduced the response by `80.13%`. Wrong-ligand
destruction did not reduce it: the wrong-ligand RMS was `2.71960e-4`, which
was `10.89%` larger than the true-pair response. Consequently T1c failed
while T1a, T1b, and the target-count gate passed.

Runtime was `5.913 s`; peak allocated CUDA memory was `627.1 MiB`.

The machine-readable result contains the complete 84-row `target_blocks`
array. The same rows are provided as
`dcst_r16_dtiod_t1_target_blocks_seed1729.csv`; each row is one exact target
block with pair count and target-macro RMS for every T1 arm. The CSV medians
reproduce the JSON `target_macro_median` values exactly.

## Interpretation

The frozen R6 privileged teacher contains a highly non-additive local protein
response that is specific to structurally supervised parameters and
contact-enriched segments. It is not identified as a ligand-specific
interaction operator: pairing the target with a different ligand strengthens
the finite difference rather than destroying it.

Training a student on this operator would faithfully distill the wrong
estimand. T2, T3, and Stage-2 affinity fitting are therefore not run.
Reopening DTIOD requires a new Stage-1 source or supervision that directly
anchors ligand substructures to contacts; changing the destruction threshold
or choosing masks after this result is forbidden.
