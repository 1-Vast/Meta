# DCST-R12 structural interaction-moment prompt decision

Date: 2026-07-28  
Decision: `STOP_PROMPT_ONLY_TEST_PROMPT_ON_STRONG_ZERO_SHOT_BACKBONE`

## Result

On ChEMBL-37 strict dual-cold development:

- B0: Spearman `0.0982`, RMSE `1.456`;
- SIMP-Priv: `0.0966`, RMSE `1.458`;
- SIMP-NoPriv: `0.0989`, RMSE `1.453`;
- SIMP-Uniform: `0.0982`, RMSE `1.455`.

The paired SIMP-Priv minus B0 effect was `-0.0017`, 95% interval
`[-0.0048, 0.0009]`. Privileged minus NoPriv was `-0.0023`,
`[-0.0054, 0.0004]`. Uniform reproduced B0's ordering as designed. Target and
ligand destruction did not show mechanism removal.

Only source admission and RMSE safety passed. Wall time was `207.851 s`; peak
allocated CUDA memory was `1296.6 MiB`. Confirmation and sealed test were not
scored.

## Diagnosis

The frozen structural moment alone is not a sufficient downstream affinity
representation. It deliberately excludes the direct trainable target–ligand
path and asks a 33k head to recover all interaction ordering from a PLINDER
teacher. The result rules out prompt-only sufficiency but does not test whether
the prompt can improve a competent ChEMBL-trained zero-shot interaction
backbone.

Historical BM0/BM1 results are not compatible backbones: they read support
affinity labels for each development target and fail protein-shuffle
specificity. The compatible local backbone is CFRI's zero-shot
ligand-to-protein cross-attention. The next route fuses a trainable CFRI
feature with the frozen privileged moment and uses matched NoPriv and uniform
prompt controls.

