# DCST-R13 privileged-prompted CFRI decision

Date: 2026-07-28  
Decision: `STOP_RAW_PRIVILEGED_PROMPT_RUN_LABEL_BLIND_TRANSPORT_AUDIT`

## Result

On the current ChEMBL-37 exact Ki/Kd strict dual-cold registry:

- B0: Spearman `0.098237`, RMSE `1.455993`;
- PP-CFRI-Priv: `0.098237`, RMSE `1.456128`;
- PP-CFRI-NoPriv: `0.098516`, RMSE `1.457321`;
- PP-CFRI-Uniform: `0.097968`, RMSE `1.455980`.

The grouped PP-CFRI-Priv minus B0 effect was exactly `0.0000` with interval
`[0.0000, 0.0000]` over 181 bootstrap units. Privileged minus NoPriv was
`-0.0004 [-0.0010, 0.0002]`; privileged minus Uniform was
`+0.0003 [0.0001, 0.0005]`, far below the frozen `0.0586` MDE.

Both target and ligand destruction reproduced the privileged prediction
exactly. The privileged interaction correction had mean absolute magnitude
`0.001869`, so Stage 2 effectively suppressed the raw Stage-1 prompt.
Only source admission and RMSE safety passed. Wall time was `254.996 s`;
peak allocated CUDA memory was `1303.9 MiB`. Confirmation and sealed test
were not scored.

## Diagnosis and decision

R12 showed that the frozen structural moment is not sufficient by itself.
R13 now shows that adding a fully trainable zero-shot direct interaction
branch does not make the unmodified moment useful. This is not evidence that
the R6 source mechanism is false: R6 remains the only source arm with a
privileged-specific held-source certificate. It is evidence that the source
mechanism is not arriving as a usable downstream reordering coordinate.

No further fusion-head variant is justified before measuring transport
support. R14 therefore audits, without downstream affinity labels, whether
PLINDER and ChEMBL overlap on the target, ligand, and frozen structural-moment
axes. Its frozen route selection determines whether the next model uses
label-blind source importance weighting, a train-domain counterfactual
centering module, or a different Stage-1 source.

