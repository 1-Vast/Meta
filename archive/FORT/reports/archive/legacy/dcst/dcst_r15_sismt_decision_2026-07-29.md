# DCST-R15 SISMT decision

Date: 2026-07-29  
Decision: `STOP_SISMT_LABEL_BLIND_SUPPORT_GATE`

## Result

SISMT used 1,613 PLINDER pairs from 767 exact targets and 4,472 ChEMBL-train
pairs from 559 targets. No ChEMBL affinity column was loaded.

The privileged representation had 143 generalized eigendirections in the
frozen support-ratio interval `[0.25, 4.0]`. Only one also met the R6
certificate-overlap threshold:

- generalized eigenvalue: `0.619565`;
- privileged mechanism overlap: `0.050474`;
- frozen threshold: `0.05`.

That borderline direction was absent from every target-block bootstrap
retained subspace: its median squared projection over 20 repetitions was
exactly `0.0`, below the frozen `0.50` stability floor. The final privileged
retained dimension was therefore zero.

NoPriv, uniform, matched random, wrong-target, and wrong-ligand controls also
retained zero stable directions. With no valid projector, partial transport
and Stage-2 affinity fitting were not run. Runtime was `15.069 s` on CUDA;
peak allocated memory was about `982 MiB`.

## Interpretation

R14 showed that the full privileged mechanism is responsive but off support.
R15 now rules out the narrower rescue hypothesis that a stable,
destruction-certified covariance direction survives in the current
source/target intersection. Lowering the overlap or stability threshold would
select an outcome-dependent, bootstrap-unstable direction and is forbidden.

This stops SISMT for the current PLINDER source. It does not affect the
separately preregistered DTIOD hypothesis, which changes the transferred
information object from an absolute mechanism state to a local mixed
finite-difference response.
