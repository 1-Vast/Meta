# DCST-R6 source-gate decision

Date: 2026-07-28  
Decision: `PASS_SOURCE_AUTHORIZE_STAGE2_DEVELOPMENT`

## Result

The preregistered source-only artifact
`dcst_r6_stage1_seed1729.json` passes all three frozen R6 gates:

- privileged joint-map mechanism: pass;
- privileged source certificate: `2/4` active two-direction bands
  (`4/8` directions);
- matched SMB-NoPriv certificate: `0/4`;
- privileged certificate strictly exceeds SMB-NoPriv: pass.

The held-source centered joint-map alignment was `0.05441` for true pairs,
`-0.01452` after exact-target derangement, and `-0.05263` after within-target
ligand derangement. The registered destruction margins were therefore
`0.06893` and `0.10704`, both above `0.05`. True joint-map cross-entropy
`4.5935` was below uniform `5.5452` and target-destroyed `6.4587`.

The active spectral bands were zero-based bands 2 and 3. Their true held-source
utilities were `0.07309` and `0.15528`, with certificate confidences `0.16513`
and `0.66323`. The no-privileged model had zero active bands under the same
architecture, initialization family, affinity labels, budget, and
certificate.

This is the first DCST route to establish both:

1. exact-target and ligand-specific structural information in the Stage-1
   student; and
2. privileged-specific affinity directions readable through the cross-stage
   interface.

Wall time was `254.396 s`, peak allocated CUDA memory was `943.0 MiB`, and no
downstream affinity label was loaded.

## Authorization

R6 may now generate the frozen 32-segment ChEMBL target feature cache using
sequence-only ESM inference and may train/score the already registered
ChEMBL-37 strict dual-cold train/development Stage 2.

This authorization does not include confirmation or sealed-test scoring.
Those remain untouched unless the complete Stage-2 development gate passes
and a separate decision explicitly authorizes them.

