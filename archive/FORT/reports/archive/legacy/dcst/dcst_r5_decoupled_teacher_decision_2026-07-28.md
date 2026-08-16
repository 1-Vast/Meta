# DCST-R5 decoupled privileged-teacher decision

Date: 2026-07-28  
Decision: `STOP_R5_ADVANCE_TO_SHARED_MECHANISM_BOTTLENECK`

## Admissible result

The corrected execution is
`dcst_r5c_stage1_seed1729.json`. The earlier artifact with the unsuffixed R5
name is invalid under
`dcst_r5_structure_only_decay_correction_2026-07-28.md`.

The correction worked at the engineering level. The teacher-frozen candidate
learned eight nonzero singular values (`0.312` to `0.653`) instead of the
invalid run's exact zeros. It therefore tested the registered hypothesis.

The hypothesis failed all frozen source-only gates:

- privileged teacher joint-mechanism gate: fail;
- privileged source certificate: `0/4` active bands;
- no-privileged source certificate: `0/4`;
- random-frozen source certificate: `0/4`;
- privileged certificate strictly exceeding both controls: fail.

The teacher's joint-map cross-entropy contained absolute target information:
`4.426` for the true pair, `6.300` after target destruction, and `4.451`
after ligand destruction, against uniform `5.545`. It did not carry
within-target ligand-reordering information. Centered alignments were
`-0.0102` true, `-0.0204` target-destroyed, and `0.0105` ligand-destroyed.
The corresponding registered margins were `0.0101` and `-0.0208`, both below
`0.05`.

All four candidate spectral-band utilities were non-positive on true held
source episodes (`-0.014`, `-0.029`, `-0.279`, `-0.131`) or failed a
destruction comparison. Wall time was `316.793 s`; peak allocated CUDA memory
was `943.4 MiB`.

## Mechanistic diagnosis

R5 exposed a broken information interface rather than insufficient parameter
count. Structural supervision trains the segment-by-Morgan joint-interaction
head, while the affinity matrix reads separate pooled target and ligand
vectors. Freezing the teacher preserves the former but does not make it
linearly readable from the latter. Joint training partially helped the
structural probe in R4 but then overwrote every affinity certificate; freezing
in R5 removed that conflict by disconnecting the useful signal entirely.

The next route must make the Stage-1 structural interaction map itself the
only affinity bottleneck. A generic residual adapter is not yet justified:
it could improve fit without proving that high-quality Stage-1 information
was transferred.

No ChEMBL affinity label was loaded or scored in R5 or its correction.

