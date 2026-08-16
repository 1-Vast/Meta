# DCST-R7 content-addressed role transport decision

Date: 2026-07-28  
Decision: `STOP_R7_REPLACE_LEARNED_ROUTER_WITH_FROZEN_ATLAS`

## Result

`dcst_r7_stage1_seed1729.json` failed every frozen source gate:

- segment-level privileged mechanism: fail;
- CART privileged certificate: `0/4`;
- CART-NoPriv certificate: `2/4`;
- privileged certificate strictly exceeding CART-NoPriv: fail.

The true centered structural alignment was `0.03125`, versus `0.00531`
target-destroyed and `-0.02444` ligand-destroyed. The ligand margin passed
(`0.05568`), but the target margin was only `0.02594`. The privileged
role-energy matrix was nearly rank one: singular values began `0.5970`,
`0.1163`, `0.0340`, and `0.0091`. No privileged band survived both
destruction controls.

No downstream affinity label was loaded. Wall time was `227.032 s`; peak
allocated CUDA memory was `945.1 MiB`.

## Router diagnosis

The learned soft router did not collapse to one role. On 87 held-source exact
targets:

- privileged mean normalized per-segment routing entropy: `0.9847`;
- privileged effective aggregate roles: `7.85/8`;
- no-privileged routing entropy: `0.9999`;
- no-privileged effective roles: `8.00/8`.

The failure is therefore role non-identifiability: nearly uniform soft roles
can rotate or exchange while the affinity matrix compensates. They are
content-conditioned in code but do not define anchored cross-protein
coordinates. The no-privileged model's 2/4 bands further proves that this
flexibility creates affinity directions without making privileged structure
load-bearing.

R7 is stopped. The next route must make role identity deterministic and
immutable across stages, using no affinity or structural label to define the
coordinate system.

