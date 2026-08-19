# Stage Q2d-1e addendum AD1 — repair of never-executed / nondeterministic branches (2026-08-19)

Frozen BEFORE the Q2d-1e ladder launch.

## Cause

The Q2d-1d ladder (runner_d) crashed at M2 truth generation with
"NameError: name 'PCA_VT' is not defined" in truth_d.py, AFTER completing
all M1 levels A-E x 3 seeds x 8 arms (M1 results in runner_d.log). The
frozen M2 definition ("block-sparse in the pre-compression 510 space,
then compressed") requires PCA_VT, which exists in the frozen label-free
feature artifact q2d1d_features.npz (key "PCA_VT", (32, 510), float64)
but was never loaded. Inspection then showed the NC1/NC2 branches of
truth_d.py were ALSO never executable: they reference A/B feature maps
that no branch assigns (UnboundLocalError), and NC1's frozen description
(I = 0, main effects only) makes the train-std normalization 0/0. These
branches were never executed in any previous stage (Q2d-1c stopped at the
oracle precheck before training). Additionally, the family-preserving
shuffle arm iterated over a raw set() of family ids, which is
PYTHONHASHSEED-dependent: the M1-A recovery run (recover_m1a.py)
reproduced 7 of 8 arms bitwise against the original runner_d.log while
family_preserving differed (e.g. seed 0 dc dz 0.683 vs 0.706), proving
process-level nondeterminism in that arm only.

## Repair (implementation-level; frozen 1d files NOT modified)

Q2d-1e now uses truth_e.py, a byte-copy of truth_d.py with documented,
label-free repairs; runner_e imports truth_e. M1/M2/M3 streams are
asserted bit-identical to truth_d (tests).

1. M2: PCA_VT loaded at module level from the frozen feature artifact.
2. NC1/NC2: A=None, B=None (no feature-conditioned map exists); NC1's
   I_raw_all is zeros per its frozen description; the runner's oracle
   diagnostic arm uses the zero-interaction bound for these mechanisms.
3. NC1 normalization: I = 0 when the train std is 0 (the frozen 0/0 is
   undefined; the description says I = 0).
4. Family-preserving shuffle: iterate sorted(set(...)) of family ids so
   fam_perm is cross-process deterministic. The arm remains an unbiased
   within-family permutation; no gate threshold changes.

No gate, budget, arm set, split, checkpoint rule, or threshold changes.
Addendum SHA frozen in commands.jsonl.
