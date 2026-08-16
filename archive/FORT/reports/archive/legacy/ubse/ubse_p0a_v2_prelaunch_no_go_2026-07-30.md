# UBSE-P0A-v2 prelaunch multi-review

Date: 2026-07-30  
Decision: `NO_GO_UBSE_P0A_V2_LONG_TRAINING_PRELAUNCH_AUDIT`

## Scope

This review asked whether the A1-C-closed P0A-v2 implementation was ready
for a roughly 5.5-hour, three-seed CUDA run. Three independent agents reviewed
the scientific closure, implementation, execution recovery, and GPU budget.
No contact-label body, affinity, event, confirmation, coordinate body, or
sealed outcome was read.

P0A-v2 remains a target-marginal proposal component. It is not pair-specific
evidence and is not required by the primary no-P0A A1-direct gate.

## Cheap checks that passed

- Safe-column closure reproduced 54,868 rows and 32,769 exact targets.
- The retained sequences contain 11,126,109 residues and 33,216 windows.
- The retained maximum A1-C ECFP4 Tanimoto is 0.495238, below 0.5.
- The A1-S retained development identity contains 59 panels.
- Focused CPU tests passed: 13/13.
- Python compilation and `git diff --check` passed.
- Fixed-shuffle calls now use the model seed at epochs 1, 2, and 4.
- The v2 output names do not collide with the accepted historical P0A
  weights, ledger, or result.
- The `drug` environment exposes PyTorch 2.6.0+cu124 and an RTX 4060 Laptop
  GPU.

These checks establish a plausible substrate, not execution readiness.

## Blocking findings

1. The purported preregistration SHA is only a nonempty command-line string;
   no real preregistration file is hashed before label access.
2. Training-label rows are checked only by count. The label projection omits
   PDB and connectivity identities, so an equal-size wrong membership could
   pass.
3. Validation labels are projected by PubMed union instead of an exact
   59-panel `(target, PubMed, scaffold, PDB, connectivity)` membership.
4. The 59-panel development role has not been independently closed against
   the whole 54,868-row training manifest on target/homology, PDB, PubMed,
   connectivity, scaffold, and ECFP4.
5. Source PDB values are upper case while A1-C PDB values are lower case.
   The current exact mask reports zero PDB exclusions. A case-folded audit
   finds 655 matching source rows; all happen also to be removed by PubMed,
   so the final row count is unchanged but the PDB-axis certificate is false.
6. Target-sequence uniqueness is checked only after `drop_duplicates`, so
   the check cannot detect one target mapped to multiple sequences.
7. A1-C homology membership is returned as a hard-coded zero instead of
   being recomputed independently on the retained manifest.
8. A completed seed writes a weight but not its validation ledger or result.
   A later crash leaves an unusable file that blocks a clean rerun. There is
   no verified per-seed or per-epoch resume contract.
9. Final and temporary run-state paths are not all checked before CUDA work,
   and a non-finite seed does not immediately prevent later seeds.
10. The ledger check counts rows but does not prove the full Cartesian grid
    over seed, evaluation epoch, retained panel, and control.
11. Telemetry can silently return zero samples; there is no persistent
    heartbeat or progress watchdog.
12. The accepted old run peaked at 7,945/8,188 MiB. Dataset shrinkage does
    not reduce the worst batch shape, and no representative AMP
    forward/backward peak-memory smoke has passed for the current free-memory
    state.

## Resource decision

The old three-seed run took 6.583 hours over 39,440 windows per seed.
Scaling by `33,216 / 39,440 = 0.8422` gives 5.50-5.54 hours for v2.
Launching before the findings above are closed would risk losing most of that
budget without a valid scientific package.

Estimated work before any future launch:

- implementation and focused tests: 1-1.5 hours;
- label-blind preflight and worst-batch CUDA smoke: 10-15 minutes;
- three-seed execution: 5.50-5.55 hours;
- acceptance aggregation: 15-30 minutes.

## Reopening conditions

P0A-v2 training may be reconsidered only after all of the following:

- real preregistration and execution-amendment files are hash-bound before
  label projection;
- all PDB identifiers are normalized and every six-axis A1-C/A1-S closure is
  independently recomputed;
- label-bearing rows are exactly equal to frozen training and validation
  manifests, not merely equal in count;
- sequence and homology checks are independent and nonconstant;
- atomic per-seed packages and verified resume behavior preserve weights,
  ledgers, RNG/protocol state, and telemetry;
- the exact ledger Cartesian product and frozen shuffle digests are tested;
- persistent heartbeat and telemetry failure gates exist;
- a representative maximum-shape CUDA AMP forward/backward smoke passes
  below the preregistered memory ceiling.

Even after those conditions pass, the 5.5-hour run should be scheduled only
if A1-direct demonstrates a nonzero, powered coupling object or another
authorized A1 branch explicitly requires the proposal. The primary no-P0A
arm remains mandatory.

