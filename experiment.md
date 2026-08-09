# Active experiment protocol

## Status

```text
S7_L2B_PHASE1_B5 ................. COMPLETE, DEVELOPMENT PASS 6/6
S7_L2B_PHASE2A_AUDIT ............. COMPLETE, AUDIT-ONLY, NOTHING TRAINED
S7_L2B_PHASE2B_PREREG_R1 ......... COMMITTED b9753db BEFORE IMPLEMENTATION
S7_L2B_PHASE2B_REAL_TRAINING ..... NOT EXECUTED, PRECONDITION FAILED
NEW MODEL TRAINING ............... NOT AUTHORIZED
SOURCE AFFINITY / DAVIS / KIBA / z  FROZEN
GIT PUSH ......................... NOT AUTHORIZED, NOTHING PUSHED
```

## Chronology

| step | commit |
|---|---|
| Phase 2A evidence closure | `0bd1702` |
| Phase 2B preregistration R1, before any implementation | `b9753db` |
| Phase 2B code and contract tests | `0a8b62e` |
| Phase 2B results and status | this commit |

The Phase 2A registration was frozen by hash only; commit `0bd1702` preserves
that evidence bundle and **does not** supply retroactive chronology for it. The
Phase 2B R1 registration **does** carry a commit timestamp that precedes its
implementation.

## Phase 2B as executed

Registered by `research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md`
(`5e6688f6…`), which supersedes `ae6d1a01…`
(`SUPERSEDED_BEFORE_EXECUTION_DESIGN_DEFECT`, never executed, kept
byte-identical; eleven defects in `PHASE2B_DESIGN_AUDIT.md`).

| stage | outcome |
|---|---|
| contract and artifact audit | `PHASE2B_CONTRACT_PASS`, 14/14 preflight items |
| census verification | matched the registration exactly |
| control materialisation | foreign-pair coverage 1.000; derangement 0 fixed points |
| synthetic trainability | **FAILED** `AP_bidir 0.3577 < 0.50` |
| real-label training | **not executed** |
| gates `R1`–`R6` | **not scored** |

```text
TERMINAL VERDICT
    PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED
```

Two defects in the preflight itself were found and fixed before the run rather
than passing silently: the gradient-reachability probe used `d.sum()`, which the
projection annihilates by construction; and the ESM-availability check demanded
states for the sealed confirmation cohort, which correctly has none.

## Phase 2B registration boundary, unchanged

The following remain frozen exactly as registered and must be carried
byte-identical into any repair: the frozen protein-only prior
`b^P = b + alpha*w_pi(GELU(W_h(LN(h))))`; the single trainable head `U` (8×1280)
and `V` (8×41), 10,568 parameters, no bias; `g(L)` = mean of the 41-D atom
features; the float64 Gram–Schmidt gauge with tolerance `1e-8`; the all-residue
gain/loss/change metric aggregated residue → unordered pair → construct →
closure component; the two-ligand foreign-pair and trained-permutation controls;
and gates `R1`–`R6`.

## Decision rule for the next stage

```text
repair the OPTIMIZATION contract only
  -> budget and sampler caps set from a measured synthetic scaling curve
  -> synthetic acceptance threshold derived from that curve, not asserted
  -> fresh synthetic teacher seed (20260905 has been observed)
  -> architecture, projection, metric, controls, gates R1-R6 unchanged

adding capacity, a PLM, attention, a GNN, a geometry branch or a
typed-interaction branch
  -> NOT an admissible response to a synthetic-precondition failure
```

No branch authorizes affinity, DAVIS, KIBA, support adaptation, production `z`,
or any modification of the frozen probability-law operator.

## Historical record

Superseded and failed protocols are recorded in `history.md`,
`report/EXPERIMENTAL_EVIDENCE_LEDGER.md` and
`report/s7_l2b_r0r/PHASE1_ARTIFACT_SUPERSESSION.json`. They are not part of the
active protocol.
