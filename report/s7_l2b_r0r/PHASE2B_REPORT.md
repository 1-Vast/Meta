# S7/L2B Phase 2B — contract repair, preflight, and a fail-closed stop

Date: 2026-08-10.
Preregistration `research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md`,
SHA-256 `5e6688f68bf214b6f44c96ef0a2b909eba99da31f1adf17cdde003addb242c96`,
**committed as `b9753db` before any Phase 2B implementation file was written**.

```text
TERMINAL VERDICT

    PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED

REAL-LABEL TRAINING   NOT EXECUTED
GATES R1-R6           NOT SCORED
BIOLOGICAL CONCLUSION NOT PERMITTED
```

The registered synthetic trainability control returned `AP_bidir = 0.3577`
against a preregistered requirement of `>= 0.50`. Under the frozen protocol that
stops the stage before the real-label run. The threshold was not lowered, no
hyperparameter was tuned against the synthetic holdout, and no second seed was
tried.

No affinity, DAVIS, KIBA or recipient value was read. The sealed confirmation
cohort was not opened. `model/`, `scripts/`, production `z`, CSMO, Band, mesh and
`A(F,z) = K(B(z)F(z))` are untouched.

---

## 1. The superseded preregistration

`PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL.md` (`ae6d1a01…`) is retained
**byte-identical** and marked `SUPERSEDED_BEFORE_EXECUTION_DESIGN_DEFECT` in
`PHASE1_ARTIFACT_SUPERSESSION.json`. It was never executed, so nothing is
withdrawn — a design was replaced before it was used.

`PHASE2B_DESIGN_AUDIT.md` (`22c132b6…`) records eleven defects with the exact
clause at fault. The four that would have silently corrupted the result:

- **D1** — `b_r(P)` was defined as the residue term of a **per-complex** additive
  decomposition of B5 pair logits. That object is fitted per `(P, L)`, so it is
  ligand-dependent; the same-protein cancellation the whole differential design
  rests on would not have held. It also existed only on held-out A, never on the
  training split.
- **D3** — the projection span `{1, b(P), c(L)·1}` contains two collinear
  columns, so the stated projector is singular.
- **D5** — the primary metric ranked only the **symmetric-difference** residues,
  i.e. the candidate set was chosen using the answer.
- **D7** — the module-participation audit required that detaching the frozen
  `h_r` change the result. `h_r` has no upstream parameters, so that is a
  mathematical no-op and could never be satisfied.

A separate integrity finding: `P1_B5_REPORT.md` no longer matches the hash
recorded for it in the Phase 1 triage (`19c9c205…` → `dbfe8b92…`). The change is
wording only, in section 3(b), and no number moved. It is recorded rather than
reverted or ignored.

## 2. What was verified before anything ran

`PHASE2B_INPUT_MANIFEST.json` — `PHASE2B_CONTRACT_PASS`, all 14 registered
preflight items.

| check | result |
|---|---|
| trainable parameters | **10,568**, names exactly `{U, V}`, no bias |
| `g(L)` atom-permutation invariance | max abs diff **0.0** |
| ligand-order swap flips `Δs` | max abs sign error **0.0** |
| `b^P` cancels in the same-protein difference | **2.05e-15** (tolerance 1e-12) |
| projection orthogonality `‖QᵀΔ‖/‖Δ‖` | **6.19e-15** (tolerance 1e-8) |
| degenerate-`b` Gram–Schmidt fallback | rank 1 constant / rank 2 generic |
| train ↔ held-out A component overlap | **0** |
| held-out A ligand-graph overlap with train | **0** |
| held-out B scaffold overlap with train | **0** |
| affinity-marked paths opened | **none**; affinity reads **0** |

Two defects in my **own preflight** were found and fixed before the run, and are
recorded because a silent pass would have been worse than a failure:

1. the gradient-reachability probe used `d.sum()` as its objective. The constant
   direction is the first column of `Q`, so the projection annihilates exactly
   that functional and the gradient was `~7e-12` — a "pass" that proved nothing.
   Replaced with a fixed generic random linear form; gradients now reach `U` and
   `V` at order `1e-1`.
2. the ESM-availability check demanded states for every construct in the corpus
   and failed on 573. Scoped to the constructs Phase 2B actually requires: **0
   missing**. Of the 573, 445 belong to the **sealed additional-PDB confirmation
   cohort** — which correctly has no states — and 131 are development constructs
   whose records were all removed by the ligand-graph disjointness filter and so
   appear in no Phase 2B split. **0** of them appear in the split.

## 3. Census — verified, not assumed

| quantity | measured | registered expectation |
|---|---:|---:|
| train constructs with pairs | **760** | ~766 |
| train closure components with pairs | **554** | ~554 |
| train eligible pairs | **226,765** | ~226,765 |
| held-out A constructs with pairs | **174** | ~175 |
| held-out A components with pairs | **112** | ~112 |
| held-out A eligible pairs | **46,818** | ~46,818 |
| held-out B eligible pairs | **30,661** | — |

Exclusions (train / held-out A / held-out B): zero symmetric difference
676 / 198 / 132; same ligand graph 1,154 / 254 / 190; scaffold not distinct
12,828 / 2,016 / 1,041.

## 4. Controls — materialised and hashed before training

`PHASE2B_CONTROL_MANIFEST.json` (`d5943258…`).

- **`R3` two-ligand foreign pair**: both ligands replaced, all four graph keys
  and all four scaffolds distinct, matched by nearest heavy-atom count with
  pooled-feature distance as tiebreak, drawn from a 7,546-ligand training pool,
  no label or score consulted. **Coverage 46,818 / 46,818 = 1.000.**
- **`R5` within-construct derangement**: 11,123 records mapped, **0 fixed
  points**, 1,043 singleton constructs correctly not permutable and counted.

## 5. The synthetic trainability control, and why it stopped the stage

The teacher is a rank-8 projected bilinear differential drawn with seed
`20260905` — it lies **exactly** in the candidate's own hypothesis class by
construction. Per pair, the synthetic gain set is the top 8 residues by
`Δδ*` and the loss set the bottom 8.

```text
required   AP_bidir >= 0.50
observed   AP_bidir  = 0.3577      chance 0.0376      FAIL
```

Four gauge-invariant diagnostics localize the shortfall. None of them is a gate
and none can rescue the result.

| diagnostic | value | reading |
|---|---:|---|
| teacher scored on its own labels | **0.99971** | the metric and evaluation code are sound |
| in-sample AP on the final-epoch sampled pairs | **0.3654** | — |
| held-out AP | **0.3577** | **no generalisation gap**; this is underfitting |
| output-level Pearson `r`(learned δ field, teacher δ field) | **0.717 mean / 0.754 median** | the class is being fitted; the ranking is not sharp enough |
| parameter movement | `U` **1.808**, `V` **0.426** | the head trained; it did not sit still |

So the failure is **not** the hypothesis class, **not** the evaluation code, and
**not** generalisation. It is the registered **optimization budget**: 6 epochs
over at most 8,864 sampled pairs (554 components × 2 constructs × 8 pairs) drive
the learned field to `r ≈ 0.75` of the teacher, but not into an exact top-8
ranking at `AP_bidir ≥ 0.50`.

There is a second, honest possibility that must be stated rather than argued
away: **the threshold itself may have been mis-calibrated.** I set `0.50` a
priori with no scaling curve to calibrate against. A field correlation of 0.75
recovered from a random rank-8 teacher is substantial. Distinguishing "budget too
small" from "threshold too strict" is exactly what the next registration must
settle **before** it touches real labels — and it cannot be settled by adjusting
either number now, because both were frozen and the synthetic holdout has been
seen.

## 6. What is therefore not concluded

Nothing about biology. In particular this run does **not** show that the frozen
sequence + 2-D ligand representation lacks a ligand-conditioned residue
correction. `R1`–`R6` were never scored. The Phase 2A finding stands unchanged
and unaffected: the MONN labels are ligand-conditioned at the residue level
(`ΔJ = +0.258 [LCB +0.234]`, chemistry association `ρ = +0.322`), and B5's
residue marginal is not.

## 7. The sole authorized next action

Preregister **one** repair of the Phase 2B **optimization contract** — not the
biology, not the architecture, not the gates on real labels. That registration
must, before any real label is read:

1. fix the training budget and sampler caps on a **measured synthetic scaling
   curve** (AP as a function of epochs × pairs), computed on synthetic data only;
2. derive its synthetic acceptance threshold **from that curve** rather than by
   assertion, and state the derivation;
3. use a **fresh synthetic teacher seed**, because seed `20260905` has now been
   observed;
4. leave the architecture, the projection, the metric, the controls and gates
   `R1`–`R6` byte-identical to `PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md`.

This is not authorized by this run and has not been performed. Adding capacity,
another PLM, attention, a GNN, a geometry branch or a typed-interaction branch is
**not** an admissible response to this failure and is not proposed.

## 8. Frozen boundaries still in force

Real ChEMBL/BindingDB affinity training; DAVIS, KIBA and recipient labels; the
sealed independent structural confirmation; any additional PLM, GNN, attention
stack, geometry branch, typed-interaction branch, affinity head, PU loss,
knowledge graph or parallel SSL module; ESM2 fine-tuning; hyperparameter search;
repeated seeds for selection; few-shot section adaptation; admission of any
statistic into production `z`; CSMO, Band, mesh, positive ridge and
`A(F,z) = K(B(z)F(z))`; P2–P4; pushing to GitHub.

Held-out A has been read during Phase 1, Phase 2A and Phase 2B design. It is
**development evidence, not independent confirmation.**

## 9. Artifacts

| file | SHA-256 (first 16) |
|---|---|
| `PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md` | `5e6688f68bf214b6` |
| `PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL.md` (superseded) | `ae6d1a0186bb37af` |
| `PHASE2B_DESIGN_AUDIT.md` | `22c132b66d1ee5d5` |
| `p2b_residue_residual.py` | `54483ed62ac39ced` |
| `p2b_run.py` | `483f2ccfef91f24f` |
| `tests/test_s7_l2b_phase2b.py` | `60dee0e7768d6467` |
| `PHASE2B_INPUT_MANIFEST.json` | `e78cb60489d93740` |
| `PHASE2B_CONTROL_MANIFEST.json` | `d5943258471e0377` |
| `PHASE2B_SYNTHETIC_AUDIT.json` | `ffe75ef15cea0863` |
| `PHASE2B_TRAINING_TRACE.json` | `bec0186cd80f4771` |
| `PHASE2B_GATE.json` | `a3939066891c9749` |

Device: CPU, chosen **before any result was seen** so the registered bit-exact
determinism check is achievable rather than a gamble on cuBLAS reduction order;
the head is 10,568 parameters and every heavy tensor is frozen. Python 3.11.15,
torch 2.6.0+cu124, numpy 1.26.4, RDKit 2023.09.6, CUDA 12.4 available but unused.
Seeds: parameters `20260901`, sampler `20260902`, bootstrap `20260903`, controls
`20260904`, synthetic `20260905`.

Regression: **100 passed** (75 pre-existing, verified, plus 25 new Phase 2B
contract tests).
