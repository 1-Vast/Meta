# Q2d terminal summary — synthetic interaction qualification chain (2026-08-19)

Authoritative record of the bounded Q2d chain. Leaf evidence: each
stage's frozen artifacts, runner logs, and the frozen adjudicators.
This document moves no gate and modifies no historical artifact.

## Current chain state

| Stage | Verdict | Evidence |
|---|---|---|
| Q2d-1b feature interaction | superseded by 1c | prereg 25b8b912... lineage |
| Q2d-1c identifiable interaction | STOP at oracle precheck (protein-cold not recoverable in closed form on the then-truth) | stageQ2d1c report |
| Q2d-1d span-restricted | **GATE FAIL** (M1:A: correct dc dz 0.5616, sp 0.128, best negative 0.612) | Q2D1D_GATE.json |
| Q2d-1e span-init + L2 1e-3 | **GATE FAIL** (all M1 levels A-E: correct dc dz 0.585/0.588/0.589/0.410/0.489, sp 0.191/0.150/0.195/-0.234/0.065; best negative beats correct at C/D/E) | Q2D1E_GATE.json |
| span-param diagnostic A=V_train.G | **GATE FAIL** (M1 A-E: correct dc dz 0.669/0.544/0.549/0.630/0.508; sp up to 0.399; best negatives still competitive; repro all True; cens 165/165; NC1/NC2 fail as required) | Q2D1E_SPANPARAM_GATE.json |
| Terminal verdict on the low-rank bilinear learner family | **CLOSED (FAIL)** — basic optimization failure under the frozen budget | this document |

## Adjudication verification checklist (executed, all confirmed)

- Frozen adjudicators used unmodified: adjudicate_d.py (1d, M1:A) and
  adjudicate_e.py (1e, M1 A-E + NC1/NC2 must fail). Adjudicator exit
  codes recorded (1 = FAIL in both).
- Preregistration SHAs verified against the files and the ladder JSONs:
  1d baf4bb72..., 1e 61bc0cc5..., diag 772d6d46... (PREREGISTRATION_
  SHA256.txt files match).
- Launch commands recorded in commands.jsonl (runner_d 13:49Z PID
  20348; runner_e 07:13Z PID 54828; runner_diag 17:12Z PID 55308).
- Protocol identical across 1d/1e/diag: seeds [0,1,2]; 8 arms (correct,
  ligand_only, additive_only, shuffled_protein,
  family_preserving_shuffle, random_protein, no_interaction_head,
  oracle_diagnostic); correct = 8 restarts, negatives = 1 restart;
  6000 steps, batch 1024, AdamW lr 5e-3, wd 1e-4; checkpoint = best
  monitor loss (1e adds the same L2 term to the monitor).
- Value-level reproduction: 1e ladder repro B/C/D/E all True (bitwise
  re-runs of level-A config); censored-count assertions D/E = 165 pass.
- Negative controls: every negative arm fails the gate in every M1
  level (adjudicator every_negative_fails = True), NC1/NC2 fail as
  required. Note: NC1 dz/sp are NaN because the NC1 truth interaction
  is exactly zero (predictions constant); the frozen adjudicator treats
  NaN as failing, which is the required outcome; the NaN appears
  verbatim in Q2D1E_GATE.json.

## What each successor fixed, and what it did not

- **1c vs 1b**: PCA-32 protein features, resolved-only ligand pool
  (157), feature-smoothed double centring. Fixed the measurement
  setup; the closed-form oracle was still unrecoverable on
  protein-cold -> STOP before training (no gate movement).
- **1d vs 1c**: truth protein map drawn entirely inside the train-row
  feature span (A_t = V_train @ C). The oracle precheck then PASSED
  (pc 0.700-0.920, lc 1.0, dc 0.753-0.893) — the truth is recoverable
  in principle. The TRAINED learner failed: 33-43% of the learned
  protein-map energy drifted into the unidentifiable null space, and
  double-cold dz fell to 0.50-0.57 while a family-preserving shuffle
  beat the correct arm. Fix target: parameter identifiability.
- **1e vs 1d**: span-initialized A (V_train @ G at init) + L2 1e-3 on
  the factor maps. The fix WORKED on its own metric: trained null
  energy dropped to 0.07-0.13 (preflight, committed). The gate still
  FAILED at every M1 level (correct dc dz 0.410-0.589, sp <= 0.195,
  best negatives competitive or better). Conclusion: null-space drift
  was real but NOT the only (or main) cold-failure cause.
- **Diagnostic vs 1e**: A = V_train @ G as a PARAMETERIZATION (gradient
  cannot create null components; null energy 0 by construction). Its
  only purpose is to distinguish "null-space redundancy/drift" from
  "basic optimization failure" under the frozen budget. No further
  successor exists or may be created after it.

## Failure-mode classification (Q2d-1e, recorded; diagnostic pending)

1. Implementation bug: excluded as the primary cause. Reproduction
   checks bitwise-pass; censored assertions pass; the AD1 truth-module
   repairs (PCA_VT load, NC1/NC2 executability, deterministic family
   order) are documented and leave M1/M2/M3 streams bit-identical.
2. Measurement-definition failure: the gate metrics are defined and
   the oracle reaches 0.76-1.0 on double-cold, so the gate is not
   vacuous.
3. Optimization failure: prime candidate — the 6000-step AdamW run
   fits the train grid (train dz ~0.99) but does not find the
   closed-form-achievable cold solution; at SNR = 1 (noise sd equals
   signal sd) the trained maps stay close to the init lottery.
   DECIDED BY THE DIAGNOSTIC.
4. Cold-generalization failure: the in-span, train-perfect maps still
   degrade on cold surfaces (level D correct 0.410 below ligand_only
   0.489) — consistent with overfitting train-only interaction
   structure at SNR 1.
5. Biological signal: synthetic only; no biological claim is made or
   falsified by any Q2d stage. Core Task 1 remains UNRESOLVED.

## Downstream authorization (frozen)

- Diagnostic PASS -> recorded conclusion limited to: the original
  learner failed mainly from unidentifiable parameterization /
  null-space drift; the span-parameterized variant is retained as a
  synthetic candidate; NOT promoted beyond the synthetic gate.
- Diagnostic FAIL -> the low-rank bilinear synthetic learner family is
  CLOSED; no further synthetic successors; Core Task 1 UNRESOLVED (not
  a biological falsification).
- Q2d-2 (representation matrix), Q2d-3, and Saifudeen B1 remain
  unauthorized (they required Q2d-1e PASS).
- The P-line (practical few-shot performance) is not blocked by this
  chain: arms 1-7 implemented and structure-tested; real training
  proceeds under the frozen promotion budget (arm 3 incumbent first).
- Core Task 1 validation now follows the CIIP funnel (short four-arm
  screen; real matched WT/variant or ligand-identical protein-pair
  data; KiRHub census -> DATA BLOCKER; Davis census -> INSUFFICIENT
  ALONE, combined-panel census next).

## Terminal verdict (2026-08-19, 19:20Z)

**The low-rank bilinear synthetic learner family is CLOSED (final
verdict FAIL).** The frozen span-param diagnostic (A = V_train @ G)
ran to completion — censored assertions passed (D/E = 165), all four
value-level reproduction checks passed bitwise, NC1/NC2 failed as
required — and the frozen adjudicator returned GATE FAIL at every
ladder level: correct-arm double-cold dz = 0.669 (A), 0.544 (B),
0.549 (C), 0.630 (D), 0.508 (E), M2 0.620, M3 0.516, all below the
never-moved 0.70 gate, with best negatives competitive (A 0.662).

Interpretation (frozen): the span parameterization improved the
cold surfaces (A: dz 0.585 -> 0.669, sp 0.191 -> 0.344; D: dz 0.410 ->
0.630 vs 1e), confirming that null-space drift was a real contributor,
but the learner still fails the gate under the frozen budget — the
residual failure is BASIC OPTIMIZATION/ESTIMATION FAILURE of the
low-rank bilinear learner at SNR 1 (train-fit perfect, cold surfaces
do not reach the closed-form-achievable oracle). No further synthetic
successors will be created; Q2d-2/Q2d-3/Saifudeen-B1 stay
unauthorized; Core Task 1 remains UNRESOLVED (this is not a
biological falsification). The Q2d chain is TERMINALLY ARCHIVED.
