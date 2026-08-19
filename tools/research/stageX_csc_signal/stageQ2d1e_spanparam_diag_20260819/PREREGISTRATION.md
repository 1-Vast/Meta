# Stage Q2d-1e span-param diagnostic preregistration (2026-08-19)

FROZEN BEFORE the Q2d-1e verdict exists (1e is still not launched; this
document is written while Q2d-1d's ladder is running). Per the
re-adjudication this is the ONLY finite diagnostic allowed after Q2d-1e,
and it may not become an open-ended tuning cycle.

## Purpose

Q2d-1e fixes two suspected failure modes of Q2d-1d: (1) the learned protein
map drifts into the train-row null space (33-43% of energy) and (2) A/B
scale ambiguity. If 1e still fails the frozen gate, this diagnostic
distinguishes "null-space redundancy / drift" from "basic optimization
failure": it EXPLICITLY parameterizes the protein map as A = V_train @ G
(G trainable, V_train the frozen feature-only right-singular basis used for
truth generation), so the map can never leave the identifiable span.

## Spec

- Everything identical to Q2d-1e (prereg SHA 61bc0cc5...): same truth
  (stageQ2d1d streams), same features/splits/arms/init policy/minibatch
  order/checkpoint rule/budget (6,000 steps, batch 1024, AdamW lr 5e-3,
  wd 1e-4, correct 8 restarts, negatives 1) and same L2 1e-3 on factor
  maps.
- ONLY change: the model's protein map is A = V_train @ G, implemented as
  fixed V_train^T projection (frozen, requires_grad=False) followed by a
  trainable 28x4 linear map G. The ligand map B and all other parameters
  are unchanged. Span-initialization of Q2d-1e is replaced by span
  PARAMETERIZATION (gradients cannot create null components).
- Ladder: M1 levels A-E x 3 seeds x 8 arms; M2/M3/NC1/NC2 level A;
  value-level reproduction checks; censored_count assertions in D/E.
- Gate: identical (double-cold, median over seeds: correct dz >= 0.70,
  sp >= 0.30, gap vs ligand_only >= 0.05; every negative arm fails;
  correct dz > best negative + 0.03; PASS requires all M1 levels A-E,
  NC1/NC2 must fail).
- Reporting: span energy of the trained map (must be ~1.0 by
  construction), A/B effective rank and scale, train fit, per-surface dz.

## Interpretation (frozen)

- 1e PASS -> learner family passes; this diagnostic never runs.
- 1e FAIL + diagnostic PASS -> the failure was null-space drift/scale
  ambiguity; the span-parameterized variant is the surviving form of the
  low-rank bilinear learner (recorded as such, NOT promoted beyond the
  synthetic gate).
- 1e FAIL + diagnostic FAIL -> basic optimization failure of the low-rank
  bilinear learner under the frozen budget; final verdict FAIL, learner
  closed, no further synthetic successors.
