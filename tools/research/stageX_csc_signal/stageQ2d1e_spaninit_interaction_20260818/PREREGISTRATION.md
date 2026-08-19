# Stage Q2d-1e preregistration - span-initialized factors + mild L2 (2026-08-19)

Corrected successor of Q2d-1d. Frozen BEFORE the Q2d-1d ladder finishes and
before its gate verdict exists, so the fix is not tuned to the verdict.
SHA-256 in PREREGISTRATION_SHA256.txt. Synthetic only; old stages remain
read-only. Baseline unchanged: pipeline qualification FAILED at Q2 (X0c);
biological conclusion UNRESOLVED; B1/B2/C/D NOT AUTHORIZED.

## Change vs Q2d-1d (with evidence)

Q2d-1d attribution diagnostic (frozen in _attr_correct.py, run beside the
ladder): the trained correct arm fits the training set essentially
perfectly (train dz 0.989-0.993, sp 0.977-0.978 across restarts) while
carrying 33-43% of the learned protein map energy in the train-row null
space (the 4 unidentifiable directions of the rank-28 train-row feature
span). Train loss is insensitive to those directions, so the frozen
best-train-loss checkpoint rule cannot select against the drift, and cold
surfaces degrade to dz ~0.45-0.72 depending on the init lottery, while the
closed-form oracle (min-norm least squares) reaches 0.753-0.893.

Fix (two parts, both frozen):
1. SPAN-INITIALIZED PROTEIN MAP. The learner s protein map A is
   initialized as V_train @ G with V_train the frozen feature-only
   right-singular basis of the train-row feature submatrix (identical to
   the basis used for truth generation; unsupervised, label-free, already
   frozen in Q2d-1d s feature artifact) and G xavier-uniform. The truth
   map lies entirely in that span, so a learner that never develops null
   components can represent it; gradients in the null directions are
   (nearly) zero, so the trained map stays inside the span and cold
   predictions use only identifiable structure. No truth or label
   information enters the initialization.
2. MILD L2 ON FACTOR MAPS. The training loss adds
   lambda * (||A||_F^2 + ||B||_F^2) with lambda = 1e-3 (frozen) in every
   arm and level (the same term enters the monitor loss used for
   checkpointing). This counters residual null-space drift and the
   A/B scale ambiguity. Applied identically to all arms.

Everything else is carried over unchanged from Q2d-1d (prereg SHA
baf4bb72...): features (PCA-32 protein, ECFP-projected 48-dim ligand, 157
resolved ligands), feature-smoothed double centring, mechanisms
M1/M2/M3/NC1/NC2 with span-restricted truth maps, splits and surfaces,
all 8 arms with identical rows/init policy/minibatch order/checkpoint
rule (correct 8 restarts, negatives 1), budget (6000 steps, batch 1024,
AdamW lr 5e-3, wd 1e-4), ladder A-E with value-level reproduction,
frozen censored_count assertion in D/E, the never-moved gate (double-cold,
median over 3 seeds: correct dz >= 0.70, sp >= 0.30, gap vs ligand_only
>= 0.05; every negative arm fails; correct dz > best negative + 0.03),
and the oracle precheck before training (already PASSED on this exact
truth: pc 0.700-0.920, lc 1.0, dc 0.753-0.893 for all 3 seeds).

## Adjudication and reporting (unchanged)

The gate is evaluated per ladder level; PASS requires the gate to hold at
every level A-E of M1 (NC1/NC2 must fail; M2/M3 reported). The oracle is
diagnostic only and never influences training, checkpoints, or selection.

## Downstream authorization (unchanged)

Q2d-2 (representation matrix incl. KLIFS per-position ESM kept as
[position, embedding]) only after Q2d-1e PASS; Q2d-3 after Q2d-2 PASS;
Saifudeen B1 requested ONLY after all synthetic gates PASS;
functional-inhibition endpoint never called pK/Ki/Kd/DTA.

## Governance (unchanged)

SHA-256 seeds; no Python hash(); artifacts carry schema / prereg SHA /
input SHA; commands.jsonl appended; restricted data never committed;
research code only under tools/research/stageX_csc_signal/; model/ and
production scripts/ untouched before B1; implementation and report commits
separated.
