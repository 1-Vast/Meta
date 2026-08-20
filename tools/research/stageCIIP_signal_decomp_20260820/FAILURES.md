# FAILURES / defects log (CIIP-S1)

## F1 (fixed before any verdict; seed-1 rerun)

- Defect: Form-2 (and Form-1 cell lists) included non-finite target cells when
  the estimand is T2 (T2 = c - cross-fitted profile contains NaN at ligand
  columns not covered by any same-parent train sibling). A single NaN cell
  made the batch MSE NaN, poisoning all model weights -> all T2 arm cells
  evaluated as n=0 in the first seed-1 run.
- Detection: SEED1_RESULT.json T2 aggregates n=0; direct probe showed
  predictions all-NaN after 3 epochs.
- Fix: cells filtered to np.isfinite(target) in both trainers (metric-side
  finite filtering was already in the frozen contract B.7).
- Consequence: the first seed-1 run was DISCARDED entirely; seed 1 re-run
  from the fixed code before any adjudication. No S1 result was adjudicated
  from the broken run. T1/T0/T3 cells were unaffected numerically by the bug
  (their targets are fully finite), but the whole run was rerun anyway so
  that all cells come from one frozen code state.
- No threshold or rule was changed by this fix.

## F2 (job-level; no scientific impact)

- The first background seed-1 job stalled at ~0 CPU (torch default thread
  oversubscription in the sandboxed job). Killed; relaunched with
  torch.set_num_threads(8) + OMP/MKL_NUM_THREADS=8 and per-arm progress
  prints. The stall produced no results and touched no outputs.

## F3-F6 (independent code review 2026-08-20; fixed per ADDENDUM_ADD2)

- F3: T3 rank-loss indexing bug (ordv indexed the full prediction vector
  instead of the selected cells -> wrong ligands compared). Voided all T3
  cells of the diagnostic run; corrected via rank_loss helper (unit test
  asserts corrected semantics differ from the buggy one).
- F4: T0m arms missing (driver trained only T0 and derived severity;
  prereg B.5 requires dedicated per-pair-scalar T0m arms). Corrected; A
  contrast now F1f-T0m vs F2-T0m.
- F5: A-contrast used Pearson contributions; prereg freezes cross-pair
  Spearman. Corrected (Spearman difference, paired parent-cluster bootstrap).
- F6: Form-2 T3 arms trained with MSE on c (T1 duplicates) instead of the
  rank loss. Corrected to the same logistic rank loss as Form-1.
- F7: adjudicator NULL-ALL screen covered only T1. Corrected to every
  estimand vs its floor (T0m: Spearman vs F7f-T0m floor).
- The diagnostic run (pre-correction) is archived as
  SEED1_RESULT_DIAG_T3BUG.json and is NEVER adjudicated.

## F8 (self-review after ADD-2 run; fixed per ADDENDUM_ADD3)

- The Form-2 T3 rank loss introduced by the ADD-2 compliance fix had a
  FLIPPED margin (h[a]-h[b] instead of h[b]-h[a]): every Form-2 T3 arm was
  trained toward anti-ranked output (negative median Spearman across all
  Form-2 T3 cells while Form-1 T3 was positive). Fixed; unit test
  test_form2_t3_rank_sign added (17/17 PASS). The ADD-2 seed-1 run is
  archived as SEED1_RESULT_ADD2_T3SIGNBUG.json (never adjudicated); all
  seeds re-run from the final code state. T3 is not a gating estimand for
  any frozen contrast.

## All other checks

- Erasure asserts, structure tests (16/16 after compliance tests), smoke
  gate (re-verified on corrected code), S0 audits: no failures.
- "Mean of empty slice" RuntimeWarnings from F9 profile construction are
  expected (ligand columns uncovered by any same-parent train sibling) and
  are handled by finite-mask filtering in training and metrics; not a crash,
  not a defect.
