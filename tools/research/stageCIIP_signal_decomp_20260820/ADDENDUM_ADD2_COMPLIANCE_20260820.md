# S1 ADDENDUM ADD-2 — protocol-compliance corrections (frozen 2026-08-20)

Frozen BEFORE the corrected seed-1 run. Triggered by an independent read-only
code review (2026-08-20) that found the frozen implementation deviated from
PREREGISTRATION.md B.4/B.5/B.8 in five places. These are compliance fixes to
make the code match the ALREADY-FROZEN protocol; no threshold, rule, arm, or
estimand definition is changed.

## Deviations found and corrected

1. T3 rank-loss indexing (train_form1). ordv sorted the batch-selected cells
   but margins indexed the full per-pair prediction vector h directly, so the
   loss compared the WRONG ligands (first-len(sel) ligands instead of the
   selected ones). Corrected: margins = h[sel][ordv[1:]] - h[sel][ordv[:-1]]
   (helper rank_loss; unit-tested against the buggy semantics).
   Affected cells: F1-T3 (Form-1) — all previous T3 results void.
2. Form-2 T3 arms trained with MSE on c (i.e. were T1 duplicates) instead of
   the preregistered pairwise logistic rank loss. Corrected: same
   adjacent-in-target-sorted-order logistic margins on batch cells grouped by
   pair. Affected cells: all Form-2 T3 arms — previous T3 cells void.
3. T0m was never trained as its own arm: the prereg B.5 requires per-arm
   T0m targets (per-pair scalar mean_l d, B.4); the driver instead derived a
   severity summary from T0-trained arms. Corrected: dedicated T0m arms
   (F1f/F2/F2w/F3/F4/F7f/F8f), one training example per pair, input =
   (protein features, panel-mean ECFP4 of the pair's ligand set) per the
   frozen "(protein-features, ECFP4)" input contract; F7f-T0m zeroes the
   protein block; F8f-T0m drops the ligand block. Checkpoint = best val MSE
   on the scalar target. The A-proposition primary contrast now uses
   F1f-T0m vs F2-T0m; the T0-derived severity contrast is retained and
   reported as SECONDARY (A_secondary_T0derived), never gating.
4. The A-contrast used Pearson contributions; prereg B.8 freezes cross-pair
   SPEARMAN. Corrected: severity_contrast computes Spearman_A - Spearman_B
   over the 9 test pairs with a paired parent-cluster bootstrap of the
   difference (same parent resample for both arms, 2000 keyed draws) and
   LOPO sign stability (sign of the difference excluding each parent).
5. NULL-ALL screen only covered T1; prereg B.8 requires every estimand cell
   vs its floor. Corrected in adjudicate_s1.py: T0/T1/T2/T3 cells vs their
   F7f floor cell (+0.03 R2); T0m cells vs the F7f-T0m Spearman floor
   (+0.10 Spearman tolerance; floor = F7f-T0m Spearman if defined, else 0).

## Frozen implementation freedoms clarified here (not covered by prereg)

- T3 rank loss operates on batch cells of one pair, margins between adjacent
  cells in target-sorted order (both forms). This is the realization of
  "pairwise logistic rank loss" frozen in B.4; it is now identical in
  Form-1 and Form-2.
- T0m ligand input = mean ECFP4 over the pair's ligand panel (the panels are
  shared, so this varies only slightly across pairs by design).

## Disclosure

- One full seed-1 run was executed with the pre-correction code
  (T2 already fixed, T3/T0m not). Its results are archived as
  SEED1_RESULT_DIAG_T3BUG.json and are NEVER adjudicated. Its T0/T1/T2 cells
  were computed by code paths untouched by this addendum, so the corrected
  run is expected to reproduce them numerically; the corrected run recomputes
  everything from one frozen code state regardless.
- No corrected-run result had been computed at the freeze time of this
  addendum (smoke at 2 pairs / 5 epochs only, used to verify the corrected
  code executes; its numbers are not arm results).

## Tests

tests/test_s1_structure.py extended: test_rank_loss_indices (asserts the
correct semantics AND that it differs from the buggy one), test_t0m_arm_
trains_scalar, test_severity_contrast_spearman. 16/16 PASS at freeze.
