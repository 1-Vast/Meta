# S1 ADDENDUM ADD-3 — Form-2 T3 rank-loss sign defect (frozen 2026-08-20)

Frozen BEFORE the final seed runs. Found by self-review immediately after the
ADD-2 corrected seed-1 run: the Form-2 T3 rank loss introduced by ADD-2 had a
flipped margin (h[a] - h[b] instead of h[b] - h[a], where b-positions hold
the LARGER target), so every Form-2 T3 arm was trained toward ANTI-ranked
output — visible as negative median Spearman on all Form-2 T3 cells while the
correctly-signed Form-1 T3 cell was positive.

- Fix: margin = h[b] - h[a] in the Form-2 T3 block (s1run.py), matching the
  Form-1 rank_loss helper exactly.
- Unit test added (test_form2_t3_rank_sign) asserting the aligned loss is the
  smaller one; suite 17/17 PASS at freeze.
- Impact: all Form-2 T3 cells of the ADD-2 seed-1 run are void. T3 is not a
  gating estimand for any frozen contrast (B/A/C_total/C_sharp/Deploy run on
  T1/T0m/T2), and Form-1 T3 was unaffected (helper was correct and unit
  tested). Nevertheless, per the one-frozen-code-state discipline, the
  ADD-2 seed-1 run is archived (SEED1_RESULT_ADD2_T3SIGNBUG.json) and ALL
  seeds are re-run from the final code state.
- No threshold, rule, arm, or estimand definition changes. The ADD-2
  compliance addendum remains in force; this addendum corrects a defect IN
  the ADD-2 implementation itself.
- Disclosure: at freeze time, the ADD-2 seed-1 T1/T0m/T2 cells had been seen
  (primary contrasts adjudicated UNRESOLVED except Deploy PASS
  power-labeled). Those cells' code paths are untouched by this sign fix;
  the re-run is expected to reproduce them numerically. The T3 fix is forced
  by the defect, not by any T3 result threshold.
