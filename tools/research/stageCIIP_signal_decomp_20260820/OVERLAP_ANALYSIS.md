# OVERLAP ANALYSIS — CIIP-S1 vs concurrent CIIP-2 (governance record)

Date: 2026-08-20. Required by execution prompt section 0.3 and plan
compatibility note (Section 9.5) because CIIP-2 HAS produced training results.

## 1. State of the concurrent work (verified now, read-only)

- tools/research/stageCIIP2_olr_potential_20260820/ is COMPLETE and TERMINAL:
  Phase-4 single-seed smoke executed under frozen prereg a7b17e8a + ADD-1
  aa8d06af + ADD-2 91e2cb3a; gate (b) FAILED on real data (permutation
  control C-perm cell R2 0.1818 ranked above the best correct arm A4-cfoie
  0.1732); Phase 5 NOT executed per the frozen gate chain; final verdict
  UNRESOLVED (power); PHASE_REPORT.md written; ledger/task/history synced.
- No background jobs running (job list verified: all completed).
- No process is writing to any stageCIIP* directory (process census checked).
- CIIP-2 goal is closed; no pending or planned CIIP-2 training exists.

## 2. Overlap assessment

- Overlapping object: the T1 (centered contrast c) estimand on the same
  49-pair covered subset, same frozen DATA1A/DATA2X2 inputs.
- Non-overlapping objects: CIIP-S1 adds the T0/T0m estimands (proposition A,
  measured nowhere before), the T2 cross-fitted parent-residual estimand
  (C_sharp), T3 rank form, the Form-2 regression probes (F2-erased FIT, F3
  full-sequence deployable, F4 KLIFS), the F9 sibling-profile ceiling as a
  GATE (plan 9.5 gap (b)), and negative-control semantics distinct from
  CIIP-2's (same-parent wrong-mutation at feature level, ligand-label
  permutation as evaluative control).
- CIIP-2's terminal finding (no arm above the shared ligand pattern at the
  available power; instrument planted-effect ceiling ~ +0.03 R2 at n=39 train
  pairs) is INPUT to S1's power table, not a conflict: S1's B/A/C_sharp
  propositions were explicitly not tested by CIIP-2 (plan 9.5 gap (a),(c)).

## 3. Adjudication

Governance instruction received 2026-08-20 (user directive accompanying the
S1 execution prompt): the concurrent programme continues; the agent selects
the optimal direction per actual conditions and re-audits the concurrent code
for bugs before proceeding. Adjudication:

- CIIP-2 is terminal with no running or planned training; the §0.3 bar
  ("must not train in parallel") cannot be violated by S1 — there is nothing
  to run in parallel.
- Plan 9.5 explicitly designates S1 as a compatible strict prequel whose
  F9/sibling-ceiling arm closes CIIP-2's design gap (b).
- DECISION: CIIP-S1 proceeds. No CIIP-2 artifact is modified; CIIP-2 results
  enter S1 only as (i) prior power context and (ii) the audit facts in its
  report Section 2.4 (independent cross-validation of S0 diagnostics).

## 4. CIIP-2 code and training re-audit (user-mandated bug re-check)

Read-only re-review of olr.py, runner.py, gen_erased.py; frozen test suite
re-executed (PYTHONDONTWRITEBYTECODE=1; no writes to the frozen directory):

- tests/test_structure.py: 11/11 PASS (20 s). Recorded count discrepancy
  with the earlier session note ("12/12") is a test-count consolidation
  during that session's own edits; current file has 11 test functions, all
  green.
- crossfit_nuisance: parent-grouped folds; each pair's m_hat uses only fold
  models that excluded the pair's own parent from fitting — conservative,
  no label leak into any nuisance.
- gain_weights: train-parent WT rows only; the 40-step bisection is weight
  preprocessing (clip normalization), not a model solve; no closed-form
  model anywhere in the deployed path.
- train_arm: train-only loss; checkpoint selection on val only (ADD-2 full
  target with m_hat subtraction for A4/A5/C-perm); test labels touched only
  at frozen final evaluation; C-perm permutes its own arm's target
  (no cross-arm contamination).
- permute_within_pair: true derangement loop; paired_parent_bootstrap:
  parent-cluster resampling, [5,95] percentiles consistent with the frozen
  90% CI rule; bootstrap mean never used as a point estimate.
- Full-batch-per-epoch training (not minibatch) is CIIP-2's own frozen
  budget (AM-5), not a defect.
- gen_erased.py: erase_at replaces the single mutated position in BOTH
  sequences; pos asserted against the Q0B-verified pair table.
- No defect found that would change CIIP-2's terminal verdict
  (UNRESOLVED (power); no authorization).

## 5. S1 handling of CIIP-2 numbers

- CIIP-2 Section 2.4 facts (variance decomposition 134.8 / 89.7 %^2; sibling
  LOSO ceiling 0.293; family prior -0.021; ligand-pattern baseline 0.1313;
  out-of-range census 23.0%) are treated as independent secondary evidence
  to be re-derived by S0 from primary data under the S1 preregistration
  before any gate use.
