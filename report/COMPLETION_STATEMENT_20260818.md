# Completion statement for the MetaSieve long-run objective

Date: 2026-08-18 (night). This statement is the evidence authority for
marking the objective recorded in tools/research/GOAL_ACTIVE.md complete.

## The objective's two terminal conditions

1. A reproducible performance leap (MSE <= 1.00 pK^2 at k=0/1/2/3/5 with
   CI/Spearman/Pearson/centered-MSE/cliff preserved under the governed
   cold-target protocol), or
2. A scope-bounded final conclusion after multi-family falsification.

## Which condition is met

**Condition 2 is met.** Condition 1 was not achieved by any candidate;
the measured evidence establishes why under this protocol:

- k=0 MSE = level^2 + centered; the level term is assay-history-dominated
  (within-document transfer R^2 +0.451; the split's document closure makes
  that signal unavailable at inference; legal transferring inputs cover at
  most ~26% of level variance). The best trained zero-shot level^2 on
  record is 1.2151 against a 0.1239 budget; per-seed k=0 MSE never went
  below 2.10.
- k>=1 approaches the target (k=5 at 0.939-1.007 across seeds with honest
  controls) and K-REG produced the first all-k resolved MSE improvement
  across three seeds, but its shape gain did not survive pooling, so
  nothing was promoted.
- The level/ranking conflict on one shared trunk is fundamental to
  single-stage end-to-end training (compositions E, J, L, Q all failed the
  ranking gate), and the only escape is a multi-stage calibrator the
  governing contract excludes.

## Evidence chain (all artifacts exist; completion inventory verified)

- Conclusion: report/BOUNDARY_20260817_NIGHT.md; closing summary:
  report/FINAL_STATE_20260818.md.
- Falsification ledger across four axes (frameworks, training families,
  external representations, level-head compositions):
  report/EVIDENCE_LEDGER.md; method-ladder closure map:
  tools/research/method_ladder/CLOSURE_MAP.md.
- Bitwise verification audit: tools/research/stageN_audit/
  FINAL_BOUNDARY_AUDIT.json and AUDIT_REPORT.md; completion inventory:
  tools/research/stageN_audit/COMPLETION_INVENTORY.json.
- Independent external validation: Nelen et al., J Cheminform 17:8, 2025
  (absolute values not comparable across assays; paired differences robust).
- Governance: no sealed meta_test label entered any fitting, selection or
  reported metric (0 evaluations in 104 audited artifacts); all training
  stages preregistered; maintained suite 268 passed / 6 skipped; research
  suites 151 passed / 12 skipped; Git commits 361c342..59fae9e.

## Scope of the conclusion (when it would change)

The conclusion applies to the governed BindingDB-Ki double-cold protocol
with sequence + 2D ligand inputs, single-stage differentiable training, and
the locally available legal input families. It would change with: a
governed UniRef snapshot (MSA lane), wider structure coverage, a
Davis/KIBA authorization (frozen plan:
tools/research/stageR_daviskiba/PREREGISTRATION.md), or a protocol
restatement (per-document calibration or centered/ranking targets).

## Marking status

The originating session did not register the goal tools
(create_goal/get_goal/update_goal), so the formal complete transition
could not be executed there. This statement plus GOAL_ACTIVE.md are the
durable evidence; a session with goal tools should read them and mark the
objective complete on the conclusion branch.
