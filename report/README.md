# Reports

Read the current record in this order:

1. `../task.md` — active five-phase execution contract and gates.
2. `POST_COMPLETION_REVIEW_20260818.md` — interpretation and governance review.
3. `BOUNDARY_20260817_NIGHT.md` — measured BindingDB development boundary.
4. `CURRENT_MODEL_EVIDENCE.md` — compact model and mechanism evidence.
5. `EVIDENCE_LEDGER.md` — stage decisions and leaf artifact locations.
6. `FINAL_STATE_20260818.md` and `COMPLETION_STATEMENT_20260818.md` — closure of
   the preceding method-ladder cycle, not the active plan.

The verification layer under `tools/research/stageN_audit/` is generated, never
hand-edited: `AUDIT_REPORT.md`, `FINAL_BOUNDARY_AUDIT.json` and
`COMPLETION_INVENTORY.json` are produced by `final_audit.py` and
`completion_inventory.py` from the artifacts themselves.

Every document above is scoped to **BindingDB-Ki double-cold development
evidence**. None of them generalizes to other DTA datasets or to architectures
that were not run, and none converts an empirical model failure into an
information-theoretic bound.

`BOUNDARY_20260816.md`, the A2 plan/handoff, and R-series plans are historical
inputs retained because tests and exact operator definitions cite them. They are
not active execution authorities. Literature context remains
in `LITERATURE_RESEARCH_SYNTHESIS_20260815.md` and
`LITERATURE_R14_20260816.md`.

`meta_fewshot/` retains formal R-series evidence. Pre-R0 experiments are
summarized once in `meta_fewshot/LEGACY_PRE_R0_SUMMARY.md`; Git preserves
removed detail without duplicate archive trees.
