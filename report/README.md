# Reports

For a complete file-by-file classification, read [`REPORT_CATALOG.md`](REPORT_CATALOG.md).
This index is intentionally separate from the evidence ledger so navigation
does not alter scientific claims.

Use the following authority order. Older narrative reports remain for provenance;
they are not alternate task contracts.

1. `../task.md` — active performance/mechanism tracks and promotion gates.
2. `EVIDENCE_LEDGER.md` — compact stage index and numerical evidence pointers.
3. `CURRENT_MODEL_EVIDENCE.md` — retained model and comparator evidence.
4. `metasieve_research_programme.md` — latest theory-only programme review;
   design input, not an experiment result.
5. `LITERATURE_R15_20260819.md` — latest literature review.
6. `STAGE_X_ROUND1_REVIEW_20260818.md` and
   `measurement_pipeline_qualification.md` — Stage X measurement history.
7. Leaf reports under `tools/research/<stage>/` — authoritative results for
   that stage, each bound to its preregistration hash.

For the current mechanism status, read the two latest leaf reports directly:
`tools/research/stageCIIP_potential_bridge/CONTROL_REPORT.md` records
`ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED`, and
`tools/research/stageCIIP_context_propagation_20260820/CONTEXT_PROPAGATION_REPORT.md`
records representation-level context propagation with context-only prediction
not evaluated. Neither report authorizes CIIP-1B, the BindingDB bridge, or
production integration. Practical performance work belongs to Main Line P and
must be interpreted separately from the mechanism track.

The remaining root documents are historical or terminal summaries:
`BOUNDARY_20260816.md`, `BOUNDARY_20260817_NIGHT.md`,
`COMPLETION_STATEMENT_20260818.md`, `COMPLETION_STATEMENT_CORE_TASK1_20260817.md`,
`FINAL_STATE_20260818.md`, `LITERATURE_R14_20260816.md`,
`NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md`, and
`POST_COMPLETION_REVIEW_20260818.md`. Do not use them to override `task.md`
or the ledger.

The verification layer under `tools/research/stageN_audit/` is generated, never
hand-edited: `AUDIT_REPORT.md`, `FINAL_BOUNDARY_AUDIT.json` and
`COMPLETION_INVENTORY.json` are produced by `final_audit.py` and
`completion_inventory.py` from the artifacts themselves.

Every document above is scoped to **BindingDB-Ki double-cold development
evidence**. None of them generalizes to other DTA datasets or to architectures
that were not run, and none converts an empirical model failure into an
information-theoretic bound.

The four superseded, unreferenced narratives were removed after being merged
into `CONSOLIDATED_RESEARCH_RECORD.md`: the A2 handoff, the long Core Task 1
brief, the pre-R14 literature synthesis, and the earlier protein-conditioned
investigation. The remaining historical files stay because tests, audit
manifests, or leaf reports cite them.

`meta_fewshot/` retains formal R-series evidence. Pre-R0 experiments are
summarized once in `meta_fewshot/LEGACY_PRE_R0_SUMMARY.md`; Git preserves
removed detail without duplicate archive trees.
