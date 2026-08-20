# Report Catalog

This catalog is the entry point for `report/`. Legacy paths remain stable
because reports are cited by `task.md`, `history.md`, stage artifacts, scripts,
and tests. Organization is therefore logical first; physical moves require a
repository-wide reference update and a separate review.

## Read First

1. `../task.md` — active research contract and stage gates.
2. `EVIDENCE_LEDGER.md` — chronological evidence and claim status.
3. `CURRENT_MODEL_EVIDENCE.md` — current model behavior and verified limits.
4. `README.md` — report navigation and retention policy.

## Authority And Governance

These files define the current state or make it auditable. Do not delete,
rename, or move them without updating every reference and the evidence ledger.

| File | Role |
| --- | --- |
| `EVIDENCE_LEDGER.md` | Evidence index, stage decisions, and claim provenance |
| `CURRENT_MODEL_EVIDENCE.md` | Current accepted baseline and limitations |
| `CORE_TASK1_UNRESOLVED_TERMINAL_20260817.md` | Core Task 1 terminal interpretation |
| `CORE_TASK1_UNRESOLVED_TERMINAL_20260817.json` | Machine-readable terminal decision |
| `COMPLETION_EVIDENCE_MANIFEST_CORE_TASK1.json` | Hash manifest for completion evidence |
| `RECORD_AUDIT.json` | Repository/report audit record |

## Current Research Inputs

These are active design inputs or recent audits, not automatically accepted
results. They must remain separate from authority files.

| File | Role |
| --- | --- |
| `metasieve_research_programme.md` | Current research-program design and proposed tracks |
| `LITERATURE_R15_20260819.md` | Latest literature synthesis |
| `measurement_pipeline_qualification.md` | Measurement/data qualification findings |
| `STAGE_X_ROUND1_REVIEW_20260818.md` | Stage X round-1 review |
| `CONSOLIDATED_RESEARCH_RECORD.md` | Consolidated replacement for superseded root narratives |

## Historical Summaries

These documents are retained for traceability and are not active model
specifications.

| Files | Role |
| --- | --- |
| `BOUNDARY_20260816.md`, `BOUNDARY_20260817_NIGHT.md` | Earlier reachable-boundary statements |
| `COMPLETION_STATEMENT_20260818.md`, `COMPLETION_STATEMENT_CORE_TASK1_20260817.md` | Completion statements |
| `FINAL_STATE_20260818.md`, `POST_COMPLETION_REVIEW_20260818.md` | Post-cycle summaries |
| `NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md` | Superseded plan retained because audit records cite it |
| `LITERATURE_R14_20260816.md` | Earlier literature review retained because a leaf report cites it |

## Stage Evidence

- `mechanism/` contains mechanism-track stage directories and their
  preregistrations, results, and reports.
- `meta_fewshot/` contains few-shot/performance-track stage directories.
- A stage directory is authoritative only through its own frozen
  preregistration, result artifact, and report; root summaries are navigation
  documents, not replacements for leaf evidence.

## Cleanup Rules

- Never delete a report solely because it is old; first run `rg` over the
  repository and update all citations.
- Do not remove untracked files owned by another agent.
- Do not delete active-stage outputs or caches while a training process is
  running. After the process exits, cache cleanup requires explicit path
  validation and a recorded command.
- New reports should be placed in the relevant stage directory. Add a root
  summary only when it is a genuine authority or navigation document.
- Keep raw restricted datasets outside Git; retain manifests, hashes, and
  semantic audits only.

The following four unreferenced root narratives were removed after
consolidation: `AGENT_HANDOFF_A2_MOMENT.md`,
`core_task_1_protein_conditioned_signal.md`,
`LITERATURE_RESEARCH_SYNTHESIS_20260815.md`, and
`protein_conditioned_signal_investigation.md`.
