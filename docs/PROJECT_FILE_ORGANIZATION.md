# Project file organization

Updated 2026-08-19 for the split performance/mechanism programme.

| Directory | Purpose | Authority |
|---|---|---|
| `model/` | Active and retained comparator models | `model/README.md` |
| `scripts/` | Governed data, training, evaluation and R-series analysis | `scripts/README.md` |
| `dataset/` | Local governed source and processed data | manifests inside each corpus |
| `report/` | Compact evidence authority and R0-R14 leaf artifacts | `report/README.md` |
| `archive/` | Recovery index only; no duplicate source tree | `archive/README.md` |
| `tools/research/` | Unadmitted research implementations, probes and stage artifacts | stage preregistration |
| `tools/tests/` | Active contracts for retained code | `tools/tests/README.md` |
| `tools/runtime/` | Ignored third-party executables and local helpers | tool-local metadata |

Root `main.py` orchestrates admitted model functionality through `scripts/`.
Research is intentionally absent from that command surface until promotion.

## Model surface

`bpsf.py`, `qpsmp_meta.py`, `encoders.py` and `interaction_grammar.py` retain the
core interaction path. `similarity_grammar.py` is the A0/Tanimoto comparator.
`level_shape.py` and `reltransport.py` reproduce the formal R-series Pareto and
transport arms. `cartesian.py` remains an algebraically tested optional module
but is not reachable from the current coordinate-free DTA dataset.

Closed `relative_grammar`, `locality_grammar`, `shape_direct`, frozen mathematical
operator and old config modules were removed with their dedicated trainers/tests.

## Report surface

At root, `task.md` is the active execution authority and `history.md` is the
chronology. `report/` keeps compact evidence authorities, historical definitions
still cited by tests/operators, and formal `stageR0`-`stageR14` leaf artifacts.
A closed proposal is not an active plan even when retained for reproducibility.

## Research surface

`tools/research/<stage_name>/` holds unadmitted probes; its tests live in
`tools/research/<stage_name>/tests/` and never in `tools/tests/`, which is
reserved for contracts on admitted or repository-wide behaviour.

| stage | what it decided |
|---|---|
| `a2_readiness/` | consolidated to `SUPERSEDED.md`; retains the literature ledger and the CPC centering probes |
| `a2_readiness_v2/` | governance incident, noise/leakage audit, causal attention audit, Stage P's frozen design. Its A2 verdict is superseded |
| `a2_exact_probe/` | **the superseding stage.** `FINAL_DECISION.md` is the single consolidated A2 record; `STAGE_R_EXACT_A2.md` closes A2 on its own operator; `STAGE_L_LIGAND_SAR.md` measures the ligand-side direction against Tanimoto |
| `stageS_sar_field/` | rejected global protein-conditioned SAR field; retains the decisive shuffled-protein and counterfactual-degeneracy evidence |
| `stageT_mmp/` | historical true-MMP probe; retained for scoped evidence and correction history |
| `stageP_practical_fewshot/` | practical target-cold performance track; frozen three-layer evaluation and baseline bake-off |
| `stageX_csc_signal/` | current mechanism qualification track; Q2d synthetic stages and matched-variant instruments |

All other `stage*` directories are historical evidence. Their source is not
production code and must not be imported by `main.py`.

Frozen-feature caches (50 MB) were deleted after recording their sha256 in
`a2_exact_probe/FEATURE_CACHE_MANIFEST.json`; the `*.meta.json` provenance
sidecars are retained beside it. Regenerate with the `extract_*` commands.
Only the `.json` results and the decision reports are evidence.

## Test surface

Test counts are reported by each completed stage and must not be copied forward
as a timeless repository total. Run the maintained and active-stage suites for
the current count.

## Recovery

Git is the archive. `archive/README.md` names recovery commits for deleted source,
reports and historical theory. Do not recreate physical archive mirrors in the
working tree. A restored artifact must be copied into a temporary branch, audited,
and promoted through a new preregistered stage rather than silently returned.

## Retention rule

For each new stage retain only its preregistration, machine result, decision
report, necessary prediction rows and loadable admitted checkpoint. Delete
generated caches, duplicate smokes, progress logs and failed checkpoints after
the stage process has exited and after recording the verdict. Never delete a
leaf artifact cited by `task.md`, `history.md` or `EVIDENCE_LEDGER.md`. Update
`history.md`, `task.md`,
`CURRENT_MODEL_EVIDENCE.md` and the ledger.

When research passes, move its reusable modules into `model/`, its commands into
`scripts/`, and its contracts into `tools/tests/`; then remove the original
`tools/research/` implementation. Git history is the recovery mechanism.
