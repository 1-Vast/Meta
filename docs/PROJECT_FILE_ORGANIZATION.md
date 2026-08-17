# Project file organization

Updated 2026-08-16 (A2-readiness v2 stage; the 2026-08-16 consolidation layout
is otherwise unchanged).

| Directory | Purpose | Authority |
|---|---|---|
| `model/` | Active and retained comparator models | `model/README.md` |
| `scripts/` | Governed data, training, evaluation and R-series analysis | `scripts/README.md` |
| `dataset/` | Local governed source and processed data | manifests inside each corpus |
| `report/` | Compact evidence authority and R0-R14 leaf artifacts | `report/README.md` |
| `archive/` | Recovery index only; no duplicate source tree | `archive/README.md` |
| `tools/research/` | Unadmitted research implementations and probes | stage preregistration |
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

At root, keep only high-density authorities: `CURRENT_MODEL_EVIDENCE.md`,
`EVIDENCE_LEDGER.md`, `BOUNDARY_20260816.md`, the next plan and report indexes.
Under `report/meta_fewshot/`, only formal `stageR0`-`stageR14`, independent
`stageM0`, and the legacy summary remain. Pre-R0 smokes and exploratory chains
were deleted after consolidation.

## Research surface

`tools/research/<stage_name>/` holds unadmitted probes; its tests live in
`tools/research/<stage_name>/tests/` and never in `tools/tests/`, which is
reserved for contracts on admitted or repository-wide behaviour.

| stage | what it decided |
|---|---|
| `a2_readiness/` | consolidated to `SUPERSEDED.md`; retains the literature ledger and the CPC centering probes |
| `a2_readiness_v2/` | governance incident, noise/leakage audit, causal attention audit, Stage P's frozen design. Its A2 verdict is superseded |
| `a2_exact_probe/` | **the superseding stage.** `FINAL_DECISION.md` is the single consolidated A2 record; `STAGE_R_EXACT_A2.md` closes A2 on its own operator; `STAGE_L_LIGAND_SAR.md` measures the ligand-side direction against Tanimoto |

Frozen-feature caches (50 MB) were deleted after recording their sha256 in
`a2_exact_probe/FEATURE_CACHE_MANIFEST.json`; the `*.meta.json` provenance
sidecars are retained beside it. Regenerate with the `extract_*` commands.
Only the `.json` results and the decision reports are evidence.

## Test surface

`tools/tests/` — 249 collected (6 skipped by default: 3 `slow`, 3
`research_gate`; enable with `RUN_SLOW=1` / `RUN_RESEARCH_GATES=1`).
`tools/research/*/tests/` — 53 stage probes. Full sweep: **302 passed**.

## Recovery

Git is the archive. `archive/README.md` names recovery commits for deleted source,
reports and historical theory. Do not recreate physical archive mirrors in the
working tree. A restored artifact must be copied into a temporary branch, audited,
and promoted through a new preregistered stage rather than silently returned.

## Retention rule

For each new stage retain only its preregistration, machine result, decision
report, necessary prediction rows and loadable admitted checkpoint. Delete
duplicate smokes, progress logs and failed checkpoints after recording the
verdict. Update `history.md`, `task.md`, `CURRENT_MODEL_EVIDENCE.md` and the ledger.

When research passes, move its reusable modules into `model/`, its commands into
`scripts/`, and its contracts into `tools/tests/`; then remove the original
`tools/research/` implementation. Git history is the recovery mechanism.
