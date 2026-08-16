# MetaSieve Repository Organization

Last audited: 2026-08-16

This document is the preflight map for Codex, Claude, and human reviewers. It
classifies files by runtime authority and prevents research prototypes,
historical theory, and failed experiment artifacts from being mistaken for the
active model.

## Folder map

| Area | Contents | Entry point | Policy |
|---|---|---|---|
| `model/` | Active neural modules plus tested compatibility operators | `model/README.md` | Runtime authority |
| `scripts/` | Active DTA, governed data, optional geometry commands | `scripts/README.md` | CLI/workflow authority |
| `tests/` | 82 maintained pytest modules, three tiers | `tests/README.md` | Must pass before promotion |
| `test/` | Manual sealing and audit commands | `test/README.md` | Not a second pytest tree |
| `research/` | Five code-only prototype families | `research/README.md` | Never production authority |
| `report/` | 336 `RESULT.json` files and narrative indexes | `report/README.md` | Immutable evidence |
| `dataset/` | ~69 GB raw/processed/sealed/cache assets | `dataset/README.md` | Local governed data |
| `contracts/` | Stable graph, context, mechanism, source schemas | `contracts/README.md` | Schema authority |
| `config/` | Defaults and remote execution helpers | `config/README.md` | Not result authority |
| `docs/` | Ownership and cleanup decisions only | `docs/README.md` | Current navigation |
| `archive/` | Historical theory, implementations, reports | `archive/README.md` | Curated provenance |
| `LLM/` | Local external-review gateway helpers | `LLM/README.md` | Credentials stay ignored |
| `weights/` | Pointer to remote model weights | `weights/README.md` | No local payload |
| `tools/` | Ignored downloaded binaries | none | Reinstallable local tools |

## Authority levels

### A0: active runtime authority

- `main.py`: retained command dispatcher.
- `model/interaction_grammar.py`: active contact-grammar trunk and
  transferability-gated transport (`--arch grammar`).
- `model/similarity_grammar.py`: grammar trunk with the chemistry-grounded
  Tanimoto transport (`--arch similarity` / `--arch similarity_only`). The
  frozen Stage R3/R4 `similarity_only` checkpoints are the **retained
  incumbent A0**.
- `model/level_shape.py`: level-shape factorized predictor (Stage R3/R4
  Innovation A), trained by `scripts/train_level_shape.py`.
- `model/reltransport.py`: the R5-R8 relative-transport candidate (bilinear
  relative potential + linear rho gate + attention-pooled level), trained by
  `scripts/train_reltransport.py`. The model family was **closed** on
  2026-08-16 under its preregistered gates (see `report/` R7/R8); the module
  and its 23-gate suite remain for evidence.
- `model/bpsf.py`: retained bipartite protein-ligand interaction trunk and
  readout, used by the `--arch bpsf` control arm.
- `model/encoders.py`: protein-slot and 2D ligand encoders.
- `model/qpsmp_meta.py`: retained zero-shot head and label-locked residual
  kernel, plus the shared `QPSMPMetaOutput` contract.
- `scripts/qpsmp_data.py`: active episode sampling/materialization contract
  (includes the 2026-08-16 physical `meta_test` seal,
  `include_meta_test=False`).
- `scripts/train_qpsmp.py`: active single-stage episodic trainer
  (`--eval-meta-test`/`--include-meta-test` opt-in, default off).
- `scripts/train_grammar_shape.py`: shape-first routed training applied to
  the incumbent grammar trunk (Stage R10 candidate; Innovation B test with
  zero architecture change).
- `scripts/evaluate_qpsmp.py`: governed nested-k evaluation and controls
  (`--split` default `meta_val`; `meta_test` needs explicit opt-in).
- `scripts/evaluate_checkpoint_nested.py`: fixed-bank checkpoint re-evaluation.
- `scripts/stageR6_compare_arms.py`: three-kind paired arm comparison with
  novelty tiers, activity-cliff sign accuracy, wrong-ligand control, and
  component bootstraps (used by R6-R10).
- `scripts/stageR9_pair_audit.py`: no-training pair-level CI decomposition
  (stratum counts/sign accuracy/margins/loss + component bootstrap).
- `contracts/ligand_graph.py`: active atom/bond dimensional contract.
- `tests/test_bpsf.py`, `tests/test_qpsmp_meta.py`,
  `tests/test_train_qpsmp_meta.py`, `tests/test_term_synthetic.py`,
  `tests/test_interaction_grammar_synthetic.py`,
  `tests/test_qpsmp_nested.py`, `tests/test_level_shape.py`,
  `tests/test_reltransport_synthetic.py`,
  `tests/test_stage0_contract_fixes.py`: minimum active regression suite.
- `scripts/stage0_diagnostics.py`, `scripts/stage0_reference_baselines.py`,
  `scripts/stage0_throughput.py`, `scripts/stage0_trunk_capacity.py`:
  read-only audit instruments (gradient coverage, label-only reference
  predictors, cost, synthetic trunk capacity).
- `report/CURRENT_MODEL_EVIDENCE.md`: consolidated experiment decision ledger.

The governed double-cold protocol
`dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1/` (assignment
sha256 frozen) is the development/confirmation population of the R5-R13
cycle; `meta_test` is sealed physically by default.

**Two different `meta_test` populations exist in this repository and must
never be conflated.** The pre-double-cold `bindingdb_ki_main_v0` `meta_test`
(42 targets / 6-7 components) was legitimately reported and consumed by
Stages 4/6/7. The double-cold `bindingdb_ki_double_cold_v1` `meta_test`
(22 targets / 10 components) has never been opened. Any `meta_test` number in
a Stage 4/6/7 report belongs to the older, consumed population and is not
evidence about the sealed split.
`python -m scripts.audit_research_record` classifies every `RESULT.json`
into double-cold-sealed (explicit R5 seal record), double-cold-sealed
(implicit, pre-R5 artifacts recording only `meta_val`), or older-protocol,
and fails if any double-cold artifact reports `meta_test` metrics.

### A1: optional active extensions

- `model/cartesian.py`: optional sparse Cartesian encoder. The active BindingDB
  path supplies no common-frame coordinates.
- `model/bpsf.py`: source-only contact/distance supervision is co-located with
  the pair field it consumes.
- `scripts/pretrain_pair_geometry.py` and `scripts/evaluate_pair_geometry.py`:
  geometry-head commands, not part of current BindingDB training.
- `scripts/build_ligand_bank.py`, `scripts/build_protein_bank.py`, and
  `scripts/precompute_qpsmp_compact_bank.py`: active bank builders.

### A2: governed data pipeline

- Acquisition: `scripts/acquire_*.py`, `scripts/source_affinity/`, and
  `scripts/structure_sources/`.
- Canonicalization: `scripts/preprocess_dataset.py`,
  `scripts/data_contract.py`, and `scripts/data_protocols.py`.
- Governance: `scripts/govern_affinity_homology.py`,
  `scripts/govern_structure_homology.py`, and `scripts/governed_sources.py`.
- Verification/sealing: `scripts/verify_*.py`,
  `scripts/seal_compiled_dataset.py`, and `scripts/project_status.py`.
- Structure-specific builders: `scripts/build_holo_complex_index.py`,
  `scripts/build_structure_supervision.py`, and `scripts/cache_*structure*.py`.

These files must not be mixed into model-architecture changes. A stage that
changes sampling or data construction is a data-protocol experiment.

### R0: current research prototypes

- `research/meta_fewshot/`: earlier BPSF/meta-section candidates, corpus
  builders, and oracle gates.
- `research/crossed_interaction/`: C/Q observables, SAR-delta studies,
  localization diagnostics, and rectangle/low-rank experiments.
- `research/correspondence_router/`: structure correspondence panels.
- `research/e0_identifiability/`: tensor-basis/directional identifiability.
- `research/source_affinity/`: source-assay SAR-delta experiments.

Research modules are imported by many tests. Do not move them until imports and
artifact paths are migrated in one dedicated cleanup change.

### R1: experiment results

`report/meta_fewshot/` contains the complete model lineage. Its directory name
is part of provenance and must remain stable.

Current retained incumbent (the comparator every R5-R13 gate is measured
against):

- `stageR3R4_level_shape_20260815/A0_incumbent_seed202608{15,16,17}/` —
  `similarity_only` grammar checkpoints, 1200 steps (double-cold k=0
  2.149, CI 0.580, three seeds).

The R5-R13 cycle (2026-08-16):

- `stageR5_reltransport_20260816/` — contract repairs + 23 structural gates
  (all pass).
- `stageR6_reltransport_screening_20260816/` — three screening falsifications
  (R6a/R6b/R6c archives retained; probes `probe_w0.1`, `probe_w0.3`,
  `smoke_v2` are generated diagnostics).
- `stageR7_reltransport_3seed_20260816/` — formal three-seed run; admission
  refused.
- `stageR8_stronger_shape_20260816/` — stronger shape signal; family closed
  for the double-cold zero-shot target; best-ever shape 0.896 and k=5 cliff
  sign 0.768 retained as findings.
- `stageR9_cliffweight_20260816/` — cliff pair-weight dose response and the
  pair-level CI audit.
- `stageR10_variance_20260816/` — shape-variance reduction; falsified.
- `stageR11_grammar_shape_20260816/` — shape-first routing on the incumbent
  trunk (zero architecture change); falsified.
- `stageR12_margin_20260816/` — margin-ranking shape objective; falsified as
  the actionable lever. `REPORT.md`/`RESULT.json` backfilled 2026-08-16 from
  the retained comparison artifact.
- `stageR13_shape_direct_20260816/` — direct interaction-head shape;
  gate-blocked at Stage 1 (18 gates: 16 pass, 2 `xfail`). No real-data run.
- `report/BOUNDARY_20260816.md` — the consolidated reachable-boundary
  statement and the correctly defined k=0 Pareto frontier.

Every result directory should contain `RESULT.json`; a training run may also
contain `checkpoint.pt`, `progress.jsonl`, and prediction manifests. Never
delete a failed result merely because its source branch was reverted.

### H0: historical and compatibility code

- `model/bands.py`, `model/config.py`, `model/mathematical.py`,
  `model/meta_operator.py`, and `model/runtime.py`: older operator diagnostics
  still covered by compatibility tests, not the active DTA path.
- `archive/legacy/retired_qpsmp/`: retired model implementation, moved out of
  the active package surface on 2026-08-15.
- `test/`: older test surface; maintained tests live in `tests/`.
- `archive/`: historical theory, curated FORT source/reports, and retired
  research. It is not current authority. Embedded Git databases, downloaded
  third-party trees, credentials, caches, and reproducible quarantine outputs
  are not historical evidence and are removed.

Do not delete compatibility modules during performance work. Retire them only
after their tests and imports are explicitly removed in a separate change.

## Data classification

### Active BindingDB task assets

- `dataset/processed/meta_fewshot/bindingdb_ki_main_v0/`
- `dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank/`
- `dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank/`
- `dataset/processed/meta_fewshot/bindingdb_ki_main_v0_ligand_bank_compact/`

### Other governed/experimental assets

- `dataset/processed/open_structures/`: optional structure-supervision corpora.
- `dataset/processed/source_affinity/`: ChEMBL/source-affinity work.
- `dataset/processed/correspondence_router/` and `crossed_interaction/`:
  research panels.
- `dataset/sealed/`: sealed DAVIS assets.
- `dataset/raw/`: downloaded source data; never edit during training.

An `.incomplete_*` directory is not an authorized model input.

## Documentation authority

Read in this order:

1. `README.md`: repository entry and scope.
2. `docs/PROJECT_FILE_ORGANIZATION.md`: file ownership and active surface.
3. `PROJECT_SUMMARY.md`: current model/result summary.
4. `task.md`: experimental constraints.
5. `report/CURRENT_MODEL_EVIDENCE.md`: current evidence and rejected stages.
6. `report/*/RESULT.json`: numerical authority.
7. `archive/theory/`: historical mathematical provenance only.

`history.md` and dated reports are provenance, not current architecture
authority.

## Archive retention and cleanup

Retain:

- unique source, reports, manifests, decision records, and failed-result
  evidence needed to reconstruct scientific conclusions;
- `archive/theory/` while its existing migration remains under review;
- curated FORT model/research/report material and its historical `history.md`;
- failed checkpoints only when a retained result or audit explicitly depends
  on the checkpoint rather than its recorded metrics.

Delete:

- nested `.git/` databases inside snapshots;
- archived `.env` or other credential-bearing local configuration;
- `tmp/` copies of third-party repositories, downloaded datasets, and smoke
  products when they are not the sole evidence for a result;
- `ephemeral_quarantine/`, `__pycache__/`, `.pytest_cache/`, and other
  reproducible caches;
- byte-identical duplicate artifacts after the retained copy and all path
  references have been recorded.

The 2026-08-15 cleanup removed the archived plaintext `.env`, the 1.33 GB FORT
temporary tree, the retired-research quarantine, and root runtime caches. It
also consolidated the root `legacy/` implementation under `archive/legacy/`.
The nested `archive/FORT/.git/` was also removed after exact-path relocation.
Only an ACL-protected, empty-for-auditing pytest temporary tree remains under
`.tmp_pytest_final_cleanup_20260815/`; it is not retained scientific evidence.
See `history.md` for the exact status.


## 2026-08-16 audit notes (R5-R13 cycle)

- **Active surface additions**: `model/reltransport.py`,
  `model/shape_direct.py`, `scripts/train_reltransport.py`,
  `scripts/train_level_shape.py`, `scripts/train_grammar_shape.py`,
  `scripts/train_shape_direct.py`, `scripts/stageR6_compare_arms.py`,
  `scripts/stageR9_pair_audit.py`, `scripts/stage_smoke.py`,
  `scripts/run_stage.py`, `scripts/audit_research_record.py`,
  `tests/test_reltransport_synthetic.py` (23 gates),
  `tests/test_shape_direct_synthetic.py` (18 gates: 16 pass, 2 `xfail`),
  `tests/test_stage0_contract_fixes.py`, `tests/test_stage_runner.py`,
  `tests/test_research_record.py`.
- **Evidence**: R5-R13 stage directories under `report/meta_fewshot/`
  (RESULT/PREDICTIONS/checkpoints/PREREGISTRATION are immutable;
  `stageR6_*/{R6a,R6b,R6c}_archive/` hold eliminated screening runs;
  `probe_w0.1`, `probe_w0.3`, `smoke_v2` are **generated** diagnostics).
- **Record-integrity check**: `python -m scripts.audit_research_record`
  recomputes the arm table, the k=0 Pareto frontier, the `meta_test` seal
  classification and every recorded checkpoint sha256 from the artifacts.
  As of 2026-08-16: 78 checkpoint hashes verified, 0 mismatched;
  75 explicitly sealed + 15 implicitly sealed double-cold artifacts,
  138 older-protocol artifacts, 0 seal violations.
- **Generated/recoverable**: `__pycache__/`, `.pytest_cache/`, the
  `.tmp_pytest_final_cleanup_20260815/` tree, `tools/downloads/`,
  `tools/cdhit/`, `tools/mmseqs2/` (reinstallable local tools).
- **Import dependency**: `research/` modules are imported by
  `scripts/cache_exact_structure_proteins.py`,
  `scripts/cache_r0b_exact_geometry.py`, and several `tests/test_*_cq_*.py`
  modules — `research/` stays in place until a dedicated migration change.
- **No deletions executed in this pass**: every model module is referenced
  by scripts, tests, or checkpoint loaders; all stage artifacts are retained.

## Safe cleanup backlog

### Can be cleaned without changing repository architecture

- Ensure local API-key files under `LLM/` remain ignored and uncommitted.
- Generate experiment indexes from result metadata without renaming artifacts.

### Requires a dedicated import migration

- Split `scripts/` into `data/`, `training/`, and `evaluation/` packages.
- Split the 81 test modules by active/runtime versus research lineage.
- Move compatibility operator modules out of `model/`.

Completed on 2026-08-15: `research/meta_fewshot/results/` was consolidated
under `report/meta_fewshot/bpsf_v2_research/`. Old narrative reports, runtime
logs, a stale implementation contract, and raw theory-input transcripts were
summarized into `report/CURRENT_MODEL_EVIDENCE.md` and then deleted.

### Prohibited during model optimization

- Moving or deleting active BindingDB banks.
- Renaming result directories or modifying historical `RESULT.json` files.
- Bulk-moving `research/` modules while tests import them.
- Restoring or deleting the already-relocated `archive/theory/` tree before
  its current user-owned migration is reconciled.
- Deleting failed checkpoints before their metrics are indexed.

## Codex/Claude preflight checklist

Before proposing a new architecture stage:

1. Run `git status --short` and preserve unrelated changes.
2. Read the A0 files and the current decision ledger.
3. Confirm failed experimental branches are absent from the active source.
4. Confirm `representation_warmup_fraction == 0.0` and zero means zero steps.
5. Run the focused active regression suite in the `drug` environment.
6. Record the frozen three-seed Stage E metrics.
7. Write one hypothesis, one changed variable, and quantitative gates.
8. Do not move files during the model experiment.
9. Preserve all run outputs under a new, uniquely named result directory.
10. On failure, revert only stage-specific source/test changes and update the
    decision ledger.

The complete maintained suite is 82 pytest modules. The minimum verification
command is:

```powershell
D:\anaconda\envs\drug\python.exe -m pytest tests -q
```

## Regression-suite tiers

| tier | how to run | contents | 2026-08-16 outcome | wall |
|---|---|---|---|---:|
| default | `pytest tests -q` | every structural, algebraic, gradient, contract, data and record-integrity gate, including all of the closed families' *non-training* gates | 412 passed, 9 skipped | 103 s |
| research gates | `RUN_RESEARCH_GATES=1 pytest tests -q` | adds the synthetic **training** gates of the two closed families (`test_reltransport_synthetic.py`, `test_shape_direct_synthetic.py`) | 416 passed, 3 skipped, 2 xfailed | 410 s |
| slow smokes | `RUN_SLOW=1 pytest tests -q` | adds the full-corpus subprocess smokes in `test_stage0_contract_fixes.py` | — | — |

The two `xfail`s are the recorded R13 gate verdicts, not failures.

The research-gate tier exists because those six tests train small synthetic
models on CPU and cost ~314 s, and because both families are **already
decided** (R8 closed the relative-transport family; R13 gate-blocked the
direct-shape family). Their measured verdicts are preserved in
`report/meta_fewshot/stageR8_stronger_shape_20260816/` and
`report/meta_fewshot/stageR13_shape_direct_20260816/RESULT.json`, so the
evidence does not depend on re-running them. **A new or reopened model family
must run its own gates in the research tier before any real-data training** —
skipping them by default is a cost decision about settled questions, not a
lowering of the admission bar.
