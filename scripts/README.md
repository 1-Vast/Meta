# Script Surface

Use the root `main.py` for public orchestration. It exposes data preparation
and verification while keeping failed research runners out of the command
surface.

`scripts/` contains the active trainer/evaluator plus governed data and optional
geometry workflows.

## Active DTA

- `qpsmp_data.py`: episode sampling and materialization.
- `train_qpsmp.py`: single-stage episodic training. `--arch {bpsf,grammar}`
  selects the architecture; `--lr-schedule cosine` and
  `--val-targets-per-component` control the schedule and the selection bank
  without touching the frozen protocol test bank.
- `evaluate_qpsmp.py`: governed nested-k multi-seed evaluation with paired
  component bootstrap; saves one checkpoint per model seed.
- `evaluate_checkpoint_nested.py`: fixed-bank checkpoint re-evaluation.
  `--targets-per-component 999999` produces the wide bank over all eligible
  meta-test targets; omit it for the frozen 6-episode protocol bank.

### Stage 0 audit instruments (read-only)

- `stage0_diagnostics.py`: per-k gradient coverage, dead-branch census,
  activation scales, parameter census, protein-swap sensitivity.
- `stage0_reference_baselines.py`: label-only reference predictors
  (global mean, ligand prior, support mean, oracle target mean).
- `stage0_throughput.py`: per-architecture seconds/step and peak memory.
- `stage0_trunk_capacity.py`: synthetic protein-by-ligand bilinear capacity
  probe used to separate expressivity failures from optimization failures.
- `build_ligand_bank.py`, `build_protein_bank.py`,
  `precompute_qpsmp_compact_bank.py`: active feature-bank builders.

## Canonical data and sealing

- `data_contract.py`, `data_protocols.py`, `preprocess_dataset.py`
- `build_ligand_bank.py`, `build_protein_bank.py`, `data.py`
- `seal_compiled_dataset.py`, `governed_sources.py`, `runtime.py`

## Open structure geometry

- `acquire_open_structures.py`, `build_holo_complex_index.py`
- `build_structure_supervision.py`, `split_structure_corpus.py`
- `govern_structure_homology.py`, `cache_structure_proteins.py`
- `pretrain_pair_geometry.py`, `evaluate_pair_geometry.py`

These implement P1A/P1B. They establish partner-specific contact/distance
geometry, not affinity energetics.

The research-only R0-C confirmation path is implemented by
`acquire_r0c_structures.py`, `build_r0c_holo_index.py`,
`build_r0c_ligand_bank.py`, `cache_exact_structure_proteins.py` and
`cache_r0b_exact_geometry.py --contract r0c`. It produced a verified fresh
geometry panel, but its downstream exact-pair residual failed admission and is
not a production data path.

## Release-pinned affinity data

- `acquire_source_release.py`, `source_affinity/`
- `build_affinity_energy_corpus.py`, `verify_affinity_corpus.py`
- `audit_affinity_corpus.py`, `census_source_affinity.py`
- `govern_affinity_homology.py`, `build_e0_input_manifest.py`
- `audit_e0_input_feasibility.py`

These implement D0-C/D1 data governance. They do not authorize model training.

## Archived unadmitted interfaces

The F6I component-statistic wrapper was archived after consolidation because
it did not complete external biological admission.

Failed P1C/P1R*/F0R implementations and terminal research branches are retained
under `archive/retired_research_20260811/`; `history.md` records their binding
decisions.

## Stage R cycle (2026-08-15)

- `audit_stage10_claims.py` recomputes eleven binding audit findings against the
  retained Stage 9/10 artifacts and checkpoints.
- `stageR0_retrieval_falsification.py` runs the strict retrieval falsification:
  ten strata, five protein counterfactual arms, and leave-one-protein-component
  -out nested selection so that tuning is separated from inference.
- `build_double_cold_split.py` builds the governed two-axis split
  (`bindingdb_ki_double_cold_v1`). Label-blind; nothing in it reads `pK`. It
  never modifies the corpus, and `QPSMPData(..., split_directory=...)` verifies
  the frozen assignment hash on load.
- `stageR2_representation_discriminator.py` compares frozen ligand and protein
  representations under the double-cold protocol on continuity, low-similarity
  performance and protein-neighbour specificity rather than overlap-heavy MSE.
- `train_level_shape.py` implements counterfactual level-shape gradient-routed
  episodic training: single stage, one backward pass, one optimizer step.
- `stageR3_compare_arms.py` scores every arm on one identical bank with
  identical controls and produces component-level paired bootstraps.
