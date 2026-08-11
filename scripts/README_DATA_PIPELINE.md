# Verified Script Surface

Use the root `main.py` for public orchestration. It exposes data preparation
and verification while keeping failed research runners out of the command
surface.

`scripts/` contains workflows whose data or geometry contracts passed their
registered gates.

## Canonical data and sealing

- `data_contract.py`, `data_protocols.py`, `preprocess_dataset.py`
- `build_ligand_bank.py`, `build_protein_bank.py`, `data.py`
- `seal_compiled_dataset.py`, `governed_sources.py`, `runtime.py`

## Open structure geometry

- `acquire_open_structures.py`, `build_holo_complex_index.py`
- `build_structure_supervision.py`, `split_structure_corpus.py`
- `govern_structure_homology.py`, `cache_structure_proteins.py`
- `pretrain_mechanistic_bridge.py`, `evaluate_mechanism_gate.py`

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
