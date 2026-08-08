# Verified Script Surface

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

## Release-pinned affinity data

- `acquire_source_release.py`, `source_affinity/`
- `build_affinity_energy_corpus.py`, `verify_affinity_corpus.py`
- `audit_affinity_corpus.py`, `census_source_affinity.py`
- `govern_affinity_homology.py`, `build_e0_input_manifest.py`
- `audit_e0_input_feasibility.py`

These implement D0-C/D1 data governance. They do not authorize model training.

## Component-statistic interface

- `evaluate_component_statistic.py` applies the internal gauge-separated
  component algebra to precomputed biological surfaces and support residuals.
- It rejects query-label fields and does not construct, train or validate the
  biological surface.

Failed P1C/P1R*/F0R implementations were removed after their evidence was
consolidated in `history.md` and
`EVIDENCE_CONSOLIDATION_AND_FAILURE_TRIAGE.md`. Terminal synthetic, structural
and source-affinity research implementations were removed from the active tree
and remain recoverable from Git history.
