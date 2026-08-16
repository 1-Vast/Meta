# Invalid A2S Artifacts

The following files were generated before the registry-to-feature row
alignment correction and are quarantined from all scientific decisions:

- `a2s_pki_smoke_seed1729.json`
- `a2s_pki_seed1729.json`
- `a2s_pkd_seed1729.json`

Reason: their filtered parquet DataFrame index was incorrectly used as the
feature row id.  `ligand_features.npz` is aligned to the complete registry
global row order.  Use only artifacts whose input record reports
`feature_row_alignment=global_registry_source_row` and
`feature_conn_sha_verified=true`.
