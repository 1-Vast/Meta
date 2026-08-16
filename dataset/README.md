# Local Data Assets

- `raw/`: downloaded sources; never edited by training.
- `processed/`: governed corpora, protein banks, ligand banks, and research
  panels. Each usable asset requires a manifest.
- `sealed/`: immutable evaluation assets.
- `cache/`: reproducible local caches, never numerical authority.
- `episode.py`, `splitseal.py`: retained dataset contract utilities.

The active Cold Target task uses the four
`processed/meta_fewshot/bindingdb_ki_main_v0*` assets listed in
`docs/PROJECT_FILE_ORGANIZATION.md`. Large data are local and ignored by Git.
