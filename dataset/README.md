# Local Data Assets

- `raw/`: downloaded sources; never edited by training.
- `processed/`: governed corpora, protein banks, ligand banks, and research
  panels. Each usable asset requires a manifest.
- `sealed/`: immutable evaluation assets.
- `cache/`: reproducible local caches, never numerical authority.
- `episode.py`, `splitseal.py`: retained dataset contract utilities.

The active Cold Target task uses the four
`processed/meta_fewshot/bindingdb_ki_main_v0*` assets listed in
`../docs/PROJECT_FILE_ORGANIZATION.md`. Large data are local and ignored by Git.

## Availability and authority

Large raw and processed datasets are not redistributed through Git. A local
asset is usable only when its governed manifest records source identity,
processing contract, split metadata and content hashes. Numerical claims must
name the exact manifest and must never treat `cache/` as authoritative input.

The active BindingDB Ki task requires the governed corpus, protein bank, compact
ligand bank and split manifest under `processed/meta_fewshot/`. The sealed
meta-test assets remain physically isolated until the current task contract
authorizes their use. Source licensing and redistribution restrictions continue
to apply to every downloaded raw asset.

`raw/` is immutable source material; `processed/` is reproducible governed data;
`sealed/` is evaluation-only; and `cache/` may be deleted and rebuilt. Dataset
code belongs here only when it defines data contracts. Executable preparation
and verification workflows belong in `scripts/`.
