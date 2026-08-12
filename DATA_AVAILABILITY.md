# Data availability

This repository tracks source code, tests, frozen theory, compact reports and
dataset manifests. It does not redistribute large third-party releases,
feature banks, tensors, checkpoints, prediction dumps or local caches.

Tracked data is limited to:

- provenance, license, checksum and manifest files;
- compact `RESULT.json` evidence;
- governance warnings such as `DO_NOT_USE_FOR_STRICT_EVALUATION.txt`.

Not redistributed:

- ChEMBL, DAVIS, KIBA, BindingDB and structure payloads;
- BioLiP2/RCSB/mmCIF/CCD downloads;
- ESM and ligand feature banks;
- model checkpoints, `.pt`, `.npz`, `.jsonl.gz` and row-level experiment dumps;
- local credentials, SSH keys and machine-specific remote configs.

Primary status and scientific authorization are in `PROJECT_SUMMARY.md`,
`project_state.json`, `task.md` and `history.md`.
