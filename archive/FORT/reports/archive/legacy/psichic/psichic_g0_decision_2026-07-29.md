# PSICHIC-G0 decision

Date: 2026-07-29  
Decision: `STOP_PSICHIC_G0_MEMBERSHIP_OR_LINEAGE_INADEQUATE`

## Result

The public code and weight inventory is usable:

- repository `huankoh/PSICHIC` was pinned at commit
  `cc445aa28d044c6212023705208b1a7704a00622`;
- GitHub detects Apache-2.0 and the root `LICENSE` is present;
- five weight blobs are listed with nonzero sizes, including the
  `23,977,503`-byte multitask `model.pt`.

Those facts pass the release-availability gate, but they do not identify the
data that produced PSICHIC-XL.

The recursive repository tree contains PDBbind-2020 benchmark CSVs but no
`large_scale_interaction_dataset` train, validation, test, or membership
file. The linked public Google Drive root exposes a
`LargeScaleInteractionDataset` folder. Its metadata page currently lists
only:

- `test.csv`;
- `finetuning.ipynb`;
- `degree.pt`.

The README says a `train.csv` was used to train PSICHIC-XL, but that training
file is not visible in either public inventory. Thus the complete model
membership is not reachable under G0.

The documented XL row schema separates protein sequence, ligand SMILES,
regression affinity, binary interaction, and multiclass functional-effect
fields. It does not expose `source_dataset`, `source_record_id`, PDB ID,
BindingDB ID, ChEMBL ID, Papyrus ID, or ExCAPE ID. The documented `key` is
only a unique pair key and is not asserted to be a source-record identifier.
Consequently, even a future copy of `train.csv` could not close source overlap
without an additional row-level mapping.

## Firewall

Only the two pinned README files, the GitHub recursive path/size inventory,
and two Google Drive folder metadata pages were read. Zero dataset rows, zero
dataset file bodies, zero outcome values, and zero model-weight bytes were
read. Development, confirmation, and sealed-test data remained untouched.

## Consequence

Do not download or run PSICHIC weights for pair-specific attribution in this
program. Reopening requires:

1. a stable, checksummed public PSICHIC-XL train/test membership release; and
2. a row-level mapping to source dataset and stable source record identifiers.

Authoritative machine result: `reports/active/psichic_g0.json`.

