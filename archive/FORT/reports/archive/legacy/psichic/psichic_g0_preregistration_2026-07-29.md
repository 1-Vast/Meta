# PSICHIC-G0 availability and overlap preregistration

Date: 2026-07-29  
Status: frozen before recursive repository inventory

## Role

PSICHIC is not treated as real contact evidence. G0 asks whether its public
release is auditable enough to justify a later pair-specific attribution
test on firewalled PLINDER/KLIFS complexes.

Allowed reads:

- GitHub repository metadata, license, pinned commit, recursive path/size
  inventory, README, and dataset README files;
- linked dataset metadata/pages without downloading data rows;
- no trained weight bytes, dataset row, prediction, affinity, development,
  confirmation, or sealed value.

## Frozen gates

All must pass before any weights or datasets are downloaded:

1. explicit code license, stable commit, and downloadable trained weight files
   with known paths and sizes;
2. complete PSICHIC-XL training membership is reachable as a manifest with
   sequence, ligand SMILES, task label type, and split;
3. every training row exposes its source dataset or stable source record ID,
   allowing PDBbind/ExCAPE/Papyrus/ChEMBL/BindingDB overlap closure;
4. the release distinguishes affinity, binary-interaction, and functional
   labels so a label-blind membership projection is possible;
5. no weight, data row, outcome, protected feature/label, or sealed value is
   loaded.

Failure returns `STOP_PSICHIC_G0_MEMBERSHIP_OR_LINEAGE_INADEQUATE`.

Pass returns `REQUEST_PSICHIC_G1_OVERLAP_DOWNLOAD`; it does not authorize
model inference. A later inference gate must compare true ligand/target
attribution against wrong ligand, within-family wrong target, ligand-only,
target-only, and matched-random fingerprints.

