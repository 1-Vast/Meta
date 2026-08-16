# PBCNet2.0-D0 remote source preregistration

Date: 2026-07-29  
Status: frozen before remote CSV-header read

## Scope

PBCNet2.0 now has an open Zenodo record (`18299525`) linked by the 2026 Nature
Chemical Biology paper. D0 determines whether the public package exposes
enough identity and provenance metadata to justify downloading its pair
manifest. D0 does not download the 5.44 GB structure archive or read any CSV
data row.

Allowed remote reads:

- Zenodo record metadata, file names, sizes, checksums, access state, and
  declared license;
- GitHub repository metadata, root file names, commit pin, and declared
  license files;
- only the first newline-terminated header of
  `training_data_8_6.csv`.

The header reader must stop at the first newline. It may not decode, log,
hash, or write any subsequent row.

## Frozen gates

All must pass before a D1 manifest download:

1. Zenodo record is open and declares an explicit data license; the pair
   manifest and source-structure archive each have a stable size and checksum;
2. the header provides reconstructible identifiers or paths for one protein
   complex and two ligand/complex members;
3. any outcome/difference column is separable from identity columns so D1 can
   project a label-blind membership table;
4. the header exposes a BindingDB row/source identifier, document/assay
   locator, or another stable lineage key sufficient for overlap and
   provenance closure;
5. no CSV data row, outcome, ChEMBL development/confirmation label, or sealed
   test is read.

Repository code-license presence is reported separately. A README claim
without a repository license file does not authorize code reuse, but does not
by itself invalidate a data-only audit whose Zenodo license is explicit.

Failure returns `STOP_PBCNET2_D0_IDENTITY_OR_PROVENANCE_INADEQUATE`.

Pass returns `REQUEST_PBCNET2_D1_MANIFEST_DOWNLOAD`; it does not authorize the
5.44 GB structure archive or model weights.

