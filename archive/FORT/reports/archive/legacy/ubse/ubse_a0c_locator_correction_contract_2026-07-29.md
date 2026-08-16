# UBSE-A0C metadata-only locator correction contract

Status: frozen before generating the corrected manifest

## Defect being corrected

The A0 manifest overloaded BioLiP annotation column 07
(`ligand_serial`) as the mmCIF `_atom_site.auth_seq_id`. The official BioLiP
schema defines column 20, not column 07, as the ligand residue sequence
number. Column 07 is used to form the BioLiP-native ligand filename.

Local coordinate evidence reproduces the distinction:

- `4ci1|B|EF2|B|1`: BioLiP column 07 is `1`, while the mmCIF EF2 residue has
  `auth_seq_id=1429`;
- `7bqu|A|EF2|A|1`: BioLiP column 07 is `1`, while the mmCIF EF2 residue has
  `auth_seq_id=501`.

## Frozen inputs

- `dataset/public/biolip2/BioLiP.txt.gz`
  - SHA-256:
    `c92229bbc8c55c3bd84a9813c3e278ba62f4cfa44e6315cc98d9bf63ed64b6ec`
- `dataset/public/biolip2/processed/ubse_a0_3d_event_sources.parquet`
  - SHA-256:
    `3d4dfcf425b3323d192a7ad50fe6873f6a84c2832a227c68bcdbdc858da8067d`
- `reports/active/ubse_a0_3d_event_source.json`
  - SHA-256:
    `30494dfa0bec15096f359e83e15155d66cff977a76cb72ee19580b53df35db58`

The old manifest and result are read-only frozen inputs and must retain their
hashes.

## Byte-level firewall

The raw gzip scanner may decode only BioLiP columns:

`0, 1, 4, 5, 6, 18, 19, 20`.

Columns `13, 14, 15, 16`, which contain affinity fields, must be skipped as
bytes and never decoded. No coordinate file, development/confirmation label,
or sealed label may be read.

## Strict join and corrected schema

The raw identity key is:

`(target_key, pubmed, pdb_id, receptor_chain, ligand_ccd, ligand_chain, column07)`.

Every old manifest row must match exactly one raw row. The corrected manifest
must replace the overloaded fields with:

- `biolip_filename_serial` from BioLiP column 07;
- `mmcif_auth_seq_id` from BioLiP column 20;
- `biolip_file_instance_id`;
- `mmcif_coordinate_instance_id`.

Missing raw matches, multiple raw matches, a non-scalar auth sequence ID, or
a duplicate corrected coordinate instance is a hard failure. Missing
`mmcif_auth_seq_id` must never be filled from the filename serial.

## Frozen outputs

- `dataset/public/biolip2/processed/ubse_a0c_3d_event_sources_v2.parquet`;
- `reports/active/ubse_a0c_locator_correction.json`.

## Gates

All of the following must pass:

1. all three frozen input hashes match;
2. the byte firewall is intact;
3. all 3,467 rows have a unique raw join;
4. all 3,467 rows have a scalar, nonempty `mmcif_auth_seq_id`;
5. corrected mmCIF coordinate instance IDs have zero duplicates;
6. the A0 scale thresholds remain satisfied;
7. PDB IDs remain disjoint across fit, validation, and audit;
8. both old artifacts retain their hashes.

A0C passing freezes the corrected metadata manifest only. Remote coordinate
availability remains unexecuted, so A0 remains `WAIT`; A1 event extraction is
not unlocked.
