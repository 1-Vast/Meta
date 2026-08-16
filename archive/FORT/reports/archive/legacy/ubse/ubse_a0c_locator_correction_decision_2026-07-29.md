# UBSE-A0C locator correction decision

Decision:
`FREEZE_UBSE_A0C_LOCATOR_V2_KEEP_A0_REMOTE_WAIT`

## Outcome

The metadata-only correction passed all eight frozen A0C gates. The corrected
manifest cleanly separates the BioLiP filename serial from the mmCIF
`auth_seq_id`. It does not unlock coordinate/event extraction: remote
coordinate availability was not executed, so A0 remains `WAIT`.

## Corrected identity

The official BioLiP schema defines annotation column 07 as the ligand serial
used in the native ligand filename, while column 20 is the ligand residue
sequence number corresponding to `_atom_site.auth_seq_id`.

Local evidence reproduces the distinction:

| locator | filename serial | mmCIF `auth_seq_id` | local coordinate evidence |
|---|---:|---:|---|
| `4ci1|B|EF2|B` | 1 | 1429 | `pdb_00004ci1_xyz-enrich.cif` |
| `7bqu|A|EF2|A` | 1 | 501 | `pdb_00007bqu_xyz-enrich.cif.gz` |

Only 39 of the 3,467 corrected rows have equal values for these two fields.
Equality is therefore incidental and cannot be used as a fallback.

## Frozen execution results

| check | result |
|---|---:|
| BioLiP physical rows scanned | 989,058 |
| malformed rows | 0 |
| rows for the 2,833 selected PDB entries | 11,609 |
| unique relevant raw keys | 11,609 |
| old manifest rows | 3,467 |
| unique matches | 3,467 |
| missing matches | 0 |
| ambiguous matches | 0 |
| scalar `mmcif_auth_seq_id` rows | 3,467 |
| duplicate corrected coordinate instances | 0 |
| fit / validation / audit rows | 3,130 / 140 / 197 |
| fit / validation / audit panels | 1,260 / 64 / 88 |
| cross-role PDB overlap | 0 |

The corrected manifest SHA-256 is:

`adc72f142e515c47ea18d20d7af08f6a434a30202a1483dcb362115062a068d5`

The old manifest and old A0 result retained their exact pre-run hashes.

## Firewall

The byte scanner decoded only BioLiP columns
`0, 1, 4, 5, 6, 18, 19, 20`. Affinity columns `13, 14, 15, 16` were skipped
as bytes and never decoded. No coordinate file, development/confirmation
label, or sealed label was read.

## Gate disposition

All A0C gates passed:

- frozen input identity;
- byte-level affinity firewall;
- strict one-to-one raw join;
- complete scalar mmCIF residue identity;
- unique corrected coordinate instance;
- unchanged scale;
- unchanged role closure;
- old-artifact immutability.

For the original A0 decision:

- A0-1 is superseded by the three-input A0C hash contract;
- the old A0-2 and A0-3 passes are withdrawn and replaced by the corrected
  A0C validations;
- A0-4, A0-5, and A0-7 remain valid and were revalidated;
- A0-6 remains unexecuted.

## Artifacts

- corrected manifest:
  `dataset/public/biolip2/processed/ubse_a0c_3d_event_sources_v2.parquet`;
- machine-readable result:
  `reports/active/ubse_a0c_locator_correction.json`;
- implementation:
  `research/ubse_a0c_locator_correction.py`;
- tests:
  `tests/test_ubse_a0c_locator_correction.py`.
