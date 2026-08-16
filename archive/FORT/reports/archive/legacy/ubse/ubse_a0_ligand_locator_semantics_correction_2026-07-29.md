# UBSE-A0 ligand-locator semantics correction

Date: 2026-07-29  
Status: frozen before external coordinate acquisition or A1 event extraction  
Binding effect: withdraw the claimed direct
`ligand_serial == mmCIF auth_seq_id` mapping and require an additive A0C
metadata rerun

## Trigger

The original A0 contract interpreted BioLiP annotation column 7,
`ligand_serial`, as the ligand residue's mmCIF `_atom_site.auth_seq_id`.
That interpretation is false.

Official BioLiP documentation assigns separate meanings:

- column 7 is the ligand serial number used to distinguish repeated
  instances of the same CCD in one chain and to construct the BioLiP ligand
  filename;
- column 20 is the ligand residue sequence number corresponding to mmCIF
  `_atom_site.auth_seq_id`.

The original source manifest therefore remains a valid BioLiP file-instance
manifest, but its `coordinate_instance_id` must not be interpreted as an
exact mmCIF ligand-residue locator.

## Evidence available before correction execution

The immutable raw annotation is
`dataset/public/biolip2/BioLiP.txt.gz`, SHA-256
`c92229bbc8c55c3bd84a9813c3e278ba62f4cfa44e6315cc98d9bf63ed64b6ec`.
An outcome-free scan joined all 3,467 A0 rows uniquely back to the raw
annotation using exact target sequence, PubMed, PDB, receptor chain, ligand
CCD, ligand chain, and column-7 serial. Column 20 was nonambiguous for every
row. Only 39 of 3,467 rows had equal column-7 and column-20 values, directly
falsifying the old equality assumption.

Two local coordinate fixtures provide concrete checks:

| A0 row | BioLiP column 7 | mmCIF `auth_seq_id` / BioLiP column 20 |
| --- | ---: | ---: |
| `4ci1|B|EF2|B|1` | 1 | 1429 |
| `7bqu|A|EF2|A|1` | 1 | 501 |

No affinity value, confirmation label, sealed outcome, or A1 event was read
to discover or verify this mismatch.

## Binding A0C schema

A0C must preserve the original manifest and create a new metadata artifact.
The corrected schema must use:

- `biolip_filename_serial`: BioLiP column 7;
- `mmcif_auth_seq_id`: BioLiP column 20;
- `biolip_file_instance_id`:
  `(pdb_id, receptor_chain, ligand_ccd, ligand_chain,
  biolip_filename_serial)`;
- `mmcif_coordinate_instance_id`:
  `(pdb_id, receptor_chain, ligand_ccd, ligand_chain,
  mmcif_auth_seq_id)`.

Missing column-20 values may never be imputed from column 7. The raw
annotation scanner may decode only the identity fields needed for the join
and must byte-skip all four affinity fields.

## Corrected gates

Before coordinate or event extraction:

1. bind the original A0 manifest hash and the raw BioLiP annotation hash;
2. require an exact, unique raw-row join for all 3,467 rows;
3. require nonempty `mmcif_auth_seq_id` for every selected row and zero
   duplicate corrected coordinate-instance IDs;
4. recompute the original scale and PDB-role-closure gates;
5. retain the original affinity/contact/confirmation/sealed firewall;
6. separately require at least 95% RCSB coordinate availability overall and
   in every role.

Until A0C is executed, the binding status is:
`INVALIDATE_UBSE_A0_AUTH_SEQ_LOCATOR_PENDING_A0C`.

If the corrected local gates pass while outbound coordinate access remains
unavailable, the allowed status is
`WAIT_UBSE_A0C_EXTERNAL_COORDINATE_FETCH`. Neither status authorizes A1.

## Surviving and withdrawn evidence

The original counts, BioLiP filename identities, PDB-entry addresses, scale,
PDB role separation, and firewall evidence survive subject to A0C
recomputation. The following statements are withdrawn:

- `ligand_serial` is an mmCIF `auth_seq_id`;
- the original five-field `coordinate_instance_id` directly locates an
  mmCIF ligand residue;
- A0-2/A0-3 passed at the exact mmCIF-residue level;
- the old manifest was ready for coordinate parsing without correction.

Affinity, Stage-2, confirmation, and sealed access remain locked.
