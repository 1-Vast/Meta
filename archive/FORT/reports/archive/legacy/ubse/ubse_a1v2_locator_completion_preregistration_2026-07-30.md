# UBSE-A1-v2 locator-completion preregistration

Date: 2026-07-30  
Status: frozen before full A1-R/A1-C raw-identity resolution  
Scope: identity metadata only; no network, coordinate body, contact, event,
affinity, confirmation outcome, or sealed-test access

## 1. Question and decision boundary

Can every frozen A1-R and A1-C metadata instance be mapped uniquely from its
BioLiP filename identity to the distinct mmCIF `auth_seq_id`, while freezing
one official RCSB URL per selected PDB?

A pass freezes locators and the URL request set. It does not pass the
extractor half of `SR-5`, establish coordinate completeness, authorize a
coordinate GET, read an event, validate a functional-group checkerboard,
admit the current P0A, or unlock affinity/Stage-2/confirmation scoring.

## 2. Frozen inputs

| Input | SHA-256 |
|---|---|
| `dataset/public/biolip2/BioLiP.txt.gz` | `c92229bbc8c55c3bd84a9813c3e278ba62f4cfa44e6315cc98d9bf63ed64b6ec` |
| `dataset/public/biolip2/processed/ubse_a1v2_a1r_metadata_manifest.parquet` | `e6eca57f15975340540d2c2c0afd4e2775f7026ca8af8375c9b3ef3e4299fee9` |
| `dataset/public/biolip2/processed/ubse_a1v2_a1c_metadata_manifest.parquet` | `6e9b2e24246db1b0853f6f6714bb2f6f2cc9e9bfd22fbbb0d665e586e67add24` |
| source-role preregistration | `1222a869ce5d9d6db72b0957a11f7478b179b29ddad1863ff3adc6cfaff1a924` |

Frozen implementation:

- `research/ubse_a1v2_locator_completion.py`:
  `58eef30a5058db53452b188e309cb7e0dc81d46bd6d14202033ee028fd18c4ed`;
- `tests/test_ubse_a1v2_locator_completion.py`:
  `48be395b23e81cf6cb4e5d315568aa196212019f47ef4bf1a7b636d746b0df67`;
- synthetic verification before freeze: `5 passed`.

Any implementation change after this freeze requires an amendment before a
new full execution.

## 3. Byte and join firewall

The raw scanner may decode only zero-based BioLiP columns:

```text
0  pdb_id
1  receptor_chain
4  ligand_ccd
5  ligand_chain
6  biolip_filename_serial
18 pubmed
19 mmcif_auth_seq_id
20 sequence
```

It must byte-skip affinity columns `13-16`. It may not decode either
binding-residue field, any contact/event field, coordinate body, PLINDER
interaction field, affinity value, confirmation outcome, or sealed label.

The strict join key is:

```text
target_key
pubmed
pdb_id
receptor_chain
ligand_ccd
ligand_chain
biolip_filename_serial
```

Every frozen row must have exactly one raw match and a nonempty scalar signed
integer `mmcif_auth_seq_id`. Filename serial and `auth_seq_id` remain separate
fields. Both BioLiP file-instance and mmCIF coordinate-instance identities
must be globally unique within the combined locator manifests.

## 4. Frozen topology

| Role | Metadata rows | Targets | Expected unique PDB URLs |
|---|---:|---:|---:|
| A1-R | 459 | 153 | 421 |
| A1-C primary | 512 | 512 | 512 |
| A1-C reserve | 64 | 64 | 64 |
| Total URL request set | - | - | 997 |

All 153 A1-R units must contain exactly three uniquely resolved instances.
A1-C primary and reserves retain their frozen order. The expected result is
zero locator failures; if that expectation fails, this execution stops and
the already registered reserve-only replacement rule must be implemented in
an amendment without reading coordinates or events.

Every selected URL must equal:

```text
https://files.rcsb.org/download/{lowercase_pdb_id}.cif.gz
```

There must be no PDB overlap across `a1r`, `a1c_primary`, and `a1c_reserve`.

## 5. Frozen outputs

- `dataset/public/biolip2/processed/ubse_a1v2_a1r_locator_manifest.parquet`
- `dataset/public/biolip2/processed/ubse_a1v2_a1c_locator_manifest.parquet`
- `reports/active/ubse_a1v2_selected_urls.parquet`
- `reports/active/ubse_a1v2_locator_completion.json`

Outputs are create-once. The runner must refuse to overwrite any existing
output or stale temporary output.

## 6. Gates

1. `L0-1 Frozen inputs`: all three input hashes match.
2. `L0-2 Byte firewall`: exactly the eight allowed identity columns are
   decoded and affinity columns `13-16` are disjoint.
3. `L0-3 Strict unique locator join`: all 1,035 frozen instance rows resolve
   once; no missing, ambiguous, nonscalar, duplicate file-instance, or
   duplicate coordinate-instance row exists.
4. `L0-4 Role scale`: 153 complete three-instance A1-R units, 512 A1-C
   primary targets, and 64 ordered reserve targets remain.
5. `L0-5 URL identity and closure`: exactly 997 unique official URLs and zero
   cross-role PDB overlap.

Pass decision:

```text
FREEZE_A1V2_LOCATORS_AND_SELECTED_URLS_KEEP_EXTRACTOR_LOCKED
```

Failure decision:

```text
STOP_UBSE_A1V2_LOCATOR_COMPLETION_INVALID
```

Regardless of outcome:

```text
SR5_locator_and_extractor_readiness = false
```

The next permissible operation after a pass is a separately frozen HEAD-only
availability audit over the exact selected-URL hash. A coordinate-body GET
remains forbidden.
