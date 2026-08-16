# UBSE-A0 three-dimensional event source decision

Date: 2026-07-29  
Decision: `WAIT_UBSE_A0_EXTERNAL_COORDINATE_FETCH`

**Superseded locator interpretation:** this decision correctly records the
BioLiP file-instance counts and PDB-role closure, but incorrectly called
column-7 `ligand_serial` an exact mmCIF `auth_seq_id`. The binding correction
is `ubse_a0_ligand_locator_semantics_correction_2026-07-29.md`; A1 remains
locked pending the additive A0C rerun.

## Outcome

The coordinate source is locally addressable at adequate scale and with
complete role closure, but the current execution sandbox denies outbound
socket access. Six of seven frozen A0 gates pass. Remote coordinate
availability remains unmeasured, so A1 event extraction is not yet
authorized.

This is an environment dependency, not a biological or model failure.

## Local source result

| Quantity | Observed |
| --- | ---: |
| coordinate instances | 3,467 |
| frozen panels | 1,412 |
| unique PDB entries | 2,833 |
| complete locator rows | 3,467/3,467 |
| duplicate coordinate locators | 0 |

Role-specific counts:

| Role | Rows | Panels | Unique PDB entries |
| --- | ---: | ---: | ---: |
| fit | 3,130 | 1,260 | 2,496 |
| validation | 140 | 64 | 140 |
| audit | 197 | 88 | 197 |

Fit, validation, and audit have zero pairwise PDB overlap. The original
G0PB audit also remains unique by homology component, scaffold, and PubMed.

## Gate result

| Gate | Result |
| --- | --- |
| A0-1 frozen input identity | pass |
| A0-2 complete locator | pass |
| A0-3 unique coordinate instance | pass |
| A0-4 source scale | pass |
| A0-5 PDB role closure | pass |
| A0-6 remote availability at least 95% | not executed |
| A0-7 firewall | pass |

The manifest was built from exactly nine metadata columns. No binding-residue
field, affinity field or value, development/confirmation data, or sealed
outcome was loaded.

## Coordinate contract

Every row has an official RCSB PDBx/mmCIF address:

```text
https://files.rcsb.org/download/{pdb_id}.cif.gz
```

and an exact `(receptor auth chain, ligand CCD, ligand auth chain, ligand auth
sequence id)` locator. BioLiP-native receptor and ligand filenames are also
recorded.

The official BioLiP download page states that its structures are freely
available and documents the receptor/ligand naming convention. The official
RCSB file-service documentation provides current HTTPS mmCIF download URLs.
The local runtime nevertheless rejected an HTTPS HEAD request at socket
creation, before any remote file could be assessed.

## Consequence

- The immutable source manifest is ready for a network-enabled fetch.
- A future run must use `--remote-check`, record per-role coverage, and pass
  at least 95% overall and within every role.
- Coordinate bytes, PLIP/ProLIF events, and a 3D teacher remain absent.
- `REQUEST_UBSE_A1_EVENT_EXTRACTION_PREREGISTRATION` is not issued.
- Affinity, confirmation, and sealed outcomes remain locked.

## Artifacts

- Contract:
  `reports/active/ubse_a0_3d_event_source_execution_contract_2026-07-29.md`
- Machine result:
  `reports/active/ubse_a0_3d_event_source.json`
- Source manifest:
  `dataset/public/biolip2/processed/ubse_a0_3d_event_sources.parquet`
- Implementation:
  `research/ubse_a0_3d_event_source.py`
- Tests:
  `tests/test_ubse_a0_3d_event_source.py`
