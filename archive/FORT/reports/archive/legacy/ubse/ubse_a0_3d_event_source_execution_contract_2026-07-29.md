# UBSE-A0 three-dimensional event source execution contract

Date: 2026-07-29  
Status: frozen before coordinate download, coordinate parsing, or event extraction

**Binding correction:** the `ligand_serial == auth_seq_id` interpretation
below was later falsified before coordinate acquisition. BioLiP column 7 is
a filename ordinal; column 20 is the mmCIF residue sequence number. See
`ubse_a0_ligand_locator_semantics_correction_2026-07-29.md`. The original
contract is retained as the pre-correction record and must not be executed.

## Purpose

UBSE-G1 showed that frozen sequence plus two-dimensional ligand features did
not transport a ligand-conditioned residue-position residual. UBSE-A0
prepares the genuinely new information source required by the successor:
experimental complex coordinates from which residue-by-functional-group-by-
interaction-type events can be extracted.

A0 is only a source-addressability and firewall step. It does not establish
event quality, a deployable student, affinity relevance, or predictive gain.

## Frozen inputs

- Closed BioLiP registry:
  `dataset/public/biolip2/processed/closed_registry.parquet`
- Required source SHA-256:
  `7905e4edf88073f564baa4b2d4fb50d496432bc4e15e97cccbfa0766b1b0638d`
- G0PB panel manifest:
  `dataset/public/biolip2/processed/ubse_g0pb_panels.parquet`
- Required manifest SHA-256:
  `4fea01e332eb3c60e41d76d5062d33cc95b13bc2e96b01df226532f78fe1b371`

Allowed registry columns:

- `target_key`, `pubmed`, `scaffold`, `conn`;
- `pdb_id`, `receptor_chain`;
- `ligand_ccd`, `ligand_chain`, `ligand_serial`.

Forbidden:

- every affinity field or value;
- binding-residue labels during A0 manifest construction;
- development/confirmation features or labels;
- sealed outcomes.

## Frozen split identity

Reuse `research.ubse_g1.build_frozen_split` exactly:

- original G0PB audit: 88 panels;
- deterministic seed-1730 validation: 64 panels;
- closed residual fit: 1,260 pre-contact panels.

No coordinate-bearing PDB identifier may occur in more than one role.

## Coordinate address

For every exact registry row, construct the official RCSB PDBx/mmCIF URL:

```text
https://files.rcsb.org/download/{pdb_id}.cif.gz
```

The row locator inside the entry is:

```text
receptor = auth_asym_id == receptor_chain
ligand   = auth_comp_id == ligand_ccd
           and auth_asym_id == ligand_chain
           and auth_seq_id == ligand_serial
```

The BioLiP-native names are retained as a second locator:

```text
receptor/{pdb_id}{receptor_chain}.pdb
ligand/{pdb_id}_{ligand_ccd}_{ligand_chain}_{ligand_serial}.pdb
```

Official BioLiP documentation defines those names, while RCSB provides the
current mmCIF file service. A0 does not assume that every historical entry is
still current; remote availability must be measured.

## Frozen output

- Row manifest:
  `dataset/public/biolip2/processed/ubse_a0_3d_event_sources.parquet`
- Machine result:
  `reports/active/ubse_a0_3d_event_source.json`

The row manifest contains no label or affinity value.

## Frozen gates

All must pass before coordinate/event extraction:

1. **A0-1 identity:** the two input hashes match.
2. **A0-2 complete locator:** every selected row has nonempty PDB, receptor
   chain, ligand CCD, ligand chain, and ligand serial fields.
3. **A0-3 unique instance:** no duplicate
   `(pdb_id, receptor_chain, ligand_ccd, ligand_chain, ligand_serial)` row.
4. **A0-4 scale:** at least 3,000 complex rows, 2,500 unique PDB entries,
   1,100 fit panels, 50 validation panels, and all 88 audit panels.
5. **A0-5 role closure:** fit, validation, and audit have zero PDB overlap.
6. **A0-6 remote availability:** at least 95% of unique PDB URLs return a
   current coordinate file, with at least 95% coverage separately in each
   role.
7. **A0-7 firewall:** no forbidden column or outcome is loaded.

Pass:
`REQUEST_UBSE_A1_EVENT_EXTRACTION_PREREGISTRATION`.

Failure:
`STOP_UBSE_A0_3D_EVENT_SOURCE_INADEQUATE`.

If A0-1 through A0-5 and A0-7 pass but A0-6 cannot be executed because the
runtime denies outbound downloads, return:
`WAIT_UBSE_A0_EXTERNAL_COORDINATE_FETCH`.

That status is an environment dependency, not a scientific pass or failure.

## Downstream constraints

If coordinates become available, A1 must be frozen before parsing them. A1
must:

- extract residue-by-functional-group-by-event labels with an auditable
  PLIP/ProLIF or equivalent contract;
- use holo coordinates only in the teacher;
- give the student only deployment-available monomer/predicted structure and
  ligand inputs;
- close homology, scaffold, PubMed/PDB source, and model-training membership;
- compare against target-marginal, pair-contact-burden, exact additive, wrong
  ligand, wrong protein, and event-shuffle controls; and
- keep affinity, confirmation, and sealed outcomes locked.
