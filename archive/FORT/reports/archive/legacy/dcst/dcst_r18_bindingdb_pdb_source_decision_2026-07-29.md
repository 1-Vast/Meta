# DCST-R18 BindingDB PDB-linked source decision

Date: 2026-07-29  
Decision: `STOP_BINDINGDB_PDB_SOURCE_INADEQUATE`

## Result

The frozen, outcome-blind R18 gate failed four of five requirements. The
firewall was clean, but the retained native BindingDB source was too narrow
in independent ligand chemistry and exact ChEMBL-train support to justify
reading numeric affinity or downloading PDB coordinates.

| Gate | Frozen requirement | Observed | Result |
| --- | --- | --- | --- |
| G1 source scale | >=2,000 rows, 150 accessions, 1,500 ligands, 1,500 PDB IDs, 300 DOIs | 1,267 rows, 159 accessions, 340 ligands, 6,561 PDB IDs, 163 DOIs | Fail |
| G2 exact TRAIN bridge | >=500 pairs, 100 targets, 80 homology components | 348 pairs, 58 targets, 55 components | Fail |
| G3 target support | >=20% of ChEMBL-train targets | 71/559 = 12.7013% | Fail |
| G4 ligand support | >=20% of the fixed 20,000-ligand sample at Tanimoto >=0.40 | 7.1450% | Fail |
| G5 firewall | No protected-entity or ChEMBL-provenance violation | Passed; 307 rows excluded before support measurement | Pass |

The large number of PDB identifiers does not repair the information problem:
many structures repeat a small ligand set. The median maximum source
Tanimoto for the fixed ChEMBL sample was only `0.252336`; its 90th percentile
was `0.367816`, still below the support threshold.

## Outcome and sealed-data firewall

- Ki/Kd numeric bytes were not decoded, parsed, copied, hashed, logged, or
  written.
- PDB coordinates were not downloaded.
- Confirmation metadata was used only for the preregistered firewall.
- Confirmation features and labels were not loaded or scored.
- The sealed test was not consumed.

The allowed registry contains no row blocked by the protected axes because
all 307 blocked rows were removed before the reported source/support counts.
The exclusion axes included 225 exact accessions, 263 sequence-homologous
rows, 52 exact ligand connectivities, 97 exact scaffolds, and 59 ligand
near-neighbours; axes can overlap.

## Engineering incident

The first formal invocation produced no result after about 73 minutes and was
terminated. Record-boundary detection repeatedly rescanned accumulated
multi-line article records, causing quadratic behavior. It was replaced by
an equivalent incremental quote/tab state machine and a multi-line regression
test. The fixed formal scan completed in `125.753 s`; four focused tests
passed. This changed parsing complexity only, not the safe-field projection,
input, seed, firewall, or frozen gates.

## Consequence

R18 does not authorize a BindingDB Stage-1 dataset. Numeric BindingDB affinity
and PDB-coordinate acquisition remain forbidden for this route. Reopening
requires a materially larger independent source with at least the frozen
ligand, bridge, target-support, DOI, and firewall thresholds; repeated
structures of the same chemistry are not sufficient.

Authoritative machine result:
`reports/active/dcst_r18_bindingdb_pdb_source_seed1729.json`.
