# S1 — Licence and Provenance Audit

Every source is assigned exactly one scientific role. A source may not be used
outside its role, and no Tier-C source may support an affinity claim.

## Structural sources

| source | release | licence | redistributable | role | exposure |
|---|---|---|---|---|---|
| RCSB PDB coordinates, `pilot20k` | acquired 2026-08-06, 14,169 mmCIF, 14,906 holo records | **CC0-1.0** | yes | `STRUCTURAL_TRAINING_ONLY` | **P1B trained on this.** All 10,468 PDB ids treated as EXPOSED, including P1B's val and test partitions, because those were consumed as P1B's own evaluation |
| RCSB PDB coordinates, independent set | RCSB Search API, X-ray, `<= 2.5 A`, `>= 1` non-polymer entity, exactly 1 protein entity, released `>= 2024-01-01`; 15,003 candidates, **exposed overlap 0** | **CC0-1.0** | yes | `STRUCTURAL_INDEPENDENT_TEST` | none by construction |
| RCSB Chemical Component Dictionary | per-component `.cif` from `files.rcsb.org/ligands/` | **CC0-1.0** | yes | ligand structure resolution | none |
| BioLiP2 annotations | used historically to select biologically relevant ligands in `pilot20k` | academic use per upstream readme; **not redistributed** | no | `ANNOTATION_ONLY` | the S2 teacher never reads BioLiP; it is computed from raw coordinates, so BioLiP's absence-as-negative convention cannot enter |

**PLINDER — not used.** Its annotations and derived files carry a separate
licence from its GPL-2.0 software, and the programme's own policy forbids pulling
the full 400k systems. The independent test set is built from raw RCSB CC0
coordinates and locally computed annotations instead, which needs no third-party
annotation licence at all. Recorded as a deliberate substitution, not an
oversight.

**PDBbind — not used.** Redistribution terms were not verified as compatible, and
policy forbids silent use of PDBbind-derived artifacts. No PDBbind-derived
benchmark (including any repackaged split) is present in this programme.

**PoseBusters, CSAR — not used in S1–S4.** They are pose-quality and
structure-affinity resources relevant only to the conditional S7 pose-aware
branch, which S4 has not authorised.

## Affinity sources (not opened in S1–S4)

| source | release | licence | role | opened? |
|---|---|---|---|---|
| ChEMBL 37 static release | tracked release manifest | CC BY-SA 3.0 | `AFFINITY_SOURCE_CALIBRATION_ONLY`, S8 only | **no** |
| BindingDB curated articles | `202608` | CC BY 3.0 | `AFFINITY_REPLICATION` | read in XP4 only; not in S1–S4 |
| Metz 2011 | publisher supplement, SHA-256 pinned | publisher supplementary data | `CONSUMED_DEVELOPMENT_DIAGNOSTIC_ONLY` | consumed by XP1/XP2/XP5 |
| Klaeger 2017 | publisher supplement, SHA-256 pinned | publisher supplementary data | `BOUNDED_SECONDARY_EVIDENCE` | consumed by XP1/XP2-F |
| NIMH PDSP Ki | full CSV export, SHA-256 pinned | free public NIMH resource | `BOUNDED_SECONDARY_EVIDENCE` | consumed by XP1-C |
| DAVIS, recipient labels | — | — | **PROHIBITED** | **0 reads, all stages** |

## Tier C — auxiliary only, never affinity evidence

Yamanishi, BIOSNAP, STITCH, DrugBank, SuperTarget, TTD, DGIdb, DrugCentral,
Open Targets, CTD, ChemProt, Hetionet, DRKG, PharmKG.

**None were downloaded or used.** They are reachable only at S9, which requires a
structurally admitted channel to exist first. Absence of an edge in any of these
resources is never treated as a biological non-interaction.

## Label-read counters for this programme

```text
DAVIS                 0
recipient             0
ChEMBL37 affinity     0
PKIS2                 0
Anastassiadis         0
any affinity value    0   (S1-S4 are entirely label-free)
```

## Attribution

RCSB PDB coordinates and the Chemical Component Dictionary are released under
CC0-1.0 by the RCSB Protein Data Bank. BindingDB is distributed under CC BY 3.0.
ChEMBL is distributed under CC BY-SA 3.0. The NIMH PDSP Ki database is a public
resource of the NIMH Psychoactive Drug Screening Program. KLIFS is open academic
access. Publisher supplementary data (Metz 2011, Klaeger 2017) are used under
their respective journal terms and are **not redistributed** through this
repository; only checksums, manifests and derived statistics are tracked.
