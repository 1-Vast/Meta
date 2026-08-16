# OpenMut `OMUT-X0` decision

**Verdict:** `OMUT_X0_EVIDENCE_INADEQUATE_STOP`

**Date:** 2026-07-28.
**Preregistration:** `reports/active/omut_x0_preregistration.md`,
SHA-256 `92e9868a015d7fd997fe3aafbf03bc1f857a2a649bb919b73ccce6b7dcd0c58f`.
The two implementation amendments are transport-only: A1 added release-keyed,
projected ChEMBL page caching after a manually terminated unobservable run; A2
kept the D1 UniProt completeness requirement strict while explicitly excluding
additional inactive ChEMBL accessions after the first cached run stopped before
producing a result.
**Runner:** `research/omut_x0.py`.
**Result:** `reports/active/omut_x0.json`.
**Predecessor:** `OMUT-D1`, verdict `OMUT_D1_TOPOLOGY_ADEQUATE`.

## 1. Gate result

Ten of eleven frozen gates pass. The sole failure is
`X0_PRIMARY_TOPOLOGY_ADEQUATE`.

The execution completed the 119,801-row ChEMBL variant census, all 143 deferred
WT queries, all 2,622 requested assay records, and all 1,381 requested document
records. The no-outcome firewall passes: zero numeric affinity, relation,
censoring, unit, temperature, pH, observed-difference, empirical-variance,
reliability, or empirical-MDE value was materialized. DAVIS-Complete,
PLATINUM outcomes, target-conditioned confirmation, and sealed test data were
not read.

## 2. BindingDB result

The SHA-256-verified BindingDB 202607 Articles archive was streamed in full:
93,712 rows. Of the 62 D1 candidates:

- 48 have at least one uniquely sequence-verified mutant construct;
- 12 fail the canonical reference-residue check;
- one has ambiguous sequence resolution;
- one has no sequence-consistent resolution.

The reconstruction found 92 WT-versus-mutant evidence components. Thirty-seven
are exact `Ki` or `Kd` components with at least four shared ligands. Every
mutant construct and reconstructed WT construct in these components is exact.

None is primary-eligible because the frozen BindingDB Articles schema has no
source-native assay identifier. A shared ligand, endpoint, construct, or paper
does not establish that the WT and mutant values were measured in the same
assay context. Inventing an assay key from document or row proximity would
relax the preregistered evidence requirement after seeing the result.

## 3. ChEMBL result

The full ChEMBL 37 variant census contains 119,801 projected activities and 279
`Ki`/`Kd` variant groups capable of the frozen `k=4` test. Deferred WT matching
recovered 5,792 projected WT activity rows and produced 241 components with at
least one shared ligand; 225 have at least four shared ligands. The relevant
assay and document queries are complete.

Sequence verification is asymmetric:

- 227 of 241 mutant sides have explicit `variant_sequence` records consistent
  with exactly the named single substitution;
- zero of 241 WT sides have an explicit assay-level sequence equal to the
  canonical reference.

This is not a missing-query artifact: every one of the 2,622 requested assay
records was retrieved. ChEMBL exposes explicit variant constructs but the
corresponding non-variant assay records do not expose an exact WT construct
sequence. Canonical target identity alone cannot distinguish full-length WT,
domain constructs, truncations, isoforms, complexes, or additional construct
changes.

Seven of 87 additional ChEMBL accessions are inactive UniProt records without a
sequence and are explicitly unresolved: `B1LRJ1`, `K7XJL6`, `O90781`,
`Q0ZMF1`, `Q9WJQ2`, `Q9WKE8`, and `Q9YQ12`. They are excluded rather than
silently imputed and do not affect completeness of the 17 frozen D1 accessions.

## 4. Scientific decision

The primary evidence topology contains:

```text
exact same-assay Ki/Kd components with >=4 shared ligands = 0
broad families represented                                 = 0
largest-accession share                                    = 0.0
```

The apparent D1 information was real row topology but not a trainable
source-resolved comparison design. The specific bottleneck is missing assay and
WT-construct evidence, not model capacity, protein embeddings, ligand
embeddings, CUDA availability, epochs, rank, seed count, or loss design.

`OMUT-I0`, `OMUT-C0`, `OMUT-M0`, and every downstream affinity-training stage
remain blocked. No neural model or real-outcome training is authorized. The
program remains category 3.

## 5. Reopening conditions

Only new information can reopen the stopped route:

1. a source-resolved reconstruction from primary papers or supplements that
   binds both sides to explicit WT and single-mutant construct sequences and a
   common assay/context, then independently reproduces at least 25 exact
   `Ki`/`Kd` components with at least four shared ligands, at least six broad
   families, and no accession share above 50%; or
2. the registered prospective two-site shared-ligand measurement program,
   beginning with the A0 reliability/variance pilot.

Paper count, database duplication, inferred canonical WT identity, assay-text
similarity, synthetic labels, teacher expansion, larger models, additional
seeds, and GPU training do not satisfy these conditions.
