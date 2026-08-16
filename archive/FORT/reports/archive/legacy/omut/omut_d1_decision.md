# OpenMut `OMUT-D1` decision

**Verdict:** `OMUT_D1_TOPOLOGY_ADEQUATE`

**Date:** 2026-07-28.
**Preregistration:** `reports/active/omut_d1_preregistration.md`,
SHA-256 `3c3abe92f0b9e92f95325bd3357df8428665bb54948fcbf4e6c995d132c0ba8f`
(amendment A1: the ChEMBL `type` alias is tolerated present but dropped unread, never
interpreted; a second correction, made before acceptance, scopes documents to the
specific shared-ligand measurements rather than the whole accession).
**Runner:** `research/omut_d1.py`. **Result:** `reports/active/omut_d1.json`.
**Predecessors:** `OMUT-D0` (complete) and `OMUT-F0`
(`OMUT_F0_DAVIS_ROLE_FROZEN__PRESERVE_CONFIRMATION`).

## 1. Gates

All nine gates pass, including `D1_NO_I0_SCOPE_LEAK` (zero keys naming censoring,
replicate covariance, reliability, observed variance, or empirical MDE) and
`D1_PROJECTOR_BLIND` (no forbidden BindingDB affinity column or unrequested ChEMBL
field was ever stored).

## 2. Topology found

BindingDB's local, SHA-256-verified 202607 Articles archive was streamed in full:
**93,712 rows**, 89,158 with a ligand identity, 1,602 excluded as multi-substitution
(sensitivity data, correctly out of scope for the primary graph).

After grouping candidate single-substitution rows by UniProt accession, and requiring
`>=4` ligands shared between the mutant row set and the WT row set of the same
accession:

| threshold `k` | components | distinct accessions | distinct documents (shared-ligand-scoped) |
| ---: | ---: | ---: | ---: |
| 4 | **62** | 17 | 98 |
| 8 | 12 | 8 | 69 |
| 16 | 8 | 5 | 55 |

The document count is deliberately **not** the union of every document touching an
accession — that was the first thing tried, gave 297 at `k=4`, and was corrected
before acceptance because it repeats exactly the over-counting this programme's
promiscuous-ligand-block audit found and rejected
([[document-overlap-binding-constraint]]): a document about some other ligand of the
same accession has no bearing on a specific shared-ligand rectangle. The reported 98
counts only documents attached to the WT or mutant measurement of a ligand that is
actually in the shared set, and is still an upper bound (it does not confirm both
measurements are in the *same* document — that is `OMUT-X0` work).

**Graph.** An accession-accession graph (edge = any shared ligand identity, restricted
to the 17 accessions with a `k=4` component) has **11 connected components** and a
largest component of size 6. The `k=4` incidence matrix (62 components × 360 shared
ligands) has effective rank **25**.

**ChEMBL variant sample** (bounded, 5,000 rows of 119,801 total, explicitly a sample
not a census): 88 distinct targets, 510 distinct mutation strings, 161 distinct
documents after collapse. This reports variant-side topology only; it does not build a
WT-comparison graph, which needs a second query and is `OMUT-X0` work.

**Projected MDE80** at the frozen assumed paired SD `0.10`, `n=62`:
**0.034**, just above the `0.03` material-effect threshold this programme has used
throughout. This is optimistic and projected, not empirical — no real paired
difference was computed, because no affinity value was read.

## 3. What this topology actually looks like — read before treating it as good news

62 crossed the frozen adequacy bar of 25, so the formal verdict is
`OMUT_D1_TOPOLOGY_ADEQUATE`. Two things temper that immediately:

**It is concentrated in two or three disease areas, not broadly multi-family.**
The dominant accession is `P04585` (HIV-1 Gag-Pol polyprotein, protease/RT region) with
6 of the 62 components and up to 33 shared-ligand-scoped documents each — HIV
drug-resistance is one of the most intensively studied mutation-affinity areas in
existence, for reasons unrelated to this programme's question. `P00519` (ABL1, the
BCR-ABL kinase behind imatinib resistance and the T315I gatekeeper mutation) and
`P00533` (EGFR) are the other dense accessions. This is the same shape this programme
has met before in KirHub (kinase-only) and the ChEMBL promiscuous-ligand block: real
signal, narrowly sourced. Family-level diversity (KLIFS/UniProt family mapping) is not
computed at D1 and is needed before "17 accessions" can be read as "17 independent
biological questions."

**The accession-only grouping is a stated simplification that specifically risks the
dominant cluster.** BindingDB's HIV Gag-Pol entries carry different residue-range
windows across papers (`[489-587]`, `[501-599]`, `[588-1127]`, ...), reflecting
different numbering conventions or construct boundaries for the protease/RT domains.
D1 groups by accession only, not by range, so a "shared ligand" match between a WT row
and a mutant row of `P04585` is not yet verified to be the same construct window.
This is exactly the risk the preregistration flagged in advance (section 2.1,
"grouping key") rather than one discovered after the fact, but it bears specifically on
the largest cluster, so it is restated here rather than left in a footnote.

**Every mutation token is a regex candidate, not a verified construct.** The detector
finds bracketed `[range,WT#POS#MUT]` notation in `Target Name`; it is not checked
against a reference sequence.

None of this makes the `ADEQUATE` verdict wrong — the gates and the frozen bar were set
before this run, and 62 clears 25 honestly. It means the verdict authorises exactly what
it says: building the evidence-bound registry, not skipping to a mechanism claim.

## 4. Claim boundary

D1 reports label-free topology only. It is **not** evidence about substitution
semantics, target-specific ligand reordering, or dual-cold transfer. The programme
verdict is unchanged:

> **3** - current data do not identify the substitution-geometry or tau-teacher
> mechanism; new prospective measurement conditions or a newly recovered,
> source-resolved public substrate are required.

## 5. What this unlocks

`OMUT-D1` is complete and adequate. It unlocks **`OMUT-X0`**: build the evidence-bound
BindingDB/ChEMBL registry, restricted to the topology found here plus whatever the
bounded ChEMBL sample and a corresponding non-Davis PLATINUM/MdrDB status contribute.
`OMUT-X0` must, at minimum:

1. verify each of the 62 candidate mutation tokens against a reference sequence for its
   accession (UniProt), resolving the accession-only grouping simplification, with
   particular attention to `P04585`'s multiple residue-range conventions;
2. map accessions to KLIFS/UniProt broad family to test whether the 17-accession,
   11-component-graph topology is truly multi-family or is, like KirHub, effectively
   one or two disease-area clusters wearing several accession numbers;
3. run the full-census ChEMBL variant query (119,801 activities, not the 5,000-row D1
   sample) and attempt the WT-comparison query that D1 explicitly deferred;
4. re-apply document-level and per-ligand provenance collapse at full resolution before
   any component count is used for an `OMUT-I0` power claim.

DAVIS-Complete remains excluded by the `OMUT-F0` policy. PLATINUM remains
`blocked_rights`. Neither contributes to this topology.
