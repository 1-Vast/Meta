# UBSE-G0P same-scaffold panel topology preregistration

Date: 2026-07-29  
Status: frozen before formal filtered component and holdout enumeration

## Purpose

UBSE-G0R showed that BioLiP reindexed binding-residue labels have a
cross-publication ligand-connectivity signal, but its wrong-ligand control
used the same scaffold only when one was available. A sequence-plus-ligand
student could therefore pass by encoding target and scaffold identities
additively.

G0P asks a stricter, label-blind question:

> Does BioLiP contain enough independent same-target, same-PubMed,
> same-scaffold multi-ligand panels to train and audit a within-panel centered
> residue-contact student under homology, scaffold, and PubMed closure?

No binding-residue value, affinity field, coordinate, or model output is used
to construct the panels or split.

## Evidence known before freeze

An exploratory multi-agent count, not a confirmatory result, found 3,738 raw
`(exact-sequence target, PubMed, scaffold)` multi-ligand panels spanning 2,598
targets, 1,053 scaffold tokens, and 2,229 PubMed IDs. The largest scaffold
contributed 18.75% of panels. No degree-filtered component, conflict-free
packing, or train/audit residual was calculated before this freeze.

Thresholds below reuse the existing program-wide floors (88 independent
mechanism units, 423 predictive-scale units, and 5% concentration); they are
not tuned to the exploratory marginal count.

## Frozen source and fields

- Source:
  `dataset/public/biolip2/processed/closed_registry.parquet`
- Required SHA-256:
  `7905e4edf88073f564baa4b2d4fb50d496432bc4e15e97cccbfa0766b1b0638d`
- Allowed BioLiP columns:
  `target_key`, `sequence`, `accession`, `pdb_id`, `pubmed`, `conn`,
  `scaffold`
- Allowed ChEMBL projection:
  TRAIN-only `target`, `accession`, `hcluster`, `dual_cold_split`
- Forbidden:
  binding-residue content, all affinity fields/values, assay outcomes,
  coordinates, development/confirmation features or labels, and sealed data.

## Frozen panel definition

1. Require nonempty exact sequence hash, sequence, PDB, PubMed, ligand
   connectivity, and scaffold token.
2. Calculate scaffold degree as the number of exact-sequence targets in the
   complete closed registry. Remove every scaffold with degree greater than
   50 targets. `ACYCLIC`, empty, or synthetic empty-scaffold tokens are also
   removed.
3. Collapse duplicate ligand observations within
   `(target_key, PubMed, scaffold, conn)` to the lexicographically lowest PDB.
4. A panel is `(target_key, PubMed, scaffold)` with at least two distinct
   ligand connectivities and at least two PDB entries.
5. Panel contrast count is `n_ligands choose 2`, capped only for reporting;
   no panel is selected by a contact or affinity result.
6. Cluster the panel target sequences using the existing conservative
   4-mer-containment union-find (`k=4`, containment at least 0.40 against the
   shorter sequence). The implementation and constants are imported from
   `tools/chembl_dualcold_registry.py`.

## Frozen conflict graph and packing

Each panel has three conflict resources:

- target homology component;
- scaffold token;
- PubMed ID.

Two panels conflict when any resource is shared. Connected components are
computed over this transitive conflict graph.

The deterministic conflict-free packing orders panels by:

1. ascending sum of their homology, scaffold, and PubMed panel degrees;
2. ascending maximum of those three degrees;
3. SHA-256 of `seed|target_key|pubmed|scaffold`;
4. lexical keys.

It greedily retains a panel only if none of its three resources has already
been used.

The first 88 packed panels form the audit manifest. All source panels sharing
an audit homology component, scaffold, or PubMed are removed from the
potential training substrate. This makes the audit jointly homology-,
scaffold-, and provenance-cold. The audit is fixed before any contact-label
read.

## Frozen gates

All gates must pass:

1. **P1 filtered scale:** at least 1,000 panels, 500 exact targets, 300
   homology components, 300 scaffolds, 500 PubMed IDs, and 2,000 within-panel
   ligand-pair contrasts.
2. **P2 independent topology:** at least 88 conflict components and the
   largest component contains at most 20% of panels.
3. **P3 resources and packing:** the optimistic resource ceiling
   `min(homology components, scaffolds, PubMeds)` is at least 423 and the
   deterministic conflict-free packing contains at least 88 panels.
4. **P4 concentration:** maximum homology-component, scaffold, and PubMed
   panel shares are each at most 5%.
5. **P5 residual training substrate:** after closing the first 88 audit panels
   over all three resources, at least 800 panels, 300 exact targets, 200
   homology components, 200 scaffolds, 400 PubMed IDs, and 1,600 pair
   contrasts remain.
6. **P6 target-domain support:** at least 20% of filtered panel targets have
   an exact accession on a ChEMBL TRAIN target; this is metadata support only.
7. **P7 firewall:** zero binding-residue value, affinity field/value,
   coordinate, development/confirmation feature/label, or sealed outcome is
   read.

Pass:
`REQUEST_UBSE-G1_CENTERED_CONTACT_STUDENT_PREREGISTRATION`.

Failure:
`STOP_UBSE_SAME_SCAFFOLD_PANEL_TOPOLOGY_INADEQUATE`.

Passing authorizes only a separately frozen, affinity-blind residue-contact
student pilot. That pilot must compare a cross-interaction model with a
same-capacity additive two-tower exact null and must evaluate only centered
within-panel residuals on the frozen 88-panel audit manifest. It does not
authorize coordinates, affinity fitting, or a causal/world-model claim.

## Compute

G0P is CPU-only identity parsing, conservative sequence clustering, graph
components, deterministic packing, and metadata support. CUDA is reserved for
a later admitted student.
