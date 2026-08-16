# UBSE-S0P same-scaffold panel topology preregistration

Date: 2026-07-29  
Status: frozen before formal panel/component enumeration

## Purpose

UBSE-G0R established that BioLiP binding-residue lists are reproducible across
PubMed/PDB records, but its hard negative used a same-scaffold ligand only
when one was available. An additive target-plus-ligand identity model can
still pass ordinary pair retrieval and wrong-partner destruction without
learning an interaction term.

S0P therefore asks a stricter, affinity-blind prerequisite:

> Can the closed BioLiP source form enough multi-ligand panels that hold exact
> target sequence, PubMed, and Murcko scaffold fixed, while admitting a
> homology x scaffold x PubMed conflict-free fit/audit split?

Only such panels may supervise the subsequent within-panel centered
residue-contact student and its exact additive null.

## Frozen source and columns

- Source:
  `dataset/public/biolip2/processed/closed_registry.parquet`
- Required SHA-256:
  `7905e4edf88073f564baa4b2d4fb50d496432bc4e15e97cccbfa0766b1b0638d`
- Allowed columns:
  `target_key`, `sequence`, `pubmed`, `scaffold`, `conn`, `pdb_id`,
  `binding_residues_reindexed`
- No affinity presence/value, coordinate, development/confirmation feature or
  label, or sealed outcome may be read.

## Frozen panel and closure

1. Reject rows with an empty identity field or an empty/malformed reindexed
   residue list.
2. Collapse duplicate `(target_key, PubMed, scaffold, connectivity)` rows to
   the lexical first PDB.
3. Define a panel as `(target_key, PubMed, scaffold)` with at least two
   ligand connectivities and at least two PDB entries.
4. Cluster panel target sequences with the existing conservative full-sequence
   4-mer containment connected-components rule at `0.40` against the shorter
   sequence.
5. A generic scaffold is one occurring on more than 50 exact target
   sequences among candidate panels. Remove all panels carrying such a
   scaffold and report the removed scaffold hash and rows.
6. Form a conflict hypergraph in which panels are joined when they share a
   target homology component, scaffold, or PubMed. Connected components are
   the independent uncertainty and split units.
7. Assign whole components to five deterministic balanced folds by descending
   component size, then component name; each next component enters the
   currently smallest fold, ties by fold index. Fold 0 is the audit fold and
   folds 1-4 are fit-only.
8. Separately compute a deterministic conflict-free packing. Sort panels by
   SHA-256 of `seed=1729|panel_id`; accept a panel only if its homology
   component, scaffold, and PubMed have not appeared in an accepted panel.

This gate deliberately uses non-isomeric connectivity only to enumerate
panels. A pass does not resolve stereo, charge, tautomer, construct, assembly,
or ligand atom mapping; those remain required before any coordinate-event
teacher.

## Frozen gates

All must pass:

1. **P1 raw panel scale:** at least 2,000 panels, 1,000 exact targets, 300
   scaffolds, and 1,000 PubMed IDs.
2. **P2 post-generic scale:** at least 1,500 panels, 800 exact targets, 300
   scaffolds, and 800 PubMed IDs remain.
3. **P3 independent topology:** at least 88 conflict components and the
   largest contains at most 20% of retained panels.
4. **P4 split support:** fit has at least 1,000 panels/600 targets and audit
   has at least 200 panels/100 targets/80 scaffolds/100 PubMed IDs.
5. **P5 packing:** the optimistic resource ceiling and deterministic feasible
   packing are at least 423 and 88, respectively.
6. **P6 concentration:** largest homology, scaffold, and PubMed shares are
   each at most 5%.
7. **P7 firewall:** zero forbidden data access.

Pass:
`REQUEST_UBSE_S1_CENTERED_ADDITIVE_NULL_STUDENT_PREREGISTRATION`.

Failure:
`STOP_UBSE_S0P_SAME_SCAFFOLD_INDEPENDENT_TOPOLOGY_INADEQUATE`.

A pass authorizes only frozen protein/ligand feature construction and a
separate GPU student preregistration. It does not authorize coordinate
downloads, affinity access, or Stage-2 training.

## Compute

Parsing, RDKit-free panel enumeration, 4-mer homology clustering, union-find,
and packing are CPU-only. CUDA is reserved for the downstream student if this
gate passes.
