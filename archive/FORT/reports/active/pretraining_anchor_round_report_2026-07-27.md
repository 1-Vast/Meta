# Pretraining-anchor round report before reopening exploration

Date: 2026-07-27

## Executive judgment

The proposed route

`build SPD dataset -> pretrain a selection-aware family-selectivity anchor -> transfer to dual-cold affinity`

is stopped before dataset expansion and neural training. The stop is not caused by GPU memory,
model size or insufficient willingness to train. The proposed supervision fails the cheapest
identifiability test: it does not contribute stable information beyond generic ligand promiscuity
across the six target classes.

This report was written before reopening the next exploration step, as required. A subsequent
exact-accession correction has now also stopped MMP-X on the current sources.

## What was actually tested

SAFSA-G0 used the Novartis SPD systematic panel only, with:

- 144 unmerged, single-effect-gene assay groups;
- missing cells kept missing and censored tested negatives retained;
- DrugCentral/ChEMBL-supported cells excluded from the primary matrix;
- each query gene removed completely before constructing its ligand score;
- own-family evidence compared with whole-family removal, wrong-family evidence, selection coverage,
  and global ligand promiscuity;
- equal family weighting and hierarchical family/gene bootstrap;
- a mandatory sensitivity analysis excluding known DrugCentral MOA cells.

The runner and four unit tests passed. The result is reproducible in
`reports/active/safsa_g0.json`.

## Why the route failed

1. **The apparent family signal is taxonomically heterogeneous.** GPCRs and nuclear receptors are
   positive, while enzyme, ion-channel and transporter families are negative. The family-bootstrap
   lower bound is -0.063.
2. **The incremental information is almost zero.** Own-family AUPRC exceeds the global
   ligand-promiscuity baseline by only +0.0047, with a wide interval crossing zero.
3. **Known on-target annotations are not the main explanation.** Removing DrugCentral MOA cells
   leaves essentially the same failed result.
4. **Selection alone is not the explanation either.** Own-family activity strongly beats
   same-family testing coverage, but that only shows real activity information exists; it does not
   show the information is a transferable family-selectivity anchor.
5. **Scaling cannot resolve a target-definition error.** More parameters or memory could fit the
   heterogeneous labels more closely, but the intended general anchor is not identified.

## What remains valid

- SPD is an admissible systematic source that retains inactive results and is useful for auditing
  measurement selection.
- The broad B0 ligand-potency signal remains reproducible.
- A successful pretraining target must encode a local, conditional interaction effect rather than
  another absolute ligand carrier or an unstable six-class average.
- Full training can move to a larger machine after mechanism proof; local compute is only a
  mechanism-test budget.

## Reopened exploration and MMP-X correction

The first reopened route was **MMP-X**, a local-interaction-primitive proposal:

`small chemical edit x protein family -> directional within-target activity change`.

Its original label-blind F0 passed chemical support with 460 transformations, 2,184
transformation--KLIFS-family units, 97 kinase families and all three sources. Before labels were
read, an exact protein-accession firewall showed that 39 nominal cross-source repeated units collapse
to **8 units / 8 transformations / 5 families** when exact target reuse is forbidden. There is no
KIRHub--Papyrus accession-disjoint edge. Verdict:
`MMPX_F0_ACCESSION_FIREWALL_INSUFFICIENT_STOP`.

Therefore both currently proposed pretraining signals are stopped for distinct reasons:

1. SAFSA: supervision is not generally identifiable beyond ligand promiscuity;
2. MMP-X on the current corpus: the mechanism may be valid, but independent protein-cold repeated
   labels do not exist at adequate scale.

Only after recording these two failures is candidate 2/3 opened. Its data gate must require
systematic negatives, multiple proteins, raw provenance and exact-target-disjoint repeated evidence
before any model is trained.

Firewall status: `confirmation_labels_read=true` is a historical project fact; no current-run
confirmation labels were opened; the confirmation partition remains permanently quarantined;
`sealed_test_consumed=false`.
