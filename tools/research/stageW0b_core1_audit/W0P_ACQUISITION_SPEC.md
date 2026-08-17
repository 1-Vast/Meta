# W0-P acquisition specification (blocking external asset)

The Core Task 1 W0/W0-P decision is NO-GO on current local assets. The
blocking asset is a standard positive-control panel. This file is the exact
acquisition contract for the next agent; it is not a research result.

## Required panel

A point-mutant / resistance / ortholog affinity panel with matched ligands.
Preferred sources, each independently usable:

1. EGFR gatekeeper panel: WT, T790M, L858R, T790M/L858R with a common
   inhibitor set (e.g., published biochemical Kd/Ki tables with per-mutant
   values).
2. ABL1 panel: WT, T315I, and other clinical resistance variants with common
   inhibitors (imatinib, dasatinib, nilotinib, ponatinib).
3. ALK / BRAF / KIT gatekeeper or resistance variants with a common panel.
4. Stanford HIVdb or equivalent protease/RT resistance variants with matched
   inhibitor IC50/Ki values.
5. Same-experiment ortholog or point-mutant panels from any licensed public
   supplementary table.

## Minimum acceptance (to be re-preregistered before reading)

- >= **20** mutation pairs (WT vs mutant or mutant A vs mutant B);
- >= **100** matched ligand rows, with >= **5** shared ligands per pair on
  median;
- censored / detection-floor fraction <= **0.25**;
- mutation position known for every pair;
- ligand structures available as isomeric SMILES;
- labels in one comparable unit system per panel (pKi or pKd), never mixed.

## Provenance contract

For every file: source DOI/URL, version/date, license, SHA-256, row/column
semantics, censoring field, original supplementary table name, and a local
parsing script under `tools/research/<new_stage>/`.

## What happens after acquisition

1. Build the panel audit (`pairs`, `rows`, censoring, per-pair sign
   consistency).
2. Preregister the W0-P pass rule from the panel's effective sample size.
3. Run correct-position vs random-position vs BLOSUM-matched unrelated vs
   global pooled vs random-protein controls, leave-one-pair-out, >=3 seeds.
4. Only if W0-P passes, re-census the single-platform affinity panels after
   censored exclusion and decide W1.
