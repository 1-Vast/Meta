# DCST-R17 KLIFS exact-overlap bridge preregistration

Date: 2026-07-29  
Status: exploratory counts observed; firewall gate frozen before formal audit

## Motivation

PLINDER transfer failed because target and ligand support were sparse, and
DTIOD failed because its active-Morgan intervention was not anchored to a
ligand-specific structural interaction. KLIFS supplies a different
high-quality source object:

- 85 aligned kinase-pocket positions;
- seven observed interaction-fingerprint types;
- exact UniProt and ligand connectivity;
- crystallographic quality and conformation metadata.

A preliminary metadata-only exact join, observed before this
preregistration, found 350 ChEMBL-train target-ligand pairs with at least one
matching KLIFS complex, spanning 66 targets and 60 ChEMBL homology
components. These are planning values, not confirmatory results.

## Firewall

The formal source registry starts from valid human KLIFS mechanisms. Exclude
every source row with any of:

- exact UniProt overlap with ChEMBL development or confirmation;
- maximum 4-mer containment at least `0.40` against an available protected
  85-position KLIFS pocket;
- exact protected ligand connectivity or Murcko scaffold;
- Morgan Tanimoto at least `0.95` to a protected ligand.

Only ChEMBL target, accession, connectivity, scaffold, endpoint, homology
component, and split metadata may be read. No affinity column may be
requested or loaded. Confirmation metadata may define the firewall but may
not fit a model or support statistic.

## Three information objects

1. `KLIFS mechanism source`: all allowed complex-level `85 x 7` IFP labels.
2. `Rectangular mechanism core`: the iterative bipartite 2-core over exact
   `(accession, connectivity)` edges. Both target and ligand degrees must be
   at least two. This is the only subset allowed to claim separation from
   target-only and ligand-only marginals.
3. `Exact affinity bridge`: allowed KLIFS pairs joined to ChEMBL TRAIN on
   exact accession and exact ligand connectivity. ChEMBL provides the
   affinity observation only after this label-blind gate passes.

Repeated structures of one exact pair remain one biological pair. Their IFP
agreement is reported; structures are not counted as independent affinity
observations.

## Frozen admission gate

R17 admits construction of a new Stage-1 dataset only if all hold:

1. at least 1,500 allowed valid complexes, 100 accessions, and 800 scaffolds;
2. rectangular 2-core has at least 500 exact edges, 50 accessions, 100
   recurring ligands, and 40 connected components;
3. exact ChEMBL-train bridge has at least 300 pairs, 60 targets, and 50
   homology components;
4. as corrected before execution in
   `dcst_r17_target_support_correction_2026-07-29.md`, at least 20 ChEMBL
   development targets have an available KLIFS pocket and at least 20% have
   allowed-source aligned KLIFS-85 pocket identity at least `0.25`, while
   the `0.40` 4-mer homology firewall remains intact;
5. at least 20% of a deterministic 20,000-ligand ChEMBL-train sample has
   allowed-source Morgan Tanimoto at least `0.40`;
6. median within-pair IFP Jaccard for repeated structures is at least `0.70`.

Failure returns `STOP_KLIFS_BRIDGE_INADEQUATE`. It does not permit relaxing
the firewall or reviving the old MNI-0 model.

## Conditional model direction

If admitted, preregister `Rectangularly Identified Fingerprint Bridge`:

- Stage 1a predicts a double-residualized IFP on the rectangular core, so
  target-only and ligand-only means are removed before supervision;
- Stage 1b calibrates only the surviving mechanism directions on the exact
  ChEMBL-train affinity bridge;
- Stage 2 trains on the full ChEMBL strict dual-cold train registry with an
  exact-null structural residual;
- scratch, ligand-only IFP, target-only IFP, shuffled IFP, and unresidualized
  MNI-0 controls are mandatory.

This architecture is not implemented until the present gate passes.
