# KIRHub ASPIRE candidate-2 aligned-site information-oracle preregistration

Date frozen: 2026-07-26, after candidate 1 was closed and before any ASPIRE outcome was computed.

## Pre-result scope amendment from the user's additional literature synthesis

Before any ASPIRE result was computed, the user supplied a compositional active-site x
ligand-fragment proposal. It is incorporated into the same candidate rather than counted as a new
candidate. Candidate 2 is therefore the **CAPIT/ASPIRE active-site compositional route**, with P0
below as its mandatory information-source gate.

If and only if P0 passes, the next preregistration may replace pooled ESM and whole-molecule Morgan
neighbourhoods with (a) aligned position-level active-site identity/physicochemical features and
(b) count-Morgan/BRICS/pharmacophore fragment features in a low-rank inductive interaction. That
follow-up must use panel-specific relative-effect losses and strict ligand-cold evaluation; it may
not use a Transformer, free target/ligand IDs, raw-label pooling across panels, or a posterior.
P0 failure closes the route before those modules are built because even a shared-ligand oracle
would then show no active-site information beyond taxonomy.

## Candidate ledger and failure-root reconstruction

This document proposes **candidate 2 -- ASPIRE-P0 (Aligned-Site Profile-Reordering Oracle)**.
The reopened-round ledger is now **2/3**. Candidate 1 changed supervision but failed because frozen
pooled ESM did not beat matched random protein or coarse taxonomy. Candidate 2 changes the
information source, not model capacity: it uses the structurally aligned 85-residue KLIFS catalytic
site sequence. The 2026-07-22 public KLIFS snapshot is already hashed and licensed for open
academic/industry use in the local manifest.

ASPIRE is not a support posterior or deployable DTA model. It is an intentionally favorable
information upper bound: when a query target is held out by full-sequence homology component, can
aligned-site similarity select training targets whose measured rankings of the *same held chemical
components* resemble the query target better than a KLIFS-group mean? If the answer is no even
with this transductive shared-panel advantage, an aligned-pocket predictor is not justified.

This differs from:

- WTPAIR, which used pooled ESM and strict ligand-cold prediction;
- SPKOP, which used a protein x ligand support kernel to predict unmeasured held ligands;
- MIF-NK, which learned crystallographic interaction-fingerprint labels. ASPIRE reads no KLIFS
  ligand, structure, interaction, affinity, or conformation label—only the aligned pocket sequence
  and taxonomy.

## Frozen coverage and protocol

- 353/358 eligible KIRHub genes have a valid 85-residue KLIFS pocket; the five invalid empty-pocket
  atypical kinases are excluded before outcomes.
- Eight main KLIFS groups contain at least two strict homology components and at least two distinct
  pockets. Across these groups there are 8,876 same-group, cross-homology target pairs; pocket
  identity has median 0.447.
- Query targets are held out by the existing 324 full-sequence homology components.
- Evaluation is repeated over the existing five strict chemical-component folds. The held ligand
  identities are not used to select target neighbours, but their measurements on training targets
  are deliberately available to form this oracle. Therefore this is not ligand-cold performance
  and cannot confirm DTA generalization.
- Only 5--95% non-saturated cells are scored. A target-fold/ligand-fold profile needs at least five
  jointly observed ligands.
- Similarity is fixed aligned-position identity over 85 residues. The eight most similar
  same-KLIFS-group training targets from distinct held homology components are weighted by squared
  nonnegative identity.
- No outcome-dependent residue selection, pocket weighting, neighbour-count search, or model
  fitting is allowed.

## Frozen arms

1. `global_centroid_oracle`: mean profile of all training targets on the held ligands.
2. `group_centroid_oracle`: mean profile of same-group training targets on the held ligands.
3. `aligned_pocket_oracle`: fixed k=8 aligned-pocket weighted mean.
4. `esm_oracle`: k=8 pooled-ESM cosine weighted mean, restricted to the same group.
5. `pocket_group_shuffle_oracle`: fixed within-group permutation of pocket sequences before k=8
   selection.
6. `random_target_oracle`: fixed eight same-group training targets with equal weights.

All oracle arms use held-ligand labels from training targets. This deliberate advantage must be
prominent in every interpretation.

## Frozen metric and gate

Primary metric is target Spearman ranking, averaged within strict homology component and
bootstrapped over components. Required gain is
`max(+0.030, component-level MDE80 at paired SD 0.10)`.

ASPIRE-P0 passes only if:

1. aligned pocket minus group centroid reaches the required gain with LCB95 > 0;
2. aligned pocket beats within-group pocket shuffle and random-target oracle with LCB95 > 0;
3. aligned pocket is not worse than the ESM oracle (LCB95 for pocket minus ESM >= 0);
4. every fold pair has at least 50 target profiles and 40 independent components.

A pass would establish only that aligned pocket composition carries target-specific profile
information under an easier shared-ligand oracle. It would authorize a separately frozen strict
ligand-cold operator using this information. A failure closes aligned-pocket sequence modeling on
this substrate; it cannot be rescued by BLOSUM tuning, learned residue attention, pocket graphs,
docking, or structure Transformers.

## Scientific basis

KLIFS defines and structurally aligns a consistent 85-residue catalytic site spanning the front
cleft, gate area, and back cleft (Kooistra et al., *Nucleic Acids Research* 2016,
doi:10.1093/nar/gkv1082). Its later review records use of this binding-site alignment as protein
descriptors for kinase profiling models (Kanev et al., *Nucleic Acids Research* 2021,
doi:10.1093/nar/gkaa895). These sources establish biological plausibility and representation
provenance, not predictive validity; the preregistered oracle decides validity here.
