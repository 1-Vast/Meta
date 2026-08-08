# PKIS typed-mechanism external pilot

Status: exploratory external preregistration, frozen before any model score was
computed. Dataset dimensions and a small number of raw label rows were inspected
for schema feasibility. No PKIS2 prediction, contrast, hyperparameter or Gate
result was inspected before this document was frozen.

This pilot is deliberately **not** `E-AFF-X1`, does not alter the ChEMBL X0-B
authorization state, does not read DAVIS or recipient labels, and cannot admit a
statistic to `model/`. Its purpose is to determine whether a dense, controlled
kinase selectivity matrix supports a transferable, biologically typed pair
statistic worth registering on an independent affinity source.

## Frozen questions

1. Does a dense selectivity panel contain structured target-by-compound signal
   after target and compound additive effects are removed?
2. Can a five-channel physicochemical interaction statistic learned on PKIS1
   predict that residual on PKIS2 for unseen targets and unseen ligand
   scaffolds?
3. Does the correct protein pocket beat both the same additive prediction and a
   capacity-matched deranged-pocket arm?
4. Does the complete zero-shot prediction improve raw, location-sensitive error
   over population-only, ligand-only and protein-only arms?

The first three questions concern interaction. The fourth concerns transferable
location. They are reported separately.

## Sources and immutable roles

- PKIS1 continuous inhibition matrix: source/development only.
- PKIS2 continuous inhibition matrix: external transfer only.
- KLIFS human kinase information: protein-pocket sequence, family and group;
  the fixed pocket has 85 aligned positions.
- PKIS SMILES and ECFP4: ligand structure only.

PKIS2 labels may not select mappings, descriptors, regularization, thresholds,
target strata, ligand strata or derangements. Hyperparameters are chosen using
PKIS1 only.

## Inclusion

- Human wild-type protein kinases with one unambiguous KLIFS 85-position pocket.
- Assay target names that map unambiguously to one HGNC gene.
- Molecules parsed by RDKit and represented by a generic Murcko scaffold.
- Main external ligand stratum: PKIS2 generic Murcko scaffold absent from
  PKIS1. Exact normalized SMILES overlap must also be zero.
- Main external target stratum: PKIS2 HGNC target absent from PKIS1.
- Harder family-cold stratum: the KLIFS family is absent from PKIS1.

Mutants, phosphorylation/autoinhibition states, explicit protein complexes,
separate kinase domains, lipid kinases without a KLIFS pocket and ambiguous
aliases are excluded rather than collapsed onto a wild-type pocket.

## Labels and estimands

Let rows be ligands and columns be targets. For each panel,

\[
R = H_L Y H_P,
\qquad
R_{lt}=Y_{lt}-\bar Y_{l\cdot}-\bar Y_{\cdot t}+\bar Y.
\]

`R` is used only for the interaction estimand. Double-centering is performed
within each panel; PKIS2 label-derived means are evaluation instruments and are
never model inputs. The location estimand is scored on the untouched `Y` scale.

## Biological statistic

Each KLIFS pocket residue is encoded by the six SiteAlign properties: size,
H-bond donor count, H-bond acceptor count, charge, aromaticity and aliphaticity.
Each ligand is encoded by pharmacophore atom counts in four topological shells
relative to its Murcko core: core, one bond, two bonds and three-or-more bonds.

Five fixed interaction channels are constructed:

1. hydrogen-bond complementarity: residue donor x ligand acceptor plus residue
   acceptor x ligand donor;
2. ionic complementarity: positive x negative plus negative x positive;
3. aromatic packing: aromatic x aromatic;
4. hydrophobic packing: aliphatic x hydrophobic;
5. steric accommodation: inverse residue size x ligand heavy-atom occupancy.

The feature tensor has `5 x 4 x 85` coordinates. A source-only ridge coefficient
section reduces it to five signed channel scores. Robust source quantiles and a
`tanh` map bound those scores in `[-1,1]`; rescaling yields five candidate
coordinates in `[0,1]`. These are pair-local and contain no target, assay or
ligand identifier.

## Three-module limit

The proposed architecture has only three novel modules:

1. `TypedPocketStatistic`: the analytic five-channel tensor and five bounded
   channel scores;
2. `CrossedNuisanceProjection`: explicit additive-versus-interaction separation;
3. `LawInterface`: a deterministic map from bounded biological coordinates and
   context to the existing ordered-anchor simplex interface.

The frozen Band, CSMO, anchors and law-valued operator are unchanged. This pilot
tests modules 1-2 and contract-tests module 3; it does not train or alter the
operator.

## Models and source-only selection

- Population: PKIS1 grand mean.
- Ligand: ridge from ECFP4 plus scalar RDKit descriptors to source ligand mean.
- Protein: ridge from the flattened KLIFS SiteAlign pocket to source target mean.
- Additive: population + ligand + protein.
- Correct: additive + typed interaction statistic.
- Deranged: the same additive prediction, but the interaction tensor uses a
  deterministic different PKIS2 pocket from the same evaluation stratum.

Regularization is selected from `{0.01, 0.1, 1, 10, 100, 1000}` using PKIS1
only. Protein models use KLIFS-group leave-one-group-out predictions. Ligand
models use five deterministic Murcko-scaffold folds. The interaction model uses
joint held-out protein-group and ligand-scaffold validation. Ties select the
largest regularization value.

## Frozen metrics and controls

Metrics are target-macro MSE, MAE and Spearman on raw activity, and target-macro
MSE/correlation on the double-centered interaction residual. Confidence
intervals use a deterministic target-cluster bootstrap (`10,000` replicates,
seed `20260808`).

The pilot interaction signal is **observed** only if all hold on the exact
target-cold, scaffold-cold PKIS2 stratum:

1. correct residual MSE is below zero/additive residual MSE with 95% upper bound
   of the paired MSE difference below zero;
2. correct residual MSE is below deranged residual MSE with the same bound;
3. correct residual Pearson correlation has a positive 95% lower bound;
4. the ligand-only raw positive control beats population-only with a positive
   lower bound in MSE reduction.

The raw location signal is reported as observed only if Correct beats Population,
Ligand, Protein, Additive and Deranged in target-macro raw MSE with a positive
95% lower bound in paired MSE reduction. No result is promoted if fewer than 20
external targets or 100 external scaffold-cold ligands remain.

Family-cold is a separately reported stress stratum and cannot rescue a failed
main Gate. Per-group effects, target counts, exact/scaffold overlap, mapping
coverage and all exclusions are mandatory.

## Stop rules

- If the typed statistic fails correct-versus-deranged, do not connect it to
  biological `z`; report representation insufficiency on this proxy.
- If interaction passes but raw location fails, retain it only as evidence for
  Claim B and do not call it an end-to-end DTA solution.
- If raw location and interaction both pass, the next legal step is an
  independent Ki/Kd source preregistration plus a sealed novel-family Gate. It
  is still not production admission.
- Orientation, many-body terms, support attention and RFSA remain frozen.

