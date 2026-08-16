# KIRHub WTPAIR candidate-1 preregistration

Date frozen: 2026-07-26, before any WTPAIR outcome was computed.

## Candidate ledger

The preceding KIRHub H0, SPKOP A1, strict-firewall correction, KLIFS coverage audit, and
literature/route-overlap audit are prerequisites and consume no slot in this reopened exploration
round, per the user's instruction. This document proposes one new candidate; the autonomous
candidate ledger is therefore **1/3**. No second or third candidate is implied by this document.

**Candidate 1 -- WTPAIR (within-taxonomy protein-pair x ligand-pair ordinal interaction
regression).** WTPAIR changes the supervision rather than the posterior, support mechanism, kernel,
protein encoder, or model capacity. It directly fits the mixed difference

`D(t,u;i,j) = [y(t,i)-y(t,j)] - [y(u,i)-y(u,j)]`

for target pairs `t,u` from the same KLIFS group but different full-sequence homology components,
and ligand pairs `i,j` from the training chemical components. The fixed score class is a small
bilinear ridge on differences of frozen ESM-2 and Morgan features. Same-group pairing removes the
coarse group-centroid shortcut; the double difference removes target intercept and global ligand
potency.

This is substantively different from:

- SPKOP A1, which transported absolute training-target ordinal profiles through a separable
  protein/ligand nearest-neighbour kernel;
- the previously blocked CLBOR idea, whose supervision was ordinary cellwise ordinal residual
  prediction;
- MIF-NK, whose label was a crystallographic interaction fingerprint and failed to establish
  target conditioning. WTPAIR uses direct KIRHub functional rank rearrangements and no KLIFS
  interaction labels or structures.

## Frozen substrate and firewall

- KIRHub Table S4 WT single-concentration matrix; only 5--95% non-saturated cells.
- Target axis: the already hashed 324 full-sequence 4-mer-containment connected components.
- Ligand axis: the already hashed 79 connected components joining exact parent identity, equal
  Bemis--Murcko scaffold, or Morgan radius-2 Tanimoto >= 0.50.
- Outer evaluation: the existing deterministic 5 target-component folds x 5 chemical-component
  folds. All held target components and held chemical components are absent from fitting.
- Training ordinals are recomputed using only the training ligand columns.
- No target-ID lookup, query labels, confirmation labels, or sealed labels enter fitting.
- The single KIRHub source cannot isolate assay/document/publication components; this remains a
  within-source mechanism experiment, not external affinity confirmation.
- KLIFS 85-residue pocket coverage (353/358 eligible genes) is a feasibility observation only and
  is not used by candidate 1.

## Fixed representations and capacity

- Protein: the existing frozen `facebook/esm2_t33_650M_UR50D` pooled vectors. Per outer target fold,
  fit a train-target-only whitened PCA with 16 components, L2-normalize, then subtract the
  train-target KLIFS-group centroid. No protein parameter is trained.
- Ligand: Morgan radius-2 2048-bit vectors. Per outer ligand fold, fit a train-ligand-only whitened
  PCA with 16 components, L2-normalize, and center on the training ligands.
- Interaction score: `s(t,l)=h_t^T W x_l`, 256 coefficients, no bias or main effects.
- Ridge penalty: fixed `alpha=10`, LSQR solver, no outcome-dependent hyperparameter search.
- WTPAIR fit: at most 20,000 deterministic seed-1729 mixed-difference examples per outer fold,
  sampled only from four-observed-cell rectangles. Target pairs must be same-group and
  cross-homology-component.
- Prediction is the train-only ligand-global B0 score plus `s(t,l)`.
- Exactly one seed is authorized. No Transformer, posterior, support adaptation, structural graph,
  second seed, or mechanism revision is authorized by this preregistration.

## Frozen arms

All arms use the identical 25 outer fold pairs and evaluation cells.

1. `ligand_only`: train-target global ordinal ligand profile transferred to held ligands by the
   existing Morgan k=8 interpolation.
2. `group_centroid`: same interpolation applied to the train-target mean profile for the known
   KLIFS group.
3. `cellwise_bilinear`: the same 256-coefficient bilinear ridge fit to ordinary training-cell
   residuals from the ligand-global profile. This isolates the effect of double-difference
   supervision.
4. `wtpair_true`: candidate 1.
5. `wtpair_group_shuffle`: a fixed within-KLIFS-group permutation of frozen protein vectors,
   refitted with the same WTPAIR procedure.
6. `wtpair_random_protein`: fixed matched 16-dimensional random protein features, refitted with
   the same WTPAIR procedure.

## Metric, power, and success gate

Primary metric is target Spearman ranking on held chemical components, first averaged within each
full-sequence homology component and then bootstrapped over components. Targets require at least
five valid held ligands. The prospective effect threshold is
`max(+0.030, component-level MDE80 at paired SD 0.10)`.

Candidate 1 passes only if all conditions hold:

1. `wtpair_true - ligand_only` reaches the threshold and has component-bootstrap LCB95 > 0;
2. `wtpair_true - group_centroid` reaches the threshold and has LCB95 > 0;
3. `wtpair_true - cellwise_bilinear` has LCB95 > 0, proving that changed supervision—not merely
   trainable bilinear capacity—contributes;
4. true minus within-group protein shuffle and true minus random protein both have LCB95 > 0;
5. ordinal RMSE is no worse than `1.02 x` the best non-destroyed baseline;
6. every fold has at least 50 evaluable target profiles and 40 independent homology components.

Failure closes WTPAIR without a larger encoder, extra seed, altered ridge penalty, pocket
augmentation, posterior, or support mechanism. A pass authorizes only review of a separately
preregistered follow-up; it does not authorize confirmation or multi-seed training.

## Scientific basis

KLIFS defines a consistent 85-residue catalytic binding site spanning front cleft, gate area, and
back cleft (Kooistra et al., *Nucleic Acids Research* 2016,
doi:10.1093/nar/gkv1082). Published kinome analyses report that aligned pocket differences,
including but not limited to the gatekeeper, are associated with kinase-inhibitor selectivity
(Müller et al., *Journal of Medicinal Chemistry* 2018, doi:10.1021/acs.jmedchem.8b00699).
These sources motivate the feasibility audit but do not establish WTPAIR's result. The experiment
must pass the frozen destruction and taxonomy controls above.
