# UBSE-A1 Coupling Identifiability and Confirmation Correction

**Date:** 2026-07-30  
**Status:** binding correction to the A1 design review  
**Decision:** `REVISE_UBSE_A1_BEFORE_PREREGISTRATION`

## Why A1-v1 is not yet executable

P0A has passed as a target-marginal proposal, and PLINDER residue-event
strings have correctly stopped as a replacement for the full teacher.
However, the current A1 unbalanced-OT design does not yet isolate
residue-functional-group coupling.

The v1 exact null is:

\[
\lambda^{(0)}_{igk}=c_k m_{ik}q_{gk}.
\]

An unbalanced transport plan may change both row and column marginals. A plan
of the form

\[
\pi_{igk}=r_{ik}u_{gk}
\]

contains no residue-functional-group coupling, but can still improve residue
AP, ligand-dependent residue rankings, wrong-ligand controls, and the
existing residue-only rectangle by changing pair-conditioned marginals.
A soft two-way purification penalty does not make this shortcut impossible.

Therefore A1-v1 cannot be preregistered or interpreted as a coupling model.

## Required exact nulls

A1-v2 must add a pair-conditioned rank-one null:

\[
\lambda^{(\mathrm{rank1})}_{tligk}
=c_{tlk}r_{tlik}u_{tlgk},
\]

where `r`, `u`, and `c` use the same deployment inputs and parameter budget
available to the full model. The full model must beat this null, not only the
target-only `m x q` null.

The comparison set must include:

- target-marginal `m x q x c`;
- pair-conditioned rank-one `r x u x c`;
- dustbin-only;
- fixed-margin balanced OT;
- fixed-mass partial OT;
- the proposed coupling model.

The coupling stage must either preserve the frozen predicted real row/column
marginals exactly, or use an explicit zero-marginal interaction residual.
Any trainable mass/burden term remains in the null and is frozen before
coupling training.

## Functional-group checkerboard

The v1 rectangle collapses over functional groups and proves at most
ligand-conditioned residue redistribution. A1-v2 needs a within-complex,
within-event typed checkerboard. For one target, ligand, and event type:

\[
Y_{ig}=Y_{jh}=1,\qquad
Y_{ih}=Y_{jg}=0.
\]

Its primary contrast is:

\[
\Omega =
\log\lambda_{ig}+\log\lambda_{jh}
-\log\lambda_{ih}-\log\lambda_{jg}.
\]

Every rank-one row-column marginal model has `Omega = 0`. A positive
directional certificate on this topology is therefore the binding evidence
for residue-functional-group placement.

Before GPU training, the source gate must separately count:

- fit functional-group checkerboards;
- development functional-group checkerboards;
- fresh confirmation checkerboards;
- supported event types;
- distinct target, ligand, PubMed, PDB, homology, and scaffold units.

Thresholds must be frozen after a label-blind power calculation and before
typed-event values are read.

## Explicit dustbin contract

The v1 report mentions a trainable dustbin but does not define its augmented
mass, score, mask, or normalization. A1-v2 must state:

- augmented row and column support;
- real and dustbin marginal masses;
- whether partial transported mass is fixed or predicted by the null;
- dustbin compatibility scores;
- masks for impossible functional-group/event cells;
- real-real normalization used by the loss and metrics;
- the exact parameter setting that returns
  `lambda = c * r outer u` cell by cell.

If zero compatibility does not exactly recover the frozen rank-one null, the
model fails its construction gate before training.

## Fresh confirmation requirement

The 88 G0PB/G1 audit panels have already had binary contact results read and
reported. They may remain development/audit evidence but cannot serve as the
paper-level independent confirmation set for A1.

A1-v2 must freeze a fresh role whose binary, typed, and functional-group
labels have never been read for model choice. The fresh role must be closed
against fit/development on:

- exact target and protein homology;
- exact scaffold and ligand-neighbour similarity;
- PubMed, PDB, assay/source lineage;
- pocket/structure similarity;
- pose/template neighbourhood;
- predicted-monomer model training cutoff, template use, and held-PDB
  membership.

The importance of ligand/pocket/template-aware split closure is consistent
with [CleanSplit](https://doi.org/10.1038/s42256-025-01124-5).

## Teacher reliability wording

Cross-PDB/cross-PubMed agreement from one PLIP implementation is
cross-deposition repeatability, not independent extractor reliability.
A1-R must either:

- add a second extractor or a frozen manual chemistry audit; or
- narrow its claim to cross-deposition repeatability.

Repeated ligand identity must close stereochemistry, formal charge,
tautomer/protomer handling, receptor construct, and unresolved residues.

## Proposal gate correction

`K = min(256, L)` selects every residue for proteins with `L <= 256`, making
recall uninformative. A1-v2 must report:

- the fraction of targets with `L <= 256`;
- random and fit-propensity recall at the identical proposal budget;
- enrichment over those controls;
- all teacher positives outside the proposal as end-to-end false negatives.

P0A top-k success does not waive this gate.

## Mandatory baselines

The equal-budget comparison must include:

- a faithful retrained LINKER-form model;
- monomer-augmented LINKER using the same deployment inputs;
- MONN/DrugBAN-style dense bilinear or cross-attention without OT;
- fixed-margin balanced OT and fixed-mass partial OT;
- pair-conditioned rank-one and dustbin-only nulls;
- binary contact, residue-event, residue-FG, and full typed-tensor ablations;
- a cofold/TankBind-like external structure upper bound with membership and
  destruction controls.

If affinity is later authorized, its baselines must include B0,
target-plus-ligand additive, typed pair-burden only, LINKER representation,
and the coupling-only residual.

Relevant prior scope:

- [LINKER](https://pubs.acs.org/doi/10.1021/acs.jcim.6c00527)
  already covers PLIP typed residue-FG labels, sequence/SMILES prediction,
  and downstream affinity use.
- [MONN](https://doi.org/10.1016/j.cels.2019.08.002) and
  [PLANET](https://doi.org/10.1021/acs.jcim.1c01475) cover contact/interaction
  auxiliary supervision with affinity.
- [TankBind](https://arxiv.org/abs/2206.14382),
  [DynamicBind](https://www.nature.com/articles/s41467-024-45461-2), and
  [NeuralPLexer](https://www.nature.com/articles/s42256-024-00792-z) cover
  monomer/apo protein plus ligand prediction of distances, poses, or complex
  states.

Typed maps, two-stage training, predicted monomers, PLIP distillation, OT,
and ANOVA purification are not individually novel.

## Surviving contribution

The defensible combination, if every corrected gate passes, is:

> a membership-closed typed-event teacher certified across deposition and
> extractor/manual review, distilled to deployment-side predicted monomer
> and 2D ligand inputs through coupling-only partial transport, with pair
> burden and the model's own pair-conditioned row/column marginals as exact
> nulls, a true residue-FG checkerboard certificate, and gain on a fresh
> strict dual-cold confirmation set.

This is a contribution about identifiable coupling plus data closure. It is
not yet an established high-innovation affinity model.

## Binding next decision

A1 remains locked until:

1. A0 remote coordinates and instance extraction pass;
2. a fresh untouched confirmation role exists;
3. the functional-group checkerboard has sufficient power;
4. the rank-one/dustbin construction exactly closes;
5. the revised baseline and membership audits are frozen.

Failure of any pretraining-source condition gives:

`STOP_UBSE_A1V2_COUPLING_OR_CONFIRMATION_NOT_IDENTIFIABLE`

Only after all construction gates pass may an A1-v2 preregistration be
written.
