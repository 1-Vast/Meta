# PCIC-A0 provenance-separated cycle-closing acquisition design

Date: 2026-07-29
Status: frozen design only; no provider contact, order, assay, or outcome
acquisition authorized

## Purpose

PCIC-O0-P showed that the public ChEMBL measurement graph has only 37 pKi
and 65 pKd joint homology-scaffold-lineage components, with 99.4690% and
95.5129% of cells in their largest components. PCIC-A0 specifies the minimum
new measurement object that would change this topology rather than
reparameterize it.

The design is a sequence of complete target-by-ligand mini-blocks measured in
new operationally independent provenance lineages. It is not a model result
and does not authorize wet-lab activity.

## A0-R reliability tranche

Before information-directed acquisition, measure a complete randomized
reliability panel:

- 12 exact wild-type human single-protein constructs;
- exactly two targets from each of six broad protein families;
- 16 shared ligands from at least 16 parent connectivities and at least 12
  Murcko scaffolds;
- two operationally independent sites, with separate reagent lots,
  operators, instruments, raw-file systems, and analysis lineages;
- one preregistered exact biochemical endpoint, Ki or Kd;
- complete `12 x 16 x 2 = 384` cells before technical replication;
- all inactive, censored, BQL/AQL, failed, and out-of-range measurements
  retained;
- blinded randomized plate position and processing order.

This tranche estimates cross-site covariance, missingness, censoring, and
mixed-difference variance. Its independent mechanism count is six family
blocks, not 384 cells or the number of rectangles. At paired SD 0.10 its
optimistic MDE80 is approximately 0.114; it cannot confirm a 0.03 gain.

## A1-M mechanism tranche

Only after A0-R freezes the empirical covariance and missingness envelope may
an information-directed mechanism tranche be selected.

One independent block contains:

- two exact constructs;
- four ligand parents from four chemical-neighbour components;
- the complete `2 x 4` target-ligand rectangle;
- the same complete block at two independent sites;
- one endpoint and one frozen assay/context per site.

Each block therefore contains 16 biological-site cells before technical
replication and supplies three within-site target-by-ligand interaction
degrees of freedom. Site repetition estimates reliability; it does not
double biological n.

The primary mechanism floor is 88 joint-independent blocks, requiring at
least 1,408 biological-site cells. Across selected blocks:

- no target or sequence-homology component is reused;
- no ligand parent, Murcko scaffold, or chemical-neighbour component is
  reused;
- no site/campaign/raw-analysis lineage is counted twice as an independent
  provenance unit;
- family share is at most 20%;
- every block has a preassigned nonzero inclusion probability;
- at least 20% of the budget is uniform-random complete blocks;
- a separate future confirmation lineage is never used adaptively.

The 88-block count is only an optimistic mechanism floor. It is not strict
dual-cold predictive validation.

## Candidate coordinates and information criterion

All candidate coordinates are frozen without outcomes and evaluated one at a
time:

- entity-weighted whitened pooled ESM-2 x ligand-feature direct operator;
- composition-only;
- family-only;
- ligand-scaffold-only;
- matched random PSD;
- any future RDIB/PD-MVR coordinate that first passes its own public-data G0.

For a candidate complete block `b`, form the nuisance-orthogonal design
contribution `M_b = X_b^T P_(Z_b)^perp X_b`. Selection compares:

- uniform randomized complete blocks;
- D-optimal pseudo-log-determinant gain in the already identified subspace;
- E-optimal minimum nonzero eigenvalue gain;
- a balanced hybrid maximizing worst-family, worst-scaffold, and
  worst-provenance information.

Ridge or diagonal jitter may be used only for numerical reporting and may not
turn a null direction into information. Primary D-optimality uses the
pseudo-determinant on the existing nonzero subspace; E-optimality reports the
smallest admitted singular value directly.

The information-directed arm must beat the uniform arm on all of:

- admitted rank;
- minimum nonzero eigenvalue;
- information participation ratio by independent block;
- empirical MDE under the A0-R covariance;
- family/scaffold/provenance balance.

If it does not, uniform complete blocks control the next tranche.

## GPU and CPU boundary

Use `D:\anaconda\envs\drug`.

- CUDA float64: coordinate whitening, batched block contributions,
  pseudo-log-determinants, minimum-eigenvalue calculations, randomized
  design simulations, covariance propagation, and power curves.
- CPU only: immutable ID normalization, homology/scaffold/provenance
  union-find, conflict-graph construction, CP-SAT or local-swap packing, and
  manifest serialization.

The combinatorial packing step is discrete and has no useful GPU equivalent.
No trainable tensor computation may silently fall back to CPU.

## Stop rules

Stop acquisition planning before provider contact if:

1. a complete A0-R roster cannot be constructed under all identity and
   chemical firewalls;
2. two operationally independent lineages cannot be contracted;
3. inactive/censored/failed outcomes would not be retained;
4. the endpoint or assay context cannot be held fixed;
5. fewer than 88 conflict-free A1-M candidate blocks exist;
6. the predicted information matrix remains rank deficient for every
   semantic coordinate;
7. candidate selection concentrates more than 20% of information in one
   family, scaffold family, or provenance lineage;
8. the uniform random arm is equivalent to the information-directed arm
   within the frozen Monte Carlo uncertainty.

## Predictive-scale boundary

Strict prediction planning remains approximately 423 independent
multi-family components with at least 40 scaffold-diverse query ligands per
target and a new sealed provenance lineage. PCIC-A0 does not pretend that the
1,408-cell mechanism tranche reaches this scale.

No RFQ, provider contact, compound order, assay scheduling, or outcome access
is part of this artifact. Those actions require explicit user authorization,
funding, reagent availability, and a source-specific amendment.

