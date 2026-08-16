# New exploration candidate 1/3: MMP-X F0 preregistration

Date: 2026-07-26  
Candidate count: agent-proposed candidate 1 of at most 3 in the newly reopened round.  
Role: label-blind matched-molecular-transformation support and power audit.

## Why this is substantively new

FACTOR-C/U learned absolute ligand carrier geometry and failed despite 50,000 unlabeled scaffolds.
MMP-X changes the statistical unit and future supervision to a local chemical edit:

`(small substituent A -> small substituent B) × protein family -> within-target rank change`.

It does not optimize an absolute molecule embedding, posterior, support kernel or Transformer.
Published kinase analysis reports that protein characteristics influence activity-cliff propensity
and identifies MMPs whose cliff/non-cliff behavior varies across kinases
(Scientific Reports 2024, DOI `10.1038/s41598-024-59501-w`). Large activity-cliff benchmarks also
show that classical fingerprints can outperform deep models, supporting an audit-first local-edit
baseline rather than immediate architecture growth.

## Frozen MMP definition

Use RDKit single-cut matched molecular pairs. For each cut:

- retained common core has at least 10 heavy atoms;
- exchanged substituent has at most 5 heavy atoms, excluding the dummy attachment;
- both molecules differ at exactly this one attachment;
- transformation orientation is lexical on the two canonical dummy-labelled substituents;
- duplicate molecule-pair/core/transformation records collapse.

This follows the size-restricted MMP convention used in the cited kinase-cliff analysis. No
similarity threshold or after-result edit vocabulary is allowed.

## Label-blind sources

Use only the already frozen structure/provenance projections from KIRHub2026, Reinecke2024 and
Papyrus-Christmann2016. Build paired observation units when both members of an MMP are present for
the same source, endpoint-specific environment and accession. KIRHub's dense design is treated as
one percent-inhibition environment; endpoints are never pooled.

Map accessions to the official raw KLIFS registry. For this F0 power lower bound, collapse targets
by KLIFS family, a conservative coarser unit than exact accession. Unmapped accessions cannot count
toward independent-family gates.

No numerical activity, inhibition, affinity or confirmation label may be read in F0.

## Frozen gates

All must pass:

1. at least 100 unique size-restricted transformations;
2. at least 25 transformations occur in paired observations from at least two sources;
3. at least 100 distinct `(transformation, KLIFS-family, source, environment)` units;
4. at least 30 distinct KLIFS-family units and all three sources contribute paired observations;
5. at least 25 transformations occur in at least two document/environments;
6. the bipartite transformation--family graph has a largest component containing all three sources
   and at least 80% of eligible units;
7. equal-source weighting gives no source more than 0.40;
8. standardized paired-effect MDE80 at SD 0.10 is <=0.05 using independent
   transformation--family units;
9. current-run numerical labels and confirmation labels unread; sealed test unconsumed.

Pass: `MMPX_F0_PASS_AUTHORIZE_DIRECTIONAL_LABEL_AUDIT`.  
Fail: `MMPX_F0_INSUFFICIENT_REPEATED_TRANSFORMATIONS_STOP`.

On pass, MMP-X1 must be separately preregistered before reading numerical development labels and
must compare true protein family against transformation-only, family shuffle and document-private
controls. On failure, do not relax core/substituent sizes or merge transformation identities.
