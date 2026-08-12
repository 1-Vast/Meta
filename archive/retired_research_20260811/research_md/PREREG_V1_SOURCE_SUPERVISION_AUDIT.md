# MetaSieve v1 source-supervision audit preregistration

Date frozen: 2026-08-11

## Question

Before adding a partner-discrimination or within-target auxiliary loss, test
whether the existing exact-Ki source data contain the measured comparisons
needed to define that loss without treating an unmeasured pair as a non-binder.

This is a census, not model training. It may read exact Ki values from
`meta_train` and `meta_val`; it must not admit a `meta_test` source row. The
consumed main-v0 test targets remain unavailable to v1 development.

## Eligible supervision

Two supervision units are counted separately for `meta_train` and `meta_val`:

1. **within-panel ligand groups**: a `(panel_id, target_id)` with at least two
   distinct measured ligands after within-panel replicate median;
2. **measured partner groups**: a `(panel_id, ligand_id)` measured against at
   least two distinct targets in different frozen CD-HIT40 groups.

Only exact positive uncensored Ki rows already admitted to the main-v0 corpus
are eligible. A panel is the frozen BindingDB `document + Ki + normalized
protocol` assay proxy. No random, graph-distant or absent database edge is
called a negative or assigned an affinity.

## Report

Persist row counts, eligible group counts, target/ligand/family coverage,
replicate counts and the distribution of absolute measured pKi differences.
The report must include input hashes and an explicit count of forbidden
meta-test source-row intersections, which must be zero.

## Decision

- A within-panel auxiliary loss is implementable only if source contains
  repeated eligible panels across multiple targets and CD-HIT40 groups.
- A measured partner loss is implementable only if source contains repeated
  eligible `(panel, ligand)` groups across multiple protein families.
- Coverage is a feasibility result only. It does not authorize v1 training or
  establish biological specificity.

