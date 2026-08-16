# SIMA-DTA FSA-D0 Preregistration Draft

**Stage:** FSA-D0 only.
**Permission:** outcome-safe topology and power audit.
**Prohibited:** model training; new numeric affinity reads; panel_davis
target-conditioned confirmation; ChEMBL confirmation; sealed outcomes.

## Aim

Determine whether a strict k=5 support-query scaffold-cold, target/homology-cold
few-shot DTA evaluation has enough independent, role-closed support for an
architecture test. This is a feasibility audit, not development scoring.

## Inputs

Allowed inputs are frozen split/role metadata, target/homology labels, scaffold
and chemical-neighbour assignments, endpoint/assay/document/provenance metadata,
cached non-outcome protein and ligand features, and TRAIN-only cross-fitted B0
predictions. pKd and pKi are separate strata. B0 prediction artifacts must
already be role-closed; FSA-D0 does not fit or recalibrate B0.

## Per-Target Outputs

For every candidate meta-test target and endpoint stratum, record:

1. support availability at k=1,3,5,10;
2. query depth after each support removal;
3. target and homology closure;
4. support-query scaffold closure;
5. query-to-meta-training scaffold and chemical-neighbour closure;
6. endpoint, assay, document, and provenance structure;
7. independent target/homology component and family concentration;
8. label-free support-feature or Jacobian-rank proxy;
9. TRAIN-only MDE and predicted component-paired arm power.

## Frozen Procedure

1. Freeze candidate target roster and all role labels before any support labels
   are opened.
2. Construct only metadata-derived support/query candidates. Enforce support-
   query scaffold disjointness and transitive query-to-training closure.
3. Compute counts, depths, component graph, family shares, and rank proxy.
4. Estimate MDE and paired-arm power from TRAIN-only cross-fitted episode
   residuals. Do not score a meta-test outcome.
5. Freeze the resulting adequacy thresholds, effect direction, support selection
   tie-breaks, metric order, and inference unit before FSA-B0.

## PASS

FSA-D0 passes only when strict k=5 retains prespecified adequate values for:
target count, independent target/homology components, support scaffolds, query
scaffolds, query depth, adaptation-rank proxy, and expected component-paired
power. The numerical floors are intentionally unset in this draft and must be
computed and frozen from TRAIN-only evidence before a PASS claim.

## STOP

Stop before M1 implementation if any strict k=5 adequacy floor fails, metadata
roles are incomplete, endpoint strata cannot be separated, closure is violated,
the rank proxy is degenerate, or power is inadequate. Do not rescue a STOP with
larger models, altered splits, lower k, larger q, more epochs, seeds, or
support-label access.

## Cheapest Decisive Computation

The first computation is CPU-only: construct the frozen candidate roster,
transitive target/homology and scaffold/chemical-neighbour closure graph, then
tabulate k=5 support and query scaffold depth per endpoint. It reads no new
numeric affinity. It stops the program immediately if strict k=5 has
insufficient independent components or query scaffolds, before Jacobians or any
GPU work.

## Planned Inference

Use target/homology-component paired contrasts and component bootstrap. Rows,
ligands, supports, seeds, and folds are descriptive, not biological n. Freeze
RMSE, MAE, within-target Spearman, ligand-reordering accuracy, concordance,
full-minus-calibration, correct-minus-wrong, and query-span-minus-random
thresholds from TRAIN-only episodes before development access.

