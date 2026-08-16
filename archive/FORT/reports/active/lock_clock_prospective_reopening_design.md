# Prospective reopening design for substitution geometry

## Purpose

The next experiment is a reliability and mechanism-identification pilot, not a performance benchmark.
Its purpose is to create target-by-ligand factorial contrasts that are unavailable in the current
observational graph and to estimate family-level variance before any predictive model is authorized.

## Minimum A0 panel

- 12 targets, exactly two from each of at least 6 protein families;
- at least 16 shared ligands measured against every target;
- scaffold diversity fixed before outcomes, with parent, Bemis-Murcko scaffold, and chemical-neighbour
  identifiers;
- two operationally independent sites using separate reagent lots, operators, instruments, raw-data
  systems, and analysis lineage;
- one endpoint chosen in advance, pKi or pKd, never a mixture;
- complete `12 x 16 x 2 = 384` target-ligand-site cells before technical replication;
- inactive, below-quantification, above-quantification, and failed curves retained with censoring state;
- full-panel inclusion probability one, randomized blinded plate position and processing order.

The 16 ligands yield 120 ligand-pair contrasts per target, but these are not independent power units.
Targets are clustered by the six families, and site replication does not create new biological units.

## Firewalls

Before measurement, freeze target accession, full-sequence homology component, pocket/family role,
binding-profile sentinel, ligand parent, scaffold, chemical-neighbour component, assay protocol,
document/raw-file lineage, site, operator, plate, and provenance block. Query-ligand selection cannot
use expected activity. Any incomplete design must retain a known nonzero inclusion probability and use
the registered inverse-probability analysis.

## A0 estimands

For targets `t,t'` and ligands `d,d'`, the basic reordering contrast is:

```text
Delta(t,t';d,d') = y(t,d) - y(t,d') - y(t',d) + y(t',d')
```

Primary A0 checks are cross-site repeatability of these contrasts, retained-inactive agreement,
family-clustered coordinate-control contrasts, and assay-noise estimates. Each candidate coordinate
enters the same fixed estimator separately. Required controls remain aligned identity, composition,
pooled ESM-2, position/sequence shuffles, BLOSUM permutation, random PSD, matched wrong target, and
within-family wrong target.

True CLOCK is eligible only if positional structure embeddings are frozen from an independent source,
the structure map is trained leave-family-out from the A0 targets, and a parameter-matched
structure-shuffled control is included.

## Power boundary

With paired family/component SD `0.10`, a two-sided 5% test with 80% power has:

```text
MDE80(n) = 2.8016 * 0.10 / sqrt(n)
```

Six independent families give MDE about `0.114`; therefore the 12-target A0 panel cannot establish a
`0.03` prediction gain. At least 88 independent components are needed at SD `0.10`, before allowing
for provenance clustering, multiplicity, missingness, or model fitting. A0 can only estimate variance,
verify two-site reliability, and decide whether expansion to roughly 70-155 mechanism components is
warranted. Strict predictive validation remains blocked until the broader approximately
423-component, scaffold-diverse, provenance-independent substrate is available.

## Reopening rule

Reopen a minimal model stage only if both sites independently show reliable target-specific mixed
differences, the fixed coordinate beats every destruction/capacity control under family-clustered
inference, the measured MDE is adequate for the registered effect, and no assay/document/provenance
sentinel explains the gain. Otherwise stop the coordinate permanently.

