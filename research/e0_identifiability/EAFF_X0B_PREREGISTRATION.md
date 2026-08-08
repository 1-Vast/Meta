# E-AFF-X0-B Crossed Design Re-Registration

Status: registered before computing any packing, cluster size, design effect or
verdict. X0-B is strictly label-blind. It reads no affinity value, trains
nothing, and does not test interaction.

X0-B replaces the X0 independence unit and its effective sample-size
calculation. It does **not** relax any requirement:

| Frozen quantity | X0 | X0-B |
|---|---|---|
| interaction RMS / assay noise | 0.5 | 0.5 unchanged |
| variance ratio under H1 | 1.25 | 1.25 unchanged |
| alpha (one-sided) | 0.05 | 0.05 unchanged |
| power | 0.80 | 0.80 unchanged |
| required effective units | 245 | 245 unchanged |
| affinity Gate margins | +0.03 / +0.03 | unchanged, out of scope |
| **independence unit** | panel/closure component | **cell-disjoint rectangle** |
| **effective n** | count of components | **design-effect adjusted count** |

## Why The Unit Changes

E-AFF-X0-FEAS established that the X0 unit cannot be produced by crossing. A
rectangle needs two proteins inside one document-keyed panel, and D1 closure
unions every pair of targets sharing a document, so both proteins of every
rectangle always already lie in one closure component. Effective components
therefore collapse onto closure components, whose universe is `245` in total
against a requirement of `245` per endpoint. The best ceiling over every
governed population is Ki `97` and Kd `56`. The X0 stop was a specification
consequence, not a measurement.

Collapsing a component to one unit also asserts that every double difference
inside it is the same random quantity. That is false: distinct rectangles carry
distinct interaction terms. The dependence is in the noise, not in the estimand.

## Model And Estimand

Within a panel `p`, for target `t` and ligand `l`:

```text
y[p,t,l] = mu[p] + alpha[p,t] + beta[p,l] + delta[t,l] + eps[p,t,l]
```

The double difference over a rectangle `(t1,t2,l1,l2)` is

```text
DD = y[t1,l1] - y[t1,l2] - y[t2,l1] + y[t2,l2]
```

`mu`, `alpha` and `beta` cancel exactly, including every panel-level offset and
every effect additive in target or in ligand. Under exchangeable interaction
`Var(delta) = tau^2` and measurement noise `Var(eps) = sigma^2`,

```text
Var(DD) = 4*tau^2 + 4*sigma^2
```

so `H0: tau^2 = 0` is `Var(DD) = 4*sigma^2`, and the design target `tau/sigma =
0.5` is exactly the frozen variance ratio `1.25`. The effect size is unchanged.

## The Unit

Two rectangles have independent `DD` when they share no measured cell:

- `mu`, `alpha`, `beta` cancel regardless of sharing, so reusing a target or a
  ligand does not reintroduce a main effect;
- distinct cells carry distinct `eps`;
- distinct `(t,l)` cells carry distinct `delta` entries, which are independent
  under exchangeability.

The registered unit is therefore the **cell-disjoint rectangle**: within a
panel, no measured `(target, ligand)` cell is used by two packed rectangles.

Packing is constructed to be auditable rather than maximal. Targets in a panel
are greedily matched into disjoint pairs by shared-ligand count; each matched
pair then splits its common ligands into disjoint consecutive pairs. Each target
belongs to one pair and each ligand pair is used once, so cell-disjointness
holds by construction. The count is a lower bound on the maximum packing and is
reported as such.

A strictly target- and ligand-disjoint packing is also computed as a
conservative comparator. It is not the registered unit; it is reported so the
choice is visible.

Because this construction emits several rectangles per matched target pair, the
packing also reports distinct target pairs, distinct targets, distinct ligands
and units per target pair. If interaction magnitude carries target-pair
structure, reused pairs are correlated, so that ratio is part of what `rho`
below must absorb and is required output rather than commentary.

## Dependence And Effective n

Cell-disjointness removes shared-measurement and shared-interaction dependence.
It does not remove shared-assay-run dependence: a panel can compress or stretch
its dynamic range, and plate structure can act on all of its rectangles at once.
That residual dependence is modelled as an intra-cluster correlation `rho`, with
the cluster taken as the **closure component**, which is strictly coarser than
the panel and already unions panels sharing documents or homology.

For clusters of sizes `m[g]`, with `N = sum(m[g])`:

```text
m_A  = sum(m[g]^2) / N
DEFF = 1 + (m_A - 1) * rho
n_eff = N / DEFF
```

`n_eff` is capped by the cluster count: as `m` grows, `n_eff -> G / rho`. The
number of clusters, not the number of rectangles, is therefore the binding
resource, and `rho <= G / 245` is necessary for feasibility at any packing size.

Because oversized clusters inflate `m_A`, X0-B also reports the best per-cluster
cap `c`, maximising `n_eff` over `m[g] -> min(m[g], c)`. Capping is a design
choice available to a later X1; it is reported, not applied to the data here.

`rho` is a property of measured values and cannot be estimated label-blind. X0-B
therefore reports the whole feasibility curve rather than assuming a point:

- `n_eff(rho)` on the preregistered grid `{0, 0.01, 0.02, 0.05, 0.10, 0.20,
  0.50, 1.00}`, at `rho = 1` reproducing the X0 collapse;
- the breakeven `rho*` at which `n_eff = 245`, maximised over the cap grid;
- the hard bound `G / 245`.

## Verdicts

Per endpoint:

- `X0B_DESIGN_INSUFFICIENT_<ENDPOINT>`: the packed unit count `N` is below `245`
  at every cap, so no value of `rho` can reach the requirement.
- `X0B_CONDITIONAL_DESIGN_SUPPORTED_<ENDPOINT>`: `N >= 245` and `rho* > 0`. The
  design meets the frozen requirement if and only if the intra-cluster
  correlation is at most `rho*`.

A conditional verdict authorizes only a separately registered X1 for that
endpoint, and that X1 must, in this order: estimate `rho` with its uncertainty
from exact-assay replicates and the fitted variance components; compare the
**upper** confidence bound of `rho` against the preregistered `rho*`; and
abstain without testing interaction if the bound is exceeded. Estimating `rho`
first and testing second is a preregistered sequence, not a data-dependent
choice of threshold.

X0-B authorizes nothing else. X2, angular or many-body bases, RFSA, DAVIS,
production integration and P2-P4 remain frozen. A conditional verdict is not
evidence that protein-by-ligand interaction exists.

## Population

X0-B runs on the already-governed E0-Core panel geometry published by X0, so no
new acquisition, canonicalization or governance is required. X0-FEAS showed the
cluster count `G` is the binding resource and that broader governed populations
could raise its ceiling to Ki `97` and Kd `56`. Extending the population is a
separate registration and is deliberately out of X0-B's scope; X0-B reports what
the current design supports and what `G` would have to become.

## Required Outputs

- packed unit counts under both packings, per endpoint;
- cluster count, cluster size distribution and `m_A`;
- `n_eff(rho)` grid, best cap per `rho`, breakeven `rho*`, and `G / 245`;
- per-endpoint verdict;
- input and output SHA-256 values, and explicit
  `affinity_value_fields_selected=0`, `davis_label_reads=0`,
  `recipient_label_reads=0`.
