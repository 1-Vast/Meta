# E-AFF-X0-B Crossed Design Re-Registration Result

## Verdict

```text
X0B_CONDITIONAL_DESIGN_SUPPORTED_KI
X0B_CONDITIONAL_DESIGN_SUPPORTED_KD
```

Under the cell-disjoint rectangle unit, the already-governed ChEMBL37 source
does contain a design that reaches the frozen requirement of `245` effective
units — **conditional on intra-cluster correlation being small**. The frozen
effect size, alpha, power, `245` requirement and the `+0.03` affinity Gate
margins are unchanged.

## Label Firewall

X0-B read only the label-blind panel geometry X0 published, after checking that
X0's own firewall record was clean. Zero affinity fields, zero SQLite access,
zero DAVIS or recipient reads, no training.

## The New Model Contains X0 As Its Total-Correlation Corner

With cluster sizes `m[g]`, `N = sum(m[g])` and `m_A = sum(m[g]^2)/N`:

```text
DEFF  = 1 + (m_A - 1) * rho
n_eff = N / DEFF
```

At `rho = 1` the best design collapses to one unit per cluster and `n_eff`
becomes exactly `36` for Ki and `12` for Kd — the published X0 counts. X0 was
therefore not a separate finding but the corner of this model at total
within-component correlation. Every other value of `rho` was unmodelled.

## Design

| Quantity | Ki | Kd |
|---|---:|---:|
| cell-disjoint units | 11,168 | 1,041 |
| target- and ligand-disjoint units (comparator) | 705 | 62 |
| X0 effective components (`rho = 1`) | 36 | 12 |
| clusters `G` | 36 | 12 |
| largest cluster | 5,381 (48.2%) | 417 (40.1%) |
| distinct target pairs | 205 | 49 |
| distinct targets | 224 | 73 |
| distinct ligands | 19,062 | 1,256 |
| units per target pair | 54.5 | 21.2 |

Cell-disjointness was verified programmatically for every packed panel.

## Feasibility Curve

Best per-cluster cap at each `rho`, maximising `n_eff`:

**Ki**

| rho | cap | units | DEFF | n_eff | reaches 245 |
|---:|---:|---:|---:|---:|---|
| 0.00 | none | 11,168 | 1.000 | 11,168.0 | yes |
| 0.01 | 200 | 2,285 | 2.411 | 947.6 | yes |
| 0.02 | 80 | 1,406 | 2.189 | 642.3 | yes |
| 0.05 | 50 | 1,081 | 2.933 | 368.6 | yes |
| 0.10 | 32 | 827 | 3.596 | 230.0 | no |
| 0.20 | 20 | 601 | 4.358 | 137.9 | no |
| 0.50 | 12 | 419 | 6.351 | 66.0 | no |
| 1.00 | 1 | 36 | 1.000 | 36.0 | no |

**Kd**

| rho | cap | units | DEFF | n_eff | reaches 245 |
|---:|---:|---:|---:|---:|---|
| 0.00 | 500 | 1,041 | 1.000 | 1,041.0 | yes |
| 0.01 | 200 | 822 | 2.507 | 327.8 | yes |
| 0.02 | 80 | 470 | 2.132 | 220.5 | no |
| 0.05 | 50 | 380 | 2.971 | 127.9 | no |
| 0.10 | 32 | 290 | 3.741 | 77.5 | no |
| 0.20 | 20 | 201 | 4.429 | 45.4 | no |
| 0.50 | 12 | 135 | 6.230 | 21.7 | no |
| 1.00 | 1 | 12 | 1.000 | 12.0 | no |

## Decision Thresholds For X1

| Endpoint | Breakeven `rho*` | Cap at `rho*` | Hard cluster bound `G/245` |
|---|---:|---:|---:|
| Ki | **0.0915** | 32 | 0.1469 |
| Kd | **0.0164** | 125 | 0.0490 |

`n_eff` is bounded by `G / rho` regardless of packing size, so the cluster count
is the binding resource and no amount of extra rectangles can rescue a large
`rho`. Ki tolerates roughly a 9% intra-cluster correlation; Kd tolerates under
2%, which is demanding.

A separately registered X1 for either endpoint must estimate `rho` with its
uncertainty from exact-assay replicates and fitted variance components, compare
the **upper** confidence bound against the `rho*` above, and abstain without
testing interaction if it is exceeded. That order is preregistered.

## Honest Limitations

1. `rho` is a property of measured values and is not estimable label-blind. Both
   verdicts are conditional design statements, not evidence that
   protein-by-ligand interaction exists.
2. The protein axis is thin. Ki's 11,168 units rest on only **205 distinct
   target pairs** and Kd's 1,041 on **49**; ligand diversity (19,062 and 1,256)
   is what makes the counts large. Any interaction variance the design detects
   would generalise over ligands far better than over proteins, and reused
   target pairs are part of what `rho` must absorb.
3. Packing counts are auditable lower bounds from a greedy target matching, not
   maximum packings.
4. Per-cluster caps are reported as a design option; no data was subsampled.
5. Cell-disjointness removes shared-measurement and shared-interaction
   dependence. It does not remove shared assay-run scale or plate structure,
   which is exactly what `rho` represents.

## If `rho` Turns Out Too Large

The lever is `G`, not more rectangles. E-AFF-X0-FEAS measured the ceiling on
cluster count across governed populations at Ki `97` and Kd `56`, which would
raise the hard bound `G/245` to `0.396` and `0.229`. Realising that requires a
separately registered population extension — lifting the E0-Core `>=20` compound
task contract for census purposes and recanonicalising ligand keys over the
broader corpus. That is out of X0-B's scope and is not authorized here.

## Scope

X0-B re-registers the crossed-census unit and its effective sample size only. It
authorizes a separately registered X1 for Ki and Kd under the `rho*` gate above,
and nothing else. X2, angular and many-body bases, RFSA, DAVIS, production
integration and P2-P4 remain frozen.
