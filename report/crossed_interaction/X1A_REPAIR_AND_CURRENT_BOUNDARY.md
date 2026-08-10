# X1A repair and current boundary

> **SUPERSEDED CURRENT ENTRY.** The repair described here has now executed and
> returned `X1A_R_DEPENDENCE_PRECONDITION_FAILED`. Use
> `X1A_R_FINAL_SYNTHESIS.md` for current status. The text below preserves the
> pre-execution rationale.

## Current verdict

```text
X1A_ICC_PRECONDITION_NOT_ESTABLISHED
X1A-R_DIRECT_DD_DEPENDENCE_REGISTERED_NOT_EXECUTED
X1B_NOT_AUTHORIZED
X2_NOT_AUTHORIZED_NOT_TRAINED
```

The historical amended X1A artifact remains intact. Its registered Gates
reported PASS, but independent review found that the resulting ICC cannot
support the planned X1B power claim.

## Why the amended ICC is not usable

All 310 targets in the label-blind census are confined to one dependency
cluster. A global target-dummy fit therefore absorbs the cluster indicator.
Ligands observed in only one cell are also fitted exactly. The remaining
nonzero cluster component is produced after changing weights through successive
cell, panel and cluster means; it is not the intended random cluster effect.

More importantly, X1B would test `q = DD^2 - v_noise`, whereas X1A measured
dependence of signed measurement residuals. Symmetric positive and negative
interactions can have near-zero mean residual dependence while `DD^2` remains
strongly clustered. The final Gate artifact also used 2,000 rather than the
registered 10,000 bootstrap draws.

Consequently, the old PASS is retained as historical development evidence but
its authorization effect is withdrawn.

## Completed repair

The original X0-B greedy packing was restored without reading affinity values.
The materialized manifest reproduces the published design exactly:

| Endpoint | Rectangles | Clusters | Frozen cap | Selected |
|---|---:|---:|---:|---:|
| Ki | 11,168 | 36 | 32 | 827 |
| Kd | 1,041 | 12 | 125 | 605 |

`rectangles.jsonl` SHA-256:

```text
22f3e738f4dbc7b53ca9ef23e995e2a398cbca280a9cdde12c546be21500d0a5
```

Each row fixes its four measurement cells, signs, endpoint, panel, dependency
cluster and label-blind cap membership. This closes selection ambiguity before
any direct DD value is computed.

## Registered next stage

X1A-R will calculate

```text
DD = y(P1,La) - y(P1,Lb) - y(P2,La) + y(P2,Lb)
q  = DD^2 - v_noise
```

without fitting target or ligand effects. The same conservative dependence UCB
must control both the rho Gate and effective sample size. Ki and Kd remain
separate; Kd also requires a small-G sensitivity because it has 12 clusters.

For consistency with X0-B's cell-interaction scale, X1B must use `D=DD/2` and
`v_D=sum(v_cell)/4`, or divide the DD-scale second moment by four. Otherwise
the biological RMS margin is off by a factor of two.

No active stage may train a model. GPU optimization becomes relevant only if
X1A-R and X1B later pass and a separate X2 contract is frozen.
