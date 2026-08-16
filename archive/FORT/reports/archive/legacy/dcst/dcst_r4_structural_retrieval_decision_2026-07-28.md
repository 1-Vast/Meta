# DCST-R4 structural-retrieval decision

Date: 2026-07-28  
Decision: `STOP_R4__STRUCTURE_AFFINITY_GRADIENT_CONFLICT`

## Result

R4 made the structure representation nearly destruction-complete:

```text
true centered alignment          +0.04355
wrong-target centered alignment  -0.00534  (difference 0.04888)
wrong-ligand centered alignment  -0.01643  (difference 0.05998)
```

The ligand margin passed, while the target margin missed 0.05 by 0.00112.
Correct structural cross-entropy was 4.414 versus the uniform 5.545.

Affinity transfer failed in the opposite direction: the privileged model
certified zero of four bands and `NoPriv` certified one. The stronger
structural objective therefore improved pair-specific structure while
destroying the affinity residual directions in the shared encoder.

## Interpretation

R4 identifies a multi-task gradient conflict. It does not justify weakening
the structural threshold or the affinity certificate. The successor must
separate structure representation learning from affinity fitting and test
whether a frozen structurally identified basis supports affinity directions.

