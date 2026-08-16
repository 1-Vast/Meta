# DCST-R5 structure-only decay correction

Date: 2026-07-28  
Status: frozen after diagnosing the first R5 execution and before rerunning it

## Invalidated execution

`dcst_r5_stage1_seed1729.json` from the first R5 execution is not an admissible
test of the registered decoupled-teacher hypothesis. Its candidate source
matrix had eight exactly zero singular values.

The failure was traced to the implementation of the structure-only phase:

```python
loss = out["g"].sum() * 0.0
```

Although this term is numerically zero, it attaches the affinity-only target
and ligand output path to autograd with explicit zero gradients. Adam's
coupled weight decay then shrank those parameters during all 4,000 teacher
steps. On 64 multi-ligand source episodes, the saved teacher had:

- mean target-vector norm: `0.0`;
- mean ligand-vector norm: `0.0`;
- within-target bilinear-design RMS: `0.0`.

The same-initialization frozen control retained target and ligand norms near
`5.65` and nonzero bilinear-design variation. Consequently, the registered
teacher-frozen affinity matrix could not receive a usable gradient and stayed
at its exact zero initialization. This is an optimizer-path defect, not
negative evidence about decoupled transfer.

## Surgical correction

The zero anchor is moved from `out["g"]` to
`out["joint_interaction_logits"]`. The structure-only computation is therefore
disconnected from `theta` and from target/ligand output layers used only by
the affinity path. A regression test requires their gradients to remain
`None`, not merely numerical zero.

No data, split, label, seed, objective weight, optimization budget,
certificate rule, mechanism threshold, or control is changed. The original
R5 preregistration remains binding. The corrected execution must use a new
artifact name and the first artifact must never be treated as an outcome.

No ChEMBL affinity label was loaded or scored while diagnosing or correcting
this defect.
