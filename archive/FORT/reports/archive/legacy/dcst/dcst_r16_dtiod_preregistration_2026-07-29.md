# DCST-R16 DTIOD preregistration

Date: 2026-07-29  
Status: frozen before implementation and execution

## Claim

`Discrete Tangent Interaction Operator Distillation` tests whether local
segment-by-substructure responses transfer better than the absolute R6
mechanism state. It is independent of SISMT and cannot inherit a pass from it.

## Discrete intervention

The frozen R6 mechanism distribution `P(t,d)` is the teacher function.

- A protein intervention replaces one raw ESM segment by the mean of the
  other 31 segments; pooled ESM is unchanged.
- A ligand intervention clears one active Morgan bit; physicochemical
  descriptors are unchanged.

For segment `r` and active bit `a`:

```text
Omega_ra =
    P(t,d)
  - P(t_without_r,d)
  - P(t,d_without_a)
  + P(t_without_r,d_without_a).
```

This finite difference is computed on CUDA. It is not a raw embedding
Jacobian.

## Mask selection

For PLINDER rows, the semantic segment is the highest-mass segment in the
privileged `32 x 8` structural target. The semantic ligand bit is the active
Morgan bit receiving maximum R6 substructure attention from that segment.
Three additional active bits are chosen deterministically by descending
attention.

Controls use:

- random segments with the same four-bit budget;
- random active bits with the same segment budget;
- wrong target;
- within-target wrong ligand;
- R6-NoPriv;
- a matched random-frozen teacher.

For ChEMBL label-blind support, segments and bits are selected only by frozen
teacher attention; no structural or affinity outcome is available or used.

## Gate T1: teacher tangent semantics

On PLINDER development:

1. Priv semantic `Omega` RMS must be at least `1.20x` NoPriv;
2. it must be at least `1.20x` each matched random-mask RMS;
3. wrong target and wrong ligand must each remove at least 30% of the
   correct-mask response;
4. at least 20 exact targets must contribute.

Failure returns `STOP_DTIOD_T1_NO_PRIVILEGED_TANGENT`.

## Gate T2: held-intervention student

Only if T1 passes, train a sequence/ligand student on PLINDER TRAIN using

```text
L = KL(P_teacher || P_student)
  + lambda_tangent * MSE(Omega_teacher, Omega_student),
```

with `lambda_tangent=1.0`, four mask pairs per row, and 4,000 CUDA steps.
No affinity outcome enters this loss. PLINDER development targets and one
deterministically held active-bit intervention per row are excluded from
training.

T2 passes only if held-target, held-mask flattened `Omega` has:

- samplewise median cosine at least `0.50`;
- global Pearson correlation at least `0.50`;
- Priv-student performance at least `0.10` above a map-only student trained
  with the same budget.

## Gate T3: tangent support

Build a TRAIN-only SVD basis retaining 95% source tangent variance, capped at
32 dimensions. Compare source-versus-ChEMBL-train domain AUC and the same
partial-transport ESS/coverage audit used by R15.

T3 passes only if:

1. tangent domain AUC is at least `0.05` lower than the R14 absolute
   privileged centered-moment AUC `0.941130`;
2. at transported mass at least `0.20`, source-target and ChEMBL-target ESS
   are both at least `88`;
3. tangent target coverage is at least 10 percentage points above the
   absolute-state coverage under the identical selector;
4. T1 attribution remains positive after projection.

Failure returns `STOP_DTIOD_TANGENT_NOT_MORE_TRANSPORTABLE`.

## Conditional Stage 2

Only after T1-T3 pass, freeze the student, tangent basis, and B0. Fit an
exact-null linear tangent readout with matched NoPriv, map-only, random,
wrong-target, and wrong-ligand controls. The existing `0.0586` MDE,
positive grouped-bootstrap LCB95, RMSE safety, and 70% destruction-removal
requirements remain binding. Confirmation and sealed test remain untouched.

