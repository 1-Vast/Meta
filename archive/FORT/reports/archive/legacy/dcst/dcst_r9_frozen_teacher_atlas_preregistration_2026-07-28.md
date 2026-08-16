# DCST-R9 frozen-teacher atlas preregistration

Date: 2026-07-28  
Status: frozen before implementation and training

## Hypothesis and cross-stage interface

R9 composes the two previously established source properties without adding a
new encoder:

1. load the train-fitted segment interaction representation from the accepted
   R6 source checkpoint, whose held-source structural mechanism passed;
2. replace its absolute `32 × 8` affinity matrix with R8's exact frozen
   source-train ESM atlas;
3. freeze every segment, ligand-token, substructure-attention, and structural
   head parameter;
4. reset an `8 × 8` atlas energy matrix to the exact null and train only that
   matrix for 4,000 source-train affinity steps.

This is the registered **Frozen-Teacher Atlas (FTA)** interface:

```text
R6 identifiable segment measure
    -> immutable nearest-centroid atlas transport
    -> trainable exact-null atlas energy
    -> held-source destructive spectral certificate
```

The source checkpoint is
`dcst_r6_stage1_seed1729.pt`; its representation weights were trained only on
source train data. The checkpoint's source certificate selected the route but
is not copied into R9: R9 recomputes its own atlas certificate on the held
source split.

## Matched no-privileged control

FTA-NoPriv loads the matched R6 no-privileged representation state from the
same checkpoint, applies the identical frozen atlas, freezes its full
representation, resets its atlas matrix, and receives the identical 4,000
affinity steps. It differs only in whether the upstream R6 representation saw
privileged structural labels.

R8 joint-atlas and R6 absolute-position results remain historical interface
controls. R9 adds no trainable adapter or router.

## Frozen source gate

Before another ChEMBL training run:

1. the loaded privileged representation reproduces the complete registered
   segment structural mechanism pass;
2. FTA certifies at least one atlas role×type spectral band;
3. FTA certifies strictly more bands than FTA-NoPriv.

The atlas hash, R6 checkpoint hash, compatible state-key set, and downstream
label non-access must be reported. Failure stops R9.

## Conditional Stage 2

On source pass, Stage 2 uses the same frozen FTA transfer representation and
the existing role atlas. A downstream atlas residual may train as in R6; all
controls, MDE, paired bootstrap, destruction, RMSE, negative-transfer,
confirmation, and sealed policies are unchanged.

No confirmation access is authorized by a source pass.
