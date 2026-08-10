# Current task

## Goal

Build a trainable few-shot drug-target affinity model that uses large open
datasets without confusing dataset size with biological identifiability.

## Established

```text
BINDINGDB_QUOTIENT_TRAINING_READY
REAL_OPEN_DATA_LINEAR_WITNESS_EXECUTED
POPULATION_SHARED_288D_AFFINITY_DIRECTION_NOT_IDENTIFIED
FEWSHOT_TARGET_SECTION_NOT_YET_TESTED
```

## Next registered design target

Prepare a separate preregistration for a single low-dimensional
target-coefficient model on frozen T-BASIS coordinates:

```text
q_t(P,L) = w0^T phi(P,L) + a_t^T U^T phi(P,L),  d <= 5
```

Use dense profiling panels for ordinal source pretraining and endpoint-specific
BindingDB/Klaeger data for quantitative constraints. Evaluate unseen targets at
`k=1/2/3/5`; adaptation must stay in the support-observable row space and must
abstain off coverage.

## Do not do

- do not rescue the failed shared linear witness by lowering a Gate;
- do not pool Ki, Kd, Kdapp, inhibition and displacement as one label;
- do not add PLMs, GNN branches, cross-attention, KG, pose or typed energy heads;
- do not open Davis/recipient confirmation labels;
- do not modify `A(F,z)=K(B(z)F(z))`, CSMO, Band or production `z`.
