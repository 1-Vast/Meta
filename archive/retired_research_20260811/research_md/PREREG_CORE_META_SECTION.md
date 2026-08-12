# Core Meta-Section preregistration

Status: design frozen before any real-data Meta-Section fit. Synthetic labels are
open for implementation validation. BindingDB affinity training and frozen
evaluation labels remain closed until execution authority is updated.

## Question and fixed model

The first experiment asks only whether a source-learned family
`m(P,L)=U^T phi(P,L)`, with `d in {1,2,3,4,5}`, permits useful unseen-target
adaptation beyond `d=0`. `phi` remains the audited 288D T-BASIS. For support
matrix `M_S` and population residual `r_S`, the only task state is

```text
a = M_S^T (M_S M_S^T + lambda I)^-1 r_S, lambda > 0.
```

No attention, Q-PMA, support encoder, MAML inner loop, new PLM, GNN branch or
high-dimensional task embedding is admitted in this experiment. Query-attention
is a later replacement hypothesis only if the linear section passes support and
biology controls but fails a preregistered query-conditioning discriminator.

## Ordered gates

1. Synthetic positive control must recover a planted low-dimensional family,
   propagate query loss to `U`, remain support-order invariant, keep section
   rank at most `k`, reject foreign/permuted support, handle off-row queries and
   reproduce analytic measurement covariance by Monte Carlo.
2. M0 uses source-only pseudo-held-out targets. Hyperparameters and stopping are
   selected without frozen evaluation targets. Arms are `d=0`, correct support,
   zero section, foreign-target support and within-task permuted labels.
3. M1 is run once on the frozen dependency-closed evaluation split only after a
   matching 21,473-cell label/feature artifact exists and authority permits the
   label reads. It adds correct, wrong-protein and ligand-only biology arms.

## Episodes and controls

Targets are sampled uniformly, not in proportion to cell count. A target may
contribute every legal `k in {1,2,3,5}` with at least three distinct query
ligands. Support and query ligand identities are disjoint; scaffold-disjoint
episodes are reported separately. Query labels are held outside the model input
and are visible only to the outer loss/evaluator.

The primary contrast is target-macro loss reduction of correct support over the
maximum of zero, foreign and permuted support. Biological admission additionally
requires correct protein to beat both wrong protein and ligand-only on identical
query sets. Primary uncertainty uses target bootstrap; dependency-component
sensitivity is mandatory. Results are reported for every `k` and `d`, with rank,
condition number, coverage and abstention rate.

`d=0..5` are all reported. `d` and ridge are selected on source-only validation
with ties resolved toward smaller `d` and then stronger ridge. Zero section and
support-free population prediction are the same arm under this contract.
Permuted labels are undefined at `k=1` and are reported as NA, not as a duplicate
control. Support draws, seeds and ligand rows never increase the inference-unit
count.

## Stop rules

- Synthetic failure: `META_SECTION_IMPLEMENTATION_NOT_VALIDATED`; do not train.
- Correct support fails a positive lower confidence bound over any support
  control: `META_SECTION_NOT_IDENTIFIED`.
- Support passes but correct protein does not beat both biology controls:
  `BIOLOGICAL_COORDINATE_NOT_AFFINITY_IDENTIFIED`.
- All core gates pass: continue to reliability propagation and compact `z`;
  do not yet alter the frozen law-valued operator.

## Required real-data materialization precondition

The frozen O1 split contains 21,473 structural cells, while the existing
T-BASIS feature bank covers the older 12,457-cell quotient corpus. A read-only
availability audit found exact numeric Ki rows for only 18,509 O1 cells; after
the frozen dependency assignment, numeric `k=5` evaluation supply is 29 rather
than 33. O1 must therefore be reopened as an estimand-matched numeric corpus
Gate before M0/M1.

Register replicate aggregation without selecting on value, persist the complete
cell/target-to-component mapping, then re-run the unchanged `>=30` and MDE
thresholds. If that passes, create matching correct/wrong-protein/ligand-only
feature rows, hashes, endpoint Ki and an access ledger. No positional or partial
join to the old feature bank is allowed.
