# SIMA-DTA Architecture Blueprint

**Status:** design only. No affinity training, confirmation, or sealed access is
authorized.

## Objective

For an unseen target t, use a scaffold-diverse support set
S_t = {(d_i, y_ti)} to predict held-out scaffold-cold query affinity. The
scientific estimand is support-identified ligand reordering beyond calibration,
not generic target offset or ligand potency.

## Frozen Boundary

B0 is frozen or cross-fitted. Residue-token embeddings are frozen and cached.
At meta-test, updates are limited to a_t, b_t, and c_t. Protein and ligand
encoders never receive meta-test gradient updates. The first exact nested null
is c_t=0; the calibration comparator is a_t + b_t B0(t,d).

## M1: Support-Conditioned Hybrid Adapter

1. Protein tokens enter bidirectional Mamba blocks for linear long-context
   propagation.
2. Sparse/local/landmark attention performs content retrieval.
3. Fixed segment boundaries carry a support-memory recap.
4. A permutation-invariant support encoder consumes ligand representation,
   B0(t,d_i), and y_ti-B0(t,d_i), yielding c_t in R^q for q in {4,8}.
5. FiLM or low-rank adapter parameters are conditioned on c_t.

The prediction is:

```text
y_hat_cal  = a_t + b_t B0(t,d)
y_hat_full = y_hat_cal + phi_tilde(t,d,c_t)^T U c_t
```

A frozen ridge residualizes phi_tilde against [1,B0] on support. This prevents
the interaction channel from receiving target intercept or affine-scale credit.

Required equal-budget arms: pure Transformer, pure Mamba, hybrid,
protein-free, random task-code, B0, intercept-only, and intercept-plus-scale.
Mamba survives only on a TRAIN-frozen material effect or an explicitly frozen
memory/throughput noninferiority margin.

## M2: Query-Span Support Design

For candidate support d, calculate a label-free Jacobian j_td with respect to
c_t. Choose S to minimize average query posterior variance:

```text
mean_q j_tq^T (lambda I + sum_d_in_S j_td j_td^T)^-1 j_tq
```

Hard constraints are scaffold diversity, support-query scaffold disjointness,
chemical-neighbour caps, a non-label tie break, and selection before support
labels are exposed. Comparator policies are uniform random, scaffold-diverse
random, k-center, and uncertainty.

## M3: Counterfactual Support Identifiability

Each training episode makes four evaluations: correct support, chemistry-matched
wrong-target support, label-permuted support, and calibration-only adaptation.
The hinge objective enforces margins for correct minus wrong, true labels minus
permuted labels, and full minus calibration performance. Protein-free is a
mandatory external control.

The curriculum uses only TRAIN-side stop-gradient statistics: support
specificity, interaction gain, Jacobian effective rank, reliability, and
family/target concentration. A uniform warm-up, nonzero uniform sampling
floor, and concentration caps prevent early noisy estimates from deleting
tasks.

## Admission Logic

M1 must beat calibration. M1+M2 must beat M1 under random support. M1+M2+M3
must beat M1+M2. The complete model stops if wrong support, protein-free, or
calibration performs approximately the same; if gains occur only for chemical
neighbours; if reordering does not improve; if one family/scaffold/document/
provenance dominates; or if Mamba has no performance or efficiency value.

## Planning Cost

For 1,024 cached residue tokens, q<=8, and batch size 8: Transformer-only
4-7 GiB, Mamba-only 3-6 GiB, hybrid 5-8 GiB. A 1,000-episode GPU smoke is
estimated at 1-3 hours on the RTX 4060 Laptop GPU after token caching. These
are implementation estimates, not benchmark claims.

