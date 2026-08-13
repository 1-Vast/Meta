# QPSMP-BPSF Architecture Refactor

Date: 2026-08-13

## Outcome

The active Cold Target few-shot model is now a single trainable architecture:

1. frozen/cached protein PLM states and a trainable ligand GNN;
2. a bond-aware bipartite atom-residue pair field;
3. independent endpoint and meta-section latent readouts;
4. a learned quotient-preserving support-set operator;
5. scalar endpoint prediction with an empirical-Bayes level correction.

The active model contains no pooled-interaction mode, one-pass atom-residue
mode, analytic ridge adapter, or V1 deployment branch. Those implementations
are isolated under `legacy/retired_qpsmp` and are not imported by `model`,
`scripts`, or the root CLI.

## Capacity

The default model has 1,660,167 trainable parameters for a 640-dimensional
protein bank. Its default configuration is:

- hidden width: 128;
- pair width: 64;
- pair blocks: 3;
- endpoint/section latent slots: 16 each;
- section dimension: 32;
- support-operator width: 128 with 2 attention blocks;
- ligand GNN layers: 3.

This is a compact small-to-medium model. Pair-field activation memory, rather
than parameter memory, is the practical RTX 4060 constraint. The implementation
therefore retains active-atom compaction, AMP, and pair chunking.

## Structural Repairs

- Pair blocks now use molecular adjacency, local residue propagation, and
  pair-to-token feedback before refreshing the pair field.
- Endpoint and section latents are independent. Section-only meta-training can
  update the full section latent branch while preserving the endpoint branch.
- The learned support state is constructed from centered residual evidence and
  the centered support row span. Constant residuals and `k=1` produce exact
  zero SAR.
- Query section coordinates are centered by the support center. This restores
  joint section-coordinate translation invariance and prevents the SAR channel
  from carrying an arbitrary section gauge offset.
- The support operator is permutation invariant and contains no analytic solve.

## Verification

Focused regression:

```text
21 passed
```

Coverage includes mask behavior, atom/residue emptiness, coordinate gauge,
support permutation, row-span membership, exact null, `k=1`, CUDA AMP, channel
decomposition, and gradients through the pair and support operators.

The final CUDA implementation smoke used one optimizer step, one episode per
step, support size 2, query size 2, and seven development episodes:

```text
AMP                         true
peak allocated CUDA memory 182.40 MiB
validation MSE             1.688635 pK^2
full MSE                   2.131852 pK^2
level-only MSE             2.364837 pK^2
SAR-cut MSE                2.072423 pK^2
SAR gain                  -0.059430 pK^2
```

This smoke validates execution only. The negative SAR gain means it does not
authorize a performance, target-specificity, biological, or safety claim.

## SCI Design Assessment

The architecture is now adequate as an SCI research model candidate, narrowly
framed as a bipartite pair-section representation coupled to an amortized,
quotient-preserving few-shot operator. Its novelty is not attention, a GNN, or
parameter count. The defensible contribution is the structural restriction of
target adaptation to centered support-observable directions while retaining a
shared scalar endpoint potential.

It is not yet adequate as an SCI result. Publication-level evidence still
requires:

- at least three frozen seeds and nested `k={1,2,3,5}` episodes;
- target-to-component macro CI/RMSE and component bootstrap intervals;
- full versus level, SAR-cut, simple learned adapter, and analytic comparator;
- correct versus permuted, foreign-support, and wrong-protein controls;
- source-only pair-geometry pretraining with receptor, ligand, scaffold, and
  structure-entry dependency closure;
- a fresh dependency-closed confirmation cohort.

Passing only endpoint RMSE would not identify the meta-section innovation.
The learned SAR term must improve endpoint and ranking risk while also beating
the support/protein controls. Until then, G2/G3 and biological claims remain
fail closed.
