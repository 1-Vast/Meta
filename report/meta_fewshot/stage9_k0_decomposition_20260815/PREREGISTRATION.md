# Stage A preregistration: what dominates cold-target k=0 error

No training. No model is modified. This stage only measures.

## Question

k=0 has now resisted a trunk redesign, a 4x budget increase, a capacity
increase, and a sequence-locality prior. Before another architecture, decompose
the k=0 error into components and find which one dominates:

1. **target-level calibration** — getting a target's mean affinity wrong;
2. **within-target shape** — getting the ligand ordering wrong inside a target;
3. **ligand representation** — how much a ligand-only predictor can do;
4. **protein representation / retrieval** — how much protein identity buys;
5. **dataset/target-prior bias** — how much a constant already achieves.

## Exact decomposition

For each target `t` with query labels `y` and predictions `p`, MSE splits
exactly:

```text
MSE(t) = ( mean(p - y) )^2   +   var(p - y)
         \_____________/         \________/
          calibration            shape
```

Reported separately, then aggregated equal-component-then-equal-target. This is
the decisive split: a calibration-dominated error and a shape-dominated error
need opposite interventions.

## Estimators, all leakage-safe

The retrieval index contains **meta_train records only**. No `meta_val` or
`meta_test` label, target, normalisation statistic or checkpoint choice enters
it. No ridge, no closed-form adaptation, no inner loop, no test-time gradient.

| estimator | uses ligand | uses protein | description |
|---|---|---|---|
| `global_mean` | no | no | meta_train mean pK |
| `ligand_prior` | yes | no | meta_train mean pK of that exact ligand, else global |
| `ligand_neighbor` | yes | no | Tanimoto-weighted kNN over meta_train ligands |
| `protein_neighbor_esm` | no | yes | ESM-pooled cosine kNN over meta_train targets, target mean pK |
| `protein_neighbor_kmer` | no | yes | 3-mer profile cosine kNN over meta_train targets |
| `dual_neighbor` | yes | yes | meta_train cells weighted by protein x ligand similarity |
| `model_f0` | yes | yes | the accepted checkpoint's zero-shot endpoint |
| `model_f0_oracle_level` | yes | yes | `model_f0` re-centred on the true target mean — isolates shape from calibration |
| `target_oracle` | — | — | true query mean; the level ceiling |

`model_f0_oracle_level` and `target_oracle` use query labels and are **diagnostic
upper bounds only**, never deployable.

## Conditioning variables

Performance is reported against: maximum protein similarity to any meta_train
target (target distance), maximum ligand Tanimoto to any meta_train ligand
(ligand novelty), assay size (cells per target), label range within target, and
homology component.

## Protein-sensitivity probe

Replace each target's protein with a cross-component meta_train donor and
measure how far `model_f0` moves. Already measured once (0.438 pK on the
grammar trunk); repeated here on the accepted `similarity_only` checkpoints for
the current active model.

## Population

`meta_val`, complete eligible bank, `evaluation_seed=73101`, k=0 episodes.
`meta_test` is not touched in this stage.

## What each outcome would imply, fixed before results

| finding | implied Stage C intervention |
|---|---|
| calibration >> shape, and `protein_neighbor` beats `global_mean` | train-only protein-conditioned calibrator |
| calibration >> shape, and `protein_neighbor` ~ `global_mean` | target-prior bias; no protein intervention can help k=0 |
| shape >> calibration, and `ligand_neighbor` strong | ligand representation is the bottleneck |
| `dual_neighbor` beats both single-modality retrievers | dual residual memory |
| `model_f0` ~ `global_mean` and protein swap moves little | cross-modal alignment failure |

## Stopping rule

This stage cannot fail; it selects the next hypothesis. Stage C proceeds with
**exactly one** intervention, chosen by the table above, and only if the
corresponding frozen oracle shows headroom of at least 5% k=0 MSE against the
accepted baseline. If no estimator shows that headroom, the honest conclusion is
that k=0 is near the achievable limit for this corpus and no intervention is
launched.
