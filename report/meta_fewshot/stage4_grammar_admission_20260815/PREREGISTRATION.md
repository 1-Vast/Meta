# Stage 4 preregistration: three-seed governed admission run

Written before any Stage 4 run. Configuration is **exactly** the Stage 3 arm C
configuration; the only change is three model seeds and the governed
nested-manifest evaluator.

## Configuration

```text
arch                     grammar (InteractionGrammarModel)
hidden_dim 384, task_dim 64, pair_dim 192, pair_latents 48, pair_heads 8
ligand_layers 5                       -> 7,294,171 trainable parameters
steps 2000, episodes_per_step 4       -> 8,000 episodes per seed
learning_rate 6e-4, backbone_lr_scale 1.0, cosine schedule (5% warmup, 10% floor)
binding_loss_weight 1.0, warmup fraction 0.0
validation bank 2 targets per component, val_interval 250
model seeds 20260812, 20260813, 20260814
```

Evaluation: `scripts/evaluate_qpsmp.py` on an outcome-redacted nested manifest
over **all eligible meta-test targets**, nested k = {0,1,2,3,5}, one common
query set per target, shared across seeds. Metrics: MSE, RMSE, concordance
index (pairwise sign accuracy), Spearman, and paired component bootstrap with
9,999 draws.

Arms: `full`, `zero_shot`, `level_only`, `sar_cut` (adaptation cut),
`permuted_state` (magnitude-matched at k=1, cyclic otherwise),
`foreign_code_state` (matched wrong support), `wrong_protein_state`.

## Admission gates

| # | requirement |
|---|---|
| 1 | k=0 does not materially regress against the retained baseline, and preferably improves |
| 2 | k=1 gain is query-specific: `full` beats `sar_cut` at k=1 |
| 3 | k=2,3,5 improve against the retained baseline in every seed |
| 4 | `full` beats `sar_cut` with a positive paired-component bootstrap lower bound |
| 5 | correct support beats permuted support and matched wrong support |
| 6 | concordance index and Spearman improve together with MSE |
| 7 | checkpoint re-evaluation reproduces the reported numbers |
| 8 | no dead trainable branch at any k |

Gate 4's bootstrap lower bound is the binding statistic. Point estimates alone
do not admit.

## Comparators

Retained baseline, frozen protocol bank (three-seed mean):
k0 2.12, k1 1.70, k2 1.34, k3 1.27, k5 1.25.

Retained baseline, wide bank (seed 20260812):
k0 3.589, k1 2.386, k2 1.993, k3 1.831, k5 1.581.

Label-only references, wide bank: global mean 3.441; ligand prior 3.119;
support mean 2.346 / 2.180 / 1.918 / 1.523 at k=1/2/3/5; oracle target-mean
level ceiling 1.100.
