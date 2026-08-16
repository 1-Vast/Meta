# Stage 3 preregistration: budget, schedule and diagnosed capacity

## Stage 2 outcome carried forward

On the wide bank (42 episodes, all eligible meta-test targets), at matched seed
and 60-step budget:

| arm | k0 | k1 | k2 | k3 | k5 | s/step | peak MB |
|---|---:|---:|---:|---:|---:|---:|---:|
| `bpsf` | 4.085 | 2.631 | 2.146 | 1.866 | 1.611 | 5.57 | 6053 |
| `grammar` | 3.996 | **2.204** | **1.829** | **1.616** | **1.372** | 4.70 | 1698 |

Passed: R1 (k0 no worse), R1b (cost), R3 (k=1 query-specific channel live and
useful: `sar_cut` 2.538 against `full` 2.204), R4 (label permutation worse than
correct at every k), R6 (k>=2 not worse).

**Failed: R1 secondary and R2.** Zero-shot spread across queries is 0.107 pK
against a 0.20 pK threshold and a 0.93 pK label spread, and the wide-bank
wrong-protein zero-shot gap is -0.019, not >= 0.05. The zero-shot endpoint is
still close to a constant and protein specificity is **not** established.

Decision: the few-shot mechanism is **not rejected** — it dominates the control
on every preregistered MSE comparison and beats the label-only support-mean
reference at every k, which the control does not. The zero-shot claim **is**
rejected at this budget. R1-secondary and R2 are carried forward as blocking
gates for any zero-shot or protein-specificity claim in Stage 4.

## Stage 3 changes, each with its diagnosed role

| # | change | diagnosed role | evidence |
|---|---|---|---|
| 1 | ligand-bank LRU holds all shards | episode materialization cost 1,158 ms/episode, ~80% of every training step, from reloading 0.62 MB npz shards | `profile`: 1,158 ms -> 15.2 ms; `THROUGHPUT.json`: 5.57 -> 0.75 s/step (bpsf), 0.205 s/step (grammar) |
| 2 | steps 60 -> 2000 | zero-shot endpoint is degenerate after 240 episodes | Stage 0 spread 0.065 pK; Stage 2 spread 0.107 pK |
| 3 | cosine schedule, 5% warmup, 10% floor, `backbone_lr_scale` 1.0, lr 6e-4 | the retained 0.25 backbone scale was tuned for a warm-started trunk; both trunks here train from scratch and the encoders are the modules with no learned signal | Stage 0 gradient table; trunk-capacity sweep stable band |
| 4 | capacity 1.82M -> 7.29M (hidden 384, embed 192, 48 contact types, 5 ligand layers) | ligand tower carries no affinity information (`ligand_only` 2.082 versus 2.063 constant); protein path is one linear map | Stage 0 sections 3.1, 3.2 |
| 5 | validation bank 1 -> 2 targets per component | 400-step probe improved a 6-episode validation score while worsening every test metric | Stage 0 section 4.1 |

The frozen protocol test bank (`eval_targets_per_component=1`,
`evaluation_seed=73101`) is **unchanged**. `val_targets_per_component` is a new,
separate field so widening selection cannot touch it.

## Arms (all 2000 steps, seed 20260812, identical schedule and banks)

| arm | architecture | size | isolates |
|---|---|---|---|
| A | `bpsf` | retained 3.79M | budget + schedule on the retained architecture |
| B | `grammar` | 1.82M | architecture, at matched size-class and budget (A vs B) |
| C | `grammar` | 7.29M | diagnosed capacity (B vs C) |

## Gates

Proceed to Stage 4 with arm C only if, on the wide bank:

1. C beats A at every k in {1,2,3,5};
2. C beats the label-only support-mean reference at every k in {1,2,3,5}
   (2.346 / 2.180 / 1.918 / 1.523);
3. `full` beats `sar_cut` at k=1 (query-specific adaptation is real);
4. magnitude-matched permuted support is worse than correct support at every k;
5. no zero-gradient trainable tensor at any k.

Zero-shot and protein-specificity claims additionally require the carried-forward
R1-secondary (spread > 0.20 pK) and R2 (wrong-protein zero-shot gap >= 0.05).
Failing those does not block a few-shot claim, but it does block any statement
that the model performs protein-conditioned zero-shot recognition.
