# Stage A preregistration: target-task inner/outer-loop meta-learning

Frozen 2026-08-17, before any arm was trained. Framework inspiration only —
this is **not** a reproduction of AdaMBind, and no result here may be described
as reproducing it.

*(Date note: the Stage 0 audit ran on 2026-08-16 local and its artifact carries
that date; the local clock rolled to 2026-08-17 during this stage. Both dates
are correct for what they label.)*

Single seed. **Development evidence, directional screening only.** Not a
performance claim, not a confirmation, and `meta_test` is not read.

## Governing change

The user instruction of 2026-08-17 supersedes the repository's prior
prohibition on inner loops and deployment-time support adaptation. Inner/outer
loops and differentiable support adaptation are now permitted. Still
prohibited, and unchanged: ridge regression, pseudoinverses, closed-form
solvers, query-label adaptation at inference, and any other solver-based
shortcut. Support **inputs and labels** may drive adaptation at inference;
query labels may not, at any point.

## Stage 0 audit result

Authority: `AUDIT_DATAFLOW.json`, produced by `audit_dataflow.py`. Every
statement below is a measured number, not a reading of the source.

### The ten-step flow, as it exists today

| # | step | where | status |
|---|---|---|---|
| 1 | governed cell | `QPSMPData._governed_cells` | 5,643 meta_train / 1,411 meta_val cells; 768 meta_test withheld |
| 2 | target-task sampling | `draw_episode` | component-uniform, then target-uniform |
| 3 | support/query construction | `_unique_ligand_order` | **0 violations in 400 draws** |
| 4 | normalization | `training_label_scale` | fitted on meta_train only; verified equal to the meta_train mean |
| 5 | zero-shot forward | `encode` | 1,798,833 trainable parameters |
| 6 | support inner-loop update | — | **does not exist** |
| 7 | adapted query forward | — | **does not exist** |
| 8 | query outer loss | `train` | present, but non-episodic: one ordinary backward per episode |
| 9 | task scoring/selection | — | **does not exist**; sampling is uniform |
| 10 | checkpoint selection | `train` | **reads meta_val labels** — see below |

### Confirmed: one target per task, unique and disjoint ligands

400 random draws across k ∈ {0..5}: `multi_target` 0, `cell_overlap` 0,
`ligand_overlap` 0, `support_duplicate_ligand` 0, `query_duplicate_ligand` 0.
`materialize` additionally raises on overlap, so the property is enforced and
not merely observed.

The frozen evaluation banks are nested: `support_is_nested_prefix` true,
`query_panel_identical_across_k` true, `support_query_ligand_disjoint` true.

### Which parameters produce the zero-shot affinity

`endpoint = ligand_value(L) + protein_value(P) + interaction(P, L)`, where
`interaction = interaction_head(cat(embed, section)) + contact_weight(occupancy)`.

| branch | trainable |
|---|---:|
| `ligand_encoder` | 684,292 |
| `protein_encoder` | 433,152 |
| `grammar` | 417,264 |
| `embed` | 170,784 |
| `ligand_head` / `protein_head` | 37,249 each |
| `interaction_head` | 14,017 |
| `section` | 4,608 |
| `contact_weight` | 24 |
| `transport` | 2 |

### Adaptable subset, and why this one

**Chosen scope: `interaction_head.2.weight` + `interaction_head.2.bias` — 97
parameters, 0.0054% of trainable.**

It is the smallest subset that can change *ordering* rather than only level.
The final readout is a linear functional of a 96-wide hidden vector: moving its
weight rotates the functional and reorders ligands, while its single bias is
exactly a scalar level shift. That split is not incidental — it is the
instrument for the required k=1 question, because the same run yields a
weight-only and a bias-only counterfactual with no extra training.

Rejected for this screen: `contact_weight` (24 params) cannot express
ligand-specific shape beyond the 24 contact-type occupancies;
`interaction_head_full` (14,017) is 145x larger and would confound "adaptation
helps" with "more task-specific capacity helps"; full-backbone second-order
MAML is excluded until partial adaptation is shown structurally incapable.

### Does the current trainer touch meta_val?

**Training gradient: no.** Episodes are drawn from `meta_train` only, and the
label scale is fitted on meta_train (verified numerically).

**Checkpoint selection: yes.** `train()` evaluates meta_val banks every
`val_interval` steps and keeps the state with the best admission score. This is
disclosed rather than repaired, for one reason: repairing it would make `A0`
stop reproducing the accepted baseline, and the experiment's validity rests on
the three arms being matched. The consequence is stated explicitly and travels
with every number in the report — **every meta_val figure here is an optimistic
development estimate, not a held-out one.** Because the rule is identical
across arms it cannot manufacture a between-arm difference, which is the only
quantity this screen decides.

The A2 task selector is a different matter and is held to the strict rule: it
reads `meta_train` only, enforced by a test.

### Conflicts with functional inner-loop adaptation

Four were checked; three are benign and one is load-bearing.

1. **`transport` already performs label-based few-shot.** `SimilarityTransport`
   shifts predictions by Tanimoto-weighted support residuals. The inner loop
   fits the same support labels, so the two mechanisms partially substitute:
   as adaptation absorbs the support residual, `locked = support_y -
   support_zero` shrinks and transport does less. This is a genuine confound
   for attributing any gain, and it is why the no-adaptation control and the
   pre-vs-post decomposition are both mandatory reports rather than extras.
2. **AMP.** Inner-loop gradients under `float16` autocast can underflow.
   Adaptation runs in `float32` outside autocast. Recorded, not negotiable.
3. **`grad_clip` and the outer optimizer** act on the outer gradient only; the
   inner update is a plain SGD step with its own step size and no clipping.
4. **Warmup (`representation_warmup_fraction`) defaults to 0.0** in the
   incumbent configuration, so no phase-A branch interferes. Left at 0.0.

## The experiment

### Arms, matched in everything but the named change

| arm | training | task sampling |
|---|---|---|
| `A0` | current accepted baseline, no inner loop | uniform |
| `A1` | inner/outer loop, 97-parameter scope | uniform |
| `A2` | identical to `A1` | adaptive task selection |

Matched: seed, steps, architecture, capacity, initialization policy, optimizer,
learning-rate schedule, query panels, episode banks, evaluation code, label
scale, and the checkpoint-selection rule. The evaluator refuses to score arms
whose recorded configs differ outside the declared fields.

### Frozen settings

- seed **20260815**, one seed;
- **1,200** outer steps, matching the Stage P budget so `A0` is comparable to
  the recorded `A0repro`;
- architecture `similarity_only`;
- k cycles over {0, 1, 2, 3, 5} exactly as the incumbent trainer does;
- inner steps **1** during training; 0/1/2/3 swept at evaluation;
- inner step size selected on **meta_train only** by train-only cross-fitting
  (see below), then frozen;
- first-order meta-gradient (`create_graph=False`);
- adaptation in float32, outside autocast.

### Inner step size selection, without meta_val

A held-out objective is needed, so it is constructed inside `meta_train`:
meta_train targets are split by homology component into fit/held folds; for
each candidate step size the model adapts on a target's support and is scored
on *that same target's* held-out ligands. `meta_val` is not read. The chosen
value and the full sweep table are recorded in `RESULT.json` before any arm
runs.

### k = 0 must be exact

With k = 0 there is no inner update and the prediction must equal the zero-shot
path bitwise. This is a test, not an intention.

### The outer objective

```
outer = pre_adaptation_query_loss + post_adaptation_query_loss
```

The pre-adaptation term protects k = 0 from being traded away for few-shot
gain; the post-adaptation term is what trains fast adaptation. Ranking and
cliff auxiliaries stay at their incumbent weights and remain secondary to MSE
for this screen.

### A2 task value

Computed on `meta_train` candidates only, per outer step:

```
value = z(-post_adaptation_query_loss) + z(cosine(g_support, g_query))
```

`z` is standardization across the candidate batch of that step. Both terms are
computed under `torch.no_grad()` on a fresh candidate batch and do not
contribute gradient. Missing gradients and zero-norm gradients yield cosine
0.0 and are counted, never silently dropped. The top `episodes_per_step` of
`3 x episodes_per_step` candidates are selected. Deterministic, parameter-free.
A learned bi-level selector is out of scope unless this simple rule produces
credible evidence.

## Evaluation

`meta_val` only, k ∈ {0, 1, 2, 3, 5}, all 41 targets / 19 components, 2 draws
per target. Reported per k: MSE, RMSE, CI, Spearman, Pearson, R².

Additionally: pre- vs post-adaptation, improvement by inner-step count, the
task-weight distribution and effective number of selected tasks, the
support/query gradient-cosine distribution, peak GPU memory, wall time, and
both optimization steps and total forward/backward count.

### Counterfactuals for k > 0, all mandatory

correct support labels; permuted support labels; matched-wrong support labels;
wrong-target protein; no-adaptation control.

### k = 1: shape or level?

The 97-parameter scope splits exactly. Re-running the adapted readout with the
bias update alone gives the pure level shift; with the weight update alone
gives the pure shape change. A k=1 gain that reproduces under bias-only is a
scalar recalibration and is reported as such, not as task-specific adaptation.

## Screening decision, frozen

`PROMISING` requires all six:

1. `A1` beats `A0` at several k > 0 without materially degrading k = 0
   (k=0 MSE increase < 5%);
2. `A2` beats `A1` at k > 0, and not merely by changing task frequency;
3. MSE and ranking move in the same direction;
4. correct support beats permuted and matched-wrong support;
5. improvement grows with support size;
6. no gain traceable to query-label leakage, meta_val training, changed
   evaluation banks, or extra data.

Uncertainty on one seed is a component-paired bootstrap over the 19 meta_val
components, 9,999 draws. A single-seed interval describes episode and target
sampling, **not** retraining noise — the recorded same-config retraining spread
is 0.058 pK in k=0 MSE (`BOUNDARY_20260816.md`), and no difference below that
is claimed as real regardless of its interval.

If the inner loop fails, the report must name which of the four causes applies
— parameter scope, optimization instability, task construction, or absent
transferable support signal — and must not add architectural complexity to
rescue it. If selection fails, the uniform result stands and the selector is
rejected.

## Data

Exactly one dataset: the governed BindingDB-Ki double-cold protocol
`bindingdb_ki_double_cold_v1`. No Davis, KIBA, ChEMBL, structural corpus, or
any other affinity source; no cross-dataset support, retrieval, normalization,
label, or checkpoint. If this candidate later passes, Davis and KIBA must be
trained independently from scratch in separate experiments.
