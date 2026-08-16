# A2S-CMAL failure handoff

Date: 2026-08-01

## Status

There is no key positive breakthrough yet. Do not read recipient labels, run the
formal five-seed recipient experiment, move the implementation into `model/`,
or publish it as the core model. All results below are source-only diagnostics.

## Execution truth

The current executable path is:

`main.py a2s-cmal -> research.a2s_cmal`

The repository's `model/` directory is a legacy Bayesian/gradient-adaptation
stack used by old scripts. It is not imported by the current CMAL command and
must not be audited as if it were the model described here.

## Scientific objective

Learn a protein-conditioned, transferable, non-closed-form meta-adaptation
operator from abundant source targets. Given only `k={1,3,5}` measured
support pairs for a strictly unseen target, it must use the support labels to
produce target-specific and query-dependent compound-ranking improvements over
a frozen support-free DTA base. Calibration, similarity lookup, interpolation,
and analytic posterior updates are baselines, not acceptable final mechanisms.

## Sealed data contract

- Episode package: `dataset/formal_training/a2s_cmal_episodes.v3`
- Content SHA-256:
  `2df5831bc8a51df93dc54531302327716fcca8900ec43f1aa37f16ed2fb9485a`
- 30,123 label-blind episodes: 23,127 meta-train, 2,595
  meta-validation, 3,456 meta-test, and 945 sealed recipient-test.
- 206 declared source targets and 63 recipient targets. One source target has
  no constructible episode, so 205 source targets occur in the episode table.
- Source splits preserve intact 40%-identity homology components.
- Support/query parent, document, measurement, and ordered-time overlaps are
  zero. All wrong-target arms preserve role/split/protocol/k and use a
  different target.
- Every label is joined by a frozen `measurement_uid`. Local runs load source
  labels only. Recipient labels have not been read in any diagnostic reported
  here.

## Current operator

The support-free base is a target-balanced pooled ligand prior plus a learned
protein-ligand residual network. The base is frozen before adapter training.

For each support measurement, the adapter computes the measured residual

`r_i = y_i - f_base(p, x_i)`

and forms pair-plus-residual tokens. A support self-attention block produces a
set state. Each query performs cross-attention to that state. The current
output is residual anchored:

`delta_q = scale_phi(p, x_q, S) * sum_i attention_phi(q,i) * r_i`.

Therefore `delta_q` is exactly zero when every measured support residual is
zero. This removed the earlier support-independent second-head shortcut while
remaining learned, protein-conditioned, query-dependent, and non-analytic.

The online counterfactual set has four arms:

1. random wrong-target support;
2. protein-hard wrong-target support;
3. support-chemistry-matched wrong-target support;
4. same-compound label swap, which keeps the correct support compounds and
   f0 values but transplants a wrong-target residual vector.

The current contrast uses frozen-base-anchored ranking gains. Correct support
is scored by `L_rank(base)-L_rank(correct)`. A wrong support is scored by
`relu(L_rank(base)-L_rank(wrong))`, so making a wrong arm worse than the base
cannot continue to improve the contrastive objective.

The latest diagnostic also separates meta-train homology components into a
fixed 1:1 base-pretraining/meta-adapter split: 54 components/79 targets for the
base and 54 components/84 targets for the adapter, with zero overlap. Sampling
is target-balanced. The base is trained only from the k=5 view because nested
k=1/3/5 episodes repeat the same query set.

## Confirmed implementation facts

- The original model had 1,420,290 parameters: 1,106,049 in the base phase and
  314,241 in the adapter phase.
- All adapter modules had nonzero gradients after the zero-initialized output
  layer opened, and their relative parameter changes were about 10%-32%.
- Tensor arm ordering is correct: arm-major concatenation, repeated recipient
  protein/query tensors, and reshape back to `[arms, batch, query]` agree.
- All 29,178 source episodes have correctly aligned negative episode IDs,
  targets, support budgets, parent IDs, and measurement IDs.
- Thus the failure is not explained by a disconnected gradient graph or a
  wrong-arm materialization bug.

## Quantified data and objective problems

1. A chemistry-only classifier that ignores labels and protein identifies the
   correct arm from support-query ECFP similarity with 51.6% meta-train and
   54.0% meta-validation top-1 accuracy; chance is 25%. The original three
   wrong-target arms therefore contain a strong chemistry shortcut.
2. Mean source residuals are similar across train/validation, but their scale
   shifts. The SD of episode support-residual means is 0.770 in train and
   1.105 in validation (+43.5%). Query within-episode residual SD is 0.951 vs
   1.138 (+19.7%).
3. Correct support and query compounds are chemically distant: sampled nearest
   ECFP Tanimoto is about 0.223 in both train and validation.
4. Ordered meta-validation contains only 18 evaluable targets in 12 homology
   components. Point estimates are noisy; a gate based only on `mean > 0` is
   insufficient evidence.

## Failure chronology

### Original operator, shared source tasks, raw-loss InfoNCE

Report: `a2s_cmal_v3_gradient_audit_seed1729.json`, base 300 steps, adapter 500.

- Training counterfactual loss: about 1.386 -> 1.068.
- Training wrong-minus-correct ranking gap: 0 -> +0.2229.
- Source meta-validation adapted minus frozen base:
  RMSE +0.09846 (worse), CI -0.00988, Spearman -0.03085,
  NDCG@10 -0.01262.
- Correct-support CI advantage over random/protein-hard/chemical arms:
  -0.00053/-0.00191/-0.00464.

This proved that a growing relative gap could be obtained without transferable
positive adaptation.

### Minimal controls

- Original model at adapter step 100 already decreased CI/Spearman/NDCG, so
  failure was not only late 500-step overfitting.
- Setting the counterfactual weight to zero still decreased absolute ranking,
  showing that the positive adapter path itself could bypass support labels.
- Residual-anchoring the final delta changed RMSE from slightly worse to
  -0.0138 (better) and made CI slightly positive, but Spearman/NDCG and one
  protein-hard comparison remained negative.
- Base-anchored capped-gain InfoNCE removed the incentive to destroy wrong
  arms, but did not produce stable correct-support specificity at steps 100 or
  300.

### Same-compound label-swap counterfactual

At adapter step 100, all source meta-validation point estimates became
positive:

- absolute CI +0.00062, Spearman +0.00187, NDCG@10 +0.00178;
- all 12 correct-vs-four-counterfactual point estimates were positive.

However, almost every component-bootstrap 95% interval crossed zero, and the
previously untouched source meta-test did not confirm the direction:

- absolute CI -0.00036, Spearman -0.00110, NDCG@10 -0.00020;
- several label-swap and protein-hard comparisons were negative.

At step 300, validation overfit and both gates failed. This was a directional
repair, not a key positive breakthrough.

### Component-disjoint base/meta-adapter training

Latest report:
`a2s_cmal_v7_disjoint_source_meta300_seed1729.json`.

At adapter step 100:

- source validation: CI +0.00457, Spearman +0.00888, NDCG@10 -0.00218;
- source holdout: CI -0.00304, Spearman -0.00820, NDCG@10 -0.00336.

At adapter step 300:

- source validation absolute gain: CI -0.01011, Spearman -0.03235,
  NDCG@10 -0.01551;
- source holdout absolute gain: CI +0.00303, Spearman +0.00898,
  NDCG@10 -0.00122;
- all four source-holdout correct-support specificity comparisons were positive
  for CI and Spearman, but absolute NDCG still decreased;
- source-holdout RMSE changed from 1.6658 to 1.7231 (worse).

This branch learned target-specific relative behavior but not a stable,
beneficial ranking correction. Its sign reverses between validation and
holdout, and it overfits by step 300.

## Open technical question

Why can the operator distinguish correct from wrong support after training, yet
fail to make correct-support predictions consistently better than the frozen
base on unseen source targets, especially for NDCG@10?

The leading possibilities are:

1. the attention operator is learning task/arm recognition rather than a
   transferable residual field over query chemistry;
2. k<=5 chemically distant supports do not identify a query-local ranking
   correction in the present ECFP/descriptor representation;
3. the positive MSE plus global pairwise ranking objective emphasizes target
   intercept and broad ordering, while NDCG@10 requires top-tail behavior;
4. the fixed 1:1 base/meta split makes the base weak and causes unstable
   residual geometry; an out-of-fold base may be necessary instead;
5. same-compound residual transplantation may be a useful destructive control
   but an imperfect training negative because donor residual slots do not
   correspond to the same ligand chemistry;
6. evaluation has too few homology components for point-estimate gating, and
   the apparent gains may be noise.

## Non-negotiable constraints for the next analysis

- Do not read or infer from recipient labels.
- Do not run recipient evaluation or formal five-seed training.
- Do not propose calibration, kNN retrieval, interpolation, ridge/GP/Bayesian
  posterior, or any other fixed/closed-form update as the final model. Those are
  baselines only.
- Do not recommend a broad hyperparameter sweep. Prefer one structural change
  and one falsifying source-only experiment.
- Correct support must improve absolute unseen-target ranking over the frozen
  base and beat chemistry-controlled label-swap support. Merely harming wrong
  arms is not success.
- Any positive claim must survive both source validation and source holdout,
  use paired homology-component uncertainty, and report k=1/3/5 separately.
- Do not move code into `model/` or publish it until a robust positive source
  mechanism gate is achieved.

## Requested analysis output

1. Audit the current equations and code path for any remaining objective,
   gradient, masking, normalization, or evaluation defect.
2. Rank the remaining causal hypotheses using the evidence above.
3. Decide whether query-specific ranking improvement is identifiable from this
   episode construction at k<=5; distinguish statistical non-identifiability
   from an inadequate neural operator.
4. Propose the smallest defensible architecture/objective repair. It must remain
   learned, protein-conditioned, query-dependent, and non-closed-form.
5. Specify one source-only experiment that can falsify the proposed repair,
   including controls, exact metrics, component-level uncertainty, and a stop
   rule. No recipient outcomes may be used for model selection.
6. State whether component-disjoint training should be retained, replaced by
   out-of-fold base predictions, or removed, and explain how to avoid changing
   representation geometry across folds.
7. Evaluate whether the label-swap arm is constructed correctly; if not,
   provide an implementable chemistry-preserving alternative that also works at
   k=1.
8. Give a line-level patch plan for `research/a2s_cmal.py`, but do not claim that
   the mechanism works without source holdout evidence.

## Recommended chat attachments

Attach only these three files first:

1. `research/a2s_cmal.py` -- the actual current executable model;
2. `reports/active/a2s_cmal_v7_disjoint_source_meta300_seed1729.json` -- the
   latest full parameter/gradient/metric report;
3. this handoff file -- data contract, controlled failures, and constraints.

If the chat accepts a fourth attachment, add `research/a2s_cmal_data.py` so the
agent can audit episode construction. Do not attach the legacy `model/` folder
unless explicitly asking for a historical comparison.
