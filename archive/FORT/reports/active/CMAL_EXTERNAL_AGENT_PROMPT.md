You are reviewing a failed but mechanically functional few-shot drug-target
affinity meta-learning system. Treat this as a research diagnosis, not as a
request to optimize a benchmark number.

I attached three files:

1. `research/a2s_cmal.py`: the actual executable model. The command path is
   `main.py a2s-cmal -> research.a2s_cmal`. Do not assume the repository's
   legacy `model/` directory is current; it is not imported by this command.
2. `a2s_cmal_v7_disjoint_source_meta300_seed1729.json`: the latest complete
   source-only run, including parameter inventory, gradients, parameter deltas,
   training logs, source validation, and source holdout results.
3. `CMAL_FAILURE_HANDOFF.md`: the sealed data contract, controlled failure
   chronology, quantitative shortcut/residual diagnostics, and non-negotiable
   constraints.

The scientific objective is strict. We need a protein-conditioned,
transferable, learned, non-closed-form adaptation operator. From only
`k={1,3,5}` recipient support measurements for a strictly unseen target, it
must produce target-specific and query-dependent compound-ranking improvements
over a frozen support-free DTA base. Calibration, interpolation, similarity
retrieval, ridge/GP/Bayesian posterior updates, and other fixed analytic updates
are baselines only and cannot be proposed as the final solution.

No recipient labels have been read. Do not request them, infer from them, or
recommend recipient evaluation. All diagnosis and model selection must remain
source-only.

The current model is mechanically active: gradients reach all adapter modules,
parameters change materially, arm tensor ordering has been audited, and all
wrong-support episode IDs and measurement IDs are aligned. The central failure
is scientific, not a dead computation graph. Earlier training reduced the
counterfactual loss and increased a wrong-minus-correct ranking gap while
worsening unseen-source RMSE, CI, Spearman, and NDCG@10.

Several causal repairs have already been tested one at a time:

- The final delta was forced to equal a learned query-specific scale times an
  attention-weighted measured support residual. This prevents a second
  support-independent query head and makes delta exactly zero when support
  residuals are zero.
- Raw-loss InfoNCE was replaced by frozen-base-anchored gain contrast. A wrong
  arm stops receiving reward once it is no better than the base, so the model
  cannot win merely by making wrong arms arbitrarily bad.
- A same-compound label-swap counterfactual was added. It keeps correct support
  compounds and f0 fixed while transplanting a wrong-target residual, forcing
  the model to use the compound-label relationship instead of identifying an
  arm from chemistry alone.
- Source meta-train homology components were split 1:1 into disjoint base and
  adapter tasks. Sampling is target-balanced, and the base no longer sees the
  same nested query three times through k=1/3/5.

These changes removed known shortcuts but did not yield robust transfer. At the
latest 300-step checkpoint, source holdout CI and Spearman improved by about
`+0.00303` and `+0.00898`, and all correct-vs-counterfactual holdout CI and
Spearman point estimates were positive. However, holdout NDCG@10 decreased by
`-0.00122`, holdout RMSE worsened from `1.6658` to `1.7231`, and source
validation simultaneously deteriorated by CI `-0.01011`, Spearman `-0.03235`,
and NDCG@10 `-0.01551`. The sign reversal across source splits means there is no
key positive breakthrough and no basis for publication or recipient testing.

Important data diagnostics are already known:

- Negative episode materialization and fused arm ordering are correct.
- A chemistry-only arm classifier, using no labels or protein, identifies the
  correct arm at 51.6% on meta-train and 54.0% on meta-validation versus 25%
  chance. The label-swap arm was added to address this shortcut.
- Support/query chemistry is weakly local: sampled nearest ECFP Tanimoto is
  approximately 0.223.
- Residual means are stable, but unseen-source scale is larger: support
  residual-mean SD is 1.105 vs 0.770 (+43.5%), and query within-episode residual
  SD is 1.138 vs 0.951 (+19.7%).
- Ordered meta-validation has only 18 evaluable targets in 12 homology
  components, so small point-estimate gates are unreliable.

Please do the following:

1. Audit the attached implementation line by line for any remaining objective,
   sign, masking, normalization, residual-transplant, attention, sampling, or
   evaluation defect. Explicitly distinguish confirmed defects from plausible
   research hypotheses.
2. Explain why the operator can learn correct-vs-wrong support specificity yet
   fail to make correct-support predictions absolutely better than the frozen
   base, especially for NDCG@10.
3. Decide whether query-specific ranking correction is statistically
   identifiable from these k<=5, chemically distant supports, or whether the
   present architecture simply has the wrong inductive bias. State what
   source-only evidence would distinguish those cases.
4. Evaluate the current same-compound label-swap construction. If transplanting
   a wrong-target residual vector by slot is invalid or creates a new shortcut,
   propose a chemistry-preserving alternative that is implementable at k=1 and
   does not use query labels for mining.
5. Decide whether the fixed component-disjoint base/adapter split is sound,
   should be replaced by out-of-fold base predictions, or should be removed.
   If you recommend out-of-fold predictions, specify how one adapter can consume
   fold-specific base representations without representation-coordinate drift.
6. Propose exactly one minimal architecture/objective repair and one decisive
   source-only experiment to falsify it. Do not suggest a broad hyperparameter
   sweep. Include the frozen controls, k=1/3/5 reporting, paired homology-
   component uncertainty, and an explicit stop rule.
7. Provide a line-level patch plan for `research/a2s_cmal.py`. The final operator
   must remain learned, protein-conditioned, query-dependent, support-label-
   dependent, and non-closed-form.

Success cannot mean only that wrong supports become worse or that a training
contrastive loss falls. It requires correct support to improve absolute source-
holdout ranking over the identical frozen base and to beat chemistry-controlled
label-swap support. Do not claim success from one small positive point estimate;
require paired component-level uncertainty and consistency across source
validation and source holdout. If the estimand is not identifiable with this
data, say so clearly and identify the smallest data/protocol change needed to
test it without touching recipient labels.
