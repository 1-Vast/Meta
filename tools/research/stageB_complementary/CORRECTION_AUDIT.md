# Correction audit of Stage A

Stage A's artifacts are **not modified**. `tools/research/stageA_innerloop/`
remains the authoritative record of what was run and reported on 2026-08-17.
This document records eight defects in that stage's *analysis and
interpretation*, the corrected treatment, and where each correction is
implemented and tested.

Two of the eight (1, 2) change measured numbers and therefore require a re-run;
Stage B supplies it. Five (3, 4, 5, 6, 7) are interpretation errors correctable
from the existing evidence. One (8) is a statistical framing error.

| # | defect | class | corrected in | test |
|---|---|---|---|---|
| 1 | matched-wrong anchored to post-adaptation | **measurement** | `arms.py`, `evaluate_stageb.py` | `test_matched_wrong_anchored_to_pre_adaptation_is_arm_independent` |
| 2 | A0 counterfactuals never summarized | **measurement** | `evaluate_stageb.py` | `test_every_arm_reports_every_control` |
| 3 | k=1 conditioning formula omitted the bias | formula | `residual.conditioning_alpha` | `test_alpha_includes_the_adapted_bias` |
| 4 | "weight-only = pure shape" | interpretation | `residual.centered_shape` | `test_a_weight_only_update_is_not_pure_shape` |
| 5 | wrong-protein framed as an inner-loop test | interpretation | this document, `REPORT.md` | — |
| 6 | "meta_test unread/untouched" | governance wording | this document | `test_meta_test_status_is_the_audited_string` |
| 7 | "A1 improves every metric at every k" | factual | this document | `test_k5_ci_regression_is_recorded` |
| 8 | component bootstrap read as seed uncertainty | statistics | this document, `REPORT.md` | — |

---

## 1. Matched-wrong support was anchored to the wrong prediction

**Defect.** `stageA_innerloop/evaluate.py:196` builds the magnitude-matched
wrong label as `support_y - 2 * (support_y - base_output["post_adaptation_support"])`.
`base_output` is the arm's *own operating condition*: `steps0` for A0 (where
post equals pre) but `steps1` for A1. So A0's corruption was anchored to its
pre-adaptation prediction while A1's was anchored to its **post**-adaptation
prediction.

**Why it matters.** Adaptation moves the support prediction *toward* the label,
so A1's `support_y - post_adaptation_support` is smaller than A0's. A1 therefore
received a **less corrupted** control. The reported "correct beats matched-wrong"
gap for A1 was measured against a weaker adversary than A0's, and the two arms'
gaps are not comparable.

**Direction of the bias.** Against A1: a weaker corruption makes the control
*less* bad, shrinking A1's apparent advantage. Correcting it should widen A1's
gap, not narrow it. The defect did not manufacture the Stage A result, but it
makes the A0/A1 comparison invalid.

**Correction.** Anchor to `pre_adaptation_support` — the shared-initialization
prediction — for every arm. Corruption magnitude then depends only on the
episode and the shared trunk, so all arms face an identical adversary.

## 2. A0's counterfactuals were computed but never reported

**Defect.** `evaluate.py` records `permuted_support` and `matched_wrong_support`
rows for **all three** arms, but the summary block at line 353 iterates
`for arm in ("A1", "A2")` and scores each control against the `"steps1"`
condition. A0 was omitted, and `steps1` is not A0's operating condition anyway.

**Why it matters.** Without A0's controls there is no way to separate the
label dependence the **incumbent Tanimoto transport already has** from the
*incremental* dependence the inner loop adds. Stage A's headline that "support
labels genuinely matter" is a statement about the combined system, not about the
inner loop. The incumbent transport is explicitly label-driven by construction,
so a large correct-vs-permuted gap was expected for A0 too.

**Correction.** Report every control for every arm, at each arm's own operating
condition, and report the **difference of differences** — the incremental
label dependence attributable to adaptation.

## 3. The k=1 conditioning formula omitted the adapted bias

**Defect.** `CONDITIONING.json` reports `alpha = 2 * lr * ||h||^2`. Stage A
adapts the weight **and** the bias. The bias gradient of a squared error is
`2 (p - y) * 1`, contributing exactly `1` alongside `||h||^2`.

**Correction.** `alpha = 2 * lr * (||h_support||^2 + 1)`. The published A0 alpha
of 1.514 is an **underestimate**; the corrected value is larger by exactly
`2 * lr = 0.2`, giving **1.714**, and A1's 0.241 becomes **0.441**. The
qualitative conclusion (A0 overshoots, A1 does not) survives, and A0 moves
closer to the oscillation boundary rather than further from it.

**Second, larger error in the same claim.** Stage A said the support contraction
"predicts" the query MSE and that the alternating sweep followed "exactly". It
does not. One inner step moves a query by

```
delta_f(q) = -2 * lr * r_support * (h_query . h_support + 1)
```

which depends on the query-support inner product and therefore differs for every
query in the panel. A single scalar governs the **support** residual only.

**Corrected claim.** A0's alpha exceeds 1, so each step overshoots and flips the
sign of the support residual; the query panel inherits a sign-alternating
correction whose per-query magnitude varies. That is *consistent* with the
observed alternating sweep (1.5352 / 2.8626 / 1.4349 / 2.2049) and is offered as
a **mechanistic hypothesis**, not as a quantitative prediction of query MSE.

## 4. "Weight-only equals pure shape" is false

**Defect.** Stage A's REPORT.md labels the bias-only update "pure level" and the
weight-only update "pure shape". The first is right; the second is not. A
weight-only step moves query `q` by `-2 lr r_s (h_q . h_s)`, whose within-episode
mean is `-2 lr r_s (mean_q h_q . h_s)` — zero only if the mean query activation
happens to be orthogonal to the support activation.

**Correction.** Shape is the query correction **after removing its own
within-episode mean** (`residual.centered_shape`). Level and shape are measured
from the correction itself, not from which parameter produced it. Stage A's
"81% of the k=1 gain is a level shift" conclusion is unaffected in direction —
it was derived from the bias-only arm, which really is pure level — but the
"weight-only = shape" half of that decomposition is withdrawn.

## 5. Wrong-protein is a full-system perturbation

**Defect.** Stage A reported the wrong-protein gap under "protein specificity" in
a section about the inner loop, implying it tested the adapter.

**Correction.** Substituting the donor protein changes the residue encoder, the
contact grammar, `embed`, the zero-shot endpoint, the support residual and hence
the transport — everything except the ligand graphs. It is a whole-system
counterfactual. The finding stands as measured (A1's gap is resolved at k=2,3,5;
A2's is not) but it licenses **no claim about the adapter specifically**.

## 6. `meta_test` status

**Defect.** Stage A's REPORT.md says "`meta_test` was not read" and
"unread/untouched".

**Correction.** The audited status is: **logical exclusion after parsing, with an
open process-isolation incident.** `QPSMPData` decompresses and parses every
`meta_test` label on construction and then discards it; the labels are in
process memory transiently. No `meta_test` label entered any fitting, selection
or reported metric in Stage A (`included: false`, `evaluated: false`, 768 cells
withheld). "Untouched" and "unread" are both stronger than the evidence
supports and must not be used.

## 7. A1 does **not** improve every metric at every k

**Defect.** Stage A's REPORT.md and RESULT.json state that A1 "beats the
accepted baseline at every k on every metric".

**Correction — this is false.** At k=5, **A1's CI is 0.62951 against A0's
0.63140**, a regression of 0.0019. The paired contrast is −0.0019
[−0.0193, +0.0144], unresolved, but the point estimate is negative and the claim
of universal improvement is withdrawn. Spearman at k=5 still favours A1
(0.3549 vs 0.3500), so the two ranking metrics disagree at the largest support
size — which is itself worth reporting rather than smoothing over.

Corrected statement: *A1's point estimates favour it at every k on MSE, Pearson
and Spearman, and on CI at k=0,1,2,3 but not at k=5.*

## 8. Component bootstrap is within-checkpoint uncertainty

**Defect.** Stage A reported component-paired bootstrap intervals alongside a
single-seed result in a way that invites reading them as the total uncertainty.

**Correction.** The bootstrap resamples the 19 `meta_val` components for **one
trained checkpoint**. It describes target and episode sampling variability
*conditional on that checkpoint*. It contains no information about retraining
variance, which the record measures separately at **0.058 pK² in k=0 MSE** for
an identical configuration. A single-seed interval excluding zero would still
not establish a reproducible effect. Only matched multi-seed runs can, and the
admission rule requires a hierarchical seed/component interval.

---

## What survives the corrections

Nothing here overturns Stage A's two decisions. The task selector remains
rejected (defects 1-8 do not touch the A2-vs-A1 contrast, which was uniformly
negative across every metric and k). The inner loop remains an unresolved weak
positive, now with one metric-and-k cell (CI at k=5) explicitly negative.

What the corrections **do** change is that Stage A's support-label counterfactual
cannot be attributed to the inner loop, because the incumbent's own transport is
label-driven and was never measured as a control. Stage B measures it.
