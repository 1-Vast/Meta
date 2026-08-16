# A2S-DTA Independent Review And V2 Source-Gate Decision

Date: 2026-08-01  
Decision: `NO_GO_INFORMATION_NOT_ADMITTED`  
Scope: source-only development evidence; locked-source and recipient labels remain sealed.

## Part 1 - Final Scientific Problem

**FACT:** The current ChEMBL construction provides passive support sets with
`k={1,3,5}`. It does not provide verified project IDs, real DMTA cycles, or the
candidate pool that existed when each compound was selected.

**INFERENCE:** The current task asks a few chemically distant, often
cross-assay measurements to identify a continuous target-specific residual
field. That adaptation object is not credibly identifiable under the observed
support/query geometry.

**HYPOTHESIS:** A defensible transferable object must be finite and auditable:
either a bounded rank intervention supported by assigned label contrasts, or a
same-assay transformation-response rule with a structural unsupported state.
Neither object has yet passed its information-admission test.

The corrected research question is:

> Under what predeclared chemical and assay relations do correctly assigned
> support labels contain incremental information about held-out query order,
> beyond an identical support-free base, and can a source-learned bounded
> operator use that information without recipient optimization or an analytic
> posterior update?

## Part 2 - What Failed And Why

### Established CMAL result

**FACT:** CMAL is mechanically active, but both source gates failed. Source
validation changed by CI `-0.01011`, Spearman `-0.03235`, and NDCG@10
`-0.01551`. The repeatedly inspected source holdout changed by CI `+0.00303`,
Spearman `+0.00898`, and NDCG@10 `-0.00122`, while RMSE worsened.

**INFERENCE:** Correct support sometimes beating wrong support only establishes
arm specificity. Because the adapted result did not reliably beat the same
frozen base, it is not evidence of beneficial target adaptation.

### V1 information-gate limitation

**FACT:** The v1 provenance union created one 380-target component. It occupied
97.0% of fit-role OOF held-out rows. Real label and assignment gaps had no
positive lower bound, while synthetic controls passed.

**INFERENCE:** The v1 null was non-confirmatory because the statistical design
had inadequate independent-component balance.

### V2 correction and result

**FACT:** `research/a2s_source_lock_v2.py` assigns homology components before
provenance handling, then quarantines rows carrying document or assay tokens
that cross roles. It reads no affinity column.

**FACT:** The v2 lock retained 68,782 of 185,591 source pKi rows (37.1%). The
retained fit/probe/locked roles contain 222/110/107 homology components and
36,609/16,068/16,105 rows. Target, homology, document, and assay overlap across
roles is zero. Five fit OOF folds contain 7,322/7,322/7,322/7,322/7,321 rows,
so the maximum fold share is 20.001%.

**FACT:** The diagnostic used separately fitted, equal-capacity G0 and G1 ridge
probes. The sole difference is whether the five support-label features are
masked in both training and evaluation. This corrects the older evaluation-only
masking design.

| k | Probe components | Delta_label rank-loss 95% CI | Delta_assign rank-loss 95% CI | Synthetic Delta_label |
| -: | -: | ---: | ---: | ---: |
| 1 | 69 | `[-0.00426, +0.01271]` | undefined | `+0.31635` |
| 3 | 45 | `[-0.00865, +0.01185]` | `[-0.00087, +0.01366]` | `+0.22976` |
| 5 | 36 | `[-0.00867, +0.01665]` | `[-0.01862, +0.00359]` | `+0.38702` |

**FACT:** A leave-query-out, high-data target-specific source oracle has
positive rank-loss headroom at all k. Its component-bootstrap lower bounds are
`+0.02647`, `+0.06979`, and `+0.07097` for k=1/3/5.

**INFERENCE:** Query ranking is improvable in principle, and the diagnostic can
recover a strong injected label channel. What remains unestablished is that the
current passive support labels provide stable, assignment-specific incremental
information. The v2 result therefore remains `NO_GO_INFORMATION_NOT_ADMITTED`.

### Root-cause ranking

1. **INFERENCE - Task and episode misdefinition.** Support/query chemistry and
   assay context do not reliably instantiate a transferable local relation.
2. **INFERENCE - Adaptation-object error.** A continuous residual field is too
   rich for at most `k-1` within-support contrasts.
3. **FACT/INFERENCE - Evaluation contamination in the old splits.** Previously
   inspected validation and holdout roles are development evidence only.
4. **INFERENCE - Shortcut susceptibility.** Earlier arm classifiers recovered
   CMAL arm identity from chemistry; matched counterfactuals remain mandatory.
5. **INFERENCE - Operator misspecification.** CMAL makes unrestricted continuous
   perturbations without structural abstention.
6. **INFERENCE - Objective mismatch.** MSE and all-pairs ranking do not directly
   optimize top-tail utility, but all major metrics worsened, so this is not the
   primary explanation.
7. **INFERENCE - Numerical optimization is not primary.** Adapter gradients and
   parameter changes were present.

## Part 3 - Identifiability Boundary

For a local residual model

\[
r_t(x)=b_t+\phi(x)^\top z_t+\epsilon,
\]

centering the support observations gives

\[
H_k e_S = H_k\Phi_S z_t + H_k\epsilon,
\qquad \operatorname{rank}(H_k\Phi_S)\leq k-1.
\]

**FACT:** k=1 has zero within-support contrasts; k=3 has at most two; k=5 has
at most four. Assay nuisance may consume additional dimensions.

**INFERENCE:** k=1 may identify an offset, sign, magnitude, or one diagnostic
response mode, but cannot test compound-label assignment. k=3/5 may identify a
very low-cardinality state only when support chemistry instantiates a
source-learned relation. None can nonparametrically identify scaffold-cold SAR.

## Part 4 - Surviving Hypotheses

### 1. Evidence-gated discrete rank intervention

**HYPOTHESIS:** Learn a support-conditioned `STOP` or bounded pair-order edit of
the frozen ranking, with zero/null evidence forcing `STOP`.

**Reviewer verdict:** `REJECT_AS_CORE_INNOVATION` at present. It is a useful
mechanism probe, but is currently a constrained support-conditioned meta-LTR.
It cannot be implemented unless balanced k=3/5 label and assignment gates pass
and a B=1/2 action oracle has practical headroom. Adjacent-position swaps imply
a fixed-library transductive claim and are not distractor/subset invariant.

### 2. Evidence-activated finite SAR grammar

**HYPOTHESIS:** In same-assay, MMP-connected episodes, infer one of
`INCREASE/DECREASE/NULL/UNSUPPORTED` for a source-frozen transformation rule and
execute at most one exact depth-1 rule.

**Reviewer verdict:** `CONDITIONAL_NOT_ADMITTED`. This is the more biologically
specific direction, but it must first pass a label-free coverage/power census
and a local k=3/5 assignment-information gate. It must beat copy-support-sign,
Matsy/MMP/MMS, hierarchical mixed effects, graph-difference, categorical
empirical Bayes, and a learned kernel. k=1 must structurally abstain.

### Rejected campaign direction

**INFERENCE:** ChEMBL publication order cannot identify a campaign-state or
selection-outcome mechanism. RetroDMTA can support retrospective prioritization,
but causal policy claims require project IDs, contemporaneous candidate pools,
selection actions, assay/batch metadata, and preferably randomized exploration
or a defensible sequential ignorability design.

## Part 5 - Decision And Exact Next Experiment

1. Freeze CMAL and do not promote any current code to `model/`.
2. Do not implement STOP/SWAP. The global v2 information gate failed.
3. Run one label-free same-assay/MMP coverage and component-power census using
   the v2 retained rows. Predeclare the MMP vocabulary, context radius, minimum
   cross-component frequency, support construction, and MDE before labels.
4. If coverage is sufficient, run one local k=3/5 G0/G1 gate with a true
   chemistry-fixed assignment derangement and component bootstrap.
5. If the local gate fails while synthetic controls pass, abandon the current
   A2S episode construction. Do not add adapter capacity.
6. If the local gate passes, compare a finite grammar first with classical
   MMP/Matsy and mixed-effects baselines. A learned model is admitted only if
   the finite rule state, protein routing, and abstention are load-bearing.
7. Keep locked-source labels sealed until the local protocol, MDE, seeds,
   baselines, controls, and stop rules are frozen. Recipient labels remain
   sealed until a one-time locked-source evaluation passes.

**FACT:** There is no key positive mechanism breakthrough and no authorization
to commit, push, or publish the current `model/` folder. The most valuable
result so far is a balanced negative information-admission result and a much
narrower falsifiable next question.

