# External Research Prompt: A2S-DTA After The Balanced V2 Information Gate

## Attachments to deliver

Preferred six-file package:

1. `A2S_POST_REVIEW_V2_GATE_DECISION_2026-08-01.md`
2. `a2s_source_information_gate_v2_2026-08-01.json`
3. `a2s_source_information_gate_lock_v2_2026-08-01.json`
4. `a2s_information_gate.py`
5. `a2s_source_lock_v2.py`
6. `A2S_META_ADAPTATION_RESEARCH_HANDOFF_2026-08-01.md`

If the chat accepts only three attachments, send items 1-3. If it accepts four,
add item 4. Do not attach recipient outcome files or the old `model/` folder.

## Prompt to paste into the chat

```text
You are an independent principal-investigator committee reviewing an
abundant-to-scarce drug-target affinity (A2S-DTA) project. Do not assume that
CMAL, CAMP, ECA, STOP/SWAP, finite SAR grammar, or any prior proposal is
correct. Do not write code. The objective is to decide the one next
source-only falsification experiment after a balanced information gate failed.

SCIENTIFIC OBJECTIVE

The desired contribution is a genuinely learned transferable adaptation
mechanism. It must learn from abundant source targets how to use k={1,3,5}
correctly assigned measurements from a strictly unseen target to improve
target-specific, query-dependent compound ranking over the identical frozen
support-free base. Calibration, retrieval, interpolation, kernel regression,
closed-form Bayesian updates, ordinary fine-tuning, and a larger DTA backbone
are baselines, not the final contribution.

EVIDENCE DISCIPLINE

Tag every substantive statement as FACT, INFERENCE, or HYPOTHESIS. Do not claim
information-theoretic impossibility from a finite negative probe. Do not treat
correct-support > wrong-support as success unless correct support also improves
over the identical frozen base. Do not use or request recipient labels.

ESTABLISHED FACTS

1. Current ChEMBL episodes use publication/document ordering, not verified
   project IDs or DMTA cycles. About 81% of recipient queries are scaffold-cold,
   support/query scaffold overlap is about 8.7%, nearest support-query ECFP
   Tanimoto is about 0.223, and exact assay-context overlap is rare.

2. CMAL is mechanically active but failed scientifically. Source validation
   changed by CI -0.01011, Spearman -0.03235, NDCG@10 -0.01551. The repeatedly
   inspected source holdout changed by CI +0.00303, Spearman +0.00898,
   NDCG@10 -0.00122, while RMSE worsened. No recipient labels were opened.

3. The old v1 information gate was non-confirmatory because one provenance
   component occupied 97.0% of OOF held-out rows.

4. The attached v2 label-free lock assigns homology components to roles first,
   then quarantines rows carrying document or assay tokens shared across roles.
   It retained 68,782/185,591 source pKi rows. Retained fit/probe/locked roles
   contain 222/110/107 homology components and 36,609/16,068/16,105 rows.
   Target, homology, document, and assay overlap is zero. Five fit OOF folds
   contain 7,322/7,322/7,322/7,322/7,321 rows.

5. G0 and G1 are separately fitted equal-capacity ridge probes. Their only
   difference is that five support-label features are masked in both training
   and evaluation for G0. G1 uses correctly assigned residual features. For
   k>=3, the assignment control cyclically deranges residuals among the same
   support compounds. k=1 has no assignment permutation.

6. The balanced v2 development result remains
   NO_GO_INFORMATION_NOT_ADMITTED. Component-bootstrap 95% intervals for
   Delta_label rank loss are:
   k1 [-0.00426, +0.01271],
   k3 [-0.00865, +0.01185],
   k5 [-0.00867, +0.01665].
   Delta_assign intervals are:
   k3 [-0.00087, +0.01366],
   k5 [-0.01862, +0.00359].
   No required lower bound is positive.

7. Synthetic label-channel controls pass strongly at k1/k3/k5 with rank-loss
   gaps +0.31635/+0.22976/+0.38702.

8. A leave-query-out high-data target-specific source oracle has positive
   headroom. Rank-loss 95% lower bounds are +0.02647/+0.06979/+0.07097 for
   k1/k3/k5. Therefore query ranking is improvable, but the current passive
   support labels have not shown stable incremental information.

9. Locked-source and recipient labels remain sealed. The v2 probe role is
   development-only. No current mechanism is admitted, and no model code may be
   promoted or published.

IDENTIFIABILITY BOUNDARY

For centered support residuals, rank(H_k Phi_S) <= k-1. Thus k=1 provides no
within-support contrast, k=3 at most two contrasts, and k=5 at most four before
assay nuisance. An arbitrary target-specific residual function is not
identifiable. A surviving object must be finite, low-cardinality, auditable,
and structurally able to abstain.

TWO SURVIVING HYPOTHESES

A. Evidence-gated discrete rank intervention: support evidence triggers STOP
or at most B=1/2 bounded rank edits. Current verdict: useful probe but rejected
as core innovation because it is close to support-conditioned meta-learning to
rank. It cannot be implemented before a positive label/assignment gate and a
positive legal-action oracle.

B. Evidence-activated finite SAR grammar: in same-assay MMP-connected episodes,
support residual contrasts activate one of
INCREASE/DECREASE/NULL/UNSUPPORTED for a fixed source-derived transformation;
the exact executor applies at most one depth-1 rule. Current verdict:
conditional and not admitted. k=1 must abstain.

Campaign-state adaptation is rejected on current ChEMBL because publication
order does not identify project decisions or candidate pools.

YOUR TASK

Design exactly one next source-only falsification experiment. The default
candidate is a label-free same-assay/MMP coverage and component-power census,
followed only if feasible by a local k=3/5 label and assignment information
gate. Challenge this default if the evidence justifies doing so.

Your response must:

1. Audit the v2 lock and probe line by line for remaining leakage, selection,
   underfitting, or invalid comparison. Pay special attention to row quarantine,
   the change in estimand, separately fitted G0/G1 preprocessing, OOF residuals,
   target-balanced weighting, component bootstrap, and the leave-query-out
   high-data oracle.

2. Decide whether the balanced global null is strong enough to abandon the
   global passive episode construction now, or whether one predeclared local
   stratum is still justified.

3. If same-assay/MMP is justified, define a precise label-free MMP coverage
   census: canonicalization, single-cut rule, attachment environment, maximum
   substituent size, stereochemistry, assay-token semantics, duplicate handling,
   minimum cross-target and cross-component rule frequency, and the statistical
   unit. State how support/query rows are constructed without outcome selection.

4. Define the minimum detectable effect and minimum component count needed
   before opening any local source labels. Do not invent a universal threshold;
   tie it to component power, assay uncertainty, and ranking utility.

5. If coverage passes, specify one local k=3/5 G0/G1 test. G0 and G1 must be
   nested equal-capacity models; complete encoder/head residuals must be OOF at
   the homology/provenance component level; G1 must use correct assignment;
   the control must derange residuals among the identical support compounds.

6. Include synthetic effect-dose controls, label-noise dose response,
   copy-support-sign, Matsy/MMP/MMS, hierarchical mixed effects,
   graph-difference, categorical empirical Bayes, learned kernel, frozen base,
   ridge, KRR, MAML, ANIL, MetaDTA, current CMAL, CNP/MetaFun, and an
   equal-capacity generic reranker where scientifically applicable.

7. Define leakage controls: protein shuffle, target shuffle, assay-matched wrong
   support, chemistry-fixed assignment permutation, residual-null, k=1 sign and
   norm controls, query permutation, distractor insertion, query subset, and
   library-size shift.

8. Give explicit STOP conditions for insufficient MMP coverage, insufficient
   power, failed local Delta_label, failed Delta_assign, absent protein
   contribution, equality with classical MMP/Matsy baselines, and instability
   under candidate-set perturbations.

9. State the exact evidence required before implementing either surviving
   mechanism, before opening the locked-source role once, and before opening
   recipient labels once.

10. Conduct a hostile novelty review against AdaMBind, MetaDTA, FS-CAP,
    CNP/ANP, MetaFun, meta-learning to rank, PiRank/selective prediction,
    Matsy, MMP/MMS, graph-difference methods, hierarchical mixed effects, and
    learned-kernel few-shot methods. Do not claim novelty merely because an
    output is discrete or a neural update is non-closed-form.

OUTPUT FORMAT

Part 1 - Critical audit findings, ordered by severity.
Part 2 - Decision: STOP_GLOBAL_NOW, RUN_ONE_LOCAL_MMP_GATE, or REQUIRE_NEW_DATA.
Part 3 - Exact label-free coverage/power protocol, if admitted.
Part 4 - Exact local source-only gate, if coverage passes.
Part 5 - Hostile novelty verdict for both surviving hypotheses.
Part 6 - A binary decision tree with all kill criteria.

Do not write code. Do not propose broad hyperparameter sweeps. Do not force a
positive mechanism. A justified rejection of the current episode construction
is a valid and potentially more valuable result.
```

