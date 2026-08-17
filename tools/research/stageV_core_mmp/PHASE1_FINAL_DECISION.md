# Phase 1 final decision — transferable protein-conditioned interaction signal

Date: 2026-08-17. Machine-readable authority:
`PHASE1_FINAL_DECISION.json`. This document closes the tested mechanisms; it
does **not** reopen, retrain or relabel any of them, and it moves no threshold.

## Question

With protein sequence-derived features and ligand 2D features only, can a
model learn a protein-conditioned local interaction signal that transfers to
unseen protein components and aligns with affinity-difference / SAR ordering?

## Verdict

**BOUNDED NEGATIVE under the current BindingDB-Ki double-cold protocol.**
The strongest admissible estimator — the core/context-strictly-matched MMP
crossed double difference — is **not estimable** on this corpus: the primary
protein-component-cold repeated-key surface does not exist (32 rows / 4
components; 0 exact keys shared with the development-validation split), and
the target × transformation interaction variance is **not identifiable above
the defensible noise envelope**.

## Evidence chain, in dependency order

1. **Stage S** (`stageS_sar_field/REPORT.md`): whole-molecule pair global FiLM
   field. Correct protein +0.0065 Pearson over ligand-only, unresolved; a
   shuffled protein reproduced the gain; the protein path moved predictions by
   62% of label spread with +0.093 alignment to truth. Rejected.
2. **Stage P**: centered protein-counterfactual objective produced a
   reproducible but truth-unaligned protein response (+0.022 alignment).
   Rejected as the Stage S failure mode in another training objective.
3. **Stage T** (`stageT_mmp/`): first true-MMP run used a **core-blind** key.
   The pooled-protein discriminator failed; the forensic correction
   (`CORRECTION_20260817_CORE_KEY.md`) measured 40.4% of training D rows with
   disjoint core sets, median across-core nuisance 0.269 pK, and **withdrew**
   Stage T's global closure claim.
4. **Stage U** (`stageU_mmp_interaction/`): core-inclusive key, frozen before
   any U0 statistic. U0 admission failed on degree concentration: one target in
   one component = 29.63% of primary observations vs the 25% cap. No U1/U2.
5. **Stage V** (`stageV_core_mmp/`): corrected successor, all Stage U
   thresholds inherited verbatim plus the missing controls. V0 fails on the
   same concentration caps; V0b: primary `internal_repeated` = **32 rows / 4
   components / EIU 4** (<100 not evaluable), internal rich keys = 0; V1:
   `MS_effect` 0.452 vs preregistered noise 0.858,
   `theta = -0.406 [-0.704, -0.073]` (resolved negative).
6. **Direct pair-level noise audit** (`PAIR_LEVEL_NOISE_AUDIT.json`): the
   preregistered cell-level reference is likely inflated. Direct repeated
   same-panel MMP deltas: 88/42,534 pairs, 40 zero-range curation duplicates;
   disagreeing-only variance **0.303 [0.200, 0.427]**. Against that
   conservative reference, V1 cross-component is **+0.391 [-0.327, +0.368]**
   — unresolved. Only the downward-biased all-group reference gives a positive
   lower bound. Therefore the interaction variance is **not identifiable above
   the defensible noise envelope**.
7. **Development-validation structure census** (`METAVAL_STRUCTURE_CENSUS.json`):
   7,209 same-panel MMP observations / 2,757 potential D rows / 19 components,
   but **0 exact keys shared with meta_train**. The double-cold split makes the
   repeated-key primary surface unsuppliable on both development splits.

## What this decision does and does not close

**Closed as tested and negative / not estimable:**

* whole-molecule global FiLM protein SAR fields;
* centered correct-vs-wrong protein training on the incumbent trunk;
* core-blind MMP pooled-protein discriminator;
* core-inclusive exact-MMP protein × transformation interaction on this corpus
  and protocol, at the measurement/identifiability gates, without any neural
  operator trained because the frozen gates forbid it.

**Not closed, and not claimed:**

* protein-conditioned local interaction as a biological possibility;
* other datasets (Davis/KIBA remain promotion-gated);
* MSA/coevolution priors (no governed UniRef snapshot exists locally; the
  route remains blocked on an external asset, not falsified);
* looser-but-principled transformation equivalence classes (not tested; they
  would need their own cancellation analysis and could only be a screen, not
  the strict positive gate);
* any architecture that was never instantiated.

## Verification

* Stage V preregistration SHA-256
  `c567f66066c301fefe293048a4643fe4f65158077c3540ce1bbb0beb5d5844d4`.
* Stage U + Stage V research suites: **55 passed** (`RUN_SLOW=1`).
* Maintained suite `python main.py verify tests`: **310 passed / 6 skipped**.
* `meta_test` evaluated: **0**. Sealed confirmation split was never mounted.
* No neural model, checkpoint or production change was produced for this
  decision.
