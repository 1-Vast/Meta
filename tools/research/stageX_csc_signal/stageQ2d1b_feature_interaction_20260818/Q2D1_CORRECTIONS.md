# Q2d-1 corrections (2026-08-18) — measurement-chain forensics

Old stage: stageQ2d_bilinear_qualification_20260818 (read-only, unmodified).
Evidence: tests/test_q2d1_forensics.py in stageQ2d1b_feature_interaction_20260818
(5 passed; each pins one defect against the old source and artifacts).

## Verified defects

1. PHASE D/E NEVER ENABLED CENSORING. q2.generate defaults to
   censoring='noclamp'; the Q2d-1 ladder called generate without the argument
   and Phase D reused that lat via a shallow copy. determinante=True on every
   cell, n_censored=0, all bounds zero. (The X0c full grid has the same
   property; only its floor-clamp negative control passed censoring='floor_clamp'.)
2. PHASE C WAS 70% OBSERVED (frozen MCAR seed) - confirmed as described; the
   defect list records it for completeness.
3. CLOSED-FORM 'TRAIN HOLDOUT' WAS INSIDE THE FIT. The diagnostic fitted the
   rank-4 SVD on the FULL train submatrix and then evaluated on a random
   subset of those same cells; the reported dz 0.951 / sp 0.927 is an in-fit
   reconstruction, not a holdout bound.
4. HALF-COLD SVD FAILURE IS BY CONSTRUCTION. Val ligands are unseen in train,
   their columns are all-NaN in the fit matrix, they are filled with a
   constant, and the rank-4 reconstruction therefore has zero variation
   across unseen-ligand columns (asserted: std < 1e-9). The dz 0.50 on that
   surface measures the diagnostic's design, not an information limit.
5. PHASE A WAS NOT INTERACTION-ONLY. Truth included mu + random pmain/lmain
   with MAIN_SD=1.0 (same scale as tau*=1.0); the learner carried per-row
   and per-ligand ID biases; checkpoint selection ran on val cells whose
   ligands are unseen in train.

## Withdrawn claims (from Q2D1_REPORT.md; the old report file itself is
## left unmodified as a historical artifact)

- WITHDRAWN: 'the learner and optimization are NOT the failure' - the 0.95
  bound was in-fit, and Phase A's confounds mean exact_bilinear 0.493 does
  not isolate the optimizer.
- WITHDRAWN: 'sigmoid, missingness, censoring and main-effect competition are
  each ruled out as the killer step' - censoring was never enabled in D/E,
  so the censoring step of the ladder was never executed.
- WITHDRAWN (weakened): 'unseen-ligand recovery is information-theoretically
  impossible for any learner' - the half-cold failure was built into the SVD
  diagnostic, and the parameter-count argument bounds only the specific
  random-factor truth on the fixed ECFP4 map, not transferable feature-
  conditioned mechanisms.

## Retained (unchanged, still supported by old artifacts)

- exact_bilinear failed the frozen gate on every ladder level under the old
  (confounded) protocol; per-level medians remain as recorded.
- Negative arms failed in every phase (additive_only/no_interaction_head
  dz 0.000; ligand_only/shuffled/random at chance).
- oracle_latent 0.662-0.699 under the old protocol; it did not clear 0.70.
- The old ID-random-factor generator cannot support ligand-cold transfer -
  retained as a statement about that generator, not about the mechanism.

## What the corrected successor (this stage) changes

Phase A becomes interaction-only (y = tau*I(P,L) + noise, no mu/pm/lm, no
learner ID biases); truth factors become feature-conditioned (per-position
pocket physicochemical z-scales x compact ligand substructure descriptors);
censoring is passed explicitly with a censored_count > 0 assertion; the
closed-form oracle fits on train only and reports protein-cold, ligand-cold
and double-cold surfaces separately; checkpoint rule is best-train-loss
(no unseen-ligand selection); each ladder level must reproduce the previous
level value-exactly with the added ingredient disabled.
