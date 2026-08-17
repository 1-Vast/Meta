# Final state of the MetaSieve research programme (2026-08-18)

This document is the closing summary of the long-run objective recorded in
tools/research/GOAL_ACTIVE.md. It ties together the three authorities:
BOUNDARY_20260817_NIGHT.md (the conclusion), EVIDENCE_LEDGER.md (the
falsification record) and tools/research/stageN_audit/AUDIT_REPORT.md
(the verification).

## Terminal condition reached

The objective had two terminal conditions: a reproducible performance leap,
or a scope-bounded final conclusion after multi-family falsification.
**The second condition is met** and verified end-to-end:

- The k=0 MSE <= 1.00 pK^2 target is not reachable under the governed
  BindingDB-Ki double-cold protocol with the locally available legal inputs:
  the level term is assay-history-dominated (within-document transfer
  R^2 +0.451, zero across documents by the split's document closure), the
  legal transferring inputs cover at most ~26% of level variance, the best
  trained zero-shot level^2 on record is 1.2151 (Stage L, ranking-degraded)
  against a 0.1239 budget, and per-seed k=0 MSE never went below 2.10.
- k>=1 is a different story: k=5 sits at 0.939-1.007 across seeds with
  honest controls, and K-REG produced the first all-k resolved MSE
  improvement across three seeds — but its shape gain did not survive
  pooling, so nothing was promoted.
- No candidate passed the full promotion gates; meta_test was never opened
  (logical exclusion after parsing; 768 cells withheld; 0 evaluations in
  104 audited artifacts); nothing moved to model/ or scripts/.

## Falsification ledger (all leak-free, preregistered stages)

Frameworks: analytic/legacy operators; BPSF/CIPF; contact-grammar;
moment-form and inner/outer-loop meta-learning (A/B); centered-objective
protein conditioning (P); panel-set level head (E); pairwise learned
transport (F); ESM-650M input lane (G/G2); live ESM-150M LoRA (I);
assay-aware level head (J); support-gated level head (L); contrastive
coembedding (K/K2). Training methods: ranking-loss substitutions (R9-R14),
orthogonal level/shape routing, paired level alignment, episodic InfoNCE,
coembedding regression alignment, LoRA conditioning.

External representations: ESM-150M/650M frozen; ESM-150M LoRA-tuned;
ChemBERTa-77M ligands; structure/pocket priors (209/387 coverage);
panel composition; assay covariates (endpoint, counts, journal/publisher,
documents); protein function annotations (ProteinKG25 GO bags, 313/387
matched, Stage P0). Every locally available legal family measured; none
breaks the k=0 level wall.

## Remaining open lanes (all externally blocked, not untested by choice)

- MSA / conservation: requires a governed UniRef snapshot; none exists
  locally and network acquisition is not performed in this environment.
- Larger structure coverage: 209/387 targets have a homologous holo
  structure locally; a fuller structure lane needs external data.
- Davis/KIBA independent replication: authorized only after a candidate
  passes the promotion gates; none did, so it was not run (the datasets
  remain sealed under dataset/sealed/).

## Governance state

- meta_test: sealed (logical exclusion after parsing), fail-closed,
  0 evaluations, written authorization required and never issued.
- Query labels: loss-and-metric-only in every stage.
- No closed-form solvers, no cross-dataset support, no transductive
  calibration, no multi-stage pretraining disguised as one run.
- All 7 training stages preregistered before results existed; single-seed
  screens gated by stop rules; the one promising lane went through a
  preregistered three-seed confirmation.
- Verification at close: maintained suite via the sanctioned entrypoint
  (python main.py verify tests) 268 passed / 6 skipped; research suites with
  the corpus slow tests enabled (RUN_SLOW=1) 147 passed; the final boundary
  audit (tools/research/stageN_audit/) re-derives every load-bearing number
  bitwise from the raw evaluation rows.

