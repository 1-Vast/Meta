# Final state of the MetaSieve research programme (2026-08-18)

This document is the closing summary of the long-run objective recorded in
tools/research/GOAL_ACTIVE.md. It ties together the three authorities:
BOUNDARY_20260817_NIGHT.md (the conclusion), EVIDENCE_LEDGER.md (the
falsification record) and tools/research/stageN_audit/AUDIT_REPORT.md
(the verification).

**Interpretation authority: `POST_COMPLETION_REVIEW_20260818.md`.** Every
statement below is scoped to **BindingDB-Ki double-cold development
evidence**. Nothing here generalizes to other DTA datasets, to untested
architectures, or to what is knowable in principle: the empirical failures
recorded in this cycle are failures of the candidates that were run, not
information-theoretic bounds.

## Terminal condition reached

The objective had two terminal conditions: a reproducible performance leap,
or a scope-bounded final conclusion after multi-family falsification.
**The second condition is met** and verified end-to-end:

- **No tested candidate reached the k=0 MSE <= 1.00 pK^2 target** under the
  governed BindingDB-Ki double-cold development protocol with the locally
  available legal inputs: the level term is assay-history-dominated
  (within-document transfer R^2 +0.451, zero across documents by the split's
  document closure), the tested governed probes explain up to 25.9% of level
  variance, the best trained zero-shot level^2 on record is 1.2151 (Stage L,
  ranking-degraded) against a 0.1239 budget, and per-seed k=0 MSE never went
  below 2.10.
- **The target is arithmetically possible.** k=0 MSE = level^2 + centered, and
  the measured centered term is 0.8648, so an oracle level predictor would put
  k=0 at about 0.865 — below 1.00. The 0.1239 budget assumes today's centered
  error is fixed. What no tested model did was move both terms in the same
  direction at once. Nothing measured here is an information-theoretic bound.
- k>=1 is a different story: k=5 sits at 0.939-1.007 across seeds with
  honest controls, and K-REG produced the first all-k resolved MSE
  improvement across three seeds — but its shape gain did not survive
  pooling, so nothing was promoted.
- No candidate passed the full promotion gates; no sealed meta_test label
  entered any fitting, selection or reported metric (logical exclusion after
  parsing; 768 cells withheld; 0 evaluations in 106 audited artifacts);
  nothing moved to model/ or scripts/.

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

## Governance state (after the 2026-08-18 post-completion repair)

- meta_test, **two surfaces, stated separately**:
  - every recorded artifact was produced on the all-label corpus, giving
    **logical exclusion after parsing**, fail-closed, written authorization
    required and never issued, **0 evaluations in 106 audited artifacts**;
  - a **physically isolated** surface now exists and is mountable
    (`scripts/build_governed_split_views.py` ->
    `dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views`,
    `QPSMPData(split_view=...)`, `--split-view` on the maintained trainer and
    evaluator). The meta_test label artifact lives out of tree; a development
    process has no filesystem path to it, proved by a file-access spy in
    `tools/tests/test_physical_meta_test_seal.py`. The isolated surface is
    element-for-element identical to the corpus construction, so it changes no
    recorded number — and no recorded artifact was produced on it.
- Checkpoint selection: the maintained trainer now defaults to
  `--selection internal` (fit/internal-validation partition of meta_train by
  homology component, `scripts/internal_validation.py`), so meta_val is not
  read during training. The legacy meta_val rule survives only as the
  disclosed diagnostic Stage B used to measure it at ~0.62 pK^2 optimism at
  k=0. Every figure in this document predates that change and therefore
  carries the meta_val-selection optimism.
- Query labels: loss-and-metric-only in every stage.
- No closed-form solvers, no cross-dataset support, no transductive
  calibration, no multi-stage pretraining disguised as one run.
- **11 retained trained stages**, discovered from the filesystem rather than
  from a hand-maintained list, all preregistered; one mtime ordering
  exception (stageI_lm) is disclosed in `stageN_audit/AUDIT_REPORT.md`. The
  earlier "7" and "8" counts were two stale hard-coded lists, both of which
  also omitted Stages A, B and P_cpc.
- Verification at close (2026-08-17T04:24Z / 04:26Z, conda env `drug`,
  Python 3.11.15, torch 2.6.0+cu124, CUDA available):
  - maintained suite, sanctioned entrypoint `python main.py verify tests`:
    **310 passed / 6 skipped**, exit 0;
  - complete research suite, `RUN_SLOW=1 pytest tools/research -q`:
    **255 passed / 2 skipped**, exit 0;
  - the final boundary audit (`tools/research/stageN_audit/final_audit.py`)
    re-derives every load-bearing number bitwise from the raw evaluation rows
    and regenerates `FINAL_BOUNDARY_AUDIT.json` and `AUDIT_REPORT.md`.
  - Historical counts reconciled, neither deleted: **147** was
    `RUN_SLOW=1 pytest tools/research/stageA_innerloop
    tools/research/stageB_complementary` (147 collected, all passing); **151**
    was that same pair plus `tools/tests/test_research_record.py` **without**
    `RUN_SLOW` (151 passed / 12 skipped); **135** was the pair alone without
    `RUN_SLOW`. All three were subsets; none was the complete research suite.

