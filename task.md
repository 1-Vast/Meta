# Current task contract

Updated: 2026-08-18 (night). Status: **final state — bounded conclusion
established across the full mechanism and covariate space.** Stage L
(support-gated assay-aware level head) closed the last composition: the gate
preserved k>=1 MSE but ordering degraded with resolved intervals (k=2/3/5
Spearman), because the zero-shot level objective and within-target ordering
conflict on the shared trunk. Record landmarks: K-REG = first all-k resolved
MSE improvement across 3 seeds (not confirmed on centered); L = best k=0
calibration ever (MSE 2.0997, level^2 1.2151) but ranking-degraded. No
candidate passed all promotion gates; meta_test stays sealed; nothing moved
to model/ or scripts/. Stage M0 (ChemBERTa-77M ligand embeddings,
`tools/research/stageM_chemberta/`) closed the last locally testable legal
input family (ordering r +0.147 below occupancy; level probe = grand mean).
Authority: `report/BOUNDARY_20260817_NIGHT.md` (final),
`tools/research/stageL_gated/REPORT.md`,
`tools/research/stageM_chemberta/REPORT.md`. Closing summary:
`report/FINAL_STATE_20260818.md`; verification:
`tools/research/stageN_audit/AUDIT_REPORT.md`.

## Objective

Produce a trainable, medium-scale model with excellent zero-shot and k=1/2/3/5
cold-target DTA performance. At least one central innovation must reside in the
training method, and every claimed few-shot gain must depend on correct support
labels and the recipient protein rather than scalar calibration or ligand recall.

## Immutable rules

- Data/evaluation: use the governed BindingDB Ki double-cold protocol. An episode
  contains one target; support/query ligands are unique. Keep current meta_test
  sealed until a frozen candidate passes all development gates. The seal is
  logical exclusion after parsing plus written authorization, not physical
  isolation; see tools/research/a2_readiness_v2/SPLIT_ISOLATION_SPEC.md.
- Learning: ordinary end-to-end forward/backward training. **Inner/outer loops
  and differentiable support adaptation are permitted as of 2026-08-17** (user
  instruction; supersedes the previous blanket prohibition on inner loops and
  deployment-time support adaptation). Still prohibited and unchanged: ridge
  regression, analytic solvers, pseudoinverses, closed-form shortcuts,
  query-label adaptation at inference, and multi-stage pretrain/finetune
  regimes disguised as one run. Adaptation at inference may read support
  **inputs and labels**; query labels remain loss-and-metric-only at every k.
  A model may jointly consume multiple governed task types, but training stays
  one ordinary single-stage optimization process.
- Datasets: one public supervised DTA source per experiment. The governed
  BindingDB-Ki double-cold protocol is the only one currently authorized. No
  merging with Davis, KIBA, ChEMBL or structural corpora, and no cross-dataset
  support, retrieval, normalization statistics, labels or checkpoints. If a
  candidate passes here, Davis and KIBA must be trained independently from
  scratch in separate experiments.
- Information: query labels are loss-only. No test labels, component leakage,
  donor-label leakage or query-panel transductive centering.
- Geometry: no atomic protein-ligand 3D claim without a legal common-frame pose.
  Current DTA coverage is 0/17,717, so Cartesian code is not a performance source.
- Evidence: smoke tests only find bugs. Promotion needs matched budgets, multiple
  fixed seeds, nested k, component-paired bootstrap and clean counterfactuals.
- Scope: one model mechanism and one training mechanism at a time. Stop when a
  preregistered family gate fails; do not rescue it by unplanned complexity.

## Established baseline and boundary

- Leak-free references (2026-08-17): T2 (incumbent recipe retrain, internal
  checkpoint selection) k=0 2.5961 / k=1 1.7712 / k=2 1.3245 / k=3 1.2197 /
  k=5 0.9859, three-seed k=0 band 2.458-2.981 and k=5 band 0.946-1.007;
  ESM-650M lane G k=0 2.239-2.790, k=5 0.944-0.987 (not confirmed).
- The k=0 <= 1.00 target is protocol-conditioned: level is assay-history
  dominated (within-document transfer R^2 +0.451; 0% across documents), the
  legal transferring inputs cover <=26% of level variance, and the best
  trained level^2 (1.52) alone exceeds the whole 1.00 budget. Full ledger:
  `report/BOUNDARY_20260817_NIGHT.md`.
- A0 is the retained incumbent and zero-shot ordering reference.
- B3 and C2 join A0 on the k0 MSE/CI Pareto frontier; none dominates.
- Fixed Morgan/Tanimoto residual weighting is the comparator to beat at k>=2.
- R7-R14 closed free query gates, direct shape heads and ranking-loss substitution.
- The 0.782 cliff-sign record is meta_val development evidence on dominated C1,
  not a confirmation or overall-performance claim.

See `report/CURRENT_MODEL_EVIDENCE.md` and `report/BOUNDARY_20260816.md`.

## Active plan

`report/NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md` is **closed**. Its exact
operator was implemented and run on real episodes
(`tools/research/a2_exact_probe/FINAL_DECISION.md`): 19 structural gates pass,
then the performance gates fail. Authority: `A2_EXACT_meta_val.json`.

*Corrected 2026-08-16 — the earlier "both falsification controls fail
inverted" wording here described the superseded pre-repair artifact and is
withdrawn.* Against the repaired deterministic nested banks, `beats_tanimoto`
fails with resolved intervals at all four k (-0.042 / -0.173 / -0.217 /
-0.252), while `degrades_under_wrong_protein` (+0.040 to +0.008) and
`degrades_under_label_shuffle` (-0.146 to -0.129) are **unresolved** — every
interval crosses zero. The operator fails its controls **by being inert, not by
being inverted**: query spread at k=5 is 0.0027 pK against Tanimoto's 0.2865
and a random feature's 0.4521. `SUPERSEDED_A2_EXACT_meta_val.json` retains the
pre-repair numbers as evidence only; it is not a comparator.

Stage P has run: two matched arms, three seeds, 1,200 steps
(`tools/research/stageP_cpc/`). P1 fails at -0.0066 [-0.0545, +0.0417]. Correct
and wrong protein give identical within-target ordering at every k in both
arms. The centered objective excluded the level branch exactly as designed
(gradient 8.1e-07) and made the protein response reproducible across seeds
(+0.316) but unaligned with truth (+0.022). The stop rule was applied; the
admission stage did not run.

## Stage A: target-task inner/outer-loop meta-learning (opened 2026-08-17)

`tools/research/stageA_innerloop/`. AdaMBind-**inspired** framework only; no
reproduction claim is made or permitted. Single seed, development evidence,
directional screening. Preregistration frozen before any arm trained.

Three matched arms at 1,200 steps, seed 20260815: `A0` no inner loop, `A1`
inner/outer loop with uniform task sampling, `A2` the same plus adaptive task
selection scored by post-adaptation query loss and support/query gradient
cosine. One code path serves all three, so `--inner-steps 0` reproduces the
accepted recipe bitwise and `A0` is matched by construction rather than by
inspection.

Adaptable scope: `interaction_head.2.{weight,bias}` — **97 parameters, 0.0054%
of the 1,798,833 trainable**. The smallest subset that can reorder ligands
rather than only shift their level, and the split is the instrument for the
k=1 shape-versus-level question. Inner step size 0.1 at 1 step, selected on
`meta_train` component folds against the frozen A0 checkpoint
(`INNER_LR_SELECTION.json`); `meta_val` was not read for it.

Stage 0 audit (`AUDIT_DATAFLOW.json`) confirmed the episode contract with 0
violations in 400 draws and the nested banks across k, and established one
disclosure that travels with every number: **the trainer selects checkpoints on
`meta_val` labels.** The rule is identical in all three arms, so it cannot
manufacture a between-arm difference, but it makes every reported `meta_val`
figure an optimistic development estimate rather than a held-out one.

### Result: NOT PROMISING on the conjunction (4 of 6 gates)

Authority: `tools/research/stageA_innerloop/RESULT.json` and `REPORT.md`.

- **`A1` (inner loop): weak positive, unresolved.** Better than `A0` at every k
  on every metric (k=0 MSE 2.0579 vs 2.0753; MSE gains 0.1111 / 0.0980 / 0.0443
  / 0.0206 at k=1/2/3/5) with **every interval crossing zero**. Only k=1 and k=2
  exceed the 0.058 retraining spread. Retained as a three-seed candidate.
- **`A2` (task selector): REJECTED.** Worse than `A1` on MSE, Pearson, Spearman
  and CI at every k. Cause identified, not noise: it selected candidates with
  support/query gradient cosine +0.9897 against a +0.6555 population mean — the
  tasks where support already predicts query. Gradient agreement measures
  redundancy, not informativeness. A learned bi-level selector is **not**
  authorized; the preregistration allowed it only on credible evidence.
- **Support labels are genuinely load-bearing** — every wrong-support
  counterfactual is resolved (permuted +0.36 to +0.43 pK² at k≥2; matched-wrong
  +0.81 to +2.05 across k).
- **`A1` produced the first resolved wrong-protein gap in the record**:
  +0.0188 [+0.0052, +0.0327] at k=2, +0.0177 [+0.0079, +0.0282] at k=3,
  +0.0085 [+0.0037, +0.0137] at k=5. Small, but every prior wrong-protein gap in
  R0-R14 and Stage P crossed zero.
- **k=1 is a level shift, not shape**: a bias-only update recovers 81% of the
  k=1 gain and the shape residual is unresolved.
- **The inner loop is free in encoder cost** (6,480 forwards for both `A0` and
  `A1`) because the 97-parameter scope is downstream of the encoder.
- **An inner loop cannot be bolted onto a trained model.** Measured
  `alpha = 2·lr·‖h‖²`: `A0` 1.514, overshooting on 100% of episodes, which
  predicts and matches its alternating sweep (1.5352 / 2.8626 / 1.4349 /
  2.2049); `A1` 0.241, stable. Training with the loop is what conditions it.

## Stage B: complementary meta-adaptation (2026-08-17) — REJECTED

`tools/research/stageB_complementary/`. Corrected Stage A's eight analysis
defects (`CORRECTION_AUDIT.md`, Stage A artifacts preserved), then ran four
matched arms — `T` transport-only, `M` meta-only, `H` naive hybrid, `C`
complementary — with **checkpoint selection that never reads `meta_val`**
(227 fit / 31 internal-validation components, partitioned by homology
component).

**Verdict: the AdaMBind-inspired framework is not admitted to production.**
Nothing promoted. Three stop conditions fired.

Two measurements matter beyond this stage:

1. **`meta_val` checkpoint selection is worth ~0.62 pK² at k=0.** `Tleak`
   (identical arm, same fit components, `meta_val` selection) scores 2.1246
   against the leak-free 2.7425, resolved at k=2/3/5. That is **93% of the gap**
   to Stage A's recorded 2.0753, **5.6× the largest mechanism effect** in the
   cycle and **10.7× the retraining spread**. Every `meta_val` number from the
   standard trainer, including the recorded incumbent band, is optimistic by
   about this margin.
2. **The ligand representation collapses within a target** — mean pairwise
   cosine 0.997 between a target's query-ligand readout activations, already
   present at `embed`. A weight update then moves every query by the same
   amount (a level shift), and against a mean-zero target it moves them by
   nothing. This is why `H`'s and `M`'s corrections are 99.7% level and why `C`
   is inert, and it links Stage R, Stage P and Stage L2 to one upstream cause.

`C` is the only arm on record to improve MSE and ranking together with resolved
intervals (Spearman +0.049 to +0.092), but **none of it comes from the
adapter**: the `C`−`T` ranking contrast is bitwise identical at k=0 and k=1,
where the meta term is exactly zero. It is a zero-shot trunk effect from the
training objective, not a few-shot mechanism.

## Stage C: the feasibility boundary (2026-08-17) — measurement only

`tools/research/stageC_level_shape/BOUNDARY.md`. No training. Baseline is the
leak-free Stage B `T` checkpoint.

**MSE decomposes exactly into `level² + centered`, and at k=0 it is 68% level.**
With a *perfect* level predictor every k lands below 1.00 (0.876 / 0.876 / 0.807
/ 0.798 / 0.734). Support labels are already the level mechanism — level² falls
1.87 → 0.28 across k while the shape term barely moves. Every "few-shot gain"
this project has measured is target-level calibration arriving through labels.

**But target level is not predictable from anything tested.** Against a
calibrated-constant reference (1.3471), the ESM linear probe is 4.85× worse, kNN
3.81×, sequence length 1.88×, the `meta_train` grand mean 1.61×, the full
trained model 1.27×, and the best ESM MLP probe 1.21×. The decay sweep selected
the largest value, driving the fit toward a constant. Since `centered ≥ 0`:

> **k=0 MSE ≥ 1.6357 with the best legitimate level predictor, even with perfect
> ordering.** Reaching 1.00 needs a 13.2× reduction in level error.

**The k=0 ≤ 1.00 target is therefore NOT reachable with the current inputs.**
The missing information is *zero-shot target-level affinity calibration for
unseen homology components*.

**Stage B's collapse claim was too strong.** Cosine 0.997 was insufficient
evidence. A frozen probe finds a **resolved** within-target ordering signal in
`occupancy` — meta_val r = **+0.2182 [+0.0751, +0.3670]** — the only
representation that survives out of component (`embed` +0.007, `readout_hidden`
+0.026, `section` +0.060, all unresolved). The model has this signal and does
not use it: `occupancy` reaches the endpoint only through `contact_weight`, a
`Linear(24→1)` of the same capacity as the probe, yet the endpoint orders no
better than a constant (ratio 1.0277). Worth ≈ 0.04 pK² if fully exploited —
real, resolved, and far too small to change the k=0 verdict.

## Next authorized work

The constraint is **information, not optimization**. No training schedule, no
adapter, no architecture growth. In priority order:

1. **external protein representation for target level** — run the preregistered
   M0/MSA lane, or structure-derived pocket descriptors, reported as external
   data, and test directly against the calibrated-constant reference (1.3471);
2. **assay/library covariate test** — regress target level on assay covariates
   within `meta_train`. If level is partly a property of testing history rather
   than the protein, the zero-shot target itself must be restated;
3. **the `occupancy` shape lever** — separate the level and shape paths so the
   24 contact-type parameters are trained on within-target centered supervision
   rather than level-dominated MSE. It can improve centered MSE, CI and
   Spearman; it cannot reach MSE ≤ 1.00 and must not be reported as if it could;
4. **remove `meta_val` checkpoint selection** from the production trainer and
   re-establish the incumbent band under a leak-free rule (Stage B measured the
   leak at ~0.62 pK² at k=0).

## Method-ladder cycle (opened 2026-08-16, CLOSED 2026-08-18)

The 2026-08-17/18 cycle superseded the ladder. Every one of the eight named
families now has a measured successor stage with a verdict — see
`tools/research/method_ladder/CLOSURE_MAP.md`. The shared M3 discriminator
harness (`tools/research/method_ladder/_shared/`, 25 structural tests) is
built and reusable; no family remains pending.

Eight named method families are tested sequentially under
`tools/research/method_ladder/<family>/`, each through the ladder
M0 primary-source reconstruction -> M1 input/identifiability gate ->
M2 structural/synthetic -> M3 frozen or low-cost discriminator ->
M4 one-seed training screen -> M5 three-seed admission. A failure stops that
family at its rung; no post-hoc rescue variants. Every family ends in exactly
one verdict: `REJECTED_BY_INPUT_CONTRACT`, `REJECTED_BY_STRUCTURAL_GATE`,
`REJECTED_BY_FROZEN_DISCRIMINATOR`, `REJECTED_BY_TRAINING_SCREEN`,
`ADMITTED_TO_FULL_EVALUATION`, or `ADMITTED_COMPONENT`.

Families: (1) multimodal representation collapse + basis reallocation;
(2) Gradient Blending / OGM; (3) Disentangled Gradient Learning;
(4) attention MIL / Set Transformer / adaptive pooling; (5) DrugBAN-style
bilinear interaction; (6) FS-CAP-style episodic scale; (7) AdaMBind-style task
valuation and label-noise robustness; (8) MMP-cliff transformation learning.

The standing measured lead into this cycle is the fusion/pooling localization:
Phase 3 put the loss of protein-differential at the `atom_context` fusion and
atom pooling inside `ContactGrammar` (150x attenuation from `context` to
`mean_state`, ~3,400x in the protein-token Jacobian). **Stage P weakens it**:
an over-driven centered objective with 4.6x amplified gradient into `embed`
extracted alignment of only +0.022 against a +0.10 threshold, so the audit
established attenuation but not *information*. Families 1 and 4 therefore both
begin by testing whether pre-fusion `atom_context` carries truth-aligned
protein-ligand information at all. If that shared gate fails, learned pooling
and basis reallocation are rejected rather than widened.

Standing constraints added 2026-08-16:

- protein counterfactuals must be **centered** to speak to ordering; every
  uncentered wrong-protein number in R0-R14 measures target level only;
- no pocket, contact or "biologically localized" language: the protein path is
  exactly invariant to residue-slot permutation;
- `include_meta_test` is fail-closed and opening it requires a written
  authorization recorded in the artifact.

## Independent M0 lane

`report/meta_fewshot/stageM0_msa_probe_20260816/PREREGISTRATION.md` is an independent
protein-calibration diagnostic. It may run only after recording the local MMseqs2
executable and a governed UniRef database snapshot. It may not consume meta_test,
share model-selection decisions with A2, or be represented as DTA improvement.

## Stage D: panel-context level + orthogonal level/shape (opened 2026-08-17)

`tools/research/stageD_level_panel/`. Preregistration frozen before any arm
trained. The D0 diagnostics re-audited the Stage C boundary and answered the
five governing re-examination questions (authority: `D0_REPORT.md` plus
`D0_AUDIT_DECOMPOSITION.json`, `D0_LEVEL_IDENTIFIABILITY.json`,
`D0_LEVEL_ANATOMY.json`, `D0_OCCUPANCY_STRATA.json`):

1. the level/shape decomposition is per episode (per draw), not per canonical
   target; the drawn-panel part of "level" is small (0.013-0.034 pK^2);
2. the calibrated constant reads meta_val labels (disclosed REFERENCE); the
   meta_train-only constant is 2.15-2.17, which the tested features do beat;
3-5. level is a joint property of protein, assay and the tested ligand panel:
   within meta_train, component identity explains 46% in-fold but transfers
   -1.1%; document identity 70% in-fold but 6.8% out; **panel composition is
   the best transferring covariate (23.9%)** vs protein sequence 11.9%;
   the occupancy ordering signal survives scaffold novelty (r +0.154) and low
   ligand recall (r +0.221); best legal level predictor on meta_val remains
   ESM-650M linear (1.6875), 13.6x the 0.1239 budget.

Stage E candidate: two innovations only — I1 panel-set level readout
(framework), I2 orthogonal level/shape routing (training module). Four arms,
one code path, leak-free internal checkpoint selection (Stage B partition).
Gates G1-G6 and stop rules S1-S4 are frozen in `PREREGISTRATION.md`.
GPU verification (torch.cuda.is_available, model/batch devices, nvidia-smi
utilization) runs before every arm and refuses to train without it.

## Stage F: pairwise learned interaction transport (opened 2026-08-17)

`tools/research/stageF_pairwise/`. Preregistration frozen. Candidate: the
transport kernel becomes a learned pairwise (query, support) edge operator over
embed-space with the fixed Tanimoto kernel kept as an additive anchor, trained
with pairwise signed-gap supervision (query labels loss-only). Arms F / F-ABS
against the frozen Stage E T2 baseline; same seed, budget, partition and
leak-free selection. Gates G1-G5 and stop rules S1-S4 in
`PREREGISTRATION.md`. This is the mechanism family that can consume the
Stage L pairwise signed-gap direction (r +0.270, orthogonal to Tanimoto),
which every prior moment-form or fixed-kernel transport could not.

## Stage G / G2: ESM-650M residue-input lane (2026-08-17, night)

`tools/research/stageG_esm650/`. Stage G single-seed screen: the incumbent
recipe retrained on the local ESM-2 650M protein bank (1280-dim pooled +
128-slot residues; provenance recorded). First arm in the record to improve
MSE, level^2, centered MSE, Spearman, Pearson, CI and cliff sign at EVERY k
against frozen T2; k=0 centered MSE resolves (-0.0396 [-0.0772, -0.0018]);
k=0/k=5 MSE gains unresolved. Controls clean. Gate G2 (resolved MSE gain at
k in {2,3,5}) failed, so the lane stopped at the screen; the continuation is
the newly preregistered multi-seed confirmation (Stage G2,
`PREREGISTRATION_G2.md`): 3 fixed seeds for both arms, pooled component
bootstrap, then freeze and a single meta_test open only if every G2 gate
passes. Davis/KIBA remain a separate later stage.

## Stage H / I: pocket-prior and live-LM lanes (2026-08-17, night)

- Stage H0 (structure/pocket, `tools/research/stageH_pocket/`): local MMseqs2
  coverage audit — 209/387 targets have a homologous holo structure (>=30%
  identity, >=50% query coverage; 152 at >=90%); pocket descriptors extracted
  and probed. Level identifiability: pocket MLP 2.4398 vs 2.6179 constant,
  shuffled control 2.4941 — REJECTED at the identifiability gate; no Stage H
  training authorized.
- Stage I (live ESM-2 150M LoRA lane, `tools/research/stageI_lm/`,
  preregistered): REJECTED — G2 fails (no resolved MSE gain at any k; k=2
  -0.0459, k=3 -0.0383, k=5 -0.0182, all intervals cross zero). Two resolved
  ranking gains recorded as observations (k=2 Spearman vs the frozen live
  control, k=3 Pearson vs T2). Engineering note: chunked LoRA backward on long
  proteins silently OOM-killed the process; the first-chunk gradient bound is
  documented for any future LM lane.

## Stage J: assay-aware level head (2026-08-18, closed)

`tools/research/stageJ_assay/`. D0c measured journal/publisher provenance as
the strongest single legal level covariate (1.619 vs 2.155 constant; shuffle
2.522; 100% meta_val share). Three preregistered arms (J / J-NOPAIR /
J-NOJRNL) vs frozen T2. REJECTED: k=0 level^2 1.73 -> 1.30 (best on record,
unresolved) but k=2/3 ranking degrades with resolved intervals; the learned
zero-shot level head cannot beat support-label calibration at k>=1 without
degrading ordering. See `REPORT.md`.

## Required stage artifacts

Each stage keeps only `PREREGISTRATION.md`, `RESULT.json`, `REPORT.md`, necessary
prediction rows and a loadable admitted checkpoint. Delete duplicate smokes,
progress logs and failed checkpoints after consolidating their verdict. Update
`history.md`, this file, `report/CURRENT_MODEL_EVIDENCE.md` and
`report/EVIDENCE_LEDGER.md` after every decision.

## Directory lifecycle

- `tools/research/` is the project's `/research` workspace. It contains only
  unadmitted hypotheses, probes and stage-specific experimental code.
- Research that passes its preregistered gates must be moved, not copied:
  reusable model code goes to `model/`, executable workflows go to `scripts/`,
  and maintained contracts go to `tools/tests/`. Delete the research copy after
  promotion so that every implementation has one owner.
- `tools/runtime/` contains ignored local executables, downloads and inspection
  helpers. It is never an evidence or source-code authority.
- Root `main.py` is the sole high-level command dispatcher. It may orchestrate
  admitted capabilities from `model/` through `scripts/`, but it must never
  import or expose `tools/research/` experiments.
- Do not recreate root `test/`, `tests/`, `research/` or `LLM/` directories.
  Tests, research prototypes and local tooling are consolidated under `tools/`.
