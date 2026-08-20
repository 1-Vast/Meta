# CIIP-S1 Signal Decomposition & Estimand Ladder — PREREGISTRATION (frozen)

Stage directory: tools/research/stageCIIP_signal_decomp_20260820/
Frozen (SHA-256 pinned below) on: 2026-08-20, BEFORE any S0 computation beyond
input-SHA verification, and BEFORE any S1 fitting.
Scientific authority: report/research_ideas/CIIP_SUCCESSOR_STAGE_RESEARCH_PLAN_20260820.md
Section 9 (verbatim extraction below). Conflicts: plan document wins over this
file's completions; implementation freedom not covered by the plan is frozen in
this file or in dated addenda only, always before seeing the corresponding
results.

## A. Verbatim extraction of plan Section 9 (9.0-9.4)

## 9. The chosen next step: CIIP-S1 stage specification (preregistration
## draft — freeze and SHA-256 before ANY computation)

### 9.0 Authority and prohibitions

- Writes only to a new directory tools/research/stageCIIP_signal_decomp_<date>/.
- Never modifies: CIIP-1A artifacts, context-propagation artifacts, model/,
  scripts/, dataset/, meta_test (never read at all).
- Duong-Ly test-pair labels: final evaluation of frozen arms only.
- No ridge/closed-form/pseudoinverse; no test-time gradients; end-to-end
  gradient training only; keyed SHA-256 rng streams; no Python hash().
- B-branches labelled: "mutation-free but not necessarily biologically causal".

### 9.1 Stage S0 — read-only audits (before any fitting)

1. Data coverage: pair/parent/position tables; ESM-length exclusion map (P9).
2. Parent overlap: split x parent incidence; WT-row sharing graph.
3. Mutation coordinate: re-verify all 65 mappings against Q0B audit; alias
   ledger check.
4. Ligand overlap: per-pair common-ligand matrix; per-ligand row coverage.
5. Assay semantics: endpoint namespace check (% inhibition only);
   out-of-range cell census (E6 re-derived); concentration metadata from
   source supplement (read-only).
6. Censoring: confirm absence of censoring annotations (blocks interval-
   censored formulations — recorded as a data limitation).
7. Re-derive Section-4 diagnostics under the preregistration (train+val
   only) with keyed rng; discrepancies vs E1-E6 reported.
8. Power computation (Q9.3) and frozen MDE table per contrast.
9. ESM erasure-inference pass (F2): rebuild erased sequences, assert
   WT-erased == variant-erased strings, inference on CPU, cache under the new
   stage directory with SHA-256; assert max erased-input delta <= 1e-5.
10. Leakage audit output: any channel found must be neutralized or recorded
    as a claim restriction BEFORE S1 starts.

S0 exit gate: audits 1-8 complete, no unresolved leakage channel, erasure
cache valid. Only then S1.

### 9.2 Stage S1 — estimand ladder and representation arms

- Targets T0, T0m, T1, T2, T3 (Section 5); T2 uses cross-fitted F9 (train
  parents only for val/test residuals).
- Arms F1-F9 (+F10 only if licensing verified; expected absent).
- Form-1 potential: frozen architecture (rank 8, hidden 64, AdamW 1e-3,
  wd 1e-4, 200 epochs, batch 512, grad clip 10) for F1/F6/F5/F7/F8 on T1/T3.
- Form-2 probes: same optimizer/budget; small MLP (hidden 64) on
  (protein-features, ECFP4) for F2/F3/F4/F9 on T0/T0m/T1/T2/T3.
- Controls (mandatory): ligand-label permutation (keyed, within pair);
  same-parent wrong-mutation (variant window from another covered pair of
  the same parent — reuse of frozen construction); family-preserving shuffle;
  random window; protein-invariant shift. All controls run through the
  identical metric pipeline.
- CPU smoke (<= 30 min, 2 pairs, 5 epochs) must pass structure tests; then
  single seed; multi-seed (3) ONLY if single-seed structure + all negative
  controls behave as frozen (floors at floor, erasure null at zero).
- Metrics per pair per arm: nonconstant flag/rate, var_true, var_pred,
  centered R2, OLS slope, scale_ratio = sqrt(var_pred/var_true), MSE,
  dead-zone (10 units) sign accuracy, Spearman (undefined if constant —
  never 0), per-parent aggregates; paired effects with parent-cluster
  bootstrap (2000 draws, keyed) + leave-one-parent-out sign stability.
- Frozen adjudication rules per proposition (lo2.5 > 0 on the primary
  paired contrast + LOPO sign stability; point thresholds frozen at freeze
  time after the S0 power table, never after results).

### 9.3 Adjudication outcomes and consequences

| outcome | definition (primary contrasts) | authorizes | closes |
|---|---|---|---|
| S1-PASS-B | B-arm beats F7/F8 floors on T1 | Direction-2 preregistration drafting | nothing |
| S1-PASS-A | F1 beats F2-erased on T0/T0m | proposition-A reporting line (mutation overall effect; still not ligand-conditioned) | nothing |
| S1-PASS-Csharp | F1 beats F2-erased on T2 | Direction-3 scoping review | nothing |
| S1-NULL-ALL | all cells at floor | — | Duong-Ly interaction route (both levels); final boundary document update |
| S1-UNRESOLVED-Csharp | C_sharp CI crosses 0 while B passes | Direction-2 only | Direction-3 recorded as underpowered-closed unless new data arrives |

In EVERY outcome: production model/scripts unchanged; CIIP-1B, BindingDB
bridge, and any deployable-representation claim remain NOT AUTHORIZED until
the relevant Direction-2/3 stage passes its own preregistered gates.

### 9.4 What CIIP-S1 licenses on PASS and forbids on FAIL

- PASS-B licenses ONLY: drafting a Direction-2 preregistration (mutation-free
  family conditioner, parent-disjoint). It does not license touching
  BindingDB, production code, or calling the parent-level component
  "mutation-specific".
- PASS-A licenses reporting a mutation overall-effect result; forbids
  ligand-conditioned framing.
- PASS-Csharp licenses Direction-3 scoping; forbids any deployment claim
  (oracle-bound).
- NULL-ALL forbids further Duong-Ly interaction work and routes the
  programme back to the performance track with the boundary updated.

---


## B. Frozen stage specifics (implementation freedoms, frozen now)

### B.1 Environment
- Python: /d/anaconda/envs/drug/python.exe (conda env drug; torch 2.6.0+cu124).
- ESM erasure inference: CPU only (plan 9.1 item 9); model
  facebook/esm2_t30_150M_UR50D via stageX x0_i2.load_esm; truncation
  max_length=1022 (identical to frozen X0c cache rule).
- S1 training: CPU (small probes); seeds frozen at {1, 2, 3}; single-seed
  gate before multi-seed.

### B.2 Frozen inputs (SHA-256 re-verified 2026-08-20 before freezing)
- tools/research/stageCIIP_potential_bridge/DATA1A.json
  sha256=1c2b92dfc9f3f84676a25a15f331cefc178005afa2207125d6dc09b771defc68
- tools/research/stageCIIP_potential_bridge/DATA1A.npz
  sha256=40f69509ff88005d05cd2e65873d08caa8e10c1b44da9e9aff256c37866f52f1
- tools/research/stageCIIP_potential_bridge/DATA2X2.json
  sha256=8caf2d854788c337729aff1eb8f7fd951ca8e2237f1b0c61e470dcde6a5c6ff0
- tools/research/stageCIIP_potential_bridge/DATA2X2.npz
  sha256=cb967441d7f43f2b624f313cf29455edc6dc2636f593ec4b4f3ec8e8c0cf76ff
- tools/research/stageX_csc_signal/stageX0c_measurement_qualification_20260818/q1_esm_cache.npz
  sha256=c8b59e33c3011ca898c51d02957f44f1865c89799ce9eb67c5327d8ff2d545db
- tools/research/stageX_csc_signal/stageX0c_measurement_qualification_20260818/Q0B_MAPPING_AUDIT.json
  sha256=f5907e286466cd83eb8e5c31fc7895c31d680e24c473e6140c0b347ef745928d
  (audit reference; read-only cross-check)
- meta_test (BindingDB sealed set): NEVER READ in this stage.
- Duong-Ly test-pair labels: final evaluation of frozen arms only; never used
  for fitting, normalization, feature construction, model selection, retrieval,
  or checkpoint selection (checkpoint selection uses val pairs only).

### B.3 Random streams (no Python hash(); keyed SHA-256 streams via
x0_common.stable_rng, imported read-only from the frozen X0 stage)
- Namespace prefix: "S1".
- Streams: S1.winperm.<arm>.<seed> (random-window positions),
  S1.famshuf.<arm>.<seed> (within-parent row shuffle),
  S1.ligperm.<arm>.<seed> (within-pair ligand-label permutation),
  S1.boot.<contrast>.<seed> (parent-cluster bootstrap draws),
  S1.power.<id>.<seed> (S0 power simulation),
  S1.diag.<id> (S0 diagnostic baselines),
  S1.order.<arm>.<epoch>.<seed> (batch order),
  S1.init.<arm>.<seed> (parameter init).
- Bootstrap: 2000 draws, parent-cluster resampling, 2.5/97.5 percentiles;
  bootstrap mean is never a point estimate.

### B.4 Targets (all from the same legal labels; endpoint = percent inhibition)
- T0  d_vl = y_var - y_wt (per pair, per ligand).
- T0m mean_l(d) (per-pair scalar; mutation severity).
- T1  c_vl = d_vl - mean_l(d).
- T2  c_vl - profile_parent(i) where profile is the cross-fitted F9
  parent-profile: train pair -> mean of c over same-parent TRAIN pairs
  excluding itself (leave-pair-out); if the pair is the only train pair of its
  parent, profile = 0 vector and the pair is flagged siblingless (pair excluded
  from T2 paired contrasts, counted in coverage); val/test pair -> mean of c
  over same-parent TRAIN pairs only (val/test labels never enter any
  nuisance). T2 for an eval pair with no same-parent train pair is undefined
  and flagged (S0 verifies none exist on test).
- T3  within-pair ranking of c (pairwise logistic rank loss over the pair's
  ligands; evaluated with the B.7 metric contract).

### B.5 Arms and forms (frozen matrix)
Form-1 potential g = s(Pv,L) - s(Pw,L), rank 8, hidden 64, AdamW 1e-3,
wd 1e-4, 200 epochs, batch 512, grad clip 10 (frozen budget; checkpoint =
best val full-target MSE on val pairs):
- F1 oracle mutation-site radius-6 window-mean ESM (640-d) on T1, T3.
- F5 family-preserving shuffle (within-parent keyed permutation of pair
  protein features) on T1, T3.
- F6 random local window (keyed non-site position, same radius) on T1, T3;
  interpreted as context-propagation value measurement, not a null.
- F7 ligand-only (protein features zeroed) on T1, T3 — structural zero
  contrast; asserted in tests (CIIP-2 audit fact (a)).
- F8 protein-invariant shift (ligand features zeroed; per-pair constant) on
  T1, T3 — structural zero contrast; asserted in tests.
Form-2 probe MLP (Linear(in,64) -> GELU -> Linear(64,1)), same optimizer and
budget, checkpoint = best val MSE on the SAME estimand being trained:
- protein features per arm:
  F1f  concat(wt radius-6 window-mean ESM, var-minus-wt window) (1280-d);
       targets T0, T0m, T1, T2, T3. (Oracle; A and C_sharp tests.)
  F2   erased full-sequence pooled ESM state (640-d; identical WT/MT by
       construction — asserted); secondary F2w = erased-site radius-6 window
       mean (640-d). Targets T0, T0m, T1, T2, T3. Counterfactual arm: B
       branch; the X token position encodes the mutation COORDINATE
       (acknowledged, non-deployable leak; B label applies).
  F3   concat(full-seq pooled WT ESM, full-seq pooled variant ESM) (1280-d);
       deployable. Targets T0, T0m, T1, T2, T3.
  F4   concat(KLIFS pocket one-hot WT row, variant row) (3400-d); deployable,
       structurally limited. Targets T0, T0m, T1, T2, T3.
  F7f  ligand-only (protein zeroed). Floor/reference = shared ligand pattern.
       Targets T0, T0m, T1, T2, T3.
  F8f  protein-only (ligand zeroed; per-pair constant output). Zero floor.
       Targets T0, T0m.
  F9   parent-profile predictor (cross-fitted as B.4; NOT trained: the
       prediction IS the profile). Ceiling reference on T1; residualizer for
       T2. Reported on T1 and T3.
- F10 external features: NOT AVAILABLE (UniRef snapshot absent locally);
  recorded absent; no arm added.
- Ligand featurization: ECFP4 2048-d (z1 lig); no test-label-derived
  statistics anywhere.

### B.6 Mandatory negative controls (same metric pipeline as all arms)
- Ligand-label permutation within pair (keyed): permute the target vector
  within each eval pair, recompute metrics on frozen arms (evaluative
  control), plus a training control C-perm for Form-1 F1.
- Same-parent wrong-mutation: at evaluation, replace the eval pair's variant
  protein features with those of a different covered same-parent train pair
  (keyed choice); defined only when a same-parent covered train pair exists
  (S0 verifies for all 9 test pairs).
- Family-preserving shuffle (F5), random window (F6), protein-invariant
  shift (F8/F8f), erasure null (Form-1 on erased inputs is structurally
  zero — asserted in tests, not trained).

### B.7 Metric contract (per pair per arm per estimand)
nonconstant flag; N_nonconstant/N_total; var_true; var_pred;
scale_ratio = sqrt(var_pred/var_true); centered MSE; centered R2 (centering
within pair over ligands; for T0m the cross-pair correlation replaces
per-pair R2); OLS slope (pred on true within pair); dead-zone sign accuracy
(|true| > 10 units); Spearman within pair (undefined if pred or true
constant — recorded NA, NEVER 0); N_rank_evaluable/N_total; per-parent
aggregates. Paired contrasts: per-pair metric differences, parent-cluster
bootstrap (2000 draws, keyed), lo2.5/hi97.5, LOPO sign stability (sign of
the aggregate when excluding each parent in turn).

### B.8 Adjudication thresholds — PLACEHOLDER until S0 power table
Point thresholds for S1-PASS-B / S1-PASS-A / S1-PASS-Csharp / S1-NULL-ALL /
S1-UNRESOLVED-* are frozen in a dated addendum (S1_ADDENDUM_THRESHOLDS.md)
AFTER the S0 power table and BEFORE any S1 fitting. Directional rule
(frozen): lo2.5 > 0 on the primary paired contrast AND LOPO sign stable.
Primary contrasts (frozen):
- B:   F2 vs F7f on T1 (Form-2; per-pair centered R2, paired).
- A:   F1f vs F2 on T0m (cross-pair Spearman over the 9 test pairs, paired,
      parent-cluster bootstrap) — secondary: T0 per-pair centered R2.
- C_total: F1 (Form-1) vs F9 on T1 (per-pair centered R2, paired).
- C_sharp: F1f vs F2 on T2 (per-pair centered R2, paired).
- Deployability: F3 vs F9 on T1 and F4 vs F9 on T1 (paired).
- NULL-ALL: every trained arm cell at or below its F7f/F8 floor cell within
  the frozen tolerance from the thresholds addendum.

### B.9 Execution ladder (no skipping)
1. CPU smoke: 2 pairs, 5 epochs, <= 30 min; ALL structure tests green
   (identity=0, antisymmetry, erasure null strictly zero, F7/F8 at floor).
2. Single seed (seed 1): all arms + all controls on the frozen matrix.
3. Multi-seed {1,2,3}: ONLY if single-seed structure + all negative controls
   behave as frozen (floors at floor, erasure null at zero, permutation
   without systematic gain).
4. Any negative-control failure: stop that branch, record; never tune to
   pass.

### B.10 Deliverables
PREREGISTRATION.md + PREREGISTRATION_SHA256.txt; OVERLAP_ANALYSIS.md (CIIP-2
terminal-state adjudication); S0 addendum with frozen thresholds;
S0_AUDIT.json/md; erased-ESM cache + SHA; RESULT.json; REPORT.md (verdict
table + authorization block); commands.jsonl; structure and data-contract
tests; SHA256SUMS; FAILURES.md (or explicit no-failure statement); append-only
syncs to history.md, task.md, report/EVIDENCE_LEDGER.md.

Thresholds in B.8 are PLACEHOLDERS by design; the addendum replacing them is
frozen before S1 fitting and may not be edited after any S1 result is seen.
