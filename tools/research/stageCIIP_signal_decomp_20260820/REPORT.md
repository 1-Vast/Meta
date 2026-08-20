# CIIP-S1 Signal Decomposition & Estimand Ladder — Stage Report (2026-08-20)

Stage dir: tools/research/stageCIIP_signal_decomp_20260820/
Preregistration: SHA 1fb7133b... (frozen before any computation)
Thresholds addendum: SHA beefb620... (frozen after S0 power table, before any
S1 fitting). Erasure cache SHA ae05ae2d... (49 pairs, string-equality and
embedding-delta asserts green, delta exactly 0.0).
Scientific authority: CIIP_SUCCESSOR_STAGE_RESEARCH_PLAN_20260820.md Section 9.

Endpoint: functional percent inhibition throughout. Never relabeled as
Ki/Kd/pK/DDG. No mutation-coordinate oracle enters any deployment path.

## 1. Governance

- OVERLAP_ANALYSIS.md: CIIP-2 was found COMPLETE and TERMINAL (verdict
  UNRESOLVED (power); gate (b) failed; Phase 5 not executed; goal closed).
  No concurrent training exists; plan 9.5 designates S1 a compatible strict
  prequel. Governance instruction (2026-08-20) selects S1 as the next
  direction. CIIP-2 code re-audited read-only: no defect that would change
  its terminal verdict; frozen tests re-run 11/11 PASS.
- Writes confined to this directory; model/, scripts/, dataset/, other
  stageCIIP_* dirs untouched; meta_test never read.

## 2. S0 read-only audit (S0_AUDIT.json/md; all items PASS)

- Coverage: 49/65 covered re-derived bit-exactly; 16 excluded pairs ALL
  pos > 1020; covered splits 32/8/9 match DATA1A.
- Parent overlap: one WT row per parent; all 9 test pairs F9-definable
  (sibling counts 1-5); test clusters {ABL1:2, KIT:2, EGFR:1, FGFR4:1,
  RET:2, TEK:1}.
- Mutation coordinates: all 65 match Q0B (0 mismatches); alias ledger clean.
- Ligand overlap: panels 179-183, pairwise common 175-183 -> F7f floor
  mandatory.
- Assay semantics: % inhibition only; raw panel out-of-range 23.02% (full
  panel) reproduces the prior census; WT panel ceiling-loaded (median WT
  mean 89.4; 87/183 ligands WT mean > 90; 0 < 10; mid-zone 96).
  Concentration metadata: NOT present as a data column in the local
  supplement copy -> recorded limitation.
- Censoring: no annotations -> interval-censored forms not identifiable
  (recorded limitation).
- Plan Section-4 diagnostics re-derived (train+val 40): ALL match (energy
  10.05%; sib 0.4425; WT-resid 0.391; diff-parent 0.041; LOPO Spearman
  0.559; LOPO per-pair median R2 0.326; ligand-global 0.060; parent-
  residualized -0.280).
- Power: sigma_R2 0.599, sigma_Spear 0.218; MDE(80%): delta-R2 0.566,
  delta-Spearman 0.208 at n=9 pairs / 6 clusters. Consequence: only ~0.5
  R2-unit effects are confirmable; every outcome carries the power label.
- Frozen-input integrity: radius-6 window features recomputed from the ESM
  cache reproduce DATA2X2.npz EXACTLY (max abs diff 0.0).
- Erasure cache: ERASED_ESM_S1.npz (49 pairs; WT-erased string == MT-erased
  string asserted per pair; max |embedding delta| == 0.0 <= 1e-5).

## 3. S1 design (frozen)

Estimands: T0 (d), T0m (mean_l d), T1 (c), T2 (c - cross-fitted F9
parent-profile), T3 (within-pair rank form).
Arms: Form-1 potential (frozen architecture/budget) F1/F5/F6/C-perm on
T1/T3; Form-2 probe MLP F1f/F2/F2w/F3/F4/F7f on T0/T1/T2/T3 and F8f on T0;
F9 profile ceiling; F10 recorded absent. Controls: within-pair ligand-label
permutation (evaluative + C-perm training control), same-parent
wrong-mutation at evaluation, family-preserving shuffle (F5), random window
(F6), protein-invariant shift (F8), erasure null (structural, asserted).
Metric contract: per-pair nonconstant/var/scale_ratio/centered MSE/R2/OLS
slope/dead-zone-10 sign accuracy/Spearman (undefined -> NA, never 0);
paired contrasts = parent-mean of per-pair deltas, parent-cluster bootstrap
2000 draws keyed, LOPO sign stability. NOTE (frozen before results): the
contrast statistic aggregates per-pair R2 deltas by parent MEANS; per-pair
MEDIANS are reported as secondary because the mean-based statistic is
noise-dominated at n=9 (S0 sanity: F9 mean per-pair cr2 0.018 vs median
0.274 while the ligand floor means 0.058 vs median 0.125).
Execution ladder: structure tests 13/13 PASS -> CPU smoke (2 pairs, 5
epochs, 17 s) PASS -> single seed 1 -> multi-seed {1,2,3} only if the
single-seed structure and negative controls behave as frozen.

## 4. Results

(RESULTS_SECTION — filled from SEED1/RESULT.json after the frozen gate)

## 5. Verdict table (frozen rules; S1_ADDENDUM_THRESHOLDS_20260820.md)

(VERDICT_SECTION)

## 6. Authorization status (plan Section 9.3/9.4)

(AUTH_SECTION)

## 7. Limitations

- n=9 test pairs / 6 parents: MDE(80%) delta-R2 = 0.566; all claims
  power-labeled.
- Endpoint is single-dose functional % inhibition; no concentration column;
  no censoring; ceiling-loaded WT panel; no replicates (mutation-specific
  variance noise-inclusive).
- Claims restricted to the 49-pair ESM-covered subset (16 pairs excluded,
  pos > 1020).
- F2 (erased context) contains the mutation COORDINATE via the X position:
  counterfactual-only; every B-branch statement carries "mutation-free but
  not necessarily biologically causal".
- Form-2 arms are information-ceiling diagnostics, never production
  mechanisms.
