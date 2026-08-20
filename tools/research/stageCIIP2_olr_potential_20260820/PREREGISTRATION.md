# CIIP-2 Stage: OLR-Potential (Orthogonalized Ligand-Routed Interaction) — Preregistration

Stage id: stageCIIP2_olr_potential_20260820
Issued: 2026-08-20, before any training in this stage
Authority: CIIP-2 research report (report/research_ideas/ciip/CIIP2_RESEARCH_REPORT_20260820.md)
Production code: UNTOUCHED (no model/ or scripts/ changes; successor stage only)
This prereg is frozen at issue; amendments require a new dated addendum file and
never edit frozen rules.

## 1. Object of study

One deployable scalar field s_theta(P, L) over (construct sequence, ligand),
parameterized as a ligand-conditioned residue router over frozen ESM-2
residue states, with structural orthogonalization of main effects and a
cross-fitted residual objective:

    s(P,L)      = alpha(P,L)^T beta(L)                    [C1 router]
    OID         = panel/protein-axis centering layers     [C2]
    target      = cross-fitted residual of c_p(L)         [C3]
    weights     = assay-gain w_L from train WT rows       [5.3]

Every scientific contrast is a finite difference of s_theta. Mutation
coordinates appear ONLY in control arms and the optional teacher, never in
the deployed path.

## 2. Frozen data contract

- Inputs (SHA-pinned): stageCIIP_potential_bridge/DATA1A.json,
  DATA1A.npz, DATA2X2.json, DATA2X2.npz, and the X0c ESM residue cache
  q1_esm_cache.npz (640-dim residue states). Reference SHAs are inherited
  from DATA1A['inputs_sha256'] and DATA2X2['inputs_sha256'].
- Endpoint: Duong-Ly raw % inhibition; target c_p(L) as frozen in CIIP-1A.
- Coverage: 65 admitted pairs, 49 ESM-covered. All primary analyses run on
  the covered subset (leakage parity with CIIP-1A controls).
- Ligand features: ECFP4 2048 (frozen X0 asset), reused unchanged.
- New derived data (this stage, written once, SHA-pinned):
  (a) mutation-erased residue states (X-substitution) for all covered
      constructs, via the same ESM pipeline as the cache (exactness rule:
      WT_erased == MT_erased elementwise per pair);
  (b) per-ligand assay-gain weights from WT rows of TRAIN parents only;
  (c) parent-disjoint split tables (CIIP-1B).

## 3. Splits (frozen construction rules)

- S1..S3 (CIIP-1A scope): pair-level 60/20/20 stratified per parent, seeds
  derived by stable RNG from ["stageCIIP2", "split", k] for k in 1..3.
  S1 replicates the CIIP-1A split (parity anchor) by reusing DATA1A's
  pair_split.
- SPB (CIIP-1B scope, primary deployment claim): parent-disjoint over the
  12 multi-pair covered parents (ABL1, ALK, CKIT/KIT, CMET/MET, EGFR,
  FGFR3, FGFR4, FLT3, LRRK2, PDGFRA, RET, TIE2/TEK): rank parents by
  covered-pair count (desc, name asc), greedily move parents to TEST until
  test holds >= 20% of covered pairs and >= 5 parents; remaining parents
  split 60/40 train/val by pair count parity. Parents with a single pair
  stay in train. Rule frozen before inspection of results.
- Ligand side: shared-ligand primary (all 183); scaffold-cold secondary
  (Bemis-Murcko scaffold clustering, 5 clusters by frozen seed; report
  hold-one-scaffold-out table; never used for selection).

## 4. Arms (single-variable ladder; all trained identically unless noted)

| arm | description | tests |
|---|---|---|
| A0-prior | train-mean centered ligand profile m(L) (no training; analytic prior) | corrected ligand-only baseline (fixes CIIP-1A design gap) |
| A1-bilinear | s = alpha(P)^T psi(L) on mean-pooled full-sequence ESM (no router) | pooling vs routing |
| A2-router | LCRR (C1) | routing value |
| A3-oid | A2 + OID centering layers (C2) | orthogonality value |
| A4-cfoie | A3 + cross-fitted residual target (C3) | nuisance orthogonality value |
| A5-gain | A4 + assay-gain weights | assay-link modeling value |
| A6-teacher | A5 + PMSTD distillation (C4; optional, separately gated) | privileged recoverability |
| C-perm | A5 trained with within-pair permuted ligand labels | train-time permutation control |
| C-randprot | A5 with random pair partner (parent-preserving) | protein-conditioning control |
| C-erased | A5 evaluated on mutation-erased sequences | mutation-information control |
| C-wrongmut | A5 predictions assigned to same-parent sibling pairs | mutation-specificity control |
| C-famprior | family-mean profile prior on SPB | trivial transfer control |
| C-free | free pairwise head (non-deployable) | integrability ceiling |

A1..A6 are the falsification ladder; C-* are controls. The DEPLOYED
candidate is A5 (A6's student if A6 is admitted by its own gates).

## 5. Model specification (frozen)

- Residue states: frozen ESM-2 640-d (cache); optional +1 surprise channel
  u_i = log p(a_i | masked) via ESM-2 masked scoring — OFF by default, ON
  only as a separate labeled arm A5s, never mixed into primary claims.
- Router: heads H in {1,4}; per head: q = MLP_phi(z_L) (2048->64), keys
  W_k h_i (640->64), temperature 1/sqrt(64), softmax over residues.
- alpha(P,L) = concat_h sum_i a_i W_v,h h_i in R^{H*r}, r in {4,8,16}
  (default r=8, H=1).
- beta(L) = MLP_psi(z_L) (2048->r).
- b_P, b_L, mu: linear heads on mean-pooled states / ligand MLP; used ONLY
  by the raw-endpoint auxiliary loss.
- Training: AdamW lr 1e-3, wd 1e-4, batch = whole pair panels (grouped so
  that OID protein-axis uses only same-split constructs), epochs <= 300,
  early stop on val weighted centered MSE (patience 40), seed set
  {11,22,33,44,55} for Phase 5; Phase 4 smoke uses seed 11 only.
- Capacity cap: total trainable params <= 2.0 M (asserted in tests).
- No test-time adaptation of any kind; no closed-form solvers; inputs to
  the deployed path are (sequence, ligand) only (asserted in tests).

## 6. Objectives

- Primary: weighted centered MSE on residual targets
  L = sum_p sum_L w_L (r_p(L) - [D_P s_tilde](L))^2, r = c - m_hat^{(-p)}.
- Nuisance m_hat: MLP (2048->64->1), K=3 folds cross-fitted BY PARENT on
  train pairs only; the fold model that excluded pair p's parent predicts
  r_p even within train (cross-fitted residuals for training targets).
  Val/test residuals always use folds fit without their parents.
- Auxiliary (never for interaction claims): raw-endpoint MSE training
  b_P + b_L + mu with s detached.
- Assay-gain weights: w_L proportional to ybar_W(L) (100 - ybar_W(L)) over
  train-parent WT rows, clipped to [0.25, 4] after normalization to mean 1.
- Teacher (A6): distill centered contrasts only, lambda_d in {0.5} frozen;
  teacher = A2 architecture + oracle site-index channel.

## 7. Metrics and reporting (frozen; all mandatory per arm)

centered MSE; raw MSE (auxiliary, separate table); R2; explained variance;
Spearman; Pearson (undefined stays undefined; no NaN->0); dead-zone sign
accuracy (|c|>10); slope (truth-on-pred regression, diagnostic);
variance recovery; scale recovery; nonconstant coverage; rank-evaluable
denominator; per-parent table; mutation-class table (frozen lists below);
assay saturation strata (WT-mean bands <50 / 50-90 / >90); scaffold strata
(secondary); parent bootstrap 2000 draws (seed 20260821); paired
correct-vs-control bootstrap; 5-seed x 3-split CIIP-1A table; SPB
parent-disjoint primary table.

Frozen mutation classes (Duong-Ly labels):
- gatekeeper: ABL1 T315I, KIT T670I, EGFR T790M
- P-loop / JMD / other hotspots (reported as "annotated-hotspot"): ABL1
  Y253H, ABL1 E255K, ABL1 Q252H, KIT V560D, EGFR L858R, EGFR d746-750*,
  RET M918T, RET C634R/W, FLT3 D835*, FGFR3 R248C, MET D1228N/H, MET
  Y1230C/H/D, ALK F1174L, LRRK2 G2019S
- all remaining: "other"

## 8. Gates and verdict rules (frozen)

Phase 4 smoke (seed 11, S1): PASS iff
 (a) pipeline end-to-end, all metrics finite or explicitly undefined;
 (b) C-perm R2 <= 0.02 or <= 10% of A5 R2 (permutation destroys signal);
 (c) A0-prior R2 in [0.08, 0.18] on covered test pairs (baseline sanity,
     audit-predicted ~0.13);
 (d) nonconstant coverage in [6/9, 9/9] for A5;
 (e) structural tests green (Section 10).
Instrument qualification (before any interpretation of real-data results):
 planted synthetic interaction field (rank 4) on real sequences/ligands
 with assay-gain transform and noise matched to observed variance
 components (parent-shared 134.8, mutation-specific 89.7); A5 must recover
 planted residual Delta R2 >= 0.25 over A0-prior at alpha=0.1 via parent
 bootstrap; else the instrument is declared underpowered for the claim and
 the programme stops honestly at "UNRESOLVED (power)".

Phase 5 CIIP-1A scope: report only; no mechanism claim.
Phase 5 SPB primary verdict (the deployment claim):
 SUPPORTED iff (i) A5 residual R2 > A0-prior R2 with paired parent
   bootstrap 90% CI excluding 0 (pooled over seeds), AND (ii) C-perm
   destroyed per Phase-4 rule, AND (iii) C-erased degrades A5 by >= 30%
   residual R2, AND (iv) C-wrongmut correct-assignment advantage > 0 in
   >= 4/6 test parents (or all available), AND (v) leave-one-parent-out
   sign stability of the pooled Delta R2.
 NOT_SUPPORTED iff (i) fails with CI excluding 0 in the negative direction
   AND (ii) passes.
 UNRESOLVED otherwise (includes CI straddling 0).
Secondary scaffold-cold: reported, never selection-bearing.

## 9. Stop rules (from CIIP-2 report Section 7; any one terminates)

All ten stop rules are inherited verbatim from the CIIP-2 research report
Section 7 and are restated in the stage runner's output header for
auditability.

## 10. Structural tests (run before training; no labels touched)

T1 antisymmetry: s(P,L) finite-difference contrast flips sign under
   construct swap (max |D + D_swapped| < 1e-5).
T2 cycle-zero: random A->B->C->A ligand cycles sum to < 1e-5.
T3 centering: panel-centering layer zeroes prediction row means exactly
   (< 1e-6) and protein-axis centering preserves them.
T4 coordinate-free contract: deployed forward signature accepts ONLY
   (residue_states, ligand_features); no site index parameter exists in
   the A1-A5 modules (asserted by signature inspection).
T5 capacity: trainable params <= 2.0 M.
T6 permutation wiring: C-perm label map is a within-pair derangement;
   fold assignment respects parent groups (no parent spans folds).
T7 leakage guards: gain weights computed from train WT rows only (asserted
   by input mask); val/test residuals use parent-excluded folds.
T8 determinism: fixed seed reproduces identical losses over 3 epochs.
T9 erased-input equality: WT_erased == MT_erased states per pair (<1e-5).
T10 no closed forms: no numpy.linalg pinv/solve/cholesky/ridge calls in
   the deployed modules (grep assertion).

## 11. Compute budget

Model <= 2M params; 9k cells; ~seconds/epoch on one GPU; Phase 4 < 10 min;
Phase 5 (14 runs x 5 seeds x 3 splits + SPB) < 3 h GPU; memory < 2 GB.
ESM erased-sequence forwards: ~67 sequences, minutes, cached once.

## 12. Authorized successors

Upon Phase-4 PASS: Phase-5 multi-seed evaluation (this prereg, no new
document needed). Phase-6 (BindingDB Ki bridge, Tanimoto increment, SCNP
few-shot) requires a NEW preregistration and remains blocked until the SPB
verdict is SUPPORTED.
