# Stage CIIP-1A 2x2 root-cause diagnostic preregistration — ORACLE-COVERED SUBSET (2026-08-19)

Frozen BEFORE any 2x2 computation. Authorized scope (per independent
review, 2026-08-19): root-cause diagnosis of the Stage-1 potential
collapse ONLY. NOT authorized: CIIP-1B, BindingDB Potential Bridge,
production integration, or any claim that this stage is CIIP-1A PASS.

## Naming

This stage is the **oracle-covered subset diagnostic**. Coverage bias
verdict (frozen in DATA2X2.json): the 16 uncovered pairs are 4 whole
families (ALK, MET, LRRK2, TEK-Y1108F; all mutation positions > 1020,
the Q1-frozen ESM window bound), with higher target variance and
informative fraction than the covered set (var_true median 243 vs 169;
informative frac 0.40 vs 0.35). The covered subset keeps 49/65 pairs
(32 train / 8 val / 9 test), 18/22 parents, and 9 test pairs across 6
parents. All cells use EXACTLY this matched 49-pair subset with the
original train/val/test assignment — no cell may use 65 pairs.

## Representation factor

- KLIFS: frozen pocket one-hot rows (1700 = 85x20), identical to Stage 1.
- oracle_local_esm: ESM-2-150M hidden-state mean over the radius-6
  window at the VERIFIED mutation coordinate (Q1-qualified
  construction), 640-dim, per pair for BOTH WT and variant rows (WT
  window taken at the same pair site). ORACLE-localized: it consumes
  the mutation annotation; it is a positive-control representation and
  can NEVER be the deployable protein representation (a
  mutation-annotation-free pocket/residue router would need separate
  qualification).

## Objective factor

- joint: L = L_contrast + 1.0 * L_abs (unchanged from Stage 1).
- centered-only: L = L_contrast ONLY. b_P/b_L receive no gradient
  (recorded in the gradient-coverage report). This is a CIIP
  qualification objective; a centered-only result is NOT a BindingDB
  absolute-DTA result.

L_abs details (joint cells): KLIFS cell — train-pair rows excluding
val/test-pair rows (Stage-1 rule); oracle_local_esm cell — VARIANT-row
cells only, each using its own pair's site window (WT rows have no
unique site; frozen rule).

## Frozen protocol (all four cells identical except the two factors)

- Model: f = b_P + b_L + s, s = alpha^T psi (rank 8, hidden 64, the
  SAME potential formula); d_p = 1700 (KLIFS) or 640 (oracle ESM);
  d_l = 2048 ECFP4; xavier init, torch.manual_seed(1).
- Training: AdamW lr 1e-3, wd 1e-4, 200 epochs, batch 512, grad clip
  10.0, end-to-end gradients only. Seed 1 (single-seed root-cause
  run). Keyed rng streams: "stageCIIP2x2|order|{cell}|{seed}|epoch|{ep}",
  "stageCIIP2x2|abs|{cell}|{seed}|epoch|{ep}".
- Checkpoint: best val contrast MSE over the 8 covered val pairs.
- Metrics (covered test pairs, per pair): R2 (PRIMARY; = 1 - MSE/var_true,
  defined also for constant predictions, = 0), Spearman/Pearson (null
  when the prediction is constant — never treated as zero), dead-zone
  sign accuracy (dead zone = 10 % inhibition units), centered MSE,
  OLS slope, var_pred, nonconstant flag, parent, n_lig.
- Variance decomposition per cell (at checkpoint): var(s), var(b_P),
  var(b_L), var(f) on the cell's abs sample; gradient coverage of
  every parameter at step 1.

## Diagnostic gates (interpretability; NOT a CIIP-1A PASS)

- A cell is COLLAPSED iff any of: N_prediction_nonconstant < 5/9,
  N_rank_evaluable < 5/9, parent coverage < 4/6, pair-mean R2 <= 0.02.
- Effects on pair-mean R2 (primary) and Spearman (secondary),
  parent-cluster bootstrap 2000 draws, keyed rng:
  rep_main_joint = R2(ESM,joint) - R2(KLIFS,joint)
  rep_main_centered = R2(ESM,centered) - R2(KLIFS,centered)
  obj_main_klifs = R2(KLIFS,centered) - R2(KLIFS,joint)
  obj_main_esm = R2(ESM,centered) - R2(ESM,joint)
  interaction = obj_main_esm - obj_main_klifs
- An effect is ESTABLISHED iff bootstrap 2.5% lower bound > 0 AND
  |point estimate| >= 0.05; ABSENT iff |point| < 0.02; otherwise
  AMBIGUOUS. No single-cell success may be attributed to either factor
  alone; conclusions come from the effect contrasts only.
- Few-pair safeguard: a parent contributes a median only with >= 1
  evaluable pair; an effect is flagged as few-pair-driven if removing
  any single parent flips the sign of the point estimate.

## Promotion rule (frozen)

After this diagnostic: request authorization for the control arms
(family-preserving shuffle, random-window control, ligand-invariant
shift, ligand-only, free-pairwise diagnostic) on the SAME matched
inputs and budget, and only then the 3-seed CIIP-1A successor — each
step separately authorized. CIIP-1B requires new authorization
regardless of any outcome here. Biological protein-conditioned signal
stays UNRESOLVED until a deployable-representation successor is
qualified.
