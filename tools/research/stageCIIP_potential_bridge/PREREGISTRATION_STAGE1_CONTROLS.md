# Stage CIIP-1A CONTROL-ARM preregistration — oracle annotation-shortcut audit (2026-08-19)

Frozen BEFORE any control training. Purpose: determine whether the
oracle local ESM nonconstant response depends on real mutation-centered
protein information, or on mutation annotation / generic local-context
/ protein-main-effect / ligand-only shortcuts. This stage CANNOT
produce CIIP-1A PASS, cannot authorize CIIP-1B, BindingDB Potential
Bridge, production integration, or a deployable representation.

## Data and protocol (identical across ALL arms)

- Matched 49-pair covered subset from DATA2X2.json; ORIGINAL
  train/val/test assignment (32/8/9); 9 test pairs, 6 parents.
- Centered target c_vl = d_vl - mean_L(d_vl); endpoint stays
  percent inhibition.
- Representation: ESM-2-150M radius-6 window-mean, 640-dim, from the
  frozen q1_esm_cache (oracle_local_esm construction), ECFP4 ligands.
- Model: f = b_P + b_L + s, s = alpha^T psi (rank 8, hidden 64, the
  SAME potential); d_p = 640.
- Objective for ALL arms: centered-only (L = L_contrast ONLY). Frozen
  justification: the 2x2 showed the objective effect on ESM is absent
  (-0.017, CI crosses 0), so centered-only does not handicap the
  correct arm, and it removes the joint-cell objective-sampling
  confound entirely. This is a NEW preregistration; the 2x2 prereg is
  untouched.
- Training: AdamW 1e-3, wd 1e-4, 200 epochs, batch 512, grad clip
  10.0, seed 1, end-to-end gradients only. Checkpoint = best val
  contrast MSE on the 8 covered val pairs. Keyed streams:
  "stageCIIPcontrols|order|{arm}|1|epoch|{ep}",
  "stageCIIPcontrols|rows|{arm}|1",
  "stageCIIPcontrols|winpos|{arm}|1".
- Checkpoint/metrics/split/budget identical across arms; arms differ
  ONLY in the protein window fed to the same potential.

## Arms

1. oracle_local_esm_correct — WT and variant windows at the VERIFIED
   mutation site (the positive-control arm).
2. family_preserving_shuffle — each pair's VARIANT window is replaced
   by the variant window of another COVERED pair OF THE SAME parent
   (keyed permutation within parent; WT windows fixed). Preserves
   family/pocket distribution; breaks mutation-window pairing.
3. random_local_window — for each pair, ONE keyed random position p on
   the construct sequence (|p - true_site| > 6; valid window range),
   used for BOTH the WT and variant rows (matched non-mutation window,
   same 640-dim, same ESM model). Tests the mutation-annotation
   shortcut: if this reproduces the correct arm, the signal is generic
   local context, not the mutation site.
4. ligand_only — protein input ZERO for both rows (centered contrast
   identically zero by antisymmetry); the strict zero floor. Spearman
   stays undefined (never treated as zero correlation); the floor is
   compared via centered MSE / R2 / sign / slope.
5. ligand_invariant_shift — per-pair ligand-invariant scalar (the
   centered contrast is exactly zero by construction); excludes the
   protein-main-effect shortcut. Implemented as the constant-zero
   predictor through the same metric pipeline.
6. random_protein — each pair's VARIANT window replaced by the variant
   window of a keyed random covered pair (across parents); all
   ligands/labels/rows/splits unchanged.
7. free_pairwise — DIAGNOSTIC ONLY: g = h([win_wt,win_var,L]) -
   h([win_var,win_wt,L]) on correct windows; conditional expressivity
   ceiling; never a mechanism; report parent coverage and few-pair
   flags.

## Metrics (per test pair, per arm)

var_true, var_pred, scale_ratio = sqrt(var_pred/var_true), centered
MSE, centered R2, OLS slope, dead-zone sign accuracy (dead zone 10 %
inhibition units), Spearman (null if constant; never 0), nonconstant
flag, parent, n_lig. Aggregates must carry N_prediction_nonconstant/N_
total and N_rank_evaluable/N_total denominators.

## Comparisons (frozen; parent-cluster bootstrap, 2000 draws, keyed)

correct vs family_preserving_shuffle, correct vs random_local_window
(PRIMARY annotation-shortcut contrast), correct vs random_protein,
correct vs ligand_only, correct vs ligand_invariant_shift (both zero
floors), free_pairwise vs correct (report-only). Each effect reports
observed_pair_mean_effect, observed_parent_mean_effect, bootstrap CI
{lo2.5, hi97.5}, leave-one-parent-out sign stability (all effects);
bootstrap mean is never a point estimate. The annotation-shortcut
audit also reports feature-level input norms: per-pair ||correct
window delta|| vs ||random matched window delta|| with a paired
bootstrap.

## Verdict rules (frozen)

- ORACLE_LOCAL_SIGNAL_SUPPORTED iff ALL: (1) correct nonconstant
  coverage exceeds random_local_window coverage; (2) correct R2 and
  sign-accuracy exceed the ligand_only floor with bootstrap lo2.5 > 0
  (sign-acc point >= 0.55); (3) correct - family_preserving_shuffle and
  correct - random_protein have bootstrap lo2.5 > 0 on pair-mean R2;
  (4) correct - random_local_window has bootstrap lo2.5 > 0 on pair-
  mean R2 AND sign-acc point >= 0.05; (5) the correct-vs-random gains
  are leave-one-parent-out sign stable and not driven by a single pair.
- ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED iff the correct-vs-random and
  correct-vs-family gaps are absent (|observed pair| < 0.02 or
  lo2.5 <= 0) AND correct nonconstant coverage does not exceed the
  random window's — the response is a generic local-context /
  annotation shortcut; biological signal stays UNRESOLVED.
- Otherwise ORACLE_LOCAL_SIGNAL_UNRESOLVED.
- In EVERY case: deployable protein representation NOT VALIDATED;
  CIIP-1A NOT AUTHORIZED (a new preregistered successor must pass);
  CIIP-1B / BindingDB Bridge / production integration NOT AUTHORIZED.
