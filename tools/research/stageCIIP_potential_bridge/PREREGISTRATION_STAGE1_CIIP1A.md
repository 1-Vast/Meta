# Stage CIIP-1A preregistration — Duong-Ly within-parent representation capacity (2026-08-19)

Frozen BEFORE any CIIP-1A computation. This stage measures whether a
unified integrable potential s_theta(P,L) (alpha(P)^T psi(L), rank 8,
hid 64) can explain the centered mutation effect on the Duong-Ly panel
(within-parent capacity ONLY; cross-parent transfer is CIIP-1B, a
separate stage). No training starts before the Q2d terminal archival
(Stage 0) is complete. Endpoint stays "percent inhibition" everywhere;
never pK/Ki/Kd.

## Data contract (frozen)

- Rows/ligands/features: x0_i1.load_features() — 97 Duong-Ly assay rows,
  183 compounds; protein = KLIFS pocket one-hot (97, 1700 = 85x20);
  ligand = ECFP4 (183, 2048). Labels = % inhibition from duongly_mmc3
  Table S2 (raw values kept; observed range -12.5..191.3 recorded; no
  windsorization at this stage because targets are centered contrasts).
- Pair table: Q0B_MAPPING_AUDIT.json duongly_variant_records —
  ADMITTED single-point pairs only (65); multi-mutant/deletion/
  insertion/unknown rows excluded from the single-mutant estimand.
- Pair-level split (capacity surface): pairs stratified by parent
  60/20/20 (train/val/test of PAIRS; parents may appear in every
  split). Seeds {1,2,3}; SHA-256 keyed rng (x0_common.stable_rng);
  no Python hash().
- Target per pair (WT row w, variant v, ligand l):
  d_vl = y[v,l] - y[w,l] over ligands with finite labels in both rows;
  c_vl = d_vl - mean_l(d_vl)  (centered mutation contrast).

## Arms (identical budget; single seed screen first, 3 seeds on promotion)

1. unified_local  — UnifiedPotential on the full pocket one-hot (THE arm)
2. unified_global — GlobalPotential (20-dim pocket composition; global
   compression diagnostic)
3. ligand_only    — s=0 and b_P=0 (predicts zero centered contrast)
4. no_interaction — s=0, b_P+b_L heads present
5. family_shuffle — variant rows' protein features permuted within the
   SAME parent (WT fixed)
6. random_protein — protein features permuted across all rows
7. ligand_invariant_shift — s replaced by a per-row scalar (centered
   contrast = 0 by construction)
8. free_pairwise   — DIAGNOSTIC ONLY: antisymmetric h([P_wt,P_v,L]) -
   h([P_v,P_wt,L]); not integrable, never a production mechanism
9. free_ligand_pair — DIAGNOSTIC ONLY: ligand-pair predictor (used on
   the ligand-contrast estimand, reported, not gated)

All gradient-trained end-to-end (AdamW lr 1e-3, wd 1e-4, epochs 200,
batch 512); no ridge/pseudoinverse/closed-form/test-time gradients.

## Loss (frozen)

L = L_contrast + 1.0 * L_abs
- L_contrast = mean over pairs,ligands of (c_hat - c)^2 with c_hat from
  the SAME s_theta (centered_mutation_effect).
- L_abs = MSE(f(P,L), y_raw) with explicit b_P/b_L heads (nuisance
  separation; b_P/b_L cancel in the centered contrast by construction).
- free_pairwise trains on c directly (L_contrast only).

## Metrics and gates (frozen; cluster = parent, 2000-draw paired bootstrap)

- pair-level Spearman(c_hat, c) per pair; pair-mean Spearman/Pearson;
  dead-zone sign accuracy (dead zone = 10 % inhibition units); MSE;
  scale = OLS slope of c on c_hat (no intercept).
- CIIP-1A PASS: unified_local pair-mean Spearman >= 0.30 AND dead-zone
  sign accuracy >= 0.65 AND (unified_local - family_shuffle) >= 0.05
  AND (unified_local - ligand_only) >= 0.05 (bootstrap 2.5% lower
  bounds), on the TEST pairs (3-seed median).
- free_pairwise interpretation: if free_pairwise - unified_local > 0.10
  (pair-mean Spearman), report "pairwise signal exists but the
  integrable potential's expression is insufficient" — never promote
  free_pairwise to a mechanism.
- No threshold may be moved after results exist.

## Deliverables

RESULT.json, REPORT.md, commands.jsonl, tests (structure + data
contract), SHA-256 pins; history/task/ledger sync. Stage 2 (CIIP-1B,
held-out parents) is authorized ONLY after an interpretable 1A result.
