# Stage Q2d-1c preregistration - identifiable feature-conditioned interaction (2026-08-19)

Corrected successor of Q2d-1b (STOP at oracle precheck; see
stageQ2d1b_feature_interaction_20260818/Q2D1B_REPORT.md). Frozen BEFORE any
Q2d-1c computation; SHA-256 recorded in PREREGISTRATION_SHA256.txt.
Synthetic-only; old stages remain read-only. Baseline unchanged: pipeline
qualification FAILED at Q2 (X0c); biological conclusion UNRESOLVED;
B1/B2/C/D NOT AUTHORIZED.

## Changes vs Q2d-1b (with evidence)

1. PROTEIN FEATURE COMPRESSION. Q2d-1b oracle failure attributed (diagnostic
   _diag_centring.py): the 510-dim protein feature map is unidentified on the
   row complement (train row span rank <= 80). Fix: protein truth features
   are the first 32 principal components of the 510-dim per-position
   physicochemical descriptor, computed over ALL 97 rows BEFORE truth
   generation (unsupervised, label-free feature compression; projection
   matrix stored with SHA-256). 32 <= train-row rank, so the factor map
   A (32x4) is identified on the train span and cold rows project into the
   same space.
2. LIGAND POOL RESTRICTED. Unresolved-SMILES ligands carry hash-fallback
   ECFP bits and dominated the cold scaffold clusters. Fix: the ligand pool
   is the 157 SMILES-resolved ligands only. Ligand truth features: ECFP4
   (2048-bit) projected through a FROZEN sparse random matrix W_L (2048x48,
   density 0.1, SHA-256-seeded, stored), standardized with TRAIN-ligand-only
   per-dimension mean/std. Scaffold clusters (Bemis-Murcko) over resolved
   ligands; largest clusters held out until >= 25% of the pool is cold.
3. FEATURE-SMOOTHED DOUBLE CENTRING. ID centring injects structure a feature
   bilinear cannot represent (Q2d-1b diagnostic). Fix: per-row and per-column
   means of I_raw over TRAIN cells are removed through feature-linear
   projections fit on train (row means regressed on P_t over train rows;
   column means regressed on L_t over train ligands), so every subtracted
   offset is representable by the learner. Then sd over TRAIN cells = tau*
   (= 1.0). The double-centring mandate is honored (row and column structure
   removed) while identifiability is preserved; the ligand-only/protein-only
   arms and the additive-only / shuffled / random arms remain the empirical
   guards that no ID structure carries the signal.

## Carried over unchanged from Q2d-1b (same SHA-anchored definitions)

- Mechanisms M1 (matched linear low-rank bilinear, A 32x4, B 48x4, QR-
  orthonormal columns, scales [1.0, 0.8, 0.6, 0.4], condition 2.5),
  M2 (sparse/local group: factor k on pocket positions 21(k-1)+1..21k of the
  PRE-COMPRESSION 85-position descriptor, block-sparse A 510x4 compressed to
  32 dims by the frozen PCA), M3 (shallow nonlinear: tanh(bilinear/sqrt(4))),
  NC1 (main-effect-only), NC2 (ID-random-factor).
- Noise N(0,1) per cell, frozen stream per (mechanism, truth seed).
- Truth weights from stable_rng("stageQ2d1c","truth",mechanism,seed) -
  disjoint from model-initialization streams (tested).
- Phase A interaction-only: y = tau * I(P,L) + noise. NO mu, NO pmain/lmain.
- Splits: protein-component cold = ALL rows whose kinase family is in
  {Tec, FGFR, LRRK, STE7, Src} (families of the five old eval parents);
  ligand-scaffold cold = frozen largest Bemis-Murcko clusters >= 25% of the
  resolved pool; surfaces protein-cold (cold rows x train ligands, 800),
  ligand-cold (train rows x cold ligands, 800), DOUBLE-COLD (cold rows x
  cold ligands, 800) = PRIMARY gate surface; train = train rows x train
  ligands; checkpoint = best TRAIN loss at 300-step monitors (no val set).
- Arms: correct, ligand_only, additive_only, shuffled_protein,
  family_preserving_shuffle, random_protein (capacity-matched Gaussian),
  no_interaction_head (scale frozen 0, requires_grad=False),
  oracle_diagnostic (P = P_t A_t, diagnostic only). All arms share identical
  rows, initialization policy, minibatch order (rng keyed by truth seed,
  phase, restart - NOT by arm) and checkpoint rule.
- Budget: 6,000 steps, batch 1024, AdamW lr 5e-3, wd 1e-4, correct arm 8
  restarts (init 0-7), negative arms 1 restart, xavier + 0.5-scaled init.
- Ladder: A identity endpoint interaction-only; B + sigmoid %-scale loss;
  C + 70% observed MCAR (frozen seed, counts recorded); D + true interval
  censoring (logit of rounded %, bounds outside 0/100) with frozen assertion
  censored_count > 0 aggregated over the level s 3 seeds; E + main effects
  (MAIN_SD=1.0) and shared-encoder competition. Value-level reproduction:
  each level L in B..E reruns with all added ingredients disabled and
  asserts bitwise equality of eval interaction outputs against level A s
  stored seed-0 outputs.
- Gate (never moved): on DOUBLE-COLD, median over 3 truth seeds: correct
  arm Spearman >= 0.30, dead-zone sign accuracy >= 0.70, gap vs ligand_only
  >= 0.05; every negative arm fails the same gate; correct dz exceeds the
  best negative-arm dz by >= 0.03; protein-cold and ligand-cold surfaces
  reported with correct beating all controls.
- Feature-space oracle precheck BEFORE any training: rank-4 ALS (train-only,
  closed form, diagnostic only) must reach dz >= 0.70 on ALL THREE surfaces
  for every truth seed of M1; otherwise STOP and no representation
  comparison is authorized.

## Downstream authorization (unchanged)

Q2d-2 (representation matrix incl. KLIFS per-position ESM kept as
[position, embedding], no premature pooling; correct must beat ALL protein
controls on double-cold) only after Q2d-1c PASS. Q2d-3 (endpoint /
within-protein ligand delta / WT-mutant same-ligand delta ablations with
gradient coverage and conflict reporting; contrast targets supervise
deployment quantities, exchange-antisymmetric, identity-zero) only after
Q2d-2 PASS. Saifudeen B1 (same-study primary; single mutants; exact matched
WT; compatible construct/substrate; same ligand; interval-censored +
responsive-window parallel estimands; held-out parent / pocket component /
ligand scaffold; correct arm improves its own prediction; cluster bootstrap
lower bound above ligand-only/shuffled/matched-wrong/random controls;
functional-inhibition endpoint never called pK/Ki/Kd/DTA) requested ONLY
after all synthetic gates PASS.

## Governance (unchanged)

SHA-256 seeds; no Python hash(); artifacts carry schema / prereg SHA /
input SHA; commands.jsonl appended; restricted data never committed;
research code only under tools/research/stageX_csc_signal/; model/ and
production scripts/ untouched before B1; implementation and report commits
separated.
