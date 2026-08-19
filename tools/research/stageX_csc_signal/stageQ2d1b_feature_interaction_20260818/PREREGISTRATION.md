# Stage Q2d-1b preregistration - feature-conditioned interaction qualification (2026-08-18)

Frozen BEFORE any Q2d-1b computation (except the Q2d-1 forensic tests, which
only read old artifacts). SHA-256 recorded in PREREGISTRATION_SHA256.txt.
Synthetic-only; no external biological labels. Old stages (X0c, Q2c, Q2d-1)
are read-only inputs. Baseline: pipeline qualification FAILED at Q2 (X0c);
biological conclusion UNRESOLVED; B1/B2/C/D NOT AUTHORIZED.

## Mechanism under test

A transferable statistical interaction field of local protein
physicochemical features x ligand substructure descriptors, learned as a
gradient-trained low-rank bilinear map. No atomic-contact, conformational or
binding-free-energy claims.

## Truth features (frozen)

- Protein truth features P_t (per row, 510-dim): KLIFS-aligned per-position
  physicochemical descriptor = 85 pocket positions x 6 frozen amino-acid
  scales (hydrophobicity, volume, charge at pH7, polarity, H-bond donor
  count, H-bond acceptor count). Each scale is z-scored over the 20 amino
  acids (constant normalization, no data governance). Rows without a KLIFS
  pocket (LRRK2) get zero features and a frozen row mask.
- Ligand truth features L_t (per ligand, 64-dim): ECFP4 (2048-bit, Morgan
  radius 2) projected through a FROZEN sparse random matrix W_L (2048x64,
  density 0.1, SHA-256-seeded, drawn once, stored), then standardized with
  TRAIN-ligand-only per-dimension mean/std (train-governed; recorded in the
  artifact). Ligands with all-zero ECFP4 (unresolved SMILES) get zero
  features and a frozen ligand mask.

## Truth mechanisms (each: effective rank, condition number, parameter count)

- M1 matched linear low-rank bilinear (pipeline minimal positive control):
  I_raw = (P_t A_t) . (L_t B_t); A_t: 510x4, B_t: 64x4; columns
  orthonormalized (QR) with column scales [1.0, 0.8, 0.6, 0.4]; condition
  number 2.5; truth parameters 2,296; observations per parameter
  = |train cells| / 2,296 (recorded).
- M2 sparse/local group interaction (local-mechanism positive control):
  4 factors; factor k acts on pocket positions 21(k-1)+1..21k (contiguous
  quarters of the 85 positions) only: A_t block-sparse 510x4 (126 nonzero
  per column), same column scales; B_t dense 64x4 orthonormal; truth
  parameters 760; obs/param recorded.
- M3 shallow nonlinear feature-conditioned (mild model mismatch):
  I_raw = tanh( ((P_t A_t) . (L_t B_t)) / sqrt(4) ) with the M1 weights.
- NC1 main-effect-only negative control: z = mu + pm + lm + noise, I = 0.
- NC2 ID-random-factor negative control: the old q2.generate truth
  (one-hot pocket x random U,V) - the matched feature learner must FAIL on
  cold surfaces, reproducing the old structural result in the clean
  protocol.

## Standardization and noise (frozen)

- I = I_raw double-centred: subtract per-row mean and per-column mean
  computed on TRAIN cells only, then scaled so sd over TRAIN cells = tau*
  (= 1.0). Eval surfaces are never used for centring or scaling.
- noise per cell ~ N(0, 1.0), frozen stream per (mechanism, truth seed).
- Phase A truth: y = tau * I(P,L) + noise. NO mu, NO pmain/lmain.
- Truth weights are drawn from stable_rng keyed by (stageQ2d1b, truth,
  mechanism, seed) - disjoint from every model-initialization stream
  (asserted by a test).

## Splits and evaluation surfaces (frozen)

- protein-component cold: hold out ALL rows whose kinase family is in the
  frozen set {BTK, FGFR, LRRK, MAPK, SRC-family} (families of the five old
  eval parents). train rows = everything else.
- ligand-scaffold cold: Bemis-Murcko scaffolds from resolved SMILES
  (unresolved ligands form their own cluster); hold out the frozen list of
  largest scaffold clusters until >= 25% of ligands are cold (list recorded
  in the artifact). train ligands = the rest.
- Surfaces: protein-cold cells = cold rows x train ligands (frozen subsample
  800); ligand-cold cells = train rows x cold ligands (800); DOUBLE-COLD =
  cold rows x cold ligands (frozen subsample 800) - the PRIMARY gate
  surface. train cells = train rows x train ligands.
- Checkpoint rule: best TRAIN loss at 300-step monitors. There is NO val
  set; no checkpoint ever sees a cold-surface label.

## Arms (all share identical rows, init policy, minibatch order, checkpoint)

correct (P_t), ligand_only (P=0), additive_only (no interaction
parameters), shuffled_protein (row-permuted P_t), family_preserving_shuffle
(permute rows within family), random_protein (frozen Gaussian 510-dim,
capacity-matched), no_interaction_head (inter_scale frozen 0,
requires_grad=False), oracle_diagnostic (P = P_t A_t, 4-dim; diagnostic
only, never deployment). Minibatch order rng keyed by (truth seed, phase,
restart) - identical across arms. Budget: 6,000 steps, batch 1024, AdamW
lr 5e-3, wd 1e-4, correct arm 8 restarts (init seeds 0-7), negative arms
1 restart, xavier + 0.5-scaled Wp/Wl init.

## Ladder A-E (each level reproduces the previous level value-exactly)

- A: identity endpoint, interaction-only truth, learner without mu/pm/lm
  biases: y_hat = scale * ((p A) . (l B) + inter_bias); loss MSE(z).
- B: + sigmoid endpoint: y% = 100*sigmoid(z); loss on % scale.
- C: + missingness: 70% observed cells, frozen MCAR seed; observed/missing
  counts recorded.
- D: + censoring: true interval censoring (logit of rounded %, bounds
  outside 0/100); asserts censored_count > 0 (frozen assertion).
- E: + main effects and competition: truth gains mu + pm + lm (MAIN_SD=1.0);
  learner gains shared linear encoders for pm/lm with per-row/per-ligand ID
  biases.
- Value-level reproduction: for each level L in B..E, the runner re-runs
  level L with the added ingredient disabled and asserts bitwise equality
  (tolerance 0) of the eval interaction outputs against level A's stored
  seed-0 outputs.

## Frozen gate (never moved retroactively)

On the DOUBLE-COLD surface, median over 3 truth seeds: correct arm
Spearman >= 0.30 AND dead-zone sign accuracy >= 0.70 AND gap vs ligand_only
>= 0.05 (frozen thresholds from the original gate). Every negative arm must
FAIL the same gate. Correct dz must exceed the best negative-arm dz by at
least 0.03. Protein-cold and ligand-cold surfaces are reported and the
correct arm must beat all controls there (reported, not gated).

## Feature-space oracle precheck (before ANY training; closed form is
## diagnostic only and is never a deployment candidate)

For M1/M2/M3 and each truth seed: rank-4 alternating least squares fit on
TRAIN cells of I using P_t/L_t (closed form, train-only), evaluated on all
three cold surfaces. Frozen rule: if the M1 oracle does not reach dz >= 0.70
on all three surfaces, the stage STOPS before training and no representation
comparison is authorized (rule 3.7).

## Downstream authorization

Q2d-2 (representation matrix: aligned one-hot, z-scales, KLIFS per-position
local ESM kept as [position, embedding] with no premature global pooling,
ESM+z-scales, global pooled ESM negative, residue permutation, shuffled,
random capacity-matched - same head/data/budget; correct must beat ALL
protein controls on double-cold) runs ONLY after Q2d-1b PASS. Q2d-3
(endpoint loss / within-protein ligand delta / WT-mutant same-ligand delta
ablations with gradient-coverage and gradient-conflict reporting; contrast
targets must supervise deployment quantities, preserve exchange
antisymmetry and identity-zero) runs ONLY after Q2d-2 PASS. Saifudeen B1
(same-study primary; single mutants; exact matched WT; compatible
construct/substrate; same ligand; interval-censored + responsive-window
parallel estimands; held-out parent, pocket component, ligand scaffold;
correct arm must improve its own prediction, not damage wrong arms; cluster
bootstrap lower bound above ligand-only/shuffled/matched-wrong/random
controls; functional-inhibition endpoint never called pK/Ki/Kd/DTA) is
requested ONLY after all synthetic gates PASS.

## Governance

SHA-256 seeds; no Python hash(); artifacts carry schema / prereg SHA / input
SHA; commands.jsonl appended; restricted data never committed; research code
only under tools/research/stageX_csc_signal/; model/ and production
scripts/ untouched before B1; implementation and report commits separated.
