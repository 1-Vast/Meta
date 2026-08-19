# Stage Q2d-1d preregistration - span-restricted identifiable interaction (2026-08-19)

Corrected successor of Q2d-1c (STOP at oracle precheck; see
stageQ2d1c_identifiable_interaction_20260818/Q2D1C_REPORT.md). Frozen BEFORE
any Q2d-1d computation; SHA-256 in PREREGISTRATION_SHA256.txt. Synthetic
only; old stages remain read-only. Baseline unchanged: pipeline qualification
FAILED at Q2 (X0c); biological conclusion UNRESOLVED; B1/B2/C/D NOT
AUTHORIZED.

## Change vs Q2d-1c (with evidence)

Q2d-1c proved that the protein factor map has an unidentifiable component:
the train-row feature submatrix has rank 28 < 32 (PCA-32 over-bases the 80
train rows) and ~8.8% of the drawn A's energy lies in that null space; no
train-only estimator can recover it, which pulled the closed-form oracle to
dz 0.587-0.621 on protein-cold/double-cold for seed 1 (true weights: 0.968+).
The ligand map is fully identified (rank 48/48).

Fix: the truth protein factor map is drawn ENTIRELY WITHIN the train-row
span of the protein features. A_t = V_train^T @ C, where V_train is the
rank-28 right-singular basis of the train-row feature submatrix (computed
from features only, before truth generation, frozen with SHA-256), and
C (28x4) is Gaussian with QR-orthonormal columns scaled [1.0, 0.8, 0.6, 0.4]
(condition 2.5). The mechanism remains feature-conditioned (the protein
factor is a linear function of the pocket physicochemical features through
the observable subspace); the truth is now entirely identifiable from the
training design. M2 likewise projects its compressed block-sparse map onto
that span. Everything else is carried over unchanged.

## Carried over unchanged from Q2d-1c

- Protein truth features: first 32 principal components of the 510-dim
  per-position physicochemical descriptor (frozen projection over all 97
  rows, label-free, SHA-256).
- Ligand pool: 157 SMILES-resolved ligands. Ligand truth features: ECFP4
  (2048) through a FROZEN sparse random matrix W_L (2048x48, density 0.1,
  SHA-256-seeded), train-ligand-only standardization.
- Feature-smoothed double centring (row/column means removed through
  feature-linear projections fit on train; offsets are representable by the
  learner), sd over train cells = tau* = 1.0.
- Mechanisms: M1 matched linear low-rank bilinear (span-restricted A_t as
  above, B 48x4 QR scales [1,1,1,1]); M2 block-sparse local-group map in the
  pre-compression 510 space, PCA-compressed and span-projected; M3 shallow
  nonlinear tanh(bilinear/sqrt(4)); NC1 main-effect-only; NC2
  ID-random-factor. Noise N(0,1) per cell, frozen per (mechanism, seed).
- Truth weight streams stable_rng("stageQ2d1d", "truth", ...) disjoint from
  model-init streams (tested). Phase A interaction-only: y = tau*I + noise,
  no mu/pmain/lmain.
- Splits: protein-component cold = ALL rows whose kinase family is in
  {Tec, FGFR, LRRK, STE7, Src}; ligand-scaffold cold = frozen largest
  Bemis-Murcko clusters >= 25% of the resolved pool; surfaces protein-cold /
  ligand-cold / DOUBLE-COLD (800 each, primary gate surface); train =
  train rows x train ligands; checkpoint = best TRAIN loss at 300-step
  monitors (no val set).
- Arms: correct, ligand_only, additive_only, shuffled_protein,
  family_preserving_shuffle, random_protein, no_interaction_head,
  oracle_diagnostic; identical rows / init policy / minibatch order (rng
  keyed by truth seed, phase, restart - not by arm) / checkpoint rule.
- Budget: 6,000 steps, batch 1024, AdamW lr 5e-3, wd 1e-4, correct 8
  restarts, negatives 1 restart, xavier + 0.5-scaled init.
- Ladder: A identity endpoint interaction-only; B + sigmoid %-scale loss;
  C + 70% observed MCAR (frozen seed); D + true interval censoring (logit
  of rounded %, bounds outside 0/100) with frozen assertion censored_count
  > 0 aggregated over the level s 3 seeds; E + main effects (MAIN_SD=1.0)
  with shared-encoder competition. Value-level reproduction: each level L in
  B..E reruns with all added ingredients disabled and asserts bitwise
  equality of eval interaction outputs against level A s stored seed-0
  outputs.
- Gate (never moved): on DOUBLE-COLD, median over 3 truth seeds: correct
  arm Spearman >= 0.30, dead-zone sign accuracy >= 0.70, gap vs ligand_only
  >= 0.05; every negative arm fails the same gate; correct dz exceeds the
  best negative-arm dz by >= 0.03; protein-cold and ligand-cold reported
  with correct beating all controls.
- Feature-space oracle precheck BEFORE any training: closed-form train-only
  rank-4 fit (truncated SVD on the complete train grid + feature least
  squares; diagnostic only) must reach dz >= 0.70 on ALL THREE surfaces for
  every truth seed of M1; otherwise STOP and no representation comparison.

## Downstream authorization (unchanged)

Q2d-2 (representation matrix incl. KLIFS per-position ESM kept as
[position, embedding]; correct must beat ALL protein controls on
double-cold) only after Q2d-1d PASS. Q2d-3 (endpoint / within-protein
ligand delta / WT-mutant same-ligand delta ablations with gradient coverage
and conflict reporting; contrast targets supervise deployment quantities,
exchange-antisymmetric, identity-zero) only after Q2d-2 PASS. Saifudeen B1
(same-study primary; single mutants; exact matched WT; compatible
construct/substrate; same ligand; interval-censored + responsive-window
parallel estimands; held-out parent / pocket component / ligand scaffold;
correct arm improves its own prediction; cluster bootstrap lower bound above
ligand-only/shuffled/matched-wrong/random controls; endpoint never called
pK/Ki/Kd/DTA) requested ONLY after all synthetic gates PASS.

## Governance (unchanged)

SHA-256 seeds; no Python hash(); artifacts carry schema / prereg SHA /
input SHA; commands.jsonl appended; restricted data never committed;
research code only under tools/research/stageX_csc_signal/; model/ and
production scripts/ untouched before B1; implementation and report commits
separated.
