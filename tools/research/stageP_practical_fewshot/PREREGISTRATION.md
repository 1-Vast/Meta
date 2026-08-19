# Stage P preregistration — three-layer evaluation hierarchy (2026-08-19)

Frozen BEFORE any P-line computation. This document only re-organizes
evaluation into three independent layers (P1 / P2 / M1) per the 2026-08-19
re-adjudication. It authorizes no single worst-case split to carry all
conclusions. It does not move any gate in the existing frozen Q2d chain
(stageX_csc_signal), does not modify old stage artifacts, and does not
touch model/ or scripts/.

## 1. Layer definitions

### P1 — Practical few-shot (new-target lead optimization)
- Allowed information: target amino-acid sequence; homology/MSA/conservation
  features; predicted monomer structure / pocket prior; assay context when
  present in the source data (never guessed); support rows from the SAME
  target, related chemical series allowed (no scaffold split); external
  pretrained representations (ablated and reported separately).
- Forbidden: query labels in adaptation or model selection; test-set
  selection of structure or hyperparameters; target-ID shortcuts; assay
  context invented when absent.
- Support sizes: k in {5, 10, 20, 40} primary; k in {1, 2, 3} extreme
  stress; k=0 reported for continuity (defined below).
- Split: protein(target)-cold. Component-balanced target split from the
  governed BindingDB-Ki corpus; support/query cells from the same target,
  ligand-unique.
- Metrics: MSE, RMSE, CI (concordance index), Spearman, Pearson, centered
  MSE; stratification: support-query Tanimoto bands, Bemis-Murcko scaffold
  novelty, activity-cliff presence (|d(pKi)| large on nearest neighbour).
- Arms (bake-off): current admitted baseline (QPSMP/BPSF), ligand-only,
  fixed Morgan/Tanimoto kNN, ordinary support fine-tuning, first-order
  MAML, AdaMBind-style task sampling+loss, CNP (support encoder), FS-CAP-
  style support encoder, ActFound-style within-task pairwise. Adoption of
  AdaMBind/MAML is decided by measured results only.

### P2 — Novel-target screening
- Allowed: same information set as P1; screening labels if present in the
  source (binary active/inactive or %inhibition).
- Forbidden: same as P1.
- Support sizes: k=0 and k=5.
- Split: protein(target)-cold; ligand novelty stratification (scaffold
  held-out targets' ligands vs training ligands).
- Metrics: regression MSE/CI/Spearman always; if screening labels exist:
  EF(1%,5%,10%), BEDROC(alpha=20), PR-AUC reported in addition, never
  instead of regression metrics.

### M1 — Fundamental stress test (mechanism)
- Split: protein-component-cold AND ligand-scaffold-cold (double-cold);
  k in {0,1,2,3,5}.
- Arms: correct protein; ligand-only; shuffled protein; family-preserving
  shuffle; similarity-matched wrong protein; residue permutation;
  capacity-matched random protein; no-interaction; component-level
  bootstrap (paired, target/component).
- Claims: M1 decides ONLY mechanism claims. M1 failure never negates P1/P2
  performance value, and P1/P2 gains are never attributed to
  protein-conditioned interaction without an M1 (or preregistered real
  positive-control B1) pass.

## 2. Cross-layer invariants

- One frozen episode bank (target, component, support cell ids, query cell
  ids per k and per draw) shared by every arm; k=0 episodes use the SAME
  rows with zero support and the SAME model code path, the few-shot module
  removed only if that removal is itself the explicit experimental
  variable.
- Fairness: same split, same support/query rows, same training budget
  (steps x batch), same backbone or explicitly reported parameter delta,
  same checkpoint rule (frozen per bake-off), single-seed screening then
  >= 3 fixed seeds, paired target/component bootstrap with intervals.
- Seeds: SHA-256 stable seeds; Python hash() banned; RNG streams keyed per
  (stage, phase, seed, restart) and shared across arms.

## 3. Dataset governance (Phase 3 ordering)

- BindingDB Ki, Davis Kd, KIBA composite, PKIS/PKIS2, and any other legal
  DTA set are trained independently or with dataset-specific output heads;
  each dataset keeps its own train/validation/test governance and
  normalization. Cross-dataset pretraining reports: shared representation,
  dataset-specific head, external-data ablation, zero-shot transfer,
  fine-tuned transfer. One dataset's test result never selects another
  dataset's structure or hyperparameters.
- Saifudeen functional inhibition panel: protein-conditioned selectivity
  positive control only; never called pK/Ki/Kd DTA.
- KIBA is a composite score, never described as physical affinity.
- Davis censored/blank values are interval-censored, never exact constants.

## 4. Promotion / stop rules (unchanged from re-adjudication)

- P-line promotion: paired improvement vs matched baselines; MSE/RMSE
  improved; CI/Spearman not significantly degraded; improvement not
  explainable by exact ligand recall; >= 3 fixed seeds direction-
  consistent; component-level intervals support the main claim.
- M-line promotion: correct protein beats ligand-only and ALL protein
  controls; improvement occurs in the correct branch; holds on double-cold
  or the preregistered real positive control; not a target-ID/family/assay
  shortcut; label permutation and residue permutation destroy the effect.
- Stop: Q2d-1e + the single finite diagnostic fail -> close the low-rank
  bilinear learner; local representation indistinguishable from
  shuffled/matched-wrong -> close that representation family; MSE-only
  gains with ranking/centered-error degradation -> calibration shortcut;
  practical few-shot valid but mechanism gate failed -> keep the
  performance model without interaction claims.

## 5. Execution order

1. P-line data census: BindingDB-Ki governed corpus (exists), Davis Kd
   (dataset/raw/dta/davis.tab), KIBA (dataset/raw/dta/kiba.tab), PKIS2
   supplements (stageX downloads), Saifudeen (stageX downloads).
2. Freeze P1 bake-off preregistration (episode bank, budget, checkpoint
   rule, metrics, promotion gates) as the next document.
3. Q2d-1d/1e adjudication continues in stageX_csc_signal; representation
   comparison (Phase 5) and any interaction-trunk change wait for the Q2d
   verdict. Training innovation (Phase 6) must show positive effect with
   the interaction trunk fixed in its ablation.
