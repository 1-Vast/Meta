# Stage U preregistration — protein × local chemical transformation interaction

Frozen **2026-08-19, before any Stage U census statistic, any interaction
variance test and any trained-arm evaluation metric was read.** The only
external knowledge allowed into this file is the already-admitted evidence
listed in `task.md` / `history.md` / Stage S / Stage T artifacts, and the
provenance counts recorded below. Nothing in this file may be changed after the
first U0 admission statistic or the first U2 evaluation metric is read; a
post-hoc threshold change voids the stage.

## 0. Estimand and why this is not Stage S and not Stage T

Stage S rejected whole-molecule pairs + global protein FiLM; Stage T built a
first true-MMP census and a pooled-protein discriminator but never produced a
promotion decision, and its exact transformation key omitted the shared core.
This stage tests the final admissible form of the Phase-1 question.

    tau      = (shared core, R_a -> R_b, attachment context,
                stereochemistry, charge change)
    delta_y(t,tau) = mu_tau + delta(t,tau) + noise
    D(tau,t1,t2)   = delta_y(t1,tau) - delta_y(t2,tau)

`D` cancels the target affinity level **and** the generic transformation effect
`mu_tau` algebraically. A model scored on `D` cannot score on target level, on
generic medicinal chemistry, or on a target-identity key. The question is
whether protein sequence-derived region tokens can predict `D` on unseen protein
components better than zero, better than a global protein summary, better than a
cross-component shuffled protein, and better than a similarity-matched wrong
protein taken from the same held-out population.

This is an admissibility test, not a DTA performance claim.

## 1. Governance

* Governed BindingDB-Ki corpus `main_v0`, double-cold split `v1`, mounted
  through the **physically isolated split view**
  (`dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1_views`). The
  sealed confirmation split artifact is not present on that surface.
* The development-validation split is used for nothing: not training, not
  selection, not thresholds, not early stopping, not a census. A parsed-AST
  test fails the suite if its name appears as a string constant in any stage
  module.
* Frozen `scripts/internal_validation.py` partition of `meta_train`
  (`PARTITION_SEED = 20260818`): 227 fit / 31 internal-validation components.
  Training reads fit only; reporting reads internal-validation only.
* Stable SHA-256 seeds throughout (`scripts.qpsmp_data.stable_seed`). Python
  `hash()` is forbidden and tested for by AST parse.
* Labels never touch feature construction, split construction, MMP definition,
  transformation-key granularity or hyperparameter selection. Rules below are
  label-blind.

### Disclosed isolation exception for the U0 provenance step

U0 must read pre-aggregation measurements, which live in a single all-label
artifact (`dataset/processed/crossed_interaction/bindingdb_202608/
exact_labels.jsonl.gz`, SHA-256 `6a44a151…c823dcc`, matching the corpus
manifest's `labels_sha256`; and `metadata_projection.jsonl.gz`, matching
`projection_sha256`). The allow-list of `source_row_id`s is derived from the
physically isolated split view, `meta_train` only; the artifact is then streamed
and a row is retained in the same pass only if its id is allow-listed.
Non-allow-listed values are never bound to any retained structure. This is
logical exclusion after parsing, strictly weaker than the model path's physical
isolation, and it is disclosed as such in `U0_RELIABILITY.json`. The sealed and
development-validation rows are excluded by the allow-list.

## 2. U0 — measurement reliability and transformation-graph audit

### 2.1 Aggregation verification

Quote and verify the corpus manifest cleaning block (verbatim in artifact):
`within_panel_aggregation: "median"`, then
`cross_panel_pair_aggregation: "equal-panel median"`, endpoint
`"exact positive uncensored Ki"`, transform `pKi = 9 - log10(Ki[nM])`,
`panel_id_role: "BindingDB assay proxy"`.

Provenance levels:
* L1: rows sharing `(panel_id, assay protocol_sha256)` — same-document repeat;
* L2: rows sharing `panel_id`, differing protocol;
* L3: rows of one `(target, ligand)` across `panel_id`s.

Report counts, robust dispersion and classical dispersion. L2 is expected to be
non-identifiable; if so, the artifact says so instead of inventing a number.

**Same-panel delta-pK reliability.** For two cells sharing a governed panel, the
panel-specific offset cancels in the difference. The frozen point estimate is

    sigma2_same = 2 * within_assay_variance

where `within_assay_variance` is the L1 pooled residual variance around group
medians, and the clustered bootstrap 95% interval over repeated-measure groups
is reported next to it. Cross-panel `sigma2_cross = 2 * (within + between)`
using L3 residual variance, with the same disclosure.

**Not a benchmark MSE floor.** Repeated-measure variance on the selected subset
is reported only as supervision reliability. It is never quoted as an MSE lower
bound.

### 2.2 MMP construction

* `rdkit.Chem.rdMMPA.FragmentMol`, Hussain–Rea, single cut (`minCuts=1,
  maxCuts=1`, max cut bonds 30, default acyclic single-bond pattern), isomeric
  SMILES throughout. No ad-hoc SMILES string surgery.
* Core = larger fragment by heavy-atom count; ties by canonical SMILES sort.
* Pair admissibility: same target, same Ki endpoint, identical core isomeric
  SMILES (including `[*:1]`), different R groups, matching attachment context.
* Attachment environment saved per transformation: element, aromaticity,
  ring membership, degree, formal charge and hybridization of the core atom
  bearing the cut dummy. Stereochemistry of core and both R groups is retained
  in isomeric SMILES and flagged separately.
* Canonical direction = canonical SMILES sort of `(R_a, R_b)`, structure only;
  `delta_y = y(ligand carrying R_b) - y(ligand carrying R_a)`. The inverse
  transformation has reversed R order and must map to `-delta_y` (tested).
* Deduplication of a repeated `(target, core, R_a, R_b)` keeps the lower cell
  index, never a label value.

**Transformation key granularities, frozen:**

1. `exact key` = SHA-256 of
   `core_isomeric_smiles | repr(full attachment context) |
   R_a_isomeric | >> | R_b_isomeric`;
2. `coarse key` = SHA-256 of
   `core_stereo_stripped | element | aromatic |
   R_a_stereo_stripped | >> | R_b_stereo_stripped`.

Both are structure-only functions; no label may enter either.

### 2.3 Primary bank and strata

* Primary bank: same target, exact positive uncensored Ki endpoint, and the two
  cells share at least one identical governed `panel_id`.
* `S1_same_panel_single`: shared panel, both `panel_count == 1`.
* `S2_same_panel_multi`: shared panel, at least one `panel_count > 1`.
* `S3_cross_panel`: no shared panel. Built and reported only as a separate weak
  noise stratum; never pooled into primary statistics.

### 2.4 Bipartite evidence graph and diagnostics

Construct `G = (protein targets, exact transformation keys, observations)`.
Report: observations, targets, components; exact/coarse key counts; target and
component degree per key; repeated keys spanning >=3 targets and >=3
components; connected components and degree concentration; target ×
transformation incidence effective rank, stable rank, condition number and
singular-value spectrum; fit/internal scaffold and transformation overlap;
deployment coverage C_k for k={1,2,3,5} on the frozen nested episode banks
(query labels not read). All of these are explicitly labeled
**empirical sufficient-richness diagnostics**, not a proof of persistent
excitation.

### 2.5 Frozen U0 admission gate — verbatim

1. same-panel fit observations >= **2,000**;
2. fit targets >= **50**;
3. at least **30** exact keys each spanning >= **3 targets** and
   >= **3 components**;
4. internal same-panel observations >= **300**;
5. internal components >= **10**;
6. no domination by a few high-degree nodes:
   * top-1 exact-key observation share <= **0.05**;
   * top-10 exact-key observation share <= **0.20**;
   * top-1 target observation share <= **0.25**;
   * top-5 target observation share <= **0.75**;
   * top-1 component observation share <= **0.25**;
   * top-5 component observation share <= **0.75**.

Thresholds 1–3 and 6 are evaluated on the exact key over same-panel fit
observations. If **any** threshold fails, stop immediately, write a negative
`REPORT.md`, and train no neural model.

## 3. U1 — interaction variance (runs only if U0 passes)

Uses fit components only for the gate; internal components are reported as
consistency only and never gate training.

1. Aggregate same-panel `delta_y` per `(exact key, target)` with the median
   (expected to be one observation because the exact key contains the core).
2. Restrict to exact keys with >=2 target effects (fit).
3. Per key compute `SS_tau = sum_t (dy_t - mean_tau)^2`, `df = k_tau - 1`;
   pool `MS_effect = sum SS_tau / sum df`.
4. Noise reference `sigma2_noise` = U0 `sigma2_same` for one observation; if a
   key has multiple observations for one target the median correction is
   applied, but this is expected to be rare.
5. Estimand `theta = MS_effect - sigma2_noise`.

**Hierarchical bootstrap (frozen).** 2,000 draws, seed 20260820.
(a) Interaction bootstrap: resample exact keys with replacement and protein
components with replacement; each (key,target) effect's multiplicity is the
product of its key draw count and its component draw count; recompute the
weighted pooled between-key MS. (b) Noise bootstrap: resample U0 L1 repeated
measurement groups with replacement (group = the repeated rows of one cell,
panel, protocol), recompute pooled residual variance, and double it. (c) The
two independent draws are paired by position, `theta_b = MS_b - sigma2_noise_b`.

**Frozen U1 gate.** The 2.5th percentile of `theta_b` must be **> 0**.
Report additionally the ratio `MS_effect / sigma2_noise`, the number of keys,
effects and components, and an internal-consistency figure on the internal
partition (descriptive only). If U1 fails, stop: write the negative report and
train no neural model. No larger or deeper network is a rescue.

## 4. U2 — minimum transformation-conditioned protein-region operator

Runs only if U0 and U1 both pass.

### 4.1 Representation contract

* Edit branch: a structured, label-blind transformation token built from the
  saved core / R_a / R_b descriptors, attachment environment, charge and stereo
  flags (counts + 256-bit folded Morgan fingerprints of core, R_a, R_b).
* Protein branch: the governed ESM-2 150M protein bank's 128 **ordered**
  residue/region tokens (640-d), projected to model width, with fixed
  sinusoidal slot position encoding and mask handling.
* The transformation token actively queries protein-region tokens through
  exactly two multi-head cross-attention layers (4 heads, model width 128,
  FFN 256, dropout 0.1). The query is the edit token; keys/values are the
  ordered region tokens. No pooled protein summary, no target embedding, no
  target index, no component ID, no assay/document ID reaches the local
  operator.
* Output `R(tau, p)` is a scalar. Prediction is
  `D_hat(tau,p1,p2) = R(tau,p1) - R(tau,p2)`, which is antisymmetric under
  p1/p2 exchange, exactly zero at p1=p2, and protein-cycle consistent for every
  parameter setting.

This is an adaptation inspired by System Identification, Conditional Similarity
and branch/trunk Operator Learning. It is **not** a reproduction of Meta-SysId,
CSN, DeepONet, ANP or PEARL, and no artifact may claim otherwise.

### 4.2 Matched arms — same fixed budget

| arm | description |
|---|---|
| A `A_zero` | constant response; `D_hat` identically 0 |
| B `B_global` | global ESM pooled protein summary + edit token (Stage S negative reference; global pooling is permitted for this arm only) |
| C `C_local` | **candidate**: transformation-conditioned local protein-region operator |
| D `D_local_shuffled` | same local operator, trained on stable cross-component shuffled protein tokens |
| E `E_mean_tokens` | same local operator, but the protein tokens are the target-independent masked mean of fit-component region tokens; structurally `D_hat = 0` |
| F `F_label_shuffled` | same local operator and correct protein, trained on within-transformation permuted `D` labels |

Wrong protein is **evaluation only**. No correct-vs-wrong training loss exists
in any arm. Training is ordinary AdamW forward/backward: no ridge, no
pseudoinverse, no closed form, no target-ID embedding, no query-label
adaptation.

### 4.3 Data banks and splits (all label-blind)

* `delta_y` effects and `D` rows are formed only within one exact key.
* Train: cross-component `D` rows whose both targets are fit components
  (cross-component = the two targets have different CD-HIT40 components).
* `fit_unsampled`: a frozen 10% sample of fit cross-component `D` rows, by
  stable seed, never used for training.
* `internal_repeated`: internal `D` rows whose exact key is also present in the
  fit bank (transformation keys repeatable). **Primary evaluation surface.**
* `internal_disjoint`: internal `D` rows whose coarse key is absent from the
  fit bank (transformation-family/scaffold cold). Evaluated if >=100 rows.
* `internal_all`: all internal same-panel `D` rows.
* All related protein components are on one side (fit vs internal); every
  instance of a coarse transformation family is held together for the disjoint
  surface.

### 4.4 Wrong and shuffled protein controls

* Candidate matching population is **the same population as the recipient**:
  internal recipients draw only from internal targets, fit-unsampled recipients
  draw only from fit targets.
* Different CD-HIT40 component, most similar admissible target by cosine on
  frozen ESM-2 150M pooled embeddings, ties broken by target id sort. Only the
  protein input changes; tau and labels are fixed.
* Shuffled protein (paired substitution and arm D): stable offset permutation
  over the same population's targets requiring different CD-HIT40 components.
* Arm D is trained on a stable cross-component permutation over fit targets
  only.

### 4.5 Training hyperparameters (frozen)

Seed `20260821` (single screen); confirmation seeds `20260822`, `20260823`.
Steps `4000`; batch `256`; AdamW lr `3e-4`, weight decay `1e-4`, cosine to 0;
Huber loss delta `1.0`; gradient clip `5.0`; row sampling weights
`1 / sqrt(deg_exact_key)` inside the fit train bank. No checkpoint selection of
any kind: every arm trains the fixed number of steps and the final parameters
are evaluated, because the reporting population is the internal partition and
selecting on it would leak. If CUDA is available it is used; device is
recorded.

### 4.6 Frozen U2 success conditions — verbatim

All are evaluated on the primary `internal_repeated` bank unless stated; paired
contrasts use identical rows. The two-way cluster bootstrap samples components
with replacement and transformation keys with replacement independently; a
row's multiplicity is the mean of its two endpoint component draw counts times
its key draw count. 2,000 draws, seed `20260820`. Effective independent units =
`min(#components, #keys)` are reported with every interval.

1. C_local minus A_zero Pearson >= **+0.05** AND C_local minus B_global Pearson
   >= **+0.05**;
2. both differences have clustered-bootstrap 95% lower bounds **> 0**;
3. C_local correct-input minus C_local shuffled-input Pearson >= **+0.05** with
   lower bound **> 0**;
4. C_local correct-input minus C_local matched-wrong-input Pearson >= **+0.05**
   with lower bound **> 0**;
5. C_local vs D_local_shuffled: MSE delta < 0 AND Spearman delta > 0 AND
   CI delta > 0 AND sign-accuracy delta > 0;
6. C_local minus F_label_shuffled Pearson >= **+0.05** with lower bound > 0
   (label shuffle destroys the effect);
7. protein-induced shift (C correct input minus C shuffled input) has Pearson
   correlation with truth >= **+0.10** (Stage P/S failure mode);
8. on `fit_unsampled`, C correct-input minus C shuffled-input Pearson
   >= **+0.05** with lower bound > 0 (a shuffled protein may not reproduce the
   fit-unsampled gain);
9. leave-one-out influence: no single transformation key and no single protein
   component accounts for more than **50%** of the C-vs-shuffled-input Pearson
   effect;
10. on `internal_disjoint` (if >=100 rows), C correct-input minus C
    shuffled-input Pearson must be **>= 0** (transformation/scaffold-cold does
    not reverse); if the surface has <100 rows the gate is recorded
    `not_evaluable` and the route cannot pass.

A single seed may reject; a single seed may not confirm. Only if the single
screen passes all ten conditions do three fixed seeds (20260821/22/23) run all
six arms, and confirmation additionally requires conditions 1–10 to hold in
every one of the three seeds.

## 5. Stop rules

* U0 fail -> stop before U1/U2; negative report; no neural model.
* U1 fail -> stop before U2; negative report; no neural model.
* U2 single-seed fail -> **the current BindingDB protocol route for
  protein-conditioned SAR interaction is formally closed**; no rescue by ANP,
  PEARL, MSA, Cartesian, conformer, larger models or moved thresholds.
* U2 single-seed pass but multi-seed fail -> recorded as NOT CONFIRMED; nothing
  promoted.
* U2 multi-seed pass -> no production integration; write only the next-stage
  proposal described in the active contract (support evidence as target-state
  identification, protein sequence as prior, ANP-style attention over
  transformation observations, k=1 conservative, k>=2 update).

## 6. Verification required before training

Structural tests must pass before any arm trains:
MMP canonical direction/inverse; attachment/stereo/charge preservation;
deterministic banks; independence from `PYTHONHASHSEED`; physical meta-test
seal; no cross split/component/transformation leakage; antisymmetry/identity/
cycle consistency; gradient coverage; no target-ID bypass in the local
operator; query labels never enter inputs, splits or checkpoint selection.
Commands, environment, raw rows and machine-readable JSON are retained in this
stage directory.

## 7. Artifacts

`U0_RELIABILITY.json`, `U0_CENSUS.json`, `U0_PROVENANCE_meta_train.jsonl.gz`,
`U1_VARIANCE.json`, `runs/<arm>/RUN.json` and `*.rows.json`, `RESULT.json`,
`CONFIRMATION.json` (if run), `REPORT.md`, `COMMANDS.md`, `ENVIRONMENT.json`.
The preregistration SHA-256 is recorded inside `RESULT.json`.
