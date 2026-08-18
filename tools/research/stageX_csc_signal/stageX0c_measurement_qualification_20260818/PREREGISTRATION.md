# Stage X0c preregistration — measurement-pipeline qualification (corrected successor)

Frozen 2026-08-18, BEFORE reading any successor evaluation output and BEFORE
running any candidate model. This stage succeeds the original Stage X0,
which is ruled INVALID INSTRUMENT (see
../X0_INVALID_INSTRUMENT_VERDICT.md): the load-bearing distance-ratio gate is
a measurement-definition failure and cannot be repaired in place. The
original preregistration, thresholds, results and failed implementations are
preserved unchanged.

Inherited (re-verified per asset: hash, license, semantics): X0-D downloads
(Duong-Ly, Anastassiadis, Davis, PKIS2), the corrected WT-mutant pair table
and its diagnostic distance QC (../X0_PAIR_TABLE.json, ../X0_I2.json — kept
as diagnostics only, NOT as capability gates), KLIFS/UniProt/PubChem fetch
manifests (../x0_fetch_manifest.json).

## 0. Ordered gates (no gate before its predecessors pass)

Q0 variant-coordinate layer -> Q1 representation capability (probe
selectivity) -> Q2 fully synthetic planted harness -> Q3 biological panel
qualification (Saifudeen 2026) -> B1 same-study WT->mutant positive control
-> B2 localization -> C cold-protein interaction test -> D final DTA.
B1 requires Q0+Q1+Q2+Q3; C requires B1+B2. Stage S-W remains read-only.
Research code lives under tools/research/stageX_csc_signal/ only.

## 1. Q0 — variant-coordinate layer

### Q0-A external validation (ProteinGym)
- Source: ProteinGym DMS_substitutions (official GitHub repository, MIT
  license); downloads recorded with URL, access date, SHA-256; raw files not
  committed.
- Sample: up to 50,000 records sampled with SHA-256-derived seed
  'stageX0c/q0a/sample' (stable across processes), restricted to records
  whose mutant annotation parses as one or more single substitutions.
- Checks per record: (a) old-residue agreement — every annotated WT residue
  equals the reference sequence residue at the annotated coordinate;
  (b) mutation application — our generated mutant sequence equals the
  dataset-provided mutated_sequence; (c) multi-mutation handling — all
  substitutions of a ':'-joined mutant verified together;
  (d) deterministic serialization — canonical JSON of a VariantRecord is
  byte-stable across processes.
- FROZEN THRESHOLD: both (a) and (b) agreement >= 99.5% on the sample.
  Q0-A fails below 99.5%.

### Q0-B historical/construct mapping
Verified item-by-item with first-hand sources (recorded in the alias ledger
with DOI/URL and access date):
- BRAF V599E -> V600E: must document the 3-nt sequence difference between
  the historical reference (Davies 2002, Nature 417:949, doi:10.1038/nature00766)
  and the current canonical sequence (UniProt P15056, VAR_018629 Val600Glu),
  and must NOT be generalized to other proteins or assumed to be a
  start-methionine shift.
- PDGFRalpha: accession, species, isoform and construct offset verified per
  row (S1 lists GenBank NP_006197 with UniProt Q9DE49 = Danio rerio; human
  canonical is P16234). D842V construct range 668-1210 exceeds the canonical
  length (1089) and must be quarantined, not silently repaired.
- Duong-Ly all 76 variant rows typed with admission status + reason.
- KLIFS pocket numbering: 85 aligned positions validated against the KLIFS
  numbering scheme (KLIFS documentation + Kooistra et al. 2016); the
  gatekeeper position must map to KLIFS pocket index 45 for the known
  gatekeeper mutations (EGFR T790M, KIT T670I, RET V804L/M, FGFR1 V561M,
  SRC T341M, PDGFRA T674I, MAPK14 T106M, ALK L1196M, FGFR4 V550E/L).
- Davis WT/mutant/phosphorylation-state constructs: census from Davis
  MOESM3 with per-row construct semantics.
- FROZEN RULE: every listed item must have either a verified mapping or an
  explicit quarantine with reason; no silent substitution/offset/reference
  modification anywhere; hard fail on any old-residue mismatch.

## 2. Q1 — representation capability via probe + control-task selectivity

Distance ratios (../X0_I2.json) are diagnostics only. The primary Q1
instrument is restricted-capacity probe selectivity on variant-relevant
labels, Hewitt-Liang style but rebuilt for protein variants.

- Units: 65 admitted Duong-Ly single-point WT-mutant pairs (../X0_PAIR_TABLE.json).
- Representations: global_pooled_esm, pair_centered_local_esm,
  local_onehot_window, klifs_pocket, residue_identity_context,
  edit_descriptor, composition, random (SHA-256-seeded, dim-matched),
  family_id, parent_id, mutation_position_only, substitution_type_only.
- Probe inputs: the PAIR representation x(pair) = concat(x_WT, x_mutant) (or
  the pair-defined edit vector where that is the representation's definition).
- Tasks:
  T-A 'pocket membership': binary label 1 if the mutation maps into the
    85 KLIFS aligned pocket positions (from Q0-B mapping), else 0.
  T-B 'substitution physicochemical class': 6 frozen classes of the (old,new)
    substitution by amino-acid physicochemical group (aliphatic, aromatic,
    polar-uncharged, positive, negative, special). Frozen class map in
    q1.py; labels come from the pair table only.
- Control tasks: same label MARGINALS as the real task, assigned by a stable
  seed within parent clusters (so the control keeps the parent-split
  structure); a probe passing a control task only memorizes marginals.
- Probes: linear (least squares via gradient descent, no closed-form
  deployment) and capacity-restricted MLP (hidden 8, Tanh). Probe capacity
  curve: linear vs MLP-8; random-label curve reported.
- Splits: (a) leave-one-parent-out (21 folds, cluster unit = parent);
  (b) substitution-type held-out for T-B (two frozen classes held out);
  (c) family held-out using KLIFS family groups. Evaluation labels never
  select probes; hyperparameters fixed in advance.
- Metrics: per-representation task accuracy, control accuracy,
  selectivity = task - control, calibration (Brier / log-loss vs random)
  per-parent accuracies, parent-cluster bootstrap (2000 draws, seed
  20260820), corruption controls (label shuffle, feature shuffle).
- FROZEN Q1 PASS: at least one representation outside
  {edit_descriptor, random} achieves selectivity >= 0.10 (point estimate)
  with cluster-bootstrap 2.5% lower bound > 0 on task T-A under the
  leave-one-parent-out split; AND the edit_descriptor's measured selectivity
  is reported honestly (no presumption of zero). Q1 result is information
  readability only; it certifies nothing about downstream causal use.

## 3. Q2 — fully synthetic planted-signal harness

Only the real observation graph, missingness, censoring topology,
protein/pocket groups, ligand/scaffold groups, degree/imbalance structure are
reused. Real endpoint values never enter the labels.

- Generative process: z(p,l) = mu + a(p) + b(l) + tau* * I(p,l) + eps,
  all draws from SHA-256-derived RNGs; y% = 100*sigmoid(z) (primary
  realization mirrors the real Duong-Ly pattern: continuous values, no hard
  floor, off-scale values reported); an emulated floor-clamp realization
  exists ONLY for censoring-machinery checks and the floor-imputation
  control. Censoring applied after full latent generation.
- Grid: tau* = SD(tau*I)/noise_sd in {0, 0.125, 0.25, 0.5, 1.0, 2.0} x
  rank R in {1, 4, 16} (dense locality), plus a sparse-pocket locality
  variant (3 driving pocket positions) at (tau*=1.0, R=4). The grid must
  bracket the graph's approximate detection threshold
  tau*_det ~ sqrt(R*(n_p+n_l)/N_obs); the computed value is recorded.
- Interaction truth forms: low-rank dense (primary), sparse aligned-pocket,
  local mutation x ligand-fragment, and a nonlinear-but-low-capacity form
  (elementwise tanh of the dense bilinear). All forms are functions of the
  protein pocket features and ligand ECFP features, so neither modality
  alone recovers them; rank, locality, SNR, variance, direction balance,
  and correlation with main effects are recorded per realization.
- Arms (identical optimizer AdamW, width 32, budget, initialization policy
  via fixed seed, early stopping on train/val only): ligand_only (protein
  input zeroed), correct_protein, shuffled_protein, family_preserving_shuffle,
  random_protein, no_interaction_head, free_target_id (NON-TRANSFERABLE
  upper bound, never a candidate), oracle_protein (the true latent protein
  factors of the planted interaction — the achievable ceiling).
- Model emits protein main, ligand main and interaction components
  explicitly; the fitted interaction component is evaluated against the
  planted interaction truth. Additionally an ANOVA projection operator Pi
  (least squares fit of mu + a(p) + b(l) on TRAINING cells only, applied
  identically to prediction and truth) cross-checks the explicit head.
- Split: strict train/val/eval blocks by parent component AND ligand
  scaffold; eval = eval_parents x eval_ligands cells; no parent construct,
  compound or scaffold crosses blocks; all hyperparameters/early stopping
  use train/val only.
- Metrics: interaction Spearman (primary), Pearson, interaction MSE,
  dead-zone sign accuracy, scale/slope recovery, intercept, tau* recovery
  curve, rank/locality recovery, per-seed results (3 frozen seeds),
  component-cluster bootstrap (2000 draws, seed 20260820).
- FROZEN Q2 PASS (at tau*=1.0, R=4, dense, on held-out eval): correct arm
  interaction Spearman >= 0.30 AND dead-zone sign accuracy >= 0.70 AND
  sign_accuracy(correct) - sign_accuracy(ligand_only) >= 0.05, on the
  median of the 3 seeds; AND every negative control fails by construction:
  tau*=0 shows no signal in any arm; no_interaction_head fails; label
  permutation destroys recovery; protein permutation destroys
  protein-conditioned recovery; main-effect-only data never classified as
  interaction PASS; floor imputation shows its bias. Other grid points are
  sensitivity only and never select the gate.

## 4. I6 — production dataflow integrity suite

Assertions frozen (each mapped to a test): contrast antisymmetry; identity
pair strictly zero; reference-term sign flip; ligand/protein/row/column ID
mapping; endpoint direction and units; old-residue consistency; stable seed
across OS processes; every regularizer finite nonzero gradient; gradient
coverage for every trainable branch; dead-branch capture; permutation
controls destroy target information; matched arms share row ids, splits,
labels, masks and censoring bounds exactly; eval labels never enter inputs,
reference or normalization; normalization statistics from train only;
interval bounds ordered and directionally correct after every transform;
sign-only target direction correct; checkpoint selection never reads eval
labels; cluster bootstrap resamples parent/pocket components (row bootstrap
forbidden for load-bearing intervals); duplicate ligands / same-parent
mutants / same-scaffold compounds not counted as independent; train/eval CSC
orientation identical; generator round-trip; planted truth bitwise
recomputable; restricted data never packaged into commits.

## 5. Q3 — Saifudeen 2026 panel qualification

First-hand audit of Saifudeen et al. (2026) Nat. Biotechnol.
(doi:10.1038/s41587-026-03090-8) and the KIRHub supplement: 92 clinical
inhibitors (86 approved), 409 WT, 349 variants/fusions (311 mutants + 38
fusions), 1 uM, per-kinase Km ATP, duplicate measurements, Supplementary
Tables 2/4/11/13. Deliverables: pairability census (exact matched WT,
same construct background, same substrate, same ATP protocol, same endpoint,
single/multi mutant, fusion, no matched WT, ambiguous construct,
censored/saturated, responsive-window status), license manifest (expected
CC BY-NC-ND 4.0: local use only; no repackaged derivative matrices in Git;
commit downloaders, code, hashes, schemas and value-free summaries only).
Saifudeen is a functional-inhibition positive control, never called pK/pIC50.

## 6. Stop rules and reporting

- Q0-A < 99.5% => stage result UNRESOLVED(coordinate-mapping), stop before Q1.
- Q1 gate fails => report representation incapability for the failing
  representation set; the distance QC is descriptive only.
- Q2 gate fails at the frozen point => the harness is unqualified; no
  real-data negative may be interpreted; stop before Q3 panel experiments
  (Q3 census work is data governance and may proceed).
- No threshold is edited after results are read. Single seeds are for
  debugging only; load-bearing conclusions use the frozen multi-seed +
  cluster-bootstrap protocol. meta_test remains sealed; ridge/pseudoinverse
  are audit tools only and never the deployment mechanism; no test-time
  query gradients.

## 7. Frozen seeds

- split seed 20260818; bootstrap seed 20260820; Q0-A sampling seed from
  SHA-256('stageX0c/q0a/sample'); Q2 seeds SHA-256('stageX0c/q2/'+tau*+'/'+rank)
  with seed index 0/1/2; all other RNGs SHA-256-derived from named strings.
Python hash() is forbidden everywhere.

## 8. Artifacts (all machine-readable, with schema version, code commit,
## input SHAs and this preregistration SHA)

VARIANT_RECORD_SCHEMA.json, Q0A_PROTEINGYM_VALIDATION.json, Q0B_MAPPING_AUDIT.json,
Q0B_ALIAS_LEDGER.md, Q0B_KLIFS_CENSUS.json, Q1_SELECTIVITY.json, Q2_PLANTED.json,
I6_TEST_REPORT.json, Q3_SAIFUDEEN_CENSUS.json, RESULT.json, REPORT.md, commands.jsonl.
