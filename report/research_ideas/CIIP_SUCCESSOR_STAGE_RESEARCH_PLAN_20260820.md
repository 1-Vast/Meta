# CIIP Successor-Stage Research Plan (DRAFT, pre-preregistration)

Date: 2026-08-20
Author role: independent senior research scientist (analysis only)
Status: DRAFT. No training, no production change, no stage directory has been
created by this document. Execution requires (1) approval of this plan,
(2) freezing the preregistration in Section 9 and recording its SHA-256,
(3) the S0 audit passing. Frozen CIIP-1A artifacts are read-only inputs.

Scope discipline: endpoint is functional percent inhibition everywhere; never
Ki, Kd, pK, or DeltaDeltaG. The Duong-Ly meta-test analogue (the sealed
BindingDB meta_test) is never read. Duong-Ly test-pair labels are used only
for final evaluation of frozen arms, never for fitting, normalization,
feature construction, model selection, or retrieval.

---

## 1. Executive summary

CIIP-1A returned ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED. The new read-only,
train/val-only diagnostics in Section 4 change the interpretation of that
failure:

1. The centered estimand c is dominated by a **parent-level (pocket-level)
   ligand-response component**: same-parent cross-mutation c profiles agree at
   median Spearman 0.442 (0.406 after WT-profile residualization) vs 0.036
   across parents. A trivial, non-learned parent-profile predictor
   (leave-one-pair-out mean of same-parent c) already reaches median Spearman
   0.579 / R2 0.326 on train pairs — above every CIIP-1A arm on test
   (Spearman ~0.33, R2 <= 0.13).
2. The **mutation-specific residual** (c minus parent profile) carries ~2/3 of
   the energy but is idiosyncratic per mutation (no positive cross-mutation
   sharing) and is the only component the correct-vs-random-window contrast can
   attribute. At 9 test pairs / 6 parent clusters, the stage had no power for
   it (bootstrap width ~0.49 R2).
3. The random-window negative control is **miscalibrated for mutation
   attribution**: contextual propagation (non-site delta 0.0575) plus
   parent-level information make "correct ≈ random" the EXPECTED outcome under
   a parent-level-signal model even if a mutation-specific increment existed.

Conclusion: CIIP-1A failed the attribution test it ran, but that test was
near-unidentifiable as designed. The sharp, still-open question decomposes
into three propositions (Section 5): A (mutation overall effect), B
(mutation-erased protein context), C (protein-conditioned ligand ranking).
The highest-information, lowest-risk next step is **CIIP-S1: signal
decomposition and estimand ladder** (Section 9): a read-only-plus-small-
training stage that measures A, B, and C on the same legal Duong-Ly data and
the existing ESM cache with a mutation-erasure counterfactual arm, before any
model-class decision.

---

## 2. What is proven vs not evaluated (ledger)

### Proven / established (frozen artifacts)

- P1. KLIFS one-hot mutation bottleneck: 38/65 pairs have ΔP = 0; constant
  predictions structurally forced by antisymmetry (STAGE1_COLLAPSE_AUDIT).
- P2. Absolute-vs-contrast objective dominance at stage-1 (R_g ~1081,
  C_g ~ -0.016: dominance proven, destructive conflict NOT established);
  val contrast MSE never left the zero floor.
- P3. Oracle local ESM restores nonconstant outputs (9/9 covered test pairs,
  6/6 parents) but pair-mean R2 <= 0.05: restored variation != restored
  signal (2x2 diagnostic).
- P4. Objective effect on ESM is absent (-0.017, CI crosses 0): the
  centered-only objective does not handicap the ESM arm (2x2 diagnostic).
- P5. Formal verdict ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED under frozen rules;
  correct-random R2 = -0.1217, parent bootstrap [-0.4569, +0.0327]; all
  correct-vs-floor/control CIs cross zero (CONTROL_ADJUDICATION.json).
- P6. ESM mutation sensitivity at the site is a representation fact:
  mean site delta 4.0111, radius-6 1.1605, correct-vs-random delta 0.5310 vs
  0.0267, 49/49 pairs larger at the site. NOT ligand-conditioned evidence.
- P7. Contextual propagation exists: non-site context delta 0.05749,
  full-sequence 0.07392; erasure control exact (49/49 identical inputs, max
  delta 0.0). A random window is not a mutation-information-free null.
- P8. Q1 (frozen): pair_centered_local_esm probe selectivity +0.189
  [0.033, 0.363] under leave-one-parent-out — capability on a control task;
  klifs_pocket -0.086 ns. Edit descriptors (position-only +0.110,
  substitution-type-only +0.209) pass as SHORTCUTS, not biological signal.
- P9. Coverage selection bias is structural: all 16 uncovered pairs fail at
  pos > 1020 (ESM length cap) and are exactly ALK/MET/LRRK2/TEK C-terminal
  mutations; the missing set is higher-variance (DATA2X2 coverage_bias).

### Established this session (exploratory, train+val covered pairs only,
### read-only, NOT preregistered; must be re-derived inside CIIP-S1)

- E1. Energy split of d: mutation main effect = 10.5% of mean d^2 energy;
  ligand-specific centered part = 89.5%. Median |mean_l d| = 3.29 units
  (IQR 1.73-5.85, max 12.97) vs median sd(c) ~13 units.
- E2. Same-parent cross-mutation c consistency: median Spearman 0.442
  (IQR 0.337-0.670, n=60 pair-combos); WT-profile-residualized 0.406;
  different-parent baseline 0.036 (IQR [-0.080, 0.173]).
- E3. Parent-profile LOPO ceiling (non-learned): median Spearman 0.579
  (IQR 0.412-0.748, n=29), R2 median 0.326; WT-residualized target 0.428.
- E4. Ligand-global main effect on c (leave-pair-out per-ligand mean):
  R2 = 0.0602 — the shared structure is parent-specific, not ligand-global.
- E5. Parent-residualized cross-mutation consistency: median -0.280
  (mechanically biased negative by leave-one-out construction, but no strong
  positive sharing survives): the mutation-specific residual is idiosyncratic.
- E6. Label-range audit (full matrix, frozen DATA1A label_stats): 4087/17751
  cells (23.0%) outside [0,100]; range -12.5..191.3.

### Not evaluated / unresolved

- N1. Context-only predictive value of ESM states (the propagation audit used
  no labels): does mutation-erased context predict c? NOT EVALUATED.
- N2. Mutation overall effect (proposition A, target d or mean_l d):
  never fitted under the CIIP contract.
- N3. Mutation-specific ligand-conditioned increment (parent-residualized c):
  never isolated; CIIP-1A's correct-vs-random contrast is an underpowered,
  miscalibrated proxy for it.
- N4. Any parent-disjoint (cold-target-family) evaluation of any CIIP object.
- N5. Transfer of anything measured here to Ki/Kd/pK or BindingDB: NOT
  EVALUATED and NOT AUTHORIZED.
- N6. Multi-seed stability of any ESM arm (single seed everywhere).
- N7. Whether a deployable (mutation-coordinate-free) representation can
  capture the parent-level component (E3 ceiling was non-learned; no trained
  deployable arm has been measured against it).

---

## 3. Root-cause ranking of the CIIP-1A failure

Ranked by contribution to the observed null, with evidence grade.

### R1 — Estimand/attribution mismatch: the stage tested the hardest, smallest
### component with a contrast that cannot see the dominant one. (PRIMARY)

The centered target c mixes a parent-level component (predictable WITHOUT
reading the mutation; E2/E3) and a mutation-specific residual (E5: no
cross-mutation sharing). The frozen primary contrast (correct window minus
random window) can only register the mutation-specific increment, because both
arms share the parent and both windows encode parent information. The
parent-level component — the only component with a demonstrated achievable
ceiling (E3: Spearman ~0.58) — was never attributed to any arm because no
parent-profile reference arm existed. Grade: NEW direct evidence (E2-E5),
consistent with P5/P7.

### R2 — Miscalibrated primary negative control (random local window)

Random windows (a) carry propagated mutation deltas (P7: 0.0575 non-site) and
(b) trivially encode parent identity (they are windows of the same parent's
sequence). Under R1, "correct ≈ random" is the expected outcome even if a
mutation-specific increment exists. The control is valid for the narrow claim
"no site-dependence detected"; it is invalid as "no protein-conditioned
signal". Grade: PROVEN direction (P7), magnitude of miscalibration
quantifiable only via the erasure arm (N1).

### R3 — Statistical power for the attributable contrast

9 test pairs, 6 parent clusters, single seed; the paired arm-difference CI
width on R2 is ~0.49 (P5). Per-pair Spearman noise (~1/sqrt(180) ≈ 0.075)
permits detecting a per-pair effect, but the BETWEEN-ARM attribution at the
parent-cluster level is severely underpowered for increments < ~0.2 R2.
Grade: PROVEN (CI widths in CONTROL_ADJUDICATION.json).

### R4 — Split semantics: pair-cold, not parent/protein-cold

All 6 test parents appear in train (verified from DATA1A/DATA2X2: covered
train/val parents include ABL1, KIT, EGFR, FGFR4, RET, TEK). The stage
measured within-family mutation-pair generalization. This is not a flaw of
the stage's own claim, but it (a) enables parent-profile riding (R1), and
(b) says nothing about cold-target transfer (N4). Grade: structural fact.

### R5 — Objective dominance (historical stage-1 cause; NOT the control-stage
### cause)

PROVEN for the joint-loss stage (P2); the control stage ran centered-only for
all arms and P4 shows the objective does not explain the ESM arm's null.
Retained here to prevent re-blending objectives in any successor. Grade:
PROVEN, closed as an explanation of the current verdict.

### R6 — Representation bottleneck (KLIFS one-hot): PROVEN (P1), superseded by
### ESM; does not explain the ESM-arm null. Closed.

### R7 — Endpoint semantics (functional % inhibition)

Single-concentration functional readout; 23% of cells outside [0,100] (E6);
bounded-interpretation values treated as unbounded linear targets; ATP
competition and construct effects survive centering. This plausibly depresses
absolute R2 and inflates heavy-tail MSE, but it acts IDENTICALLY on all arms,
so it cannot explain the attribution failure by itself; it does cap achievable
R2 and distorts linear metrics. Grade: plausible contributor, bounded role.

### R8 — Coverage selection bias (P9)

The covered 49 pairs exclude C-terminal mutations of large kinases; claims are
restricted to the covered population. Not a cause of the null (all arms share
coverage); a generality limit. Grade: structural fact.

### R9 — Oracle representation

The tested representation needs mutation coordinates; it was never a
deployment candidate. This constrains what a PASS could have authorized; it
did not cause the FAIL. Grade: design boundary.

### Explicitly rejected root causes

- "ESM encodes nothing about the mutation": falsified by P6.
- "Protein-conditioned interaction is biologically absent": not tested; the
  stage is estimator-bound, and E2/E3 show a large parent-conditioned
  component exists in the LABELS regardless of representation.
- "Potential capacity insufficient": not implicated (P3, collapse audit §6).

---

## 4. New read-only diagnostics (this session; train+val covered pairs only)

Method notes: all quantities use ONLY the 40 covered train+val pairs
(DATA1A/DATA2X2); no test label was read; no model was fitted; keyed rng
(numpy default_rng with fixed seeds) for the different-parent baseline.
These numbers are exploratory and must be re-derived inside CIIP-S1's
preregistered S0/S1 audit before gate use.

| diagnostic | value | interpretation |
|---|---|---|
| main-effect energy share of d | 10.5% | centering discards only ~10% of contrast energy |
| median |mean_l d| | 3.29 units | mutation overall effect is small vs c sd ~13 |
| same-parent cross-mutation corr(c_i, c_j) | median Spearman 0.442 | parent-level component dominates c |
| same-parent, WT-residualized | 0.406 | not a shared-WT-profile artifact |
| different-parent baseline | 0.036 | component is parent-specific |
| parent-profile LOPO predictor | Spearman 0.579 / R2 0.326 | non-learned ceiling above every CIIP-1A arm |
| ligand-global leave-pair-out | R2 0.060 | ligand-only shortcuts are weak on c |
| parent-residualized cross-mutation | median -0.28 (LOO-biased) | mutation-specific residual is idiosyncratic |

Consequences for design:

- D1. Any future claim of "protein-conditioned interaction" on Duong-Ly MUST
  state which component it attributes: parent-level (large, achievable,
  mutation-coordinate-free) or mutation-specific (small, hard, oracle-bound).
- D2. The correct null for site-attribution is the **mutation-erasure
  counterfactual** (identical inputs at the site; validated at delta=0 in P7),
  not a random window.
- D3. The correct ceiling reference for c prediction is the cross-fitted
  parent-profile baseline, not the zero floor alone.
- D4. The mutation overall effect (proposition A) is a separate, small,
  possibly more mutation-specific-predictable target that the current
  estimand deliberately discards.

---

## 5. Three propositions, kept distinct

- A. Mutation-site representation predicts the mutation OVERALL effect
  (target: d, or mean_l d; includes the ligand-invariant shift).
- B. Mutation-ERASED protein context provides protein-state information
  (parent/pocket state) predictive of ligand-specific response.
  Every B branch is labelled: "mutation-free but not necessarily
  biologically causal".
- C. Protein representation changes ligand ranking / ligand-specific response.
  C_total: predicts c. C_sharp: predicts the parent-residualized c (the
  mutation-specific increment). CIIP-1A partially tested C_total with an
  oracle representation and a miscalibrated attribution contrast; it did not
  test C_sharp at usable power, and did not test A or B at all.

### Experiment matrix (representations x estimands)

Estimand ladder (all within-parent pairs; targets from the same legal labels):

| id | target | proposition | notes |
|---|---|---|---|
| T0 | d (raw WT->variant difference) | A (+C_total) | includes main effect |
| T0m | mean_l d (per-pair scalar) | A (pure) | mutation severity scalar |
| T1 | c (centered contrast) | C_total | current CIIP-1A estimand |
| T2 | c - parent-profile (cross-fitted) | C_sharp | mutation-specific residual |
| T3 | ligand ranking within pair (rank loss on c) | C_total | scale-free robust form |

Representation arms (features; identical potential/training where applicable):

| id | representation | deployable? | proposition addressed |
|---|---|---|---|
| F1 | mutation-site local ESM (oracle, radius-6) | NO | A, C |
| F2 | mutation-erased local/context ESM (site replaced by X in BOTH rows; identical inputs) | n/a (counterfactual) | strict zero-mutation-info null; B (in regression form) |
| F3 | full-sequence ESM mean (cached) | YES | B, C_total |
| F4 | KLIFS pocket one-hot | YES (structurally limited) | B, C_total |
| F5 | family-preserving random representation (within-parent keyed permutation) | control | shortcut audit |
| F6 | random local window (non-site, propagation-carrying) | control | context-propagation value; NOT a mutation-free null |
| F7 | ligand-only (protein zeroed) | control | zero floor |
| F8 | protein-invariant constant shift | control | excludes protein-main-effect shortcut by construction |
| F9 | parent-profile predictor (cross-fitted mean of same-parent c) | ceiling reference (uses labels of other pairs) | ceiling for T1; residualizer for T2 |
| F10 | optional external features (conservation/MSA/PSS) | conditional | only if data+license verified; UniRef snapshot currently absent locally |

Model forms: Form-1 = antisymmetric integrable potential
g = s(Pv,L) - s(Pw,L) (as frozen; required for any potential-compatible
claim). Form-2 = direct regression g(P_features, L) -> target (information
ceiling; F2/F9 only make sense here, since erased Form-1 inputs are
identical and yield g = 0 by construction — that structural zero is itself
the strict-null evidence). Mixing forms is explicit: Form-2 arms are
ceiling/diagnostic, never production mechanisms.

Decisive contrasts (pre-registered, paired, parent-cluster bootstrap):

- A-test: F1 vs F2-erased on T0/T0m (Form-2): does site information predict
  the overall mutation effect beyond erased context?
- B-test: F2-erased (context) vs F7/F8 floors on T1 (Form-2): does
  mutation-free context beat the floors? Flag: "mutation-free but not
  necessarily biologically causal".
- C_total-test: F1 (Form-1) vs F9 ceiling on T1: how much of the achievable
  parent-level signal does the potential capture?
- C_sharp-test: F1 vs F2-erased on T2 (Form-2): the mutation-specific
  increment — the sharp residue of the CIIP question.
- Deployability-test: F3/F4 vs F9 on T1 (Form-2): can a mutation-coordinate-
  free representation reach the parent-profile ceiling?

---

## 6. Estimator attribution review (each: main-cause? / evidence / cheap
## validation / fix if confirmed / still "interaction"?)

### Q1. Did centering over-remove transferable signal?

1. Main cause? NO for energy (keeps 89.5%), YES for predictability structure:
   the removed 10.5% (mean d) is the component most likely to be
   mutation-specific-predictable (it is a direct function of the edit).
2. Evidence: E1 (10.5% energy; median |mean d| 3.29 vs sd(c) 13).
3. Cheap validation: S0 computes T0m predictability from edit descriptors and
   erased context (train-only, read-only linear probes).
4. Fix: add T0/T0m to the ladder (planned); never merge with T1 into one
   scalar target.
5. Still interaction? The main effect is a protein main effect of the
   mutation; reporting it as "protein main effect of edit" is legitimate, but
   it is NOT ligand-conditioned interaction unless shown ligand-dependent.

### Q2. Is percent inhibition a valid linear difference target?

1. Main cause? Partial for absolute metric values; not for the attribution
   null (acts on all arms equally).
2. Evidence: E6 (23% out of [0,100]); functional assay semantics
   (Duong-Ly; Cheng-Prusoff logic says single-concentration inhibition mixes
   potency and mechanism).
3. Cheap validation: S0 reports heavy-tail diagnostics (skew/kurtosis,
   dead-zone composition) and rank-vs-linear metric agreement on existing
   frozen arm outputs (read-only re-scoring of CONTROL_RESULT.json).
4. Fix: T3 rank estimand as co-primary; robust regression (Huber) as a
   sensitivity arm; never relabel endpoint.
5. Still interaction? Yes — rank-based within-panel contrasts are still
   ligand-conditioned if they pass the same controls.

### Q3. logit(p)/rank/within-panel contrast/interval-censoring?

1. Main cause? The wrong transform alone cannot explain the null (arms share
   the target), but the transform determines the achievable ceiling.
2. Evidence: E6 shows logit(y/100) is UNDEFINED for 23% of cells (y>100 or
   y<0); logit requires clipping with explicit justification — undesirable.
   No censoring annotations exist in Duong-Ly Table S2 -> interval-censored
   formulation is not identifiable from current data.
3. Cheap validation: recompute frozen-arm metrics under rank target
   (read-only); report Spearman stability.
4. Fix: adopt T3 (rank) co-primary; keep raw-centered MSE as scale reporter;
   skip logit (invalid domain) and skip interval-censoring (no bounds data).
5. Still interaction? Yes, unchanged semantics.

### Q4. Should mutation overall effect and ligand-specific interaction be two heads?

1. Main cause? Not of the null, but the single-head centered design hid the
   decomposition that made R1 invisible.
2. Evidence: E1-E3; eSIG-Net (peer-reviewed) uses a joint objective exactly
   to avoid a discrepancy-only collapse.
3. Cheap validation: re-score frozen arms on T0m (read-only) to see whether
   any arm accidentally predicts the main effect.
4. Fix: two explicit heads (T0m head + T1 head) with disjoint reporting;
   the T1 head is the interaction claim; the T0m head is a nuisance/severity
   reporter.
5. Still interaction? Only the T1/T2 head is interaction; the T0m head is
   main effect — keep names separate.

### Q5. Predict protein main effect first, then residual interaction?

1. Main cause? The ABSENCE of a cross-fitted parent-profile nuisance made R1
   invisible; presence of the absolute head in stage-1 caused dominance (P2).
2. Evidence: E3 (parent profile predicts 1/3 of c variance); P2.
3. Cheap validation: S0 computes cross-fitted F9 and the T2 residuals
   (train-only, read-only) — done exploratorily (E5), redo preregistered.
4. Fix: cross-fitted parent-profile nuisance INSIDE the estimand definition
   (T2), fit on train parents only for val/test evaluation; never let the
   nuisance see the evaluated pair's labels.
5. Still interaction? The T2 residual is the sharpest interaction claim; the
   removed part is reported separately as family-level signal.

### Q6. Absolute vs centered gradient competition?

1. Main cause? PROVEN for stage-1 (P2); NOT operative in the control stage
   (all arms centered-only; P4).
2. Evidence: R_g ~1081, C_g ~ -0.02 (dominance without proven opposition);
   P4 (objective effect on ESM absent).
3. Cheap validation: none needed — closed by P4; if any future arm uses a
   joint loss, re-run the step-1 gradient audit first (read-only).
4. Fix: keep objectives on disjoint parameter subsets or sequential stages;
   if a joint loss is ever reintroduced, pre-register a gradient-dominance
   gate (R_g < 10) as a continuation condition.
5. Still interaction? N/A (optimization issue).

### Q7. Coverage selection bias in the 49-pair subset?

1. Main cause? No for the null (shared across arms); YES for generality.
2. Evidence: P9 (all missing pairs pos > 1020; ALK/MET/LRRK2/TEK; higher
   variance).
3. Cheap validation: S0 coverage table (position distribution, kinase length,
   pocket proximity) — read-only.
4. Fix: claims restricted to "covered population"; optional long-construct
   arm (chunked ESM with position bookkeeping) as a coverage-repair
   diagnostic, clearly labelled non-comparable budget.
5. Still interaction? N/A.

### Q8. Is random-window an over-strong/incorrect negative control?

1. Main cause? YES, jointly with R1, for the ATTRIBUTION failure.
2. Evidence: P7 (propagation), E2/E3 (parent-level component accessible from
   any window of the same parent).
3. Cheap validation: none further needed for direction; magnitude via F2
   erasure arm (N1).
4. Fix: replace the primary site-attribution null with F2 (erasure); keep F6
   only as a "context-propagation value" measurement, renamed and
   reinterpreted.
5. Still interaction? The erasure-null design preserves the interaction
   semantics of the correct arm.

### Q9. Statistical power?

1. Main cause? YES for the attributable contrast (R3).
2. Evidence: CI widths in CONTROL_ADJUDICATION.json; 6 clusters; single seed.
3. Cheap validation: S0 power computation: per-pair Spearman se ~0.075;
   paired-arm differences; parent-cluster bootstrap behavior under the
   observed parent-effect heterogeneity (read-only simulation from train
   labels, keyed rng).
4. Fix: power statement frozen in preregistration; UNRESOLVED is a
   pre-registered legitimate outcome for C_sharp; no post-hoc metric
   shopping; multi-seed only after single-seed structure+controls pass.
5. Still interaction? N/A.

### Q10. Hidden leakage via same parent / same ligand / same assay?

1. Main cause? For a cold-target claim: disqualifying. For the stage's own
   within-family claim: part of the design, but it enables parent-profile
   riding (R1) and shared-WT-row correlations (quantified: WT-residualized
   consistency stays 0.406, so shared WT is NOT the driver).
2. Evidence: split tables (R4); E2 WT-residualization; single-platform panel
   (assay identity constant -> assay batch is not a WITHIN-panel confound but
   blocks cross-assay claims).
3. Cheap validation: S0 leakage audit: WT-row sharing graph, ligand overlap
   matrix, parent x split incidence — read-only.
4. Fix: parent-disjoint evaluation lane (CIIP-1B semantics) remains gated;
   within-panel stage must report parent-cluster CIs and per-parent effects
   (already required).
5. Still interaction? Yes, but claims stay within-panel until the
   parent-disjoint lane is separately authorized and passed.

### Q11. Does Duong-Ly functional inhibition transfer to Ki/Kd/pK DTA?

1. Main cause? Not a cause of the null; a hard boundary on interpretation.
2. Evidence: endpoint definitions (functional % inhibition at fixed
   concentration vs equilibrium affinity); Cheng-Prusoff dependence on
   ATP/substrate; Nelen 2025 (matched pairs robust across assays; absolute
   values not comparable); Kalliokoski 2013 (mixed IC50 comparability
   limits); the project's own BOUNDARY_20260817_NIGHT (cold-target level wall
   on BindingDB-Ki).
3. Cheap validation: none that would authorize bridging; any bridge requires
   its own qualification stage (assay-semantics mapping, per-endpoint
   analysis) — currently NOT AUTHORIZED.
4. Fix: none within this stage; keep endpoints in separate namespaces
   forever.
5. Still interaction? The Duong-Ly findings are interaction claims about a
   functional kinase panel only.

---

## 7. Literature synthesis and applicability boundaries

Ledger notation: [PR] peer-reviewed; [PP] preprint/theoretical; suitability
for CURRENT data (Duong-Ly 97x183 % inhibition, ESM cache, no structures, no
MSA, no external labels) marked YES/PARTIAL/NO.

### 7.1 Main effect vs interaction decomposition; matched pairs; assay noise

- [PR] Pahikkala et al. 2015, KronRLS, Brief. Bioinform. 16(2):325 —
  pairwise kernel DTA with explicit cold-start settings S1-S4; canonical
  split-semantics reference. Suitability: YES (split design).
  https://academic.oup.com/bib/article/16/2/325/246489
- [PR] He et al. 2017, SimBoost, J. Cheminform. 9:24 — pair-feature gradient
  boosting; continuous interaction scores; cold-start variants. Suitability:
  PARTIAL (feature design).
  https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0209-z
- [PR] Nelen et al. 2025, J. Cheminform. 17:8 — matched pairs robust to
  inter-assay variability; absolute values rarely comparable across assays.
  Direct external validation of within-panel contrasts. Suitability: YES
  (estimand choice). https://pubmed.ncbi.nlm.nih.gov/39833966/
- [PR] Kalliokoski et al. 2013, PLoS ONE 8:e61007 — comparability limits of
  mixed public IC50. Suitability: YES (caveat source).
  https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0061007
- [PR] Cheng & Prusoff 1973, Biochem. Pharmacol. 22:3099-3108 — Ki-IC50
  relation depends on substrate concentration: single-concentration %
  inhibition is mechanism-mixed. Suitability: YES (assay semantics).
  https://doi.org/10.1016/0006-2952(73)90196-2
- [PR] Haibe-Kains et al. 2013, Nature 504:389-393 — inconsistency in large
  pharmacogenomic studies; batch/platform confounds. Suitability: YES
  (batch-leakage analogy).
  https://www.nature.com/articles/nature12831
- [PR] Duong-Ly et al. 2016, Cell Reports 14:772-781 — the panel itself:
  functional HotSpot assay, % remaining activity at fixed concentration;
  mutant kinase opportunities. Suitability: primary data source.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4740242/
- [PR] Anastassiadis et al. 2011, Nat. Biotechnol. 29:1039-1045 — comprehensive
  kinase selectivity, same functional style; census says replication-only here.
  https://www.nature.com/articles/nbt.2017
- [PR] Davis et al. 2011, Nat. Biotechnol. 29:1046-1051 — KINOMEscan Kd
  panel; different endpoint; pooling with % inhibition forbidden.
  https://www.nature.com/articles/nbt.1990

### 7.2 PLM mutation effects; context propagation; control tasks

- [PR] Rives et al. 2021, PNAS 118:e2016239118 — ESM-1b; contextual residue
  representations emerge from scale. Basis of all ESM features here.
  https://www.pnas.org/doi/10.1073/pnas.2016239118
- [PR] Meier et al. 2021, NeurIPS 34 — ESM-1v zero-shot mutation scoring
  (masked marginals). NOTE: scoring needs the mutation coordinate; predicts
  function-fitness, NOT ligand-conditioned response. Suitability: PARTIAL
  (proposition A flank only).
  https://proceedings.neurips.cc/paper/2021/hash/f51338d736f95dd42427296047067604-Abstract.html
- [PR] Notin et al. 2023, NeurIPS Datasets & Benchmarks — ProteinGym;
  mutation-effect benchmarks and their substitution-level scope. Suitability:
  PARTIAL (external benchmark; not ligand-conditioned).
  https://proceedings.neurips.cc/paper_files/paper/2023/hash/cac723e5ff29f65e3fcbb0739ae91bee-Abstract-Datasets_and_Benchmarks.html
- [PR] Riesselman et al. 2018, Nat. Methods 15:816-822 — DeepSequence VAE;
  mutation effects from MSAs. Needs MSA/UniRef snapshot (absent locally).
  Suitability: NO (blocked input).
  https://www.nature.com/articles/s41592-018-0138-4
- [PR] Frazer et al. 2021, Nature 599:91-95 — EVE; disease-variant effect
  from MSA. Suitability: NO (MSA + different endpoint).
  https://www.nature.com/articles/s41586-021-04043-8
- [PR] Notin et al. 2022, ICML — Tranception; autoregressive + retrieval
  mutation scoring. Suitability: PARTIAL (A flank).
  https://proceedings.mlr.press/v162/notin22a.html
- [PR] Hewitt & Liang 2019, EMNLP — control tasks and selectivity: a probe
  must beat a control task, not just decode. This IS the Q1 design's source;
  extend the same logic to predictivity (representation contains info !=
  downstream increment). Suitability: YES (methodology).
  https://aclanthology.org/D19-1275/
- [PR] Belinkov 2022, Comput. Linguist. 48:207-219 — probing survey:
  representational information vs causal use. Suitability: YES (frame for
  "delta exists but no increment").
  https://aclanthology.org/2022.cl-1.7/
- [PR] eSIG-Net 2026, Nat. Methods — residue-level WT/MT ESM at the mutation
  site with joint discrepancy + interaction objectives (PPI endpoint).
  Precedent for the two-head design; endpoint not equated with ours.
  Suitability: PARTIAL (design precedent).
  https://www.nature.com/articles/s41592-026-03086-x

### 7.3 Mutation effects on binding; kinase resistance

- [PR] Li et al. 2021, Commun. Biol. 4:1262 — PremPLI: mutation-induced
  protein-ligand affinity change (ΔΔG-style) from complex structures.
  Closest published task; REQUIRES co-complex coordinates we do not have.
  Suitability: NO (structure input), boundary reference.
  https://www.nature.com/articles/s42003-021-02826-3
- [PR] AI for kinase-inhibitor resistance (review), 2025, Wiley JIM4 —
  method census for resistance-mutation prediction; most require structure or
  large labelled cohorts. Suitability: PARTIAL (landscape).
  https://onlinelibrary.wiley.com/doi/full/10.1002/jim4.70021
- [PR] Kanev et al. 2021, NAR 49:D562 — KLIFS 85-residue pocket definition
  (our F4 features). Suitability: YES (feature semantics).
  https://academic.oup.com/nar/article/49/D1/D562/5932301

### 7.4 DTA architectures and their cold-start failure modes

- [PR] Öztürk et al. 2018, Bioinformatics 34:i821 — DeepDTA.
  https://academic.oup.com/bioinformatics/article/34/17/i821/5093245
- [PR] Tsubaki et al. 2019, Bioinformatics 35:309 — GNN CPI.
  https://academic.oup.com/bioinformatics/article/35/2/309/5050024
- [PR] Nguyen et al. 2021, Bioinformatics 37:1140 — GraphDTA.
  https://academic.oup.com/bioinformatics/article/37/8/1140/5942970
- [PR] Huang et al. 2021, Bioinformatics 37:i151 — MolTrans.
  https://academic.oup.com/bioinformatics/article/37/Supplement_1/i151/6319694
- [PR] CS-DTA — entity-disjoint and similarity-controlled cold-start
  evaluation; documents overestimation under random splits.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC13161074/
- [PR] Chen et al. 2019, J. Chem. Inf. Model. 59 — hidden bias in DUD-E:
  analogue proof that benchmark bias masquerades as model skill (VS domain).
  https://pubs.acs.org/doi/10.1021/acs.jcim.9b00554
- Internal [PR-equivalent, frozen]: BOUNDARY_20260817_NIGHT.md — the
  BindingDB-Ki double-cold level wall across all tested representation
  families; any bridge claim must confront it.

### 7.5 Shortcut learning; annotation artifacts; counterfactual evaluation

- [PR] Geirhos et al. 2020, Nat. Mach. Intell. 2:665-673 — shortcut learning
  taxonomy; exactly the family_shuffle/annotation-shortcut phenomenon.
  https://www.nature.com/articles/s42256-020-00257-z
- [PR] Geirhos et al. 2019, ICLR — texture bias (shortcut) in CNNs.
  https://openreview.net/forum?id=Bygh9j09KX
- [PR] Gururangan et al. 2018, NAACL — annotation artifacts in NLI;
  hypothesis-only baselines == our ligand-only floor logic.
  https://aclanthology.org/N18-2017/
- [PR] Poliak et al. 2018, EMNLP — hypothesis-only NLI baselines.
  https://aclanthology.org/D18-1024/
- [PR] Kaushik et al. 2020, ACL — counterfactually-augmented data: edit the
  input to flip the label-relevant factor, keep the rest fixed; the exact
  isomorph of the mutation-erasure counterfactual (F2).
  https://aclanthology.org/2020.acl-main.711/

### 7.6 Conditional effect estimation; invariance; nuisance removal

- [PR] Shalit, Johansson, Sontag 2017, ICML — counterfactual regression with
  balanced representations (CFR); nuisance balance before effect estimation.
  Suitability: YES (design logic). http://proceedings.mlr.press/v70/shalit17a.html
- [PR] Künzel et al. 2019, PNAS 116:4156-4165 — meta-learners (T/S/X) for
  heterogeneous effects; the two-head decomposition is the T/X-learner
  isomorph. Suitability: YES. https://www.pnas.org/doi/10.1073/pnas.1804597115
- [PR] Nie & Wager 2021, Biometrika 108:299-319 — R-learner: residualize
  outcome AND treatment on nuisances, then estimate the effect; isomorph of
  "parent-profile first, residual interaction second" (Q5).
  Suitability: YES. https://academic.oup.com/biomet/article/108/2/299/5911092
- [PR] Chernozhukov et al. 2018, Econometrics J. 21:C1-C68 — double ML;
  cross-fitting discipline (our F9/T2 construction).
  Suitability: YES. https://academic.oup.com/ectj/article/21/1/C1/5056401
- [PP] Arjovsky et al. 2019, IRM — invariance across environments; powerful
  but needs multiple environments and known failure modes at our n.
  Suitability: PARTIAL (conceptual only at this scale).
  https://arxiv.org/abs/1907.02893
- [PR] Gulrajani & Lopez-Paz 2021, ICLR — domain-generalization methods often
  fail to beat ERM under fair evaluation; warns against importing DG
  machinery at 49 pairs. Suitability: YES (negative calibration).
  https://openreview.net/forum?id=lQdXeXDoWtI
- [PR] Veitch et al. 2021, ICML — counterfactual invariance to spurious
  correlations; formalizes when an edited-input null is the right null
  (supports F2). Suitability: YES (theory).
  https://arxiv.org/abs/2102.09916
- [PR] Crawford et al. 2017, PLoS Genet. 13:e1006869 — marginal epistasis
  test: detect interaction WITHOUT joint interaction features; an alternative
  detection-first strategy if estimation is underpowered.
  Suitability: PARTIAL. https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006869

### 7.7 Meta-learning and few-shot bioactivity

- [PR] Finn et al. 2017, ICML — MAML. http://proceedings.mlr.press/v70/finn17a.html
- [PR] Raghu et al. 2020, ICLR — ANIL: MAML's gain is mostly feature reuse;
  warns that adaptation-head-only designs may not add interaction knowledge.
  Suitability: YES (interpretation discipline).
  https://openreview.net/forum?id=rkgMkCEtPB
- [PR] Garnelo et al. 2018, ICML — Conditional Neural Processes; support-
  conditioned prediction. Suitability: YES (few-shot head design, P-line).
  http://proceedings.mlr.press/v80/garnelo18a.html
- [PR] Garnelo et al. 2018 — Neural Processes (latent).
  https://arxiv.org/abs/1807.01613
- [PR] Feng et al. 2024, Nat. Mach. Intell. — ActFound: pairwise meta-learning
  bioactivity foundation model trained on within-assay pairwise differences
  (ChEMBL); the closest deployed analogue of our within-panel contrast.
  Requires large external assay corpus for pretraining — not available
  locally; reuse report: https://www.nature.com/articles/s42256-026-01187-y ;
  code: https://github.com/BFeng14/ActFound . Suitability: PARTIAL (design
  precedent; data not present).

### 7.8 Conditional energy models and structured prediction

- [PR] Belanger & McCallum 2016, ICML — Structured Prediction Energy
  Networks: energy-based structured outputs; conceptual basis for the
  integrable potential s(P,L). Suitability: YES (conceptual).
  http://proceedings.mlr.press/v48/belanger16.html
- [PR] LeCun et al. 2006 — A Tutorial on Energy-Based Learning.
  http://yann.lecun.com/exdb/publis/pdf/lecun-06.pdf

Boundary summary: everything that needs co-complex coordinates (PremPLI
class), MSAs (DeepSequence/EVE class), or large external assay corpora
(ActFound pretraining class) is OUT of scope for the next stage. Everything
imported must survive the 49-pair / 18-parent scale and the parent-cluster
bootstrap.

---

## 8. Three candidate directions (prioritized)

### Direction 1 (PRIORITY 1, chosen): CIIP-S1 — signal decomposition and
### estimand ladder on Duong-Ly (details in Section 9)

- Hypothesis: the Duong-Ly centered response decomposes into a parent-level
  component (predictable without mutation coordinates; ceiling ~Spearman
  0.58 by E3) and a mutation-specific residual (idiosyncratic; E5); CIIP-1A
  failed because its attribution contrast could only see the latter.
- Protein info: F1-F4 (+F2 erasure counterfactual); ligand: ECFP4 (existing).
- Training: minimal — Form-2 ridge-free linear/MLP probes are EXCLUDED by
  contract (no closed form); small gradient-trained probes and the frozen
  rank-8 potential only; CPU.
- Splits: original 32/8/9 pair split for comparability; PLUS parent-disjoint
  re-split as a SECONDARY lane (report-only against the primary).
- Controls: ligand-label permutation; same-parent wrong-mutation;
  family-preserving shuffle; random window (reinterpreted);
  protein-invariant shift; F2 erasure null.
- Mutation annotation needed: for F1/F2 arms only (never for deployment
  claims).
- Transferable to BindingDB cold-target few-shot? NOT from this stage; it
  only authorizes the Direction-2 preregistration.
- Extra structure/external data: none.
- Risks: C_sharp underpowered (pre-registered UNRESOLVED outcome); erasure
  arm needs one new ESM inference pass (feasibility already demonstrated by
  the propagation audit).
- Success: B-arm beats floors (lo2.5 > 0) AND/OR A-arm separation AND/OR
  C_sharp separation; each with parent-cluster bootstrap + LOPO stability.
- Failure: every ladder cell at floor -> close the Duong-Ly interaction
  route entirely (both parent-level and mutation-specific).

### Direction 2 (PRIORITY 2, gated by Direction 1 B/C_total pass):
### mutation-free family-level conditioner probe

- Hypothesis: a deployable, mutation-coordinate-free protein representation
  (full-sequence ESM / KLIFS pocket) captures the parent-level ligand-response
  component well enough to serve as a transferable family conditioner.
- Protein info: F3/F4 only; ligand: ECFP4.
- Training: small potential/probe on T1; parent-profile cross-fitted
  nuisance reported separately.
- Splits: parent-disjoint (leave-one-family-out) primary — this is the first
  lane that speaks to cold-target transfer.
- Controls: F5/F7/F8 + ligand permutation; the F9 ceiling defines the target
  performance band.
- Mutation annotation: NO (deployment), YES (only in audit arms).
- Transferable to BindingDB? POTENTIALLY — a family-level conditioner is the
  only object that could legally enter a future BindingDB bridge
  preregistration; that bridge is a separate stage with its own assay-
  semantics qualification (Q11) and remains NOT AUTHORIZED now.
- Extra data: none.
- Risks: family-level conditioning may be exactly the "family key" shortcut
  Stage S warned about; parent-disjoint power at 18 parents is limited.
- Success: deployable arm within a pre-registered fraction (e.g., >= 70%) of
  the F9 ceiling on held-out families, beating all controls.
- Failure: closes the mutation-free interaction route; the parent-level
  component is then label structure without a deployable reader.

### Direction 3 (PRIORITY 3, gated by Direction 1 C_sharp pass):
### erasure-counterfactual mutation-specific channel

- Hypothesis: a site-specific delta (correct site vs erased context) predicts
  the parent-residualized response T2 — the only claim that would preserve a
  mutation-specific interaction on this panel.
- Protein info: F1 vs F2 paired; ligand: ECFP4.
- Training: Form-2 two-head (T0m + T2) probe; potential-form only if the
  probe passes.
- Splits: pair-level primary; parent-disjoint secondary.
- Controls: full Section-9 matrix plus same-parent wrong-mutation as the
  decisive negative.
- Mutation annotation: YES (fundamental to the claim; deployment would need
  variant calling — acceptable for a mutation-panel product, not for general
  DTA).
- Transferable to BindingDB? NO — mutation panels do not exist there.
- Extra data: none for the probe; if it passes, a census for external
  mutation panels (licensing-checked) becomes a separate decision.
- Risks: highest; smallest signal component; n=49 pairs.
- Success: pre-registered C_sharp separation with bootstrap + LOPO.
- Failure: closes the mutation-specific ligand-conditioned route on Duong-Ly
  (final, with the estimand explicitly documented as exhausted).

---

## 9. The chosen next step: CIIP-S1 stage specification (preregistration
## draft — freeze and SHA-256 before ANY computation)

### 9.0 Authority and prohibitions

- Writes only to a new directory tools/research/stageCIIP_signal_decomp_<date>/.
- Never modifies: CIIP-1A artifacts, context-propagation artifacts, model/,
  scripts/, dataset/, meta_test (never read at all).
- Duong-Ly test-pair labels: final evaluation of frozen arms only.
- No ridge/closed-form/pseudoinverse; no test-time gradients; end-to-end
  gradient training only; keyed SHA-256 rng streams; no Python hash().
- B-branches labelled: "mutation-free but not necessarily biologically causal".

### 9.1 Stage S0 — read-only audits (before any fitting)

1. Data coverage: pair/parent/position tables; ESM-length exclusion map (P9).
2. Parent overlap: split x parent incidence; WT-row sharing graph.
3. Mutation coordinate: re-verify all 65 mappings against Q0B audit; alias
   ledger check.
4. Ligand overlap: per-pair common-ligand matrix; per-ligand row coverage.
5. Assay semantics: endpoint namespace check (% inhibition only);
   out-of-range cell census (E6 re-derived); concentration metadata from
   source supplement (read-only).
6. Censoring: confirm absence of censoring annotations (blocks interval-
   censored formulations — recorded as a data limitation).
7. Re-derive Section-4 diagnostics under the preregistration (train+val
   only) with keyed rng; discrepancies vs E1-E6 reported.
8. Power computation (Q9.3) and frozen MDE table per contrast.
9. ESM erasure-inference pass (F2): rebuild erased sequences, assert
   WT-erased == variant-erased strings, inference on CPU, cache under the new
   stage directory with SHA-256; assert max erased-input delta <= 1e-5.
10. Leakage audit output: any channel found must be neutralized or recorded
    as a claim restriction BEFORE S1 starts.

S0 exit gate: audits 1-8 complete, no unresolved leakage channel, erasure
cache valid. Only then S1.

### 9.2 Stage S1 — estimand ladder and representation arms

- Targets T0, T0m, T1, T2, T3 (Section 5); T2 uses cross-fitted F9 (train
  parents only for val/test residuals).
- Arms F1-F9 (+F10 only if licensing verified; expected absent).
- Form-1 potential: frozen architecture (rank 8, hidden 64, AdamW 1e-3,
  wd 1e-4, 200 epochs, batch 512, grad clip 10) for F1/F6/F5/F7/F8 on T1/T3.
- Form-2 probes: same optimizer/budget; small MLP (hidden 64) on
  (protein-features, ECFP4) for F2/F3/F4/F9 on T0/T0m/T1/T2/T3.
- Controls (mandatory): ligand-label permutation (keyed, within pair);
  same-parent wrong-mutation (variant window from another covered pair of
  the same parent — reuse of frozen construction); family-preserving shuffle;
  random window; protein-invariant shift. All controls run through the
  identical metric pipeline.
- CPU smoke (<= 30 min, 2 pairs, 5 epochs) must pass structure tests; then
  single seed; multi-seed (3) ONLY if single-seed structure + all negative
  controls behave as frozen (floors at floor, erasure null at zero).
- Metrics per pair per arm: nonconstant flag/rate, var_true, var_pred,
  centered R2, OLS slope, scale_ratio = sqrt(var_pred/var_true), MSE,
  dead-zone (10 units) sign accuracy, Spearman (undefined if constant —
  never 0), per-parent aggregates; paired effects with parent-cluster
  bootstrap (2000 draws, keyed) + leave-one-parent-out sign stability.
- Frozen adjudication rules per proposition (lo2.5 > 0 on the primary
  paired contrast + LOPO sign stability; point thresholds frozen at freeze
  time after the S0 power table, never after results).

### 9.3 Adjudication outcomes and consequences

| outcome | definition (primary contrasts) | authorizes | closes |
|---|---|---|---|
| S1-PASS-B | B-arm beats F7/F8 floors on T1 | Direction-2 preregistration drafting | nothing |
| S1-PASS-A | F1 beats F2-erased on T0/T0m | proposition-A reporting line (mutation overall effect; still not ligand-conditioned) | nothing |
| S1-PASS-Csharp | F1 beats F2-erased on T2 | Direction-3 scoping review | nothing |
| S1-NULL-ALL | all cells at floor | — | Duong-Ly interaction route (both levels); final boundary document update |
| S1-UNRESOLVED-Csharp | C_sharp CI crosses 0 while B passes | Direction-2 only | Direction-3 recorded as underpowered-closed unless new data arrives |

In EVERY outcome: production model/scripts unchanged; CIIP-1B, BindingDB
bridge, and any deployable-representation claim remain NOT AUTHORIZED until
the relevant Direction-2/3 stage passes its own preregistered gates.

### 9.4 What CIIP-S1 licenses on PASS and forbids on FAIL

- PASS-B licenses ONLY: drafting a Direction-2 preregistration (mutation-free
  family conditioner, parent-disjoint). It does not license touching
  BindingDB, production code, or calling the parent-level component
  "mutation-specific".
- PASS-A licenses reporting a mutation overall-effect result; forbids
  ligand-conditioned framing.
- PASS-Csharp licenses Direction-3 scoping; forbids any deployment claim
  (oracle-bound).
- NULL-ALL forbids further Duong-Ly interaction work and routes the
  programme back to the performance track with the boundary updated.

---

## 9.5 Relationship to the concurrent CIIP-2 (OLR-Potential) work

After this plan's diagnostics were completed, a concurrent work product was
found in the tree (created 2026-08-20 09:16-09:27, no results yet):
report/research_ideas/ciip/CIIP2_RESEARCH_REPORT_20260820.md and
tools/research/stageCIIP2_olr_potential_20260820/PREREGISTRATION.md
(+ ADDENDUM_ADD1_SPB.md). Independent review notes:

1. CONVERGENCE (independent read-only diagnostics agree): their sibling-LOPO
   ceiling R2 0.293 / correlation 0.386 vs my E2/E3 (0.442 / R2 0.326);
   their panel-shared ligand baseline R2 0.131 on test vs my E4 (0.060
   leave-pair-out on train+val; same quantity family); their variance
   components (parent-shared ~60% of interaction signal) vs my E3/E5; their
   out-of-range census 23.0% == my E6. Two independent derivations reaching
   the same decomposition materially strengthens R1-R3 of Section 3.
2. Their audit adds two facts this plan did not compute: the CIIP-1A
   random-window arm is fully explained by the train-mean ligand profile
   (R2 0.1291 vs 0.1313), and the WT panel is ceiling-loaded (median WT
   mean 90.9%; sensitivity concentrated in ~84 mid-zone ligands) — a
   one-sided compression that reinforces Q2/Q3 (Section 6) and motivates
   their assay-gain weighting; and the Manning-family-mean prior does NOT
   transfer to held-out parents (R2 -0.021), which directly raises the bar
   for Direction 2 of Section 8: parent-level interaction structure exists
   but does not transfer by family CATEGORY; only sequence-feature transfer
   remains open (their SPB lane).
3. DESIGN GAPS in the concurrent preregistration that this plan covers and
   recommends adding before its training starts:
   (a) no proposition-A lane (T0/T0m) — the mutation overall effect remains
       unmeasured everywhere;
   (b) its residual target removes the LIGAND-side pattern only; the
       dominant PARENT-side component has no within-family ceiling arm on
       the pair-level lane (its C-famprior is family-category-level, which
       its own audit shows fails); add a sibling-profile-prior arm on S1
       lanes (their own 2.4d shows this ceiling was legally reachable yet
       no CIIP-1A arm approached it — the same check must gate A5);
   (c) C-erased is an evaluation-time control on a trained model; the
       cleaner B-measurement (FIT on erased inputs, Form-2) is absent;
   (d) scale: 5 seeds x router + distillation + 13 arms is far from
       minimal; the Section-9.1 S0 audits + estimand ladder are a cheaper
       prequel that can invalidate the premise before the mechanism spend.
4. Compatibility: nothing in CIIP-2 conflicts with CIIP-S1; S1 is strictly
   a prequel. If governance elects to run CIIP-2 first, this plan's S0
   audits (9.1) and arms (b)/(c) above should be prepended by addendum,
   and CIIP-2's A-ladder results must be reported against the
   sibling-profile ceiling, not only against floors.

## 10. Production constraints (restated as binding rules)

1. No oracle mutation-coordinate ESM into production model/.
2. No CIIP potential into BindingDB; no CIIP-1B; no BindingDB Bridge start.
3. No biological-mechanism claim from correct-vs-random-window differences.
4. No success claim from single seed / single parent / few pairs.
5. No bigger backbones/training budget to mask non-identifiability — the
   next stage SHRINKS the question rather than growing the model.
6. % inhibition never relabelled Ki/Kd/pK.
7. Context-propagation magnitude never cited as predictive value.
8. Failure never reported as "biological absence of protein-conditioned
   signal".

## 11. Open assumptions register (must be checked inside S0)

- A1. E3/E5 diagnostics replicate under preregistered keyed rng on train+val.
- A2. All 9 test pairs retain same-parent train support for F9 (verified
  exploratorily: parents ABL1/KIT/EGFR/FGFR4/RET/TEK all present in train).
- A3. Erasure inference reproduces delta=0 on all 49 pairs under the new
  cache (P7 precedent).
- A4. The parent-profile component is not an assay artifact: it survives
  WT-residualization (E2) but its biological meaning (pocket-level ligand
  pharmacology vs platform correlation) remains open — hence the B-branch
  flag "mutation-free but not necessarily biologically causal".
- A5. Direction-2's eventual BindingDB relevance is unproven and gated by
  the frozen BOUNDARY document.

---

## 12. Core Task 1 completion contract (what this plan does and does not finish)

Core Task 1 — "prove the model genuinely exploits protein-ligand conditional
interaction, on deployable information, in service of cold-target DTA" — is a
FOUR-GATE chain. This plan (CIIP-S1) executes gate CT1-a only. Completion of
Core Task 1 requires all four gates, or a rigorous falsification closure at
any one of them.

| gate | claim to be adjudicated | adjudicating stage | status |
|---|---|---|---|
| CT1-a signal attribution | an attributable protein-conditioned ligand-specific response exists on the surrogate panel, decomposed into parent-level vs mutation-specific components, surviving the matched controls | CIIP-S1 (this plan; Duong-Ly functional % inhibition) | DEFINED here; NOT yet executed |
| CT1-b deployable reader | a mutation-coordinate-free representation (sequence/legal priors) reads the signal under parent-disjoint splits | Direction 2 (gated by S1-PASS-B); note family-category prior already fails (R2 = -0.021) | NOT AUTHORIZED |
| CT1-c model-level attribution | the FINAL model's predictions are caused by the interaction term at inference (erasure/permutation/ablation counterfactuals on the native task), not by a family key or assay batch | control matrix of this plan re-applied to the production candidate | NOT AUTHORIZED |
| CT1-d native-task transfer | the interaction exists and is exploited on BindingDB Ki/Kd/pK cold-target (assay-semantics qualification + governed double-cold + beats Tanimoto/ligand-only incumbents) | separate bridge qualification stage (Q11) | NOT AUTHORIZED |

Explicit non-completions of this plan, stated to prevent overclaiming:

1. CIIP-S1 PASS-B does NOT complete Core Task 1: it is within-panel,
   functional-endpoint, and carries the "mutation-free but not necessarily
   biologically causal" flag; the parent-level component may be
   platform-correlated structure rather than transferable binding
   pharmacology.
2. CIIP-S1 PASS-Csharp is oracle-bound and therefore CANNOT satisfy the
   deployable-information constraint of the final task.
3. CIIP-S1 NULL-ALL closes the Duong-Ly route; for Core Task 1 that is a
   rigorous falsification step, not a failure of process — but Core Task 1
   then remains unproven under current legal data and must be recorded as
   such, with the boundary document updated.
4. No outcome of this plan authorizes CT1-b/c/d work, production changes,
   CIIP-1B, or the BindingDB bridge; each requires its own preregistered
   stage descending from the S1 verdict table (Section 9.3).

Positive completion of Core Task 1 is the conjunction CT1-a ∧ CT1-b ∧
CT1-c ∧ CT1-d. Negative completion (a legitimate terminal state per the
programme goal) is a pre-registered falsification at any gate, propagated to
the boundary document with the exact estimand, data, and controls recorded.
