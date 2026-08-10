# MetaSieve experimental evidence ledger

Updated: 2026-08-10.

This is the canonical human-readable summary of completed experiments. It
does not replace immutable manifests, raw JSON results, or the full historical
ledger in `history.md`.

## Evidence precedence

When records disagree, use this order:

1. immutable input/output manifests and raw machine-readable Gate artifacts;
2. `project_state.json` for current authorization and freeze state;
3. this ledger for the current interpretation of completed experiments;
4. `history.md` for the complete chronological failure and decision record;
5. prose supplied from outside the repository, which is not evidence until its
   commits, manifests, predictions and label-access audit are recovered.

Historical test counts and status statements describe the repository at the
time of each experiment. The current consolidated regression count is 207.

## Current conclusion

```text
MATHEMATICAL_OPERATOR_IMPLEMENTED_AND_CONTRACT_TESTED
GEOMETRY_IDENTIFIED
PAIR_COMPATIBILITY_IDENTIFIED
FROZEN_ESM2_B5_DEVELOPMENT_GATE_PASS_6_OF_6
EXACT_RESIDUE_LOCALISATION_IDENTIFIED_IN_DEVELOPMENT
TEACHER_LIGAND_CONDITIONALITY_IDENTIFIED_IN_LABELS
B5_RESIDUE_MARGINAL_IS_GENERIC_POCKET
B5_LIGAND_DEPENDENCE_CONFINED_TO_THE_COUPLING_TERM
TEACHER_EDGE_COUPLING_NOT_IDENTIFIED
EXACT_RESIDUE_ATOM_COUPLING_NOT_IDENTIFIED
LABEL_SEMANTICS_NOT_AMBIGUOUS
X1A_R_DEPENDENCE_PRECONDITION_FAILED
X1B_NOT_RUN_PRECONDITION_FAILED
CYCLE_QUOTIENT_ALGEBRAICALLY_AVAILABLE_BUT_DEPENDENCY_NOT_REPAIRED
CQ_R0_BINDINGDB_SOURCE_CENSUS_REGISTERED_NOT_EXECUTED
AFFINITY_ENERGETICS_NOT_IDENTIFIED
K_SHOT_SECTION_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

Phase 2A resolved the label-side structural attribution question. S2R validated
a bounded, gauge-free binary ordinal estimator on fresh synthetic components.
S3R transferred that estimator to real MONN residue-differential labels and
returned `REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED`. The labels are
ligand-conditioned, but the current frozen ESM2 plus mean-pooled 41-D ligand
basis did not recover the registered direction beyond the controls. A crossed
affinity existence test is now being repaired at its dependence precondition;
no trainable affinity model is authorized.

### X1A authorization correction — direct-DD dependence required (2026-08-10)

The amended X1A Gate historically reported `X1_ICC_PRECONDITION_PASSED`, but
its target/ligand fixed effects absorb cluster-exclusive targets and singleton
ligands, and its signed-residual ICC does not estimate dependence of X1B's
`q=DD^2-v_noise`. The final artifact also used 2,000 rather than the registered
10,000 bootstrap draws. Its historical numbers remain evidence of that amended
development run, but its X1B authorization is withdrawn.

The exact X0-B packing is now materialized label-blind: Ki 11,168 rectangles in
36 clusters and Kd 1,041 in 12; caps select 827 and 605. Current verdict is
`X1A_ICC_PRECONDITION_NOT_ESTABLISHED`. X1A-R subsequently executed the repaired
direct-DD dependence audit and failed for both endpoints. X1B was not run and
X2 is not authorized.

### X1A-R repaired dependence result (2026-08-10)

The repaired audit used exact-assay-aligned rectangles and the planned
statistic `Z=(DD/2)^2-v_D,U`, with no target/ligand nuisance fit. Ki returned
`rho_U=0.120406` and effective `n=200.43`; Kd returned conservative
leave-one-component `rho_U=0.101078` and effective `n=61.05`. Both exceed their
rho limits and fall below the frozen 245 effective-unit requirement. This is a
source-design dependence/information failure, not a training failure.

### Cycle-quotient feasibility (2026-08-10)

The proposed quotient is algebraically available: raw panel cycle dimensions
are 29,677 Ki and 3,279 Kd. It does not repair independence. Exact-assay cycle
dimension is zero because no assay spans multiple targets, and the largest
dependency component holds 48.9%/46.1% of panel quotient dimension. This
authorizes only the registered BindingDB curated-article metadata census, not
numeric affinity reads or training.

### S7/L2B Phase 2B-S3R — real structural direction not identified (2026-08-10)

S2R first returned `BINARY_ORDINAL_IDENTIFIABILITY_REPAIRED`: the sealed
synthetic seed reached held-out component-macro `AP_bidir = 0.662021`. S3R then
used the same gauge-free direct-W estimator on the real structural labels.

The primary panel comprised 46,818 ligand pairs and 112 closure components.
Candidate AP was `0.035880`, chance `0.025472`, frozen B5 differential
`0.031582`, foreign-pair `0.035735`, context-corrupted `0.032336`, and the
trained permuted-label learner `0.037125`.

The candidate exceeded chance by `+0.010408 [LCB +0.006920]`, below the
registered `+0.05` margin. It failed R2-R5 as well. Module participation,
train-only label access, unit norm, stream equality and exact prediction replay
passed, so this is not a numerical-training verdict. It is scoped to frozen
ESM2 residue states, the global mean of 41-D ligand atom features and this
ordinal estimator. Heldout-B and R6 were not opened; affinity reads were zero.

The parent R6 was superseded before execution because pairwise differences do
not identify absolute amplitude, feature origin or difference-null directions.
No raw residual was added to B5. See `PHASE2B_S3R_GATE.json` and
`PHASE2B_S3R_EVIDENCE_CONSOLIDATION.md`.

## Retained PASS evidence

| Stage | Verdict | Verified result | Scope |
|---|---|---|---|
| P0 | PASS | Canonical DTA, sealing, label-firewall and frozen-operator contracts pass regression. | Software/data contract only. |
| P1A | PASS | 14,906 governed holo complexes, 14,906 receptor sequences and 2,869 chemotypes after 739 protected-homology exclusions. | Open structural corpus; no affinity claim. |
| P1B | `PASS_GEOMETRY_IDENTIFIED` | On 1,477 controlled test complexes, correct contact AUPRC 0.43885 versus wrong-protein 0.05149 and wrong-ligand 0.23895. Correct distance MAE 1.97541 A versus wrong-protein 2.66531 A. | Correct-protein contact/distance geometry only. |
| D0-C | PASS | Official ChEMBL37 SQLite archive verified by SHA-256; 343,562 canonical Ki/Kd rows and 41,619 tasks constructed deterministically. | Immutable source corpus provenance; no model trained. |
| D1 | PASS | 3,817 governed tasks, 697 targets and 253 homology/document closure components; fixed fold sizes 1467/588/588/587/587 with zero closure crossing. | Independence governance; no affinity identification. |
| E0S/E0R | Synthetic diagnostic PASS | Teacher statistic was reconstructable from frozen geometry (maximum error 2.19e-7). Objective mismatch was found; a Moore-Penrose witness gave train RMSE 3.18e-8 and holdout RMSE 0.01269, with correct CI 0.99737. | Synthetic realization and numerical diagnosis only. |
| T-BASIS-R0 | Structural PASS | Fixed 288D radial chemistry basis was held-out recoverable and partner-dependent: reconstruction gain 0.5312 [0.4433, 0.5962], partner gain 0.1561 [0.1070, 0.2007]; all 288 coordinates active. | Structural statistic only; affinity reads zero. |

The retained P1B raw Gate is
`report/mechanism_refactor/p1b_gate_pilot20k_seed17_v4/gate_report.json`.
The retained D0/D1 report is
`report/mechanism_refactor/p1r2b_d0_chembl37_v1/STAGE_REPORT.md`.
Later synthetic and research implementations were removed after consolidation;
their results remain in `history.md` and their source trees are recoverable from
Git history.

## Mechanism-to-affinity experiments

| Stage | Main result | Verdict and interpretation |
|---|---|---|
| P1C | Ligand CI 0.71110; correct mechanism 0.60629; deranged 0.60715. Correct-minus-ligand -0.10481 and correct-minus-deranged -0.00086. | FAIL. The legacy readout had no usable correct-protein affinity increment. |
| P1R0 | Legacy statistic was atom-order dependent; PCA32 retained 48.78% of correct-vs-deranged energy. | Readout contract defect confirmed; PCA compression alone did not explain failure. |
| P1R1 | Invariant 288D MIF: ligand 0.71110, correct 0.71874, deranged 0.69516. Correct-minus-ligand +0.00764 [-0.00179, 0.01685]; correct-minus-deranged +0.02359 [0.00990, 0.03761]. | FAIL below frozen +0.03 Gate. Protein contrast recovered, but not ligand-baseline affinity value. |
| P1R2A | Variance decomposition: ligand 77.58%, protein 20.76%, non-additive interaction 1.67%. Correct-minus-ligand -0.00596; correct-minus-deranged +0.02373. | FAIL. The interaction residual remained partner-sensitive but was not affinity-incremental. |
| P1R2B0 | Source OOF correct-minus-ligand: Ridge -0.00263, spline -0.01980, MLP -0.03026. Metaval Ridge +0.00441; nonlinear arms remained non-positive. | FAIL. Greater nonlinear capacity after global pooling did not rescue affinity semantics. |
| P1R2B1 | Strongest metaval MLP correct-minus-deranged +0.03895 [0.01791, 0.06010], but correct-minus-ligand +0.01915 [-0.00296, 0.04096]. | FAIL. Compatibility/wrong-protein penalty without a stable positive correct-protein increment. |
| E-AFF-P0/H0A/H0C | No population-shared radial affinity direction; task-local headroom was ligand/series SAR; centered radial interaction residual did not recover partner affinity. | Negative evidence against the tested radial affinity mappings, not against all pair-local biology. |

These experiments establish the recurring failure mode:

```text
correct protein changes geometry or compatibility
    !=
correct protein adds transferable affinity-ranking information
```

## Data and estimand feasibility

| Stage | Status | Meaning |
|---|---|---|
| F0R | Historical failure closed | Live ChEMBL API rehydration could not reproduce legacy JSONL bytes. This does not invalidate ChEMBL37; D0 replaced the live API with a release-pinned dump. |
| E-AFF-X0 | `STOP_SOURCE_INTERACTION_UNDERDETERMINED` | The original component-as-unit independence requirement was unattainable. X0-B later registered cell-disjoint DD rectangles plus capped cluster dependence; X1A/X1A-R evaluate only that replacement design and do not retroactively convert X0 to PASS. |
| E-AFF-L0/L0R | NOT RUN scientifically | The positive control failed. The 195-task/78-component audit therefore did not test protein-specific affinity location. |
| XP1/XP2 | Development evidence only | Consumed kinase panels contained interaction signal, but did not meet the `k<=5`, double-held-out and fresh external-admission requirements. |
| XP3/XP4/XP5 | Public-data boundary | Low-noise panels had too few independent protein groups; broader BindingDB panels had estimated assay noise (0.7774) above interaction SD (0.4058), and the tested radial basis did not generalize. |

A blocked or underdetermined estimand is not recorded as a failed biological
model. New data may reopen it only under a separately frozen acquisition,
independence and power contract.

## Structural self-supervision evidence

### S7/L2B Phase 1 frozen-ESM2 B5 development Gate

Phase 0 repaired atom correspondence, tie-aware AP, sealed per-pair predictions,
control manifests, and publication/document closure before B5 scoring. The
frozen ESM2 B5 run then achieved macro-AP `0.069601` on held-out A and passed
all six preregistered contrasts: against prevalence, ligand-only, wrong protein,
motif shuffle, wrong ligand, and B4. The smallest 95% lower bound was
`B5-B4 = +0.040392`.

Marginal decomposition localizes the gain to residue prediction. Residue AP
improved by `+0.177221` over B4 (LCB `+0.160064`); atom AP changed by
`-0.009949` with an interval spanning zero. Wrong-ligand residue AP was
`0.245328` versus B5 `0.265114`, so about 92.5% of the residue signal survives
ligand substitution. The admitted conclusion is therefore development-level
exact-residue localisation of a mostly generic pocket. Exact residue-atom
coupling, ligand-conditioned pocket choice, affinity direction, transfer, and
production `z` admission are not established.

MONN supplies no viable time-forward confirmation panel: two additional entries
survive the 2019 cutoff and none survives a 2024 cutoff. Held-out B remains
nested development evidence, not independent replication.

### S7/L2B Phase 2A audit-only attribution (2026-08-10)

Registered by `research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2A.md`, SHA-256
`4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e`, frozen
before any metric, with amendments 01-03 each frozen before the phase it
governs. **Not committed** — chronology rests on hashes embedded in every output
artifact, which is weaker than a git commit timestamp. Nothing was trained; no
affinity value was read.

```text
TERMINAL VERDICT   LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING
NEXT ACTION        preregister one ligand-conditioned residue residual head
```

**The Phase 1 wording is corrected, not withdrawn.** Phase 1's "92.5% survives a
wrong ligand" was a statement about **B5** and was being carried as if it also
described the **corpus**. Its control was an arbitrary foreign ligand. Phase 2A
used a real alternative ligand of the same exact construct against a noise floor
measured from the data — two crystals of the same construct with the same ligand.

| component-macro Jaccard of residue masks | value |
|---|---:|
| replicate: same construct, same ligand, different crystal | 0.6361 |
| alternative scaffold-distinct ligand | 0.4165 |
| **ΔJ, paired over 292 closure components** | **+0.2580 [LCB +0.2344]** |

Registered minimum meaningful effect was 0.05. The variation is
chemistry-associated (Spearman `ρ = +0.3221 [LCB +0.2987]`, ligand-permutation
control median `p = 0.03`, 84.8% of constructs positive), and 80.4% of
scaffold-distinct pairs change the mask meaningfully. Both registered criteria
pass, so `TEACHER_GENERIC_POCKET_ONLY` is refuted on this corpus.

Decomposing the sealed logits on the complete uniformly weighted mask
(orthogonality `1.18e-9` against a `1e-8` tolerance; weighted ALS agreed with
double centering to `1.07e-14`):

| arm | full | residue marg. | atom marg. | additive | coupling |
|---|---:|---:|---:|---:|---:|
| B5 | 0.06975 | 0.04045 | 0.00514 | 0.03983 | 0.01133 |
| B4 | 0.02323 | 0.01619 | 0.00550 | 0.01510 | 0.00637 |
| BX5 wrong ligand | 0.01969 | 0.03595 | 0.00326 | 0.02769 | 0.00346 |
| BP5 wrong protein | 0.00464 | 0.00461 | 0.00492 | 0.00563 | 0.00355 |
| BL ligand-only | 0.00573 | 0.00330 | 0.00573 | 0.00573 | 0.00305 |

A wrong ligand retains 89% of B5's residue marginal but only 31% of its coupling
term: B5's ligand dependence is confined to the pair term, and that term is
small. `B5 coupling − rewiring null = +0.00601 [LCB +0.00461]` and
`B5 coupling − BX5 coupling = +0.00787 [+0.00620]` — both clearly above zero,
both **below** the preregistered 0.01 practical margin, and not rounded up. The
teacher's own edge coupling fails too (median `z = +0.413` against a threshold
of 2.0, 63.4% above their own null), reproducing I-2 under a stricter rewiring
specification with zero degree-preservation violations.

The well-posed label-fitted additive ceiling is **0.3889**; B5 reaches **17.9%**
of it, and its residue marginal reaches 19.8% of the true residue-margin ceiling
of 0.2043. **The bottleneck is the residue marginal, not the coupling.** The
registered logistic Rasch null is reported but flagged `rasch_converged: false`
— the design is completely separated at 0.07% positive density — and is
explicitly **not** used as a ceiling.

Label semantics are audited and **not ambiguous**: water-mediated edges are 8.2%
(threshold 20%) and removing them strengthens the teacher result
(ΔJ 0.258 → 0.278); a local dense-distance comparator built on 1,909 complexes
(median mapped identity 1.000) places 88.1% of PLIP positives within 5.0 Å of a
ligand heavy atom. A second interaction-annotation tool does not exist locally,
so that one comparison remains **UNRESOLVED**. PU learning and a soft teacher
remain unauthorized.

Contract: all seven Phase 0 checks passed over 26 hashed artifacts, including
proof that the B4-family and B5-family offset tables are identical key-for-key —
which is what validates the Phase 1 marginal decomposition. Census:
`DATA_IDENTIFIABLE` with 1,093 multi-scaffold constructs across 779 closure
components and 323,410 within-construct scaffold-distinct pairs.

### S7/L2B Phase 2B — contract repaired, then stopped fail-closed (2026-08-10)

Registered by `PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL_R1.md` (`5e6688f6…`),
committed `b9753db` **before** any implementation existed. It supersedes
`ae6d1a01…`, which was **never executed** and is kept byte-identical under
`SUPERSEDED_BEFORE_EXECUTION_DESIGN_DEFECT`; its eleven defects are itemised in
`PHASE2B_DESIGN_AUDIT.md`. Since that document produced no result, nothing is
withdrawn — a design was replaced before use.

```text
TERMINAL VERDICT  PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED
```

Every artifact and numerical precondition passed: exactly 10,568 trainable
parameters with no bias, `g(L)` atom-permutation invariant at `0.0`, ligand-order
swap sign-exact at `0.0`, the protein-only prior cancelling in the same-protein
difference at `2.05e-15`, projection orthogonality `6.19e-15`, zero train/held-out
component overlap, zero held-out ligand-graph overlap. The census matched the
registration: 226,765 training and 46,818 held-out A eligible pairs, foreign-pair
control coverage 1.000, derangement with 0 fixed points.

The stage then stopped at its own synthetic trainability control:
`AP_bidir = 0.3577` against a preregistered `>= 0.50`. The threshold was not
lowered, nothing was tuned against the synthetic holdout, and no second seed was
tried. The real-label run was not executed and gates `R1`–`R6` were not scored,
so **no biological conclusion is permitted from this stage**.

Gauge-invariant diagnostics localise the shortfall: the teacher scores `0.99971`
on its own labels, so the metric is sound; in-sample `0.3654` versus held-out
`0.3577`, so there is no generalisation gap; and the learned delta field
correlates with the teacher field at `0.754` median, so the hypothesis class is
being fitted. What failed is the registered optimization budget — and possibly
the a priori `0.50` threshold, set without a calibration curve. Neither may be
adjusted now, because both were frozen and the synthetic holdout has been seen.

### S7/L2B Phase 2B-S0 — the synthetic control was itself invalid (2026-08-10)

Registered by `PREREG_PHASE2B_S0_SYNTHETIC_FAILURE_LOCALIZATION.md`
(`81675578…`), frozen before any S0 measurement.

```text
TERMINAL VERDICT  SYNTHETIC_CONTROL_LOSS_MISALIGNED
```

S0-A and S0-B passed, establishing that the contract and the implementation are
sound: the 48-epoch nested stream was hashed with 1,680 updates / 166,300
presentations / 30,552 unique pairs counted separately; antisymmetry error 0.0;
identical-ligand differential 0.0; projection orthogonality `5.7e-15`; same-seed
replay bit-identical; and the teacher parameters copied into the production head
reproduced the teacher to `4.06e-07` relative field error, `2.2e-16` AP agreement
and `3.46e-08` on `W = UᵀV`.

S0-C then rejected the **control**. Initialising the student exactly at the
teacher, where `AP_bidir = 1.0000`, the registered loss and optimizer reduce BCE
monotonically while AP collapses — `1.0000 → 0.9732 → 0.5797 → 0.3899 → 0.4963`
at 0/1/10/100/210 updates against BCE `0.6364 → 0.3845`. The preregistered rule
fired on both conditions. A control whose own answer is destroyed by its own
training procedure cannot adjudicate a student, so the failed Phase 2B
precondition carries **no information about the candidate**.

Mechanisms, separated and ranked by an addendum diagnostic that could not change
the verdict: **scale** (BCE `0.6364` at the teacher's own scale versus `0.3408`
at the ray optimum `a* = 20.2`; AdamW moves about 4.5% of the parameter scale per
update, so a 20× scale change rewrites the direction first); **residual
misalignment** (from the ray optimum BCE still falls `0.34075 → 0.33669` while AP
degrades `1.000 → 0.891`); and **knife-edge labels** (rank-8/rank-9 gap median
`0.00222`).

**The Phase 2B report's budget attribution is downgraded to not uniquely
established.** The budget was never tested — the earlier-applicable objective
cause fired first — and it cannot explain a pipeline that discards the answer
when handed it. Sampled-pair coverage, optimizer updates and low-rank
factorization remain untested; implementation and synthetic generalization are
excluded. The original result
`PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED` is preserved
unchanged as historical evidence.

The sole next action is a separate preregistration for a repaired synthetic
control, `PREREG_PHASE2B_S1_REPAIRED_SYNTHETIC_CONTROL.md` (`4850c7d5…`),
written and hashed but **not executed**. Its methodological core is a mandatory
**alignment certificate**: a synthetic control must first prove that its own
answer survives its own training procedure before its verdict on a student means
anything.

Also recorded: `P1_B5_REPORT.md` no longer matches the hash in the Phase 1 triage
(`19c9c205…` → `dbfe8b92…`). The change is wording only, in section 3(b); no
number moved. See `PHASE1_ARTIFACT_SUPERSESSION.json`.

S0-S4 used 1,118 RCSB complexes with 621 protein clusters and 586 ligand
scaffolds, disjoint from the 10,468 P1B-exposed PDB IDs. The deterministic
six-channel 3D teacher passed rotation, translation, atom-permutation and
determinism checks at machine precision.

The aggregate mean-ESM plus ECFP Ridge probe recovered two teacher totals:

| Channel | R2 versus mean | R2 versus random | Correct minus deranged |
|---|---:|---:|---:|
| Directional H-bond | +0.268 [0.166, 0.378] | +0.366 [0.222, 0.505] | +0.037 [-0.015, 0.084] |
| Hydrophobic burial | +0.299 [0.162, 0.454] | +0.307 [0.167, 0.431] | -0.006 [-0.035, 0.024] |

Attribution localized this aggregate signal to the ligand. For H-bond totals,
ligand-only R2 was 0.266, protein-only 0.009 and joint-minus-ligand 0.0015. For
hydrophobic burial, ligand-only was 0.331, protein-only -0.017 and
joint-minus-ligand -0.032.

The valid conclusion is therefore
`AGGREGATE_ESM_ECFP_PROBE_NOT_PROTEIN_SPECIFIC`. S4 did not evaluate P1B's
atom-local, residue-local or atom-by-slot tensors and is not an upper bound on
the whole sequence-plus-2D class. No GPU training was authorized by this
aggregate result.

### ⚠ S5R-2 — VERDICT WITHDRAWN 2026-08-09 (independent audit)

The entry below is retained as the original claim. **Its terminal verdict is
withdrawn** and its "all five conditions pass" line is unit-dependent.

Corrected: `S5_DATA_OR_PREREGISTRATION_CONTRACT_FAIL_CLOSED`;
`SLOT_PROXY_POCKET_CHEMISTRY_SIGNAL_OBSERVED_AS_DEVELOPMENT_EVIDENCE`;
`EXACT_RESIDUE_PAIR_MECHANISM_NOT_IDENTIFIED`; `AFFINITY_DIRECTION_NOT_TESTED`;
`S6_NOT_AUTHORIZED`.

The registered inference unit is the **closure component**. The implementation
used homology groups and scaffolds separately. Recomputed:

| unit | n | cond 1 `A5` vs mean | cond 3 `A5−AD` |
|---|---:|---|---|
| homology | 359 | +0.447 [+0.304, +0.578] | +0.298 [+0.229, +0.371] |
| scaffold | 198 | +0.447 [+0.290, +0.597] | +0.298 [+0.237, +0.353] |
| cells | 546 | +0.447 [+0.310, +0.577] | +0.298 [+0.235, +0.362] |
| **closure components** | **54** | **+0.447 [−0.020, +0.497]** | **+0.298 [−0.076, +0.354]** |

The largest closure component holds **1,291 of 1,449 = 89.1 %** of the panel.
`A5 − AL` corrects to **+0.1383 [+0.0827, +0.2037]** (67 rows had entered arm
`AL` as zero vectors). The panel is **development-exposed** — it is the split the
P1B Gate was decided on. The learned-head comparison is **withdrawn**. Full
record: archived historical path
`report/ssl_s5/S5R2_INDEPENDENT_EVIDENCE_AUDIT.md` inside
`archive/legacy_pre_b5_20260810.zip`; `history.md` F-97.

### S5R-2 pair-local pocket-chemistry enrichment (2026-08-09) — ORIGINAL CLAIM

S5R-2 tested the pair-local contract S4 omitted, on P1B's own held-out test
split of `pilot20k_homology_split_v2` (1,449 usable complexes, 359 homology
groups, 198 scaffolds), with a **zero-parameter** coordinate: the
background-subtracted pocket residue-type composition
`e[t] = q[t] - b[t]`, pocket-weighted by the frozen P1B contact posterior.

| Arm | R2 vs mean, CORE channel APOLAR |
|---|---:|
| oracle (target) | +1.0000 |
| **A5 frozen P1B x composition, 0 parameters** | **+0.4467 [+0.2899, +0.6062]** |
| A1 ligand-only | +0.3396 |
| AL wrong ligand, correct protein | +0.3129 |
| A2 background composition | +0.1879 |
| AD wrong protein | +0.1486 |
| AS chemistry shuffle | +0.0004 |
| AG geometry shuffle | -0.0006 |
| AR capacity-matched random | -0.0023 |

All five registered conditions pass: vs mean +0.4467 [+0.2899]; minus
ligand-only +0.1070 [+0.0240]; minus deranged protein +0.2980 [+0.2289]; minus
geometry shuffle +0.4473 [+0.2907]; minus background +0.2587 [+0.1415].
Terminal verdict `P1B_PAIR_LOCAL_STRUCTURAL_MECHANISM_OBSERVED`.

**Two qualifications are part of the evidence, not footnotes.** Re-derived at
exact residue resolution (n=700, 265 groups), the CORE margin over ligand-only
falls to **+0.078 [-0.003, +0.152]** — not separable from zero — because the
128-slot target itself retains only R2 = 0.516 of the exact-residue quantity
(2.46 residues per occupied slot, 69.1 % chemically mixed). Per-slot
localisation is weak (+0.085 [-0.004, +0.174]). Partner specificity is
unambiguous (+0.356 [+0.279]); pair specificity is small (+0.114 [+0.064]).

Two registered hypotheses resolved: **H-DILUTION REFUTED** (explicit slot
composition and pocket-weighted frozen ESM are indistinguishable, -0.008
[-0.045, +0.041]), and **T-RULE applied** (`K1` pocket coverage, R2 = +0.390,
excluded as tautological with P1B's training target).

A 1,726,406-parameter learned pair map scored **+0.015 against the
zero-parameter readout's +0.394**, so no larger model is justified. The k<=5
support section is **not identified** on this channel: rank saturates at 2.84 of
a possible 5.

Full historical evidence is archived under its original `report/ssl_s5/` paths
inside `archive/legacy_pre_b5_20260810.zip`; see `history.md` F-96.

## Claims not admitted as evidence

The externally supplied S5-S9 SSL report has no matching commits, immutable
manifests, checkpoints, prediction files or result JSON locally or on the
audited remote branches. Its numerical claims remain
`EXTERNAL_CLAIM_NOT_REPRODUCED`. Independent code-contract findings from that
report were separately audited and fixed; those fixes do not validate its
experimental claims.

## Failure taxonomy

The 94 historical failure entries in `history.md` are grouped as follows:

| Range | Theme | Durable lesson |
|---|---|---|
| F-01--F-19 | Early meta-architecture and representation attempts | More attention, deeper encoders and generic adaptation did not establish support or correct-protein causality. |
| F-20--F-60 | Phase-Z data, context, acquisition and recovery | Provenance, endpoint semantics, target independence and immutable releases must be treated as model prerequisites. |
| F-61--F-85 | P1 mechanism and E0 realization | Geometry is identifiable; global summaries and tested affinity readouts are not sufficient. Synthetic objective defects must not be confused with biological failure. |
| F-86--F-94 | Affinity estimand and public-data feasibility | Some questions are underdetermined by available independent panels; a blocked Gate is not permission to relax independence. |
| F-95--F-96 | S5 blocking gate, then S5R-2 execution | A blocked stage is not a negative result. Once unblocked, the pair-local structural mechanism is observable with **zero parameters** — but a channel must be checked against the biological quantity, not only against its own representational proxy, and the proxy is where a PASS can hide a ligand shortcut. |

| F-107 | S4R-A ligand representation audit, then the S4R single-axis graph-aware repair | The mean-pooled 41-D ligand basis was measurably collapsed and replacing it doubled the above-chance gain, but the recovered signal survives a foreign ligand pair intact — a real representation bottleneck can hide a construct-level prior rather than the ligand-conditioned mechanism you were looking for. |

| F-108 | S5D ligand-steering collapse and conditional-estimand diagnostics | The registered mechanism was falsified by its own diagnostic: the estimator *does* steer on the ligand, and an estimand that cancels the pocket confound exactly still finds nothing — a control failing is not by itself evidence that the model ignored the input. |

| F-109 | C0/C1 untouched-corpus exact-correspondence audit | On 1,862 never-scored systems, within-slot AP is 0.9856 against a 0.9540 degree-preserving null — exact atom-residue correspondence is nearly a function of residue contact degree, so the geometry-gated router was closed before training. A registered rule can be wrong about provenance, and a fail-closed check is what catches it. |

| F-110 | X1A crossed-interaction ICC precondition | Historical amended PASS; later audit showed its fitted residual ICC did not identify dependence of the X1B statistic. Authorization withdrawn and superseded by X1A-R. |
| F-112 | X1A-R direct-DD dependence | Exact-assay repair failed for Ki and Kd (`n_eff=200.43/61.05 <245`). X1B was not run; no training authorized. |
| F-114 | BindingDB quotient corpus | Development optimization became executable on 12,457 Ki cells / 320 panels; strict closure has 31 components but largest share 0.8586, so population claims remain closed. |
| F-115 | BindingDB T-BASIS linear witness | First real open-affinity training completed; explained fraction 0.000709 and no correct-over-zero/foreign/deranged confidence bound. Shared 288D linear direction not observed. |

## Active authorization

The BindingDB development pipeline is executable, but its first shared linear
witness failed. No production admission or external confirmation is authorized.
The next stage must be separately preregistered and may change only the
coefficient-sharing assumption to a `d<=5` target subspace; dense profiling
modalities must remain separate from Ki/Kd calibration.

None of the following changes. S2R completed the synthetic repair, S3R completed the authorized real
structural run, S4R completed the single authorized single-axis ligand
representation repair, S5D completed a no-training diagnostic of the estimand,
and C0/C1 completed a no-training correspondence-information audit on an
untouched corpus. S4R failed at R1, S5D failed its own D1 rule and all three of
its Gates, and C1 failed C1a, so heldout-B, R6 and all downstream stages remain
closed. Heldout-B was created by none of them.

The C2 Geometry-Gated Coarse-to-Exact Correspondence Router was **not
preregistered and not trained**. Its authorizing Gate did not pass.

No repair of this estimand is eligible. The registered S4R stopping rule closes
the pose-free ligand representation route, including re-running the stage at a
larger vocabulary or radius; the S5D stopping rule closes the conditional
estimand route and forbids a fourth estimand variant on heldout-A, which has now
been consumed three times. Any further work on this question requires a
separately governed information stage with its own preregistration.

The following remain frozen:

- real ChEMBL/BindingDB affinity training;
- DAVIS, KIBA and recipient labels;
- new PLM, second protein encoder, attention stack, geometry/pose branch,
  typed-interaction branch, affinity head, PU loss, knowledge graph or parallel
  module;
- a larger ligand vocabulary or Morgan radius as a rescue of S4R;
- the C2 correspondence router at the 6.0 A / 128-slot contract, and any
  widening of the correspondence corpus or relaxation of its threshold;
- typed-interaction production integration;
- few-shot adaptation and any `k`-shot claim;
- production biological `z`;
- CSMO, Band, mesh and frozen theory;
- P2-P4.

The B5 development PASS does not authorize a source-affinity Gate. Production
admission still requires closure-component OOF
`correct-ligand >= 0.03` and `correct-deranged >= 0.03`, both with 95% lower
confidence bounds above zero, followed by a sealed transfer Gate.

## Recovery and reproducibility

Terminal-negative implementations and duplicate reports were removed after
their evidence was recorded. They remain recoverable from commits `3281780`,
`12a2765`, and `608decf`. Large releases, embedding banks and caches are not
redistributed; see `DATA_AVAILABILITY.md`.

For the current executable state run:

```powershell
conda run -n drug python -m pytest -q
```
