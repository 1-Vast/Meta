# Current task contract

Updated: 2026-08-20. This file is the single authority for active work. Detailed
chronology belongs in `history.md`; numerical evidence belongs in
`report/EVIDENCE_LEDGER.md` and each stage's frozen artifacts.

## Mission

Build a trainable model that produces a reproducible performance step-change on
cold-target zero-shot and few-shot DTA. The research target remains MSE <= 1.00
pK^2 across k in {0, 1, 2, 3, 5}, accompanied by competitive RMSE, CI and
Spearman rather than an MSE-only calibration shortcut. This is an aspiration,
not permission to weaken the split, tune on held-out labels or relabel a failed
gate.

One core innovation must be a training mechanism with a separately attributable
effect. The second intended innovation is a protein-ligand interaction
representation.

The programme has two parallel tracks. **Main Line P** may continue to test
practical cold-target few-shot performance without waiting for a mechanism
claim. **Main Line M / Core Task 1** decides whether a gain may be attributed to
transferable protein-conditioned interaction. Failure or non-identifiability on
Main Line M does not erase a valid performance result on Main Line P; it blocks
only the protein-interaction claim, the CIIP/BindingDB mechanism bridge, and
production integration of that unvalidated mechanism.

## Current execution state (authoritative, 2026-08-20)

- **Final mission:** a reproducible performance step-change in cold-target
  zero-shot and few-shot DTA with a trainable model and no ranking regression.
- **Core Task 1:** **UNRESOLVED overall.** It asks whether a deployable protein
  representation can produce ligand-dependent effects that transfer to unseen
  proteins/families and survive matched protein, label, assay and context
  controls.
- **Current CIIP-1A mechanism:** **ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED.** On the
  49-pair oracle-covered Duong-Ly subset, correct and random local ESM windows
  were both nonconstant on 9/9 test pairs; correct-minus-random pair-mean R2 was
  -0.1217 with parent bootstrap [-0.4569, +0.0327]. This is a scoped failure of
  the oracle-coordinate, low-capacity potential, not biological falsification.
- **Context propagation:** observed at representation level. Mutation erasure
  made 49/49 WT/variant inputs identical and produced maximum embedding delta
  0.0; original non-site states had mean delta norm 0.05749. Context-only
  ligand-conditioned predictive value remains **NOT EVALUATED**.
- **Authorization:** CIIP-1B, a BindingDB interaction bridge, a deployable
  protein router, and production `model/`/`scripts/` integration are **NOT
  AUTHORIZED**. The completed read-only context audit does not itself authorize
  a successor. Any new mechanism experiment needs a new preregistration and
  explicit authorization.
- **Performance track:** remains scientifically independent of the mechanism
  gate. Performance gains without a passed mechanism attribution matrix must be
  reported only as performance gains.

## Current scientific state

- The leak-free BindingDB-Ki double-cold T2 baseline has MSE
  2.5961/1.7712/1.3245/1.2197/0.9859 at k=0/1/2/3/5.
- At k=0, MSE decomposes as level^2 1.7314 plus centered MSE 0.8648. The target
  level is the dominant error, while current within-target ordering is also weak.
  An oracle target level would give approximately 0.865 MSE, so the excellence
  target is arithmetically possible even though no tested model reaches it.
- Fixed Morgan/Tanimoto residual weighting is the strongest reproducible k>=2
  query-specific comparator. It is ligand-only and is not a protein-conditioned
  meta-learning result.
- `occupancy` retains a small resolved ligand-side ordering signal. Existing
  interaction endpoints underuse it.
- Stage S rejected a global protein-conditioned SAR field. Correct protein gave
  only +0.0065 Pearson over ligand-only, with an unresolved interval, while a
  shuffled protein reproduced the gain. The live protein path behaved as a
  target key rather than transferable biology.
- Stage S also proved that a correct-vs-wrong protein hinge can cheat by making
  the wrong branch arbitrarily bad (hard-wrong MSE 12.62). Wrong-protein damage
  alone is therefore not evidence of protein-conditioned signal.
- Stage T tested a core-blind MMP key and its pooled-protein discriminator
  failed decisively (zero predictor wins 2.4x on MSE; shuffled protein improves
  it). A forensic correction
  (`tools/research/stageT_mmp/CORRECTION_20260817_CORE_KEY.md`) scoped that
  verdict: 40.4% of Stage T training D rows compared targets with disjoint
  cores, so the exact-cancellation claim does not hold for Stage T, and the
  Stage T closure of the whole family is withdrawn.
- Stage U (`tools/research/stageU_mmp_interaction/`) froze a core-inclusive key
  (`tau = shared core + R_a -> R_b + attachment environment + stereochemistry +
  charge`) and stopped at **U0**: the primary same-panel fit bank has 37,945
  observations / 243 targets / 30,463 exact keys, but one target in one
  component contributes **29.63%** of observations, above the frozen 25%
  degree-concentration cap. Per the frozen stop rule, U1/U2 were not run and no
  neural model was trained. This is an admission/balance negative, **not**
  biological absence of interaction signal.
- Stage V (`tools/research/stageV_core_mmp/`) is the corrected successor and
  closes the tested mechanism at V0/V1: the primary internal repeated-key
  surface has **32 rows / 4 components** (<100 evaluability floor), and on fit
  the between-target interaction MS is **0.452 pK² vs 0.858 pK² supervision
  noise**, with bootstrap `theta = -0.406 [-0.689, -0.058]` entirely below
  zero. Internal rich exact keys = **0**. No V2 operator was trained. The
  correct reading is **not estimable**, not "closed": `MS_effect` is an upper
  bound on protein-attributable interaction variance (sd <= 0.67 pK) even
  before any noise subtraction, and the noise reference itself is inflated by
  T0's non-identifiable technical-versus-condition split. This bounds what the
  current BindingDB-Ki protocol can *identify*; it is **not** evidence that
  protein biology fails to modulate SAR, and **not** a universal absence claim.
- No present DTA cell has a legal common-frame protein-ligand pose. Cartesian or
  complex-geometry claims are prohibited on this corpus.
- CIIP-1A has now supplied a scoped model-level negative: oracle mutation-site
  ESM is mutation-sensitive, but the tested potential did not extract
  ligand-conditioned value beyond random-window and matched controls. This does
  not supersede the broader Core Task 1 status of UNRESOLVED.

Authority for these statements:

1. `report/POST_COMPLETION_REVIEW_20260818.md`
2. `report/BOUNDARY_20260817_NIGHT.md`
3. `report/CURRENT_MODEL_EVIDENCE.md`
4. `report/EVIDENCE_LEDGER.md`
5. `tools/research/stageS_sar_field/REPORT.md`
6. `tools/research/stageT_mmp/CORRECTION_20260817_CORE_KEY.md` (scopes Stage T)
7. `tools/research/stageV_core_mmp/REPORT.md` (corrected Phase-1 test)
8. `tools/research/stageV_core_mmp/STAGE_U_GOVERNANCE_AUDIT.md` (Stage U stopped)
9. `tools/research/stageCIIP_potential_bridge/CONTROL_REPORT.md` (current CIIP-1A verdict)
10. `tools/research/stageCIIP_context_propagation_20260820/CONTEXT_PROPAGATION_REPORT.md`

## Historical five-phase mechanism programme

The phases below remain the dependency chain for **mechanism attribution**.
They do not block Main Line P performance baselines. A phase may reject with one
seed; promotion requires the stated evidence. A failed prerequisite blocks the
downstream mechanism claim and mechanism integration. Do not jump from an
unidentifiable interaction representation to a larger meta-adapter. Later dated
status blocks are chronology; they cannot override the current execution state
above.

### Phase 1 - transferable protein-conditioned interaction signal

**Question.** Using only protein sequence-derived features and ligand 2D
features, can the model learn a local interaction representation `Z(P,L)` whose
dependence on the correct protein improves true SAR ordering on unseen protein
components?

**Representation contract.** Preserve ligand atoms or functional groups,
protein residues or regions, and their local cross-interaction tokens until the
readout. Early global pooling is not allowed. A protein embedding changing under
a swap is insufficient; the change must align with affinity differentials.

**Current status: UNRESOLVED overall (2026-08-20).** BindingDB exact-MMP evidence
remains not estimable, while the separate Duong-Ly CIIP-1A oracle-local
potential now has a scoped model-level verdict of **NOT SUPPORTED**. Neither
result establishes biological absence or validates a deployable representation.
This supersedes both the earlier whole-family "CLOSED" wording and the later
statement that no model-level verdict existed.

**Stage T (`tools/research/stageT_mmp/`) — scoped to what it tested.** Its
`exact_key` omitted the shared core, it median-pooled several cores into one
target effect, and its descriptor never read the core. Forensics
(`tools/research/stageV_core_mmp/STAGE0_FORENSICS.json`): **40.4% of Stage T fit
`D` rows and 28.9% of internal rows compared targets with disjoint core sets**,
carrying an uncancelled generic residual of median **0.269 pK** (p95 1.268)
against a `D` truth sd of 0.804. So Stage T's "cancels `mu_tau` exactly" claim
is false for those rows, and its global closure of the family is **withdrawn**.
What stands: the **coarsened-key pooled-protein discriminator is rejected**, and
the defect explains its own inverted controls. Correction document:
`tools/research/stageT_mmp/CORRECTION_20260817_CORE_KEY.md`. T0 stands
unchanged; the "1,112 rich keys" figure must not be reused.

**Stage U — stopped and superseded.** Its preregistration was frozen at 17:12,
four minutes after Stage T's gate metrics were read at 17:08, so it is an
**adaptive correction, not an independent confirmation**. Its chemistry was
correct (core-inclusive key, core-consuming descriptor, interaction-variance
gate, local region operator), but four load-bearing requirements were missing
after metrics had been read: no residue-token permutation control, no
capacity-matched random protein, a `fit_unsampled` bank that retained the same
targets and keys (so it cannot detect target-key memorisation), and no
identical-initialization / identical-minibatch-order rule. Per the governing
rule it is stopped, its preregistration is **not** edited, and no Stage U number
is used as evidence. Audit:
`tools/research/stageV_core_mmp/STAGE_U_GOVERNANCE_AUDIT.md`.

**Stage V (`tools/research/stageV_core_mmp/`) — the corrected test, stopped
before training.** It inherits every Stage U threshold verbatim (frozen before
any core-inclusive number existed) and adds the four missing controls plus two
defect repairs; **no threshold was loosened**. Three frozen gates fail:

- **V0 census**: all five admission thresholds pass (37,945 fit observations,
  243 targets, **1,001** core-inclusive rich keys, 4,589 internal observations,
  25 internal components) but **two concentration caps fail** — one target,
  alone in its component, carries **29.63%** of same-panel fit observations
  against a 25% cap. This independently reproduces Stage U's U0 stop.
- **V0b evaluability** (inherited <100-row rule): the designated primary
  surface, internal `D` rows whose key repeats in fit, has **32 rows over 4
  protein components**, and **internal rich exact keys = 0**. The primary
  evaluation surface for the corrected estimand does not exist here.
- **V1 interaction variance**: within one *complete* transformation, the pooled
  between-target mean square is `MS_effect` **0.4517** against a supervision
  noise reference of 0.8576, giving `theta` **-0.4059 [-0.6889, -0.0577]** —
  resolved below zero, ratio **0.527**; internal agrees descriptively (0.2610,
  ratio 0.304). Even ignoring the noise reference entirely, `MS_effect` is an
  **upper bound** on protein-attributable interaction variance: sd <= 0.67 pK
  before subtracting any noise.

**Consequence.** Per the frozen stop rules, **no neural arm was built or
trained**. The scoped conclusion is about **estimability, not biology**: on this
corpus, with same-panel Ki supervision and full-context single-cut MMP
transformations, the protein x transformation interaction is not identifiable
above supervision noise and no adequately supported evaluation surface exists on
the withheld components. **Insufficient support is not biological absence.** It
may **not** be claimed that protein-conditioned interaction representations are
impossible, nor that Phase 1 is closed. Authority:
`tools/research/stageV_core_mmp/REPORT.md`.

**What would unblock Phase 1.** Not a larger model. Either a corpus in which a
complete transformation recurs across many protein components, or a
preregistered looser-but-principled transformation equivalence class with its
own cancellation analysis — since loosening the key reintroduces exactly the
0.269 pK residual measured above. Stereochemical and charge-changing edits
remain **unmeasured** (1 and 2 internal observations).

**Required controls.** Correct protein, shuffled protein, similarity-matched
wrong protein, protein-blind ligand-only, residue permutation, and a
capacity-matched random protein branch. All comparisons use identical rows.

**Pass gate.** Correct protein must improve signed relative-affinity prediction
over ligand-only and every wrong/shuffled control, with the component-bootstrap
lower bound above zero on the preregistered estimand. Report MMP coverage,
scaffold novelty, protein novelty, ligand similarity and activity cliffs.

**Stop.** If correct protein is statistically indistinguishable from shuffled or
matched-wrong protein, close the proposed representation family and do not enter
Phase 2 or few-shot adaptation.

### Phase 2 - protein-conditioned relative SAR

**Question.** Can the admitted representation predict
`Delta y(P,i,j) = y(P,L_i) - y(P,L_j)` across unseen protein components?

**Model contract.** The relative operator must be antisymmetric under ligand
exchange and exactly zero for an identity pair. Train directly on signed
within-target differences with a robust relative loss. Target balance and assay
or panel provenance are mandatory because raw pair counts are highly skewed.

**Counterfactual rule.** Wrong, shuffled and residue-permuted proteins are
evaluation controls by default. A protein-counterfactual training loss is
forbidden unless a preceding synthetic test proves that its bounded objective
cannot pass by worsening the wrong branch, and the real-data report separately
shows improvement of the correct branch. This rule supersedes the naive margin
loss proposed in the original five-phase note.

**Pass gate.** Relative Pearson/Spearman/CI and sign accuracy must beat the
protein-blind relative model and all protein controls, including in
high-similarity/high-affinity-gap activity cliffs. Improvement only on familiar
ligands or within seen components is insufficient.

### Phase 3 - absolute level and SAR-shape ownership

**Question.** Can zero-shot absolute affinity improve without sacrificing the
relative SAR coordinate established in Phase 2?

**Model contract.** Use an explicit decomposition
`y_hat(P,L) = mu_hat(P, assay_context) + s_hat(P,L)`. The shape branch is
centered and cannot cheaply relearn a target constant. Level and shape use
separate readouts and explicit gradient ownership.

**Required training innovation.** Test a focused trainable mechanism such as
disentangled gradient routing, gradient projection or independently supervised
heads. It must be implemented as one forward/backward optimizer step unless a
separately preregistered inner/outer loop is the experimental variable. Report
per-loss gradient norms, conflicts and parameter-update coverage so the claimed
training contribution is observable.

**Pass gate.** Against the Phase-2 model, target-level squared error decreases,
centered MSE does not increase, and Spearman/CI do not regress. A level gain
obtained by degrading ordering fails. The training mechanism must show a
positive effect in the attribution arm with the interaction trunk held fixed.

### Phase 4 - genuine few-shot target-state identification

**Question.** Given support `{(Z(P,L_i), r_i)}`, can a trainable set-conditioned
operator identify a new target's local SAR state and improve each query beyond
fixed chemical similarity?

**Contract.** k=0 is exactly the admitted zero-shot model. k=1 may make a
conservative level update plus a strongly regularized local correction; it
cannot identify a complete SAR field. At k>=2, use relative support evidence
such as residual or affinity differences. The correction must vary by query and
must be support-permutation invariant. Query labels, closed-form solvers, ridge,
pseudoinverses and deployment-time query gradients are prohibited.

Attention is admissible only as a correspondence mechanism over preserved local
interaction tokens. Generic attention or a pooled task code is not itself an
innovation. The support label must be bound to its molecular context before set
aggregation.

**Pass gate.** Beat level-only adaptation and fixed Morgan/Tanimoto residual
transport at k>=2 in MSE and ranking, with positive correct-vs-label-permuted and
correct-vs-matched-wrong support gaps. At k=1, report level and query-specific
effects separately. A support-independent or label-insensitive adapter fails.

### Phase 5 - causal attribution, confirmation and external replication

Run four matched arms:

1. incumbent trunk + incumbent training;
2. new interaction trunk + incumbent training;
3. incumbent trunk + new training mechanism;
4. new interaction trunk + new training mechanism.

This factorial design attributes representation, training and interaction
effects. Evaluate k=0/1/2/3/5 with MSE, RMSE, CI, Spearman, centered MSE,
target-level error, activity-cliff sign accuracy, scaffold novelty and protein
novelty. Required counterfactuals are correct support, label-permuted support,
matched wrong-target support, wrong/shuffled protein and ligand-only controls.

Promotion requires three fixed seeds, identical nested episode banks,
component-paired confidence intervals, no k=0 ranking regression, and a loadable
checkpoint. Only after BindingDB promotion may Davis and KIBA be trained and
evaluated as separate single-dataset experiments. Do not merge their labels or
report cross-dataset training as the BindingDB result.

## Phase status

| Phase | Status | Admission authority |
|---|---|---|
| 1 | **BOUNDED NEGATIVE (tested mechanisms closed)** — Stage S/P/T/U/V chain; Phase 1 final decision records the exact-MMP route as not estimable and interaction variance as not identifiable above the defensible noise envelope | `tools/research/stageV_core_mmp/PHASE1_FINAL_DECISION.md`, `PHASE1_FINAL_DECISION.json` |
| 2 | BLOCKED pending Phase 1 | this contract |
| 3 | BLOCKED pending Phase 2 | this contract |
| 4 | BLOCKED pending Phases 1-3 | this contract |
| 5 | BLOCKED pending a promoted candidate | this contract |

## Immutable governance

1. BindingDB experiments use the governed double-cold split and physical split
   view. Protein components and ligand identities/scaffolds must be audited.
2. Model selection, early stopping and threshold tuning use only a
   meta-train-derived internal component split. `meta_val` is a development
   evaluation, not a clean held-out claim after selection. No `meta_test` label
   may enter fitting, selection or a reported development metric; future Phase-5
   confirmation must use the repaired physically isolated split-view surface.
3. Stable SHA-256 seeds are mandatory. Python `hash()` is prohibited for any
   reproducibility-critical value.
4. Support and query are disjoint ligands from one recipient target. Query
   labels are used only in the outer training loss or final evaluation.
5. Each experiment trains on one named open dataset. Davis/KIBA are independent
   replications, not silently integrated training data.
6. No ridge, differentiable ridge, analytic solve, Cholesky, pseudoinverse or
   search for a strongest closed-form adapter.
7. Normal trainable inner/outer loops are allowed only when they are the
   preregistered mechanism. No test-time query labels and no hidden multi-stage
   pretraining claim.
8. External representations or labels require a provenance ledger, frozen-data
   baseline and an attribution arm. Literature inspiration is not evidence that
   an implementation works.
9. Common-frame 3D interaction is forbidden without legal coordinates. Ligand
   conformers and protein structures in separate frames may only produce
   independently invariant features.
10. One scientific variable changes per screening stage. A single seed may
    reject; positive promotion requires three seeds.
11. Structural tests precede training: shape/mask contracts, permutation
    invariance/equivariance, k=0 identity, k=1 semantics, label sensitivity,
    counterfactual distinctness, finite gradients and checkpoint round-trip.
12. Do not reinterpret a failed preregistered threshold after seeing results.

## Closed directions and retained comparators

- A2 moment adaptation, seven output-side query adapters, the partial MAML lane,
  adaptive task selection, shared-trunk level heads, global protein-conditioned
  SAR fields and naive wrong-protein margin training are closed on their tested
  implementations.
- These are empirical closures, not universal theorems about MAML, attention,
  Set Transformers, DrugBAN, FS-CAP, OGM or all protein language models.
- Retain the leak-free incumbent, fixed Morgan/Tanimoto transport, ligand-only
  signed SAR field, `occupancy`, and K-REG only as measured comparators or
  training diagnostics. They are not final candidates.

## File and evidence lifecycle

- `tools/research/<stage>/` is the only location for unadmitted research code,
  tests and temporary stage artifacts.
- Every stage needs a frozen `PREREGISTRATION.md`, machine-readable result and
  compact decision report. Progress logs, generated caches, duplicate smokes and
  failed checkpoints are removed after the verdict; Git is the recovery layer.
- A passed implementation moves once into `model/`; its governed workflow moves
  into `scripts/`; maintained contracts move into `tools/tests/`. Remove the
  research copy after promotion so there is one implementation.
- `main.py` orchestrates only admitted `model/` and `scripts/` functionality.
- `history.md` records decisions chronologically. `task.md` contains only the
  active contract. `report/EVIDENCE_LEDGER.md` points to numerical leaves.
- Historical documents still referenced by tests or exact operator definitions
  remain in place but are explicitly non-authoritative. Do not delete an
  evidence dependency merely because its proposal is closed.

## Immediate next action

Phase 1 is decided: `tools/research/stageV_core_mmp/PHASE1_FINAL_DECISION.md`
records the bounded negative. No training is authorized on the Stage U or
Stage V constructions; keep all Stage S/T/U/V artifacts read-only. A future
positive claim would require a new preregistered stage and at least one of:
(a) a corpus where a complete core/context-matched transformation recurs
across many protein components, (b) a governed MSA/coevolution snapshot
(currently absent), or (c) a looser-but-principled equivalence class with its
own cancellation analysis — and (c) alone may only be an identification
screen, never the strict positive gate. No threshold may be moved to revive
this route.

The remaining-lanes audit
(`tools/research/stageV_core_mmp/REMAINING_LANES_AUDIT.json`) records that no
local legal route can change this verdict; the negative branch is complete
pending only external assets or a governance change.

Completion evidence: `report/COMPLETION_STATEMENT_CORE_TASK1_20260817.md`;
goal-tool unavailable, so `tools/research/GOAL_ACTIVE.md` is the authority and
records the complete-on-negative-branch state. No further internal experiments
are to be started until one of the three reopening conditions appears.

Goal status: **COMPLETE** on the bounded-negative branch as recorded in
`tools/research/GOAL_ACTIVE.md` (round 10). No further internal work until an
external asset or governance change opens a new cycle.

## Active direction update (2026-08-17 user redirect)

The exact-MMP branch is closed. The new active goal is Stage W
(`tools/research/stageW_soft_mmp/`): build a transferable local
protein-ligand signal on independent open datasets before any meta-learning.
W0 is complete: **Davis FAIL** (only 7 rich soft families; closed for this
surface), **KIBA PASS** (2,420 rich families / 127 components / 5.5M
cross-component D rows; same-core residual median 0.667 <= 1.0;
between-component MS outside component-permuted null). Immediate next action:
freeze a W1 preregistration for KIBA-only local interaction representation
(multi-pocket residue states x ligand pharmacophore subgraphs, local
cross-attention, independent level/shape heads, full protein counterfactuals),
then execute W1 gates exactly as frozen.

Stage W W1 update: preregistration frozen
(`tools/research/stageW_soft_mmp/W1_PREREGISTRATION.md`, SHA-256
`038f4d97…49082`). W1 split admission **PASSES**: KIBA fit 103 / heldout 24
components; heldout_repeated 186,673 rows / 23 components / 2,040 families.
Next: build the KIBA ESM-2 150M protein token bank, implement the local
interaction operator, run structural tests, then train arms A-F exactly under
the frozen gates.

## Stage W0b audit update

Preregistration frozen (`stageW0b_core1_audit/PREREGISTRATION.md`, SHA-256
`ff23c408…70631f`). Verdict: **W1 GO/NO-GO = NO-GO on current local assets**.
Reasons: W0-P positive control not runnable; Davis/Metz/Klaeger censored at
71.2%/60.4%/93.5%; Davis strict MMP only 7 classes; KIBA descriptive only.
Stage W W1 is PAUSED (no training artifact, no training metric). Immediate next
action: acquire an admissible mutation/ortholog positive-control panel, then
re-census support after censored exclusion/interval-censoring before any W1
training decision.

Stage W0-P panel update: local positive-control candidates found.
`stageW0P_positive_control/W0P_PANEL.json`: 6 near-identical BindingDB
sequence pairs (1–5 mismatches, >=3 shared ligands, 32 total rows), sufficient
by frozen rule. W0-P is now runnable; the gradient-trained positive-control
test with correct/random/BLOSUM-approximate/global/random-protein controls is
the next step. W1 remains NO-GO until W0-P passes and censoring re-census
passes.

Stage W0-P result: **FAIL** (leave-one-pair-out, 3 seeds). Correct positions
sign accuracy 0.240 does not beat corrupted positions; global ESM pooled
difference reaches 0.760 but is recorded as unexplained under n=6 pairs, not a
positive. W1 remains NO-GO. Next authorized action: enlarge the point-mutant
panel or acquire a standard positive-control panel, then re-preregister the
W0-P pass rule from effective sample size.

## Core Task 1 W0/W0-P final decision

`stageW0b_core1_audit/W0B_W0P_FINAL_DECISION.md/json`: **NO-GO on current local
assets.** W0-P FAILED (leave-one-pair-out, 3 seeds; correct positions do not
beat corrupted positions); local panel cannot be enlarged beyond 7 pairs;
Davis/Metz/Klaeger censoring 60-94%. W1 remains PAUSED/NO-GO. Next action is
external W0-P panel acquisition with provenance, then censored re-census and
re-preregistration — not W1 training.

Blocker recorded: `stageW0b_core1_audit/W0P_ACQUISITION_SPEC.md`. Current
status is NO-GO until a standard W0-P panel matching the acquisition spec is
present and provenance-recorded. No W1 training is authorized meanwhile.

W0b censored re-census (round 25): broad all-pairs layers survive censoring
on Davis/Metz/Klaeger (support GO if W0-P passes); strict MMP remains small.
Future W1 screen target = all-pairs/similarity layer; strict MMP as
confirmation layer. W0-P blocker and W1 NO-GO unchanged.

W0b cross-platform residual (round 26): Metz-Klaeger +0.642 [0.482,0.777];
Metz-Davis -0.666; Klaeger-Davis -0.625. Interpreted only as cross-platform
transfer gate: direct residual sharing with Davis is closed; single-platform
signal is not killed. W0-P blocker and W1 NO-GO unchanged.

Core Task 1 terminal verdict: **UNRESOLVED on current local assets**
(`report/CORE_TASK1_UNRESOLVED_TERMINAL_20260817.md/json`). Not SOLVED; not
FALSIFIED-AS-TESTED; not biological absence. Reopening requires the frozen
W0-P acquisition spec, censored re-census and a new preregistered cycle.
No further internal experiments before external asset injection.

## Stage X new independent cycle (2026-08-17 user authorization)

New work in `tools/research/stageX_csc_signal/`; old S-W records read-only.
The entries in this dated Stage X chronology preserve the state at each
historical round. They do not override the authoritative current execution
state above; in particular, Q2c/Q2d are complete and CIIP-1A/context audit are
complete as recorded below.
X0 prereg frozen (SHA-256 `03cdc907…9683`). X0-D downloaded and audited:
Duong-Ly 2016 mutant panel, Anastassiadis 2011 complete matrix, Davis 2011 raw
supplementary tables, PKIS2 supplements.

**Round-1 independent review (2026-08-18): X0 remains active and has not
passed.** Global pooled ESM and amino-acid composition are measured as mutation
insensitive, but the claimed local-representation pass is invalid: WT and
mutant windows were extracted at different coordinates, the mutation-token
inter-protein denominator is zero, four mutation annotations do not match the
downloaded reference residue, and KLIFS is not implemented. The current I1
draft also fails at runtime and does not isolate the planted interaction. I3-I5
are incomplete and the seven tests are initial integrity smokes, not full I6.
Authority: `report/STAGE_X_ROUND1_REVIEW_20260818.md`. Required next action is
instrument correction and X0 qualification, not X0-P/X1 training.
**Round-2 governance + corrected successor (2026-08-18):** original X0 ruled
INVALID INSTRUMENT — the distance-ratio capability gate is a measurement-
definition failure, not repairable in place (verdict
`tools/research/stageX_csc_signal/X0_INVALID_INSTRUMENT_VERDICT.md`; frozen
artifacts untouched). New successor
`tools/research/stageX_csc_signal/stageX0c_measurement_qualification_20260818/`
with its own frozen preregistration (SHA-256
`7de23c81…9922cd`) and ordered gates Q0 -> Q1 -> Q2 -> Q3 -> B1 -> B2 -> C -> D.
Status: Q0-A PASS (ProteinGym 45,623 records, 100% old-residue and 100%
mutated-sequence agreement), Q0-B PASS (BRAF 3-nt alias evidence; PDGFRalpha
= P16234; D842V quarantined; 76/76 Duong-Ly variants typed; KLIFS gatekeeper
= pocket index 45; Davis census), Q1 PASS (probe selectivity:
pair_centered_local_esm +0.189 [0.033,0.363] on pocket membership under
LOO-parent), Q3 census delivered (Saifudeen 2026: 313/349 matched WT gene,
21.3% saturated, CC BY-NC-ND 4.0), I6 23 production contracts + 8 Q0 unit
tests green. Q2 planted harness FAILED its frozen gate (correct-arm Spearman 0.033 /
dead-zone sign accuracy 0.504 / gap -0.018 at tau*=1.0, rank 4; needs
0.30/0.70/0.05). Oracle arm recovers the centred interaction (dz 0.68-0.76),
so the failure is optimization/representation capacity, not information
absence. No B1/B2/C/D authorization; next step is a new preregistration for
the representation fix or a revised gate point.
**Q2c successor (2026-08-18, running):** frozen prereg SHA
1027ccde…26a5c in `stageQ2c_harness_audit_20260818/`. Review amendments:
oracle clause withdrawn (in-artifact oracle dz 0.607-0.674 < 0.70);
Q1 split Q1-A/B/C; Q3 census-only; ANOVA projection removed from diagnostics
(dead-arm Pearson 0.51-0.54 explained by full-yhat projection + inter_scale
drift in the no_interaction_head arm; interior eval cells 0/545; endpoint
distortion Pearson 0.59). Q2c-0/1/1b done, Q3b done: oracle dz 0.664 at
tau*=1.0 vs 0.733 at tau*=2.0 -> 0.70 threshold unreachable at the frozen
gate point even with oracle factors; Q2c-2 precondition unmet -> NOT STARTED;
B1/B2/C/D NOT AUTHORIZED; next = new prereg (Q2d) along the measured power
curve.
**Q2d-1b/1c/1d oracle chain (2026-08-19, running):** Q2d-1 forensics 5/5
(Phase D/E never censored, C 70%, closed-form holdout leaked, half-cold
zero-variance by construction, Phase A main effects + ID bias + unseen-ligand
checkpointing); Q2d-1b STOP (broken ALS oracle + train-only ID centring);
Q2d-1c STOP (train-row feature rank 28 < 32, 8.8% of drawn protein map
unidentifiable; true weights 0.968+ everywhere, ligand map fully
identified); Q2d-1d (prereg baf4bb72...) span-restricted truth: oracle
precheck PASS all seeds (pc 0.700-0.920 / lc 1.0 / dc 0.753-0.893), training
authorized; minibatch target-alignment defect found+fixed+regression-tested
during startup verification; ladder M1 A-E x 3 seeds x 8 arms + M2/M3/NC1/NC2
running; gate adjudication pending. Q2d-2/Q2d-3/B1 still NOT AUTHORIZED.

## Active direction update (2026-08-19, new-window takeover, user re-adjudication)

The single worst-case protocol is split into Main Line P (practical
cold-target few-shot performance; k={5,10,20,40} primary, k={0..5} stress;
mechanism gates do NOT block it; un-gated gains are performance claims only)
and Main Line M (strict mechanism: protein-component cold + scaffold cold +
double-cold + all protein counterfactuals; decides interaction claims only).
Q2d synthetic chain is bounded: adjudicate the running Q2d-1d ladder ->
frozen Q2d-1e (span-init A=V_train@G + L2 1e-3) if 1d fails -> at most ONE
limited diagnostic (explicit A=V_train@G parameterization) -> final PASS/FAIL
on the low-rank bilinear learner. No Q2d-1f/1g chain.

New stage: tools/research/stageP_practical_fewshot/ (prereg SHA
b81283c2..., P1 bake-off prereg SHA 59a90ef2..., addendum AD1 SHA
a675b0ec...). Three evaluation layers frozen: P1 practical few-shot
(protein-cold, same-series support/query allowed, k={5,10,20,40} main /
k={1,2,3} stress, k=0 continuity), P2 novel-target screening (k={0,5},
screening metrics only with frozen threshold labels), M1 fundamental stress
(double-cold + counterfactuals, k={0..5}, mechanism claims only).
P-line artifacts built and tested (p_split.py: quota-balanced whole-cluster
60/20/20 target split, 298/100/101 targets; p_bank.py: per-(target,draw,k)
episode bank, Q=8, 6,376 records, SHA-pinned; 6 tests green). Bake-off arms:
ligand-only, fixed Tanimoto, ordinary fine-tuning, first-order MAML, CNP,
FS-CAP-style, ActFound-style pairwise, AdaMBind-style (frozen addendum
required after FULL-text inspection), current admitted baseline. Adoption of
AdaMBind/MAML decided by measurement only. Multi-dataset phase: BindingDB Ki
/ Davis Kd (raw davis.tab, 379 proteins x 68 ligands, 71.2% sentinel
censored) / KIBA (raw kiba.tab, 229 x 2,068) each with dataset-specific
heads and normalization; Saifudeen = functional selectivity positive control
only, never pK/Ki/Kd DTA. Literature R15 recorded
(report/LITERATURE_R15_20260819.md; AdaMBind Nat Commun 2026 full-text
inspection pending).


## Phase 1 verdict update (2026-08-19)

Q2d-1d GATE FAIL (M1:A double-cold: correct dz 0.5616 / sp 0.1278 /
family-preserving negative 0.6121 beats correct; correct dc dz 0.432-0.618
across ladder A-E vs oracle 0.761-0.992). Ladder crashed at M2 (frozen
PCA_VT NameError); M1:A recovered with the exact frozen code path
(7/8 arms bitwise; family-preserving PYTHONHASHSEED note recorded).
Q2d-1e (span-init + L2, frozen prereg 61bc0cc5...) is RUNNING on GPU with
AD1 truth repairs frozen BEFORE launch (sha 0b405df9...; M1/M2/M3 streams
bit-identical, 9 tests). If 1e fails -> one frozen diagnostic
(A=V_train.G span-param) -> terminal PASS/FAIL on the low-rank bilinear
learner family. Q2d-2/Q2d-3/B1 remain unauthorized.

## Round update (2026-08-19, Q2d-1e running)

- Global re-adjudication executed: full read of task/history/ledger/
  record/evidence/programme + Q2d-1d report + 1e/diag preregs + P1 specs.
  Q2d-1e ladder untouched and healthy (M1 D seed 0 at last check; watcher
  auto-adjudicates on exit). 1e level-A medians (0.585/0.588/0.372) already
  below the 0.70 gate -> 1e FAIL expected; frozen span-param diagnostic
  will follow; no 1f/1g chains.
- Core Task 1 funnel reset and frozen as plan:
  tools/research/stageX_csc_signal/CORE1_FUNNEL_PLAN_20260819.md
  (structure tests -> single-seed 4-arm screen -> stop rules -> 3-seed
  full negatives + bootstrap -> real matched-pair data requirement).
- P-line fairness fixes (no training started): MAML outer-gradient
  contamination reproduced by a red test, fixed in task_fomaml_grad
  (adapted-model grads cleared before query backward) and verified against
  a toy functional reference; CNP re-adjudicated to deterministic
  Deep Sets (AD2 sha f8909ede..., decision A over latent NP), k=0
  context correction exactly 0, permutation invariance / query
  equivariance / query-label isolation tests green. 17 stage-P trainer
  tests green. Real arm-3/4/5 training remains queued behind the 1e
  ladder (GPU budget).

## Round update (2026-08-19): P1 arms 6/7 implemented

- AD3 frozen (sha 5c573132...): arm 6 = FS-CAP-style ligand-only
  Deep-Sets support encoder (no protein features anywhere; k=0 context
  correction exactly 0; total param delta vs arm 3 +96,446 recorded);
  arm 7 = ActFound-style pairwise supervision (antisymmetric +
  identity-zero by construction; eval = support-anchored differences;
  k=0 = frozen p_train label mean). Both: shared sampler rng stream,
  AdamW 3e-4 / 6000 steps, same monitor + checkpoint rule, eval with
  query labels never entering. 9 tests green; stage-P suite rerun.
- P1 bake-off implementation state: arms 1-7 implemented and tested;
  arm 8 (AdaMBind-style) remains BLOCKED on full-text inspection
  (snippet-only access in this environment; bake-off prereg requires
  full text before its addendum). No real training started (1e ladder
  holds the GPU).

## Round update (2026-08-19): Q2d-1e GATE FAIL; diagnostic + arm-3 running

- Q2d-1e ladder completed 17:09Z: censored assertions D/E=165 pass,
  value-level reproduction B/C/D/E all True; frozen adjudicator ->
  GATE FAIL at every M1 level A-E (correct dc dz 0.585/0.588/0.589/
  0.410/0.489, sp <= 0.195; best negatives beat correct at C/D/E);
  NC1/NC2 fail as required. Q2D1E_GATE.json + Q2D1E_LADDER.json committed
  (751a05b).
- Q2D_TERMINAL_SUMMARY.md written (stageX_csc_signal/): full 1b->1e
  chain, adjudication verification checklist, failure-mode classification
  (optimization failure prime candidate; diag decides), frozen downstream
  authorization, terminal slot pending the diagnostic.
- Frozen span-param diagnostic LAUNCHED 17:12Z (runner_diag.py, PID 55308,
  job uhq, watcher 3dz auto-adjudicates on exit; estimate ~8-10 GPU-h).
  This is the ONLY successor; no further synthetic stages will be created.
- P-line arm-3 ordinary FT REAL training launched on CPU (3 seeds, 4
  threads, job w61; GPU reserved for the diagnostic; ~2-3 CPU-h).
- CIIP data censuses: KiRHub -> DATA BLOCKER (stageCIIP_kirhub_census_
  20260819, prereg 319d505a...); Davis panel -> INSUFFICIENT ALONE
  (67 WT-variant pairs, median 33 common ligands, 7 held-out parents <
  10; stageCIIP_davis_census_20260819, prereg 2dd8b708...); combined
  census (Anastassiadis, Duong-Ly) next.

## Round update (2026-08-19): combined census SUFFICIENT; device audit

- Combined census (stageCIIP_combined_census_20260819, prereg
  6952ed1a...): SUFFICIENT — Duong-Ly panel ALONE admits CIIP-1A and
  CIIP-1B: 70 single-mutant WT-variant pairs, median 183 common ligands,
  12 held-out parents (>=10 frozen), 6 multi-mutant rows excluded and
  counted, endpoint = % inhibition (never relabeled pK/Ki/Kd), NA 0.23%.
  Anastassiadis = cross-endpoint replication only; cross-panel pooling
  of Kd with % inhibition FORBIDDEN.
- Device audit: _diag16.py (PID 46544, started 2026-08-18 14:16Z, CPU
  355s total, WS 2MB, no output files in the repo) — idle orphan from a
  previous window; taskkill/Stop-Process failed from this shell (no
  task instance — integrity-level mismatch), LEFT RUNNING as-is per
  'no blind kill'; recorded. Runner processes: diag PID 19360 (job 6l9,
  watcher ddf), arm-3 CPU PID 53968 (job w61) — both legitimate.
- GPU metrics mid-diagnostic (2026-08-19 17:45Z): util 32%, mem
  4219/8188 MiB, power 10.44 W, temp 47 C. Step-time measurement and
  post-run metrics will be recorded when the diagnostic ends.

## Round update (2026-08-19): unified potential authored; Stage-1 gated on Stage-0

- Unified deployable function frozen in code: f_theta(P,L) = b_P(P) +
  b_L(L) + s_theta(P,L), s = alpha(P)^T psi(L) (rank 8, hid 64), one
  module for scalar potential / protein contrast / ligand contrast /
  CIIP double contrast / zero-shot endpoint / few-shot correction.
  Free pairwise predictors exist ONLY as diagnostic arms (cycle test
  fails for them by design). 10 structure tests green (identity exact,
  antisymmetry exact, cycle ~1e-6, one-s contrasts, permutation/batch
  contracts, no target-ID input, all branches non-zero grad, stable
  seed). tools/research/stageCIIP_potential_bridge/ (commit 0958e6b).
- CIIP-1A prereg frozen (sha 31d3eeaf...): Duong-Ly within-parent
  capacity; pair-level 60/20/20 split stratified by parent; Q0B-admitted
  single-point pairs only; centered mutation contrast target; 9 arms
  incl. global-compression diagnostic + free pairwise (never a
  mechanism); gates sp>=0.30, dead-zone sign-acc>=0.65 (dead zone 10 %
  units), gaps vs family-shuffle/ligand-only >=0.05 with parent-cluster
  bootstrap; single seed screen then 3 seeds.
- Stage ordering: NO CIIP-1A execution before the Q2d Stage-0 archival
  (diagnostic in final repro checks; watcher ddf armed).

## Q2d TERMINAL VERDICT (2026-08-19 19:20Z) — family CLOSED

- span-param diagnostic: GATE FAIL at all levels (correct dc dz
  0.669/0.544/0.549/0.630/0.508 < 0.70; repro all True; cens 165/165;
  NC1/NC2 fail as required). Per the frozen interpretation: basic
  optimization failure of the low-rank bilinear learner under the
  frozen budget. Learner family CLOSED; no further synthetic
  successors; Q2d-2/Q2d-3/Saifudeen-B1 unauthorized; Core Task 1
  UNRESOLVED. Q2D_TERMINAL_SUMMARY.md is the terminal archive.
- Stage 0 (Q2d archival) is COMPLETE. Next authorized stages in order:
  Stage 1 CIIP-1A (Duong-Ly, prereg 31d3eeaf..., module + 10 structure
  tests ready) -> Stage 2 CIIP-1B -> Stage 3 BindingDB Potential
  Bridge -> Stage 4 few-shot Potential Transport. P-line arm 4-7
  screening can now use the GPU (arm-3 already done: Tanimoto ahead).

## Stage 1 collapse audit (2026-08-19): FAIL explained, not biological

- tested one-hot potential: FAIL (3/13 nonconstant; the other 10 pairs
  have IDENTICAL WT/variant KLIFS inputs -> contrast forced to 0).
- biological protein-conditioned signal: UNRESOLVED (no falsification).
- primary cause: representation (38/65 pairs dP=0; Q1 klifs ns) +
  objective competition (R_g=1081, C_g=-0.016; contrast never left the
  zero floor). Potential capacity not implicated.
- authorized successor: preregistered 2x2 {KLIFS one-hot, local ESM} x
  {joint, centered-only}; freeze new prereg before ANY training.
- CIIP-1B / BindingDB Bridge / production integration: NOT AUTHORIZED.

## 2x2 root-cause diagnostic (2026-08-19): historical run state

- Authorized: freeze matched-subset 2x2 prereg; structure/coverage tests;
  single-seed 4-cell root-cause diagnostic; controls + 3-seed only after
  interpretable results (separate authorization).
- NOT authorized: CIIP-1A PASS claim, CIIP-1B, BindingDB bridge,
  production model/scripts, oracle local ESM as deployable representation.
- Frozen gates for collapsed cells + effects; no single-cell attribution.

## 2x2 reviewer verdict #2 (2026-08-19): design limits accepted; run restarted

- Upheld: prereg SHA match, matched 49-pair subset, oracle naming, no
  production changes, no frozen-result overwrites.
- Design limits recorded (final report MUST respect them):
  1) KLIFS cells structurally capped: only 3/9 covered test pairs have
     nonzero KLIFS input -> 5/9 nonconstant gate unreachable by
     construction; KLIFS collapse = structural statement; objective
     main effect not fairly estimable on KLIFS; interaction confounded.
  2) Objective-sampling confound: KLIFS joint L_abs uses WT+variant rows
     (pool ~9.1k cells), ESM joint variant rows only (~5.8k) -> not a
     pure factorial estimate; per-cell pool/valid-label/WT-variant
     counts + step-1 R_g/C_g + epoch L_abs stats now reported.
  3) Coverage selection bias: 4 whole families missing; conclusions
     restricted to the oracle-covered subset.
  4) Allowed interpretations only: KLIFS structural deficit; whether
     oracle ESM restores nonconstant potential; whether centered-only
     relieves absolute-loss dominance. No universal objective claim,
     no CIIP-1A PASS, no deployable-representation claim.
- Run state resolved: first relaunch had crashed after cell 1 (var_dec
  shape bug 44 vs 64); fixed + reporting fields added (protocol
  unchanged); relaunched (job 4gq). Commit ad28ce9.

## 2x2 reviewer verdict #3 (2026-08-19): three P1s fixed; status corrected

- Corrected status: launch 3 (job 4gq) trained all four cells then
  crashed in effect_boot LOPO (str-vs-int index); exit code had been
  masked by the grep pipeline. State was 'execution incomplete /
  process absent' - now fixed and relaunched from scratch (job x2k,
  no pipes, PYEXIT captured). Stale watcher replaced (v2).
- P1 abs-sampling repetition FIXED (amendment A1): the abs rng stream
  now includes the minibatch index; previous launches re-consumed the
  same <=512 absolute cells ~12x per epoch (inflated L_abs). Split/
  seed/update count/loss weights/gates unchanged; prereg untouched.
- P1 point estimate FIXED (A2): observed_pair_mean_effect (frozen
  primary estimand) + observed_parent_mean_effect reported; bootstrap
  mean explicitly NOT a point estimate; status rule uses observed
  point + CI.
- P1 LOPO FIXED (A3): leave-one-parent-out sign stability now computed
  for all five effects.
- Regression tests added (test_2x2_impl.py): unique-per-minibatch and
  deterministic abs sampling, observed-point semantics, LOPO coverage,
  determinism, status rule. 36/36 green. IMPLEMENTATION_AMENDMENTS.md
  frozen (A1-A4) with code SHA f333c979... Commit d256beb.

## Control-arm stage (2026-08-19): historical launch state; completed 2026-08-20

- Goal: attribute the oracle ESM nonconstant response - real mutation-
  centered protein info vs annotation/generic-context/main-effect
  shortcuts. Only the matched 49-pair covered subset; centered-only
  objective for all arms; identical budget/checkpoint/bootstraps.
- Deliverables completed: CONTROL_RESULT.json, CONTROL_REPORT.md,
  annotation shortcut audit, per-pair metrics, gradient coverage,
  no-leakage audit, tests, SHA256SUMS.
- Final status vocabulary only: ORACLE_LOCAL_SIGNAL_SUPPORTED /
  NOT_SUPPORTED / UNRESOLVED + deployable NOT VALIDATED + CIIP-1A/1B/
  Bridge/production NOT AUTHORIZED.

## CIIP contextual-propagation audit (2026-08-20)

- Completed under frozen preregistration `cdd6e0a8...e9b5b19` in
  `tools/research/stageCIIP_context_propagation_20260820/`. All 49 mutation-
  erasure pairs used identical masked inputs and produced maximum absolute ESM
  embedding delta 0.0 (tolerance 1e-5).
- Original residue states show context propagation: mean site delta 4.0111,
  radius-6 mean 1.1605, non-site context mean 0.05749, full-sequence mean
  0.07392. A random ESM window can therefore carry nonzero edit information.
- No labels or predictor fitting were used. Context-only predictive value and a
  site-specific predictive residual are NOT EVALUATED; CIIP-1B, deployable
  router work, BindingDB bridge, and production integration remain unauthorized.

## CIIP-1A control-arm adjudication (2026-08-20)

- Formal verdict: **ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED**. Correct and random
  local windows are both nonconstant on 9/9 test pairs; correct-minus-random
  pair-mean R2 = -0.1217, parent bootstrap [-0.4569, +0.0327].
- ESM mutation sensitivity remains a representation fact, not ligand-conditioned
  interaction evidence: correct-site delta norm 0.5310 vs random 0.0267, with
  all 49 covered pairs larger at the true site.
- The result is restricted to oracle coordinate, Duong-Ly centered % inhibition,
  49 covered pairs, current low-capacity potential, and pair-level split. CIIP-1B,
  BindingDB bridge, and production integration remain unauthorized. The
  read-only contextual-propagation audit is complete and does not authorize a
  successor. No further CIIP training is currently authorized.

## CIIP-2 independent research cycle (2026-08-20, user-mandated)

- A user-issued phased mandate authorized a full independent cycle: Phase 0
  audit, Phase 1 diagnostic designs, Phase 2 candidate comparison across four
  domains, Phase 3 preregistered successor, Phase 4 gated smoke. Master report:
  `report/research_ideas/ciip/CIIP2_RESEARCH_REPORT_20260820.md`. Recommended
  mainline: OLR-Potential (ligand-conditioned residue router with mean-pool
  skip + cross-fitted ligand nuisance + assay-gain weights, unified under one
  deployable s_theta(P,L)).
- Audit additions: shared-ligand-pattern baseline R2 = 0.1313 on covered test
  pairs (random-window CIIP-1A arm was riding it); parent-shared vs
  mutation-specific variance 134.8 / 89.7 %^2; sibling LOSO ceiling 0.293;
  family prior -0.021; 99/183 WT-ceiling ligands.
- Successor stage `stageCIIP2_olr_potential_20260820` (prereg a7b17e8a...,
  ADD-1 aa8d06af..., ADD-2 91e2cb3a...): 12/12 structural tests green;
  instrument qualification INSTRUMENT_UNDERPOWERED (planted transferable field
  recovers +0.016..+0.030 R2 vs 0.25 standard, robust across estimators);
  real-data smoke gate (b) FAILED (C-perm 0.1818 ranks above all correct
  arms; best increment +0.042 R2). Phase 5 not executed per frozen chain.
- Terminal verdicts: R1 representation SUPPORTED; R2 identification NOT
  SUPPORTED at this power; R3 deployment UNRESOLVED (power); R4 few-shot bridge
  BLOCKED; R5 binding interpretation NOT CLAIMED. Overall **UNRESOLVED
  (power)**: the 49-pair panel cannot adjudicate deployable protein-conditioned
  interaction learning in either direction at the pre-registered standard.
  Successor unblock: >=100 mutation conditions across >=30 parents on one
  endpoint, or a Ki/Kd DeltaDeltaG corpus with its own prereg. Production
  `model/` and `scripts/` untouched.
