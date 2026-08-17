# Current task contract

Updated: 2026-08-17. This file is the single authority for active work. Detailed
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
representation. Few-shot adaptation is blocked until that representation is
shown to carry transferable, affinity-relevant protein information.

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

Authority for these statements:

1. `report/POST_COMPLETION_REVIEW_20260818.md`
2. `report/BOUNDARY_20260817_NIGHT.md`
3. `report/CURRENT_MODEL_EVIDENCE.md`
4. `report/EVIDENCE_LEDGER.md`
5. `tools/research/stageS_sar_field/REPORT.md`
6. `tools/research/stageT_mmp/CORRECTION_20260817_CORE_KEY.md` (scopes Stage T)
7. `tools/research/stageV_core_mmp/REPORT.md` (corrected Phase-1 test)
8. `tools/research/stageV_core_mmp/STAGE_U_GOVERNANCE_AUDIT.md` (Stage U stopped)

## Active five-phase programme

The phases are sequential scientific gates. A phase may reject with one seed;
promotion requires the stated evidence. A failed prerequisite blocks every
downstream phase. Do not jump from an unidentifiable interaction representation
to a larger meta-adapter.

### Phase 1 - transferable protein-conditioned interaction signal

**Question.** Using only protein sequence-derived features and ligand 2D
features, can the model learn a local interaction representation `Z(P,L)` whose
dependence on the correct protein improves true SAR ordering on unseen protein
components?

**Representation contract.** Preserve ligand atoms or functional groups,
protein residues or regions, and their local cross-interaction tokens until the
readout. Early global pooling is not allowed. A protein embedding changing under
a swap is insufficient; the change must align with affinity differentials.

**Status: NOT ESTIMABLE on this dataset (2026-08-17). Phase 1 is neither closed
nor in progress — it is blocked by identifiability, and no model-level verdict
exists.** This supersedes the earlier "CLOSED" wording, which rested on a
defective transformation key.

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
