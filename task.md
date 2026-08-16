# Active Cold-Target Mechanism-Meta Work Contract

## 2026-08-15 architecture amendment

The mathematical theory is a source of design inspiration, feasibility checks,
information-limit warnings, and scalar-output consistency.  It is not a frozen
neural architecture contract.  The active learned core may be deeply redesigned
when module-level tests and governed cold-target evidence justify the change.

The active candidate is the interaction-grammar trunk with a label-locked
residual transport whose per-support coefficient depends on the query. It uses
no ridge, matrix solve, inner loop, or deployment gradient update. k=0,1,2,3,5
are first-class modes of one single-stage episodically trained model. k=1 is now
structurally query-specific rather than scalar, but that channel has **not**
passed its governed controls; ligand-specific SAR remains unestablished at every
k. The retained BPSF endpoint and its scalar-only k=1 kernel remain available as
the control arm.

An optional sparse Cartesian rank-0/vector/symmetric-traceless-rank-2 encoder is
available only when real coordinate inputs are declared.  Common-frame
protein--ligand edges require a verified complex pose.  Independently framed
protein and ligand coordinates may only be fused after invariant reduction.
The current BindingDB main bank has no coordinate sidecar, so its production
path is the exact sequence+2D fallback and cannot be described as atomic 3D
recognition.

## Question

Can one shared, query-loss-trained QPSMP neural meta-potential use correctly bound support labels to
improve protein-specific Cold Target affinity prediction beyond identical-budget additive, level,
ligand-only, SAR-cut, wrong-protein, shuffled-protein, and design-nuisance controls?

## Allowed Inputs

- source-trained or legally frozen protein residue representations;
- ligand molecular graphs encoded by the shared ligand encoder;
- declared measurement context;
- `k={0,1,2,3,5}` support observations from the unseen recipient target.

Target IDs may index cached tensors but cannot enter the model. Query labels, target memory, and
recipient-specific trainable parameters are prohibited.

## Primary Module

`QPSMPBioModel` is the primary learned path. Query loss must deliver finite gradients to the protein
encoder, ligand encoder, localizer, crossed scalar head, section basis, and neural support adapter.
Closed-form/ridge adaptation is excluded from the active few-shot path.

## Frozen Invariants

- support/query rows are disjoint;
- pair inclusion and orientation are outcome-independent;
- support ordering does not change output;
- k=0 reads no support and is exactly the shared zero-shot endpoint;
- k=1 may select or reweight source-learned mechanisms using an absolute,
  label-bound support residual, but cannot be claimed to identify unrestricted SAR;
- delta and rectangle outputs are differences of the retained scalar endpoint path;
- additive, cross-zero-shot, level, and SAR-adaptation channels are reported separately;
- foreign support changes only the transient SAR state;
- validation episodes are fixed before checkpoint selection;
- component/dependency weighting and consumed-development status are explicit.

## Next Gate

Updated 2026-08-15 after the interaction-grammar stage series.

The active candidate reduced MSE by 12-18% at every k on both the frozen
protocol bank and a wide bank over all 42 eligible meta-test targets, in three
seeds, and produced the first genuinely protein-conditioned zero-shot endpoint
in this lineage (cross-component protein swap moves the output by 0.438 pK
against 0.0093 pK). Governed admission was still refused.

The next experiment must attack the three specific failures, not the MSE:

1. **Ranking.** The concordance index falls from 0.647 (zero-shot / level) to
   0.571-0.610 and Spearman from 0.372 to 0.169-0.257 whenever the
   query-specific gate is active. Any new adaptation channel must improve CI and
   Spearman, not only squared error.
2. **Support identity.** Permuting the support labels leaves `mean(r)` exactly
   unchanged, so the permutation contrast isolates the query-specific channel.
   It is currently negative at k=2, 3 and 5.
3. **Zero-shot resolution.** The endpoint spread across the queries of one
   episode is 0.087-0.186 pK against a 0.93 pK label spread.

Every stage must report both banks. The 6-episode frozen bank cannot resolve
differences below about 0.05 MSE and must not be used alone for selection or
decision; the 400-step Stage 0 probe improved a 6-episode validation score while
worsening every test metric.

### Geometry is settled, not open

`scripts/audit_geometry_coverage.py` establishes that **zero of 17,717**
BindingDB deployment cells have a common-frame protein-ligand complex. Every
Cartesian or equivariant interaction encoder (PBCNet2.0, TensorNet, PaiNN, MACE,
Equiformer, SE(3)-EGNN) therefore has no legal input on this task.
`model/cartesian.py` stays verified and unused, and the active model refuses
coordinate inputs by raising. Reopening this requires **new data**, not a new
architecture: either a second training stream on the invariant holo
contact/distance supervision, or co-folded complexes that do not yet exist here.
Do not fuse independently framed protein and ligand structures.

### Stage 5 outcome: the objective is the blocker, not the operator

`model/relative_grammar.py` (signed antisymmetric reference-query difference
operator plus a leave-one-out label-consistency credit) passed all 17 algebraic
and synthetic gates and was then **rejected on real data**: on `meta_val` its
`full` equals its `level_only` and its permutation gap is identically zero at
k>=2, the algebraic signature of flat weights and a null operator.

The decisive diagnosis came from evaluating both budgets on the same split. At
800 steps the earlier `rho` gate *improves* the concordance index at k=2,3,5; at
2000 steps it *degrades* it in 9 of 12 (seed, k) cells. Two structurally
different transports, one shrinking and one inert, both converge to level
calibration, and the ranking cost grows with optimization.

Cause: the transport is trained on squared error, whose optimum for k noisy
support residuals is shrinkage toward their mean; a level shift is constant
across queries and cannot change ranking. **The next candidate must make ranking
the primary training signal for the transport component itself.** The current
`pairwise_ranking_loss` carries weight 0.5 against a dominant MSE term and is
applied to the whole prediction rather than to the query-specific part.

Two constraints established: adding supervision to a failed mechanism did not
rescue it, and a label-consistency credit cannot be evaluated while the operator
it weights is null.

### Stage 6 correction and outcome: the bottleneck was the similarity metric

The Stage 5 statement above is **over-generalised**. A label-and-chemistry-only
audit shows the MSE optimum is not the support mean: a fixed Morgan/Tanimoto
kernel beats it by 0.19/0.21/0.25 MSE at k=2/3/5, and helps the residuals of a
frozen checkpoint that never saw the mechanism. Those two transports failed
because they found no usable similarity metric, not because squared error
forbids a query-specific channel.

`model/similarity_grammar.py` (`--arch similarity_only`) is **accepted as a
validated mechanism at k>=2**: within-checkpoint, three seeds, complete banks,
18/18 positive point estimates for MSE, CI and Spearman, 16/18 component-level
bootstrap lower bounds above zero, permutation gaps +0.40 to +0.51.

Stage 6 is frozen. The next experiment must change **one** thing:

1. resolve the unresolved cross-arm contradiction (F wins on `meta_val`, the
   incumbent wins on `meta_test`) with identical initialisation and a bank
   neither arm selected on; or
2. sharpen the kernel — `nearest_residual` beats the current soft weighting at
   every k, so `gamma` (which barely trains, 7.99 from an init of 8.0) is too
   soft.

k=0 remains the dominant error term and is untouched; do not bundle it in.

### Stage 7-8 outcome: kernel confirmed, Mac-Diff direction closed

**Stage 7 (accepted).** Freezing each incumbent `grammar` checkpoint and swapping
only the support transport at inference shows the fixed Morgan/Tanimoto kernel
beats the level baseline at k=2/3/5 in MSE, CI and Spearman with 9/9
component-level lower bounds above zero. This resolves the Stage 6 cross-arm
contradiction: it was trunk and training variance, not a mechanism conflict.
Hard nearest-support is **rejected** — its MSE intervals always cross zero and it
significantly degrades ranking. Replacement of the incumbent's learned transport
is still not established (meta_val favours the kernel, meta_test the incumbent,
nothing significant).

**Stage 8 (rejected).** A Mac-Diff-inspired sequence-derived locality prior over
the 128 ordered protein slots passed all 14 structural gates, including exact
zero-gate identity with the accepted baseline, but regressed the preregistered
k=0 target in 3/3 seeds with negative cross-arm point estimates at every k.

Because the Mac-Diff conformational sidecar and the support-conditioned
conformational router were conditional on Stage 8 passing, **neither is
authorised**. No Mac-Diff weights or inference were used. Reopening that
direction requires new evidence that protein representation — not ligand-side
signal or calibration — is the k=0 bottleneck.

### Stage 9-10 outcome: k=0 diagnosed and improved

k=0 error is **59% target-level calibration**, and the trained zero-shot
endpoint has essentially no within-target ligand discrimination (re-centring it
on the true target mean gives 0.7403 against a flat constant's 0.7430; CI 0.525
against a 0.500 coin flip). Ligand retrieval beats the protein-conditioned model
at calibration using no protein at all.

A training-free `meta_train`-only retrieval prior blended into the endpoint at
w=0.5, transport unchanged, **reduces k=0 MSE by 12.3%** with component
bootstrap +0.198 [+0.012, +0.416] and consistent direction in 3/3 seeds. k>=2
gains of 17-23% are seed-consistent but their intervals cross zero and are not
claimed.

### Stage R0 audit: what that result actually supports

Eleven binding corrections were verified by recomputation
(`report/meta_fewshot/stageR0_retrieval_falsification_20260815/`). All eleven
hold. Consequences that bind the next cycle:

* the 12.3% is **development evidence, conditional on ligand overlap**. 305 of
  624 query cells (48.9%) contain a ligand present verbatim in `meta_train`;
  restricted to the 12 exact-free targets the effect is **+0.050
  [-0.074, +0.175]**, unresolved. 6 of 10 components improve;
* selection (`beta`, source, `w`) and inference used the same population;
* Stage 9's composed 25.5% is **transductive** — it re-centres on the query
  panel. The best per-query train-only estimator gives 10.8%;
* the prior is an offline evaluator, not part of the model or checkpoint;
* **protein representation for k=0 is reopened.** Raw pooled ESM cosine spans a
  band of width 0.21 around 0.90 with a 0.024 spread across the nearest 16
  training targets, so `softmax(16*sim)` was near-uniform by construction;
  train-only centring widens that spread to 0.238. Mac-Diff locality, conformer
  routing, PBCNet2.0 and Cartesian equivariance stay closed — they were rejected
  on structural-input coverage and multi-seed training evidence, which this does
  not touch.

### Stage R0 outcome: the retrieval prior is falsified, and Stage 10 is withdrawn

All five preregistered gates fail on the identical Stage 10 population
(`report/meta_fewshot/stageR0_retrieval_falsification_20260815/REPORT.md`):

* **exact-ligand-free k=0 is -0.217 [-0.785, +0.261]** — the prior makes
  genuinely new ligands *worse*, 2.8019 -> 3.0193;
* the entire benefit is exact overlap (1.3581 -> 1.0114) and near-duplicates
  (1.3273 -> 0.9842); every low-similarity stratum regresses;
* tuning the same 200-configuration search on the population it is reported on
  is worth **0.468 MSE by itself** (2.5514 tuned against 3.0193 nested);
* 9 of 10 outer folds select a **protein-blind** source, so protein specificity
  is zero by construction — and the one protein-conditioned fold loses to its
  own shuffled control. Adding the sharper train-centred ESM retriever did not
  help, which lowers the prior on protein-representation interventions again.

**Retrieval is therefore a named baseline only. It is not part of any core
innovation, and no protein-conditioned language attaches to it.**

One exploratory signal survives and sets the target: on 62 exact-ligand-free
activity-cliff pairs (Tanimoto >= 0.6, gap >= 1.0 pK) the trained zero-shot
endpoint orders at **chance, 0.519**, while a parameter-free Morgan/Tanimoto
prior reaches 0.716. The trunk cannot read chemistry where chemistry is decisive.

The binding questions are now, in order:

1. **Build a double-cold protocol** that controls ligand identity, scaffold and
   chemical similarity as well as protein homology. 48.9% exact ligand overlap
   made every earlier development decision incapable of separating recall from
   capability.
2. **Stop the absolute-affinity objective collapsing the interaction trunk into
   a target-level constant** — the calibration/shape decomposition and the
   activity-cliff result agree that this, not any single operator, is the
   mechanism behind every failed transport.
3. **Give the trunk within-target ordering ability** and demonstrate it against
   wrong-protein, ligand-only and permuted-support controls.

### Stage R5-R6 (2026-08-16): contract repairs and the relative-transport cycle

The 2026-08-16 mandate adds two binding requirements on top of the standing
gates: the few-shot correction at k=1 must depend on the complete
(protein, support ligand, support label, query ligand) relation — a scalar
residual shift is no longer a legal k=1 mechanism — and the interaction
branch must measurably contribute to zero-shot ordering instead of
degenerating to a target-level constant. At most two core innovations are
claimed for the final candidate, and the training innovation must be one of
them, proven by a same-architecture ordinary-training ablation.

The experimental contract was repaired first (R5): evaluation wrong-protein
donors come from the same evaluation split with meta_train-only whitening;
gradient cosines are aggregated across episodes/steps/seeds; `meta_test` is
sealed physically (`QPSMPData include_meta_test=False`); every run records
config, split hash, seed, checkpoint sha256, per-target predictions, donors,
activation statistics, gradient coverage and resources.

The current candidate (R6b, in screening) is the relative-transport model:
`f0 = ligand_prior + target_level + mean_m delta(P, L, anchor_m)` with an
antisymmetric protein-conditioned relative potential shared by the
zero-shot shape and the few-shot correction
`t(q) = shrink * sum_k a(q,k) * [r_k + delta_hat(q,k) - (f0(q) - f0(k))]`
(the exact, label-free residual identity). Its Stage 1 gate suite (23 tests)
passes, including the synthetic interaction-branch gate and the private-task
abstention gate. Training is counterfactual gradient-routed shape-first:
pairwise ranking with cliff weighting plus relative supervision as the shape
signal, wrong-protein and wrong-support contrasts, one backward pass. The
R6a screening eliminated the earlier multiplicative-gate design under its
own preregistered gates (gate inert, k=0 17.4% worse than A0 at 300 steps);
the additive correction is the recorded single-variable response. No
real-data performance claim exists for the current design yet; screening
results are elimination-only and the formal gates live in the R7
preregistration.

### Stage R7 outcome (2026-08-16): admission refused; R8 preregistered

The three-seed formal run of the relative-transport design (A2) was refused
under its preregistered gates: k=0 2.420 vs the incumbent's 2.149 (-12.6%,
bootstrap -0.271 [-0.683, +0.091]), CI 0.542 vs 0.580, wrong-protein gap
negative at k=0, level-only beating full at k=1-2. Two mechanism facts are
now established: (1) the shape-first training produces the project's first
real shape gain (0.943 -> 0.895; activity-cliff sign 0.536 vs 0.512) and the
routed level readout converges to the incumbent's calibration at the full
budget (A3: 1.292 vs 1.236); (2) the query-specific rho gate is again
eval-inert (nogate gap ~0.000) while its training disturbs calibration —
the seventh query-specific channel in this project with that signature,
now under ranking-primary objectives, so the objective is no longer the
explanation. Stage R8 is preregistered: A3's configuration exactly, with
stronger shape signal (shape_variance 1.5, relative 1.0) and no
query-specific gate; if three seeds do not reach k=0 1.934, the model
family is closed for the double-cold zero-shot target. meta_test remains
sealed.

### Stage R8 outcome (2026-08-16): family closed for the double-cold zero-shot target

The stronger-shape arm (B1: shape_variance 1.5, relative 1.0, no gate)
reaches k=0 2.167 vs A0's 2.149 (-0.8%, an unresolved tie), with the best
shape term recorded in this project (0.896) and the best activity-cliff
ordering on record (k=5 cliff sign 0.768 vs A0's 0.675) — but CI regresses
(0.535 vs 0.580). Both preregistered advance gates fail, so under the R8
decision rule the model family is **closed for the double-cold zero-shot
target as a claimed core innovation**, and meta_test remains sealed. The
shape-first training is retained as the project's first measured shape
source; the open question for the next cycle is retaining that shape gain
without the CI regression and with better-than-A0 calibration (candidates:
support-conditioned calibration at k>=1 only, LambdaRank-style pair
weighting, budget scaling beyond 1200 steps).

### Stage R9 (2026-08-16): pair-level diagnosis and the cliff-weight dose response

The R9 pair audit (`stageR9_cliffweight_20260816/PAIR_AUDIT_meta_val.json`,
no training) decomposes the R8 CI regression stratum by stratum. At k=0 the
only component-resolved stratum is the **mid-similarity band
(0.4 <= Tanimoto < 0.6): +0.119 [+0.022, +0.220]** per-target sign accuracy
against A0 — the band immediately below the activity-cliff weight's
discontinuity at 0.6. Cliff pairs themselves improve (-0.049, unresolved);
low-similarity pairs are unresolved at the target level (-0.022) despite
their pair-count share; mid-gap pairs (0.5-1.0 pK) are near-resolved
(+0.120 [-0.008, +0.263]). The hypothesis: the x4 cliff weight starves the
0.4-0.6 band. The dose response (cliff_pair_weight in {1.0, 2.0, 4.0},
three seeds, 1200 steps, everything else B1) confirmed and sharpened it:
the x4 cliff weight is a **net negative for the ranking itself** — C1
(weight 1) beats B1 (weight 4) on global CI (0.562 vs 0.535) *and* on cliff
pairs (0.606 vs 0.577 pooled); the cliff-ordering ability comes from the
shape-first training, not the cliff emphasis. C2 (weight 2) gives the
family's first three-seed k=0 below A0 (2.119, unresolved) with the best
calibration of the family (1.218 vs A0 1.236) and k=5 cliff sign 0.775; no
dose passes Z1'/Z5' together. After the weight removal no stratum remains
resolved and margins stay compressed (C1 0.097 vs A0 0.121). R10 tests the
next single variable on the C1 base: `shape_variance_weight 1.5 -> 0.5`,
three seeds, 1200 steps, via the smoke-first stage runner
(`scripts/run_stage.py`).

### Stage R10-R11 (2026-08-16): two falsifications complete the variable ladder

R10 (variance 1.5 -> 0.5 on the C1 base) failed all four gates — the
variance term is not the margin-compression cause. R11 (shape-first routing
on the incumbent trunk, zero architecture change) failed H1-H3: the
incumbent's calibration lives in the interaction branch, so routing the
level away from it degrades calibration (1.236 -> 1.488) and CI
(0.580 -> 0.525). The level/shape routing trades calibration for shape on
every architecture tested. The shape-first gains remain real but
unconverted; the remaining preregistered lever is budget scaling (requires a
matched-budget A0 retrain and the learning-curve condition). meta_test
sealed.

### Stage R12-R13 (2026-08-16): margin loss and direct shape — the ladder closes

R12 (hinge margin ranking on the C2 base) moved CI by +0.003 only — the
margin compression is a symptom of the shape branch's expressivity, not the
loss form. R13 (direct interaction-head shape with difference supervision)
was gate-blocked at Stage 1: the MLP shape branch collapses under the shape
variance term on the synthetic interaction task (thresholds unmoved). The
R9-R13 ladder is now a closed, evidence-consistent chain, and the
consolidated reachable-boundary statement is `report/BOUNDARY_20260816.md`.

### Stage R13.5 (2026-08-16): record audit before any new experiment

Six defects in the record were found and repaired before proposing anything
new. Each is now recomputed from the leaf artifacts by
`python -m scripts.audit_research_record` and held in place by
`tests/test_research_record.py`:

1. **R13 gate count.** "16 gates / 15 of 16 pass" alongside two recorded
   xfails is inconsistent with itself and with the suite. The suite collects
   **18 gates: 16 pass, 2 `xfail`**. Corrected; `RESULT.json` backfilled.
2. **The k=0 "frontier" mixed models.** "MSE 2.055-2.119, CI 0.548-0.580"
   read B3's MSE with C2's and C2's CI with A0's, describing a configuration
   that does not exist. The real Pareto set over (MSE down, CI up) is three
   whole configurations: **B3 (2.055, 0.531), C2 (2.119, 0.548), A0 (2.149,
   0.580)**. No model reaches both ends; no MSE difference against A0 is
   resolved.
3. **The 0.782 cliff sign** is a double-cold **`meta_val` development
   record** on arm C1 — which is itself Pareto-dominated on both primary
   metrics. It has never been measured on `meta_test`.
4. **Two different `meta_test` populations** were being written with one
   name. Stages 4/6/7 report the consumed `bindingdb_ki_main_v0` split
   (42 targets); the double-cold split (22 targets / 10 components) is
   sealed and unopened. The audit classifies all 228 double-cold and
   older-protocol artifacts and finds **0 seal violations**.
5. **R12 had no `REPORT.md` or `RESULT.json`.** Both backfilled from the
   retained comparison artifact; nothing re-run. Its gate **M5 is recorded
   as not evaluable** — the artifact has no `D2_vs_C2` contrast, so the
   preregistered bootstrap against its stated control was never computed.
6. **Stale counts** in `docs/PROJECT_FILE_ORGANIZATION.md` (R5-R10, 394
   tests, 72/79 modules, 212 results) updated to R5-R13 and the measured
   values.
7. **Six checkpoints do not reload into the current model**, and the record
   did not say so. They are the R3R4 pre-fix arms
   (`A1_shared`/`A2_routed`/`A3_full` seed 20260815 and their `_predrift`
   copies), whose `TypedLigandChannels` was replaced by the documented
   capacity fix. This is by design — their `RESULT.json` metrics are the
   evidence, not the bytes — but it is now stated, classified by the audit,
   and covered by a test that fails if any *other* checkpoint is orphaned.

Verified in the `drug` environment: **78 recorded checkpoint sha256 hashes
recomputed, 0 mismatched**; **39 checkpoints reload strictly, 0 broken**,
covering all three frontier arms. Suite tiers:
default 413 passed / 9 skipped (105 s); `RUN_RESEARCH_GATES=1` 416 passed /
3 skipped / 2 xfailed (410 s). The six deferred tests are the synthetic
training gates of the two **closed** families, whose verdicts are immutable
evidence; a new or reopened family must run its own gates in that tier.

### Stage R14 (2026-08-16): ordering localized, then the loss-form axis closed

**Phase 2 (no training)** decomposed the within-target shape term exactly into
an ordering floor `Var(y)(1-r²)` and an amplitude excess `(sd_p - r·sd_y)²`.
Results, all on double-cold `meta_val`, three seeds, eight arms:

* the regression-dominant incumbent has the **lowest ordering floor of every
  arm in the project** (0.692, `r` 0.213); every ranking-primary routed arm is
  worse, 8/8;
* the cause is the **training method, not the architecture** — G1 is the
  incumbent architecture with zero changes, and shape-first training takes
  `r` from 0.213 to 0.134;
* the retained "first within-target shape source" claim is **corrected**:
  B1's shape gain is +0.043 worse ordering offset by 0.060 less amplitude,
  i.e. shrinkage;
* a protein-conditioned amplitude head was **falsified before
  implementation** (global rescale worsens 6/8 arms; the per-target optimal
  scale is negative in 25.2% of targets). It was dropped and **not replaced**,
  so R14 claimed one core innovation rather than two.

**Phase 3** implemented the surviving training innovation — a within-target
listwise term whose optimum coincides with the regression optimum — as a
loss-form swap at fixed weight on the incumbent. 29 gates pass; the alignment
identity was verified numerically first.

**The claim is withdrawn.** O1 fails (floor 0.6931, worse than both A0 arms),
O3 fails (CI 0.544 vs 0.570), and **O4, the necessity control, fails**:
deleting the ranking term entirely (0.6951) is indistinguishable from
replacing it with the aligned one (0.6931). The measured cause is that at
this model's operating point the aligned term supplies **1.7% of the
regression term's gradient** — the `1/T(p)` factor that creates stationarity
also damps it wherever predictions are under-dispersed. Exact
regression-compatibility and useful ranking pressure are in tension.

**Secondary finding.** The incumbent configuration retrained under an
identical setup differs by **0.058 k=0 MSE and 0.051 in `r`**. The whole
frontier spread (2.055-2.149) is 0.094, about 1.6x that. Frontier
differences are one to two retraining standard deviations and none was ever
resolved by a component bootstrap.

**The loss-form axis is closed** (R9, R10, R12, R14 all varied the
within-target ranking objective; none moved the ordering floor). The next
evidence-supported hypothesis is that `r` is bounded by what the ligand
representation carries about within-target ordering — a representation-side
**diagnostic**, not a model change. Not tuned after the fact: the ListCE
weight and shift were fixed in advance and a post-hoc sweep is refused.

meta_test remains sealed and unopened.

## PASS

Proceed beyond development only if preregistered component-level lower bounds show useful full-scalar
gain and correct-protein crossed/SAR specificity, with no support-binding, target-main, document,
panel, or query leakage.

## Stage M0: corrected MSA diagnostic proposal (2026-08-16)

The MSA direction is diagnostic only and does not replace the core training
innovation. See `report/meta_fewshot/stageM0_msa_probe_20260816/PREREGISTRATION.md`
and `report/M0_GLOBAL_RESEARCH_EXPANSION_20260816.md`.

Before execution, D0 must verify a local UniRef snapshot and MMseqs2 or
jackhmmer in the `drug` environment. M0-A uses per-seed A0 residuals and a
fixed low-capacity train-only diagnostic probe; no Ridge, solve, checkpoint
training or deployment adaptation is allowed. M0-B fixes all kernel
bandwidths and normalization inside meta_train. M0-C reports MSA-depth,
popularity, family-overlap, label-count and ligand-novelty strata.

M1 is authorized only after M0-A passes its preregistered increment,
permutation and depth-confounding gates. The training innovation remains the
single core track: calibration-preserving within-target ranking with
counterfactual label-binding controls. `meta_test` remains sealed.

## STOP

Fail closed if the learned arm only beats zero, only separates an artificially destructive wrong
arm, loses to level/additive/ligand-only controls, uses unstable validation selection, or depends on
overlapping development units. Such failure closes the frozen recipe, not the entire function class.
