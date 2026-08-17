# MetaSieve research history

This is the compact authoritative chronology after the 2026-08-16 cleanup.
The pre-consolidation full log remains recoverable from Git commit `361c342`
and earlier history. Deleted legacy implementations are indexed in
`archive/README.md`.

## Governing task

Predict BindingDB Ki for unseen protein targets at k=0/1/2/3/5 support sizes.
Episodes contain one recipient target; support and query ligands are disjoint.
The active development protocol is protein-component-hard CD-HIT40 with a
sealed double-cold meta_test (logical exclusion after parsing). No query labels, closed-form adaptation,
inner-loop optimization, test-time gradients, or fabricated protein-ligand 3D
coordinates are allowed.

## Legacy cycle, through 2026-08-14

Analytic QPSMP/LIRMS, HyperSAR, D-MEMT/DORM, CIPF/TERM and K3/ELMT were built and
falsified. Recurring failures were k=1 structural death, uniform or label-insensitive
support routing, weak protein specificity, dead confidence branches, residual/level
semantic ambiguity and insufficiently identifiable primitive slots. The useful
lessons are consolidated in `report/meta_fewshot/LEGACY_PRE_R0_SUMMARY.md`; the
working artifacts and dedicated code were removed on 2026-08-16.

## R0-R4: governed protocol and representation ladder

- R0 rejected the selected ligand-retrieval prior on exact-free cold targets and
  exposed selection/transduction caveats in the earlier 12.3% k0 result.
- R1 established the double-cold split and the meta_test seal.
- R2 showed representation, not another loss form, was the decisive question.
- R3/R4 built the level-shape family and established A0, B3 and later C2 as the
  k0 MSE/CI Pareto set rather than a single dominant model.

## R5-R8: transport and shape-first training

R5 repaired donor pools, whitening scope, gradient accounting, seals and artifact
contracts. R6-R7 tested multiplicative, additive and reliability-gated relative
transport; query-specific gates were deployment-inert or harmed calibration.
R8 showed shape-first pairwise/cliff-aware training can improve within-target
shape and activity-cliff sign, but the gain traded against CI and failed promotion.

## R9-R14: objective family closure

- R9 cliff weighting produced C2, a Pareto arm, but no decisive all-metric win.
- R10 variance objectives, R11 grammar-shape, R12 margin ranking and R13 direct
  shape readout failed their gates.
- R13's final gate record is 18 collected, 16 pass, 2 expected failures. The
  implementation was removed after consolidation; its leaf report remains.
- R14 exact decomposition found A0 has the best zero-shot ordering floor
  (shape 0.692, correlation 0.213). Ranking-primary G1 reduced correlation to
  0.134. Regression-compatible ListCE was structurally valid but contributed
  only about 1.7% of the MSE gradient and was inert. The loss-form axis is closed.

## Current evidence boundary

No trained model met the preregistered excellence threshold. Development k0 MSE
frontier is approximately 2.055-2.149 with CI 0.531-0.580. The best development
activity-cliff sign is 0.782 on Pareto-dominated C1. Fixed Tanimoto residual
transport remains the strongest reproducible query-specific k>=2 mechanism, but
it is ligand-only and cannot establish protein-conditioned meta-learning.
The double-cold meta_test population defined for the current protocol remains
sealed; older protocol and pre-authorization computations are separately
quarantined and must not be described as current confirmation.

## 2026-08-16 repository consolidation

Removed the obsolete FORT/theory/retired archive trees, the standalone research
prototype tree, pre-R0 experiment directories, closed relative/locality/direct-shape
models, their trainers and dedicated tests. Retained R0-R14 leaf evidence, current
Pareto loaders, governed data builders, A0/B3/C2 model paths and R14 diagnostics.
Recovery commits and the deletion policy are recorded in `archive/README.md` and
`docs/REPOSITORY_CLEANUP_20260816.md`.

The final layout pass consolidated experimental research, maintained tests and
ignored third-party/runtime utilities under `tools/`. The obsolete `LLM/`, root
`test/` and root `tests/` surfaces were removed. Root project/context/data
summaries were merged into `README.md`, `task.md`, `dataset/README.md` and the
report authorities. `main.py` explicitly dispatches admitted `scripts/`
workflows only; successful research must move from `tools/research/` into
`model/` and `scripts/` before entering that command surface.

## 2026-08-16: A2 readiness, and the closure of the A2 family

Two no-training audit cycles under `tools/research/`. The second was run
because the first drew causal conclusions from correlational evidence; two of
its four load-bearing claims did not survive.

**Governance first.** `QPSMPData.include_meta_test` had defaulted to `True`
since R5 while the R5 contract recorded the opposite. Six analysis scripts took
the default and parsed the sealed cells; one wrote a seal claim its own code
did not implement into seven artifacts. No recorded number is affected, and
this is demonstrated rather than argued: re-running the identical audit under
the repaired seal reproduces 105 of 105 numeric fields bit-identically. The
seal is now fail-closed, opening it requires a written authorization, and
artifacts derive their seal block from the dataset object. A second defect
surfaced during the repair — the record audit was classifying ten double-cold
stage summaries as older-protocol evidence on a corpus none of them names.

**A2 is closed.** The plan's premise was that frozen A0 contains a
protein-conditioned SAR coordinate identifiable from k<=5 labels. Probing the
internal representations A2 actually consumes — not the endpoint scalar the
first cycle measured — shows every representation's ligand-differential
survives a wrong-protein swap with cosine 0.998-1.000, against a protein-blind
reference of exactly 1.000. A trained low-rank readout gains **+0.017** from the
real protein over a capacity-matched permuted one on `embed`, and loses on
`mean_state` and `max_state`. The family is rejected, not amended: four
representations spanning the trunk were tested, including the two that retain
the most protein-differential.

**The diagnosis it rests on is now controlled.** The protein-inert ordering
holds across five distance-stratified donors (level 0.215 -> 0.342 pK, centered
0.0007 -> 0.0011 pK), with a measurement floor of exactly zero and working
shuffled-label, foreign-ligand and scrambled-protein controls.

**Two earlier claims are withdrawn.** The first cycle's "the architecture can
express what training removes" came from one random initialisation; ten of them
move in uncorrelated directions (pairwise cosine -0.003 over 1845 pairs),
unaligned with truth, producing no usable ordering. And the collapse is not in
the readout: causal interventions on the two protein channels put it at the
fusion and atom pooling inside `ContactGrammar`, where a 47% change in the
retrieved context becomes a 0.31% change in `mean_state`.

**One new positive result.** `embed` carries a protein-*independent*
transferable SAR direction, delta-affinity `r` +0.2623 [+0.1295, +0.4055] on
held-out protein components from 1,553 parameters, better than the raw ligand
encoder. It is ligand-side transfer, and it has not been compared against the
fixed Morgan/Tanimoto transport that already exploits the same continuity.

**Standing language constraint.** The protein path is exactly invariant to
residue-slot permutation, so no result from this architecture may be called
pocket-aware or biologically localized.

## 2026-08-16 (later): the exact A2 operator, and the ligand-side result

A third cycle, because the second closed A2 from the wrong object. The v2
probe rejected A2 using a zero-shot bilinear pair predictor, which shares A2's
feature space but forms no support moment and has no shrinkage in k.
`tools/research/a2_exact_probe/` implements the operator the plan specifies and
tests it on real episodes.

It passes 19 structural gates, including one A0 cannot: its k=1 correction is
genuinely query-specific, where A0's is provably a pure level shift. Then it
fails 5 of 6 preregistered performance and control gates with resolved paired
intervals. k=5 MSE 1.1765 against a two-scalar level baseline's 1.0746 and
fixed Morgan/Tanimoto's 0.9101. A wrong protein (1.0866) and shuffled support
labels (1.0820) both make it *better*: both falsification controls fail
inverted. Its query-specific content is 0.0028 pK against a 0.884 pK label
spread, while the same operator given noise features produces 0.3497 pK - so
the mechanism is intact and the moment simply carries nothing. **A2 is closed
on its own operator and its own gates.**

Stage L then tested the one positive result against the comparator it had never
been measured against. The protein-independent `embed` direction predicts the
**signed** within-target affinity gap at r +0.270 [+0.128, +0.418] on held-out
components. Tanimoto scores +0.028 on that quantity - structurally, since
similarity is symmetric and a signed gap is not - but wins on the *magnitude*
(+0.288 vs +0.161). The two predictions correlate +0.026 and `embed` is
unchanged when residualised against Tanimoto. On activity cliffs `embed` is
+0.379 where Tanimoto is -0.370. The signal is real, orthogonal, and
**pairwise** - which is exactly why the A2 moment form, which averages the
support into one vector, cannot use it.

Governance was finished first: the seal is now correctly described as logical
exclusion after parsing rather than a physical seal, `violations = 0` no longer
masks the open process-unsealed incident (the audit exits non-zero), five of
seven split-undeclared artifacts were resolved from their own leaf runs and two
recorded as unresolvable, 28 future-dated records were corrected, and a
specification for a genuine physical seal was written but not implemented or
authorized.

## 2026-08-16 (Stage P): objective-only protein conditioning is closed

Six matched runs, two arms, three seeds, 1,200 steps each. `A0repro` is the
incumbent retrained; `CPCoverdrive` is identical except that the protein
counterfactual is computed on the within-target *centered* prediction, at four
times the weight, on every episode including k=0. The configs were verified to
differ in exactly those two fields before any metric was computed.

**The primary gate fails.** Correct-protein k=0 ordering changes by
-0.0066 [-0.0545, +0.0417] against a requirement of a positive lower bound and
a mean of at least +0.05.

The decisive measurement is not the gate. In **both** arms, at **every** k,
substituting a similarity-matched wrong protein leaves within-target ordering
unchanged in the third decimal: A0repro 0.156/0.156 at k=0 and 0.334/0.334 at
k=5; CPCoverdrive 0.149/0.148 and 0.331/0.331.

The mechanism was not at fault. Gradient of the centered contrast into the
level branch is 8.1e-07 - float32 zero - so `protein_value(P)` was excluded
exactly as the algebra requires, while gradient into `embed` and
`interaction_head` rose 4.6x and 3.6x over the incumbent. The objective even
made the protein response *reproducible*: seed-to-seed cosine of the
protein-induced shift is +0.316 for CPCoverdrive against A0repro's -0.059, the
undirected signature a random initialisation shows. But its alignment with
centered truth is +0.022 against a +0.10 threshold. The model was made
consistently sensitive to the protein in a direction that carries no affinity
information.

Stage P existed to separate "no objective ever asked" from "the data does not
contain the signal". The first is now excluded. The evidence shifts toward the
second without proving it: a different architecture, budget or protein
representation could still succeed.

Closed at exactly that scope: **centered-objective training on the current
ContactGrammar, at a 1,200-step budget on the double-cold protocol with
sequence + 2D inputs, does not produce protein-conditioned within-target
ordering.** The stop rule was applied - the admission stage did not run and
CPCpos/CPCwrong/CPCrand/A3perm were not trained.

## Next authorized work

**None.** No training is authorized. `NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md`
is superseded for its A2 content, and the centered-protein-counterfactual
training innovation, though structurally valid, has lost its premise and its
target and is not authorized either.

`tools/research/a2_readiness_v2/PREREGISTRATION_V2.md` specifies **Stage P**:
whether protein-conditioned within-target ordering is learnable at all on this
protocol when an objective demands it - the one question that separates "no
objective ever asked" from "the data cannot support it". Frozen, donor rule
verified meta_train-only. Two arms (A0repro, CPCoverdrive) x 3 seeds,
~2.3 h measured; `tools/research/stageP_cpc/PREREGISTRATION.md` is
authoritative.

A second, independent lane is now evidence-supported: a **pairwise** learned
operator over (query, support) that could use the Stage L direction, where the
moment form cannot. Untested, and it must beat Tanimoto transport to matter.

M0/MSA remains an independent calibration diagnostic.

---

# 2026-08-17 — Stage A: inner/outer-loop meta-learning permitted and screened

## The governing change

The user instruction of 2026-08-17 **supersedes the repository's prior
prohibition on inner loops and deployment-time support adaptation.** Inner/outer
loops and differentiable support adaptation are now permitted. Unchanged and
still prohibited: ridge regression, pseudoinverses, closed-form solvers,
query-label adaptation at inference, and multi-stage pretrain/finetune regimes.
Support **inputs and labels** may drive adaptation at inference; query labels
remain loss-and-metric-only at every k.

Also fixed by that instruction: exactly one public supervised DTA dataset per
experiment — here the governed BindingDB-Ki double-cold protocol — with no
merging and no cross-dataset support, retrieval, normalization, labels or
checkpoints. Davis and KIBA, if ever used, must be trained independently from
scratch in separate experiments.

This is **AdaMBind-inspired framework testing, not a reproduction of AdaMBind**,
and no artifact in this stage may claim otherwise.

## Stage 0 audit (`tools/research/stageA_innerloop/AUDIT_DATAFLOW.json`)

The ten-step flow was traced and measured rather than read. Steps 6, 7 and 9
(support inner-loop update, adapted query forward, task scoring) **did not
exist** before this stage.

Confirmed by measurement: **0 violations in 400 draws** for one-target-per-task,
support/query cell disjointness, ligand disjointness and within-side ligand
uniqueness; and the frozen evaluation banks are nested — support is a true
prefix chain, the query panel is identical across k, and support/query ligands
never overlap.

Two audit findings shaped the design:

1. **The trainer selects checkpoints on `meta_val` labels.** The training
   gradient never touches meta_val and the label scale is fitted on meta_train
   (verified numerically), but `train()` keeps the state with the best meta_val
   admission score. Disclosed rather than repaired, because repairing it would
   stop `A0` from reproducing the accepted baseline. It is identical in all
   three arms, so it cannot manufacture a between-arm difference — but **every
   meta_val figure in this stage is an optimistic development estimate, not a
   held-out one.**
2. **The incumbent `transport` already performs label-based few-shot.** The
   inner loop fits the same support labels, so the two mechanisms partially
   substitute: as adaptation absorbs the support residual, transport's `locked`
   residual shrinks. This is a genuine attribution confound, which is why the
   no-adaptation control and the pre-versus-post decomposition are mandatory
   reports here rather than extras.

## Adaptable scope

`interaction_head.2.{weight,bias}` — **97 parameters, 0.0054% of the 1,798,833
trainable.** The smallest subset that can reorder ligands rather than only shift
their level, and the weight/bias split is the instrument for the k=1
shape-versus-level question: bias alone is exactly a scalar level shift.

Rejected: `contact_weight` (24 params, cannot express ligand-specific shape),
`interaction_head` in full (14,017, would confound adaptation with capacity),
and full-backbone second-order MAML (excluded until partial adaptation is shown
structurally incapable).

Inner step size **0.1 at 1 step**, chosen on meta_train component folds against
the frozen A0 checkpoint (`INNER_LR_SELECTION.json`): held-out MSE 0.9282 vs
0.9612 unadapted. `lr=1.0` diverges to 1578 pK² at 3 steps, so the grid brackets
the useful range. `meta_val` was not read for this choice.

## Two identities that make the comparison meaningful

One code path serves all three arms, so `--inner-steps 0` is the accepted
recipe. Both halves of that claim are tested rather than asserted:

* the re-implemented readout equals `InteractionGrammarModel.encode`'s own
  endpoint **bitwise** at k ∈ {0,1,3,5};
* the Stage A episode loss at `inner_steps=0` matches the production episode
  loss to **2e-5 on all seven terms** at every k, with each term verified
  non-zero so the match is not vacuous.

The second test exists because the first implementation silently dropped three
auxiliary terms (`support_match`, binding contrastive, protein contrast). Had
that survived, any A1 gain would have been "the inner loop recovers what we
deleted".

## Bugs found and repaired in this stage

1. **The inner loop could not run under `torch.no_grad()`.** Adaptation *is* a
   gradient computation, and evaluation calls it inside `no_grad`. Fixed with an
   explicit `enable_grad` scope plus detachment of the fast weights on the way
   out. Without this, support adaptation at inference was impossible.
2. **Auxiliary loss terms dropped** — see above.
3. **`hash()` reintroduced in the evaluator** for a permutation seed. Caught and
   replaced before any run; this is the same PYTHONHASHSEED defect repaired in
   Stage R and it must not return.
4. **`paired()` compared a condition against itself**, which made every
   counterfactual an exact zero that looked like a clean null. Fixed to take
   separate left/right conditions, with a regression test.
5. **Two trainers wrote to one output directory.** Killing a shell wrapper's
   Python child let the wrapper walk on to the next arm while a fresh launch used
   the same path. A fail-closed guard now refuses to start on an existing
   `progress.jsonl` without `--force`.

## Result: NOT PROMISING on the conjunction; the two mechanisms split

`tools/research/stageA_innerloop/REPORT.md` and `RESULT.json` are
authoritative. Four of six preregistered gates pass.

**The inner loop (`A1`) is a weak positive that did not resolve.** It beats the
accepted baseline at every k on every metric — k=0 MSE 2.0579 vs 2.0753,
Pearson 0.2069 vs 0.1735, and MSE gains of 0.1111 / 0.0980 / 0.0443 / 0.0206 at
k=1/2/3/5 — and **every interval crosses zero**. Only k=1 and k=2 exceed the
recorded 0.058 retraining spread. `A0`'s k=0 MSE lands inside the incumbent
band, so the baseline reproduced.

**The task selector (`A2`) is rejected.** Worse than `A1` on MSE, Pearson,
Spearman and CI at *every* k. The diagnosis is not noise: the selector
concentrated on candidates with support/query gradient cosine **+0.9897**
against a **+0.6555** population mean — that is, on the tasks where support
already predicts query. Gradient agreement measures redundancy, not
informativeness. It cut effective diversity to 5.29 of 9 candidates, cost
2.33x the encoder forwards, and erased the protein specificity `A1` had.

**Two results worth keeping regardless of the verdict.**

*Support labels genuinely matter.* Every wrong-support counterfactual is
**resolved**: permuted support costs +0.36 to +0.43 pK^2 at k>=2, and
matched-wrong support costs +0.81 to +2.05 across k. This mechanism is
label-bound, not a calibration artifact.

*`A1` produced the first resolved wrong-protein gap in this project's record* —
+0.0188 [+0.0052, +0.0327] at k=2, +0.0177 [+0.0079, +0.0282] at k=3,
+0.0085 [+0.0037, +0.0137] at k=5. Small, but every prior wrong-protein gap in
R0-R14 and Stage P crossed zero. `A2` erases it to four decimals.

**k=1 is a level shift, not shape adaptation.** The 97-parameter scope splits
exactly, so one run answers it: bias-only (a pure scalar shift) recovers **81%**
of the k=1 gain, and the shape residual is +0.0214 [-0.0913, +0.1446],
unresolved. The shape share grows with k (19% -> 41%) but never resolves.

**Why the gain shrinks with k rather than growing (gate G5 failure).** The
incumbent `SimilarityTransport` already shifts predictions by Tanimoto-weighted
support residuals. The inner loop fits the same labels, so as adaptation absorbs
the residual, transport's `locked` term shrinks and transport does less. The two
mechanisms substitute. This was preregistered as a confound and is now observed.

**Why an inner loop cannot be bolted onto a trained model.** At k=1 the step has
a closed form: `alpha = 2*lr*||h||^2`, and the post-step residual is `(1-alpha)`
times the pre-step residual. Measured: **`A0` alpha = 1.514, overshooting on
100% of episodes**; `A1` 0.241, `A2` 0.216, neither overshooting. Overshoot
flips the residual sign each step, which predicts an alternating sweep — and
`A0`'s observed sweep alternates exactly: 1.5352 / 2.8626 / 1.4349 / 2.2049 at
0/1/2/3 steps. **Training with the loop is what makes the loop
well-conditioned**, pushing alpha 6.3x down into the stable region.

**Cost.** The inner loop is *free in encoder terms*: 6,480 encoder forwards for
both `A0` and `A1`, because the 97-parameter scope sits downstream of the
encoder and re-evaluates on cached features. Only the selector costs more
(15,120, 2.33x). Zero added parameters.

## Next authorized work

**Three matched seeds of `A0` vs `A1` alone.** Not the selector, and not
additional architecture. The single-seed screen cannot separate a consistent
0.02-0.11 pK^2 improvement from retraining noise, and that is the only open
question the mechanism raises. The measured obstacles are named and should not
be "fixed" by complexity: competition with the existing Tanimoto transport
explains the shrinking gain with k, and a one-point inner loop on 97 parameters
explains why k=1 is mostly a level shift.

The A2 selector is closed. A learned bi-level task adapter is **not** authorized
— the preregistration permitted it only if the simple rule produced credible
evidence, and it produced a clean negative with an explained cause.

---

# 2026-08-17 — Stage B: complementary meta-adaptation rejected, with a measured cause

`tools/research/stageB_complementary/` — REPORT.md and RESULT.json are
authoritative. Stage A's artifacts are preserved unmodified; `CORRECTION_AUDIT.md`
records eight defects in Stage A's analysis, all corrected and each guarded by a
regression test.

## The decision

**The AdaMBind-inspired target-task meta-adaptation framework is NOT admitted to
the production model.** Nothing was promoted to `model/` or `scripts/`. Three
preregistered stop conditions fired: the complementary arm does not beat both
constituent mechanisms; the improvement is only level calibration; ranking
degrades while MSE improves.

## Finding 1: the biggest effect in the cycle is an evaluation leak

`Tleak` is the baseline with one change — checkpoints selected on `meta_val`,
the rule Stage A and the incumbent trainer use — on the *same* fit components.

**k=0 MSE 2.7425 (leak-free) vs 2.1246 (meta_val-selected): -0.6180.** Resolved
at k=2, 3 and 5. k=0 Pearson 0.0561 vs 0.1557.

Against Stage A's A0 (2.0753) the 0.6672 gap decomposes as **0.6180 from the
selection leak (93%)** and 0.0492 from the 12% smaller training set. The leak is
**5.6x the largest mechanism effect measured in this cycle** and **10.7x the
recorded retraining spread**. Every meta_val figure from the standard trainer,
including the recorded incumbent band, is optimistic by roughly this margin.

## Finding 2: the ligand representation collapses within a target

Mean pairwise cosine between the readout hidden vectors of one target's query
ligands: **T 0.99859, M 0.99784, H 0.96309, C 0.99659.** The collapse is not
created by the readout MLP — `embed`, its input, has the same relative spread
(0.033 vs 0.034).

The consequence is algebraic. A weight update `dw = -lr * sum_i c_i h_i` moves
query q by `<dw, h_q>`; with `h_i ~ h_j` that is the same number for every query,
i.e. **a level shift**. And if the target is mean-zero, `sum_i c_i h_i ~ 0` and
**the adapter produces nothing**.

That single measurement explains every arm: M's correction is 99.7% level
(0.6542 vs 0.0021 centered shape), H's is 99.7% level (0.1482 vs 0.0004), and C
— forbidden from level by construction — is inert.

## Finding 3: C improves both metrics, and it is still a rejection

C is the only arm in this project's record to improve MSE **and** ranking
together with resolved intervals (Spearman +0.049 to +0.092, resolved at every
k). But **none of it comes from the adapter**:

- C's meta correction is 0.0000 / 0.0006 / 0.0434 / 0.0389 at k=1/2/3/5;
- removing it changes nothing or slightly helps (+0.0004 / -0.0360 / -0.0159);
- **the C-minus-T ranking contrast is bitwise identical at k=0 and k=1**
  (+0.052852 / +0.092001 / +0.036328), and at k=1 the meta term is exactly zero
  and the transport is a pure level shift, which cannot reorder anything.

So C's advantage over T is a **zero-shot trunk difference from training with the
complementary objective**, not a few-shot mechanism. Recorded as a training-time
observation deserving its own test, and explicitly not as evidence for
meta-adaptation.

## Finding 4: the Stage A support-label counterfactual was never evidence for adaptation

With the corrected pre-adaptation anchor and every arm measured, the *baseline*
shows the same large gap (T matched-wrong minus correct: +2.41 at k=1, +5.89 at
k=5). The incumbent transport is label-driven by construction. The correct
statistic is incremental dependence, and for C it is **negative at every k**
(-0.21 to -0.79): the candidate depends on correct support *less* than the
baseline.

## Corrections applied to Stage A

Eight, all tested. Two changed measured numbers: the matched-wrong control was
anchored to each arm's own post-adaptation prediction (so A0 and A1 faced
differently corrupted adversaries), and A0's counterfactual rows were computed
but dropped from the summary. Six were interpretation errors, including the
withdrawn claim that A1 improved every metric at every k — **its k=5 CI is
0.62951 against A0's 0.63140** — and the conditioning formula, which omitted the
adapted bias (`alpha = 2*lr*(||h||^2 + 1)`, so A0's 1.514 becomes 1.714) and was
wrongly said to predict query MSE when the query effect is
`-2*lr*r_s*(h_q . h_s + 1)` and varies per query.

## Next authorized work

**None on this framework.** No rescue by attention blocks, full-backbone MAML,
learned selectors or extra datasets: the measured obstacle sits *upstream* of the
adapter, and no adaptation rule downstream of a collapsed representation can
repair it.

The two measurements worth carrying forward are the leak (remove meta_val
checkpoint selection from every future comparison) and the collapse (within-target
ligand cosine 0.997), which links Stage R's inert operator, Stage P's
uninformative protein response and Stage L2's weak directional signal to one
upstream cause.

---

# 2026-08-17 — Stage C: the feasibility boundary for MSE <= 1.00 pK^2

`tools/research/stageC_level_shape/BOUNDARY.md` is authoritative; leaf artifacts
are `FEASIBILITY.json`, `LEVEL_CEILING.json`, `REPRESENTATION.json`. No training
was done in this stage. Baseline is the leak-free Stage B `T` checkpoint.

## The error is mostly level, and a perfect level predictor would meet the target

MSE decomposes exactly into `level^2 + centered_MSE`:

| k | MSE | level^2 | centered | level share | MSE with perfect level |
|---|---|---|---|---|---|
| 0 | 2.7425 | 1.8664 | 0.8761 | 68% | **0.8761** |
| 1 | 1.8549 | 0.9788 | 0.8761 | 53% | **0.8761** |
| 2 | 1.3628 | 0.5560 | 0.8068 | 41% | **0.8068** |
| 3 | 1.2481 | 0.4500 | 0.7981 | 36% | **0.7981** |
| 5 | 1.0096 | 0.2754 | 0.7342 | 27% | **0.7342** |

So the <=1.00 target is arithmetically reachable, and reachable through level
calibration alone. Support labels are already the level mechanism: level^2 falls
1.87 -> 0.28 across k while the shape term barely moves. **Every "few-shot gain"
this project has measured is target-level calibration arriving through labels.**

k=0 centered MSE is 0.8761 against a within-target label variance of 0.8525 —
ratio **1.0277**, i.e. the endpoint orders *worse* than predicting each target's
own mean.

## Target level is not predictable from anything tested

Predicting an unseen protein's mean pK, meta_val, decay chosen on meta_train
component folds, meta_val read once:

| method | level MSE | vs calibrated constant |
|---|---|---|
| ESM linear probe | 6.5368 | 4.85x worse |
| ESM nearest neighbour | 5.1292 | 3.81x worse |
| sequence-length probe | 2.5283 | 1.88x worse |
| meta_train grand mean | 2.1703 | 1.61x worse |
| incumbent model | 1.7078 | 1.27x worse |
| ESM MLP probe (best) | 1.6357 | 1.21x worse |
| *calibrated constant (reference)* | *1.3471* | *1.00* |

**Not one method beats a constant.** The decay sweep selected the largest value,
driving the fit toward a constant because the features carry no cross-component
level signal. There is also a ~0.91 pK covariate shift between the two splits'
target-level distributions.

Since `MSE = level^2 + centered` and `centered >= 0`:

> **k=0 MSE >= 1.6357 with the best legitimate level predictor, even with
> perfect within-target ordering.** Reaching 1.00 needs level MSE <= 0.1239, a
> 13.2x reduction, i.e. explaining ~91% of a between-target variance that
> nothing currently explains at all.

A plausible reason this is irreducible rather than a modelling failure: a
BindingDB target's mean pKi depends on which ligands were tested against it, so
part of "target level" is a property of assay history and library composition
rather than of the protein sequence. Consistent with the measurements; not
proven by them.

## Stage B's collapse claim was too strong — `occupancy` carries real signal

Cosine 0.997 was necessary but insufficient evidence. The ligand-varying share
is 17-33% across representations, over 2-4 effective directions. A frozen linear
probe on within-target centered affinity, fitted on meta_train, read once on
meta_val:

| representation | train-fold r | meta_val r |
|---|---|---|
| embed (96) | +0.2383 | +0.0074 [-0.1550, +0.1718] |
| readout_hidden (96) | +0.2301 | +0.0256 [-0.1375, +0.1936] |
| **occupancy (24)** | +0.2029 | **+0.2182 [+0.0751, +0.3670]** |
| section (48) | +0.2174 | +0.0603 [-0.1026, +0.2250] |

**`occupancy` is the only representation whose within-target ordering signal
survives out of component, and its interval excludes zero.** The wide
representations overfit; the 24-dimensional contact-type vector holds.

The model has this signal and does not use it: `occupancy` reaches the endpoint
only through `contact_weight`, a `Linear(24->1)` of the *same capacity as the
probe*, yet the endpoint orders at ratio 1.0277. The likely cause is that those
24 parameters are trained against a total MSE that is 68% level error, so the
optimizer spends them on level. Mechanistic hypothesis, not a measured cause.

Exploiting it fully is worth `1 - r^2 ~ 4.8%` of the centered term, about
**0.04 pK^2** — real, resolved, and far too small to change the k=0 verdict.

## Boundary

- **k=0 MSE <= 1.00: not reachable** with current inputs; bounded below by 1.6357.
- **k>=1: reachable only with near-perfect ordering**; k=5 is already 1.0096.
- **Missing information: zero-shot target-level affinity calibration for unseen
  homology components.**
- **Present but unused: a resolved ordering signal in `occupancy`.**

## Next authorized work

Not a training schedule and not an adapter — the constraint is information, not
optimization. In priority order:

1. run the preregistered **M0/MSA lane** (`report/meta_fewshot/stageM0_msa_probe_20260816/`)
   or another external protein representation, reported as external data, and
   test it directly against the calibrated-constant reference on target level;
2. test whether target level regresses on **assay/library covariates** within
   meta_train — if it does, the level is partly not a protein property and the
   zero-shot target must be restated;
3. the **`occupancy` shape lever** — separate the level and shape paths so the
   24 contact-type parameters are trained on within-target centered supervision
   rather than level-dominated MSE. This can improve centered MSE, CI and
   Spearman; it cannot reach MSE <= 1.00 and must not be reported as if it could.

---

# 2026-08-17 (later) — Stage D: panel-context level + orthogonal level/shape training

Opened after the Stage C boundary. `tools/research/stageD_level_panel/`.
Preregistration frozen before any arm trained; D0 diagnostics ran first.

## D0 re-audit of Stage C (answers to the five governing questions)

1. The level/shape decomposition is per episode (per draw), not per canonical
   target; the drawn-panel share of level^2 is small (panel-sampling variance
   0.013-0.034 pK^2), so level is a between-target quantity.
2. The calibrated constant (1.3471) reads meta_val labels and is a disclosed
   REFERENCE; the meta_train-only constant is 2.15-2.17, which tested features
   do beat (best legitimate level predictor: ESM-650M linear, 1.6875).
3-5. Level is a joint protein/assay/panel property: within meta_train,
   component identity explains 46% in-fold but -1.1% held-out; document
   identity 70% in-fold but +6.8% held-out (and meta_val has zero document
   overlap); **panel composition transfers best (+23.9%)** vs protein
   sequence (+11.9%). A shuffled-panel control (5.07) confirms the panel
   association is real signal.
4. The occupancy ordering signal survives stratification: scaffold-novel
   r +0.154, low ligand recall r +0.221, meta_val overall +0.203; per-component
   values heterogeneous but positive on average.

New external representation recorded: ESM-2 650M pooled embeddings computed
locally (tools/runtime/esm2_t33_650M_pooled/, manifest with model/weights
hashes), reported as external data.

## Stage E (preregistered, running)

Two innovations: I1 panel-set level readout (framework); I2 orthogonal
level/shape routing (training). Arms T2 (leak-free baseline), T2-LEVEL
(loss-only ablation), LSP (both innovations), LSP-NOROUTE (framework-only
ablation). 1,200 steps, seed 20260815, internal checkpoint selection. Gates
G1-G6 and stop rules S1-S4 in PREREGISTRATION.md. Mandated GPU verification
(torch.cuda.is_available, model/batch devices, nvidia-smi utilization) runs
before every arm.

---

# 2026-08-17 (evening) — Stage D/E result: panel-set level head + orthogonal routing REJECTED

All four Stage E arms ran to completion (1,200 steps, seed 20260815, leak-free
internal checkpoint selection, GPU verification before each arm). Frozen
meta_val banks, read once. Authorities:
tools/research/stageD_level_panel/REPORT.md, LSP_vs_T2.contrast.json,
LSP_PANEL_SHUFFLE.json, per-arm RESULT.json.

T2 (leak-free retrain of the incumbent recipe) reproduces the band: k=0 MSE
2.5961 (level^2 1.7314, centered 0.8648), k=5 0.9859 (Spearman 0.314, CI
0.619, cliff 0.609), honest label/protein controls. T2's k=5 is already
below the 1.00 target.

LSP: k=0 MSE 2.3935 — gain -0.2026 [-0.5195, +0.0869] UNRESOLVED; k=5
+0.2269 [+0.0685, +0.4098] RESOLVED DEGRADATION; ranking (Spearman/Pearson/
CI) lower at every k (unresolved). G1 and G2 fail; stop rule S1 fires.

Attribution (ablations T2-LEVEL, LSP-NOROUTE): the panel level head carries
the k=0 level gain (level^2 -0.32) but double-fits the target mean with the
transport at k>=1; the trained head is inert (own level error 1.438 vs 1.539
panel-shuffled, versus the D0 frozen probe's 1.887 vs 5.075); orthogonal
routing improves centered only 0.915->0.878 (baseline 0.865) — level
contamination of the contact dictionary was NOT the shape bottleneck; the
level term routed into the shape path (T2-LEVEL) damages level to 2.636
while slightly helping centered (0.807).

Nothing promoted. The measured frontiers now stand as: zero-shot level is
assay-history-dominated (70% in-fold document variance; 6.8% out-of-document
transfer; 23.9% panel composition; 11.9% protein sequence) and no legal
representation family tested (ESM-150M/650M pooled, panel composition, assay
covariates, trained panel head) approaches the 0.1239 level budget; the
shape term is representation-limited (within-target information r ~ 0.22).
Structure lane audited: only 15/499 governed targets have an exact-sequence
holo structure in the local pilot20k corpus — insufficient for a full-protocol
pocket lane; MSA lane remains blocked on a governed UniRef snapshot.

---

# 2026-08-17 (night) — Stage F: pairwise learned transport REJECTED

tools/research/stageF_pairwise/. Both arms ran (1,200 steps, seed 20260815,
leak-free internal selection, GPU verification). Against the frozen T2
baseline every F-minus-T2 interval crosses zero at every k; k=5 MSE
-0.0076, k=5 Spearman -0.0190, k=5 CI -0.0124; the framework-only arm
(F-ABS) is worse at every k. Gates G2/G3 fail; nothing promoted.

This is the fifth learned-kernel family in the record that fails to beat the
fixed Morgan/Tanimoto weighting, now including the Stage L pairwise
direction as input. The accumulating ledger supports the representation
limit: the trunk's ligand-varying subspace carries r ~ 0.22 at most, and no
downstream operator can extract more than the representation contains.

---

# 2026-08-17 (night) — Stage G: ESM-650M residue-input trunk (single-seed screen)

tools/research/stageG_esm650/. The incumbent similarity_only recipe retrained
with the local ESM-2 650M protein bank (1280-dim pooled + 128-slot residues,
recorded provenance). Single seed, leak-free internal selection, GPU
verification.

First arm in the record to improve MSE, level^2, centered MSE, Spearman,
Pearson, CI and cliff sign at EVERY k against the frozen T2 baseline.
k=0 centered MSE resolves: -0.0396 [-0.0772, -0.0018]. k=0 MSE -0.2136 and
k=5 MSE -0.0417 do not resolve; ranking gains do not resolve. Controls clean
(matched-wrong/permuted/wrong-protein all above correct at every k).
Cost: 2.04M trainable (1.14x T2), 172.6 s, 516.7 MB peak.

Preregistered gates: G1 PASS, G2 FAIL (no resolved MSE gain at k in {2,3,5}),
G3 PASS, G4 PASS, G5 partially (parameter criterion as written). Per the
stop rules the lane halts at the screen; because the failure mode is
statistical power on a new input lane rather than mechanism inertness, the
continuation is the newly preregistered multi-seed confirmation
(PREREGISTRATION_G2.md) — a new stage, not a tuning step.

---

# 2026-08-17 (night) — Stage G2: ESM-650M lane NOT CONFIRMED; boundary restated

tools/research/stageG_esm650/REPORT_G2.md and G2_multiseed_contrast.json:
three fixed seeds for both arms (T2 and G), pooled component bootstrap. The
single-seed Stage G pattern did not survive: k=0 centered MSE pooled
-0.0035 [-0.0244, +0.0168] (the single-seed resolved gain was a seed
artifact); k=0 MSE -0.2078 [-0.7001, +0.2128]; k=5 MSE -0.0206. Every
interval crosses zero. G2-1 and G2-2 fail; no meta_test opened.

Per-seed bands: T2 k=0 2.458-2.981, k=5 0.946-1.007; G k=0 2.239-2.790,
k=5 0.944-0.987. The retraining spread dominates every single-seed
difference — the reason the multi-seed protocol is mandatory.

Decisive new measurement (D0b_DOC_TRANSFER.json): within-document assay
history transfers level across targets with R^2 +0.451, and the double-cold
split's document closure makes exactly that signal unavailable at
inference. The final boundary is restated in
report/BOUNDARY_20260817_NIGHT.md with the full falsification ledger.

---

# 2026-08-17 (night) — Stage H0: structure/pocket lane audited and rejected at the identifiability gate

tools/research/stageH_pocket/. Local MMseqs2 search of the 387 governed
meta_train/meta_val targets against the pilot20k holo corpus: 152 targets
have a holo homolog at >=90% identity (exact-sequence only 15), and 209/387
have one at >=30% identity with >=50% query coverage - so the earlier
15/499 exact-coverage figure understated the lane.

Pocket descriptors (pocket atoms within 6.0 A of the holo ligand: size,
volume, amino-acid and element composition, holo ligand size, mapping
identity/coverage; up to 5 complexes averaged per target) were extracted for
209 targets. Frozen SGD probes, component-fold selection, meta_val read
once (H0_POCKET_IDENTIFIABILITY.json):

- grand-mean baseline 2.6179;
- pocket MLP 2.4398, linear 2.9446, shuffled-pocket control 2.4941;
- pocket + ESM-150M MLP 2.3447.

The pocket association is ~6.8% over the constant and barely above the
shuffled control: at this resolution the structure lane carries no
level-breakthrough signal, and the preregistered identifiability threshold
was not met, so no Stage H training is authorized. The lane is recorded as
REJECTED_BY_FROZEN_DISCRIMINATOR for level calibration (pocket priors for
shape/ordering remain unmeasured and are not claimed).

---

# 2026-08-17 (night) — Stage I: live ESM-150M LoRA lane REJECTED (ranking-only)

tools/research/stageI_lm/. Arms I (LoRA trainable, r=8) and I-FROZEN (live
path, adapters frozen), 1,200 steps, seed 20260815, leak-free internal
selection, GPU verification. G2 fails: no resolved MSE gain at any k
(k=2 -0.0459, k=3 -0.0383, k=5 -0.0182, all crossing zero). S1 fires;
nothing promoted.

Resolved positives (observation, not promotion): I minus I-FROZEN k=2
Spearman +0.0744 [+0.0118, +0.1420]; I minus T2 k=3 Pearson +0.0553
[+0.0007, +0.1205]. Trained adapters slightly worsen level (k=0 level^2
+0.113 unresolved). Engineering: chunked LoRA backward on long proteins
silently OOM-killed the process; the first-chunk gradient bound fixed it and
is documented for any future LM lane. The live frozen encoder beat the
cached fp16 bank in this seed (k=0 -0.259 resolved, cause undisclosed).

---

# 2026-08-18 — Stage J: assay-aware level head + paired level alignment REJECTED

tools/research/stageJ_assay/. D0c first: journal/publisher codes parsed from
panel_ids carry level signal (probe 1.619 vs 2.155 constant, shuffle 2.522;
100% of meta_val episodes share a journal code with meta_train). Three arms
(J, J-NOPAIR, J-NOJRNL), 1,200 steps, seed 20260815, leak-free internal
selection, GPU verification.

The level head works at k=0 (level^2 1.7314 -> 1.297-1.321; MSE -0.3941
[-0.8918, +0.0312], the largest k=0 level improvement on record but
unresolved) and the journal covariates are load-bearing (1.50 without vs
1.30-1.32 with). The paired alignment term adds nothing measurable. But G3
fails with RESOLVED ranking degradation (k=2 Spearman -0.0624, k=2 CI
-0.0282, k=3 Spearman -0.0598) and G2 fails: coupling a learned zero-shot
level head to the k>=1 transport again substitutes worse calibration for
the support labels' near-optimal one and shrinks the transport's shape
residuals. REJECTED; nothing promoted.

---

# 2026-08-18 — Stage K/K2: contrastive coembedding (K-REG) — strongest mechanism, NOT confirmed

tools/research/stageK_contrastive/. Screen: episodic InfoNCE (K) fails;
positive/negative regression alignment (K-REG) passes the single-seed gates
with resolved k=2/3/5 MSE and k=0 centered/Pearson gains. Three-seed
confirmation: ALL five k MSE improvements resolve pooled (k=0 -0.1118
[-0.1851, -0.0490]; k=1 -0.0480; k=2 -0.0273; k=3 -0.0218; k=5 -0.0122,
all hi < 0) with ranking preserved and zero control inversions in any
seed — the first all-k resolved mechanism in the project's record.

But K2-2 fails: the k=0 centered gain did not survive pooling (-0.0154
[-0.0304, +0.0001]), so per the preregistered stop rule the configuration
is NOT CONFIRMED: no meta_test opened, nothing promoted to model/ or
scripts/. The mechanism is read as calibration-consistency regularization
(collapse 0.99859 -> 0.9908) rather than a new information source; k=0 MSE
stays >= 2.44 in every seed and the bounded conclusion stands.

---

# 2026-08-18 (night) — Stage L: support-gated level head REJECTED (G3)

tools/research/stageL_gated/. The gate worked for MSE (k>=1 statistically
indistinguishable from T2; k=0 MSE 2.0997, the best calibration in the
record, level^2 1.2151) but ordering degrades with RESOLVED intervals at
k=2/3/5 (Spearman -0.087/-0.075/-0.060; CI -0.040/-0.035/-0.026) and k=0
CI -0.030: training the head on k=0 episodes reshapes the SHARED trunk,
which then orders worse at k>=1 even with the head off. Three compositions
of a learned zero-shot level head with this trunk (E ungated, J ungated,
L gated) have now all failed the ranking gate — the zero-shot level
objective and the within-target ordering objective conflict on the same
representation. A separate frozen-feature calibrator has a measured
ceiling (1.62 -> k=0 ~2.5) below L's 2.10. The bounded conclusion in
report/BOUNDARY_20260817_NIGHT.md stands as the programme's final state.

---

# 2026-08-18 (night) — Stage M0: ChemBERTa ligand-side LM probes REJECTED at identifiability

tools/research/stageM_chemberta/. Frozen ChemBERTa-77M pooled ligand
embeddings (local snapshot, 600-dim, manifest recorded). Within-target
ordering probe: meta_val r +0.1472 [-0.0261, +0.3179] — crosses zero and
sits below the occupancy record (+0.218). Level probe: collapses to the
grand mean exactly (2.1547). The ligand-side language-model input family is
falsified for both level and shape. The external-representation ledger now
covers every locally testable legal family; the bounded conclusion in
report/BOUNDARY_20260817_NIGHT.md stands as final.

---

# 2026-08-18 (night) — final boundary audit: the conclusion re-derives exactly

tools/research/stageN_audit/. FINAL_BOUNDARY_AUDIT.json and AUDIT_REPORT.md
re-derive every load-bearing number of the boundary from the raw evaluation
rows with no training: T2 k=0 decomposition (2.5961 / 1.7314 / 0.8648,
66.7% level share), within-document level transfer R^2 +0.4515, and the K2
three-seed pooled MSE contrasts (-0.1118 / -0.0480 / -0.0273 / -0.0218 /
-0.0122, all hi < 0) — all matching their stored authorities bitwise.
104 result artifacts audited for the meta_test seal: none evaluated; only
the two disclosed legacy R14 artifacts record included=True. All seven
training stages are preregistered. The bounded conclusion is verified
end-to-end and stands as the programme's final state.

---

# 2026-08-18 (night) — programme close: FINAL_STATE document written

report/FINAL_STATE_20260818.md consolidates the closing summary: the
bounded conclusion (BOUNDARY_20260817_NIGHT.md), the falsification ledger
(EVIDENCE_LEDGER.md) and the verification (stageN_audit/AUDIT_REPORT.md).
Open lanes are recorded as externally blocked (MSA: no UniRef snapshot;
structure: 209/387 coverage; Davis/KIBA: promotion-gated, not run).
Governance state: meta_test 0 evaluations, query labels loss-only,
7/7 stages preregistered, 135 tests passed.

---

# 2026-08-18 (night) — Stage P0: protein function annotation (GO) family REJECTED at identifiability

tools/research/stageP_go/. ProteinKG25 (477k proteins, 47k GO terms)
sequence-matched to 313/387 governed targets (81% coverage). GO-bag level
probe: 2.2699 vs 1.4329 grand-mean baseline on the covered subset — 58%
worse than the constant; fold-selected decay 1.0 (max shrinkage). Protein
function annotations carry no transferable cross-component level signal,
consistent with component identity transferring -1.1%. The
external-representation ledger is now complete for every locally available
legal family; the bounded conclusion stands.

---

# 2026-08-18 (night) — round 9: closing document reconciled with the final ledger

report/FINAL_STATE_20260818.md updated to include Stage P0 (GO annotations,
313/387 matched, falsified) in the external-representation ledger and the
final verification tallies (maintained suite 268/6 via the sanctioned
entrypoint; research suites 147 passed with RUN_SLOW=1; boundary audit
bitwise). All authorities (task.md, history.md, EVIDENCE_LEDGER.md,
BOUNDARY_20260817_NIGHT.md, FINAL_STATE_20260818.md, GOAL_ACTIVE.md) are
mutually consistent. The programme remains in its final state; the goal
record stays active only for continuation under new external data.

---

# 2026-08-18 (night) — Stage Q: decoupled frozen-feature level head REJECTED (G3)

tools/research/stageQ_frozenhead/. Q0 probe met its threshold (joint frozen
features 1.3416, best frozen level predictor on record). The trained
composition fails G3: k=0 MSE -0.1479 (unresolved) with RESOLVED ranking
degradation at k=0/2/3 (Spearman -0.059/-0.050/-0.050). The decoupling
hypothesis is falsified: even a head with zero trunk-coupled features, gated
to k=0 only, reshapes the shared trunk through the k=0 training signal
itself. Four compositions (E, J, L, Q) have now failed this gate; the
level/ranking conflict on one shared trunk is fundamental to single-stage
end-to-end training, and the only escape (a separately trained inference
calibrator) is excluded as a multi-stage regime. Q-UNGATED confirms the
gate is necessary but insufficient (k=1 MSE 2.3095 while k=5 Spearman
0.3445/CI 0.635 are record-best). The bounded conclusion stands.

---

# 2026-08-18 (night) — literature round 2: external validation of the level wall

1. Nelen et al., Matched pairs demonstrate robustness against inter-assay
   variability. J Cheminform 17(1):8, 2025. doi:10.1186/s13321-025-00956-y.
   PMID 39833966. Abstract-level evidence, independent of this project:
   "absolute values from different assays are rarely comparable" while
   potency DIFFERENCES between matched compound pairs are more consistent
   (agreement within 0.3 pChEMBL 44-46% uncurated, 66-79% curated). This is
   the same level/ordering structure the boundary measured on BindingDB:
   the level term is assay-dominated and non-transferable, the ordering
   term transfers. Recorded as external validation of
   report/BOUNDARY_20260817_NIGHT.md.
2. Xu et al., CrossLinker: Aligning Relational and Sequential Contexts for
   DTI Prediction in Cold-Start and Few-Shot Scenarios. J Chem Inf Model
   66(7):4257-4267, 2026. doi:10.1021/acs.jcim.5c03216. PMID 41874971.
   Link-based (fine-grained) contrastive learning over sequence + a
   relational knowledge graph. Design evidence only: its relational
   modality has no local analogue with measured transferable signal (GO
   bags were falsified for level in Stage P0; component identity transfers
   -1.1%), so no local lane is reopened by this paper.

The boundary stands; the literature basis is updated accordingly.

---

# 2026-08-18 (night) — round 12: documentation closure

report/CURRENT_MODEL_EVIDENCE.md extended to cover the complete cycle
(Stage I through Q plus the audit and external validation), and
tools/research/README.md now carries a stage index of all eleven research
directories with their verdicts. All six authorities (task.md, history.md,
CURRENT_MODEL_EVIDENCE.md, EVIDENCE_LEDGER.md, BOUNDARY_20260817_NIGHT.md,
FINAL_STATE_20260818.md) plus GOAL_ACTIVE.md are mutually consistent and
complete.

---

# 2026-08-18 (night) — round 13: Davis/KIBA boundary-check plan frozen (not run)

tools/research/stageR_daviskiba/PREREGISTRATION.md: a frozen design for
re-testing the boundary's scope on the sealed Davis assets (mechanism-v2
and homology-v1 inventories only; no label read). The check would replicate
the D0 assay-history anatomy, the T2 baseline and the K-REG mechanism under
the same leak-free rules, gated on whether the level wall reproduces
(k=0 level share >= 50%, within-document transfer >> cross-document).
Status: NOT AUTHORIZED, NOT RUN — the standing governance authorizes only
the BindingDB-Ki lineage for training until a candidate passes promotion.
KIBA has no local assets and is not planned here.

---

# 2026-08-18 (night) — round 14: completion-evidence inventory passes

tools/research/stageN_audit/COMPLETION_INVENTORY.json: all 20 artifacts
referenced by the final authorities exist on disk with the expected
schemas; all 8 trained stages carry PREREGISTRATION.md + REPORT.md + JSON
artifacts; zero missing paths. The objective's terminal condition
(scope-bounded conclusion after multi-family falsification) is evidenced,
verified and inventoried.

---

# 2026-08-18 (night) — round 16: method-ladder cycle formally closed

tools/research/method_ladder/CLOSURE_MAP.md maps each of the eight named
ladder families to its measured successor stage and verdict (families 1-8
all closed by Stages A/B, E, F, J, K, Q measurements). task.md's
method-ladder section is updated from paused to CLOSED. No open research
item from the pre-session plan remains: every named family, every locally
available input family and every level-head composition has a recorded
verdict.

---

# 2026-08-18 (night) — round 17: seal-phrasing correction (maintained test caught the overstatement)

The maintained research-record test (test_research_record.py) failed on the
round-12 edit to CURRENT_MODEL_EVIDENCE.md, which overstated a
logical-exclusion seal whose labels were parsed in process.
Corrected in all three authorities (CURRENT_MODEL_EVIDENCE.md,
FINAL_STATE_20260818.md, BOUNDARY_20260817_NIGHT.md) and in GOAL_ACTIVE.md
to the supportable claim: no sealed meta_test label entered any fitting,
selection or reported metric (0 evaluations in the audited artifacts).
Full maintained suite green again: 268 passed / 6 skipped.
