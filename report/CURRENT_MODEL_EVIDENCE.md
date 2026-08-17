# Current model evidence

Authority date: 2026-08-18. This file is the compact scientific state; leaf
artifacts live under `report/meta_fewshot/stageR0_*` through `stageR14_*`.

Interpretation note: `POST_COMPLETION_REVIEW_20260818.md` controls where the
closing documents use broader language than the measured experiments. In
particular, measured probe coverage is not an information-theoretic upper
bound, repeated level/ranking conflict is not a universal theorem, and proxy
experiments do not directly close their named external method families.

Scope: every conclusion in this file is **BindingDB-Ki double-cold
development evidence**. It does not extend to other DTA datasets, to
architectures that were not run, or to what is achievable in principle.

## Protocol

The task is BindingDB Ki cold-target prediction with k=0/1/2/3/5 support. The
current split is CD-HIT40 protein-component-hard and double-cold. Support and
query belong to one recipient target. Current-protocol meta_test is sealed by logical exclusion after parsing
(fail-closed default plus written authorization); it is not physically
isolated from the read path. Older-protocol and pre-authorization calculations are quarantined and
are not confirmation evidence.

## Retained systems

- **A0**: incumbent similarity-only interaction trunk and zero-shot reference.
- **B3**: level-shape Pareto arm with the best observed development k0 MSE.
- **C2**: cliff-weighted Pareto arm.
- **Fixed Tanimoto transport**: strongest reproducible k>=2 query-specific
  comparator; inactive at k0 and degenerate at k1.

The true development k0 Pareto set is B3 (MSE 2.055, CI 0.531), C2 (2.119,
0.548), A0 (2.149, 0.580). No arm dominates both calibration and ranking.

## Positive evidence

1. Fixed Morgan/Tanimoto support weighting improves MSE and ranking together at
   k>=2 and reacts strongly to label permutation. It proves transferable ligand
   SAR continuity, not protein-conditioned meta-learning.
   **Correction (2026-08-16): every wrong-protein control in R0-R14 is computed
   on uncentered error and therefore measures target-*level* specificity, not
   ordering. R3R4's "first resolved protein specificity" (+0.4216 at k=2) is a
   level result.** The ordering version of the same control measures
   −0.0002 [−0.0015, +0.0008].
2. Shape-first training is the first learned route that measurably improves
   within-target shape and activity-cliff sign. The best cliff sign is 0.782 on
   meta_val development, but the responsible C1 arm is Pareto-dominated.
3. Attention-pooled target level can converge to incumbent calibration under a
   complete budget, but it has not produced a decisive new zero-shot frontier.

## Closed directions

- QPSMP/LIRMS analytic adaptation: rejected and incompatible with the no-solver
  research objective.
- HyperSAR, D-MEMT/DORM, CIPF/TERM and ELMT: correct-support binding or protein
  specificity was absent, k1 often collapsed to scalar calibration, and several
  routing/confidence paths were dead or deployment-inert.
- Seven query-specific learned gate variants: failed across multiplicative,
  additive and reliability forms.
- Ranking-loss substitution: RankNet variants, variance, margin, grammar shape,
  direct shape and regression-compatible ListCE failed to improve the joint
  MSE/ranking frontier. R14 showed ListCE was only ~1.7% of the MSE gradient.
- Full Cartesian/PBCNet2.0 transfer: no legal common-frame input exists for any
  of 17,717 DTA cells. The implementation is algebraically tested but unused.

## R14 diagnosis

The exact k0 decomposition separates target calibration from within-target
ordering. A0 supplies the best ordering floor observed (shape 0.692,
correlation 0.213); ranking-primary G1 reduces correlation to 0.134. Therefore
the next change must test representation identifiability, not another ranking
loss. Earlier apparent shape gains were partly shrinkage/calibration effects.

## A2 readiness v2 (2026-08-16): the A2 family is closed

`tools/research/a2_readiness_v2/`, no training, frozen A0, `meta_val` scored
once after every choice was frozen on `meta_train` folds.

- **Where the ordering comes from.** `r_full` 0.213, `r_ligand_only` 0.027,
  `r_interaction` 0.221; increment **+0.1855 [+0.0566, +0.3236]**, resolved.
  Ordering is an interaction-branch property, and it is concentrated in the two
  least novel ligand terciles (0.289/0.295 vs 0.051).
- **It is not protein-conditioned.** Across five donor strata the level moves
  0.215 → 0.342 pK while the centered ordering moves 0.0007 → 0.0011 pK, a
  260-330x ratio. The correlation contrast is a decisive null at four strata.
  Measurement floor exactly 0.
- **A2 is closed on its own operator.** The v2 rejection used a zero-shot
  bilinear pair predictor, which is not A2's operator; that was an over-reach.
  Stage R implements the exact episodic form (`z = A_phi(e0)`,
  `c_S = mean_i r_i z_i`, `delta = eta(k)<c_S, z_q>`), passes 19 structural
  gates including a genuinely query-specific k=1 correction, and then fails
  5 of 6 preregistered gates on real episodes with resolved paired intervals.
  k=5 MSE 1.1765 against a two-scalar level baseline's 1.0746 and fixed
  Tanimoto's 0.9101. **A wrong protein (1.0866) and shuffled support labels
  (1.0820) both make it better** — the two falsification controls fail
  inverted. Its query-specific content is 0.0028 pK against a 0.884 pK label
  spread; given noise features the same operator produces 0.3497 pK, so the
  mechanism works and the moment simply carries nothing.
- **A protein-independent signed-SAR direction exists in `embed`**, measured
  under a *fair* construction (Stage L2, which supersedes Stage L entirely).
  Balanced (i,j)+(j,i) pairs make a symmetric score's signed correlation
  identically zero by construction, verified at -0.0000. `embed` scores
  **+0.2119 +/- 0.0112** across three pair-sampling seeds, beating the
  capacity-matched Morgan-difference control (-0.042, which *can* carry
  direction and carries none), the protein-blind ligand encoder (+0.119) and a
  random directional head (+0.054); incremental value over the ligand encoder
  is resolved at **+0.188 [+0.052, +0.325]** with the slope fitted on
  meta_train. Zero shared ligand identities and zero shared scaffolds with
  meta_train.
- **Two Stage L claims are withdrawn.** The activity-cliff figure reverses from
  the reported +0.379 to **-0.1178** (188 pairs, 18 targets, 10 components -
  adequate, so confirmatory), and the claim that Tanimoto "points the wrong
  way" on cliffs is void: a symmetric similarity has no direction, and the
  earlier number was an artifact of unbalanced pair sampling. The signal is
  concentrated in the two least-novel ligand terciles (+0.240, +0.245) and
  nearly vanishes in the most novel (+0.052). It is therefore strongest where
  fixed Morgan/Tanimoto transport is already strong and fails where a
  directional mechanism would earn its keep.

## Stage P (2026-08-16): objective-only protein conditioning fails

Two matched arms, three seeds, 1,200 steps: `A0repro` (incumbent, uncentered
contrast at 0.5) and `CPCoverdrive` (centered contrast at 2.0 on every episode).
Configs verified to differ in exactly those two fields.

- **Primary gate P1 fails**: `r_correct(CPC) - r_correct(A0repro)` at k=0 is
  **-0.0066 [-0.0545, +0.0417]** against a requirement of a positive lower
  bound and a mean >= +0.05.
- **The decisive number is not P1.** Correct and wrong protein give the *same*
  within-target ordering in both arms at every k: A0repro k=0 0.156/0.156,
  k=5 0.334/0.334; CPCoverdrive k=0 0.149/0.148, k=5 0.331/0.331.
- **The mechanism worked.** Gradient into `protein_head` from the centered
  contrast is 8.1e-07 (float32 zero), so the level branch was excluded exactly
  as designed, while gradient into `embed` and `interaction_head` rose 4.6x and
  3.6x over the incumbent.
- **The objective made the protein response reproducible but uninformative.**
  Seed-to-seed cosine of the protein-induced shift: A0repro -0.059 (undirected,
  the random-init signature), CPCoverdrive **+0.316**. Its alignment with
  centered truth is **+0.022** against a +0.10 threshold.
- No material regression: k=0 MSE +0.018, CI -0.019, Spearman -0.035,
  calibration +0.028, all unresolved.

**Scope of the closure**: centered-objective training on the current
ContactGrammar, at this budget and protocol, does not produce
protein-conditioned within-target ordering. Not a claim that
protein-conditioned architectures are impossible.

**What it resolves**: Stage P separated "no objective ever asked" from "the data
does not contain the signal". The first is now excluded - the objective asked at
4x weight through a verified gradient route and the model complied with a
consistent protein response carrying no ordering information. The evidence
shifts toward the second without proving it.

## Unresolved questions

1. ~~Does frozen A0 encode a low-dimensional protein-conditioned SAR
   coordinate?~~ **Answered: no**, on four representations spanning the trunk.
2. ~~Can k<=5 support residuals identify that coordinate?~~ **Moot** — there is
   no coordinate to identify.
3. ~~Is protein-conditioned within-target ordering learnable at all when the
   objective explicitly demands it?~~ **Answered: not by an objective alone on
   this architecture.** Stage P asked at 4x weight on every episode and got a
   reproducible but truth-unaligned protein response and no ordering gain.
4. ~~Does the protein-independent SAR direction in `embed` add anything to
   Morgan/Tanimoto?~~ **Answered: yes, orthogonally** - it predicts the signed
   gap where Tanimoto predicts only the magnitude, and it is strongest on
   activity cliffs where Tanimoto inverts. **Open:** whether a *pairwise*
   operator (not a moment) can convert that into a deployment gain; Tanimoto
   transport still wins every episodic comparison run so far.
5. Can a training objective preserve MSE calibration while making correct
   support binding causally necessary?
6. Does MSA information explain protein-side target-level calibration residuals
   beyond frozen ESM, after depth and research-bias stratification?

## Next decision

**No training is authorized.** `NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md` is
superseded for its A2 content: the model innovation it proposes is falsified
before implementation. The centered-protein-counterfactual training innovation
is structurally valid but has lost its measured premise and its located target,
and is **not authorized**.

`tools/research/a2_readiness_v2/PREREGISTRATION_V2.md` specifies Stage P, the
one prerequisite that would separate "the objective never asked" from "the data
cannot support it". Its donor rule is verified meta_train-only by construction
(`draw_episode` selects donors from the episode's own split), it is frozen, and
it costs **~23 min per 1200-step run (measured 1.14 s/step), 2 arms x 3 seeds = 6
runs, ~2.3 h**.
`tools/research/stageP_cpc/PREREGISTRATION.md` is the authoritative arm and
cost statement; the five-arm figure previously recorded here counted the
later admission-stage controls and is superseded.

The independent M0/MSA probe is unaffected and cannot be co-trained with any of
the above.

## Claim boundary

There is no demonstrated excellent/SOTA candidate yet. All current values are
development evidence. Any future claim must report MSE/RMSE, CI, Spearman,
activity-cliff sign, novelty strata, component bootstrap, seed consistency,
resource cost and clean label/protein counterfactuals — and protein
counterfactuals must be **centered** to speak to ordering.

**No result from this architecture may be described as pocket-aware,
contact-resolved, binding-site-localized or biologically localized.** The
protein path is exactly invariant to residue-slot permutation (measured
2.4e-08 pK, proved algebraically): the atom-to-residue cross-attention reads an
unordered bag of sequence-window summaries, and no complex geometry exists for
any of 17,717 DTA cells.

## 2026-08-17 update: Stage C/D/E/F/G cycle results

Leak-free references (meta_train internal-validation checkpoint selection;
frozen meta_val banks, component-weighted, restored pK^2):

- T2 (incumbent recipe retrain, seed 20260815): k=0 2.5961, k=1 1.7712,
  k=2 1.3245, k=3 1.2197, k=5 0.9859; k=5 Spearman 0.314, CI 0.619.
  Multi-seed band (3 seeds): k=0 2.458-2.981, k=5 0.946-1.007.
- G (ESM-650M residue-input lane, 3 seeds): k=0 2.239-2.790, k=5
  0.944-0.987; pooled G-T2 intervals all cross zero - the lane is NOT
  confirmed (REPORT_G2.md).

Every mechanism arm of this cycle (LSP panel level head, T2-LEVEL,
LSP-NOROUTE, F pairwise transport, F-ABS) failed its preregistered gates;
none is promoted. The level/shape anatomy is quantified in
tools/research/stageD_level_panel/D0_REPORT.md and D0b_DOC_TRANSFER.json:
within-document assay history transfers 45% of level variance across targets,
the double-cold split makes that signal unavailable at inference, and the
tested governed probes (panel composition + protein sequence) explain up to
25.9%. No tested candidate reached k=0 MSE <= 1.00 with any representation or
framework family run to date; k=5 sits at the 1.00 boundary in development.
The target remains arithmetically possible: the measured centered term is
0.8648, so an oracle level predictor would put k=0 near 0.865. See
report/BOUNDARY_20260817_NIGHT.md.

## 2026-08-18 update: cycle completed (I through Q) - final state

- Stage I (live ESM-150M LoRA): REJECTED; two resolved ranking observations
  (k=2 Spearman vs its frozen live control, k=3 Pearson vs T2), no resolved
  MSE movement, slightly worse level.
- Stage J (assay-aware level head, journal/panel/protein): REJECTED; best
  k=0 level^2 on record (1.30) but resolved k=2/3 ranking degradation.
- Stage K/K2 (contrastive coembedding): K-REG = first ALL-k resolved MSE
  improvement across three seeds (k=0 -0.112, k=1 -0.048, k=2 -0.027, k=3
  -0.022, k=5 -0.012; ranking preserved; zero control inversions); its k=0
  centered gain did not survive pooling, so NOT CONFIRMED and nothing
  promoted.
- Stage L (support-gated assay-aware level head): REJECTED; best k=0
  calibration in the record (MSE 2.0997, level^2 1.2151) but resolved
  k=2/3/5 ranking degradation.
- Stage M0 (ChemBERTa-77M ligand LM): REJECTED at identifiability (ordering
  r +0.147 below the occupancy record; level probe = grand mean).
- Stage P0 (ProteinKG25 GO annotations, 313/387 matched): REJECTED at
  identifiability (2.27 vs 1.43 constant on the covered subset).
- Stage Q (decoupled frozen-feature level head): REJECTED; Q0 joint frozen
  probe 1.3416 (best frozen predictor on record) but the trained composition
  degrades k=0/2/3 ranking with resolved intervals - the level/ordering
  conflict on one shared trunk was reproduced across four tested compositions
  (E, J, L, Q). That is evidence about those compositions, not a theorem
  about single-stage training in general.
- Final boundary audit (stageN_audit, regenerated 2026-08-18): every
  load-bearing number re-derives bitwise from the raw rows; 106 seal
  artifacts, 0 evaluations; 11 retained trained stages discovered from the
  filesystem, all preregistered (one mtime ordering exception disclosed in
  AUDIT_REPORT.md). External validation: Nelen et al., J Cheminform 17:8, 2025.

No candidate passed all promotion gates; no sealed meta_test label entered
any fitting, selection or reported metric (0 evaluations in the audited
artifacts); nothing moved to model/ or scripts/. Closing summary:
report/FINAL_STATE_20260818.md.

## 2026-08-17 update: Phase 1 protein-conditioned interaction decision

`tools/research/stageV_core_mmp/PHASE1_FINAL_DECISION.md` is the bounded
negative for the tested protein-conditioned interaction mechanisms under the
current BindingDB-Ki double-cold protocol: Stage S global FiLM rejected;
Stage T core-blind pooled discriminator rejected with its closure withdrawn;
Stage U U0 degree-concentration fail; Stage V primary repeated-key surface
not evaluable (32 rows/4 components; 0 exact keys shared with the
development-validation split) and interaction variance not identifiable above
the defensible noise envelope (`theta = -0.406 [-0.704, -0.073]` with the
preregistered noise; direct pair-level noise 0.303 [0.200, 0.427] leaves
cross-component V1 unresolved at +0.391 [-0.327, +0.368]). This is not a
biological impossibility claim; MSA/coevolution is externally blocked, and
Davis/KIBA remain promotion-gated.
