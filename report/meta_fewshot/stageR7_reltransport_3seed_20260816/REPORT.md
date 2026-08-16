# Stage 3 (R7): three-seed formal development — admission refused

Numerical authority: `COMPARE_R7_meta_val.json` (+`.rows.jsonl`), per-arm
`RESULT.json`, gates frozen in `PREREGISTRATION.md` before the runs.
Population: governed double-cold `meta_val` (41 targets, 19 components),
`evaluation_seed=73101`, `query_size=16`, one draw, all eligible targets.
A0 = the frozen Stage R3/R4 incumbent checkpoints (seeds 20260815/16/17,
1200 steps, `similarity_only`). A1/A2/A3 = relative-transport trunk,
3 seeds each, matched 1200-step budget. Seeds are averaged inside a target
before components are resampled, so every interval is conditional on the
trained seeds. **`meta_test` was not opened.**

## Headline table (3-seed mean)

| arm | k=0 MSE | CI | Spearman | cliff sign | calib | shape | k=5 MSE |
|---|---:|---:|---:|---:|---:|---:|---:|
| A0 incumbent | **2.149** | **0.580** | **0.223** | 0.512 | 1.236 | 0.913 | 0.915 |
| A1 ordinary | 2.519 | 0.520 | 0.059 | 0.518 | 1.576 | 0.943 | 0.958 |
| A2 full method | 2.420 | 0.542 | 0.114 | **0.536** | 1.525 | **0.895** | **0.907** |
| A3 full, no rho gate | 2.197 | 0.551 | 0.136 | 0.533 | **1.292** | 0.905 | 1.092 |

Paired component bootstrap, A2 minus A0 (negative = A2 worse): k=0
**-0.271 [-0.683, +0.091]**, k=1 **-0.399 [-0.745, -0.103]** (resolved),
k=2 **-0.247 [-0.434, -0.097]** (resolved), k=3 -0.081 [-0.199, +0.022],
k=5 +0.009 [-0.045, +0.061].

## Gate outcomes

Zero-shot (k=0): **Z1 fail** (-12.6% against a +10% target); **Z2 fail**
(lower bound -0.683); **Z3 fail** (seeds 2.130 / 2.399 / 2.730 — mixed
direction); **Z4 fail** (Tanimoto<0.4 tier 2.373 vs 2.103); **Z5 fail** (CI
-0.038 against A0, Spearman 0.114 vs 0.223); **Z6 fail at k=0** (wrong-protein
2.312 is *below* full 2.420); **Z7 pass** (interaction-cut `ligand_only` is
clearly worse than full: seed 1 shows 3.266 vs 2.130 — the branch is a real
performance source at the full budget).

k=1: **F1 pass** (the correction is query-specific by construction; the
Stage 1 gate suite proves the dependence on the query ligand); **F2 fail**
(level-only beats full by -0.458, point estimate); **F3 pass** (correct
support beats the magnitude-matched wrong label: permuted gap +0.65);
**F4 pass** (support-ligand replacement moves the prediction);
**F5 fail** (MSE improves over k=0 but CI does not: 0.542 at both).

k>=2: **G1 fail** at k=2,3 (level-only beats full by -0.358/-0.225; +0.02 at
k=5); **G2 fail** (the trained rho is eval-inert: `nogate` gap 0.004..0.000,
so "full" is the Tanimoto baseline plus nothing); **G3 pass** (permuted gaps
+0.70/+0.24/+0.18/+0.38 — correct support beats permuted labels); **G4 fail**
(CI/Spearman worse than A0 at k=0-3); **G5 inconclusive** (k=5 gains exist
but are not separated from the zero-shot).

Innovation B: **T1 pass in point estimate** (A2 beats A1 by +0.099 at k=0 —
the training method is a real source for the new architecture); **T2 pass**
(shape 0.943 -> 0.895 under the method — the first real shape improvement
recorded in this project); **T3 fail** (the total A2-vs-A0 change is
negative, so no share of a positive total can be attributed); **T4 not
evaluable** given T3.

## What was learned (the mechanism, measured)

1. **The routed level readout converges to the incumbent's calibration at
   the full budget.** A3: calibration 1.292 vs A0's 1.236, at 1200 steps —
   the same convergence the R3R4 ladder needed. The 300-step screenings
   could not see this; recorded as the reason the formal budget is the
   decisive test.
2. **The shape-first training produced the project's first real shape
   improvement.** A2's shape 0.895 vs A1's 0.943 (same architecture), and
   A2's activity-cliff sign accuracy 0.536 beats A0's 0.512. The
   relative-supervision + cliff-weighted ranking signal works.
3. **The query-specific transport channel costs more than it returns —
   again.** The rho gate is eval-inert (nogate gap ~0) while its training
   gradients disturb calibration: A2 calib 1.525 vs A3's 1.292 at the same
   budget. This is the seventh query-specific channel in this project to
   show the same signature (inert at deployment, costly in training), now
   with the gate trained under ranking-primary shape objectives — so the
   objective is not the explanation this time.
4. **The zero-shot target remains unmet by every design tested.** The best
   new arm (A3) is 2.2% behind A0 at k=0; the best-ever double-cold k=0
   (B3's 2.0554) is 4.3% ahead. The 10% Z1 target has never been approached.

## Verdict

**Admission refused.** Z1, Z2, Z3, Z4, Z5 and Z6 fail at k=0; F2 and G1/G2/G4
fail. `meta_test` remains sealed and unopened. All artifacts (checkpoints,
per-target predictions, gradient diagnostics, coverage censuses, resource
records) are retained in this directory.

## Next single-variable hypothesis (for Stage R8, preregistered before any run)

The measured facts select the hypothesis: A3's training (routed +
counterfactual + relative + ranking, **without** the rho gate) reaches
near-incumbent k=0 calibration while A2's shape gain shows the shape
objectives work. The next experiment therefore holds A3's configuration
exactly and raises the shape signal only — `shape_variance_weight 1.0 -> 1.5`
and `relative_loss_weight 0.5 -> 1.0` — with the transport left as the
retained Tanimoto+key baseline (no query-specific gate). Success criterion
unchanged (Z1-Z7); if the stronger shape signal does not move k=0 MSE below
1.934 across three seeds, the double-cold k=0 target is recorded as unmet by
this model family and the family is closed.
