# Stage R8: stronger shape signal — family closed for the double-cold zero-shot target

Numerical authority: `COMPARE_R8_meta_val.json` (+`.rows.jsonl`), per-arm
`RESULT.json`. Population and contract identical to R7; gates frozen in
`PREREGISTRATION.md` before any run. A0 = frozen R3R4 incumbent checkpoints;
A3 = R7's best arm (routed + counterfactual + relative + ranking, no gate);
B1 = A3's configuration with `shape_variance_weight 1.5` and
`relative_loss_weight 1.0` — one changed variable. Three seeds, matched
1200-step budget. `meta_test` not opened.

## Result (3-seed mean)

| arm | k=0 MSE | CI | Spearman | calib | shape | k=5 cliff sign |
|---|---:|---:|---:|---:|---:|---:|
| A0 incumbent | 2.149 | **0.580** | **0.223** | **1.236** | 0.913 | 0.675 |
| A3 (R7) | 2.197 | 0.551 | 0.136 | 1.292 | 0.905 | 0.676 |
| B1 stronger shape | **2.167** | 0.535 | 0.096 | 1.271 | **0.896** | **0.768** |

B1 vs A0, paired component bootstrap: k=0 MSE -0.018 [-0.243, +0.229]
(unresolved tie); k=0 CI -0.045 (regression); B1's shape 0.896 is the best
shape term recorded anywhere in this project, and its k=5 activity-cliff
sign accuracy 0.768 is the best cliff-ordering result recorded (A0: 0.675).

## Gate outcomes

- **Z1' fail**: -0.8% against the preregistered -2% threshold (and far from
  the standing Z1 target of -10%);
- **Z5' fail**: CI regresses by 0.045 against a 0.02 tolerance;
- **S-shape pass**: shape 0.905 -> 0.896 (the stronger shape signal moves
  shape);
- **S-corr pass**: B1 beats A3 at k=0 (2.167 vs 2.197) — the shape gain is
  not paid for by worse calibration.

## Verdict

Under the preregistered decision rule (Z1' and Z5' both required to
advance), **the model family is closed for the double-cold zero-shot target
as a claimed core innovation.** `meta_test` remains sealed and unopened, and
is explicitly not opened on the strength of these results.

The positive findings are retained, not discarded: the shape-first training
method is a real, measured performance source for within-target ordering
(the first in this project: shape 0.943 -> 0.896 across the ordinary-to-full
ladder, and cliff sign 0.518 -> 0.768 at k=5); the routed level readout with
attention pooling converges to the incumbent's calibration at the full
budget (1.271 vs 1.236); the counterfactual and identifiability machinery is
verified at the gradient level. The binding constraint remains the
level/shape tradeoff at k=0: every configuration that improves shape pays
more in calibration (or CI) than it recovers in the zero-shot MSE.

## Next cycle (not preregistered here; recorded as the handoff)

The single open question for a future cycle: how to retain B1's shape gain
(0.896, cliff 0.768) without the CI regression (0.535 vs 0.580) and with a
calibration better than A0's 1.236. Candidates, in order of the evidence:
(1) a calibration head that reads the support set at k>=1 but leaves k=0
untouched; (2) listwise (LambdaRank-style) weighting of the pairwise shape
loss instead of flat RankNet weighting; (3) budget scaling beyond 1200
steps, which every full-budget trajectory in this project was still
improving at when stopped.
