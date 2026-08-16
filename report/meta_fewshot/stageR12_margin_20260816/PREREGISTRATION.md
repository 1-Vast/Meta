# Stage R12 preregistration: margin-ranking shape objective on the C2 base

Frozen before any R12 run. Population, bank, seeds, budget, selection and the
wrong-protein contract are identical to R7-R11. `meta_test` remains sealed.

## Hypothesis

The R9 audit measured the residual CI deficit as **margin compression**: C1's
mean absolute predicted margin is 0.097 vs A0's 0.121, and the small-gap /
mid-gap strata remain below A0 while no stratum is resolved. The RankNet
(softplus) shape loss gives vanishing gradients on correctly ordered pairs
with modest margins, so the model stops pushing exactly where label noise
flips the ordering. A **margin-ranking (hinge) loss**
`max(0, m - sign(dy) * dp)` keeps pushing every comparable pair until its
predicted margin exceeds `m`, directly countering the compression. Literature
basis: margin-based ranking in metric learning (Hadsell, Chopra & LeCun 2006;
Weinberger & Saul 2009) and the pairwise-loss/global-epistasis analysis of
Diaz et al. [arXiv:2305.03136] (https://arxiv.org/abs/2305.03136), which
frames pairwise losses over biological fitness/affinity exactly as a
level/shape decomposition — the structure of this project's factorized model.
This is a **loss-form change on the shape objective** (a training-innovation
candidate), not a loss-weight adjustment of the closed gate family: the
transport remains the retained Tanimoto baseline, the architecture is
untouched.

## The single variable

C2's configuration exactly (cliff_pair_weight 2.0, shape_variance 1.5,
relative 1.0, no gate), except the ranking loss form:
`ranknet` -> `margin` with `ranking_margin = 0.1` (normalized pK units).
1200 steps, 3 episodes/step, lr 6e-4 cosine, seeds 20260815/16/17. One
margin value, no sweep.

## Arms

- **A0** frozen R3R4 incumbent checkpoints;
- **C2** the R9 checkpoints (RankNet, cliff weight 2) — the single-variable
  control, already trained;
- **D2** the R12 arm (margin ranking).

## Gates

- **M1** D2 k=0 CI improves over C2's 0.548 and is no more than 0.02 below
  A0's 0.580 (target: >= 0.560);
- **M2** D2 k=0 MSE does not regress beyond C2's 2.119 (point estimate);
- **M3** D2 k=5 activity-cliff sign accuracy stays >= 0.70;
- **M4** D2's k=0 shape term stays <= A0's 0.913;
- **M5** the CI direction holds in all three seeds and the D2-vs-C2 CI
  contrast has a positive component-bootstrap lower bound.

## Failure conditions

M1 failure (no CI improvement) falsifies the margin-compression hypothesis
as the actionable lever: the compression is then a *symptom* of the shape
branch's expressivity, and the next hypothesis targets the shape readout's
capacity (documented separately, not swept). M2/M3/M4 failures record the
tradeoff; the dose is not swept after the fact.

## Resources

Three 1200-step runs, executed by `scripts/run_stage.py` (GPU smoke first);
commands recorded in the stage's `commands.jsonl`.
