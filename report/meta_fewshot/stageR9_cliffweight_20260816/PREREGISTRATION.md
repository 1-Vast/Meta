# Stage R9 preregistration: activity-cliff pair-weight dose response

Frozen before any R9 run. Population, bank, seeds, budget, selection and the
wrong-protein contract are identical to R7/R8. `meta_test` remains sealed.

## Measured premise (R8 artifacts, before any R9 run)

The stronger-shape arm B1 trades overall ordering for cliff ordering at
k=0: CI 0.535 vs A0's 0.580 while cliff sign accuracy rises 0.512 -> 0.598
(+0.086); at k=5, CI -0.024 for cliff +0.093. The shape objective's
cliff-pair weighting (x4) is the prime suspect for this distortion: it
concentrates the ranking gradient on the ~5-10% of pairs with
Tanimoto >= 0.6 and gap >= 1.0 pK, at the expense of the small-gap pairs
that decide the overall concordance index.

## Hypothesis (single variable)

Reducing the activity-cliff pair weight recovers overall CI while retaining
most of the shape gain. Everything else is exactly B1's configuration:
`--no-gate --shape-variance-weight 1.5 --relative-loss-weight 1.0`,
1200 steps, 3 episodes/step, lr 6e-4 cosine, seeds 20260815/16/17.

## Arms

- **A0** frozen R3R4 incumbent checkpoints (unchanged);
- **B1** the R8 checkpoints (cliff weight 4.0) — the dose-response control,
  already trained;
- **C1** cliff weight 1.0 (no cliff emphasis), 3 seeds, 1200 steps;
- **C2** cliff weight 2.0 (halved), 3 seeds, 1200 steps.

## Gates (identical in form to R8)

- **Z1'** arm k=0 MSE point estimate at least -2% vs A0;
- **Z5'** arm k=0 CI no more than 0.02 below A0's;
- **S-shape** arm's k=0 shape term below A3's 0.905;
- **S-corr** arm beats A3 (2.197) on k=0 MSE.

## Decision rule

If C1 or C2 passes Z1' and Z5' simultaneously, that dose is carried forward
to a three-seed re-evaluation of the standing Z1-Z7 gates. If neither passes
both, the cliff-weight axis is recorded as falsified and the next
single-variable hypothesis is the ranking pair-weighting scheme itself
(gap-weighted / LambdaRank-style), preregistered separately. No gate moves
after the fact.
