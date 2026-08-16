# Stage R8 preregistration: stronger shape signal on the A3 configuration

Frozen before any R8 run. Population, bank, seeds, budget, selection and the
wrong-protein contract are identical to R7. `meta_test` remains sealed.

## Hypothesis

R7 measured two facts: A3's configuration (routed level + counterfactual +
relative supervision + ranking, Tanimoto transport, no query-specific gate)
reaches near-incumbent calibration (1.292 vs 1.236) and k=0 2.197 (-2.2% vs
A0); and the shape objectives produce a real shape gain (A2 shape 0.895 vs
A1's 0.943). The single variable changed here is the **shape signal
strength**: `shape_variance_weight 1.0 -> 1.5` and
`relative_loss_weight 0.5 -> 1.0`, everything else identical to A3.

## Arms

- **A0** frozen R3R4 incumbent checkpoints (unchanged);
- **A3** the R7 checkpoints (shape 1.0 / relative 0.5, no gate) — the
  single-variable control, already trained;
- **B1** the R8 arm (shape 1.5 / relative 1.0, no gate), 3 seeds,
  1200 steps, 3 episodes/step, lr 6e-4 cosine.

## Gates (unchanged from R7 where applicable)

- **Z1'** B1 k=0 MSE vs A0: point estimate at least -2%, and
- **Z5'** B1 k=0 CI does not regress below A0 by more than 0.02;
- **S-shape** B1's k=0 shape term is below A3's 0.905 (the stronger shape
  signal must move shape);
- **S-corr** B1 beats A3 on k=0 MSE (the shape gain must not be paid for by
  worse calibration).

## Decision rule

If Z1' and Z5' both pass, the design advances to a full re-evaluation of
Z1-Z7 (the Z1 target of -10% is then re-tested at three seeds). If Z1'
fails (B1 k=0 above 0.98 * A0), the model family is **closed for the
double-cold zero-shot target** and recorded as such; meta_test stays sealed.
No gate moves after the fact.
