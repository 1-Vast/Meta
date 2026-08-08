# E-AFF-R0 Readout Diagnosis

Status: registered before running the invariance check, the credit simulation or
the published-effect summary. R0 reads no affinity label. It is a diagnosis of
the measuring instrument, not a new affinity experiment and not a Gate.

## Question

Every affinity result in this project — P1C, P1R1, P1R2A/B, E-AFF-P0, H0A and
H0C — was scored with within-task concordance, macro-averaged over closure
components. The frozen theory proves convergence in a different functional,

```text
d_M(F,g) = Hausdorff-W1 between K(B(z)F(z)) and K(B(z)g(z)),
```

and its scope chapter explicitly declines to provide "pairwise, listwise, or
metric ranking" guarantees or "derivation of ranking from affinity regression".

R0 asks a single question:

> What does the registered readout assign credit to, and what is it blind to?

## Parts

1. **Algebraic invariance.** Using the repository's own `metrics.concordance`
   and the real H0C task sizes, check whether within-task concordance changes
   under a per-task shift or positive rescaling of predictions or labels. A
   task-level affinity location and scale are exactly such maps.
2. **Credit simulation.** Generate `y[t,i] = level[t] + within[t,i]` and score
   oracles that see only the level, only the within-task variation, or both,
   under the registered readout and under a location-sensitive error. Sweep the
   share of variance held by the level channel. The level is the channel a
   protein sets for a chemical series.
3. **Published effect.** Summarise what the frozen geometry actually did to
   H0C's per-task scores, to separate "the feature was inert" from "the feature
   moved predictions without a consistent direction".

## Verdicts

- `READOUT_BLIND_TO_TASK_LEVEL_AFFINITY_LOCATION` if part 1 finds exact
  invariance.
- `PERFECT_LEVEL_PREDICTOR_SCORES_CHANCE_AT_EVERY_VARIANCE_SHARE` if part 2
  finds the level oracle at exactly `0.5` for every variance share.
- `READOUT_DIAGNOSIS_INCONCLUSIVE` otherwise.

## Scope Limits

R0 does **not** show that protein-specific affinity lives in the location
channel. It shows only that if it does, the registered readout assigns it zero
credit. Establishing that the channel carries correct-protein information is the
job of a separately registered Gate with its own controls, and R0 authorizes
nothing: not X1, not angular work, not RFSA, not DAVIS, not production.

R0 also does not retroactively convert any past negative result into a positive
one. Those results remain valid statements about within-task ranking.
