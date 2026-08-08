# E-AFF-L0 Protein Affinity Level Gate

Status: registered before any arm is fitted or scored. L0 is the first Gate in
this project to measure the functional the frozen theory actually controls. It
reads governed ChEMBL37 source Ki/Kd labels, which the registered E-AFF stages
are permitted to read. DAVIS and recipient labels remain prohibited.

## Question

Does the correct protein, entering only through validated pair-local geometry,
improve prediction of the **affinity level** of a chemical series beyond
ligand-only information, beyond protein sequence-family information, and beyond
a deranged protein?

This is Claim A in `THEORY_BIOLOGY_INTEGRATION.md`. It is not the crossed
interaction claim, which remains gated by X0-B.

## Correct Characterisation Of The Theory Metric

The permitted statement, used throughout this package:

> `W1(P, P + c) = |c|`, so `W1` is sensitive to common translations, while it may
> also respond to distribution-shape changes. The frozen theory controls
> law-class distances and declared Lipschitz losses but does not automatically
> derive pairwise or listwise affinity ranking.

`W1` must not be described as only a location metric. Earlier wording in this
package that called it a location metric is superseded by the statement above.

## Sign Convention Is Not Information

Ordered anchors impose a **deployment sign convention**. They do not prove that
biology contains protein-specific affinity information:

```text
sign fixed by deployment convention
    !=
protein-specific affinity information identified from biological input
```

A model built on ordered anchors that carries no protein information will place
weight without regard to the protein and will fail the Gate. The convention
removes an estimation problem; it manufactures no signal.

## Pre-Execution Amendment A1

Recorded before any arm was fitted or any outcome viewed.

The first draft of this preregistration set `m = 4`, citing the frozen few-shot
limit `d_adapt <= k`. That limit is **not engaged by L0**: every arm maps
features to `p` through a globally fitted, closure-cross-fitted map, and no arm
performs support-based adaptation, so there is no `d_adapt` to bound. Truncating
the frozen anchor lattice to its first four entries would also have removed the
top of the ladder (centres `0.75` and `0.90`), compressing the achievable
affinity range for every arm equally and blunting the Gate's sensitivity.

L0 therefore uses the deployment's own frozen `MetaSieveConfig.m = 7`: six
stochastically ordered logistic anchors at centres
`0.15/0.30/0.45/0.60/0.75/0.90` plus one broad uniform anchor. The
`d_adapt <= k` limit remains binding for any future few-shot adapter, which
remains frozen.

## Operator And Anchor Contract

Frozen by `audit_l0_operator_contract.py`, verdict
`L0_OPERATOR_AND_ANCHOR_CONTRACT_FROZEN`, artifact
`artifacts/eaff_l0_contract_v1/`:

- all seven anchors lie in the frozen band polytope;
- the six logistic anchors satisfy the declared stochastic-dominance order with
  maximum gap `0.0`, and their band mean intervals are strictly increasing from
  `[0.186, 0.303]` to `[0.845, 0.889]`;
- moving simplex weight up the ladder never lowers either endpoint of the
  emitted mean interval, maximum violation `0.0`;
- simplex assembly stays inside the polytope;
- `K` law classes are valid on the fixed mesh and Hausdorff-`W1` respects the
  frozen stability bound, maximum violation `0.0`;
- the output mesh is unchanged, `M = 32`, `h = 1/32`;
- all `258` frozen theory files hash-match, core theory intact.

The broad uniform anchor is deliberately **outside** the dominance order. It is
a width and abstention channel, and no monotonicity is claimed for mixtures that
move weight into or out of it.

## Why This Is A New Question

E-AFF-R0 established that within-task concordance, the readout behind every
previous affinity result, is exactly invariant to per-task affinity location and
scale, and that a perfect level predictor scores `0.5000` under it at every
variance share. The level channel has therefore never been measured. L0 measures
it, in the operator's own metric.

## Output Object

Every arm emits a band through the frozen operator and nothing else:

```text
A(F,z) = K(B(z)F(z)),   B(z) = [ beta_0(z) | beta_1 | ... | beta_m ],   p in Delta_m
```

- `kappa(z) = (endpoint family, assay format stratum)`; `beta_0(z)` is that
  stratum's population band, fitted on training data only.
- `beta_1 <= ... <= beta_m` are fixed anchor bands, stochastically ordered on
  `V`, calibrated once on training data, never per target, and frozen before any
  evaluation arm is scored.
- `m = 4`, so `d_adapt = 4 <= k = 5` and the few-shot dimensional limit holds
  with one residual degree of freedom.

Arms differ **only** in what drives `p`. No arm may use a target identifier, an
assay identifier, or a task identifier.

## Arms

1. `population` — `p_0 = 1`. Prior band for the assay stratum only.
2. `ligand_only` — `p` from ligand state alone.
3. `sequence_only` — `p` from the protein ESM state alone, no geometry.
4. `correct` — `p` from `z_bio`, the pair-local geometry statistic of section
   4.3 of the integration design, built from the frozen P1B bridge and the
   T-BASIS-R0 radial basis with the correct protein.
5. `deranged` — identical to arm 4 with a score-blind wrong protein.

Arm 3 is mandatory. Without it a family-level shortcut ("kinases read high")
passes as a mechanism claim.

## Nuisance Rule

No support-fitted nuisance may be shared between arms. Either every arm fits its
own nuisance from its own inputs, or no arm fits one. This repairs the H0C
defect where `local_score`, fitted on the correct protein's own support labels,
was added to both the correct and the deranged arm, making the derangement
contrast uninformative.

## Splits

- Evaluation proteins lie in closure components disjoint from every component
  used for training, anchor calibration or population-band fitting. The level
  must be predicted, not memorized.
- Ligand scaffold disjointness between fit and evaluation within each task, as
  in H0C.
- No panel consumed by P0, H0A, H0C, T-DIR-P0 or T-BASIS-R0 may be reused as
  untouched validation.
- Assay stratification enters through `kappa` so a protocol offset cannot be
  read as a protein effect.

## Metrics

Primary: the mesh-discretized **interval score** of the emitted band against the
observed `Y`, which is convex in the band, bounded on compact `V`, Lipschitz in
band sup norm — the loss class the frozen theory requires — and which penalizes
width, so abstaining everywhere cannot win.

Secondary: `W1` between the emitted law class and the observed task law;
tertiary: MAE in log units, reported for interpretability only.

Aggregation: macro over closure components, with a component bootstrap for
95% intervals, matching prior stages.

## Frozen Arms

Named `A0..A4` and fixed here. Every arm emits a law class through the same
frozen operator, uses identical outer folds, and is capacity matched: the same
conditional-CDF estimator, the same bandwidth rule, and the **same seven-
dimensional bounded input**, so no arm wins on feature count.

| Arm | Input | Coefficient |
|---|---|---|
| `A0` | none | `p = e_0`, all mass on the stratum population band |
| `A1` | 7 bounded unsupervised components of the ligand state | frozen target |
| `A2` | 7 bounded unsupervised components of the protein sequence state | frozen target |
| `A3` | the 7 `z_bio` coordinates from the correct protein | frozen target |
| `A4` | the 7 `z_bio` coordinates from a score-blind wrong protein | **A3's estimator**, evaluation only |

`A1` and `A2` receive an unsupervised PCA of their native state to seven
dimensions, min-max bounded on the training pool. The projection sees no label
and no evaluation fold.

`A4` is evaluation-only. No arm is ever trained on wrong proteins as negative
affinity examples.

The coefficient is the frozen positive-ridge target
`g_mu^*(z) = argmin_p [L_0(z, B(z)p) + (mu/2)||p||^2]`, computed exactly by
`target_from_conditional_cdf`, from a conditional CDF estimated by the
repository's Nadaraya-Watson `KernelConditional` on the training pool only.
Bandwidth is fixed by the declared rule `n_pool ** (-1/(d+4))` with `d = 7`,
identical for every arm, and is not tuned.

## Frozen Metrics

Primary theory-compatible loss: the frozen `band_loss`, coverage violation plus
`lambda_w` times mean width, which is convex in the band, bounded on compact
`V`, Lipschitz in band sup norm, and penalizes width so a vacuous wide band
cannot win.

Also reported: empirical coverage, mean interval width, location error,
calibration by endpoint, closure-component macro estimates, closure-component
bootstrap 95% intervals, and every arm separately. Ki and Kd are never pooled.

Location is the midpoint of the band-induced mean interval
`[a_max - integral U, a_max - integral L]`, reported in log-affinity units
through the inverse of the frozen training-range affine scaling.

## Frozen sigma_assay Estimator

- **Estimator.** Pooled within-cell standard deviation over `(task, ligand)`
  cells holding at least two governed activities. A task key already fixes one
  exact assay, so within-cell replicates are within-assay replicates. Cells are
  weighted by `n_c - 1`.
- **Missing replicates.** Cells with fewer than two measurements are excluded
  and are never treated as zero noise.
- **Precondition.** At least `100` replicate cells, else
  `L0_NOT_RUN_NUMERICAL_PRECONDITION_FAILED`.
- **Confidence.** A chi-square 95% interval on the pooled SD is reported. The
  margin uses the **point estimate**, so a wide interval can neither soften nor
  harden the Gate.
- **Aggregation.** Estimated per admitted endpoint, on the fitting folds only,
  and reported before any arm is scored.

```text
margin_L0 = 0.5 * sigma_assay
```

the same `0.5` design ratio the project already froze for X0.

## Gate

For each control `C` in `{A1, A2, A4}`, `A3` passes against `C` only if **all
three** hold:

1. `band_loss(C) - band_loss(A3) > 0` with a 95% closure-component-bootstrap
   lower bound above zero;
2. `location_error(C) - location_error(A3) >= margin_L0` in log-affinity units,
   with a 95% closure-component-bootstrap lower bound above zero;
3. `coverage(A3) >= coverage(C) - 0.05`, so `A3` cannot buy an advantage by
   under-covering or by manipulating interval width.

`A3` passes overall only if it passes against `A1`, `A2` **and** `A4`
simultaneously. Beating `A4` alone is a wrong-protein-penalty result and is
insufficient. Beating `A1` alone is insufficient because it may be a
protein-family shortcut. Beating `A2` is required to show value beyond
sequence-level protein information.

Also reported and required non-degenerate: the fraction of queries answered by
abstention, that is `p_0` mass on the population band. An arm that abstains on
most queries has not identified anything, whatever its score.

## Ki/Kd Direction Consistency

If both endpoints were admitted, a pass may not rest on an unexplained direction
reversal between them. The Stage 4 identifiability check admitted **Ki only**
(`C1 = 79`, `C2 = 218`, `C3 = 0.630`) and excluded **Kd** (`C1 = 10 < 30`), so
L0 runs on Ki alone and reports Kd as not identified. No cross-endpoint
aggregation is performed.

## Verdicts

- `L0_PROTEIN_LEVEL_IDENTIFIED` — all three contrasts pass with positive lower
  bounds and non-degenerate identifiability diagnostics.
- `L0_FAMILY_SHORTCUT_ONLY` — arm 4 beats 2 and 5 but not 3.
- `L0_WRONG_PROTEIN_PENALTY_ONLY` — arm 4 beats 5 but not 2 or 3.
- `L0_PROTEIN_LEVEL_NOT_IDENTIFIED` — otherwise.

## Scope Limits

A pass admits `z_bio` as a **candidate** for `z` and authorizes a sealed
transfer Gate. It does not by itself admit anything to `model/`, does not
establish Claim B, does not authorize angular or many-body bases, RFSA, DAVIS,
production integration or P2-P4, and does not license describing any coordinate
as free energy. Ligand-size normalization remains a declared gauge choice.

A failure is informative and terminal for this parameterization of `z_bio`. It
would mean the validated geometry does not carry a transferable protein-specific
affinity level, which — combined with the existing within-task ranking negatives
— would be the strongest evidence yet that the P1B representation is
structurally correct but energetically uninformative.
