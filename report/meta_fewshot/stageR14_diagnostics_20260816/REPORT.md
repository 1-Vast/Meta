# Stage R14 diagnostics: the within-target ordering coefficient is the only lever

No training. Frozen checkpoints, the identical fixed `meta_val` bank, forward
passes only. `meta_test` never read.

Numerical authority: `DISPERSION_meta_val.json`,
`DISPERSION_attribution_meta_val.json`, `DISPERSION_rows.json` (+ per-target
`.rows.jsonl`), produced by `scripts/r14_dispersion_audit.py`. The audit
reproduces the retained authority exactly — A0 k=0 MSE 2.1488, calibration
1.2358, shape 0.9130 — under the project's equal-component, equal-target
weighting.

## The decomposition

For centered predictions `p` and centered labels `y` inside one target, with
`r = corr(p, y)`, the shape term splits **exactly**:

    shape = E[(p - y)^2] = Var(y)(1 - r^2)  +  (sd_p - r·sd_y)^2
                           \___ordering___/    \____amplitude____/

The two halves behave completely differently:

* the **ordering floor** `Var(y)(1-r²)` is untouchable by any rescaling. It
  falls only if `r` rises;
* the **amplitude excess** vanishes at `sd_p = r·sd_y`. Rescaling a centered
  prediction by a positive scalar is monotone within a target, so it moves
  MSE and **cannot change CI, Spearman or sign accuracy at all**.

This matters because the R7-R13 ladder tracked CI and cliff-sign accuracy,
which are rank statistics dominated by easy large-gap pairs, and never
tracked `r`, which is the quantity that enters the MSE.

## Result 1 — every ranking-primary arm has worse k=0 ordering than plain MSE

Double-cold `meta_val`, k=0, three seeds, equal-component weighting:

| arm | training | MSE | calib | shape | **ordering floor** | amp. excess | **r** | sd_p/sd_y | CI |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| **A0** incumbent | MSE-primary | 2.149 | 1.236 | 0.913 | **0.692** | 0.221 | **0.213** | 0.180 | 0.580 |
| D1 (R10) | shape-first, var 0.5 | 2.285 | 1.358 | 0.927 | 0.706 | 0.221 | 0.156 | 0.023 | 0.552 |
| C2 (R9) | shape-first, cliff 2 | 2.119 | 1.218 | 0.901 | 0.716 | 0.185 | 0.148 | 0.142 | 0.548 |
| D2 (R12) | shape-first, margin | 2.154 | 1.245 | 0.908 | 0.725 | 0.183 | 0.129 | 0.084 | 0.551 |
| C1 (R9) | shape-first, cliff 1 | 2.235 | 1.332 | 0.903 | 0.731 | 0.172 | 0.150 | 0.132 | 0.562 |
| B1 (R8) | shape-first, var 1.5 | 2.167 | 1.271 | 0.896 | 0.735 | 0.161 | 0.096 | 0.123 | 0.535 |
| G1 (R11) | shape-first, **same trunk** | 2.405 | 1.488 | 0.917 | 0.751 | 0.166 | 0.134 | 0.204 | 0.525 |
| B3 (R3R4) | full routed method | 2.055 | 1.131 | 0.924 | 0.753 | 0.172 | 0.073 | 0.225 | 0.531 |

**The plain MSE-trained incumbent has the lowest irreducible ordering error of
every arm in the project — 8 of 8, without exception.** The ordering floor is
monotone in the opposite direction to everything the ladder was optimising.

### The attribution is clean, because G1 changes nothing but the training

G1 is the incumbent architecture trained with the shape-first routed method —
**zero architecture change**, by construction (R11). At k=0 it takes `r` from
**0.213 to 0.134**, a 37% loss of within-target ordering, on identical
weights-shape. The damage is a property of the training method, not of the
factorized trunks.

The damage is also specific to k=0. At k=5 the ordering floors converge and
G1 slightly *beats* A0 (0.631 vs 0.646, r 0.357 vs 0.348) — support
observations supply ordering that the zero-shot readout cannot.

## Result 2 — a retained claim needs correcting

The record states that shape-first training is "the project's first measured
within-target shape source", citing B1's shape 0.896 against A0's 0.913. The
decomposition shows that reading is **causally wrong**:

| | ordering floor | amplitude excess | shape |
|---|---:|---:|---:|
| A0 | 0.692 | 0.221 | 0.913 |
| B1 | 0.735 | 0.161 | 0.896 |
| difference | **+0.043 worse** | −0.060 | −0.017 |

B1's shape improvement is **shrinkage toward the target mean**, and it is
bought on top of *worse* ordering. The same holds for every shape-first arm.
The honest statement is that shape-first training is a measured **amplitude
suppression** source, not an ordering source, at k=0.

This does not retract the activity-cliff finding — C1's k=5 cliff sign 0.782
is a real measurement — but it explains it: cliff pairs are high-similarity,
large-gap pairs, and ordering them correctly is compatible with poor global
within-target correlation. C1 holds the cliff record with `r` = 0.150 against
A0's 0.213.

## Result 3 — amplitude is not a usable lever, tested two independent ways

**Test A, honest global rescale.** One scalar fitted leave-one-component-out
on 18 components, applied to the held-out 19th:

| arm | k=0 MSE | after LOCO rescale |
|---|---:|---:|
| A0 | 2.1488 | 2.1570 |
| B1 | 2.1666 | 2.1668 |
| B3 | 2.0554 | 2.0581 |
| C1 | 2.2348 | 2.2374 |
| C2 | 2.1195 | 2.1263 |
| D1 | 2.2846 | 2.2613 |
| D2 | 2.1537 | 2.1493 |
| G1 | 2.4053 | 2.4096 |

Six of eight get **worse**; the two that improve are the two most collapsed
arms (D1 at sd ratio 0.023, D2 at 0.084) and the gains are ≤0.023. A single
global amplitude is already at its optimum.

**Test B, is the per-target optimum a stable property of the target?** The
per-target *oracle* rescale would take A0 from 2.1488 to **1.9280** — 0.22
of MSE, which would by itself clear the preregistered Z1 target of 1.934.
That gap is only reachable if the optimal scale is a property of the target
rather than a fit to one panel's label noise. Transferring the scale across
seeds of the *same* target (identical target, different trained model)
recovers 2.0053, about 65% of the gap. But:

* the per-target optimal scale has median 2.00 with **IQR [0.02, 5.36]**;
* it is **negative in 25.2% of targets** — for a quarter of the population the
  model's within-target ordering is *anti-correlated* with truth, so the
  MSE-optimal action is to flip the sign. A mechanism that learns which
  targets the model gets backwards is an ordering fix, not an amplitude one;
* the estimate is numerically unstable wherever `sd_p` is small, which is
  every ranking-trained arm — C2's seed-transfer diverges to 29.8.

**A protein-conditioned amplitude head is therefore not selected.** It was the
model-side candidate in `report/LITERATURE_R14_20260816.md`; it is falsified
here, before implementation, by two independent tests.

## What this localizes

At the MSE-optimal amplitude the shape term *equals* the ordering floor, so

    achievable k=0 MSE = calibration + Var(y)·(1 - r²),   Var(y) ≈ 0.725

`r` is the only quantity that improves MSE and CI **together**. Everything the
R7-R13 ladder varied — cliff weight, variance weight, loss form, routing,
readout parameterization — moved amplitude and coarse rank statistics while
leaving `r` at or below the incumbent's 0.213.

Targets, with A0's calibration held:

| k=0 `r` | ordering floor | k=0 MSE | vs A0 |
|---:|---:|---:|---:|
| 0.213 (A0 today, at its own optimal amplitude) | 0.692 | 1.928 | −10.3% |
| 0.30 | 0.660 | 1.896 | −11.8% |
| 0.40 | 0.609 | 1.845 | −14.1% |
| 0.50 | 0.544 | 1.780 | −17.2% |

The first row is the important one: **A0's existing ordering, if its amplitude
were per-target optimal, already clears Z1.** The blocker is that the optimal
amplitude is not predictable from the protein (Result 3). So the reachable
route is to raise `r` far enough that a *near-constant* amplitude is close to
optimal — which happens automatically, because the optimal scale `r·sd_y/sd_p`
becomes both larger and less sensitive as `r` grows and `sd_p` grows with it.

## Result 4 — the alignment identity holds, and it names the mechanism

`scripts/r14_alignment_check.py`, 200 synthetic panels of 16 queries with
realistic pK statistics. The quantity is the gradient the *ranking* term
contributes **at the regression optimum** `s = y` — that is, at a prediction
that is already exactly right:

| within-target ranking term | max ‖∂L/∂s‖ at `s = y` |
|---|---:|
| regression-compatible ListCE | **1.7e-17** (machine zero) |
| RankNet softplus (R8/R9/R10 shape loss) | 1.02e-01 |
| hinge margin 0.1 (R12 shape loss) | 7.36e-02 |

**RankNet and hinge keep pushing even when the prediction is perfect.** Their
optima are scale-free monotone arrangements, so they always want wider
margins; that residual gradient is a force pulling the model away from
`s = y`, and Phase 2 measured what it costs — `r` falls in 8 of 8 arms.

The regression-compatible ListCE is exactly stationary at `s = y`, and also
at any positive multiple of the centered label (3.9e-18 at scale 1.7), which
is the property that lets the MSE term pin the amplitude while the listwise
term supplies ordering pressure *only when the ordering is wrong*.

The derivation, for labels `y`, a fixed shift `m` below the label range,
`T(x) = x - m`, and weights `w = y - m`:

    ∂ListCE/∂s_k ∝ w_k/T(s_k) - Σ_j w_j / Σ_j T(s_j)

vanishes for every `k` exactly when `T(s) ∝ w`. This was the design's first
listed failure mode ("the alignment identity does not hold for squared
error"); it is cleared.

## Consequence for the R14 design

One core innovation survives Phase 2, and it is a training innovation: an
objective whose ranking term **cannot** trade `r` away for coarse rank
statistics. The measured premise is Result 1 (8/8 arms, with a zero-
architecture-change attribution), and the mechanism is the regression-
compatible construction surveyed in `report/LITERATURE_R14_20260816.md`.

The model-side candidate is **dropped**, not replaced. Phase 2 falsified it,
and substituting an unmeasured alternative to keep a two-innovation shape
would be exactly the auxiliary decoration the mandate forbids.

`meta_test` remains sealed and unopened.
