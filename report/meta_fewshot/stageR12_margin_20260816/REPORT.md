# Stage R12: margin-ranking shape objective — falsified as the actionable lever

Numerical authority: `RESULT.json` and `COMPARE_R12_meta_val.json` in this
directory. Population: double-cold `meta_val`, 41 targets / 19 components,
three seeds, 1200 steps. `meta_test` sealed and not evaluated.

> **Record note (2026-08-16).** This stage was run and compared on 2026-08-16
> but its narrative and `RESULT.json` were never written; the verdict lived
> only in `task.md` and `EVIDENCE_LEDGER.md`. Both files are backfilled here
> from the retained comparison artifact. Nothing was re-run and no number
> changed.

## The hypothesis

The R9 pair audit measured the residual CI deficit as **margin compression**:
C1's mean absolute predicted margin is 0.097 against A0's 0.121, and no
stratum remained resolved after the cliff-weight removal. The RankNet
(softplus) shape loss has vanishing gradients on correctly ordered pairs with
modest margins, so the model stops pushing exactly where label noise flips
the ordering. A margin-ranking (hinge) loss `max(0, m - sign(dy) * dp)` keeps
pushing every comparable pair until its predicted margin exceeds `m`.

## The single variable

C2's configuration exactly (cliff_pair_weight 2.0, shape_variance 1.5,
relative 1.0, no gate), changing only `ranknet -> margin` with
`ranking_margin = 0.1`. One margin value, no sweep.

## Outcome (k=0, three-seed means)

| arm | MSE | CI | Spearman | calib | shape | cliff sign |
|---|---:|---:|---:|---:|---:|---:|
| A0 incumbent | 2.149 | 0.580 | 0.223 | 1.236 | 0.913 | 0.512 |
| C2 (RankNet, control) | **2.119** | 0.548 | 0.126 | **1.218** | **0.901** | **0.621** |
| D2 (margin) | 2.154 | 0.551 | 0.136 | 1.245 | 0.908 | 0.618 |

At k=5, D2 reaches CI 0.617 and cliff sign 0.742 against C2's 0.612 / 0.775.

## Gate outcomes

| gate | target | measured | outcome |
|---|---|---:|---|
| M1 | k=0 CI >= 0.560 and above C2's 0.548 | 0.551 | **fail** |
| M2 | k=0 MSE no worse than C2's 2.119 | 2.154 | **fail** |
| M3 | k=5 cliff sign >= 0.70 | 0.742 | pass |
| M4 | k=0 shape <= A0's 0.913 | 0.908 | pass |
| M5 | 3/3 seed CI direction + positive D2-vs-C2 component bootstrap | — | **not evaluable** |

**M5 is recorded as not evaluable, not as a pass.** The retained comparison
artifact carries a `D2_vs_A0` contrast but no `D2_vs_C2` contrast, so the
bootstrap the preregistration named against its stated control was never
computed. For what it is worth, the `D2_vs_A0` k=0 CI interval is +0.029
[-0.034, +0.099] — it crosses zero.

## Verdict

M1 fails: the margin loss moved k=0 CI by **+0.003** against its control
while k=0 MSE regressed. Under the preregistered failure condition, **margin
compression is a symptom of the shape branch's expressivity, not the loss
form.** The loss-form axis is closed. The one resolved effect in the
comparison is negative: D2's k=0 cliff sign is -0.107 [-0.201, -0.011]
against A0, so the hinge loss also costs cliff ordering relative to the
incumbent even while beating it on the same metric in absolute terms
(0.618 vs 0.512) — the interval is against A0 on 14 components, and it is the
*only* resolved cell in this stage.

## Where this leaves the ladder

R12 eliminates the loss form. Together with R10 (variance term, falsified)
and R11 (trunk routing, falsified), the remaining hypothesis the R9-R13 chain
points at is the shape readout's **expressivity and inductive bias** — which
R13 tested directly and which failed its Stage 1 gates. The consolidated
statement is `report/BOUNDARY_20260816.md`.
