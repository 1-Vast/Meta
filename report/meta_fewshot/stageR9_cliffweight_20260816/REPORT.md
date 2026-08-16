# Stage R9: pair-level diagnosis and the activity-cliff pair-weight dose response

Numerical authority: `COMPARE_R9_meta_val.json` (+`.rows.jsonl`),
`PAIR_AUDIT_meta_val.json` (A0/B1), `PAIR_AUDIT_R9_meta_val.json`
(A0/B1/C1/C2). Gates fixed in `PREREGISTRATION.md` before any run. Population
and contract identical to R7/R8; three seeds, matched 1200-step budget.
`meta_test` not opened.

## The question (posed before the audit)

R8's B1 improved the shape term (0.896, best on record) and the k=5
activity-cliff ordering (0.768 vs A0's 0.675) while the global CI regressed
(0.535 vs 0.580). Where does the CI loss live, pair by pair?

## Pair audit (no training; A0 vs B1, k=0)

| stratum | pairs | A0 sign | B1 sign | per-target B1-A0 (bootstrap; + = B1 worse) |
|---|---:|---:|---:|---|
| all | 10,824 | 0.572 | 0.541 | +0.024 [-0.035, +0.088] |
| cliff (sim>=0.6, gap>=1.0) | 570 | 0.509 | **0.577** | -0.049 [-0.190, +0.096] |
| **mid_sim (0.4-0.6)** | 2,031 | 0.599 | 0.554 | **+0.119 [+0.022, +0.220]** |
| low_sim (<0.4) | 5,964 | 0.644 | 0.584 | -0.022 [-0.117, +0.074] |
| mid_gap (0.5-1.0 pK) | 2,535 | 0.544 | 0.491 | +0.120 [-0.008, +0.263] |
| small_gap (<0.5 pK) | 3,405 | 0.522 | 0.478 | +0.021 [-0.029, +0.070] |
| top5 pairs | 1,149 | 0.483 | 0.491 | +0.014 [-0.082, +0.113] |

The only component-resolved stratum is the **mid-similarity band** — the
band immediately below the activity-cliff weight's discontinuity at
Tanimoto 0.6 — while cliff pairs themselves improved. This is the
discontinuity hypothesis, measured rather than assumed.

## Dose response (cliff_pair_weight in {4.0=B1, 2.0=C2, 1.0=C1}, 3 seeds)

| arm | w | k=0 MSE | CI | Spearman | calib | shape | k=5 MSE | k=5 cliff sign |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A0 | — | 2.149 | **0.580** | **0.223** | 1.236 | 0.913 | 0.915 | 0.675 |
| B1 | 4 | 2.167 | 0.535 | 0.096 | 1.271 | **0.896** | 1.072 | 0.768 |
| C1 | 1 | 2.235 | 0.562 | 0.159 | 1.332 | 0.903 | 1.096 | **0.782** |
| C2 | 2 | **2.119** | 0.548 | 0.126 | **1.218** | 0.901 | 1.089 | 0.775 |

Paired bootstraps vs A0 (k=0 MSE): B1 -0.018 [-0.243, +0.229]; C1 -0.086
[-0.334, +0.171]; C2 +0.029 [-0.233, +0.317] (sign convention: negative =
arm better). All unresolved.

**The dose response confirms the audit and adds two findings:**

1. The x4 cliff weight is a **net negative for the ranking itself**: C1
   beats B1 on the global CI (+0.027) *and* on cliff pairs (0.606 vs 0.577
   pooled). The cliff-ordering ability comes from the shape-first training,
   not from the cliff emphasis; the emphasis distorts the shape fit even
   for the pairs it targets.
2. The weight trades MSE against CI (w=2 minimizes k=0 MSE at 2.119 with
   the family's best calibration 1.218; w=1 maximizes CI at 0.562). The
   C1-vs-A0 audit shows **no stratum remains resolved** (mid_sim +0.068
   [-0.032, +0.183], mid_gap +0.113 [-0.005, +0.255]) — the cliff-weight
   removal closed the resolved band; the remaining CI gap is diffuse, and
   C1's margins are still compressed (0.097 vs A0's 0.121).

## Gates

- Z1' (arm k=0 >= -2% vs A0): C2 passes at -1.4% (point estimate,
  unresolved); B1 and C1 fail.
- Z5' (CI no more than 0.02 below A0): **C1 passes** (0.562 vs 0.580,
  -0.018); B1 and C2 fail.
- No single dose passes both simultaneously; per the preregistered decision
  rule the cliff-weight axis is therefore recorded as **partially
  resolved**: the CI-optimal dose (w=1) and the MSE-optimal dose (w=2)
  split the gates, and the remaining CI deficit is no longer attributable
  to the cliff weight.

## Next single-variable hypothesis (R10, preregistered before any run)

The audit narrows the remaining CI gap to the margin compression (C1 margins
0.097 vs A0 0.121; small/mid-gap strata still below A0). The leading suspect
is the shape variance term (weight 1.5), which pulls the shape branch toward
the centered labels and flattens it on uncertain pairs. R10 changes exactly
one variable on the C1 base (cliff weight 1.0): `shape_variance_weight
1.5 -> 0.5`. Gates: CI improves over C1 and is no more than 0.02 below A0;
k=0 MSE does not regress beyond C1's 2.235; k=5 cliff sign stays >= 0.70;
three seeds, component bootstrap.
