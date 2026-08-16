# Stage C: train-only retrieval prior blended into the zero-shot endpoint

> **Corrected 2026-08-15 by the Stage R0 audit.** Evidence grade of this stage is
> **development, conditional on ligand overlap** — not confirmed. Three specific
> corrections, all verified in
> `../stageR0_retrieval_falsification_20260815/AUDIT_VERIFICATION.json`:
>
> 1. `meta_val` has **10** homology components in this bank, not 11.
> 2. `w = 0` reproduces the accepted model **to within 0.03 in the transport
>    scale**, not exactly: this script uses a fixed `beta = 8` while the three
>    checkpoints learned 7.9743 / 7.9849 / 7.9897.
> 3. 305 of the 624 query cells (48.9%) contain a ligand that appears verbatim in
>    `meta_train`. On the 12 targets with **no** exact ligand overlap the k=0
>    reduction falls to **+0.050 [-0.074, +0.175]** — a quarter of the headline
>    effect, and unresolved. The 12.3% below therefore describes a population
>    that is half exact recall.
>
> All numbers below are the original artifacts and are unchanged.

Numerical authority: `BLEND_meta_val.json` (+`.rows.jsonl`),
`BOOTSTRAP_meta_val.json`. Selected by the Stage A decision table.
Population: `meta_val`, all eligible targets, nested k = 0/1/2/3/5, the three
accepted `similarity_only` checkpoints. `meta_test` untouched.

## Intervention

```text
f0'(q) = (1 - w) * f0(q) + w * retrieval(q)          w = 0.5
r_k    = y_k - f0'(L_k)
f(q)   = f0'(q) + s(n) * sum_k softmax_k(8*Tanimoto) r_k    (transport unchanged)
```

`retrieval(q)` is per-query, `meta_train`-index only, and uses no query label,
no target identity and no query-set statistic. Three sources were compared:
`ligand` (Tanimoto kNN over meta_train ligand means, beta=24), `dual` (protein
kNN then ligand kNN within neighbours), and `blend` (their mean).

No training. `w = 0` reproduces the accepted model to within the transport-scale
difference noted above (1.6123/1.3071/1.0641/1.0187/0.8959 at k=0/1/2/3/5). It is
a close numerical check, **not** an exact identity: an exact check requires
reading `transport.similarity_scale` from each checkpoint instead of the fixed
`--transport-beta 8`.

## Result, `meta_val`, three seeds pooled

| k | accepted (`w=0`) | `blend_w0.5` | reduction | CI change | Spearman change |
|---|---:|---:|---:|---:|---:|
| 0 | 1.6123 | **1.4148** | **12.3%** | +0.026 | +0.081 |
| 1 | 1.3071 | **1.0437** | 20.2% | +0.026 | +0.081 |
| 2 | 1.0641 | **0.8841** | 16.9% | -0.001 | +0.007 |
| 3 | 1.0187 | **0.8045** | 21.0% | +0.004 | +0.005 |
| 5 | 0.8959 | **0.6871** | 23.3% | -0.002 | +0.013 |

Per-seed k=0 MSE reduction: **+0.164 / +0.275 / +0.154 — positive in 3/3**.
Per-seed k=0 CI change: **+0.025 / +0.030 / +0.022 — positive in 3/3**.
Per-seed k=5 MSE reduction: +0.208 / +0.240 / +0.178 — positive in 3/3.

## Paired component-level bootstrap, 9,999 draws

| k | metric | mean | 95% CI | LB>0 |
|---|---|---:|---|---|
| 0 | MSE reduction | +0.198 | **[+0.012, +0.416]** | **yes** |
| 0 | CI | +0.026 | [-0.083, +0.113] | no |
| 2 | MSE reduction | +0.180 | [-0.019, +0.425] | no |
| 3 | MSE reduction | +0.214 | [-0.027, +0.545] | no |
| 5 | MSE reduction | +0.209 | [-0.037, +0.594] | no |

The intervals are wide because the retrieval advantage is strongly heterogeneous
across the 10 `meta_val` homology components, which is exactly what the ligand
novelty stratification predicted. Recomputed per component, **6 of 10 improve**
and 4 regress (-0.155 to -0.012).

## Decision

**Accept the k=0 claim as development evidence. Do not claim the k>=2 gains.**
The gates below were the ones preregistered for this stage; the Stage R0 audit
adds an exact-ligand-exclusion gate that this stage did not have and does not
pass (see the correction block at the top).

| gate | outcome |
|---|---|
| >= 5% k=0 MSE reduction | **pass** — 12.3% |
| no k>=2 regression | **pass** — every k improved |
| non-negative CI/Spearman | **pass** on Spearman at every k and on CI at k=0/3; CI is -0.0013 at k=2 and -0.0021 at k=5, i.e. flat, not a ranking collapse |
| positive component-bootstrap lower bound | **pass at k=0 only** (+0.012); k>=2 intervals cross zero |
| consistent direction in every seed | **pass** — 3/3 for k=0 MSE, k=0 CI and k=5 MSE |

So the supported statement is: **a training-free, `meta_train`-only retrieval
prior blended at w=0.5 into the zero-shot endpoint reduces cold-target k=0 MSE
by 12.3% on `meta_val`, with a positive component-bootstrap lower bound and
consistent direction in all three seeds, without degrading ranking** — on a
population in which 48.9% of query cells contain a ligand that appears verbatim
in `meta_train`, with hyperparameters selected on that same population. On the
exact-ligand-free subset the effect is +0.050 [-0.074, +0.175] and unresolved.

The k>=2 point estimates (17-23%) are consistent in all three seeds but do not
clear the component bootstrap on 11 components, and are **not claimed**.

## The `dual` alternative, and why it was not selected

`dual_w0.5` gives the better ranking (k=0 CI +0.056, k=3 CI +0.049
[+0.007,+0.094] and Spearman +0.119 [+0.035,+0.212], both with positive lower
bounds) and positive MSE lower bounds at k=2/3/5, but **no k=0 MSE gain**
(-0.002, and negative in 2 of 3 seeds). Since Stage C was selected to fix k=0,
`blend` is the arm that answers the question. `dual` is recorded as the better
choice if ranking rather than k=0 error becomes the objective.

## Standing caveats

* **Ligand novelty.** The CD-HIT40 split is component-hard on proteins only.
  Stage A showed retrieval collapses to the global mean for the most novel
  ligands (novelty < 0.4) and delivers its gain at novelty >= 0.6. Roughly half
  this improvement is genuine chemical-neighbourhood transfer and roughly half
  is ligand recall on a split that never controlled for ligand overlap. Both are
  protocol-legal; only the first is a capability. **This sentence must accompany
  the 12.3% figure.**
* `beta = 24`, `w = 0.5` and the `blend` source were selected on `meta_val`, the
  development split. Legal, but tuned, and `meta_test` has not been used to
  confirm them.
* The prior is a *retrieval memory*, not a learned mechanism. It adds no
  parameters and no training, and it does not make the trained trunk better —
  it compensates for the trunk's near-total lack of within-target ranking
  (Stage A: model CI 0.525 against a 0.500 coin flip).

## Resources

Training-free. One frozen forward pass per episode per checkpoint plus a
meta_train index of 7,661 ligands and 399 targets. No GPU training, no new
parameters, ~8 GB envelope untouched.
