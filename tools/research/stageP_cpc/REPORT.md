# Stage P: centered protein supervision does not produce protein-conditioned ordering

Verdict: **P1 FAILS. The family closes at the boundary the preregistration
names, and no admission stage runs.**

Authority: `STAGE_P_meta_val.json` (+ `.rows.jsonl`). Preregistration:
`PREREGISTRATION.md`, frozen before any arm was trained. Population: double-cold
`meta_val`, 41 targets / 19 components, k ∈ {0,1,2,3,5}, three matched seeds,
component-paired bootstrap (9,999 draws). `meta_test` labels were used for no
fitting, selection or reported metric.

## What was asked

Whether the current `ContactGrammar` can learn useful protein-conditioned
within-target ordering when an objective explicitly demands it. Every
wrong-protein control in R0-R14 is computed on *uncentered* error, which the
additive `protein_value(P)` branch satisfies alone, so no objective had ever
asked (DATAFLOW_AUDIT F6/F7).

Two arms, identical in every respect except the objective:

| arm | protein contrast | weight | fires on |
|---|---|---:|---|
| `A0repro` | uncentered (incumbent) | 0.5 | k>0 episodes only |
| `CPCoverdrive` | **centered** | **2.0** | **every episode** |

Verified matched: the recorded configs differ in exactly
`protein_contrast_form` and `protein_contrast_loss_weight` and nothing else.
The evaluator refuses to score arms that differ in anything more.

## Primary gate

| gate | requirement | measured | verdict |
|---|---|---|---|
| **P1** | `r_correct(CPC) − r_correct(A0repro)` at k=0, component-paired lower bound above zero **and** mean ≥ +0.05 | **−0.0066 [−0.0545, +0.0417]** | **FAIL** |

Not merely short of the threshold — the point estimate is slightly negative and
the interval spans zero. The over-driven objective did not improve
correct-protein ordering at all.

## Why this is a clean negative rather than a null result

The decisive number is not P1. It is that **correct and wrong protein give the
same ordering, in both arms, at every k**:

| arm | k | `r` correct | `r` wrong | CI correct | CI wrong |
|---|---:|---:|---:|---:|---:|
| `A0repro` | 0 | +0.156 | +0.156 | 0.545 | 0.547 |
| | 2 | +0.238 | +0.239 | 0.566 | 0.566 |
| | 5 | +0.334 | +0.334 | 0.623 | 0.622 |
| `CPCoverdrive` | 0 | +0.149 | +0.148 | 0.526 | 0.524 |
| | 2 | +0.250 | +0.251 | 0.556 | 0.557 |
| | 5 | +0.331 | +0.331 | 0.612 | 0.611 |

Substituting a similarity-matched wrong protein changes within-target ordering
in the third decimal place. Four times the weight, on every episode, on an
objective algebraically barred from being satisfied by the level branch, moved
this not at all.

**The gap decomposition is therefore uninformative and is reported as such.**
Improvement −0.0066, donor degradation +0.0079, degradation share 0.54 — but
there is no gap to decompose. The preregistered concern was an arm that widens
the gap by damaging the donor; this arm did neither.

## The mechanism worked; the signal was not there

This is the part worth keeping. The objective did exactly what it was designed
to do:

**The level branch was excluded, as intended.** Gradient of the centered
contrast into `protein_head`: **8.06e-07** for CPC and 8.50e-07 for A0repro —
float32 zero. Centering removes `protein_value(P)` identically, and it did.

**The gradient reached the interaction path, and training amplified its
sensitivity.** Gradient norms of the centered contrast, averaged over 12
episodes and 3 seeds:

| branch | `A0repro` | `CPCoverdrive` | ratio |
|---|---:|---:|---:|
| `protein_encoder` | 1.67e-01 | 3.20e-01 | 1.9× |
| `ligand_encoder` | 1.40e-01 | 3.55e-01 | 2.5× |
| `grammar` (attention) | 1.23e-01 | 2.37e-01 | 1.9× |
| `embed` | 5.11e-02 | 2.35e-01 | **4.6×** |
| `interaction_head` | 3.48e-02 | 1.25e-01 | **3.6×** |
| `protein_head` (level) | 8.50e-07 | 8.06e-07 | — |
| `transport` | 0.00e+00 | 0.00e+00 | — |

**The protein-induced shift became reproducible across seeds.** Mean pairwise
cosine of the k=0 protein-induced shift vectors, 123 pairs:

| arm | seed cosine |
|---|---:|
| `A0repro` | −0.0591 (sd 0.576) |
| **`CPCoverdrive`** | **+0.3161** (sd 0.599) |

A0repro's protein response is undirected across seeds — the same signature a
random initialisation shows (−0.003, `NOISE_AND_LEAKAGE_AUDIT.md` §7).
CPCoverdrive's is *consistent*: three independently seeded runs move the same
way under the same protein substitution.

**But the shift is not aligned with truth.** Correlation of the protein-induced
shift with the centered labels: `A0repro` −0.007, `CPCoverdrive` **+0.022**.

So the objective produced a **reproducible protein-dependent response carrying
no affinity information**. Had P1 passed, this would have been caught by P3
(alignment ≥ +0.10, measured +0.022) — the preregistered "arbitrary but
consistent movement" failure mode. Because P1 failed first, the secondary gates
were not formally evaluated; these are observations, recorded because they name
the cause.

## No compensating regression, and no gain

| quantity, CPC vs A0repro at k=0 | measured |
|---|---|
| MSE | +0.0180 [−0.0822, +0.1293] |
| CI | −0.0186 [−0.0495, +0.0067] |
| Spearman | −0.0349 [−0.1033, +0.0233] |
| calibration | +0.0280 [−0.0724, +0.1369] |

All unresolved. The change is neither harmful nor useful.

## Conclusion, at exactly the preregistered scope

> **Centered-objective training on the current `ContactGrammar`, at a 1,200-step
> budget on the double-cold protocol with sequence + 2D inputs, does not produce
> protein-conditioned within-target ordering.**

This is **not** a claim that protein-conditioned architectures are impossible,
that the data cannot support protein conditioning under any method, or anything
about `meta_test`.

### What it does resolve

Stage P was designed to separate two explanations that fit every prior
measurement equally well: **(a)** no objective ever asked, versus **(b)** the
data does not contain the signal.

**(a) is now excluded.** The objective asked, at four times the weight, on every
episode, through a gradient route verified to reach the interaction path and
verified to exclude the level branch. The model *complied* — its protein
response became reproducible across seeds — and the response carried no
ordering information.

That shifts the weight of evidence toward **(b)**, without proving it: a
different architecture, a longer budget, or richer protein features could still
succeed. What is now measured is that *this* architecture, given an explicit and
over-driven demand, has no protein-conditioned within-target ordering to give.

## Stop rule applied

The preregistration says: P1 fails → stop, do not run the admission stage, do
not add post-hoc weights, schedules or arms. `CPCpos`, `CPCwrong`, `CPCrand` and
`A3perm` are **not** run.

## Cost, as measured

| quantity | value |
|---|---|
| runs | 6 (2 arms × 3 seeds), all 1,200 steps |
| measured rate | 1.14 s/step (the 40-step calibration's 0.70 s/step excluded validation passes) |
| wall time | ≈23 min per run, ≈2.3 h total |
| peak GPU memory | 2,697 MB |
| added parameters | **0** |
| added forwards | one wrong-protein forward on the 240 k=0 steps the incumbent skips |

## Commands

```bash
# per seed in {20260815, 20260816, 20260817}
conda run -n drug python -m scripts.train_qpsmp --arch similarity_only \
  --steps 1200 --seed <SEED> \
  --split-directory dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1 \
  --output report/meta_fewshot/stageP_cpc_20260816/A0repro_seed<SEED>

conda run -n drug python -m scripts.train_qpsmp --arch similarity_only \
  --steps 1200 --seed <SEED> \
  --protein-contrast-form centered --protein-contrast-loss-weight 2.0 \
  --split-directory dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1 \
  --output report/meta_fewshot/stageP_cpc_20260816/CPCoverdrive_seed<SEED>

conda run -n drug python -m tools.research.stageP_cpc.evaluate \
  --stage report/meta_fewshot/stageP_cpc_20260816 \
  --output tools/research/stageP_cpc/STAGE_P_meta_val.json

conda run -n drug python -m pytest tools/research/stageP_cpc/tests -q   # 15 passed
```
