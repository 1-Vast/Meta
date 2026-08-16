# Stage A: what dominates cold-target k=0 error

> **Corrected 2026-08-15 by the Stage R0 audit** (verified in
> `../stageR0_retrieval_falsification_20260815/AUDIT_VERIFICATION.json`):
>
> 1. The composed `retrieval[shape|level]` row is built as
>    `shape - shape.mean() + level.mean()`, where both means are taken over the
>    **query panel of the episode**. It is therefore a **transductive diagnostic
>    upper bound**, not an ordinary per-query predictor, and the 25.5% figure is
>    **exploratory** evidence about available headroom. The single-source rows
>    (`ligand_neighbor_*`, `dual_neighbor_*`, `model_f0`) are per-query and are
>    unaffected.
> 2. "Protein retrieval is weak" holds **only for raw pooled ESM cosine**, whose
>    similarities occupy a band of width 0.21 around 0.90 with a mean spread of
>    0.024 across the 16 nearest training targets. A `softmax(16 * sim)` over that
>    band is nearly uniform, so `protein_neighbor_esm` was close to a global mean
>    by construction. Centring on a `meta_train`-only mean widens the spread to
>    0.238. Finding 3 below is corrected accordingly: it rules out *this*
>    protein retriever, not protein representation in general.
>
> All numbers below are the original artifacts and are unchanged.

No training. Numerical authority: `DECOMPOSITION_meta_val.json`
(+`.rows.jsonl`), `NOVELTY_meta_val.json`. Gates fixed in `PREREGISTRATION.md`.
Population: `meta_val`, 50 k=0 episodes, all eligible targets. Every retrieval
index contains **meta_train records only**. `meta_test` was not touched.

## Headline

**k=0 error is 59% target-level calibration and 41% within-target shape, and the
trained zero-shot endpoint contributes essentially nothing to shape.** A
train-only *transductive* composition of two retrievers beats the trained
endpoint by **25.5% MSE** (exploratory upper bound; see the correction block).
The strongest per-query train-only predictor is `ligand_neighbor_b24` at 1.6250,
a **10.8%** reduction, with concordance 0.635 against 0.525.

## The exact decomposition

`MSE(t) = (mean(p-y))^2 + var(p-y)` = calibration + shape, aggregated
equal-component-then-equal-target.

| estimator | MSE | calibration | shape | CI | Spearman |
|---|---:|---:|---:|---:|---:|
| `global_mean` | 2.4195 | 1.6765 | 0.7430 | 0.500 | — |
| `protein_neighbor_kmer` | 2.1257 | 1.3827 | 0.7430 | 0.500 | — |
| `protein_neighbor_esm` | 2.0785 | 1.3355 | 0.7430 | 0.500 | — |
| `ligand_prior` | 2.0127 | 1.2599 | 0.7528 | 0.611 | 0.355 |
| `dual_neighbor_b24` | 1.9244 | 1.2651 | **0.6593** | 0.660 | 0.392 |
| **`model_f0`** | **1.8212** | **1.0809** | 0.7403 | 0.525 | 0.075 |
| `ligand_neighbor_b24` | 1.6250 | **0.6967** | 0.9283 | 0.635 | 0.313 |
| `retrieval[dual_b24 shape + ligand_b24 level]` *(transductive)* | 1.3560 | 0.6967 | 0.6593 | 0.660 | 0.392 |
| `model_f0_oracle_level` *(oracle)* | 0.7403 | 0 | 0.7403 | 0.525 | 0.075 |
| `target_oracle` *(oracle)* | 0.7430 | 0 | 0.7430 | 0.500 | — |

## Five findings

**1. Calibration dominates.** For the accepted model, 1.081 of 1.821 (59%) is
target-level miscalibration. Knowing only the target's mean removes it.

**2. The model has no within-target ligand discrimination.**
`model_f0_oracle_level` (0.7403) is statistically indistinguishable from
`target_oracle` (0.7430) — the model's ligand ordering improves on a flat
constant by 0.4%. Its concordance is 0.525 against a 0.500 coin flip and its
Spearman is 0.075. **This is the sharpest negative result in the project:** the
zero-shot endpoint is, for ranking purposes, a per-target constant.

**3. This protein retriever is weak; ligand retrieval is strong.** ESM-pooled and
3-mer protein kNN reduce calibration only from 1.677 to 1.336/1.383, and being
ligand-blind they cannot rank at all. Sharp ligand retrieval reaches calibration
0.697 — **better than the protein-conditioned trained model, using no protein**.
*Corrected scope:* the ESM retriever's weights are near-uniform because raw
pooled cosine is compressed into a band of width 0.21 (see the correction block),
so this finding indicts the retriever, not protein representation generally.

**4. Calibration and shape have different best sources.** Sharp ligand retrieval
(beta=24) wins calibration (0.697) but has poor shape (0.928); dual
protein x ligand retrieval wins shape (0.659) but poor calibration (1.265).
Composing them — level from one, shape from the other — gives **1.356**, a
**25.5% reduction** against the model. *This composition re-centres on the query
panel and is therefore transductive: it measures how much headroom exists if the
two sources could be combined, not a deployable per-query predictor.* The best
deployable train-only predictor here is `ligand_neighbor_b24` at 1.6250 (10.8%).

**5. Protein sensitivity is modest.** Swapping in a cross-component donor
protein moves the accepted model's zero-shot output by 0.21-0.24 pK.

## The decisive caveat: ligand novelty

The CD-HIT40 split is component-hard on **proteins only**. Ligands may recur
across splits, so retrieval can partly recall a ligand's typical potency.
Stratified by mean max-Tanimoto of query ligands to `meta_train`:

| ligand novelty | targets | `global_mean` | `model_f0` | `retrieval` | retrieval CI |
|---|---:|---:|---:|---:|---:|
| [0.00, 0.40) most novel | 14 | 2.747 | 3.086 | **2.754** | 0.608 |
| [0.40, 0.60) | 3 | 3.941 | 1.332 | **1.250** | 0.519 |
| [0.60, 0.80) | 8 | 1.874 | 1.243 | **1.052** | 0.741 |
| [0.80, 1.01) near-duplicate | 25 | 1.875 | 1.908 | **1.391** | 0.702 |

* Retrieval beats the trained model in **4 of 4** buckets, so its advantage is
  not purely a recall artefact.
* But for the most novel ligands it collapses to the global mean (2.754 against
  2.747): it degrades **gracefully** rather than helpfully. Its genuine
  capability lives at novelty >= 0.6, which is 33 of 50 targets here.
* The trained model is **worse than a constant** in the most novel bucket (3.086
  against 2.747) and ties it in the largest bucket (1.908 against 1.875).

So roughly half of the retrieval gain is real chemical neighbourhood transfer
and roughly half is ligand recall on a split that never controlled for it. Both
are protocol-legal — no target ID, no query label — but only the first is a
capability, and this must be stated whenever the 25.5% is quoted.

## Decision

The preregistered table row that fires is: *dual retrieval has the strongest
frozen oracle -> build a train-only, confidence-weighted protein-target x
ligand-neighbour residual memory for k=0, abstaining when neighbours are weak.*
Headroom is 25.5%, far above the 5% bar.

One refinement the data forces: the brief expected abstention **to `f0`** when
neighbours are weak. The novelty table shows `f0` is *worse than a constant*
exactly there, so abstention should target the ligand/global prior instead, and
soft retrieval already self-abstains toward it as similarities flatten.

`beta = 24` was selected on `meta_val`, the development split. That is legal but
is a tuned hyperparameter and is reported as such.

## What this rules out

* ~~Another protein-representation intervention for k=0.~~ **Withdrawn by the
  Stage R0 audit.** What is ruled out is *raw pooled ESM cosine retrieval*, whose
  weights are near-uniform by construction. Train-only centring or whitening was
  never tested and widens the usable spread tenfold, so protein representation
  for k=0 is reopened as a Stage R2 question. Mac-Diff, conformer routing,
  PBCNet2.0 and Cartesian equivariance remain closed: they were rejected on
  structural-input coverage and multi-seed training evidence, neither of which
  this correction touches.
* Adding capacity or training budget to the interaction trunk to fix k=0 ranking:
  the trunk's shape contribution is ~0 and prior capacity/budget increases did
  not move it.

## Resources

Training-free. One pass over 50 `meta_val` episodes plus a meta_train index of
7,661 ligands and 399 targets; seconds of CPU, no GPU beyond three frozen
forward passes per episode.
