# Stage R0 preregistration: falsify the retrieval prior

Written before any Stage R0 result exists. No model is trained or modified.
Every retrieval index contains `meta_train` records only.

## Why this stage exists

The Stage 10 k=0 result (12.3%) was measured on a population where 48.9% of
query cells contain a ligand present verbatim in `meta_train`, with `beta`, the
retrieval source and `w` selected on the same population the bootstrap
resamples. On the 12 exact-free targets the effect is +0.050 [-0.074, +0.175].
Stage R0 decides whether any of that survives when exact recall is removed and
selection is separated from inference.

At k=0 the support transport is inactive, so `w = 0` is **exactly** the
checkpoint's zero-shot endpoint `f0`. The Stage 10 `beta=8` defect (finding 9)
does not touch this stage.

## Population

`meta_val`, `evaluation_seed=73101`, k=0 bank: 44 targets, 10 homology
components, 624 query cells. Baseline `f0` comes from the three accepted
`similarity_only` checkpoints; per-target values are averaged over the three
seeds before component aggregation, so **all intervals are conditional on those
three trained seeds**. `meta_test` is not touched.

## Strata, all computed per query cell

Novelty is **never** classified by a target-average similarity.

| stratum | definition |
|---|---|
| `exact_overlap` | the query ligand id occurs in `meta_train` |
| `exact_free` | it does not — **the primary population** |
| `scaffold_overlap` | its Bemis-Murcko scaffold occurs in `meta_train` |
| `scaffold_disjoint` | it does not |
| `tanimoto_lt40` | per-query max Morgan(r=2,1024) Tanimoto to `meta_train` < 0.4 |
| `tanimoto_40_60`, `tanimoto_60_80` | the two middle bands |
| `near_duplicate` | max Tanimoto >= 0.8 |
| `activity_cliff` | reported over within-target ligand **pairs** with Tanimoto >= 0.6 and abs label gap >= 1.0 pK; scored as the fraction of such pairs ordered correctly |

## Retrieval arms

Sources: `ligand` (protein-blind Tanimoto kNN over `meta_train` ligand means),
`dual[arm]` (top-16 protein neighbours, then ligand-weighted within them), and
`blend[arm]` (their mean).

Protein arms, all using the same top-16 machinery so only the selection of
neighbours differs:

| arm | neighbours |
|---|---|
| `correct` | top 16 by raw pooled ESM cosine |
| `centered` | top 16 by cosine after subtracting the `meta_train` mean vector (a train-only statistic) |
| `shuffled` | the similarity vector is permuted across `meta_train` targets |
| `random16` | 16 uniformly random `meta_train` targets, equal weight |
| `matched` | for each correct neighbour at similarity `s`, a random non-neighbour whose similarity is closest to `s` — same similarity level, different identity |

`correct` and `centered` are candidate mechanisms. `shuffled`, `random16` and
`matched` are falsifiers: if the prior works without correct protein identity,
it is a ligand memory and must not be described as protein-conditioned.

## Nested leave-one-component-out selection

For each of the 10 components `c`:

1. select `source` in {ligand, dual_correct, dual_centered, blend_correct,
   blend_centered}, `beta` in {8, 16, 24, 32}, `w` in {0, 0.25, 0.5, 0.75, 1.0}
   and `confidence` in {fixed, novelty_gated} using **only** the other nine
   components, scored by `exact_free` k=0 MSE;
2. apply that one configuration to `c` and record its predictions;
3. never look at `c` before step 2.

Only the pooled outer-fold predictions are analysed. No global winner is chosen
and then bootstrapped on the same data. `novelty_gated` sets
`w_q = w * clip((novelty_q - 0.2)/0.4, 0, 1)`, so the prior abstains toward `f0`
when the nearest training ligand is far.

## Gates, fixed now

| # | gate | threshold |
|---|---|---|
| G1 | `exact_free` k=0 MSE improvement over `f0`, outer-fold pooled | component bootstrap **lower bound > 0** |
| G2 | `tanimoto_lt40` MSE | not worse than `f0` |
| G3 | CI and Spearman on `exact_free` | change >= -0.01 |
| G4 | correct protein beats `shuffled`, `random16` and `matched` on `exact_free` | positive point estimate, and component lower bound > 0 against `shuffled` |
| G5 | the gain is not entirely exact recall | `exact_free` improvement >= 40% of the all-cell improvement |

9,999 paired component bootstrap draws, components resampled with replacement.

## Decision rule

* **All of G1-G5 pass** — the retrieval prior may enter Core Innovation A as an
  abstaining prior, and the protein arm that passed G4 is carried into Stage R2.
* **G1 passes, G4 fails** — retrieval is a real but purely *chemical* prior. It
  may still enter the model, but no protein-conditioned language is permitted
  anywhere and `dual` is dropped in favour of `ligand`.
* **G1 fails** — retrieval is retained only as a protocol-specific baseline. It
  does not become part of a core innovation, and Stage 10's result is recorded
  as an artefact of ligand overlap.

This stage cannot be rescued by adding arms after seeing results. Any arm not
listed above may only be reported as exploratory.
