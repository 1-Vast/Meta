# Stage R0: the retrieval prior is falsified

Numerical authority: `FALSIFICATION_meta_val_stage10bank.json` (+`.rows.jsonl`),
`FALSIFICATION_meta_val.json` (secondary population),
`AUDIT_VERIFICATION.json`. Gates fixed in `PREREGISTRATION.md` before any result
existed. No training; `meta_test` untouched.

## Result

**All five preregistered gates fail.** Under the decision rule that means:
*retrieval is retained only as a protocol-specific baseline. It does not become
part of a core innovation, and Stage 10's k=0 result is recorded as an artefact
of ligand overlap plus selection bias.*

Primary population = the exact Stage 10 k=0 bank (44 targets, 10 components, 624
query cells; 3 checkpoints -> 1,872 rows). At k=0 the support transport is
inactive, so `w = 0` is exactly the checkpoint endpoint `f0`.

| stratum | rows | `f0` MSE | outer-fold selected | dCI | dSpearman |
|---|---:|---:|---:|---:|---:|
| all | 1872 | 1.6123 | 1.5884 | +0.018 | +0.055 |
| **exact_free** | 957 | **2.8019** | **3.0193** | -0.038 | -0.096 |
| exact_overlap | 915 | 1.3581 | **1.0114** | +0.017 | +0.070 |
| scaffold_disjoint | 894 | 2.5523 | 2.5895 | +0.022 | +0.072 |
| scaffold_overlap | 978 | 1.3993 | 1.4657 | +0.030 | +0.045 |
| tanimoto < 0.40 | 669 | 4.4967 | 4.9126 | -0.010 | -0.030 |
| tanimoto 0.40-0.60 | 180 | 2.0658 | 2.5707 | +0.034 | +0.080 |
| tanimoto 0.60-0.80 | 102 | 1.1188 | 1.4530 | +0.118 | +0.207 |
| near-duplicate (>= 0.80) | 921 | 1.3273 | **0.9842** | +0.032 | +0.093 |

The pattern is unambiguous. Every stratum where the query ligand is genuinely
new is **worse** with the prior; the entire benefit is concentrated in the
near-duplicate and exact-overlap strata. Aggregated over all cells the effect
collapses to +0.024 [-0.496, +0.471] — the Stage 10 headline of +0.198 does not
survive honest selection.

| gate | outcome |
|---|---|
| G1 exact-free component lower bound > 0 | **fail** — -0.217 [-0.785, +0.261] |
| G2 low-Tanimoto not worse | **fail** — 4.4967 -> 4.9126 |
| G3 CI/Spearman >= -0.01 | **fail** — -0.038 / -0.096 |
| G4 correct protein beats shuffled/random/matched | **fail** — see below |
| G5 exact-free gain >= 40% of all-cell gain | **fail** — ratio -9.1 |

## Two things this stage measured that no earlier stage could

**1. The size of the selection bias.** The same 200-configuration search run on
the whole population instead of nested folds selects `ligand, beta=16, w=0.75`
and reports exact-free MSE **2.5514** against `f0`'s 2.8019 — an 8.9%
improvement. Honest leave-one-component-out selection on the identical data and
identical search space gives **3.0193**, a 7.8% *degradation*. The gap of 0.468
MSE **is** the selection bias, measured rather than assumed. The configuration
that wins globally is not stable across components: `f0`'s exact-free MSE ranges
from 0.49 to 8.31 across the eight components that have exact-free cells, and
only 3 of 8 improve under their own outer-fold configuration.

**2. Protein identity is not used at all.** In 9 of 10 outer folds the selected
source is `ligand` — protein-blind. A protein-blind predictor cannot exhibit
protein specificity, so its counterfactual contrast against shuffled, random and
similarity-matched protein neighbours is exactly zero by construction, and the
one fold that did choose a protein-conditioned source (`blend_centered`) scores
*worse* than its shuffled counterfactual (-0.113 [-0.340, 0.000]). Adding the
train-centred ESM retriever to the search space — the sharper representation
that finding 10 of the audit exposed — did not change this: sharper protein
neighbours still lose to no protein neighbours at all.

**Consequence: no protein-conditioned language may be attached to this prior**,
and `dual`/`blend` are dropped.

## The one exploratory signal worth carrying forward

Activity cliffs — within-target ligand pairs with Tanimoto >= 0.6 and an
affinity gap >= 1.0 pK — are where chemistry must be read rather than recalled.
82 such pairs exist in this bank, 62 of them exact-free.

| ordering accuracy | all cliff pairs | exact-free cliff pairs |
|---|---:|---:|
| trained endpoint `f0` | 0.6211 | **0.5192** |
| retrieval prior | 0.6965 | **0.7162** |

**On the pairs that matter most, the trained zero-shot endpoint is at chance.**
A fixed Morgan/Tanimoto prior with no parameters is not. This is a diagnostic
about the trunk, not evidence for retrieval: it was not a preregistered gate,
62 pairs is a small sample, and the same prior degrades MSE on the same cells.
Graded **exploratory**. It is the sharpest available statement of what the model
is missing, and it sets the target for Core Innovation A.

## Secondary population

Repeating everything on the k=0-only bank (50 targets, 11 components, 751 cells)
reproduces the conclusion: exact-free 2.6904 -> 3.0016, same gate outcomes. The
result is not an artefact of the nested-bank eligibility rule.

## Consequences for the rest of the cycle

1. **Stage 10's accepted k=0 claim is withdrawn.** The supported statement is now
   that a train-only retrieval prior reduces k=0 MSE **on query ligands that
   already appear in `meta_train`** (1.3581 -> 1.0114) and increases it on ligands
   that do not (2.8019 -> 3.0193).
2. **Retrieval is out of Core Innovation A.** The factorised predictor must earn
   `ligand_prior` with learned parameters, not a memory. Retrieval survives only
   as a named baseline in comparison tables.
3. **The double-cold protocol (Stage R1) is now mandatory, not optional.** 48.9%
   exact ligand overlap made the development split incapable of separating recall
   from capability, and every earlier stage inherited that.
4. **Protein representation is not the k=0 lever either.** The audit reopened it;
   R0 tested it with a sharper retriever and it lost. Stage R2 keeps a protein
   arm, but the prior on it is now low.
5. **The trunk's failure is specifically within-target ordering**, and it is
   visible at chance level on activity cliffs. That is what Innovations A and B
   must attack.

## Resources

Training-free. Three frozen forward passes per episode plus a `meta_train` index
of 7,661 ligands and 399 targets; 200 configurations x 10 outer folds evaluated
in numpy. No GPU training, no new parameters.
