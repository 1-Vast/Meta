# LOCK/CLOCK G0 preregistration

**Frozen before any new KirHub activity statistic is computed.**

## 1. Question and claim boundary

The question is whether task-matched amino-acid substitution geometry supplies target-side
information for target-specific ligand reordering that is not explained by KLIFS taxonomy, pocket
composition, exact aligned identity, pooled ESM-2, or coordinate corruption.

This is not a reproduction of LOCK-GP or CLOCK-GP. Jankowiak et al. study fixed-length local mutation
landscapes with at least about 1,800 measured sequences per landscape. CLOCK learns about 49k
parameters from positional Chroma structure embeddings across many mutation landscapes. The present
DTA data contain one wild-type sequence per target and multiple ligands, so neither the original GP
likelihood nor the learned CLOCK map is identified here.

G0 may establish only a label-free geometry result and a single-source, ligand-warm shared-panel
reordering mechanism result. It cannot establish strict dual-cold prediction, cross-family transfer,
or affinity improvement.

## 2. Fixed data

Source metadata are frozen in `manifests/lock_clock_g0_sources.json`.

The label-free stage uses the intersection of:

- human 85-position KLIFS pockets;
- the frozen KirHub pooled ESM-2 cache;
- the frozen KirHub strict homology-component registry.

Only if the label-free execution gate passes may the runner read KirHub Table S4. The KirHub stage
uses its existing target and ligand component folds, retains the previously registered 5--95 activity
window, and reports that training-target measurements at query ligands are used. It is an information
oracle, not a deployable predictor.

Seed is `1729`. No development, confirmation, Davis, ChEMBL confirmation, or sealed label may be read.

## 3. Frozen substitution geometries

Let `B` be BLOSUM50 in log-score form, including the gap token. Following the reference
implementation:

```text
L[a,b] = (B[a,b] - 0.5 * (B[a,a] + B[b,b]))
L       = L * 0.25 / abs(median(L))
C       = exp(L)
```

For aligned 85-residue pockets `x,y`:

```text
K_linear(x,y)    = mean_l C[x_l,y_l]
K_nonlinear(x,y) = product_l C[x_l,y_l]
K_LOCK(x,y)      = normalized(K_linear * (1 + K_nonlinear))
```

No exponent, kernel scale, rank, or bandwidth is fitted to activity.

The separately reported label-free position-specific arm uses
`alpha_l = 0.5 + (1 - H_l/log(21))`, with entropy `H_l` computed from pocket residues, and replaces
`C` by the Hadamard power `C ** alpha_l` at position `l`. This is called
`conservation_LOCK`, not CLOCK: it has no positional structure embedding.

The comparison geometries are:

1. fixed `LOCK`;
2. `conservation_LOCK`;
3. exact aligned one-hot identity;
4. BLOSUM-aware pocket composition;
5. frozen pooled ESM-2 cosine similarity.

No coordinate concatenation is allowed.

## 4. Frozen destruction controls

Every arm enters the same neighbour-profile estimator.

- `position_shuffled`: independently permute the 85 positions within every target, preserving its
  exact pocket composition.
- `sequence_shuffled`: derange intact pocket sequences across targets within KLIFS group.
- `blosum_permuted`: permute the 20 canonical amino-acid labels of `C`, keeping the gap fixed.
- `random_psd`: an RBF correlation matrix on seeded random amino-acid vectors, scaled to match the
  median off-diagonal correlation of `C`.
- `matched_wrong_target`: derange coordinates within `(KLIFS group, eligible-ligand quartile)`.
- `within_family_wrong_target`: derange coordinates within KLIFS family.

Fixed points are forbidden wherever a block contains at least two targets and are counted otherwise.

## 5. Stage G0-L: label-free audit

For every kernel report:

- PSD minimum eigenvalue after symmetrization;
- centered effective rank and fraction of centered kernel energy retained by 16 dimensions;
- centered-kernel alignment with family, group, composition, exact identity, and pooled ESM-2;
- fraction of centered Frobenius energy remaining after projection on family and composition kernels;
- top-8 neighbour overlap with pooled ESM-2;
- within-family non-constant pair fraction.

The label-free execution gate passes only if fixed LOCK:

- is PSD to tolerance `-1e-8`;
- has at least `0.05` residual energy after family and composition projection;
- has centered-kernel alignment at most `0.95` with pooled ESM-2;
- is non-constant for at least `0.80` of within-family target pairs.

The separate low-dimensional claim requires 16 dimensions to retain at least `0.80` of centered
kernel energy. Failure of that claim does not authorize increasing the dimension after seeing results.

## 6. Stage G0-R: fixed-estimator reordering audit

If G0-L passes, run the existing 5 x 5 target/ligand-fold shared-panel oracle. For a query target,
the estimator selects the top eight training targets under one similarity matrix and averages their
activity profiles with squared non-negative similarity weights. All arms use the same folds,
candidate set, missingness rule, ligand window, and Spearman calculation.

Two candidate sets are reported separately:

- `group`: training targets in the query KLIFS group;
- `family`: training targets in the query KLIFS family, requiring at least two candidates.

The independent unit is the strict full-sequence homology component. Target values are first averaged
within component. Contrasts use paired component bootstrap intervals. The measured MDE is reported
for every candidate set; no target- or cell-level pseudo-replication is allowed.

## 7. Frozen gates

The fixed LOCK mechanism gate requires all of:

- group-mode `LOCK - group_centroid` mean at least `max(0.03, measured MDE)` and LCB95 above zero;
- group-mode LOCK LCB95 above pooled ESM-2, aligned identity, composition, position shuffle,
  sequence shuffle, BLOSUM permutation, random PSD, and matched wrong target;
- family-mode LOCK LCB95 above the uniform family centroid and within-family wrong target;
- at least 70 paired homology components in family mode and family-mode MDE at most `0.03`.

`conservation_LOCK` is a separate hypothesis. It is credited only if its LCB95 is above both fixed
LOCK and composition. It is not allowed to rescue a failed fixed LOCK gate.

No bandwidth, exponent, top-k, activity window, fold, dimension, seed, or threshold may change after
the first formal statistic is computed.

## 8. Decisions

```text
LOCK_G0_LABEL_FREE_DEGENERATE_STOP
LOCK_G0_REORDERING_NOT_IDENTIFIED_STOP
LOCK_G0_COORDINATE_SIGNAL_SURVIVES__TRAIN_ONLY
```

Even a pass maps only to final category 2. A strict predictive stage remains blocked until a
multi-family, provenance-independent factorial substrate exists. True CLOCK remains blocked until
audited positional structure embeddings and enough independent property landscapes exist.
