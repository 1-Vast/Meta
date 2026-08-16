# LOCK/CLOCK G0 implementation amendment

**Frozen before G0-L is run.** This amendment resolves implementation details that were not explicit
in `lock_clock_g0_preregistration.md`; it does not change an arm, threshold, seed, or stop rule.

Bound preregistration SHA-256:
`0f5932ac130ecb70f127878f23d8b347794f959c1883a3924b3cc7a2a583c762`.

Bound source-manifest SHA-256:
`68d92261d0330ee57854bf09dadb28519ed6e15ae2b1d79dab114070ab58316b`.

The byte-level local inputs, including the G0-R component registry omitted from the source manifest,
are separately frozen in `manifests/lock_clock_g0_input_bindings.json`. The runner verifies every
stage-relevant file against that manifest before loading it.

The implementation identity is a deterministic bundle digest over the runner and the three imported
helper modules that supply profile weighting/selection, gene aggregation, eligibility, and fold
allocation. Every stage records the four individual file hashes and the bundle digest; any change
after G0-L is a hard stop before workbook or statistic access.

## Kernel definitions

`normalized(K)` means `D^-1/2 K D^-1/2`, where `D = diag(K)`. Both the minimum eigenvalue before
and after this normalization are reported.

Pocket composition is the 21-token frequency matrix `F`. Its registered BLOSUM-aware kernel is
`K_comp = normalize(F C F^T)`. The nuisance projection uses the joint span of the centered
21-column frequency matrix and the centered feature span of `K_comp`, not a scalar regression on
the composition kernel. This makes the adjustment conservative to the row scaling introduced by
diagonal kernel normalization.

The random-PSD control draws one eight-dimensional standard-normal vector per amino-acid/gap token
with seed `1729 + 13`. Its RBF exponent is chosen label-free so that the median off-diagonal token
correlation exactly matches the median off-diagonal correlation of `C`. This token matrix replaces
`C` in the otherwise unchanged fixed LOCK formula.

All randomizations operate on genes in lexicographic order. A block derangement first randomly orders
the members and then cyclically shifts that order. Thus blocks of size at least two have no fixed
point; singleton fixed points are retained and counted. The matched wrong-target quartile is the
deterministic rank quartile of each target's number of non-saturated KirHub measurements.

## Label-free statistics

All alignments use double-centered kernels and normalized Frobenius inner products. The centered
effective rank is `(sum lambda_i)^2 / sum lambda_i^2` over positive eigenvalues. Top-16 energy is
`sum_{i=1}^{16} lambda_i^2 / sum_i lambda_i^2`.

For the family-plus-composition residual, let `Z` contain centered KLIFS-family one-hot columns,
centered pocket-frequency columns, and the centered `K_comp` column span. With `P` the orthogonal
projector onto the numerical column span of `Z`, the reported fraction is
`||(I-P) K_centered (I-P)||_F^2 / ||K_centered||_F^2`.

Top-eight ESM overlap is the mean fraction of non-self neighbours shared with pooled ESM-2; ties are
resolved by gene order. A within-family pair is non-constant when its normalized similarity is below
`1 - 1e-12`.

## Reordering statistics

Every arm uses the same target folds, ligand folds, candidate restriction, activity window, top-eight
rule, squared non-negative weights, and target-level Spearman calculation. Within each
target-by-ligand-fold evaluation, the truth and predictions from **all** reported arms are intersected
before any Spearman value is computed. This common mask prevents arm-specific prediction coverage from
creating a paired gain. If any arm has a non-finite Spearman value on that mask, the entire
target-by-ligand-fold cell is removed for every arm. Uniform group/family centroids differ only in
using every eligible candidate with unit weight.

The five ligand-fold Spearman values are first averaged within target. Target values are then averaged
within strict full-sequence homology component. Every contrast is complete-case paired at that component
level. Its measured 80%-power MDE is
`(z_0.975 + z_0.80) * sd(paired component deltas) / sqrt(n_components)`.

Family-restricted results test only within-family weighting and cannot establish cross-family transfer.
KirHub is a single-source, ligand-warm information oracle and cannot separate assay, document,
preparation, construct, or provenance effects.

The registered gate remains a strict-component conditional bootstrap. A non-gating family-macro
bootstrap and family-level empirical MDE are also reported because components within a family reuse
candidate pools and measurements. Neither interval propagates noise from the shared training-target
profiles; all inferential statements are conditional on those observed profiles.

After G0-L passes, a non-inferential preparation stage may open Table S4 only to reproduce the
existing oracle's `min_ligands=20` eligibility, form the exposure-matched corruption stratum, and
freeze the G0-R gene/ligand order plus target-component and ligand-component folds. It computes no
reordering statistic. All coordinate and destruction kernels are generated on the frozen 372-gene
G0-L universe and are only subset to the resulting 353-gene oracle universe; conservation exponents,
random coordinates, and sequence/position corruptions are never rebuilt on the activity-selected
subset. G0-R refuses any universe, quartile, corruption digest, or fold mapping that differs from the
preparation manifest.

Both the preparation stage and G0-R recompute the complete G0-L report before workbook access. All
discrete fields must match exactly and numeric fields must agree within absolute and relative
`1e-12`; the four gate predicates are separately re-derived from the metrics. G0-R also independently
reconstructs both balanced fold maps from the frozen eligible universe rather than trusting editable
manifest fields.
