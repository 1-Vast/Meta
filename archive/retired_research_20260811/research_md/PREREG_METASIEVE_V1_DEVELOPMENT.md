# MetaSieve v1 development preregistration

Date frozen: 2026-08-11

## Scope and claim boundary

This is a source/meta-validation repair experiment for the biological axis
that failed main-v0. The 50 main-v0 meta-test targets and their 1,934 cells are
physically excluded from v1 artifacts and must never be loaded by the trainer
or evaluator. Main-v0 test is permanently consumed.

Meta-validation is a development set with 37 k=5 tasks in nine eligible
CD-HIT40 clusters. Its maximum verdict is:

```text
V1_DEVELOPMENT_CANDIDATE_SELECTED
```

It cannot establish biological specificity or authorize production migration.
A later confirmation requires a fresh, exposure-logged set with at least 18
eligible independent protein clusters and at least 30 k=5 targets, disjoint
from source, meta-validation and main-v0 test.

## Fixed data and solver

- use the frozen main-v0 exact-Ki corpus and 288D T-BASIS;
- only `meta_train` is used for optimization;
- `meta_val` is evaluated once after all architecture, losses and seeds are
  frozen here;
- k=5, `d=2`, ridge `1.0`, five seeds `20260821..20260825`;
- sample `CD-HIT40 cluster -> target -> episode`, uniformly at each level;
- use the same optimizer, step count and paired episodes for all arms;
- the only support-conditioned state remains the ridge coefficient, whose
  dimension is at most `rank(M_S) <= min(d,k)`.

Cluster-balanced sampling is a MetaSieve repair, not an AdaMBind or CleanSplit
reproduction. It is applied identically to every arm so it cannot explain an
arm contrast.

## Nested arms

### V0

The existing frozen-T-BASIS model: ligand population prior plus the low-rank
`w0^T U^T phi` pair term and ridge section. It is retrained under the common
cluster-balanced sampler.

### V1-A: shared pair prior, MSE only

```text
h(P,L)  = phi(P,L) + W2 SiLU(W1 phi(P,L)) / sqrt(32)
mu(P,L) = f_L(L) + beta^T h(P,L)
m(P,L)  = U^T h(P,L), d=2
```

The bottleneck width is 32. `beta`, `W1`, `W2` and `U` are shared source-learned
parameters; they are not adapted from a new target's support labels. The v0
`population_coordinate` is removed because the full pair head subsumes it.

This trains a pair-specific readout over frozen T-BASIS. It is not an
end-to-end ESM/GINE/bridge update and must not be described as one.

### V1-B: measured-contrast auxiliary training

V1-B has exactly the V1-A architecture and initialization schedule. It adds two
source-only squared-difference losses, each with frozen weight `0.1`:

1. within-panel ligand contrast: same panel and target, two measured ligands;
2. measured partner contrast: same panel and ligand, two measured targets in
   different CD-HIT40 groups.

For either pair `(i,j)`, the loss is
`[(mu_i-mu_j) - (y_i-y_j)]^2` in source-standardized pKi units. It uses only
observed exact Ki values. An absent database edge, random protein, graph-distant
pair or docking pose is never assigned a label. The auxiliary losses train the
shared `h` and `mu`; no extra task state exists at inference.

## Physical isolation and controls

Build separate source cells/features, meta-validation episode inputs and
meta-validation query truth. The prediction process receives support labels but
has no query-label argument and writes predictions before the evaluator opens
query truth. Loaders reject every frozen main-v0 meta-test target and cell.

Wrong-protein controls use a frozen map `meta_val target -> meta_train donor`
from a different CD-HIT40 group, selected without affinity by nearest sequence
length and amino-acid composition. The old global, length-unmatched wrong arm
is forbidden because 74/399 source targets map to held-out donors.

Evaluate a 2x2 biology factorial with identical ligand episodes:

```text
support correct / query correct
support correct / query wrong
support wrong   / query correct
support wrong   / query wrong
```

Also evaluate ligand d0, pair d0, zero section and permuted support. Wrong
features are controls only and never enter a training loss.

## Development criteria

Report target-macro and CD-HIT40-cluster-macro MSE, R2, CI, Pearson and
Spearman. Report paired MSE reductions and one-sided 95% cluster-bootstrap
lower bounds.

V1-B becomes the selected development candidate only if all hold on meta-val:

1. correct absolute MSE is no worse than V1-A;
2. its cluster-macro correct-versus-permuted reduction exceeds V1-A's;
3. its cluster-macro correct-support/query versus correct-support/wrong-query
   reduction exceeds V1-A's;
4. correct beats zero and ligand d0 at cluster level;
5. finite gradients, support permutation invariance, `rank<=min(d,k)`, and no
   query-label or main-v0-test access are verified.

These are development selection rules, not confirmatory scientific Gates.
Failure stops end-to-end frontend work; success authorizes only a fresh-data
confirmation preregistration. Q-PMA, CSMO and production `model/` remain closed.

