# CFRI Gate Z preregistration - ChEMBL-37 powered substrate

Registered on 2026-07-25 after the B0-only power audit and before any CFRI result on this substrate.

## Task and substrate

The task is drug-and-target dual-cold affinity prediction. Development targets are disjoint from
training by ChEMBL target id, UniProt accession and sequence-homology component. Development drugs
are disjoint by parent connectivity, Bemis-Murcko scaffold and exact Morgan Tanimoto >= 0.95. ChEMBL
document and assay overlap are also zero. Only exact pKi and pKd are used; pIC50 is excluded.

Registry SHA256: `357256be57210d6cb44a560809751df1a5577aa231e74d8d50b0df84424ae168`.
Ligand-feature SHA256: `14c5e8de3ec1629de0fb51ca275091a6cf51e0962542658e4b30b88e9209db30`.
Frozen target-feature SHA256: `c3ac20e8500b8428bbb776a851da496f5d524d2467f02bc46ac5c1febf57b5d0`.

The development set has 181 targets with at least four ligands in 165 independent homology
components. Confirmation labels remain unopened.

## Frozen model and training protocol

The model is `y(t,d) = b(d) + g(t,d)`. The shared B0 path uses only Morgan-1024 plus ten
train-standardized physicochemical descriptors. CFRI is the only candidate mechanism: a frozen ESM-2
target representation conditions a target-ligand residual trained against five-fold
component-cross-fitted B0 predictions. No target, ligand, family, document, assay, source or split id
is a feature.

The run uses seed 1729, 4000 steps for every B0 fit and every arm, five cross-fitting folds, Adam with
the checked-in learning rates, and at most 256 randomly sampled training rows per target episode.
The cap prevents quadratic pair-matrix memory growth; targets remain equally weighted and development
evaluation uses every ligand. The complete frozen arms are B0, T0, A0, I0, R0, CFRI,
CFRI-Tshuffle, CFRI-Lshuffle, CFRI-Tpool and CFRI-Trandom.

Loss weights remain `affinity=0.3`, `center=1.0`, `orth=1.0`, `cal=0.1`. No weight or capacity was
tuned on ChEMBL development results.

## Frozen power and decision rule

Four matched B0 retrains (1729, 2027, 4241, 5501) produced full-query development target-macro
Spearman 0.1165, 0.1145, 0.1285 and 0.1239. The exact Gate-Z B0 trainer was used. The grouped empirical
MDE80 is 0.0586, so the primary effect threshold is frozen at `max(0.03, 0.0586) = 0.0586`.

Gate Z passes only if all conditions hold:

1. CFRI minus B0 target-macro Spearman is at least 0.0586.
2. The paired homology-component bootstrap LCB95 is greater than zero.
3. The paired direction is positive in every development homology component.
4. CFRI RMSE is no more than 2% worse than B0.
5. Both target derangement and ligand derangement cause significant damage (LCB95 greater than zero).
6. At least one registered stratum has an improvement LCB95 greater than zero.

Failure stops three-seed CFRI, BERP-GS, confirmation evaluation and long training. No threshold
relaxation, capacity increase or extra development-tuned rescue is admissible.

`affinity_labels_read=true` for train/development only; `confirmation_labels_read=false`;
`sealed_test_consumed=false`.
