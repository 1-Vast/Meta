# MetaSieve-main v0 real-training report

Date: 2026-08-11

## Outcome

The preregistered target-level Gates all pass, but the result is not admitted as
independent biological specificity:

```text
registered target-level verdict  REAL_BIOLOGICAL_META_SECTION_V0_PASS
scientific admission verdict     BIOLOGICAL_SPECIFICITY_NOT_IDENTIFIED_CLUSTER_SENSITIVITY
production migration             NOT AUTHORIZED
```

Meta-effect, correct-support specificity, and the full-versus-ligand-only
contrast survive a CD-HIT-cluster bootstrap. Correct protein versus wrong
protein does not. The 33 eligible test targets occupy only six CD-HIT clusters;
one cluster contains 21 targets and drives the target-level wrong-protein
contrast. The other five cluster means favor the wrong-protein control.

## Frozen protocol

- exact positive Ki; CARA-referenced measurement cleaning and a declared
  MetaSieve equal-panel median for protein-task labels;
- one protein sequence per task;
- CD-HIT 4.8.1 at 40% identity, complete-cluster 8:1:1 assignment;
- k=5, five training seeds and five paired test support draws;
- frozen 288D T-BASIS, differentiable `d<=5` ridge section;
- selected on meta-validation only: `d=2`, ridge `1.0`;
- 285/37/33 k=5 train/validation/test tasks;
- CSMO law metrics `NA_NOT_ADMITTED`; no 28D biological z map was invented.

## Point results

All values are target-macro over the same paired episodes.

| Arm | MSE down | RMSE down | R2 up | CI up | Spearman up | Pearson up | Band loss |
|---|---:|---:|---:|---:|---:|---:|---:|
| Population d=0 | 8.711 | 2.625 | -24.132 | 0.529 | 0.079 | 0.082 | NA |
| Ligand-only section | 3.426 | 1.725 | -5.012 | 0.517 | 0.050 | 0.059 | NA |
| Full correct | **1.916** | **1.298** | **-1.244** | 0.531 | 0.086 | 0.097 | NA |
| Full zero | 7.926 | 2.500 | -21.320 | 0.523 | 0.065 | 0.060 | NA |
| Full foreign | 6.466 | 2.226 | -14.585 | 0.521 | 0.062 | 0.065 | NA |
| Full permuted | 2.047 | 1.339 | -1.467 | 0.511 | 0.036 | 0.035 | NA |
| Full wrong protein | 2.113 | 1.353 | -1.683 | **0.537** | **0.103** | **0.103** | NA |

Absolute generalization is weak: full-correct R2 remains negative and ranking
correlations are about 0.1. Wrong-protein ranking metrics also exceed the
correct-protein values. This is a mechanism-screen result, not a competitive
DTA performance claim.

## Paired Gates

The contrast is `MSE(control) - MSE(full correct)`; positive favors the full
correct arm. The target analysis follows the preregistration. Cluster results
are the required dependence sensitivity.

| Control | Target mean | Target one-sided 95% LCB | Target Gate | Cluster mean | Cluster LCB | Cluster sensitivity |
|---|---:|---:|---:|---:|---:|---:|
| Population d=0 | 6.795 | 4.179 | PASS | 2.445 | 0.694 | PASS |
| Zero support | 6.010 | 3.750 | PASS | 2.154 | 0.665 | PASS |
| Foreign support | 4.550 | 3.199 | PASS | 3.543 | 2.078 | PASS |
| Permuted support | 0.130 | 0.060 | PASS | 0.145 | 0.087 | PASS |
| Ligand-only | 1.509 | 1.046 | PASS | 0.950 | 0.452 | PASS |
| Wrong protein | 0.197 | 0.020 | PASS | **-0.081** | **-0.227** | **FAIL** |

## Interpretation and next Gate

The closed-form section is doing real support-conditioned work on this split;
the result is not explained by no support, foreign support, permuted labels, or
a capacity-matched ligand-only family. What remains unidentified is whether the
frozen protein-ligand coordinate contributes reproducibly across independent
protein families.

Do not add Q-PMA or migrate this implementation to production. The next
experiment must be preregistered on source/meta-validation data only and target
the failed axis directly: train the biological coordinate for partner
specificity, then evaluate on fresh protein clusters. The consumed main-v0 test
targets may be used only descriptively for amendments.

## Audit artifacts

- machine result: `MAIN_V0_RESULT.json`
- paired draw metrics: `draw_metrics.json`
- predictions written without query labels before scoring:
  `predictions_before_query_labels.jsonl.gz`
- result SHA256: `c119407a729f181e30cf1eb24e06368f0d229456f355d8c9d10285c286c29642`
- prediction SHA256: `cdc9edbf4ad1b2d0886bd864e38099ff48274d0e0027641cf237bc49e3db01c6`
- corpus SHA256: `8f223a276840153c99671cfa64c8195f947e66e9b1d2a455fbc5137976a4f4fa`
- feature SHA256: `3d5303871ab1428e44c06c16c58e6ba329bc3ee61691868ecb925cb849cb4311`

