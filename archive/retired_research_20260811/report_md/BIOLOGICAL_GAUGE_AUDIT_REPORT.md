# Biological gauge and partner-information audit

Date: 2026-08-11

## Terminal decision

```text
A0  NO_SINGLE_EXACT_ORTHOGONAL_GAUGE_IDENTIFIED
A1  TBASIS_SELECTIVITY_SIGNAL_NOT_IDENTIFIED
A2  NOT_RUN_GATE_CLOSED
```

The wrong/wrong recovery is not explained by one exact orthogonal change of
basis. Partial local support/query kernel alignment, scale changes and
population/section cancellation remain plausible contributors. The audit does
not use query affinity labels and cannot apportion their causal contributions.

The frozen calibrated 288D T-BASIS did not show a component-generalizing
measured-selectivity advantage under the preregistered capacity-matched ridge
probe. Therefore no partner anchor, frontend retraining, A2 localization,
Q-PMA, CSMO bridge or production migration is authorized.

## A0: gauge audit

The audit recomputed 925 sealed meta-validation episodes from five matching v0
checkpoints. They cover 37 episode-eligible targets and nine CD-HIT40 clusters.
The label-free global rotation audit separately covers all 3,150
meta-validation cells, 50 targets and 11 clusters.

Synthetic orthogonal controls preserve the ridge kernel and correction to
machine precision. Scale and shear controls correctly do not. Real
correct/wrong coordinates are far from an exact shared orthogonal gauge:

| Quantity | Median |
|---|---:|
| support Procrustes residual | 0.445 |
| held-out query transfer residual | 0.478 |
| support Gram relative error | 0.729 |
| `H_CC` versus `H_WW` relative error | 0.383 |
| leave-one-cluster global query residual | 0.801 |
| wrong/correct Gram trace ratio | 0.800 |

The wrong/wrong kernel is nevertheless closer to correct/correct than either
single-sided replacement (`0.383` versus `0.837/0.925`). This is evidence of
partial local consistency, not an exact or globally transferable gauge.

## A1: measured selectivity

The label-blind closure starts from 1,820 same-panel, same-ligand measured
cross-family groups. Repeated targets are aggregated within
`(panel, ligand, CDHIT40 family)`. Sharing a document, exact ligand, Murcko
scaffold, protein family or verified local identity at least 0.40 closes groups
into 21 components. The component MDE is 0.543, so the probe Gate opened, but
one component contains 86.43% of groups; all inference is component-macro.

The five-fold outer split holds out complete components. Standardization and
ridge selection occur inside each outer training fold. The planted linear
positive control passes (`R2=0.999999`, group Pearson `0.999999`). The real
result is negative:

| Frozen representation | Component-macro MSE | Group Pearson | Pair sign accuracy |
|---|---:|---:|---:|
| zero / ligand-only | 0.585 | NA | 0.500 |
| nuisance length/composition | 0.839 | -0.035 | 0.494 |
| ESM pooled additive | 0.843 | 0.244 | 0.622 |
| calibrated T-BASIS | 0.926 | -0.180 | 0.407 |
| rewired T-BASIS coupling null | 0.635 | -0.028 | 0.482 |

The preregistered loss reductions are both negative:

```text
rewired coupling null minus T-BASIS  -0.2909, one-sided 95% LCB -0.5220
ESM additive minus T-BASIS           -0.0833, one-sided 95% LCB -0.4519
999 fixed-hyperparameter diagnostics uncalibrated p = 1.0
```

The negative direction is robust to weighting: 15/21 components favor the
controls; group-weighted control-minus-T-BASIS differences are `-2.046` for the
single rewiring, `-2.229` for ESM-additive and `-2.562` for zero. Removing the
giant component leaves the rewiring and zero contrasts negative.

An independent implementation audit found that this is not a publication-ready
coupling-specific randomization test. The single rewiring leaves 1,870/4,726
rows unchanged and is identity for 719/1,820 groups. The label permutation does
not preserve repeated-family incidence, and its ridge hyperparameters were
selected on the observed labels. Its `p=1.0` is retained only as an uncalibrated
diagnostic, not confirmatory evidence. These defects cannot turn the observed
negative directions into a PASS; they do forbid a calibrated causal statement
about coupling.

This supports only the bounded fail-closed statement that the current audit did
not identify selectivity information recoverable from the 288D representation
by the source-only, component-held-out linear probe. It does not prove zero
mutual information, exclude a weaker/nonlinear signal, localize where
information was lost, or establish that local interaction states are absent.

## A2 and next action

A2 was conditional on an A1 PASS and was not run. Inspecting raw aggregates or
local atom-residue states now would be a post-outcome expansion of the
hypothesis family.

The next admissible action is data acquisition, not architecture search:
prospectively construct an independent assay-matched dense crossed-selectivity
cohort in which the same ligands have continuous measurements across related
targets. Preserve assay conditions, endpoint, unit, qualifier, replicate,
document, scaffold and protein-family provenance. Freeze dependency closure and
effect thresholds before labels are used. Only a fresh confirmatory A1 PASS can
authorize A2 and a separately registered low-capacity partner anchor.

## Artifacts

- `research/meta_fewshot/PREREG_BIOLOGICAL_GAUGE_AUDIT.md`
- `report/meta_fewshot/a0_gauge_audit.json`
- `report/meta_fewshot/a1_selectivity_dependency_audit.json`
- `report/meta_fewshot/a1_selectivity_probe.json`
