# T-BASIS-R0 Fixed Radial Basis Result

Date: 2026-08-07

## Verdict

```text
RADIAL_BASIS_PARTNER_RECOVERABILITY_IDENTIFIED
```

This is a structure-only, research-stage PASS. It establishes recovery of one
fixed two-body radial basis. It does not establish angular/many-body recovery,
affinity direction, universality, few-shot adaptation or admission to `z`.

## Governed Panel

- `192/64/64` train/validation/test complexes;
- 320 distinct homology groups, PDB IDs and exact protein sequences;
- all 40 T-DIR-P0 records excluded;
- validation/test are P1B-held-out splits;
- no cross-split selected scaffold overlap;
- score-blind wrong-protein map is one-to-one with reuse `0`;
- maximum wrong-protein sequence identity `0.32632`;
- affinity, DAVIS and recipient label reads `0`.

## Fixed Basis

The privileged teacher is a permutation-invariant `8 x 6 x 6 = 288`
chemogeometric tensor:

```text
ligand atom chemistry x residue chemistry x continuous Gaussian radial basis
```

It uses actual holo atom-to-slot distances only to construct the teacher. The
student uses frozen P1B five-bin distance probabilities and sequence/2D
chemistry. A single train-only shared `6 x 6` Ridge radial calibration maps
P1B bin moments into the fixed radial coordinates.

These coordinates are structural moments. They are not binding energy and do
not constitute a complete uncertainty law.

## Results

| Split | Mean MSE | Correct MSE | Deranged MSE | Reconstruction gain | Partner gain |
|---|---:|---:|---:|---:|---:|
| Train | 1.0000 | 0.5097 | 0.7462 | 0.4903 | 0.2365 |
| Validation | 1.1823 | 0.5529 | 0.7249 | 0.5324 | 0.1455 |
| Test | 1.1001 | 0.5157 | 0.6875 | 0.5312 | 0.1561 |

Test reconstruction-gain 95% CI is `[0.44328, 0.59621]`. Test partner-gain
95% CI is `[0.10695, 0.20070]`. All six preregistered conditions pass.

The uncalibrated P1B basis had test MSE `26.6533`, so the PASS depends on the
small shared radial calibration. This is expected because uniform within-bin
moments and the teacher's continuous radial functions have different scale and
shape; it also means the raw P1B probabilities are not themselves an admitted
basis.

## Interpretation Boundary

The result supports the first UMBD premise:

```text
fixed privileged 3D radial basis
        -> recoverable from sequence + 2D frozen P1B moments
        -> recovery degrades with a score-blind wrong protein
```

It does not yet show that angular or many-body information absent from the
five-bin output can be distilled. It also does not establish affinity value.

The wrong-protein arm changes both residue chemistry composition and predicted
distance distributions. Therefore `partner_gain` is a partner-conditioned
contrast, not a clean decomposition of pair geometry versus protein marginal.
A future basis stage must add marginal-preserving controls before claiming
pair-specific mechanism recovery.

## Authorization

This PASS authorizes only a separately preregistered structure-only study of
low-order angular/many-body privileged basis distillation with:

- a new sealed panel or a declared development/validation split;
- fixed analytic coordinates and gauge;
- correct, geometry-deranged and chemistry-marginal controls;
- local and global reconstruction reporting;
- abstention from affinity, DAVIS and production integration.

CDAC affinity calibration, rank-aware few-shot adaptation, biological `z`,
CSMO/Band changes and P2-P4 remain frozen.
