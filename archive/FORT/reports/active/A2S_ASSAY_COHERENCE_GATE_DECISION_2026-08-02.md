# A2S Assay-Coherence Information Gate Decision

Date: 2026-08-02
Status: source-only stop decision
Branch at execution: `research/a2s-assay-coherent-20260802`

## Question

**HYPOTHESIS.** Affinity is a conditional thermodynamic observable. Restricting
support and query measurements to the same exact ChEMBL assay may remove enough
construct/readout noise for a k-shot target adaptation state to become
identifiable in chemically distant queries.

This was an information gate, not a meta-learning model. The tested empirical
Bayes and KRR estimators are baselines only.

## Preregistered Admission Rule

At both `k=3` and `k=5`, in query strata with nearest-support Tanimoto below
0.35, at least one of the two global response bases had to satisfy all of:

1. CI gain over the frozen base: component-bootstrap lower 95% bound `> 0.005`.
2. Correct residual assignment over residual permutation: lower bound `> 0`.
3. Correct target support over wrong-target support: lower bound `> 0`.
4. A synthetic two-state positive control had to pass.

Thresholds and methods were not changed after execution.

## Mechanical And Firewall Checks

**FACT.** `D:/anaconda/envs/drug/python.exe -m pytest
tests/test_a2s_assay_coherence_gate.py -q` returned `6 passed`.

**FACT.** The run opened only source roles `fit` and `probe`:

- fit: 5,928 episodes, 143 targets;
- probe: 2,488 episodes, 74 targets, 68 homology components;
- fit/probe target overlap: 0;
- locked labels requested: false;
- recipient labels requested: false;
- CUDA device: NVIDIA GeForce RTX 4060 Laptop GPU.

## Results

**FACT.** Exact-assay grouping reduced residual dispersion but did not meet the
adaptation-information gate:

- median within-assay residual SD: 0.8275 pKi;
- median within-target residual SD: 1.2647 pKi;
- ratio: 0.6543;
- synthetic k=3 correct-vs-wrong CI: 0.3124, lower 95% 0.2290;
- synthetic k=5 correct-vs-wrong CI: 0.3963, lower 95% 0.3116;
- synthetic verdict: pass.

The decisive real-data cells were:

| k | Tanimoto | Components | Basis | CI gain mean | Gain lower 95% | Assignment lower 95% | Wrong-target lower 95% | Pass |
|---|---|---:|---|---:|---:|---:|---:|---|
| 3 | 0.00-0.20 | 5 | descriptors | +0.0022 | -0.0143 | +0.0072 | -0.0125 | no |
| 3 | 0.20-0.35 | 12 | descriptors | -0.0069 | -0.0250 | -0.0307 | -0.0350 | no |
| 3 | 0.00-0.20 | 5 | original | -0.0020 | -0.0184 | -0.0040 | -0.0001 | no |
| 3 | 0.20-0.35 | 12 | original | -0.0071 | -0.0277 | -0.0186 | -0.0274 | no |
| 5 | 0.00-0.20 | 4 | descriptors | +0.0029 | 0.0000 | -0.0230 | +0.0022 | no |
| 5 | 0.20-0.35 | 5 | descriptors | +0.0309 | +0.0015 | -0.0001 | -0.0326 | no |
| 5 | 0.00-0.20 | 4 | original | +0.0032 | -0.0476 | -0.0818 | -0.0025 | no |
| 5 | 0.20-0.35 | 5 | original | +0.0643 | +0.0223 | +0.0172 | +0.0378 | yes |

**FACT.** The machine verdict was
`ASSAY_COHERENT_GLOBAL_ADAPTATION_INFORMATION_NOT_ADMITTED`: k=3 failed,
k=5 passed only one cell, and the synthetic control passed.

**FACT.** In chemically local/all-query cells, scaled KRR remained stronger.
For example, at k=5 and Tanimoto >=0.55, scaled KRR CI gain was +0.0534
(lower 95% +0.0374), versus +0.0213 for the original-basis EB adapter.

## Interpretation

**INFERENCE.** Assay context is a real noise variable, but noise reduction is
not the missing transferable target state. The k=3 distant-query effects were
null or negative and generally failed assignment and wrong-target controls.

**INFERENCE.** The isolated k=5 success is not generalizable evidence. It is
supported by only five homology components, does not replicate at k=3, and is
surrounded by null/negative adjacent cells.

**INFERENCE.** Randomly adding more support compounds also raises each query's
nearest-support similarity. Consequently, the strict low-similarity cells have
only 4-12 components even though the exact-assay episode census is large. This
is a power limitation of the distant stratum, but it cannot be repaired by
tuning the failed estimator or relaxing the gate after seeing the result.

## Decision

1. Stop exact-assay EB/global-basis adaptation as a primary route.
2. Keep assay metadata as a conditioning/censoring variable and stratification
   control, not as the meta-adaptation object.
3. Keep scaled KRR/TRACE as the chemically local baseline.
4. Move to a new branch based on a different biological principle: relative
   free-energy responses to defined chemical perturbations and thermodynamic
   cycle consistency.
5. Do not open locked or recipient labels.

## Artifacts

- JSON content hash: `c1917620a8c71dd29adbb18a3c57bb18a579dbffc63171f75f2a29530b590c0d`
- JSON file SHA-256: `d75e6c980185a5eade634b0cb1506860682bd567e6ca44ced52bdff54a4952d3`
- Records SHA-256: `000337fce3f733e0fa5e76f4c0b408bbf8fc449c94a5074e39de7b68763de87c`

