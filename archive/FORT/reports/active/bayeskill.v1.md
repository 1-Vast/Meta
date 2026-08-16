# Wave 1 Bayesian Kill Test

Date: 2026-07-30  
Verdict: `STOP_W1_BAYES_BELOW_LIGAND_AND_GRADIENT`

## Protocol

- Endpoint: pKi only.
- Seed: 1729.
- Support: k=5.
- Train/development episodes: 64/64; 61 development homology components.
- Shared label-free query-span selector for train and development.
- Shared train-only B0 for every arm.
- Pair ambient dimension: 8; protein-conditioned subspace rank: 2.
- One training epoch after three B0 epochs.
- No flexible kernel, graph ligand encoder, structure, confirmation, or sealed data.

Before this run, the implementation fixed the rank bound, evaluation gate,
Bayesian-model-averaging variance, calibration variance, raw/applied correction
logging, support-label bypass, pair-feature normalization, landmark masked
mean, 1034-wide ligand validation, and latent/source observation boundary.
The resulting Wave 0 suite passed 21 CUDA tests.

## Development Metrics

| Arm | RMSE | MAE | Spearman | Pairwise | NLL |
| --- | ---: | ---: | ---: | ---: | ---: |
| B0 | 7.163 | 7.052 | -0.069 | 0.477 | 3.406 |
| Calibration | 1.736 | 1.500 | -0.069 | 0.477 | 2.915 |
| Ligand-only Bayesian | 1.324 | 1.070 | 0.098 | 0.538 | 2.133 |
| Gradient baseline | 1.247 | 1.012 | 0.073 | 0.527 | not modeled |
| Protein-conditioned Bayesian | 1.693 | 1.461 | -0.043 | 0.486 | 2.915 |
| Wrong support | 2.212 | 1.977 | -0.064 | 0.478 | 2.943 |
| Permuted labels | 1.792 | 1.552 | -0.099 | 0.465 | 2.921 |
| Protein-free basis | 1.737 | 1.507 | -0.069 | 0.477 | 2.918 |

Positive bootstrap gain means the protein-conditioned Bayesian arm is better
than the named control.

| Control | RMSE gain, 95% component bootstrap interval | Spearman gain, 95% interval |
| --- | ---: | ---: |
| Calibration | 0.05 [-0.03, 0.14] | 0.03 [0.01, 0.05] |
| Ligand-only | -0.36 [-0.46, -0.26] | -0.13 [-0.22, -0.04] |
| Gradient | -0.43 [-0.57, -0.26] | -0.09 [-0.19, 0.02] |
| Wrong support | 0.50 [0.26, 0.77] | 0.02 [-0.01, 0.06] |
| Permuted labels | 0.11 [-0.02, 0.28] | 0.07 [0.02, 0.14] |
| Protein-free basis | 0.05 [-0.03, 0.13] | 0.03 [0.01, 0.05] |

The correct-support destruction behaves in the expected direction, and true
labels improve ranking over permuted labels. These attribution checks do not
rescue failure against the two predictive controls. Family dominance is not
claimed because the frozen registry has no family field; homology-component
dominance is recorded in the JSON result.

## Posterior And Uncertainty

Ranking log-BF had mean 0.433, median -0.028, range [-1.157, 12.851]. Ranking
inclusion probability had mean 0.101 and median 0.049. Most episodes therefore
received little ranking mass, with the positive mean driven by a small tail.

Bayesian 50/80/95% coverage was 0.996/1.000/1.000. The corresponding mean
interval widths were 9.63/18.30/27.98. The model is severely overdispersed;
these intervals are not evidence of calibrated uncertainty.

## Compute

- Device: NVIDIA GeForce RTX 4060 Laptop GPU.
- Wall time: 20.72 seconds.
- Mean/peak utilization: 23.05%/41%.
- Mean/peak power: 12.94/24.38 W.
- Peak NVIDIA memory: 5613 MiB.
- Peak PyTorch allocated memory: 2245.8 MiB.

## Decision

Wave 1 fails because the Bayesian arm is materially worse than ligand-only and
gradient baselines in both absolute error and reordering. The calibration gain
is small and its RMSE interval includes zero. Uncertainty is also too broad.

Do not proceed to a flexible kernel, additional seeds, rank changes, wider
encoders, extra epochs, graph ligands, query-selector changes, or confirmation
as a rescue. Any reopening requires a separately preregistered diagnosis that
changes the identified failure mechanism rather than tuning this stopped arm.

Machine-readable result: [`bayeskill.v1.json`](bayeskill.v1.json).
