# MAML Stability Audit

## Evidence

The fixed-seed pKi E0 run used 64 train targets and 64 strict development
episodes across 61 homology components. The complete gradient adapter scored
RMSE 1.3088, MAE 1.0420, and within-target Spearman 0.0719. Frozen ligand-only
B0 plus support calibration scored 1.2902, 1.0288, and 0.0841. The adapter is
therefore worse by 0.0185 RMSE, 0.0132 MAE, and 0.0122 Spearman.

## Causes

1. The outer optimizer updates one target at a time, in lexical target order;
   64 tasks and one epoch are high-variance rather than a meta-batch estimate.
2. Training uses first-sorted-scaffold supports while evaluation uses the
   frozen query-span roster. This is a train/evaluation episode mismatch.
3. Zero-initializing U prevents random residual blow-up, but makes its first
   update a readout-only bootstrap; the support and protein paths learn later.
4. The MAML inner loop has no posterior covariance, evidence gate, or
   abstention. k=5 optimization noise is not separated from task variation.
5. SetEncoder consumes ligand, B0, and residual only, so it can encode a
   ligand shortcut instead of a protein-specific adaptation state.
6. The active backbone uses a single forward Mamba scan, not TrueBiMamba.

## GPU Observation

The E0 process used 5.71 GiB peak memory but only 16.7% mean and 44% peak
CUDA utilization at 14.10/17.51 W mean/peak draw. This is consistent with
small sequential target kernels and CPU orchestration. A Windows 3D graph near
89% is a different engine counter and cannot establish CUDA training use.

## Decision

Do not increase epochs, width, or seeds for this primary adapter. Retain it as
a matched MAML baseline. Replace meta-test gradients with the
registered Bayesian posterior before the next training run.
