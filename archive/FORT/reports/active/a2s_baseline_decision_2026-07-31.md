# A2S-DTA Baseline Decision (SUPERSEDED)

> This report summarizes pre-correction runs and is retained only for audit
> history.  Its metrics are invalid because filtered parquet indices were
> misaligned with the global ligand feature cache.  Use
> `a2s_alignment_correction_decision_2026-07-31.md` and the corrected JSON
> artifacts instead.

Date: 2026-07-31  
Protocol: target-disjoint A2S target-side single-cold  
Source rule: `n_eff >= 100`  
Recipient rule: `n_eff < 30`  
Seed: 1729

## Execution

The closed-form A2S runner used only TRAIN rows. Source-target labels fitted a
pooled 1,034-dimensional ridge and source-target affine adapters. Recipient
labels were used only in the declared support/query episode. The primary
support budgets were `k={1,3,5}`; query compounds were allowed to overlap
source chemistry, as required by the target-side single-cold primary track.

The runner compared B0, recipient calibration, support-compatible source
routing, random-source routing, protein-similarity routing,
chemistry-similarity routing, and a support-evidence abstention gate. No
Transformer/Mamba architecture was trained.

## Roster and compute

| Endpoint | Source targets | Recipient candidates | k=1 episodes | k=3 episodes | k=5 episodes |
| --- | ---: | ---: | ---: | ---: | ---: |
| pKi | 242 | 193 | 151 | 137 | 117 |
| pKd | 41 | 256 | 119 | 101 | 87 |

The pKi run used 31.9 seconds wall time and 48.6 MiB peak Torch memory on the
RTX 4060 Laptop GPU. The pKd run used 6.3 seconds and 48.5 MiB. These are
closed-form baseline costs, not a neural training budget.

## Results

Values are RMSE gain relative to recipient calibration; positive is better.
Intervals are target-level bootstrap 95% intervals.

| Endpoint / k | Support-compatible source | Abstention gate | Random source | Protein source | Chemistry source |
| --- | --- | --- | --- | --- | --- |
| pKi / 1 | `+0.034 [-0.004, +0.075]` | `0.000 [0.000, 0.000]` | `-0.374 [-0.558, -0.180]` | `-0.377 [-0.568, -0.209]` | `-0.350 [-0.521, -0.187]` |
| pKi / 3 | `-0.020 [-0.058, +0.012]` | `-0.033 [-0.058, -0.011]` | `-0.695 [-0.865, -0.524]` | `-0.627 [-0.788, -0.456]` | `-0.689 [-0.859, -0.518]` |
| pKi / 5 | `-0.013 [-0.044, +0.021]` | `-0.007 [-0.029, +0.015]` | `-0.696 [-0.878, -0.520]` | `-0.626 [-0.791, -0.471]` | `-0.744 [-0.918, -0.583]` |
| pKd / 1 | `+0.286 [+0.203, +0.389]` | `0.000 [0.000, 0.000]` | `+0.083 [-0.052, +0.229]` | `+0.155 [-0.004, +0.302]` | `+0.020 [-0.128, +0.177]` |
| pKd / 3 | `+0.163 [+0.101, +0.226]` | `+0.125 [+0.081, +0.171]` | `-0.198 [-0.328, -0.076]` | `-0.177 [-0.298, -0.053]` | `-0.299 [-0.423, -0.180]` |
| pKd / 5 | `+0.132 [+0.065, +0.195]` | `+0.151 [+0.101, +0.200]` | `-0.360 [-0.528, -0.209]` | `-0.235 [-0.374, -0.108]` | `-0.424 [-0.590, -0.268]` |

## Decision

1. **A2S data and evaluation chain: PASS.** The target-disjoint source-only
   protocol, primary pKi roster, support/query separation, CUDA execution, and
   target-level bootstrap all ran successfully.
2. **Universal scalar source-adapter claim: NO-GO.** It does not give a
   positive pKi transfer gain at the primary k=3/5 budgets. The pKi support
   router's negative-transfer rates were 43.0%, 48.2%, and 52.1% for k=1,3,5.
3. **Abstention as implemented: diagnostic only.** It reduced pKi harm at k=3/5
   but did not produce a positive gain. It cannot yet be claimed as the A2S
   innovation.
4. **pKd: secondary positive signal, not a pooled conclusion.** The pKd
   support router is positive at all three budgets, but endpoint differences
   require an independent replication and endpoint-specific mechanism audit.

The A2S task is not stopped. The next admissible model step is a
recipient-conditioned router trained only on source-target pseudo-recipient
episodes with cross-fitting, followed by a held-out natural pKi recipient
test. It must beat the current recipient-calibration baseline and the random
source control. No Mamba/Transformer capacity expansion is justified by this
baseline result.
