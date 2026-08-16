# A2S-DTA Alignment Correction And Decision

Date: 2026-07-31  
Status: corrected closed-form evidence; neural architecture remains blocked

## P0 correction

The first A2S runs are invalid and must not be cited:

- `reports/active/a2s_pki_smoke_seed1729.json`
- `reports/active/a2s_pki_seed1729.json`
- `reports/active/a2s_pkd_seed1729.json`

The loader filtered `registry.parquet` before assigning `source_row`.  That
created a local index while `ligand_features.npz` is aligned to all 343,211
registry rows.  The corrected loader reads the complete registry, resets the
global row id, verifies feature length and `conn_sha`, and only then filters to
TRAIN and the requested endpoint.  The regression suite has 72 passing tests.

## Corrected execution

All runs used CUDA in the `drug` environment, TRAIN rows only, source targets
with `n_eff >= 100`, recipient candidates with `n_eff < 30`, and target-side
single-cold episodes.  Source and recipient target IDs remain disjoint.  The
primary corrected fit gives every source target equal total weight while
preserving the original total source weight
(`target_macro_equal_total_weight`).  The earlier `a2s_*_corrected_seed1729.json`
artifacts are alignment-valid row-weighted diagnostics; the target-balanced
artifacts below supersede them for the primary result.

| Endpoint | Source targets | Recipient candidates | k=1 | k=3 | k=5 | Artifact |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| pKi | 242 | 193 | 151 | 137 | 117 | `a2s_pki_targetbalanced_seed1729.json` |
| pKd | 41 | 256 | 119 | 101 | 87 | `a2s_pkd_targetbalanced_seed1729.json` |

Corrected source-support RMSE gain versus recipient calibration:

| Endpoint | k=1 | k=3 | k=5 |
| --- | ---: | ---: | ---: |
| pKi | +0.137 [0.070, 0.207] | +0.112 [0.054, 0.174] | +0.087 [0.029, 0.149] |
| pKd (secondary) | +0.373 [0.262, 0.497] | +0.220 [0.140, 0.312] | +0.242 [0.149, 0.359] |

Intervals are target-level bootstrap 95% intervals. pKi and pKd are not pooled.
These are descriptive source-support controls: source selection uses declared
support labels and is not the proposed held-out router.

## Cross-fitted router decision

The target-balanced source-fold router used 726 pseudo-recipient episodes and
117,126 candidate rows; natural recipient query labels were not used for
fitting.  Corrected pKi router gain versus recipient calibration was:

| k | Router | Router-gated | Random source |
| --- | ---: | ---: | ---: |
| 1 | -1.144 [-1.378, -0.914] | -1.144 [-1.371, -0.906] | -0.121 [-0.259, 0.029] |
| 3 | -1.213 [-1.443, -0.979] | -1.213 [-1.455, -0.977] | -0.420 [-0.580, -0.272] |
| 5 | -1.315 [-1.601, -1.050] | -1.315 [-1.586, -1.056] | -0.466 [-0.616, -0.320] |

The registered router route is **NO-GO**: it fails the primary pKi k=3/5
positive-gain criterion and its gate does not abstain reliably.  The positive
source-support control cannot be promoted to a learned transfer claim because
it uses recipient support labels directly for source choice.

The artifact schema now also reports NDCG@10, benefiting-recipient rate, and
RMSE-gain AULC.  The target-balanced pKi router AULC is `-1.221`, and its mean
benefiting-recipient rate is `0.201`; the stop is therefore not caused by
omitting ranking or learning-curve metrics.

## Decision and next gate

1. Accept the target-balanced artifacts as the primary corrected A2S evidence;
   retain row-weighted corrected outputs as alignment diagnostics only.
2. Keep pKi primary and pKd secondary; do not pool endpoints.
3. Do not train Mamba/Transformer, MAML, or a larger router.
4. Before another model attempt, freeze support draws, target-balanced source
   fitting, and a genuine natural-tail document/time/source-closed roster.
5. Re-open a small router only if those protocol corrections provide adequate
   pKi k=5 query/component power.  Pseudo-tail-only gains remain insufficient.
