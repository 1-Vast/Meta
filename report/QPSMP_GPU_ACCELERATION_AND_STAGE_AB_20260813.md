# QPSMP GPU Acceleration and Stage A/B Development Results

## Scope

This report records a consumed-development engineering run. It does not
authorize G2, G3, biological-mechanism, or external-generalization claims.

## Training-system correction

The original runner decompressed padded ligand-graph shards and forwarded one
episode at a time. On the RTX 4060 system this produced disk saturation and low
GPU utilization. The corrected path:

- precomputes a label-free compact ligand bank using real atom counts;
- caches compact episodes before optimization;
- batches protein, ligand, support, and query tensors across episodes;
- executes pooled and atom-residue interaction modes on CUDA; and
- writes validation progress during training.

The compact bank contains 9,880 ligand graphs and occupies approximately 3 MB.
It was generated in 13.5 seconds. The batched implementation matches the
single-episode implementation to floating-point tolerance in both interaction
modes.

The 300-step, 32-episode-per-step training phase plus validation required about
31 seconds. The previous serial 300-step attempt ran for more than 85 minutes
without producing a result. Sampled atom-residue training reached about 46%
GPU utilization and 26 W while disk activity was no longer saturated. Final
control evaluation remains less parallel than training and dominates total
wall time.

## Stage A: zero-support representation

Both arms used the same seed, budget, support/query protocol, and consumed
development cohort (18 test episodes).

| Arm | Validation MSE | Test zero-shot MSE | Ligand-only MSE | Wrong-protein zero-shot MSE | Wrong-protein gap |
|---|---:|---:|---:|---:|---:|
| Pooled | 1.2046 | 3.0108 | 2.6114 | 3.0857 | +0.0750 |
| Atom-residue | 1.4853 | 2.5541 | 3.1939 | 3.7873 | +1.2333 |

The atom-residue arm improved the test point estimate and passed the tested
point-estimate controls, but its validation MSE was worse than the pooled arm.
With one seed, 18 episodes, and no component-bootstrap lower bound, this is a
development signal rather than an admitted G2 result.

## Stage B: learned support-span posterior

The atom-residue representation was trained for 500 steps with the learned
support-span posterior.

| Metric | MSE (pK squared) |
|---|---:|
| Full | 1.2767 |
| Level-only | 1.4325 |
| SAR-cut | 1.2609 |
| Permuted state | 1.3129 |
| Foreign state | 1.2700 |
| Wrong-protein state | 1.2875 |

The full endpoint improved over level-only by 0.1558, but the isolated SAR
contribution was harmful (`SAR-cut - full = -0.0157`) and the foreign-state arm
was better than correct state. Therefore the proposed learned meta-posterior is
trainable but is not admitted as a target-specific few-shot SAR mechanism.

## Decision

The GPU batching, compact-bank preprocessing, and atom-residue representation
remain in the core implementation. The learned support-span posterior remains
an experimental candidate behind an explicit mode selector; it is not treated
as validated performance evidence. Further model expansion stops until a
multi-seed component-paired Stage A confirmation and a positive SAR-cut/control
result are available.
