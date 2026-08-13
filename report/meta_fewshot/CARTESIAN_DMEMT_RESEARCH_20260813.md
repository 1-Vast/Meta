# Cartesian mechanism meta-learning: research, implementation, and decision record

## Decision

The old HyperSAR chain is replaced as the active research candidate by a
**Difference-constrained Mechanism-Evidence Meta-Transformer (D-MEMT)**.  The
mathematical theory is now inspiration and a feasibility/consistency audit, not
a frozen implementation specification.

This decision followed three independent read-only research tracks covering
Cartesian tensor networks, meta-learning, and data/evaluation governance, four
provided design memos, primary literature, current code, and the actual tensor
banks.  Implementation proceeded only after module-level falsification tests.

## Data finding that bounds the 3D claim

The governed BindingDB Ki task contains 499 targets, 9,880 ligands, and 17,717
cells.  Its active protein bank contains ESM2 pooled/residue tensors and masks;
its ligand bank contains 2D atom/bond graphs and masks.  Neither contains
Cartesian coordinates or a pair-specific pose.

The repository's BioLiP2/RCSB corpus contains 14,906 governed holo structures,
but the compiled supervision shards contain contacts/distance bins rather than
coordinates, and exact overlap with the BindingDB task is negligible (15 exact
protein sequences, about 205 exact ligands, two exact protein--ligand cells).
It is therefore an independent structure-supervision source, not a 3D sidecar
for the main DTA task.

Consequences:

- the active BindingDB run must use the sequence+2D fallback;
- contact/distance labels cannot be presented as Cartesian inputs;
- independent protein and ligand frames cannot create intermolecular vectors;
- only a verified common-frame complex may enable cross Cartesian edges;
- Cartesian performance is not claimed until a split-safe, hashed mapping and
  coverage manifest exist.

## Architecture

### Unified mechanism-slot contract

`BipartitePairSectionFormer` now retains all aligned latent interaction slots
before mean pooling.  Projecting their exact mean through the original heads
reconstructs the previous endpoint and section exactly.  These slots provide
the coordinate-free mechanism representation.

### D-MEMT core

For each support observation, the model computes a zero-shot residual and a
slot sensitivity.  Their product forms label-bound gradient-shaped evidence.
A support-order-free Transformer produces per-slot task prompts.  Query
transport receives aligned support--query slot differences, and scalar prompt
gates modulate the slots.  A difference-only odd reference kernel supplies an
affinity-anchor path.  k=1 is no longer structurally zero; k=0 returns exactly
the shared zero-shot prediction.

This is gradient-shaped evidence, not full-parameter MAML and not a solver.
There is no ridge, Cholesky, pseudoinverse, or deployment optimization.

### Optional Cartesian encoder

The sparse encoder maintains scalar, polar-vector, and symmetric-traceless
rank-2 states.  Edge displacement unit vectors create vector and rank-2 bases;
radial and chemical networks produce scalar coefficients.  Only scalar channel
mixing/gating is used, preserving O(3) transformation laws.  Invariant norms
are pooled into the same mechanism-slot contract used by D-MEMT.

The current main training path does not instantiate this branch because it has
no legal coordinates.  This is a deliberate general missing-geometry design,
not an implicit 3D claim.

The optional training/evaluation control path is nevertheless executable:
`EpisodeBatch` can carry common-frame coordinates, packed edges, availability,
and an explicit common-frame declaration.  The model rejects undeclared or
independently framed cross interaction and rejects every edge whose endpoints
belong to different flattened samples.  Thus future governed coordinate
sidecars can train the same D-MEMT path without weakening episode isolation.

## Module tests completed before training

- retained-slot endpoint/section exact reconstruction;
- support joint-permutation invariance;
- label-only permutation changes mechanism output;
- k=0 exact zero-shot fallback and truly empty support materialization;
- k=1 non-zero evidence, non-zero adaptation, and label gradients;
- task-state shape checks and override sensitivity;
- Cartesian O(3) equivariance including reflection;
- translation invariance, node permutation, padding invariance;
- symmetric and trace-free rank-2 state;
- missing-coordinate vector/tensor channels exactly zero;
- finite gradients through Cartesian and meta paths;
- common-frame Cartesian slots remain invariant after a joint rigid transform.
- Cartesian packed edges cannot cross support/query or batch samples;
- training-wrapper geometry fields reach the Cartesian parameters and receive
  finite gradients;
- undeclared and explicitly independent frames are rejected by the joint
  Cartesian interaction API.

## Sequential candidate results

All results below are development diagnostics, not confirmatory claims.

| Candidate | steps | best val MSE | test full MSE | SAR-cut MSE | SAR-cut − full | permutation − full |
|---|---:|---:|---:|---:|---:|---:|
| prompt scalar slot gate | 80 | 1.2046 | 1.5492 | 1.5492 | +0.00001 | −0.00041 |
| + difference-only reference anchor | 80 | 1.2358 | 1.5945 | 1.5741 | −0.02040 | −0.04094 |
| + support LOO relative supervision | 80 | 1.1523 | 1.5408 | 1.5327 | −0.00808 | −0.01340 |
| prompt-coupled relative kernel | 120 | 1.1629 | 1.5488 | 1.5422 | −0.00661 | −0.01515 |

The engineering hypothesis passed; the performance hypothesis has not yet
passed.  Stronger reference paths create visible corrections, but on the
development test they worsen MSE and correct label binding is not superior to
permutation.  This is evidence against claiming that the current mechanism is
an important performance source.

## Frozen full-evaluation protocol

The final development evaluation uses one model family, three fixed seeds,
all eligible targets, nested common-query k={0,1,2,3,5}, component-first
aggregation, and paired component bootstrap.  It reports MSE/RMSE, CI,
Spearman, label permutation, foreign mechanism prompts, wrong-protein state,
SAR-cut, level-only, and zero-shot.  k=0 uses a genuinely empty support set.

The historical meta-test has been repeatedly observed and is development-only.
No publication or SOTA claim is authorized without a newly sealed confirmation
split or an external cold-target benchmark.

## Completed three-seed nested evaluation

The frozen run trained three seeds for 120 steps and evaluated 42 targets from
six CD-HIT40 components with a single nested common-query manifest.  Peak GPU
memory was 3.96 GiB.  The table reports equal-component/equal-target means.

| k | full MSE | level-only MSE | SAR-cut MSE | CI | Spearman | SAR-cut − full (95% component CI) |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 3.3371 | 3.3371 | 3.3371 | 0.5261 | 0.0656 | 0.0000 [0.0000, 0.0000] |
| 1 | 2.7003 | 2.7411 | 2.7031 | 0.5246 | 0.0730 | 0.0028 [−0.0086, 0.0170] |
| 2 | 2.3018 | 2.3802 | 2.3100 | 0.5392 | 0.1130 | 0.0082 [−0.0049, 0.0260] |
| 3 | 1.8314 | 1.9084 | 1.8482 | 0.5457 | 0.1227 | 0.0168 [0.00002, 0.0405] |
| 5 | 1.5643 | 1.6433 | 1.5915 | 0.5559 | 0.1598 | 0.0272 [0.0064, 0.0600] |

At k=3, correct support also beats label permutation by 0.0166 pK²
(95% component CI [0.0027, 0.0309]).  At k=2 this contrast is positive but
uncertain; at k=5 its lower bound is slightly negative.  k=1 cannot have a
nontrivial permutation control because one label has no alternate binding.

Interpretation is deliberately narrow:

- D-MEMT supplies a statistically positive marginal MSE contribution at k=3
  and k=5 on this consumed development population;
- the effect is modest (about 0.9% and 1.7% relative to SAR-cut), below the
  preregistered ~5% target;
- k=1 and k=2 are inconclusive, so the requirement of excellent performance
  throughout k=1/2/3/5 is not met;
- absolute ranking remains weak (CI 0.556 and Spearman 0.160 at k=5);
- foreign/wrong-protein *prompt-only* gaps are near zero because those controls
  replace only task prompts, not the complete reference memory and level state;
- the run supports continued development, not an outstanding/SOTA claim.

Artifacts are in `report/meta_fewshot/dmemt_nested_k01235_3seed_20260813/`:
`EPISODES.json`, `PREDICTIONS.jsonl`, and `RESULT.json`.

The final repository-wide regression run completed with 272 passing tests;
the independent final review found no remaining P0/P1 issue after the packed
edge and training-control-path fixes.

## Primary sources

- PBCNet2.0: https://doi.org/10.1038/s41589-026-02241-x
- official PBCNet2.0 code: https://github.com/YuJie-0202/PBCNet2.0
- TensorNet: https://proceedings.neurips.cc/paper_files/paper/2023/hash/75c2ec5f98d7b2f50ad68033d2c07086-Abstract-Conference.html
- AdaMBind: https://www.nature.com/articles/s41467-026-70554-5
- Conditional Neural Processes: https://proceedings.mlr.press/v80/garnelo18a
- Set Transformer: https://proceedings.mlr.press/v97/lee19d.html
- CrossTransformers: https://proceedings.neurips.cc/paper/2020/hash/fa28c6cdf8dd6f41a657c3d7caa5c709-Abstract.html
