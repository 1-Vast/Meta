# DCST two-stage preregistration

Date: 2026-07-28  
Status: frozen before any DCST development score  
Runner: `research/dcst_two_stage.py`  
Model: `model/dcst.py`

## 1. Question and scope

The user has explicitly selected a two-stage model program:

1. learn or adapt an interaction representation on a constructed high-quality
   source dataset;
2. train the final predictor on the strict dual-cold affinity dataset.

The candidate innovation is **DCST — Destruction-Certified Spectral
Transfer**. Ordinary pretraining followed by fine-tuning is a baseline, not
the novelty claim. DCST asks whether the source interaction can be decomposed
into pair-dependent directions, whether those directions survive held-source
destruction tests, and whether transferring only the certified directions
prevents negative transfer while a separate Stage-2 residual learns new
interactions.

This preregistration supersedes `OMUT-X7` as the active executable program.
It does not invalidate the historical OpenMut evidence or convert mutation
geometry into a model input.

## 2. Frozen assumptions

- Stage 1 source: PLINDER 2024-06/v2 processed dual-cold
  train/development registry, using exact BindingMOAD-derived affinity,
  Morgan-plus-descriptor ligand features, and frozen ESM-2 pooled plus eight
  ordered segment features.
- Stage 2 target: ChEMBL-37 exact Ki/Kd registry, train for optimization and
  development for scoring.
- The locally reported FIRE-DTA atom caches are not present in the current
  workspace. DCST is therefore a sequence/ligand-feature model, not a 3D
  structure model.
- ChEMBL confirmation and sealed partitions are not loaded by the Stage-2
  dataset object and are never scored.
- Endpoint, source, assay, document, target-only, and ligand-only observation
  heads do not transfer.
- The PLINDER processed registry has no DOI/document field compatible with
  ChEMBL documents. Cross-source document-lineage closure is therefore an
  explicit unresolved limitation, not a passed gate.

## 3. Cross-source firewall

Before ChEMBL development affinity is loaded, source train and source
development rows are projected against ChEMBL development plus confirmation
metadata. A source row is excluded by any of:

- exact UniProt accession overlap;
- protein 4-mer containment at least 0.40 against the shorter sequence;
- exact ligand parent connectivity;
- exact Bemis–Murcko scaffold;
- Morgan radius-2 Tanimoto at least 0.95.

The label-blind projection observed 42,539 protected downstream metadata
rows. From 6,821 source train/development candidates it retains 5,130 train
rows and 447 development rows; 1,244 rows are excluded by the union. Counts by
axis are accession 829, homology 898, parent 185, scaffold 572, and chemical
near-neighbour 478. Overlaps are not additive.

Every Stage-1 normalizer, base model, cross-fitted residual and interaction
model is rebuilt on the retained Stage-1 train rows only.

## 4. Frozen model

### Stage 1

The deployable interaction branch contains:

- a compact ordered-segment encoder over frozen ESM-2 features;
- a ligand MLP over the common 1,034-dimensional ligand contract;
- a direct bilinear interaction `g(t,d)=z_t^T Theta z_d` with exact regular
  null `Theta=0`.

A separate ligand-only base is five-fold cross-fitted on target/homology
components. The interaction branch is trained on the cross-fitted residual
with rank-reversal, residual MSE, within-target centering, and base-output
orthogonality. The ligand-only base is never transferred.

### Held-source certificate

`Theta=U diag(s) V^T` is truncated to rank 8. For each component, source
development target-macro Spearman against the cross-fitted-base residual is
compared with:

- a deranged source target representation;
- a deterministic within-target ligand derangement.

For component `j`,

```text
margin_j = utility_true_j
           - max(abs(utility_target_destroyed_j),
                 abs(utility_ligand_destroyed_j)).
```

A component is inactive unless both `utility_true_j>0` and `margin_j>0`.
Positive margins map to confidence
`1-exp(-margin_j/0.05)`. No ChEMBL development value selects a component or
certificate threshold.

### Stage 2

The certified source encoders and spectral components are frozen. Each active
component receives one global gate initialized from its source confidence.
A separate interaction residual branch is initialized with the Stage-1
encoders and exact-zero bilinear operator, then trained on ChEMBL train. It
cannot overwrite the frozen source path. Its output is penalized for
within-target correlation with the transferred output. A weak gate-prior
penalty permits a harmful source direction to close.

The final standardized score is:

```text
b_ChEMBL(d)
+ sum_j gate_j * certified_source_component_j(t,d)
+ g_stage2_residual(t,d).
```

## 5. Frozen comparisons

All trained Stage-2 arms receive the same target-balanced optimizer steps:

| Arm | Meaning |
| --- | --- |
| `B0` | ChEMBL ligand-only base |
| `Scratch` | same interaction branch trained from zero on ChEMBL |
| `NaiveFT` | complete Stage-1 branch conventionally fine-tuned |
| `FrozenEncoderFT` | Stage-1 encoders frozen, full interaction operator tuned |
| `FullTransferResidual` | all top-eight source directions plus the same residual branch |
| `DCST-CertShuffle` | DCST with source certificate confidences permuted |
| `DCST` | destruction-certified spectral transfer |
| `DCST-Tshuffle` | inference-time target derangement |
| `DCST-Lshuffle` | interaction-branch ligand derangement; base unchanged |

The default accepted run uses seed 1729, 4,000 base steps per fold, 4,000
Stage-1 interaction steps, 4,000 Stage-2 steps per arm, width 32, rank 8, Adam
at `2e-3`, and the frozen loss weights in the runner. Shorter runs are
engineering smoke tests only.

## 6. Metrics and pass gates

Primary metric: ChEMBL development target-macro Spearman.  
Independent bootstrap unit: target homology component.  
Secondary: target-macro RMSE/MAE/concordance and negative-transfer rate.

The pre-existing zero-shot MDE is 0.0586. A scientific pass requires all:

1. at least one source direction receives a positive destruction certificate;
2. `DCST-B0 >= 0.0586`;
3. grouped bootstrap LCB95 for `DCST-B0` is positive;
4. grouped LCB95 is positive for both `DCST-Scratch` and `DCST-NaiveFT`;
5. DCST target-macro RMSE is no more than 2% worse than B0;
6. both target and ligand destruction remove at least 70% of the incremental
   Spearman gain.

Failure of the complete gate does not authorize rank, scale, width, epoch, or
certificate-threshold tuning on ChEMBL development. Diagnosis must use the
registered arms to decide whether the failure is source information,
cross-stage transfer, downstream residual learning, or shortcut dependence.

## 7. Claim boundary and current literature

Recent DTA work already uses pretrained molecular/protein encoders, ordinary
feature fusion, and two-stage latent objectives. Co-Diffusion freezes a
Stage-1 affinity-aligned latent space before Stage-2 diffusion; DTIAM and
ColdstartMHDTI use independently pretrained encoders; TAPB uses target
randomization for debiasing. DCST does not claim novelty for pretraining,
bilinear interaction, SVD, gating, or orthogonality separately. Its candidate
novelty is the combination of:

1. held-source destructive certification of individual interaction spectral
   directions;
2. selective cross-dataset transfer of only those directions;
3. a protected frozen source path plus a separately trained downstream
   residual that can add information without catastrophic forgetting;
4. negative-transfer and destruction gates under a simultaneous target- and
   ligand-cold split.

This remains a candidate novelty until the registered controls show that the
certificate, rather than extra capacity or initialization, is load-bearing.

Primary references inspected before freezing:

- https://arxiv.org/abs/2603.11125
- https://www.nature.com/articles/s41467-025-57828-0
- https://www.nature.com/articles/s41467-025-66915-1
- https://doi.org/10.3389/fchem.2026.1846850

