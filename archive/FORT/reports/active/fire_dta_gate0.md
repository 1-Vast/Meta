# FIRE-DTA Gate 0: Feasibility, Novelty and Engineering

## Verdict

`FIRE_DTA_GATE0_CONDITIONAL_PASS`

BFEO and BERP are algebraically coherent, independently testable and implementable on the local RTX
4060 while preserving the validated B1 posterior. They are not yet validated scientific innovations.
Long training remains locked until a leakage-safe external structural subset proves that the
bound-minus-free ensemble representation transfers beyond IPBind-style controls.

## 1. Task definition and biological objective

Predict ChEMBL-37 pKi for previously unseen targets at k=0/4/8/16. The deployment goal is reliable
compound prioritization with calibrated uncertainty and fewer harmful support updates. The independent
statistical unit is the target/sequence cluster, not an affinity row, pose, frame or compound pair.

## 2. Deployment scenario

At inference, protein sequence and ligand graph are universal. AlphaFold coordinates are nearly
universal but a validated pocket and bound pose are not. Candidate states may therefore be missing or
unreliable; FIRE-DTA must abstain or fall back rather than route around BFEO through a direct sequence or
ligand affinity head.

## 3. Current-model failure

The current char-CNN/Morgan/FiLM base lacks atom-residue and state information. B1 is mathematically
sound but adapts in a learned ligand latent. Prior residue fields, pocket distillation, pretrained
compatibility and structural priors learned auxiliary signals without improving the cold-target
estimand. The new hypothesis is permissible only because its supervision source is external complex,
apo/holo and trajectory evidence rather than scalar affinity alone.

## 4. Concentrated innovation modules

**BFEO is the required core module.** It uses one shared invariant atomistic encoder for complex, free protein and free ligand, computes
an exact channel contrast, and marginalizes multiple candidate states. The implementation has no direct
protein/ligand affinity bypass. IPBind makes single-state bound-minus-free a standard control, so BFEO's
claim is limited to ensemble marginalization and protected state information flow.

**BERP is an optional performance module.** It replaces B1's projected ligand latent with `[1,z_phys]` and performs one exact target-level
Bayesian update. Its implementation exactly matches B1 when both see the same design matrix. Its novelty
is provisional: only a gain specific to interpretable physical channels over ligand-latent B1 can make
it load-bearing.

## 5. Data and label audit

Local data contain 401 AlphaFold PDB files but no paired experimental pose or apo/holo registry. Public
PLINDER, MISATO and BigBind are reachable, but only a PLINDER subset fits a responsible local workflow.
MISATO trajectories describe bound-complex dynamics and are not binding free-energy labels. BigBind
uses ChEMBL33, so its affinities are ancestor data rather than independent evidence.

## 6. Leakage threat model

The external registry must block target identities and homologs, pocket/template neighbors, exact and
parent compounds, scaffold and ligand-similarity neighborhoods, pose neighbors, current affinity edges,
assay/document replication, report roles and sealed roles before any protected label is deserialized.
Generated poses must be selected without affinity or query-label feedback.

## 7. Strongest baselines

The required controls are current B1, pLM+ligand-GNN concatenation, ordinary atom-residue
cross-attention, complex-only invariant 3D GNN, IPBind-style single-state bound-minus-free, multi-state
mean pooling, and identical BFEO with ligand-latent B1. Parameter and compute matching is mandatory for
scored neural comparisons.

## 8. Complete information flow

```text
sequence / apo coordinates + ligand graph
-> label-blind pocket and candidate-state generation
-> shared standard invariant encoder
-> complex - free protein - free ligand channels                 [BFEO]
-> confidence-aware state marginalization
-> zero-shot prediction + physical coordinate + state variance
-> one target-level Cholesky posterior over support residuals    [BERP]
-> prediction + aleatoric/state + epistemic uncertainty
```

## 9. Mathematical and engineering checks

Eleven focused tests pass: rigid-transform invariance, state-order invariance, exact cancellation without
protein-ligand cross edges, finite gradients, CUDA BF16, B1/BERP Cholesky equivalence, exact k=0 fallback,
target-ID isolation, grouped BERP shrinkage, nested-support contraction and support-dependent posterior
change. The configured Gate-0 model has 146,606
trainable parameters and an 18-dimensional physical coordinate.

A 50-step deterministic synthetic CUDA BF16 smoke completed in 11.10 seconds on the RTX 4060 Laptop GPU,
with finite loss/gradients and 84.84 MiB peak allocated memory. This is an engineering test, not an
affinity experiment, checkpoint or scientific effect.

## 10. Training objectives and generalization

Stage A is limited to native/decoy pose, bound-free contrast, state consistency and channel supervision
on target/pocket/ligand-disjoint external structures. Only a Stage-A mechanism pass permits one-seed
train-target affinity training with censored NLL and within-target ranking. BERP episodes reuse the
current k=0/4/8/16 protocol. Target, pocket, scaffold, assay/document, template and temporal sensitivities
are required before multi-seed work.

## 11. Ablation and statistics

The frozen arm matrix and thresholds are in `fire_dta_preregistration.md`. Primary effects use paired
target-cluster bootstrap confidence intervals, empirical MDE, all-fold guards, harmful-adaptation rate,
calibration and destruction-based gain removal. Auxiliary pose success cannot substitute for target
macro RMSE and within-target Spearman.

## 12. Practical validity boundary

BFEO channels are currently representation contrasts, not thermodynamic energy components. Missing or
low-confidence states must increase applicability warnings. Boltz-2/UMA outputs may be state/feature
teachers only, never affinity truth. The full PLINDER/MISATO plan exceeds local storage and compute;
subset selection is part of the scientific design, not a convenience sample.

## 13. Go/no-go decisions

- **GO:** retain B1; keep FIRE-DTA isolated; use the new core for deterministic tests and a registered
  external subset audit.
- **NO-GO:** no 15-25M build, no BigBind affinity warm-up, no full MISATO download, no one-seed ChEMBL
  training, no multi-seed run and no sealed evaluation yet.
- **DATA READY:** D0 passed with 9,264 train systems, 1,828 exact experimental-apo-linked train systems
  and fully audited model-independent atom caches. Training and model-input execution remain stopped so
  BFEO/BERP can be adjusted without changing the frozen data.
- **Next decision:** review the module adjustment, then construct pocket/edge/model features from the
  frozen atom caches before any structural mechanism probe.

No checkpoint, outer-report or sealed label was read. `sealed_test_consumed=false`.
