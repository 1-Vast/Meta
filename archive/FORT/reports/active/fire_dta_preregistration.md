# FIRE-DTA Conditional Preregistration

## Scope

Two concentrated innovation modules are registered: BFEO is the required core mechanism and BERP is an
optional performance enhancement. Multiple contribution claims are allowed only when they are
inseparable parts of these modules. The residue/atom encoder, graph construction, optimizer, censored
likelihood, Cholesky solve, conformal calibration and state generator are standard infrastructure. No
third predictive module, learned query gate, affinity teacher, retrieval system or assay de-noising path
is authorized. BFEO must pass on its own; BERP cannot rescue a failed BFEO mechanism.

The existing `model/pipeline.py` and B1 implementation remain unchanged. FIRE-DTA is isolated until all
preceding gates pass. `sealed_test_consumed=false`.

## Gate D0: external structural registry

Acquire only the PLINDER index, frozen split metadata and a 5k-20k-system high-quality subset. Select
systems using structure-quality and uniqueness fields without reading protected affinity labels.
Quarantine before label parsing:

- checkpoint, outer-report and sealed target accessions or homologous sequence clusters;
- exact/current parent compounds and near-identical ligand neighborhoods reserved for evaluation;
- report-set pocket/template neighbors;
- ChEMBL/BindingDB affinity edges duplicated by current or protected observations;
- post-cutoff structures or labels when a temporal protocol is used.

Require at least 5,000 high-quality training systems, at least 500 independent protein/pocket clusters,
at least 1,000 linked apo/holo pairs for the bound-free mechanism check, and a protein/pocket/ligand
similarity-disjoint validation set. Otherwise return `FIRE_DTA_STRUCTURAL_REGISTRY_INSUFFICIENT`.

Do not download full MISATO MD. The 0.32 GiB QM asset or a registered trajectory subset may be acquired
only after D0 passes. BigBind labels remain ancestor-data evidence and cannot serve as independent
validation.

## Gate A: cheapest BFEO mechanism probe

Use at most 20k structures and a Gate-0 encoder below 1M trainable parameters. Compare identical data,
initialization budget and optimizer steps:

| Arm | Model |
| --- | --- |
| A0 | complex-only invariant 3D GNN |
| A1 | IPBind-style single-state bound-minus-free model |
| A2 | multi-state model with ordinary mean pooling |
| A3 | BFEO without bound-free subtraction |
| A4 | BFEO without learned/confidence-weighted marginalization |
| A5 | complete BFEO |

Destructive controls: pocket shuffle, ligand shuffle, coordinate randomization, complex/free branch
swap, state-confidence shuffle, distal-residue control and template-neighbor deletion.

Primary Gate-A quantities are target/pocket-cluster macro native-pose ranking, decoy discrimination,
apo/holo consistency and, only on safe independent labels, target-disjoint affinity. BFEO must beat A1
and A2 by at least the empirical target-cluster MDE, improve every held-out fold, and lose at least 70%
of the incremental effect under the corresponding destruction. A complex-only or single-state control
matching A5 closes BFEO as `BFEO_NOT_LOAD_BEARING`.

No free-energy decomposition claim is allowed at Gate A. Passing pose or state tasks alone authorizes
only a one-seed train-target affinity probe; it does not establish the deployment estimand.

## Gate B: one-seed BFEO cold-target probe

Seed `1729`; train targets only; current target episode and leakage firewall. Required before any
15-25M scale-up:

- k=0 target-macro RMSE gain over current B1 >= 0.03;
- within-target Spearman gain >= 0.02;
- all target folds positive, with no fold worse than 0.015 RMSE;
- pocket and ligand destruction remove >=70% of the incremental gain;
- parameter/compute-matched cross-attention, single-pose 3D GNN and IPBind controls do not match it.

Failure returns `BFEO_COLD_TARGET_FAIL_STOP`. BERP is not evaluated after BFEO failure.

## Gate C: BERP adaptation probe

BERP retains B1's single exact target-level posterior over a low-dimensional `[1, z]` design; B1 is
already a Bayesian last-layer residual model, not merely a scalar random intercept. The BERP prior has
four globally learned precision groups: intercept, BFEO channel means, BFEO state variances and the two
ensemble summaries. Reliability and structural state variance enter only as support likelihood
precision (a conservative tempered-evidence interpretation). No target-conditioned hyper-prior is
allowed because constant/random-prior controls already explained the earlier AHP/PME gain.

Freeze the passed BFEO architecture and compare:

| Arm | Adaptation |
| --- | --- |
| C0 | no support adaptation |
| C1 | current ligand-latent B1 |
| C2 | physical coordinates with ordinary ridge/MAP |
| C3 | BERP exact posterior |
| C4 | BERP with physical channels permuted |
| C5 | BERP with cross-target or label-permuted support |

BERP must add at least 0.02 RMSE and 0.02 Spearman beyond the strongest identical-BFEO control, reduce
harmful adaptation by at least 20%, have no negative fold, preserve k=0 exactly, and lose its gain under
support and physical-channel destruction. Multi-anchor or query gates are not allowed rescue paths.

## Scale and stop rule

The implemented 146,606-parameter core is an engineering object only. A 15-25M model, BigBind affinity
warm-up, full MISATO trajectories, multi-seed training and sealed evaluation are prohibited until their
preceding gates pass. A successful one-seed Gate B/C stops for review before long training, as requested.
