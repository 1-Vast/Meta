# HTL-DTA Baseline Preregistration

Date: 2026-07-31  
Status: protocol registration; pKi A2S baseline execution may proceed after
recipient roster validation

## Objective

Estimate the value of **abundant-to-scarce cross-target transfer** before
implementing a new expert-transfer architecture. For every support budget
`k in {1,3,5}`, the primary estimand is the paired target-macro change
from a transfer arm to a no-transfer arm on the same scarce recipient roster.
The resulting transfer-gain curve and its area under the curve are reported;
absolute performance alone is insufficient. Pseudo-tail and natural-tail are
separate estimands. `k=10` is conditional on query depth, `k=20` is restricted
to pseudo-scarce/medium-resource upper bounds, and `k=0` is a secondary
zero-shot diagnostic. pKi is primary; pKd is secondary replication.

## Data and firewall

The only active source is ChEMBL-37 dual-cold TRAIN/development metadata and
affinity rows after the HTL-1 admission gates. pKi and pKd are separate
endpoint strata. Confirmation, Davis, sealed, pKd_app, pIC50, and censored
inactive rows are excluded. The historical confirmation-consumption caveat in
the migration report remains binding; this preregistration does not reopen any
protected partition.

Support selection is covariate-only. Every declared cold axis closes target
homology component, scaffold, parent connectivity, document, assay, and source
lineage. A target or component, not a row, pair, query, seed, or fold, is the
primary statistical unit.

Source targets and recipient targets are target-disjoint. Only head-target
labels enter source expert fitting; recipient labels enter through the declared
support set and are never reused as source data for that recipient's query.
Homology-compatible and homology-cold source/recipient strata are scored
separately when each has sufficient power.

The preliminary topology audit uses the candidate source rule `n_eff >= 100`
and recipient rule `n_eff < 30`. It finds 242/193 pKi source/recipient targets
and 41/256 pKd source/recipient targets. Because only 12 pKi and 5 pKd
recipients share a homology component with a source, similarity-only transfer
is a control rather than an assumed solution. Final thresholds and rosters are
still frozen only after closure and power checks.

## Evaluation tracks

### Pseudo-scarce source control

Use data-rich head targets as pseudo-scarce recipients, with support sizes
`k in {1,3,5}`. For each target/seed, freeze the target-side single-cold query
and the following strengthening strata before labels are scored:

1. random chemical;
2. scaffold-cold;
3. low-similarity chemical.

`k=10` and `k=20` are included only when the recipient has enough closed query
units; `k=20` is not used for the natural scarce-tail primary claim.

Pseudo-tail success is reported as sample efficiency and cannot be called
natural-tail success.

### Natural-tail

Use genuinely low-resource targets admitted by HTL-3. Queries must be held out
by document/time/source where the source supports that ordering, and must pass
scaffold and chemical-component closure. Random deletion of rows is not an
admission rule. Targets without adequate query depth or component-aware power
are excluded before model scores are inspected.

### Chemical and homology strengthening tracks

The primary track is target-side single-cold: a query compound may have
appeared on source targets. Report separate strengthening tracks for global
drug-cold, scaffold-cold, and homology-cold recipients. Do not require their
intersection for the primary claim.

### Secondary strict stress test

The legacy strict unseen-target `k=5` roster is reported separately as Track C
only. It cannot promote a protein-conditioned architecture or rescue a failed
natural-tail result.

## Required baselines

The following no-transfer and transfer controls are frozen before any
head-to-tail transfer model:

* B0 ligand-only shared ridge/calibration baseline;
* per-target ridge and random forest where support is sufficient;
* pooled random forest/ExtraTrees;
* PCM using frozen target and ligand descriptors;
* kNN-DTA with a pre-registered similarity metric;
* shared multi-task target-conditioned model;
* standard MAML-style adaptation;
* AdaMBind-style task scheduling with the same encoder budget;
* random expert mixture;
* protein-similarity-only transfer; and
* support-compatibility-only transfer.

For pseudo-tail, the recipient is sampled from a head target but its support is
restricted to `k`; the remaining labels are query-only. The source pool is
constructed from other data-rich head targets, so the experiment measures
abundant-to-scarce transfer rather than within-target label leakage.

The random expert and similarity-only arms are mandatory controls for source
selection. A parameter-padded model is not a compute-matched control.

## Metrics and uncertainty

Report, per endpoint and track:

* target-macro RMSE and MAE;
* within-target Spearman and pairwise accuracy;
* NDCG@10 for ranking;
* sample-efficiency curve and area under that curve;
* paired transfer gain `Delta(k)` against the no-transfer arm and transfer-gain
  AULC;
* negative-transfer rate relative to B0 and calibration-only transfer;
* head-target retention;
* uncertainty coverage and risk-coverage curves.

Aggregate confidence intervals with a target-level or independent homology-
component bootstrap. If targets share a component, the component is the
resampling unit. Report the number of analyzable targets, components, and
source/document families in every table.

## Seeds and compute

Run one fixed seed for the admission smoke and baseline comparison. Run five
seeds only after the one-seed comparison passes the pre-registered effect and
diagnostic gates. All numerical paths use `D:\anaconda\envs\drug` and CUDA;
record GPU name, utilization, power, peak allocated memory, wall time, model
parameters, gradient evaluations, and episode count.

## Admission and stop rules

Before scoring, freeze the HTL-1 threshold and the HTL-3 roster only after
scaffold/document/assay/source closure and component-aware MDE calculation.
Baseline execution is admitted only if at least 50 pKi natural-tail targets
survive the declared `k=5` query depth and the corresponding component power
floor. pKd has its own power decision.

Stop the route if any of the following holds:

* a proposed transfer arm does not beat B0 and calibration-only baselines on
  both an error and a ranking metric;
* gains exist only for random-chemical pseudo-tail queries;
* random experts match learned source selection;
* most tail targets are harmed or negative-transfer rate is not below the
  no-gate control;
* head-target performance materially regresses;
* uncertainty is uncalibrated or coverage is not reported; or
* a result depends on relaxed closure, a post-hoc threshold, confirmation,
  structural input, or an unregistered loss.

No architecture claim is made from a baseline comparison alone. Only after the
baseline kill test passes may the two registered hypotheses in `task.md` be
implemented: resource-aware head-to-tail expert transfer and
negative-transfer abstention.
