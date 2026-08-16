# HTL-DTA Task Migration

Date: 2026-07-31  
Program: `HTL_DTA`  
Status: topology-first; training not yet authorized

## Decision

FORT is migrating its primary benchmark from strict unseen-target `k=5`
adaptation to **A2S-DTA: Abundant-to-Scarce Few-Shot Target Transfer**. The
strict task remains a secondary stress test and is not a capacity-rescue route.

The migrated question is whether **abundant head-target data can transfer to a
scarce recipient target** and improve prediction, ranking, uncertainty, and
experiment selection at the same support budget as no-transfer baselines, while
identifying and suppressing negative transfer. The primary output is a paired
transfer-gain curve over `k in {1,3,5}`; `k=10` is conditional on closed query
depth, `k=20` is restricted to pseudo-scarce/medium-resource upper bounds, and
`k=0` is a secondary zero-shot diagnostic. The primary track is target-side
single-cold; global drug-cold, scaffold-cold, and homology-cold are separate
strengthening tracks. This is a change in the estimand and evaluation topology, not a
relaxation of endpoint, provenance, chemical-closure, public-data, or
confirmation-firewall rules. Source targets and recipient targets are
target-disjoint: recipient labels enter only through the declared support set,
never through a source expert trained on that same target. Homology-compatible
and homology-cold transfer strata are reported separately where power permits.

## Why migrate

The existing strict route has not established a stable protein-conditioned
signal. In the corrected low-capacity probe, the protein-conditioned result was
worse than the protein-free control (`RMSE gain -0.0709`, 95% CI
`[-0.1221,-0.0200]`). The public crossed measurements also do not provide
verified protocol-comparable independent units for a strong interaction claim.
Adding another Transformer/Mamba interaction block would therefore confound
capacity with task identifiability.

The TRAIN-only ChEMBL-37 dual-cold registry does contain a real target-resource
imbalance:

| Endpoint | Targets | Median rows/target | Targets `<5` | Targets `<10` | Targets `>=100` |
| --- | ---: | ---: | ---: | ---: | ---: |
| pKi | 559 | 68 | 40 | 61 | 242 |
| pKd | 407 | 17 | 124 | 168 | 41 |

Median unique-scaffold depth is 34 for pKi and 12 for pKd. These are topology
signals, not independent experimental units; all admission decisions must use
provenance- and component-aware counts.

Under a provisional transfer topology (`source n_eff >= 100`, `recipient
n_eff < 30`), the TRAIN metadata contains 242 pKi source targets and 193 pKi
recipient targets. Only 12 pKi recipients share a homology component with a
source; 181 are homology-cold to the source pool. The corresponding pKd counts
are 41/256, with 5 same-component and 251 homology-cold recipients. The task is
therefore large enough to test abundant-to-scarce transfer, but a
protein-similarity-only route is not a sufficient prior and must be treated as
a control.

## Estimand and roles

For endpoint `e` in `{pKi, pKd}`, define `n_eff(t,e)` as the number of unique
endpoint-consistent parent-target-provenance units after document/assay/source
deduplication and chemical closure. Head, medium, tail, and extreme-tail
thresholds are frozen from the empirical `log(n_eff)` distribution before any
model score is inspected.

Two roles are kept separate:

* `pseudo_tail`: head targets whose support is covariate-only restricted to
  `k in {1,3,5,10,20}`. Query tracks are random-chemical, scaffold-cold, and
  low-similarity. This measures sample efficiency only.
* `natural_tail`: genuinely low-resource targets with document/time/source-held-
  out queries, scaffold and chemical-component closure, and enough effective
  target-level units for power. Random row deletion cannot create this role.

No pKi/pKd pooling is allowed. pKd, pKd_app, and pIC50/IC50 have separate
semantics; censored IC50 remains a separate ordinal/safety task or audit-only
source.

## Required evaluation

The primary unit is a target or independent target component, never a raw row,
query, ligand pair, or random seed. Required outputs are target-macro RMSE/MAE,
within-target Spearman, pairwise accuracy, NDCG@10, sample-efficiency AULC,
negative-transfer rate, head-target retention, and uncertainty
coverage/risk-coverage. Pseudo-tail and natural-tail results must be reported
separately.

Baseline kill tests must include ligand-only B0, per-target RF/ridge, pooled
RF/ExtraTrees, PCM, kNN-DTA, a shared multi-task model, standard MAML,
AdaMBind-style scheduling, random expert mixtures, protein-similarity-only
transfer, and support-compatibility-only transfer.

## Proposed innovations after the gates

At most two innovations will be admitted after the benchmark is established:

1. **Resource-aware head-to-tail expert transfer**: combine head-target
   adapters using protein relation, support compatibility, chemical-space
   compatibility, and source reliability; apply only a rank-1/2 residual or
   equivalent small update.
2. **Negative-transfer abstention**: a support-only leave-one-out/evidence gate
   may set transfer weight to zero when the expert mixture does not beat the
   global or calibration-only baseline.

These are hypotheses, not current claims. Scheduler details are implementation
controls, not a third innovation.

## Gates and current blocker

The next authorized computation is a read-only HTL-1 topology audit producing
endpoint-separated `n_eff`, scaffold/component depth, document/assay/source
counts, homology components, candidate strata, query depth, and effective
component power summaries. Training is blocked until HTL-1 and the source-role
audit pass. Natural-tail admission additionally requires provenance-separated
held-out queries and positive power after closure.

The audit must also preserve the governance caveat that a prior cross-document
artifact recorded `confirmation_labels_read=true`; therefore the program does
not claim that confirmation labels have never been read. Any contaminated
confirmation partition remains quarantined, and the current migration reads no
new sealed or confirmation target-conditioned labels.

## Source boundary

Only public sources with version, checksum, license/usage terms, schema, role,
endpoint semantics, lineage, and overlap records may enter active roles.
ChEMBL-37 dual-cold TRAIN/development is the current source. Papyrus is
aggregated and not an independent replicate source; BindingDB is target-shallow;
Reinecke is one development kinase campaign; and SPD is predominantly censored
inactive. Public availability alone does not make these sources admissible for
pooled inference.

## Reproducibility artifacts

The migration contract is in `task.md`, the source summary is in
`DATASET_RECORD_SUMMARY.md`, and the control-plane state is in
`manifests/state.v1.json`. The planned topology output is
`dataset/processed/htl_target_topology.v1.json`; no model checkpoint or
training result is implied by its generation.
