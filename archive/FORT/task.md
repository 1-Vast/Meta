# FORT Active Task

## Current Verdict

`DATA_NOT_READY` (data gate) / `CORE_META_MECHANISM_IDENTIFIED` (design gate)

Gate D0 is complete and returned a DATA STOP. The executed A2S-CFRA linear
router remains stopped for negative transfer, the natural-tail roster has zero
strictly admitted recipients, and no affinity-model training is authorized.

**2026-08-01 amendment.** The D0 stop is now attributed to three protocol
defects rather than to a structural absence of natural-tail data, and the
primary model is re-selected. See
`2026-08-01 A2S-SDO Preregistration And D0 Protocol Revision`. The binding next
action is a `dataset-run` rebuild of the natural-tail roster under the
within-recipient document-ordered estimand, not further historical-snapshot
acquisition. Model training remains unauthorized until that rebuild passes.

The former primary task, strict unseen-target k=5 adaptation, is retained as a
secondary stress test but is no longer the active development objective. The
current data do not identify a stable protein-conditioned signal for that task:
the corrected low-capacity probe is worse than protein-free, and the public
crossed measurements lack verified protocol-comparable independent units.

The active program is now **A2S-DTA: Abundant-to-Scarce Few-Shot Target
Transfer** (implemented under the `HTL_DTA` control plane). Its core estimand
is abundant-to-scarce target transfer: only data-rich source targets are used
for large-scale fitting; a target-disjoint recipient target contributes `k`
support labels at test time; the model predicts that recipient's query
compounds. Long-tail topology defines source and recipient roles; it is not the
claim by itself. No new architecture training is authorized until the
target-frequency topology and recipient roster pass the stages below.

No structure teacher, pose, pocket, private label, future measurement, or sealed
outcome enters the migrated program. This is a task-definition migration, not a
relaxation of public-data, endpoint, provenance, or falsification requirements.

## Scientific Question

For endpoint `e` in `{pKi, pKd}`, let `n_eff(t,e)` be the number of unique,
endpoint-consistent parent-target-provenance units for target `t`, after
document/assay/source deduplication and chemical closure. Partition targets
into head, medium, tail, and extreme-tail strata using a frozen empirical
distribution of `log(n_eff)`.

Given a target-disjoint scarce recipient with a small, legally available
support set `k in {1,3,5}`, does transfer from abundant source targets improve
over no-transfer baselines:
no-transfer baselines:

1. scaffold-cold query affinity and within-target ranking;
2. uncertainty-aware experiment prioritization; and
3. tail-target performance without increasing negative transfer?

The primary claim is a positive, target-macro **transfer gain curve** over
`k={1,3,5}`, not universal protein understanding. `k=10` is conditional on
enough closed query units; `k=20` is limited to pseudo-scarce or
medium-resource learning-curve upper bounds. `k=0` is a secondary zero-shot
diagnostic. The main task is target-side single-cold: query compounds may have
appeared on source targets. Global drug-cold, scaffold-cold, and
homology-cold variants are separate strengthening tracks, not a joint primary
requirement.

AdaMBind remains a literature and equal-budget baseline. Its random/CD-HIT
protocol does not define this task because it does not close support/query
scaffolds, chemical neighbours, documents, assays, and source provenance.

## Permissions And Data

`affinity_training_authorized = false`, `active_training_authorized = false`,
and `new_numeric_affinity_read = false`. Confirmation, Davis, and sealed rows
remain excluded. Until a new source passes D0, only outcome-blind metadata
needed for source admission and roster auditing may be read. pKi is the primary
endpoint; pKd is a secondary replication endpoint and is never pooled with
pKi. Any later numerical path must use CUDA in the `drug` environment and
record compute utilization, power, memory, wall time, and peak memory.

The migrated task has two data roles:

- `pseudo_tail`: head targets with a covariate-only support/query simulation;
- `natural_tail`: genuinely low-resource targets with a pre-registered
  document/time/source-held-out query.

For both roles, source targets and recipient targets are target-disjoint. The
source pool contains only data-rich head targets; recipient labels may enter
only through the declared support set. Transfer is reported in separate
homology-compatible and homology-cold strata where the roster has power. A
recipient target must not become a source expert for its own query, and a source
target must not be counted as an independent recipient episode.

Pseudo-tail results cannot be described as natural-tail success. A route that
passes only pseudo-tail is a sample-efficiency result, not a solution to the
real low-resource target problem.

## Public Data Source Boundary

`public_datasets_only = true`. Every source used for training, validation,
feature construction, pretraining, teacher or privileged information,
calibration, support selection, and evaluation must be a publicly accessible
dataset with a verifiable version, checksum, license or usage terms, schema,
and provenance record. Private, proprietary, unpublished, newly measured, or
sealed data are outside this program.

Public availability is necessary but not sufficient. A source may enter an
active role only after an audit records its exact version, allowed fields,
endpoint semantics, assay/document/source lineage, licensing boundary, target
and homology overlap, ligand/scaffold/chemical-neighbour overlap, and role
assignment. Public datasets must not be pooled merely to inflate sample size;
duplicate records and shared provenance are collapsed before statistical
inference, and pKi/pKd remain separate.

External-paper datasets, including AdaMBind-related BindingDB, Davis, and KIBA
files, are baseline or audit material until they pass this same public-source
and strict-closure protocol. If no audited public source supplies the
target-specific interaction information required by a proposed mechanism, the
route is stopped; prospective measurement is not an implicit fallback.

The frozen preprocessing protocol for HTL-DTA is:

1. Separate pKi, pKd, pKd_app, and pIC50/IC50 before all construction. pIC50
   and censored IC50 are audit-only or a separate ordinal safety task; they are
   never converted into pKi/pKd or pooled with them.
2. Compute `n_eff` from unique parent-target-provenance units, not raw rows.
   Record target counts, unique scaffolds, chemical components, documents,
   assays, sources, and homology components.
3. Freeze head/medium/tail thresholds from the target-count distribution before
   model scores. A candidate threshold grid is `{5,10,30,100,300}`, but the
   final thresholds require target count and power checks.
4. Build target-disjoint A2S episodes with primary support sizes `k in {1,3,5}`.
   Use `k=10` only when the recipient retains enough closed query units and
   `k=20` only for pseudo-scarce or medium-resource upper-bound curves. The
   primary query track is target-side single-cold; chemical similarity,
   scaffold-cold, and global drug-cold are reported as separate strata.
5. Build natural-tail episodes from real low-resource targets. Hold out later
   document/time/source units where available; never manufacture a natural tail
   by randomly deleting rows and then call it prospective evidence.
6. Close target homology, support/query scaffold, chemical component, document,
   assay/protocol, and source lineage wherever the track declares the axis
   cold. Retain source-row IDs and deterministic episode hashes.
7. Fit B0 and all nuisance statistics out-of-fold within the same chemical and
   provenance closure. Target-macro weighting is mandatory; row-frequency
   weighting is not an acceptable primary objective.

The existing `dataset/processed/episodes.v1.json` is a legacy strict k=5
roster. It remains available for Track C stress testing but is not the HTL-DTA
natural-tail registry.

## Active Research Modules

The task benchmark is the first contribution. No new architecture is active
until the benchmark passes topology and power gates. If a model is admitted,
it may claim at most two core innovations:

1. **Resource-aware head-to-tail expert transfer**: head-target adapters are
   combined using protein relation, support compatibility, chemical-space
   compatibility, and source reliability. The tail update is rank 1-2 or a
   comparably small residual; it is not a collection of full target networks.
2. **Negative-transfer abstention**: a support-only leave-one-out/evidence gate
   may set transfer weight to zero when the expert mixture does not beat the
   global or calibration-only baseline. A scheduler is an implementation
   control, not a third innovation.

Frozen ESM/Mamba/Transformer features, Morgan/physicochemical features, B0,
per-target calibration, kNN-DTA, RF/ridge, MAML, and AdaMBind-style scheduling
are infrastructure or controls. The former Bayesian reordering adapter remains
an archived baseline; it is not the new primary path. Parameter padding is not
a compute-matched control. Report parameters, gradient evaluations, peak
memory, wall time, and episode count for every ablation.

## Current Evidence

The local ChEMBL-37 dual-cold TRAIN registry already shows a target-frequency
tail, but the counts are not independent experiments:

| Endpoint | Targets | Median rows/target | Targets <5 | Targets <10 | Targets >=100 |
| --- | ---: | ---: | ---: | ---: | ---: |
| pKi | 559 | 68 | 40 | 61 | 242 |
| pKd | 407 | 17 | 124 | 168 | 41 |

Unique scaffold medians are 34 for pKi and 12 for pKd. These numbers justify a
topology audit and pseudo-tail benchmark, not an immediate transfer claim.
After provenance closure, the current audit still finds zero verified same-assay
crossed rectangles, a largest document-source family of about 69.2%, and
DD/noise(q90) of about 0.802 for pKi and 0.296 for pKd. Papyrus is aggregated,
BindingDB is target-shallow, and Reinecke is one development kinase campaign.

The former strict protein-conditioned route remains stopped: the IDG-RBP
correct-protein RMSE gain over protein-free is `-0.0709 [-0.1221,-0.0200]`.
This is why HTL-DTA starts with task topology and target-level baselines rather
than another interaction-capacity upgrade. Details are in
`reports/active/pcic_rr_dta_research_design_2026-07-31.md`.

## Ordered Stages

| Stage | Status | Admission |
| --- | --- | --- |
| HTL-0 source/role audit | PARTIAL / TRUE TIME MISSING | Current source/version/lineage are hashed, but publication or measurement dates are unavailable. |
| HTL-1 target-frequency topology | STOP / DATA_NOT_READY | D0's optimistic closed upper bound is 40 pKi recipients, below the frozen minimum of 50. |
| HTL-2/A2S recipient roster | PASS (single-cold) | Target-disjoint pKi roster and k={1,3,5} episodes are constructed; strengthening closure tracks remain separate. |
| HTL-3 natural-tail roster | STOP / DATA_NOT_READY | Strict roster has 0 recipients and 0 independent components; MDE80 is undefined. |
| HTL-4/A2S baseline diagnostic | CORRECTED / CROSS-FITTED ROUTER NO-GO | Global feature-row alignment and target-macro source weighting are repaired; pKi source-support control is positive, but the held-out router is negative at k=3/5. |
| HTL-5 head-to-tail transfer | BLOCKED | Requires HTL-1 through HTL-4 PASS; one seed first, then five seeds. |
| HTL-6 budgeted screening | BLOCKED | Requires predictive and uncertainty gains before acquisition-function testing. |
| Track C strict unseen-target stress test | SECONDARY/STOPPED | Reopen only as a stress test after a separately powered source passes; never use it as a capacity rescue. |
| C0 confirmation/sealed | FORBIDDEN | Separate authorization only. |

## PASS And STOP Rules

HTL-1 PASS requires endpoint-separated counts, at least 50 analyzable tail
targets for the primary pKi track where the k=5 query depth survives closure,
and a documented effective-component MDE. pKd is secondary and must pass its
own power audit. Thresholds are frozen before model scores.

HTL-3 PASS requires genuine natural-tail provenance separation, no source or
document family dominating the target-level bootstrap, adequate query depth,
and a positive power calculation after scaffold/chemical closure. A random row
deletion is not natural-tail evidence.

HTL-4 PASS requires target-macro gains over global B0, per-target RF/ridge,
kNN-DTA, PCM, standard MAML, and AdaMBind-style scheduling on both pseudo-tail
and natural-tail where the latter is admissible. Report RMSE/MAE, Spearman,
pairwise accuracy, NDCG@10, uncertainty coverage, AULC, and negative-transfer
rate. The primary statistical unit is target or independent target component,
not row, pair, query, or seed.

HTL-5 PASS additionally requires: most tail targets benefit; gains survive
scaffold-similarity strata; learned source selection beats random and
protein-similarity-only transfer; negative-transfer abstention reduces harmed
targets; head performance does not materially regress; and 95% target-bootstrap
lower bounds exceed `max(MDE80, material floor)`.

STOP for pseudo-tail-only gains, nearest-neighbour-only gains, head-target
dominance, random-source equivalence, no natural-tail power, source/document
leakage, uncalibrated uncertainty, or negative-transfer rates no better than
the no-gate baseline. Do not rescue a stop with capacity, epochs, seeds,
unregistered losses, relaxed chemical closure, confirmation, structural inputs,
or post-hoc thresholds.

## Active Artifacts

- `main.py`
- `model/protein.py`, `model/ligand.py`, and `model/interaction.py` (encoding)
- `model/posterior.py`, `model/reorder.py`, and `model/likelihood.py` (inference)
- `model/ligandbase.py` and `model/gradadapt.py` (comparison baselines)
- `scripts/preprocess.py`, `scripts/audit.py`, `scripts/episode.py`, and
  `scripts/htl_topology.py`, and `scripts/train.py`
- `research/a2s_baseline.py`
- `dataset/processed/episodes.pKi.v1.parquet`
- `dataset/processed/episodes.pKd.v1.parquet`
- `reports/active/fewshot.v1.json`
- `reports/active/stability.v1.md`
- `reports/active/bayeskill.v1.json`
- `reports/active/bayeskill.v1.md`
- `reports/active/pcic_rr_dta_research_design_2026-07-31.md`
- `reports/active/htl_dta_task_migration_2026-07-31.md`
- `reports/active/htl_baseline_preregistration_2026-07-31.md`
- `reports/active/a2s_alignment_correction_decision_2026-07-31.md`
- `reports/active/a2s_dta_master_design_2026-07-31.md`
- `reports/active/a2s_invalid_artifacts_2026-07-31.md`
- `reports/active/a2s_baseline_decision_2026-07-31.md` (superseded)
- `reports/active/a2s_crossfitted_router_preregistration_2026-07-31.md`
- `manifests/state.v1.json`
- `configs/a2s_dta_minimal.v1.json`

Planned HTL-DTA artifacts (not yet generated and not an authorization to
train):

- `dataset/processed/htl_target_topology.v1.json` (preliminary metadata-only audit; thresholds not frozen)
- `dataset/processed/htl_episodes.pKi.v1.parquet`
- `dataset/processed/htl_episodes.pKd.v1.parquet`
- `reports/active/a2s_pki_targetbalanced_seed1729.json` (pKi primary baseline)
- `reports/active/a2s_pkd_targetbalanced_seed1729.json` (pKd secondary replication)
- `reports/active/a2s_router_pki_targetbalanced_seed1729.json` (cross-fitted router NO-GO)

## Closed Routes

| Family | Established fact | Reopening condition |
| --- | --- | --- |
| A1 structure | 32 complete units vs frozen 128 minimum | A separate compliant source route |
| DCST and source transfer | Strict interaction transport was not identified | New deployable input and preregistered evidence |
| Gradient task-code primary | E0 did not exceed calibration and lacks posterior evidence | Baseline comparison only |
| Flexible-kernel posterior | Finite-rank Wave 1 failed ligand-only and gradient controls | New preregistered mechanism after failure diagnosis; no hyperparameter rescue |
| pKd primary arm | 51 strict k=5 episodes, median query depth 7 | Powered endpoint-specific audit |
| Strict unseen-target k=5 | Correct protein failed protein-free in the low-capacity probe | A new powered, provenance-comparable source; retained only as Track C |
| HTL-DTA transfer | New active task | HTL-0 through HTL-4 must pass before head-expert implementation |

## HTL-1 Preliminary Audit Result

The read-only TRAIN metadata audit generated
`dataset/processed/htl_target_topology.v1.json`. It reads no affinity or
replicate labels and excludes confirmation/development rows by parquet split
filter. It reports 201,827 TRAIN rows, 559 pKi targets, and 407 pKd targets.

Using the candidate rule `n_eff < 10`, only 5 pKi targets have both at least
six unique scaffolds and six provenance units, so that threshold cannot satisfy
the primary tail-power requirement. The unselected candidate `n_eff < 30`
provides 115 pKi targets under this optimistic upper-bound screen, but this is
not yet an episode admission: document/assay closure, natural-tail provenance,
independent component count, and MDE still need to be checked. Thresholds
remain unfrozen and no model training is authorized.

For a separate candidate transfer topology (`source n_eff >= 100`, `recipient
n_eff < 30`), pKi has 242 source targets and 193 recipient targets; only 12
recipients share a homology component with a source and 181 are homology-cold
to the source pool. pKd has 41 source targets and 256 recipients, with 5
same-component and 251 homology-cold recipients. This supports a real
abundant-to-scarce benchmark, but it also means protein-similarity-only transfer
cannot be assumed to solve it; both compatible and homology-cold strata must be
scored separately.

The initial CUDA A2S artifacts were quarantined after a global registry-to-
feature row alignment audit found that filtered parquet indices had been used
as feature row ids. The corrected loader verifies the full 343,211-row
registry and `conn_sha` before filtering, and the primary ridge gives each
source target equal total weight. Target-balanced pKi support-compatible
source routing gains versus recipient calibration are `+0.137 [0.070,0.207]`,
`+0.112 [0.054,0.174]`, and `+0.087 [0.029,0.149]` at k=1,3,5. Corrected
secondary pKd gains are `+0.373`, `+0.220`, and `+0.242`, with positive
target-bootstrap intervals; endpoints remain separate. However, the
target-balanced cross-fitted pKi router is `-1.213` and `-1.315` at k=3/5,
so the registered router route is NO-GO and no new architecture training is
justified. The complete model, literature, metric, kill-test, and PASS/STOP
specification is frozen in `reports/active/a2s_dta_master_design_2026-07-31.md`.

The completed Gate D0 audit is recorded in
`dataset/processed/a2s_natural_tail_d0.v1.json` and
`reports/active/a2s_natural_tail_d0_decision_2026-07-31.md`. Among 193 candidate
pKi recipients, an intentionally optimistic five-draw parent/document/assay-
closed upper bound retains only 40 targets; adding distinct support scaffolds
retains 34. ChEMBL-release temporal closure retains at most 22, and strict
time/source/parent/scaffold/document/assay closure retains zero. Shared source
family or homology joins all 193 candidates into one dependency component, the
strict component count is zero, and MDE80 is undefined. This is a binding DATA
STOP before any A2S-MAP, MAML, AdaMBind, or encoder training.

## 2026-07-31 Historical Multi-Source Reopening Amendment

The final-snapshot D0 result above remains a fact and is not relaxed. The route
may reopen only under a different, prospectively frozen estimand: **historical
natural scarcity at an index date**, followed by genuinely later experiments.
This amendment is registered before reading any new affinity value or fitting
any model.

### Revised natural-tail estimand

Natural scarcity is defined at a frozen historical index date `tau`, not by the
recipient's final ChEMBL-37 count. Candidate index dates are fixed at
`2018-12-31`, `2020-12-31`, and `2022-12-31`. A recipient is eligible at its
earliest qualifying date only when all conditions hold:

1. it has 5-29 independent pre-`tau` pKi parent-provenance units;
2. it has at least 10 post-`tau` query units after parent, scaffold, document,
   assay, upstream-source, and time closure;
3. publication date is after `tau` and first-seen database release is after the
   frozen support snapshot for every query unit;
4. the target is removed from the abundant source pool for its own episode;
5. the source pool contains only targets with at least 100 closed pre-`tau`
   units; and
6. a target enters once, at its earliest eligible index date, so repeated time
   points are never counted as independent recipients.

The deployment question is therefore: given only the information publicly
available for target `t` at `tau`, can abundant targets improve prediction of
experiments that appeared later? This is not random deletion or pseudo-tail.

### Frozen source order and evidence roles

1. **HIST-S small development substrate:** versioned ChEMBL-37 records with
   document publication date and first-seen release. It is used first for
   metadata admission and, only after PASS, low-cost mechanism development.
2. **HIST-L large confirmation substrate:** BindingDB fixed-archive records
   with unique upstream origin plus GtoPdb 2026.2 unique references, combined
   with ChEMBL after cross-database provenance deduplication. It is not read for
   model selection after its components are assigned.
3. **Conditional supplement:** the traceable published PDSP Ki subset may be
   separately preregistered only if HIST-L topology is close to the frozen
   power floor.
4. **Mapping only:** Papyrus/Papyrus++ may discover or align records but never
   increases the independent source, document, target, or component count.
5. Davis, KIBA, AdaMBind CSV splits, converted IC50, censored measurements, and
   randomly thinned head targets are baselines or diagnostics, not D0 evidence.

Only exact `Ki` or native pKi with relation `=` and an unambiguous single human
protein, parent structure, assay, document, and upstream provenance are
admissible. pKd remains a later endpoint-separated replication and cannot
rescue pKi.

### Unified metadata and provenance contract

Before numerical affinity access, construct a versioned registry containing:

```text
database_source, upstream_source, source_record_id, target_uniprot,
target_sequence_hash, homology_component, endpoint_type, standard_relation,
assay_id, assay_description_hash, document_id, DOI, PMID, patent_id,
publication_date, curation_date, first_seen_release, parent_inchikey,
connectivity_key, scaffold_id, organism, target_type
```

The canonical document key priority is DOI, then PMID, patent family, then a
source-specific document ID. The dependency-level provenance family joins the
underlying document/patent family, institution or campaign, and assay family.
ChEMBL, BindingDB, GtoPdb, Papyrus, and PDSP are aggregators or databases, not
automatically independent upstream sources. The same DOI/PMID/patent/assay in
multiple databases is one unit.

### Small-to-large model validation

No model fit is authorized until HIST-S metadata D0 passes. If it passes, the
small phase uses frozen 1,034-D ligand features, frozen protein covariates, and
the least expensive A2S-MAP implementation. Recipient and provenance
components are assigned by immutable hash before labels are read:

- `DEV-S` is the small mechanism-development set. Nested
  `S_t^1 subset S_t^3 subset S_t^5` supports are label-blind and fixed across
  models. Query labels may tune only the registered mechanism within DEV-S.
- `LOCK-L` contains all held-out target/provenance components from the larger
  multi-source registry. Its outcomes are not used for architecture, rank,
  source-count, threshold, seed, or loss selection.

Small-data success is only an admission to large confirmation, never an
excellence claim. It requires at k=3 and k=5: positive target-macro RMSE gain
over the equal-capacity zero-prior fallback; at least one positive ranking
gain; majority-recipient benefit; at least a 10 percentage-point reduction in
negative transfer at coverage >=0.60; and loss of the effect under source-code
shuffle. Exact MAML, equal-budget AdaMBind, pooled fine-tuning, chemistry-only,
protein-only, random/all-source, and current affine support routing remain
mandatory controls.

Only a completely frozen small-phase model may enter HIST-L. Large-phase
success requires the original Excellence PASS rules, target/provenance-
component bootstrap, five registered seeds, and direction agreement across
ChEMBL-only and unique-origin external strata. Failure on LOCK-L returns to the
Failure Ledger; LOCK-L cannot be recycled into DEV-S.

### Compute placement, generalization, and reproducibility

This workstation is authorized only for metadata audits, tests, frozen-feature
probes, and the small one-seed mechanism kill. **HIST-L training must not run on
this device.** After small-phase PASS, the local deliverable is a portable,
hash-locked job bundle for an external GPU system; execution requires a
separate compute record and does not inherit local authorization implicitly.

An innovation is not accepted from an aggregate DEV-S improvement. Before
large execution, freeze all generalization analyses:

- leave-one-database-source-out and unique-origin external-source performance;
- leave-one-homology-component/family-out performance;
- earlier-to-later index-date transfer without backfilled documents;
- scaffold-cold and middle/low chemical-similarity strata;
- leave-one-provenance/campaign-component-out influence; and
- calibration, coverage, and negative-transfer rate under each shift.

Large PASS requires the gain direction to agree across ChEMBL-only,
BindingDB-unique, and GtoPdb-unique strata where powered; no single database,
protein family, document/patent family, assay campaign, chemical series, or
seed may carry the claim. A mechanism that wins only in-distribution or on one
source is a generalization STOP even if pooled RMSE improves.

Every small or large run must emit a reproducibility certificate containing:

```text
git_commit_or_diff_hash, code_bundle_sha256, environment_lock_sha256,
python_torch_cuda_versions, data_manifest_and_file_hashes, source_versions,
registry_and_roster_hashes, split_component_ids, episode_hashes, config_hash,
seed, command, hardware, wall_time, utilization, power, peak_memory,
parameter_and_gradient_counts, checkpoint_hashes, prediction_hashes,
metric_code_hash, and complete stdout/stderr log
```

The external runner must support a deterministic smoke replay and restart from
an immutable checkpoint. All target-level predictions and chosen transfer
actions are retained so statistics can be recomputed without retraining.
Missing hashes, environment lock, component assignments, logs, or predictions
invalidate the run; a reported metric alone is never reproducible evidence.

### Power and iteration discipline

Historical D0 reports MDE80 over independent dependency components for
`sigma_delta in {0.2,0.3,0.4,0.5}`. A conservative planning SD must be selected
from valid pre-existing, non-LOCK-L target-level contrasts before model scores;
the former `0.10` convenience proxy cannot authorize training. Fifty targets
is only a topology minimum and does not imply power for a 0.05 pKi gain.

Each optimization cycle must be preregistered before execution, contain at
most two core innovations, use seed 1729 first, and end in PASS or STOP. A
failed cycle cannot be rescued by extra epochs, more capacity, a wider atlas,
post-hoc recipients, relaxed closure, or immediate reuse of LOCK-L. This
permits iterative scientific improvement without converting the small dataset
into an architecture-search test set.

### Reopening gates

1. **H0-S metadata inventory:** verify versions, checksums, license, schema,
   true dates, upstream lineage, and cross-database overlap for local ChEMBL
   and BindingDB; locate or acquire the fixed GtoPdb release metadata.
2. **H1-S historical D0:** build the affinity-blind HIST-S registry and common
   roster, dependency graph, cutflow, and MDE envelope. Fewer than 50 closed
   recipients or inadequate power is a DATA STOP.
3. **S0-S static headroom:** test rank-2/4 residual-atlas headroom without
   encoder updates. Failure stops A2S-MAP.
4. **K1-S one-seed mechanism kill:** run the frozen A2S-MAP and destruction
   controls on DEV-S. Failure starts a new preregistered hypothesis cycle, not
   hyperparameter rescue.
5. **L0-L large metadata admission:** add only BindingDB/GtoPdb unique-origin
   records, freeze LOCK-L, and recompute dependency power.
6. **L1-L confirmation:** run the frozen winning mechanism, strong baselines,
   ablations, and five seeds on separately recorded external compute. Only this
   stage can support an excellence claim; it must not execute on this device.

### H0-S executed metadata inventory

**FACT:** `D:\anaconda\envs\drug\python.exe main.py h0-inventory
--acquire-remote-metadata` completed without reading or materializing an
affinity value and without invoking CUDA. The frozen artifact is
`dataset/processed/a2s_h0_metadata_inventory.v1.json`, SHA-256
`b8812489f4d8535f28fbd8c3bce5c3d6c39539da03031e82634e958594ad7942`.

**FACT:** ChEMBL document metadata now has publication year for 9,640/9,643
documents and a first-document-release date for 9,643/9,643. This does not
close historical time: the ChEMBL-37 registry has no activity-row first-seen
field, and the official ChEMBL 24.1, 27, and 31 database snapshots are not
local. A document release cannot prove that every current activity under that
document was present at that release.

**FACT:** BindingDB 202607 has the required document/date/source identifiers
in its 640-field header and version-matched assay companion archives are
available, but the main archive alone has no generic assay lineage. Official
GtoPdb 2026.2 metadata is now frozen locally with version, license, PubMed,
patent, assay-description, and InChIKey evidence; it has no per-interaction
first-seen history. Both remain HIST-L candidates and cannot increase HIST-S.

**INFERENCE:** H0-S returns `DATA_NOT_READY`. H1-S and all model fitting remain
blocked. The sole reopening route is historical-snapshot identity recovery;
publication year plus final-snapshot document release is an invalid backfill
shortcut.

### H0-S identity-projection protocol correction

This correction is registered before extracting a historical SQLite database,
executing an activity projection, reading a numeric affinity value, or fitting
a model. The ChEMBL-37 processed registry is a median aggregate without
`activity_id`; therefore 24.1/27/31 alone cannot bridge historical rows to the
current endpoint table. H0-S now requires the same identity-only projection
from **ChEMBL 24.1, 27, 31, and 37**.

Each archive must match the publisher byte count and SHA-256 and bind the
publisher checksum, license, and attribution files. The extracted database
must pass `PRAGMA integrity_check`, certify the exact release through its
`version` table, and freeze a normalized `sqlite_master` hash. File presence is
not admission evidence.

The projection connection is immutable and query-only. A SQLite authorizer
allows only `SELECT`, the fixed `upper/trim/coalesce` functions, and an exact
column allowlist on the `main` database; all writes, attachments, pragmas,
views that expand to forbidden fields, arbitrary functions, and undeclared
column reads are denied. Numeric or textual outcomes including `value`,
`published_value`, `standard_value`, `pchembl_value`, range values, text values,
and activity comments are forbidden in every SELECT, WHERE, JOIN, expression,
aggregate, ORDER, GROUP, HAVING, subquery, output, and log. Raw assay
descriptions and sequence bodies are also forbidden; target identity uses the
source-provided sequence MD5 plus accession. `standard_type`,
`standard_relation`, and `standard_units` remain explicitly declared
outcome-blind eligibility metadata because the frozen estimand requires exact
uncensored Ki and a known unit; they never authorize reading a measurement.

Release-local `activity_id` is only a candidate identity anchor. The audit also
freezes a SHA-256 semantic fingerprint over public source, document, assay,
compound, parent, target, endpoint, relation, and unit identities. Same native
ID with incompatible document/assay/source-compound/endpoint lineage is a hard
collision and H0-S STOP. Changed native IDs may link only when the semantic
fingerprint is one-to-one in every involved release. Many-to-many fingerprints,
join fan-out, document or molecule conflicts, target/parent drift, and any
presence pattern with deletion or reappearance are quarantined.

The cross-release output uses the ordered mask `(24.1,27,31,37)`. Only masks of
form `0*1+` are persistent. `first_seen_release` means earliest observed
database availability, never experiment or publication time. A later activity
under a document already present in an earlier full document registry is
classified as direct backfill. Publication year supplies a second, independent
end-of-year condition; missing or retrospective years are not future-query
evidence.

Document-presence matching uses the union of ChEMBL document ID, normalized
DOI, exact PMID, and normalized patent aliases. These aliases validate or link
document presence but do not form a DOI-else-PMID precedence identity whose
meaning could change when a later release enriches an alias.

Every activity and document projection must reproduce in an independent
process with identical row count, canonical row hash, Parquet hash, firewall
read log, and conflict/null counts. H0-S cannot pass until all four replay
certificates and the cross-release presence audit pass. H0-S PASS would only
authorize H1-S metadata roster construction; it would not authorize affinity
training. HIST-L training and confirmation remain prohibited on this
workstation regardless of H0-S outcome.

### Dataset-processing provenance gate

**FACT:** From this point onward, every dataset acquisition, resumed transfer,
checksum, extraction, projection, replay, presence audit, roster build, split,
and inventory refresh must execute through `main.py dataset-run`. The immutable
history root is `dataset/processing_history/v1`. A data artifact that is present
on disk but is not bound to a passing contemporaneous run remains uncertified.

Before a child process starts, the runner freezes an intent and prepared record
containing the exact argument vector, absolute working directory, controlled
child environment, executable hash, input and output pre-state, expected byte
counts and SHA-256 values, per-file code and configuration bundle hashes, Git
commit/diff/status fingerprints, platform, and runner version. The child uses
`shell=False`; raw stdout and stderr are written separately. After execution,
the runner records PID, UTC start/end, monotonic duration, exit code, output
paths/sizes/hashes, code/config post-state, input-mutation checks, and log
hashes. `process_status` and `validation_status` are distinct: exit code zero
cannot override a missing output, checksum mismatch, input mutation, or code
change, and neither status constitutes an H0/H1 scientific PASS.

Each run directory is create-once and ends with a content seal. The append-only
JSONL index is hash-chained; `dataset-run verify` must validate the entire chain
and every sealed file before a downstream artifact is admitted. Existing
outputs may be resumed only when explicitly declared as both input and output;
otherwise the runner refuses overwrite. Dataset Python commands must use
`D:\anaconda\envs\drug\python.exe` and are allowlisted to `historical-project`,
`historical-presence`, and `h0-inventory`. Curl is restricted to official
ChEMBL release downloads. Training, router, baseline, arbitrary Python, shell,
upload, and HIST-L commands are policy failures and are logged rather than run.

Primary projection and independent replay are separate logged child processes,
one release per process. The same rule applies to the cross-release presence
audit and its replay. Publisher byte count and SHA-256 are mandatory expected
outputs for archive transfers. Projection/presence manifests and replay
certificates remain the scientific commit markers; runner success alone cannot
certify their internal outcome-blind or identity contracts.

Operations completed before this gate are reconstructed only as explicitly
retrospective records. Such records must name their evidence sources and
limitations and mark unavailable command, PID, timing, exit code, stdout, or
stderr as unknown; an empty reconstructed log is never contemporaneous proof.
Retrospective records preserve the corrupt 24.1 transfer/quarantine, clean 24.1
and 27 transfers, partial 31 transfers, 24.1 extraction, earlier H0 inventories,
and both interrupted no-output projection attempts, but do not certify those
operations. The relevant final bytes must still be rehashed or replayed in a
new contemporaneous run before admission.

## 2026-07-31 Meta-Learning Committee Preregistration

This section supersedes the model-design and next-action portions of the
earlier A2S-CFRA record. It does not erase, reinterpret, or reopen any earlier
failure. Evidence statements use exactly these labels:

- **FACT**: established by a local artifact, source code, or executed result;
- **LITERATURE**: established in the cited primary paper or official code;
- **INFERENCE**: a conclusion drawn from facts or literature;
- **HYPOTHESIS**: an unverified claim that requires the registered experiment.

### Committee role conclusions

**FACT:** Agents 1-3 completed their task/history, meta-learning, and
cross-domain first-round reviews independently and did not read one another's
outputs. The Agent 4 design role and Agent 5 adversarial role were performed
only after those reviews completed; they are recorded separately below.

- **FACT - Agent 1, task/history audit:** the estimand is paired, target-macro
  pKi transfer gain for target-disjoint abundant sources and scarce recipients
  at `k={1,3,5}`. The current `n_eff>=100` and `<30` roles are provisional,
  thresholds are not frozen, and no document/time/source-closed natural-tail
  roster exists. The current same-roster MAML, AdaMBind, pooled fine-tuning,
  RF/ridge, PCM, kNN-DTA, and all-source comparisons have not been executed.
- **FACT - Agent 1, failure audit:** the corrected support-compatible source
  control is positive in RMSE, but every source prediction is only an affine
  transform of the same pooled `B0`. It therefore supplies calibration but no
  general target-specific reordering. The learned linear router is strongly
  negative and its uncalibrated gate accepts every evaluated recipient.
- **LITERATURE - Agent 2, meta-learning review:** MAML/FOMAML/Reptile learn a
  global initialization; ANIL restricts the inner loop to a head; BOIL changes
  the body; Meta-SGD learns update directions; LEO, CNAPs, VERSA, HyperMAML,
  and MetaDTA already establish support-conditioned or amortized task
  parameters; ATS and AdaMBind already schedule detrimental meta-training
  tasks. None of these mechanisms alone is new for A2S-DTA.
- **INFERENCE - Agent 2:** at `k<=5`, a frozen representation plus a
  low-dimensional, support-conditioned update is better identified than a
  full-model inner loop, but it remains source-blind unless source identity and
  source competence are retained explicitly.
- **LITERATURE - Agent 3, cross-domain review:** Taskonomy and Task2Vec motivate
  directed task affinity; LogME motivates support-label model evidence;
  DSelect-k supplies sparse expert selection; SelectiveNet supplies
  risk-coverage accounting; DML supplies target-level cross-fitting discipline;
  and conservative contextual bandits motivate a baseline-relative lower-bound
  action. Their guarantees do not transfer automatically to natural scarce DTA.
- **INFERENCE - Agent 3:** the reusable cross-domain quantity is paired,
  directed transfer gain relative to an identical no-routing fallback. Sparse
  MoE, similarity routing, predictive variance, or abstention alone is not a
  mechanism innovation.
- **HYPOTHESIS - Agent 4, model-design role:** source-specific residual
  functions can be factorized into a small task-code atlas. A recipient's
  support can identify a sparse source-code prior, after which only a
  low-dimensional code and intercept need adaptation.
- **INFERENCE - Agent 5, adversarial role:** this hypothesis is admissible only
  if a cross-fitted frozen-feature probe first shows ranking headroom, natural
  query closure and power, no source/recipient leakage, superiority to an
  equal-capacity zero-prior adapter, and useful abstention at nontrivial
  coverage. It is not presently authorized for training.

### Corrected current diagnosis

- **FACT:** globally aligned, target-balanced pKi support-compatible routing
  has RMSE gains `+0.1365 [0.0702,0.2067]`,
  `+0.1117 [0.0536,0.1743]`, and `+0.0873 [0.0287,0.1493]` at
  `k=1,3,5`. This updates the objective file's older statement that simple
  pKi routing had no positive k=3/5 result.
- **FACT:** the positive control is not a learned held-out router and does not
  improve ranking consistently. At k=3, Spearman changes from `0.0937` to
  `0.0897`; at k=5 it changes from `0.0829` to `0.0776`.
- **FACT:** the cross-fitted linear router has gains `-1.1441`, `-1.2134`, and
  `-1.3146` at k=1/3/5, with negative-transfer rates `0.801`, `0.774`, and
  `0.821`. Its score is positive for all 405 evaluated episodes, so its gate
  never falls back. Only 13-14 sources are selected at each k; `CHEMBL233`
  alone receives 60/151, 66/137, and 56/117 selections.
- **FACT:** the implemented meta-training pseudo-recipients have median query
  depths 320/318/316 at k=1/3/5, while the count-defined scarce recipients
  have medians 15/13/13. The implementation did not match the pseudo query
  distribution to natural scarcity.
- **FACT:** a token-level metadata reconstruction finds support/query document
  or assay overlap in 131/151, 134/137, and 114/117 episodes at k=1/3/5.
  Scaffold overlap occurs in 55/151, 89/137, and 81/117 episodes. Exact parent
  connectivity is disjoint in this reconstruction, but natural time/source
  separation is absent.
- **FACT:** the first A2S artifacts remain invalid because local filtered row
  indices were used against a global feature cache. They must never be cited.
- **INFERENCE:** the learned router failure is consistent with three jointly
  plausible causes: pseudo-to-scarce task shift, pointwise-utility argmax
  winner's curse over hundreds of sources, and collapse onto source-depth or
  provenance shortcuts. Existing artifacts do not identify which cause is
  dominant.
- **INFERENCE:** the current affine source library cannot satisfy the required
  causal chain for compound ranking. A new hypothesis must first add a
  source-specific residual direction that can reorder compounds, then test
  whether support identifies that direction.

### Meta-learning literature comparison

| Family | **LITERATURE:** meta-knowledge and inner loop | **INFERENCE:** A2S fit and boundary | Disposition |
| --- | --- | --- | --- |
| [MAML](https://proceedings.mlr.press/v70/finn17a.html) | One global initialization; canonical inner loop updates all learner parameters and differentiates through the update. | Original regression evidence does not cover k=1/3. It preserves no source identity and has no refusal rule. | Mandatory exact baseline. |
| [FOMAML/Reptile](https://arxiv.org/abs/1803.02999) | First-order initialization learning; Reptile moves parameters toward within-task solutions. | Cheaper, but still source-blind and vulnerable to full-model overfit. | Compute baseline. |
| [ANIL](https://arxiv.org/abs/1909.09157) | Shared body plus adaptable head initialization; only the head changes. | Restriction is appropriate for k<=5, but no source selection or negative-transfer decision is learned. | Strong baseline and design constraint. |
| [BOIL](https://arxiv.org/abs/2008.08882) | Frozen head and body adaptation to force representation change. | Large-body updates from <=5 affinities are weakly identified; original evidence is classification. | Diagnostic only. |
| [Meta-SGD](https://arxiv.org/abs/1707.09835) | Initialization plus per-parameter update direction and step size. | Use only on a small adapter; a global full-model step vector does not solve routing. | Optional constrained baseline. |
| [LEO](https://openreview.net/forum?id=BJgklhAcK7) | Amortized low-dimensional latent initialization followed by latent-gradient adaptation. | Supports the hybrid low-dimensional principle, but does not select donor targets or abstain. | Mechanism precedent and baseline. |
| [ProtoNet](https://proceedings.neurips.cc/paper_files/paper/2017/hash/cb8da6767461f2812ae4290eac7cbc42-Abstract.html) | Shared embedding plus class-mean prototypes; no gradient inner loop. | Literal prototypes require classification; pKi binning would discard continuous ordering. | Metric/kNN control only. |
| [TADAM](https://proceedings.neurips.cc/paper_files/paper/2018/hash/66808e327dc79d135ba18e051673d906-Abstract.html) | Support task embedding generates constrained FiLM metric modulation. | Low-dimensional conditioning is useful prior art, not source competence or harm control. | Conditional-modulation baseline. |
| [MetaDTA](https://openreview.net/forum?id=yzlif16IASM) | Attentive Neural Process amortizes a target affinity function and uncertainty from context. | Official quick-run evidence uses much larger context and does not establish k=1/3/5; variance is not abstention. | Domain-specific amortized baseline. |
| [ATS](https://arxiv.org/abs/2110.14057) | Query loss and support/query gradient agreement schedule meta-training tasks. | It controls training allocation, not recipient-time donor selection. | Source-task curriculum control. |
| [AdaMBind](https://doi.org/10.1038/s41467-026-70554-5) | MAML-style DTA learner plus ATS-like scheduler; k=5 and larger contexts. | It has no recipient-specific donor composition or no-transfer action and no k=1/3 evidence. | Mandatory equal-budget baseline. |
| [CNAPs](https://arxiv.org/abs/1906.07697) / [HyperMAML](https://arxiv.org/abs/2205.15745) | Context networks or hypernetworks generate task-conditioned parameters/updates. | Full parameter generation is underconstrained at k<=5 and remains source-blind unless donor identity is explicit. | Low-capacity hypernetwork control. |
| [VERSA](https://arxiv.org/abs/1805.09921) | Amortized probabilistic inference over task parameters from a variable support set. | Probabilistic adaptation is useful, but raw uncertainty is not negative-transfer control. | Probabilistic adapter baseline. |

**LITERATURE:** none of the primary sources above establishes strict,
target-disjoint natural-recipient pKi at k={1,3,5}, target-macro positive
transfer, majority-recipient benefit, scaffold/provenance closure, and explicit
negative-transfer reduction. The new work must not claim novelty for
meta-learning, low-rank adaptation, task conditioning, MoE, or abstention in
isolation.

### Cross-domain mechanism transfer

| Work | **LITERATURE:** original mechanism | **INFERENCE:** transferable part | Non-transferable part / failure condition |
| --- | --- | --- | --- |
| [Taskonomy](https://doi.org/10.1109/CVPR.2018.00391) | Directed transfer affinity and a budgeted task-transfer graph. | Learn directed source-to-recipient gain, not symmetric similarity. | Original target labels are abundant; affinity changes with support budget and chemistry. |
| [Task2Vec](https://doi.org/10.1109/ICCV.2019.00653) | Fisher task embeddings recommend representations. | Shrunk support-gradient/task-code signatures. | k<=5 regression Fisher estimates can be rank-deficient and unstable. |
| [LogME](https://proceedings.mlr.press/v139/you21b.html) | Bayesian linear-model evidence ranks frozen representations, including regression. | Support-only evidence baseline and router feature. | Evidence is not paired post-adaptation transfer gain and is fragile at n<<d. |
| [DSelect-k](https://arxiv.org/abs/2106.03760) | Differentiable cardinality-constrained expert selection. | Enforce a small donor set. | Sparsity does not identify competence and can collapse early. |
| [Negative transfer gap](https://doi.org/10.1109/CVPR.2019.01155) | Defines harm relative to the same target algorithm without source data. | Use a paired, equal-capacity fallback comparison. | High-dimensional density-ratio estimation is impossible from <=5 supports. |
| [SelectiveNet](https://proceedings.mlr.press/v97/geifman19a.html) | Joint prediction/selection with explicit coverage and validation calibration. | Report risk-coverage and full-cohort fallback performance. | Pseudo-recipient calibration gives no guarantee under natural-recipient shift. |
| [DML cross-fitting](https://doi.org/10.1111/ectj.12097) | Sample splitting and orthogonal scores isolate nuisance estimates. | Group folds by target/component and exclude each pseudo-recipient from every nuisance fit. | Root-n or causal guarantees do not automatically apply to A2S. |
| [Conservative contextual bandits](https://papers.nips.cc/paper_files/paper/2017/hash/bdc4626aa1d1df8e14d80d345b2a442d-Abstract.html) | Execute an exploratory action only when a pessimistic reward bound clears a baseline. | Make no-transfer an explicit action and gate on predicted gain lower bound. | A2S is one-shot and receives no deployment reward; no regret/safety theorem transfers. |

### Candidate comparison and Agent 0 decision

| Candidate | Core mechanism | Identifiable new information | Main risk | Decision |
| --- | --- | --- | --- | --- |
| **A2S-MAP**: abundant-to-scarce meta-adapter prior | Low-rank source residual task-code atlas; support-conditioned sparse source prior; exact low-dimensional posterior; paired-gain fallback. | Source-specific ligand reordering directions and their uncertainty, absent from affine source calibrators. | Residual atlas may repeat the failed ligand-kernel route or encode chemistry/provenance shortcuts. | **SELECTED MAIN HYPOTHESIS**, conditional on static gates. |
| **LPGR**: listwise pessimistic gain router | Fit full source residual experts; learn a listwise source ranking and lower-bound no-transfer action. | Replaces pointwise utility regression with relative gain ranking. | k<=5 still selects among hundreds of high-dimensional experts; winner's curse remains and adaptation is not low-dimensional. | Reject as main; retain only as a routing-loss control. |

**INFERENCE:** A2S-MAP is selected without a vote. It is the only candidate
that directly addresses all three diagnosed causes: it gives source knowledge
a reordering degree of freedom, constrains recipient adaptation to a small
task code, and compares the source-conditioned prior to an equal-capacity
zero-source prior before transfer. LPGR changes the loss but not the
identifiability bottleneck.

### Literature-derived deferred hypotheses (2026-08-01)

These are architecture or selection controls for a future authorized run, not
new core innovations. They do not amend the two A2S-MAP hypotheses, do not
permit training while D0 is `DATA_NOT_READY`, and cannot consume recipient,
development, or confirmation labels during selection.

| Paper | Assessment for the current DTA model | Decision |
| --- | --- | --- |
| [Qiu et al., Gated Attention for Large Language Models](https://arxiv.org/abs/2505.06708) | **HYPOTHESIS - deferred architecture control:** a query-dependent, head-specific sigmoid gate after SDPA in `LandmarkAttention` may add nonlinear, sparse selection over protein segment summaries and reduce redundant landmark mixing. Its LLM attention-sink result is not evidence for protein sequences; the only admissible claim is a source-fold improvement at equal parameter/FLOP budget. | Retain as a one-factor, source-only control after D0 and ranking-headroom gates pass. Compare ungated attention to post-SDPA head gating with fixed embeddings, target/component folds, compute, and an equal-capacity control. Report target-macro RMSE, Spearman, NLL/coverage, gate sparsity, and training failures across the registered seeds. |
| [Specformer: Spectral Graph Neural Networks Meet Transformers](https://arxiv.org/abs/2303.01028) | **LITERATURE - not admitted as a hypothesis:** Specformer learns a set-to-set spectral filter over the full Laplacian eigenvalue set. The current package supplies sequence segments and ECFP features, not a frozen molecular-contact or protein-interaction graph whose edges are valid under recipient/source closure. A target homology graph is specifically split-disjoint, so it cannot provide the proposed cross-target signal without violating the estimand. | Do not implement. Reconsider only with a pre-split, label-free per-example graph and a proof that all graph construction and spectral decomposition preserve the frozen closure contract. |
| [He et al., WMAGT drug-disease association model](https://doi.org/10.1186/s12859-024-05705-w) | **LITERATURE - not admitted as a hypothesis:** its weighted multi-aggregate GCN plus graph transformer operates on drug-drug, disease-disease, and known drug-disease edges for link prediction. Those relation types are neither current inputs nor the pKi target-transfer estimand; importing them would create a domain and leakage confound. | Exclude from the active DTA model. The generic local/global aggregation pattern is not a distinct mechanism once the required heterogeneous graph is absent. |
| [TG-NAS: Generalizable Zero-Cost Proxies with Operator Description Embedding and Graph Learning](https://github.com/YeQiao/TG-NAS) | **HYPOTHESIS - deferred architecture-selection control:** a graph surrogate over a fixed catalogue of model DAGs, operator descriptions, and source-only zero-cost statistics may discard clearly weak encoder variants before expensive fitting. It is an efficiency tool, not a DTA mechanism or performance claim. | Retain only for nested source-fold architecture ranking after D0 passes. Fit and evaluate the proxy on disjoint source targets; forbid recipient/development/confirmation labels, final metrics, and manual post-hoc catalogue changes. Require held-out source-fold rank correlation and superiority to random/budget-matched selection before it can reduce the registered architecture budget. |

**FALSIFICATION REQUIREMENT:** Neither retained hypothesis may replace the
equal-capacity zero-source-prior arm or be counted as a third core innovation.
The gated-attention effect must survive the frozen target/component folds and
the TG-NAS proxy must predict out-of-fold source architecture ranking before
either is evaluated once on a sealed recipient roster. Failure returns the
corresponding item to excluded-control status.

### Unique main model: A2S-MAP

Full name: **Abundant-to-Scarce Source-Atlas Meta-Adapter Prior**.

Only two core innovations are claimed:

1. **HYPOTHESIS - uncertainty-matched source-atlas meta-initialization.**
   Abundant targets define target-specific residual functions in a shared
   low-rank task-code atlas. A scarce recipient's support posterior selects a
   sparse, uncertainty-aware barycenter of source task codes, which becomes a
   recipient-conditioned prior. The inner update changes only the recipient
   intercept and low-dimensional task code.
2. **HYPOTHESIS - cross-fitted paired-gain meta-policy.** The outer objective
   and gain head compare the source-conditioned posterior to an identical
   zero-source-prior adapter. Transfer is used only when a target-fold
   calibrated lower gain estimate clears a frozen material margin; otherwise
   prediction falls back exactly to the zero-source-prior arm.

Low rank, Bayesian updating, sparse routing, cross-fitting, and abstention are
implementation primitives. Novelty, if supported, lies in the paired mechanism
and strict A2S estimand, not any primitive alone.

### Four responsibility zones and frozen interfaces

**Zone P - preprocessing and task construction.** Inputs are immutable source
files plus publisher metadata; outputs are disjoint `D_source` and frozen
`{S_r^k,Q_r}` rosters, component/split manifests, episode hashes, and leakage
reports. It owns compound/target/document/assay/source canonicalization,
endpoint/unit/censoring eligibility, parent/scaffold/chemical-neighbour closure,
historical time closure, label-blind nested support sampling, and transforms.
Any fitted normalization, vocabulary, scaffold rule, or feature transform is
fit on source-training data only and then frozen. Recipient query labels and
query-derived selection statistics are forbidden. Zone P currently remains
`DATA_NOT_READY`; no downstream zone may fit until H0-S and H1-S pass.

**Zone M - shared representation and meta-training.** Inputs are only Zone P's
abundant-source training folds and source-fold pseudo-recipient episodes.
Outputs are the frozen ligand features/base predictor `f0`, residual basis `U`,
cross-fitted source codes `{z_h,Sigma_h}`, source reliability nuisances, prior
precision, and paired-gain meta-policy parameters. Meta-parameters are
`theta_M={f0,U,source-code metric,prior precision}`; pseudo-task parameters are
an intercept plus q-dimensional code. Outer objectives use only held-out source
pseudo-query labels. Natural-recipient identities, labels, query metrics, and
test-component summaries are forbidden, and all baselines receive identical
episodes and compute accounting.

**Zone C - transfer control and negative-transfer suppression.** Inputs are a
recipient target descriptor, its declared support features/labels, frozen
source-atlas summaries, and source-out-of-fold reliability/provenance features.
Outputs are support task state `(mu_r^0,Sigma_r^0)`, sparse donor weights
`pi_rh`, source-conditioned prior `(m_r,Sigma_r^T)`, predicted lower paired
gain, and a frozen transfer/fallback action. Learnable parameters are the
source-code metric, low-capacity routing coefficients, and cross-fitted gain
quantile head; none are updated on natural recipients. Query labels and query
metrics are forbidden. Primary evaluation is inductive and does not read query
features in Zone C; any future transductive variant requires a separate
preregistration and report. Random, wrong, chemistry-only, protein-only,
provenance-only, and always-transfer controls are mandatory.

**Zone I - recipient adaptation and inference.** Inputs are the frozen Zone M
artifacts, Zone C prior/action, recipient target, exactly k nested support
pairs, and query compound features. It adapts only one intercept plus q=4 code
coordinates through one exact posterior solve; no encoder/full expert gradient,
early stopping, recipient-specific step choice, or query-driven checkpoint is
allowed. All queries for one recipient/support draw share the same adapted
state. The same frozen model and protocol serve k={1,3,5}. Outputs are both
potential predictions (zero-source and source-prior), the chosen action,
uncertainty, and final pKi/ranking scores. Query labels enter only the final
Zone-I evaluation artifact and can never flow back to Zones P/M/C or adaptation.

The only permitted cross-zone direction is `P -> M -> C -> I`, with frozen Zone
P artifacts also supplied read-only to C/I for declared target/support/query
features. Natural query outcomes flow to the final evaluator only. Any reverse
flow, hidden recipient statistic, or query-dependent model selection is a
leakage STOP.

### Mathematical definition

Let `phi(d)` be frozen ligand features and let `f0(d)` be a target-balanced
pooled source predictor. On abundant target `h`, fit the residual model

```text
y_hi = f0(d_hi) + b_h + g_U(d_hi)^T z_h + epsilon_hi
g_U(d) = U^T phi(d),  z_h in R^q.
```

`U` is a shared residual basis and `z_h` is the source task code. The primary
probe uses `q=4`; `q in {2,4}` is the only source-fold model-selection range.
This residual term can change ligand ordering, unlike the stopped affine
source adapters.

For recipient support `S_r`, compute a zero-source-prior support posterior
over the code and intercept. In schematic Gaussian form,

```text
Sigma_r^0 = (Lambda_0 + G_S^T G_S / sigma^2)^-1
mu_r^0    = Sigma_r^0 G_S^T (y_S - f0_S - b_r) / sigma^2.
```

For each candidate source, combine recipient and source code uncertainty:

```text
s_rh = -0.5 (mu_r^0-z_h)^T (Sigma_r^0+Sigma_h)^-1 (mu_r^0-z_h)
       -0.5 logdet(Sigma_r^0+Sigma_h) + beta^T q_rh.
```

`q_rh` contains only source-out-of-fold reliability, protein relation,
support/source chemistry coverage, and provenance compatibility. Select the
top `J` candidates, with `J in {4,8}` chosen only inside grouped source folds,
and compute sparse weights `pi_rh`. The source-conditioned initialization is

```text
m_r = sum_h pi_rh z_h
Sigma_r^T = (Lambda_T + G_S^T G_S / sigma^2)^-1
mu_r^T = Sigma_r^T (Lambda_T m_r + G_S^T (y_S-f0_S-b_r)/sigma^2).
```

The proposed and fallback predictions differ only in the source prior:

```text
yhat_r^T(d) = f0(d) + b_r + g_U(d)^T mu_r^T
yhat_r^0(d) = f0(d) + b_r + g_U(d)^T mu_r^0.
```

For a held-out pseudo-recipient, define the paired source-transfer gain

```text
G_r(k) = L_Q(yhat_r^0) - L_Q(yhat_r^T).
```

The target-macro meta-objective is query loss plus a one-sided penalty on
`G_r(k)` below the frozen material margin. A separate low-capacity quantile
head predicts a lower gain quantile from support-only context, posterior
entropy, source-weight concentration, and source-fold reliability. The test
action is

```text
a_r = 1 if calibrated_lower_gain_r > delta_material else 0
yhat_r = a_r yhat_r^T + (1-a_r) yhat_r^0.
```

**FACT:** natural recipient query labels are not an input to `U`, `z_h`, source
weights, hyperparameters, the gain head, its threshold, or early stopping.
They are used once for final evaluation. **INFERENCE:** cross-fitted calibration
does not guarantee natural-recipient safety under distribution shift; this is
why coverage and natural-tail performance are empirical gates.

### Meta-knowledge, meta-training, and meta-test

**HYPOTHESIS:** A2S-MAP meta-knowledge comprises the residual basis `U`, source
task codes and uncertainty, the source-code metric, source prior precision,
and the support-only gain quantile head. It is not a single MAML initialization.

Meta-training protocol:

1. Use pKi only. Partition abundant source targets by homology/provenance
   component into five deterministic outer folds.
2. For each fold, refit `f0`, `U`, source codes, normalizers, reliability, and
   every routing nuisance on the other folds. The held-out target is never in
   its own atlas or nuisance fit.
3. Construct nested, label-blind `S^1 subset S^3 subset S^5` supports. Truncate
   or resample pseudo queries to the frozen natural query-depth distribution;
   retain all episode hashes and row IDs.
4. Use held-out pseudo query labels only for meta-objective and gain labels.
   Fit routing parameters on meta-train folds and the gain threshold on a
   disjoint calibration fold. Rotate folds; aggregate target-macro.
5. Fit the final source atlas on all abundant sources only after all choices
   are frozen. Natural scarce targets contribute no labels to this fit.

Meta-test protocol:

1. Expose only the declared recipient support labels and allowed covariates.
2. Infer `b_r, mu_r^0, Sigma_r^0`, select at most J source codes, and compute
   `mu_r^T` analytically. No encoder or full expert is updated.
3. Apply the frozen gain policy. If it abstains, return `yhat_r^0`, not a
   missing prediction.
4. Predict the untouched query once. Store both potential predictions and the
   chosen action for audit, but never feed query outcomes back into the model.

Inner-loop range: the primary model adapts one intercept plus `q=4` task-code
coordinates by one exact posterior solve; no SGD step is used. `q=2` is the
only dimension ablation. Full-model, Transformer, Mamba, protein encoder, and
ligand encoder updates are forbidden at this stage.

### Why the performance mechanism is testable

```text
source-specific residual atlas
 -> query-reordering headroom beyond affine calibration
 -> support identifies a low-dimensional recipient task code
 -> uncertainty-aware sparse source prior reduces k<=5 variance
 -> paired gain policy rejects harmful priors
 -> lower negative-transfer rate
 -> better natural-recipient RMSE and within-target ranking.
```

The expected order of intermediate improvement is frozen:

1. cross-fitted residual-atlas query R2 and ranking over affine `B0`;
2. support-code recovery and top-J positive-source recall;
3. calibration of predicted versus realized paired gain;
4. negative-transfer rate at fixed coverage;
5. natural-recipient RMSE, Spearman, pairwise accuracy, and NDCG@10.

If an earlier variable fails, a later aggregate gain cannot be attributed to
the proposed mechanism.

### Essential differences from alternatives

- **INFERENCE - fine-tuning:** pooled fine-tuning has no donor identity or
  explicit fallback and usually updates more parameters. A2S-MAP adapts at
  most five scalars after a source-conditioned prior.
- **INFERENCE - global MAML/FOMAML/ANIL/Meta-SGD:** these learn one global
  initialization or update rule. A2S-MAP retains a source task-code atlas and
  creates a different prior for each recipient.
- **INFERENCE - AdaMBind/ATS:** they schedule meta-training tasks. A2S-MAP
  selects source knowledge at recipient inference and judges paired gain
  against an equal-capacity fallback.
- **INFERENCE - ordinary MoE/DSelect-k:** ordinary gates choose experts from
  context but do not define source competence as cross-fitted paired gain.
  A2S-MAP weights source task codes, then performs a support posterior update;
  it does not average full networks.
- **INFERENCE - MetaDTA/CNAPs/HyperMAML/VERSA:** these amortize a task predictor
  or parameters globally. A2S-MAP exposes donor identity and donor uncertainty
  and can refuse donor-specific information.
- **INFERENCE - stopped A2S-CFRA:** A2S-CFRA regresses pointwise utility for
  affine experts and chooses an argmax. A2S-MAP first supplies a residual
  reordering basis, selects a small uncertainty-matched code set, optimizes the
  final adapted query loss, and uses a separately calibrated paired-gain arm.
- **INFERENCE - prior Bayesian failures:** the stopped Bayesian routes used a
  global/protein-conditioned support kernel on a different strict task and
  failed ligand/gradient controls. A2S-MAP is reopened only conditionally and
  must first show source-atlas headroom; a Bayesian solve is not evidence by
  itself.

### Data contract before any model fit

- pKi is primary. pKd is a separately constructed, separately fit, separately
  tested secondary replication and cannot rescue pKi.
- Source and recipient target IDs are disjoint. Source/pseudo-recipient folds
  are grouped by homology and provenance component.
- Natural recipients are defined by frozen resource and temporal/provenance
  rules, not random row deletion. Support precedes or is source-disjoint from
  query; exact parent connectivity and document/assay/source units do not cross.
- The common primary roster must support all k values with nested support and
  at least 10 independent query units after closure. NDCG@10 is undefined and
  not imputed when fewer than 10 remain.
- Primary target-side single-cold may allow query compounds seen on source
  targets, but exact overlap, nearest-neighbour similarity, scaffold-cold,
  global drug-cold, homology-cold, and provenance strata are reported
  separately. No broad scaffold-cold claim may use the single-cold result.
- Use at least five frozen, label-blind nested support draws per recipient;
  average draws and seeds within recipient before inference. Draws, rows,
  pairs, episodes, queries, folds, and seeds are not biological replicates.
- Construct independent bootstrap components by joining recipients that share
  a homology component or held-out provenance/source family. Freeze the graph
  before model scores.

### Minimum kill tests and ordered experiment plan

**Gate D0 - unique next action, no model fit:** build the pKi natural-tail
document/time/source/parent-closed common roster and report target/component
count, query-depth distribution, support/query and source/recipient overlap,
dependency concentration, and MDE80. STOP if fewer than 50 recipients or the
component MDE exceeds the material effect that the study can detect.

**Gate D0 result - DATA STOP:** the metadata-only audit found 40 recipients in
the most optimistic provenance-closed upper bound, 34 after support-scaffold
diversity, 22 in the ChEMBL-release temporal envelope, 11 after temporal and
scaffold/provenance closure, and zero in the strict source-family-closed
roster. True publication/measurement time is absent, all candidates form one
source/homology dependency component, and the strict MDE80 is undefined. S0
and every affinity-model training stage are therefore blocked.

**Gate S0 - frozen-feature atlas headroom:** after D0 PASS, fit source-fold
rank-2/4 residual SVD/ridge atlases on frozen 1,034-D ligand features. On
held-out pseudo-recipients with natural-matched queries, compare affine `B0`,
zero-prior code adaptation, query-oracle code upper bound, and support-selected
source prior. PASS only if the oracle and deployable support code both add
ranking headroom at k=3/5 and source-code shuffle destroys it. This is the
cheapest kill for repeating the old ligand-kernel failure.

**Gate K1 - one-seed mechanism kill:** seed 1729, frozen features, q=4, J=8,
one exact inner solve. Compare zero prior, source prior, source prior plus gain
policy, and matched random policy. Require positive pKi gain at k=3 and k=5,
source prior above zero prior, at least a 10 percentage-point negative-transfer
reduction, coverage at least 0.60, and mechanism destruction under source-code
shuffle. Failure stops the model; do not tune width, rank, J, losses, or epochs.

Conditional formal order after K1 PASS:

1. execute equal-budget recipient calibration, ridge/RF/kNN, pooled
   fine-tuning, MAML/FOMAML/ANIL, and AdaMBind;
2. execute source-routing controls and the two core ablations;
3. execute chemistry/protein/provenance destruction controls;
4. run five seeds only after the one-seed gate passes;
5. run independent pKd replication only after pKi mechanism gates pass;
6. consider encoder expansion only after natural pKi and mechanism attribution
   both pass. Encoder expansion is a new preregistration, not an automatic step.

### Strong baselines

Every baseline uses identical recipients, support/query rows, target weights,
frozen encoder inputs, and permitted label access. Record trainable parameters,
gradient evaluations, support use, episode count, wall time, peak CUDA memory,
and cross-fitting status.

1. recipient calibration-only and global `B0`;
2. recipient-only ridge and RF, kNN-DTA, and LogME/evidence model selection;
3. target-balanced pooled model plus equal-parameter fine-tuning;
4. exact MAML, FOMAML, ANIL, low-dimensional Meta-SGD, and optional Reptile;
5. MetaDTA/VERSA-style amortized adapter and constrained conditional modulation;
6. equal-budget AdaMBind and ATS-only scheduler control;
7. global zero-prior task-code adapter with the same `U`, q, and inner solve;
8. all-source, random-source, protein-only, chemistry-only, provenance-only,
   current support-compatible affine routing, and stopped linear router;
9. LPGR listwise router, ordinary sparse MoE, and Mera-style source weighting;
10. A2S-MAP without and with the paired-gain policy.

### Ablations and destructive controls

The two core innovation ablations are mandatory:

1. replace the recipient-conditioned source-code prior by the zero/global prior;
2. remove the gain policy and always transfer at identical capacity.

Destructive controls are source-code identity shuffle; source reliability
shuffle; protein relation shuffle; support chemistry shuffle; support-label
permutation; wrong-recipient support; no-label recipient; source-mean code;
chemistry-only nearest neighbour; provenance-only routing; uniform all-source
prior; random source prior; matched random abstention; and high-abstention
stress. An oracle using pseudo query labels is reported only as headroom and is
never a deployable arm.

PASS requires source-code shuffle, support-label permutation, and wrong-support
to remove the proposed gain. If chemistry-only or provenance-only matches the
full model, the mechanism is a shortcut. If always-transfer matches the gain
policy at equal coverage, the negative-transfer innovation is not identified.

### Primary metrics and statistical analysis

For recipient `r` and budget `k`, the primary paired estimand is

```text
Delta_rk = RMSE_r(strongest frozen no-routing fallback)
            - RMSE_r(A2S-MAP policy).
```

Primary budgets are k=3 and k=5, fixed in advance; k=1 is fully reported but
does not rescue failure at k=3/5. Report target-macro RMSE and MAE, within-target
Spearman, pairwise accuracy, NDCG@10, gain AULC on the common roster, median
gain, benefiting-recipient fraction, negative-transfer rate, coverage,
risk-coverage AUC, calibration of predicted gain, and leave-one-component-out
influence.

Average support draws and seeds inside recipient first. Bootstrap 5,000 times
over the frozen independent target/provenance components, paired across arms.
Rows, pairs, episodes, support draws, queries, folds, and seeds are never
resampled as independent biological units. Use the same component draws for
all primary contrasts.

### Excellence PASS/STOP rules

Full PASS requires all conditions below:

- natural-tail pKi target-macro RMSE gain at both k=3 and k=5 exceeds
  `max(MDE80, 0.05 pKi units)` and each paired component-bootstrap 95% lower
  bound is above zero;
- A2S-MAP beats the strongest no-routing baseline, pooled fine-tuning, exact
  MAML, equal-budget AdaMBind, global zero-prior adapter, random/all-source,
  protein-only, chemistry-only, and current affine support routing;
- median gain is positive and more than 50% of recipients benefit;
- the gain policy reduces negative-transfer rate by at least 10 percentage
  points versus always-transfer while retaining coverage >=0.60 and improving
  full-cohort fallback metrics;
- at k=3 or k=5, at least one prespecified ranking metric has a positive paired
  95% lower bound, and no primary ranking metric materially regresses;
- positive gain remains in the prespecified middle/low chemical-similarity or
  scaffold-cold stratum and is not confined to one homology/provenance family;
- pseudo-recipient and natural-recipient directions agree; pseudo-only success
  is a STOP, not a partial natural-tail claim;
- five-seed estimates are stable after the one-seed gate, with seeds averaged
  within target and no post-hoc seed, k, recipient, or stratum selection.

STOP for failed D0 power/closure; no residual-atlas ranking headroom; global
prior equivalence; source-code-shuffle invariance; learned routing no better
than random/all-source; gate coverage below 0.60; no negative-transfer
reduction; chemistry/provenance shortcut; majority-recipient harm; pseudo-only
gain; query-label access; or any source/recipient, parent, document, assay, or
declared chemical closure violation. A STOP cannot be rescued by capacity,
epochs, extra seeds, wider rank, a new unregistered loss, relaxed split, pKd,
confirmation, or post-hoc thresholding.

### Compute budget and execution environment

- All numerical probes must use
  `D:\anaconda\envs\drug\python.exe` and CUDA. Metadata-only D0 does not
  require CUDA and must not fit an affinity model.
- S0 frozen-feature probe: at most 1 GPU-hour, 2 GiB peak Torch allocation,
  zero neural-encoder updates, seed 1729.
- K1 one-seed kill: at most 4 GPU-hours and 4 GiB peak Torch allocation.
- Conditional five-seed round: seeds `{1729,1731,1733,1741,1753}`, at most 20
  aggregate GPU-hours. Stop the job and record failure if a bound is exceeded.
- Report GPU name, utilization, power, wall time, peak NVIDIA memory, peak
  Torch memory, parameters, and gradient evaluations. Parameter padding is not
  a compute control.

### Failure Ledger

Existing historical failures remain binding. New entries are append-only.

| Date | Scientific hypothesis / command and seed | Data role | Observation | Failure type and alternative explanation | Decision / evidence required to reopen |
| --- | --- | --- | --- | --- | --- |
| 2026-07-31 | Initial A2S pKi/pKd baseline commands; seed 1729 | TRAIN source and provisional scarce recipients | Filtered local row IDs were applied to the global ligand-feature cache. | Software/data alignment failure; all apparent model effects are uninterpretable. | **STOP** and quarantine. Reopen only with global row IDs, feature length and connectivity hash verification; corrected artifacts now satisfy this. |
| 2026-07-31 | `D:\anaconda\envs\drug\python.exe main.py a2s-router --endpoint pKi --out reports/active/a2s_router_pki_targetbalanced_seed1729.json`; seed 1729 | Source-fold pseudo-recipients -> count-defined scarce recipients | Gain `-1.144/-1.213/-1.315`; negative-transfer `0.801/0.774/0.821`; gate accepts all; source collapse. | Mechanism failure. Alternatives are pseudo-natural shift, depth/provenance shortcut, pointwise utility miscalibration, and argmax winner's curse. | **STOP** A2S-CFRA. Reopen only as a control after a closed natural roster and a source residual atlas with static headroom. |
| 2026-07-31 | Read-only `build_episodes` metadata reconstruction; no seed, no model fit | Provisional pKi recipient episodes | Most episodes share document/assay with query; many share scaffold; pseudo query median is about 24x natural. | Estimand/closure failure for a natural-tail interpretation. The existing result remains a weak single-cold upper control. | **STOP** natural-tail claims. Reopen only after D0 document/time/source/parent closure, common roster, dependency graph, and MDE PASS. |
| 2026-07-31 | `D:\anaconda\envs\drug\python.exe main.py natural-tail-audit --out dataset/processed/a2s_natural_tail_d0.v1.json`; no seed-dependent scores, no model fit | TRAIN pKi metadata only | Optimistic provenance upper bound 40 (<50); temporal/scaffold closure 11; strict source-family closure 0 recipients/0 components; MDE80 undefined. | Data topology and provenance failure. `chembl_release` is only an ingestion proxy, and shared source family joins all candidates into one component. | **DATA STOP.** Reopen only with a hashed provenance-rich pKi source containing true time/lineage metadata and at least 50 recipients that pass the identical D0 contract. |
| 2026-07-31 | Post-run adversarial code review of Gate D0; no model fit or query-label read | D0 implementation and emitted certificate | Strict diagnostics add an unregistered scaffold-cold constraint, may repeat a k=5 support set across draws, use provisional resource thresholds, and do not recompute emitted overlap. | Protocol/certificate limitation. It does not reverse the DATA stop because the separate exhaustive five-distinct-support parent/document/assay upper bound is only 40 and true time/lineage metadata are absent. | **STOP** reuse of the strict roster. Reopen D0 only after preregistering thresholds and correcting draw uniqueness, primary closure, overlap, component, and immutable-certificate checks. |
| 2026-07-31 | `D:\anaconda\envs\drug\python.exe main.py h0-inventory --acquire-remote-metadata`; no model fit or numeric label read | H0-S ChEMBL/BindingDB/GtoPdb metadata | ChEMBL publication year is available for 9,640/9,643 documents and document-release date for all, but the final registry has no activity-row first-seen field and historical 24.1/27/31 databases are absent. | Data-time identifiability failure. Document release can predate a later backfilled activity; BindingDB/GtoPdb are later HIST-L candidates, not HIST-S substitutes. | **DATA STOP.** Reopen only after checksum-verified historical ChEMBL snapshots yield an affinity-blind stable activity-identity projection for all three index dates. |
| 2026-07-31 | Multiple pre-gate curl processes wrote `chembl_24_1_sqlite.tar.gz`; exact complete command set unavailable retrospectively | H0-S ChEMBL 24.1 archive acquisition | Final size matched 3,659,492,620 bytes but SHA-256 was `46937b804020de62714fd7791de1a7b2d303652c094009e1b2bebfdf42681922`, not the publisher digest. Bytes remain quarantined under the digest-bearing `.corrupt` name. | Data acquisition/concurrency failure; same-path concurrent resume/write made the archive uninterpretable. | **STOP** use of those bytes. Reopen only through a clean staged download with official size/SHA and an immutable run record. |
| 2026-07-31 | Three ChEMBL 31 curl resume attempts; contemporaneous cells and PIDs were intentionally terminated | H0-S ChEMBL 31 archive acquisition | The surviving partial is 4,242,055,168 of 4,505,413,744 bytes and has no contemporaneous digest. | Procedural interruption to migrate execution under the dataset-processing provenance gate; not evidence of archive corruption. | **CONTINUE** only after preserving the partial in a retrospective record and resuming or restarting through a staged, hash-validated logged run. |
| 2026-07-31 | First 24.1 `historical-project` attempt; PID 42148 later stopped | H0-S primary identity projection | Interrupted before any activity/document Parquet or manifest was published. | Implementation/protocol correction: replay redundantly repeated the full integrity scan. | **STOP** that run. Reopen only with primary and replay in separate logged processes. |
| 2026-07-31 | Second 24.1 `historical-project` attempt; exact PID unavailable | H0-S primary identity projection | Interrupted after the complete-processing-log requirement arrived; projection output directory remained empty. | Provenance-gate interruption, not a scientific result. | **STOP** that run. Reopen only through `dataset-run` with prepared/log/run/seal/ledger evidence. |
| 2026-07-31 | Logged `download-chembl37`; run `20260731T125817198299Z_download-chembl37_29142fe9`; no seed/model | H0-S ChEMBL 37 archive acquisition | curl exit 56 after 1,494,429,696 bytes; schannel reported that the server closed without `close_notify`. Expected size/SHA did not pass, canonical publish was not attempted, and the partial SHA-256 is `520e2c4d46f3cbef83505b00b32f185d9e750a427c9e89aee726c7b0bef25c22`. | Network/transport interruption, not evidence of corrupt final content. | **CONTINUE** only as a logged resume whose `parent_run_id` names this failed run; official final size and SHA remain unchanged. |
| 2026-07-31 | Logged resume `20260731T130610857908Z_download-chembl37-resume-1_c8bfda37`, then logged control `20260731T131131874669Z_terminate-c8bfda37_7c2d1562` | H0-S ChEMBL 37 archive acquisition | The resumed connection stayed below about 25 KiB/s with ETA above 58 hours and was deliberately terminated; child exit 15. Partial grew to 1,502,298,112 bytes, SHA-256 `d130d58cdfe1ec45db05d2f7d02e14b2a0fe1b0c8b1768aed3c1fdcf6d0c565f`. Editing the runner to create the logged termination action also correctly triggered `code_unchanged=false`. | Operational throughput STOP plus expected reproducibility failure for that run; no canonical publish. | **CONTINUE** only with this failed run as parent and a preregistered low-speed abort/retry policy; final publisher size/SHA remain mandatory. |
| 2026-07-31 | Logged `project-chembl241-primary`; run `20260731T133335519650Z_project-chembl241-primary_f7a55599` | H0-S ChEMBL 24.1 identity projection | After the full integrity scan, code rejected the real SQLite version row `('ChEMBL_24','2018-04-23 ...')` because it inferred that the internal name must contain archive revision `24.1`. Exit 1; all three declared outputs remained missing. | Implementation/version-contract bug. ChEMBL 24.1 is a revised release archive whose SQLite internal major version is exactly 24; this is not evidence of database corruption. | **CONTINUE** only after registering exact per-release internal version names, testing 24.1 acceptance and wrong-name rejection, and rerunning the full primary projection. |
| 2026-07-31 | Concurrent `project-chembl27-primary-v1` runs `20260731T152113745605Z_project-chembl27-primary-v1_0f428010` and `20260731T152138390894Z_project-chembl27-primary-v1_3ff76781`; no seed/model | H0-S ChEMBL 27 identity projection | Two runners initially targeted the same three canonical outputs. The later run was terminated before any output existed and indexed `FAILED`; the retained run completed successfully with activity/document identity outputs. | Dataset-processing concurrency/provenance failure; prior runner lacked an atomic output reservation. | **CONTINUE** only after adding sealed `output_claim.json`, identity-aware orphan reconciliation, and a passing same-output concurrency test; then independently replay the retained projection. |
| 2026-08-01 | `project-chembl37-primary-v1` dataset-run wrapper; no seed/model | H0-S ChEMBL 37 identity projection | The outer PowerShell tool reached its 1-hour timeout while runner `36100` and child `46976` remained active. The existing run later completed with exit 0, all outputs validated, and claim/seal/ledger evidence present; duration was 3658.30 s. | Operational monitoring timeout, not a dataset or scientific failure; the child process was preserved and independently verified. | **CONTINUE** with the sealed successful run; future long operations must use a detached monitor whose timeout cannot obscure runner completion. |
| 2026-08-01 | `h0-s-metadata-inventory-v3`; run `20260731T181218233594Z_h0-s-metadata-inventory-v3_a850289a`; no seed/model | H0-S four historical ChEMBL snapshots plus cross-release presence | All four projection and replay certificates passed and the outcome firewall read/materialized zero affinity values. Presence replay passed, but the presence manifest reports `12782` hard native collisions and `identity_status=false`; H0 decision is `DATA_NOT_READY`. | Identity/provenance closure failure. Native IDs are not sufficient stable measurement identities across releases; collision resolution or a stronger immutable lineage key is required. | **DATA STOP.** Do not build H1-S roster or fit any model. Reopen only after a versioned, independently replayed collision-resolution audit makes presence identity pass without reading affinity values. |

### Current status and unique next action

### Model input contract (2026-08-01)

All preprocessing is consolidated in `scripts/preprocess.py`. It reads
only the audited ChEMBL-37 dual-cold registry, frozen strict pKi roster, and
deterministic feature caches, then emits the immutable package
`dataset/processed/a2s_validation_small.v1/`. Models may read only that
package's manifest and listed arrays/tables; raw registries and databases are
forbidden at model time. Source-only normalization and nested k={1,3,5}
support are enforced, and the package is explicitly development-only while H0
remains blocked. The preparation run and hashes are sealed in
`dataset/processing_history/v1`.

The model-facing package is also copied to the separate `dataset/ready/`
directory. `dataset/ready/manifest.json` is the authoritative index; raw
inputs remain under `dataset/public/` and processing intermediates remain under
`dataset/processed/`. Ready packages are versioned and must never be deleted
or overwritten in place.

Formal training data admission is governed by
`dataset/ready/FORMAL_TRAINING_DATA_SPEC.md`; candidate collection is tracked
in `dataset/ready/formal_training_candidates.json`. Collection does not imply
training authorization.

Current final state: **`DATA_NOT_READY`**.

The H0-S identity-collision route (resolve or quarantine the 12,782 hard native
collisions across the four sealed ChEMBL projections) is **superseded as the
unique next action** by the 2026-08-01 section below. It remains the only route
to a *database-release-time* estimand and is retained as a deferred track; it is
no longer a prerequisite for the primary scientific claim. Do not fit A2S-MAP,
rerun the stopped router, run MAML/AdaMBind, relax closure, lower the target
floor, or change an encoder until Gate D0-R passes.

## 2026-08-01 A2S-SDO Preregistration And D0 Protocol Revision

This section is registered before any affinity value is read by a model, before
any fit, and before any probe. It supersedes the model-selection portion of the
2026-07-31 Meta-Learning Committee Preregistration. It does not erase, reopen,
or reinterpret any recorded failure. Evidence labels are as previously defined.

### Corrected diagnosis: what is and is not identified

**FACT** — from `reports/active/a2s_pki_targetbalanced_seed1729.json`, the
within-target Spearman of the pooled ridge `B0` on scarce pKi recipients is
`0.0824 / 0.0937 / 0.0829` at `k=1/3/5`, and
`recipient_calibration/spearman` is **bit-identical to** `b0/spearman` at every
`k`. A per-target additive offset cannot reorder compounds; this identity is a
structural check, not a result.

**FACT** — in the same artifact every routing arm sits within `±0.02` of `B0`
on Spearman: `source_support` `0.0878/0.0897/0.0776`, `source_protein`
`0.0846/0.1001/0.1055`, `source_chemistry` `0.0769/0.0871/0.0811`,
`source_random` `0.0792/0.0694/0.0816`. Source routing does not move within-target
ranking in any direction.

**FACT** — from `reports/active/anchordelta_retrain_precheck_decision_2026-07-31.md`,
a trained support-relative antisymmetric comparator on frozen pair features
reaches RMSE `1.0775` and Spearman `0.2521` where recipient calibration reaches
`1.2402` and `-0.0026`. The wrong-protein arm reaches `1.0775` and `0.2522`;
the correct-minus-wrong-protein component-bootstrap intervals cross zero on
RMSE, Spearman, and pairwise accuracy. Wrong-target support labels leave
ranking unchanged at `0.2521` while degrading RMSE to `1.5173`.

**INFERENCE — the identified decomposition.** For scarce recipient `t` and
compound `d`, prediction decomposes into three terms whose identifiability
status is now settled by project evidence:

```text
yhat(t,d) = A(t)            target level/scale        IDENTIFIED from support labels
          + R(d)            target-independent        IDENTIFIED (comparator, Spearman ~0.25;
                            ligand potency ordering   pointwise ridge only ~0.085)
          + E(t,d)          target-conditioned        NOT IDENTIFIED in any framing tried
                            reordering residual
```

**INFERENCE — prohibited reinterpretation.** `E(t,d)` has now been the target of
at least five independent framings — DCST interaction transport, the
flexible-kernel finite-rank posterior, the IDG-RBP correct-protein probe,
BridgeFIRE/Gate-P physical structure, and protein-conditioned AnchorDelta — and
has failed each time, three of them with an explicit wrong-protein or
protein-free control that matched the correct-protein arm. The preregistered
A2S-MAP core hypothesis, `g_U(d)^T z_h` as a *source-target-specific residual
reordering direction*, is a sixth framing of the same unidentified term. It may
not be described as a new mechanism class, and its failure may not be
attributed to capacity, rank, `J`, or episode count.

**INFERENCE — the reframe.** The A2S estimand does not require transferred
knowledge to be protein-specific. A comparator meta-learned on abundant targets
and applied to a disjoint scarce recipient *is* abundant-to-scarce transfer even
when its parameters carry no protein conditioning. The program has been
discarding components that work (`R(d)`) because they fail a protein-conditioning
test that the A2S estimand does not impose. The correct no-transfer control is
therefore recipient-only calibration, **not** `B0`, because `B0` already contains
abundant-source knowledge.

### Corrected diagnosis: why Gate D0 failed

**FACT** — the D0 stop, the historical amendment, and the v4 formal package all
evaluate natural scarcity at a *global* index date `tau` drawn from the
preregistered grid `{2018-12-31, 2020-12-31, 2022-12-31}`. The v4 sealed report
`dataset/formal_training/chembl37_pki_formal.v4/reports/natural_tail_d0.json`
records `1` recipient at 2022 and `5` at 2018.

**FACT — outcome-blind metadata sweep, 2026-08-01, non-admitted diagnostic.**
A read-only sweep of `canonical/pki_measurements_exact.parquet` restricted to
`target_uid, accession, protein_class_id, compound_parent_uid,
connectivity_inchikey, document_uid, document_year, assay_context_uid,
document_src_id, assay_id` (no `pKi`, `standard_value`, or `pchembl_value`
column was read) establishes:

1. The index-date grid was never swept below 2018. Recipient counts under the
   registered band (`5-29` pre-`tau` parents, `>=10` post-`tau` parents) are
   `69/67/64/57/46/42` at `tau = 2008/2010/2012/2013/2014/2015`, against `6` at
   2018 and `2` at 2022. The frozen floor of 50 is cleared at
   `tau <= 2013`.
2. The binding constraint at late `tau` is post-`tau` query depth, not the
   scarcity band. At `tau=2018`, `254` targets sit in the `5-29` band but only
   `6` have `>=10` post-`tau` parents and `150` have zero post-`tau` rows.
   ChEMBL-37 Ki coverage thins in its own curation-lag window; a global `tau`
   inside that window cannot produce recipients.
3. `document_src_id` is a curation-channel label, not a provenance family:
   `135,775` of `157,613` corpus rows carry `src_id=1` (scientific literature).
   The 2026-07-31 D0 dependency graph, which joined all 193 candidates into one
   component through shared `src_id`, used a degenerate dependency definition.

**FACT — within-recipient document-ordered feasibility, same sweep.** Under a
per-recipient design where support is drawn from a target's earliest documents
and query from strictly later document years, with support/query parent- and
document-disjoint by construction, the corpus yields, at `k=5` support parents
and `>=10` query parents: `31` recipients for a `[5,29]` total-parent scarcity
band, `82` for `[5,50]`, `106` for `[10,60]`. For the `[5,50]` roster:
`82` recipients across `82` distinct accessions and `58` protein classes,
largest protein-class share `0.073` and top-3 share `0.171`; query depth median
`22.5`, `q25 = 16`, min `11`; support pool median `9`, min `5`; split years
`1981-2020`, median `2009`. Against the `>=100`-parent source pool of `214`
recipient-excluded targets: `0` shared `target_uid`, `0` shared `accession`,
`0` shared `assay_id`, `169` shared documents and `756` shared parent compounds
requiring closure, and `559 / 1698` recipient query connectivity keys
(`0.329`) also present on source targets.

**INFERENCE.** A natural-tail roster above the frozen floor of 50 recipients,
with adequate query depth and low dependency concentration, exists inside the
already-sealed ChEMBL-37 corpus. The D0 DATA STOP is attributable to (a) an
unswept global index-date grid placed in the curation-lag window, (b) a global
rather than per-recipient time origin, and (c) a degenerate `src_id`-based
dependency definition — not to an absence of natural-tail data.

**INFERENCE — the H0-S route is not on the critical path.** Historical database
snapshots and the 12,782-collision identity audit establish *database-release*
time. The primary scientific claim is prospective *publication*-ordered
extrapolation, whose correct clock is `document_year`, already present in
ChEMBL-37. Release time is required only for a deployment-replay estimand that
the primary claim does not make. This does not relax the recorded H0-S facts.

**Prohibited.** The sweep above is a scratchpad diagnostic. It is not an
admitted artifact, not a roster, and may not be cited as evidence in any result.
An admitted roster must be produced through `main.py dataset-run` under Gate
D0-R below.

### Revised estimand (Gate D0-R)

For scarce recipient `r`, let `D_r` be its documents ordered by
`(document_year, document_uid)` and let `tau_r` be a per-recipient split year.
Support `S_r` is drawn label-blind from documents with year `<= tau_r`; query
`Q_r` is the parent compounds appearing only in documents with year `> tau_r`
and not in any support document. Sources are targets with `>=100` closed parent
compounds, target-, accession-, document-, and parent-disjoint from every
recipient.

The deployment question is: *given only what was published about target `r` up
to `tau_r`, and abundant data on other targets, can we predict the compounds
published on `r` afterwards?* This is neither random row deletion nor a global
index date. It is a distinct estimand from the 2026-07-31 historical amendment
and replaces it as primary; the global-`tau` historical estimand is retained as
a deferred secondary track.

### Literature conclusions

**LITERATURE.** The 2026-07-31 comparison tables for MAML, FOMAML/Reptile,
ANIL, BOIL, Meta-SGD, LEO, ProtoNet, TADAM, MetaDTA, ATS, AdaMBind, CNAPs,
HyperMAML, VERSA, Taskonomy, Task2Vec, LogME, DSelect-k, negative-transfer gap,
SelectiveNet, DML cross-fitting, and conservative contextual bandits remain in
force and are not restated. Additions:

| Work | **LITERATURE:** mechanism | **INFERENCE:** transferable / non-transferable | Disposition |
| --- | --- | --- | --- |
| [ADKF-IFT](https://arxiv.org/abs/2205.02708) (ICLR 2023) | Bilevel deep-kernel GP: meta-learn the feature extractor, fit task-specific GP parameters per task, solved by the implicit function theorem; contains DKL and DKT as special cases. | Transferable: the explicit *partition* between meta-learned and task-fitted parameters, and an exact per-task solve. Non-transferable: reported strength is at **large** support sizes on FS-Mol classification; `k<=5` pKi regression is outside its evidence. | Closest architectural precedent; mandatory low-cost baseline in the amortized-adapter slot. |
| [MHNfs](https://openreview.net/pdf?id=XrMWUuEevr) (ICLR 2023) | Context-enriched molecule representations: support-set retrieval against a large context memory via modern Hopfield networks. | Transferable: enrich a query representation using the *support set* rather than a task embedding. Non-transferable: binary FS-Mol actives/inactives; no affinity scale, no provenance closure. | Retrieval-conditioning control. |
| [Conformal prediction under covariate shift](https://arxiv.org/abs/1904.06019) and audited/pseudo-calibrated variants | Reweighted exchangeability restores coverage when calibration and test distributions differ; audit models use a small labelled target sample to flag likely failures. | Transferable: the formal statement that **support and query are not exchangeable**, and that the correction is a learnable reweighting. Non-transferable: coverage guarantees require a known or estimable likelihood ratio, unavailable from `k<=5`. | Source of the core mechanism; provides the abstention-calibration control, not a guarantee. |
| [Task Singular Vectors](https://arxiv.org/abs/2412.00081) (CVPR 2025) and task-arithmetic merging | Layer-wise SVD of task matrices; low-rank task subspaces compress to ~10% while retaining ~99% accuracy, and singular-vector interaction quantifies interference. | Transferable: the diagnostic that per-task residuals are low-rank and that interference is measurable. Non-transferable: assumes a shared architecture fine-tuned per task with abundant task data; A2S recipients have `<=5` labels. | Supports demoting, not reviving, the source-atlas route. |
| [Mera, Vogt & Bajorath](https://doi.org/10.1038/s41598-025-22058-3) (Sci Rep 2025) | Meta-model learns source *instance* weights against recipient training loss; defines a negative-transfer index. | Nearest precedent for meta-learned negative-transfer control in drug design. Non-transferable: binary PKI classification with source/target chemical overlap. | Mandatory near-neighbour strong baseline. |

**LITERATURE — the exchangeability gap.** MAML, ANIL, MetaDTA, and AdaMBind
sample support and query as exchangeable draws from the same task. AdaMBind's
published protocol uses random and CD-HIT-40% task splits with support sizes 5
and 40; MetaDTA is an ICLR-2022 workshop paper using randomly drawn context.
**INFERENCE:** none of them is evaluated under a support set that precedes its
query in publication time, chemistry, and optimization stage. Under natural
scarcity that non-exchangeability is the defining feature of the task, and it is
the gap this program can occupy without claiming novelty for meta-learning,
low-rank adaptation, task conditioning, MoE, or abstention in isolation.

### Selected scientific hypothesis

**HYPOTHESIS (H-SDO).** Under document-ordered natural scarcity, the support set
`S_r` is a systematically biased sample of the recipient's query distribution
`Q_r`. That bias — in level, in scale, and in support-relative chemical position
— is a low-dimensional, target-independent, *identifiable* function of
support-only observables, and it is learnable across abundant source targets by
replaying the identical document-ordered protocol. Correcting it improves
target-macro RMSE and within-target ranking on natural scarce recipients, where
protein-conditioned reordering does not.

### Corrected meta-learning requirement: a learned adaptation operator

The primary A2S-DTA objective is not episodic packaging of a conventional DTA
model and is not learning parameters of a closed-form calibration solver. The
core object must be a trainable support-conditioned operator

```text
A_theta(S_t, x_q, p_t) -> Delta y_q
```

whose parameters are learned from abundant source-target episodes by held-out
query loss:

```text
theta* = argmin_theta E_T[ L_Q( A_theta(S_T, Q_T) ) ]
```

The model must therefore learn **how to adapt**, rather than specify an
analytic support-to-query formula. Ridge, kernel ridge, Gaussian-process or
Bayesian posterior updates, similarity interpolation, and the existing
A2S-BIR/A2S-SDO analytical forms remain permitted only as baselines, bounds,
diagnostics, or initializers. They cannot be the primary mechanism.

Valid adaptation must be query-specific: for a fixed recipient support set,
the correction must be able to differ across compounds,
`Delta(x_1,p_t) != Delta(x_2,p_t)`. Intercept-only or scale-only calibration
does not count as meta-adaptation. Candidate implementations may include a
low-capacity local residual network, neural process, conditional adapter,
hypernetwork, learned iterative optimizer, or support-conditioned modulation,
but no candidate is selected before the information and power gates pass.

The few-shot budget is binding. At `k in {1,3,5}`, adaptation must be
low-capacity, regularized, uncertainty-aware, and support-conditioned. High-
dimensional task codes, full-model fine-tuning, unconstrained hypernetworks,
large recipient-specific updates, and query-label transduction are forbidden.

Every candidate must pass all of the following before it can be called the
A2S-DTA core innovation:

1. full learned operator versus the identical model without adaptation;
2. learned operator versus random/frozen adaptation;
3. superiority to ridge residual, kernel ridge, GP/Bayesian shrinkage,
   pooled fine-tuning, and the strongest no-adaptation baseline;
4. support-label permutation and wrong-support destruction;
5. correct-target versus wrong-target degradation;
6. ligand-only versus protein-conditioned comparison;
7. stable target/component-bootstrap results at `k=1,3,5` across seeds.

Primary metrics are RMSE/MAE, Spearman, concordance, pairwise accuracy,
NDCG@10, target-macro gain, benefiting-recipient fraction,
negative-transfer rate, coverage-risk AUC, and uncertainty calibration where
applicable. The primary statistical units are targets and independent target
components, never rows, pairs, or query compounds.

The priority order is:

```text
learnable adaptation mechanism
  > identifiable target-conditioned signal
  > architecture complexity
  > model scale
```

No backbone expansion, larger Transformer/Mamba, auxiliary loss, rank,
capacity, or training duration may rescue a failed adaptation or information
gate. If the learned operator cannot beat the best closed-form baseline under
the same strict roster, the correct conclusion is that the available
few-shot information supports calibration/local interpolation only and does
not identify complex meta-adaptation.

### Demoted analytical route: A2S-SDO

The earlier A2S-SDO description is retained as an auditable hypothesis and
baseline, but it is **not** the required core meta-learning contribution. A
support-drift formula, scalar anchor, source prior, Bayesian posterior,
kernel/ridge solve, or abstention rule cannot by itself satisfy the present
definition of a learnable adaptation mechanism.

Full name: **Abundant-to-Scarce Support-Drift Operator**.

**Core mechanism (first and only core innovation): meta-learned correction of
support/query non-exchangeability.** Every prior meta-learner in this literature
treats `(S_r, Q_r)` as exchangeable draws and learns *an initialization, a task
code, or a task weight*. A2S-SDO instead makes the **shift operator from `S_r`
to `Q_r`** the meta-learned object. What is learned is not "what this target
is", but "how the compounds already published on a target relate to the compounds
published on it next".

**Second innovation (retained from A2S-MAP, unchanged in role):
cross-fitted paired-gain abstention.** A support-only quantile head predicts the
lower bound of the paired gain against an identical no-correction arm; transfer
is applied only when that bound clears a frozen material margin, otherwise the
prediction falls back exactly to the no-correction arm.

No third innovation is claimed. Low rank, comparators, cross-fitting, empirical
Bayes, and abstention are implementation primitives.

**Root problem solved.** The two failure modes of every executed A2S arm: the
pointwise pooled ridge cannot rank within a target (`Spearman ~0.085`), and the
support mean is both noisy and biased as an estimator of the query level at
`k<=5`. Neither is a protein-representation problem, so neither is solved by
encoder capacity.

**Exactly what is meta-learned.** `theta_M = {f0, c, h, q_head}`:

- `f0`: frozen target-balanced pooled source predictor (the generic ligand
  potency prior);
- `c(d, d')`: antisymmetric support-relative comparator, `c(d,d') = -c(d',d)`,
  on frozen ligand features;
- `h(u)`: the drift operator, mapping a support-only context vector `u` to a
  query-dependent correction;
- `q_head`: the cross-fitted lower-gain quantile head.

**New identifiable quantity.** The realized per-target drift
`delta_r = E[y | Q_r] - E[y | S_r]` and its query-conditional profile. This is
directly measurable on every abundant source target without any model, which is
what makes the mechanism killable before it is built.

**Why the alternatives cannot implement the same function.**

- **Fine-tuning / pooled fine-tuning:** has no representation of the support set
  as a *biased sample*; it treats support labels as unbiased supervision and has
  no explicit fallback.
- **MAML / FOMAML / Reptile:** learn one global initialization. The inner loop
  fits `S_r`; under non-exchangeability, fitting `S_r` harder moves *away* from
  `Q_r`. This yields the registered differential prediction below.
- **ANIL / Meta-SGD / LEO / BOIL:** restrict or reparameterize *which*
  parameters adapt. They still target `S_r`'s conditional mean.
- **AdaMBind / ATS:** schedule meta-training tasks by loss and support/query
  gradient similarity. Under an ordered protocol, high support/query gradient
  agreement selects precisely the *exchangeable* tasks, i.e. it schedules away
  the phenomenon A2S-SDO models.
- **MetaDTA / CNAPs / VERSA / HyperMAML / ADKF-IFT:** amortize a task predictor
  or task parameters from context. An amortized model can in principle absorb
  drift, but only if it is *trained* under the ordered protocol; none is, and
  none exposes drift as an estimand, so none can be ablated on it.
- **Ordinary MoE / DSelect-k:** route among experts. Routing does not represent
  a within-task sampling bias.
- **Stopped A2S-CFRA and demoted A2S-MAP:** both estimate a discrete or
  continuous *source-target-specific* object. A2S-SDO estimates a
  target-independent operator conditioned on observed support.

**Registered differential prediction.** MAML, ANIL, MetaDTA, and AdaMBind gains
measured under **random** support sampling will shrink or reverse under
**document-ordered** support sampling on the identical roster; A2S-SDO's gain
will not. Both protocols are run for every arm at equal budget. This is the
primary mechanism-attribution experiment and it is what distinguishes the claim
from a generic accuracy improvement.

### Mathematical formulation

Let `phi(d)` be frozen ligand features and `f0` the frozen pooled source
predictor. For recipient `r` with support `S_r = {(d_i, y_i)}_{i=1..k}` and
query compound `q`:

```text
support-anchored ranking term
  s_r(q) = (1/k) * sum_i [ y_i + c(q, d_i) ],      c(d,d') = -c(d',d)

support-only context (no protein features, no query labels)
  u(q, S_r) = [ k,
                mean_i(y_i),  sd_i(y_i),
                mean_i(y_i - f0(d_i)),
                mean_i(f0(d_i)) - f0_corpus_mean,
                min_i sim(q, d_i),  mean_i sim(q, d_i),
                support chemical diversity,
                tau_r - median document_year(S_r) ]

drift correction
  delta_hat_r(q) = h( u(q, S_r) )

prediction
  yhat_r(q) = s_r(q) + delta_hat_r(q)

no-correction fallback (identical capacity, drift term removed)
  yhat_r^0(q) = s_r(q)

paired gain on a held-out source pseudo-recipient
  G_r(k) = L_Q( yhat_r^0 ) - L_Q( yhat_r )

test action
  a_r = 1 if calibrated_lower_gain_r > delta_material else 0
  yhat_r = a_r * yhat_r + (1 - a_r) * yhat_r^0
```

**Adaptation variables and dimensionality.** At meta-test the recipient adapts
**zero free parameters**. `s_r(q)` is a closed-form aggregation over `k <= 5`
support pairs and `delta_hat_r(q)` is one forward pass of a frozen `h` on a
`<=10`-dimensional context vector. There is no inner SGD step, no per-recipient
step size, no early stopping, and no query-driven checkpoint. This is a strictly
smaller adaptation surface than A2S-MAP's intercept plus `q=4` code, and far
smaller than any MAML variant.

### Data flow and leakage contract

The four responsibility zones and the frozen `P -> M -> C -> I` direction are
retained verbatim from the 2026-07-31 preregistration, with these bindings:

- **Zone P** emits the Gate D0-R roster: per-recipient `tau_r`, label-blind
  nested `S^1 subset S^3 subset S^5`, `Q_r`, closure manifests, homology
  components, and episode hashes. `tau_r` and support membership are selected
  from document years and counts only.
- **Zone M** fits `f0`, `c`, and `h` on abundant source folds only, grouped by
  homology and provenance component, with each held-out pseudo-recipient
  excluded from every nuisance fit.
- **Zone C** computes `u`, the predicted lower gain, and the transfer action.
  It reads support labels and query *features*; it never reads query labels.
- **Zone I** emits both potential predictions, the action, and the final scores.

**Query-label leakage prevention.** Recipient query labels enter no fit, no
hyperparameter, no threshold, no seed choice, no stopping rule, and no
normalizer. They are read once, by the evaluator, after predictions are hashed
and written. `tau_r` is chosen by the frozen rule "earliest split year admitting
`k=5` support parents and `>=10` closed query parents", which reads counts and
years only. Any reverse flow is a leakage STOP.

### Meta-training and meta-test procedures

Meta-training:

1. pKi only. Partition abundant source targets into five deterministic outer
   folds grouped by homology and provenance component.
2. For each source target, enumerate document-ordered pseudo-episodes using the
   **identical** rule applied to recipients: earliest documents to support,
   strictly later documents to query, parent- and document-disjoint.
3. Truncate or resample pseudo-query sets to the frozen natural query-depth
   distribution (median `22.5`, `q25 = 16` under the candidate roster). Retain
   episode hashes and row IDs.
4. Refit `f0`, `c`, `h`, and all normalizers per fold on the other folds only.
   Fit `q_head` on a disjoint calibration fold. Rotate folds; aggregate
   target-macro.
5. Freeze everything before any natural recipient is touched.

Meta-test: expose only `S_r` and permitted covariates; compute `s_r(q)`,
`u(q,S_r)`, `delta_hat_r(q)`; apply the frozen gain policy; predict `Q_r` once;
store both potential predictions and the action.

### Ordered intermediate variables

The expected order of improvement is frozen. A later aggregate gain may not be
attributed to the mechanism if an earlier variable fails.

1. `delta_r` measured on source targets is non-zero, and its dispersion across
   targets exceeds its within-target standard error;
2. `delta_r` is predictable out-of-fold from `u` alone;
3. drift correction improves target-macro RMSE over `s_r(q)` on held-out source
   pseudo-recipients at `k=3` and `k=5`;
4. `s_r(q)` improves within-target Spearman over pointwise `f0` at `k=3/5`;
5. the gain head reduces negative-transfer rate at coverage `>=0.60`;
6. natural recipients reproduce 3-5 with the same direction.

### Baselines

Identical recipients, support/query rows, target weights, frozen inputs, and
permitted label access. Record trainable parameters, gradient evaluations,
episode count, wall time, and peak CUDA memory. **Every arm is run under both
random and document-ordered support sampling.**

1. recipient-only calibration (the primary no-transfer control) and
   recipient-only ridge/RF;
2. global `B0` / pooled source ridge (already contains source knowledge; it is a
   transfer arm, not a no-transfer control);
3. kNN-DTA and LogME/evidence selection;
4. target-balanced pooled model plus equal-parameter fine-tuning;
5. exact MAML, FOMAML, ANIL, low-dimensional Meta-SGD, optional Reptile;
6. MetaDTA/VERSA-style amortized adapter, ADKF-IFT, constrained conditional
   modulation, MHNfs-style retrieval conditioning;
7. equal-budget AdaMBind and an ATS-only scheduler control;
8. Mera-style meta-weighted source transfer;
9. demoted A2S-MAP at its registered `q in {2,4}`, `J in {4,8}`;
10. the stopped A2S-CFRA linear router and current affine support routing;
11. A2S-SDO without and with the paired-gain policy.

### Ablations and destructive controls

Core ablations (mandatory):

1. remove `delta_hat` and keep `s_r(q)` at identical capacity (kills innovation 1);
2. remove the gain policy and always transfer at identical capacity (kills
   innovation 2);
3. replace `c` by uniform anchor aggregation (isolates the comparator's ranking
   contribution from the drift contribution).

Destructive controls: support-label permutation; wrong-recipient support;
no-label recipient; support/query order reversal (query-ordered support —
must destroy the drift gain); random `tau_r`; support chemistry shuffle;
`u`-component ablation one at a time; protein-feature injection into `h` (must
**not** improve — it is the recorded unidentified term); chemistry-only nearest
neighbour; provenance-only routing; matched random abstention; high-abstention
stress. A query-oracle drift is reported as headroom only and is never a
deployable arm.

**PASS requires** support-label permutation, wrong-support, and order reversal to
remove the gain. **If protein-feature injection improves `h`,** the recorded
`E(t,d)` stop is contradicted and the run halts for a separate preregistration
rather than being reported as a win.

### Primary metrics and statistical unit

Primary paired estimand for recipient `r` and budget `k`:

```text
Delta_rk = RMSE_r(recipient-only calibration) - RMSE_r(A2S-SDO policy)
```

with the strongest frozen no-routing fallback reported as the secondary
contrast. Primary budgets `k=3` and `k=5`, fixed in advance; `k=1` is fully
reported and cannot rescue failure at `k=3/5`. Report target-macro RMSE and MAE,
within-target Spearman, pairwise accuracy, NDCG@10, gain AULC, median gain,
benefiting-recipient fraction, negative-transfer rate, coverage, risk-coverage
AUC, predicted-vs-realized gain calibration, and leave-one-component-out
influence.

**Statistical unit** is the recipient target or a frozen independent
target/provenance component. Support draws and seeds are averaged **within**
recipient first. Bootstrap 5,000 times over frozen components, paired across
arms, using the same draws for all primary contrasts. Rows, pairs, episodes,
draws, queries, folds, and seeds are never resampled as independent units.

**Power.** For a paired design at `alpha=0.05` two-sided and 80% power,
`MDE = 2.802 * sigma_delta / sqrt(n)`. At the candidate `n = 82` components,
`MDE80` is `0.062 / 0.093 / 0.124 / 0.155` for
`sigma_delta = 0.2 / 0.3 / 0.4 / 0.5`. **INFERENCE:** an 82-component roster can
detect a `~0.09` pKi target-macro gain, not a `0.05` one. The material floor is
`max(MDE80, 0.05)` and is therefore expected to bind at `~0.09`. `sigma_delta`
must be fixed from existing non-LOCK target-level contrasts before any score is
read; the retired `0.10` convenience proxy remains invalid.

### PASS / STOP criteria

Gate **D0-R** PASS requires all of: `>=50` recipients surviving full
parent/document/assay/scaffold and source closure; `>=10` closed query parents
per recipient; `k=5` nested label-blind support available for every recipient;
`>=25` independent homology/provenance components with no component exceeding
`0.15` of the roster; a computed `MDE80` under the frozen `sigma_delta`; and a
sealed `dataset-run` record. Anything less is a DATA STOP.

Gate **S0-R** (cheapest static probe) PASS requires, on source targets only:
non-zero `delta_r` dispersion exceeding within-target standard error, and
out-of-fold `R^2 > 0` for `delta_r` regressed on `u` with a positive
component-bootstrap lower bound. **STOP the entire A2S-SDO route on failure.**

Gate **K1-R** (one-seed kill, seed 1729) PASS requires: positive target-macro
RMSE gain at `k=3` **and** `k=5` over recipient-only calibration exceeding
`max(MDE80, 0.05)`; A2S-SDO above the no-drift ablation; `>=10` percentage-point
negative-transfer reduction versus always-transfer at coverage `>=0.60`; at
least one ranking metric with a positive paired 95% lower bound at `k=3` or
`k=5`; and destruction under support-label permutation and order reversal.

Full PASS additionally requires all conditions in `Excellence PASS/STOP rules`
above, plus: the registered differential prediction holds (competitor gains
shrink or reverse under document-ordered support while A2S-SDO's does not);
gains survive the middle/low chemical-similarity and scaffold-cold strata; and
the `0.329` source-seen query fraction is stratified with the gain present in
the source-unseen stratum.

STOP for: D0-R power or closure failure; `S0-R` failure; drift-ablation
equivalence; order-reversal invariance; protein-injection improvement; coverage
below `0.60`; no negative-transfer reduction; majority-recipient harm;
pseudo-only gain; query-label access; or any declared closure violation. No STOP
may be rescued by capacity, epochs, seeds, rank, unregistered losses, relaxed
splits, pKd, or post-hoc thresholds.

### Compute budget

- Gate D0-R: metadata only, CPU, no CUDA, no affinity read, via
  `main.py dataset-run`.
- Gate S0-R: `<=20` CPU-minutes, no GPU, frozen features, seed 1729.
- Gate K1-R: `<=2` GPU-hours, `<=4 GiB` peak Torch allocation, seed 1729.
- Conditional five-seed round: seeds `{1729,1731,1733,1741,1753}`, `<=10`
  aggregate GPU-hours.
- All numerical work uses `D:\anaconda\envs\drug\python.exe` and CUDA and
  records GPU name, utilization, power, wall time, peak memory, parameters, and
  gradient evaluations.

### Failure Ledger additions

| Date | Scientific hypothesis / command and seed | Data role | Observation | Failure type and alternative explanation | Decision / evidence required to reopen |
| --- | --- | --- | --- | --- | --- |
| 2026-08-01 | Audit of the Gate D0 global index-date grid `{2018,2020,2022}`; outcome-blind metadata sweep, no seed, no model fit | v4 formal pKi corpus metadata | Recipient counts under the registered band are `69/67/64/57` at `tau=2008/2010/2012/2013` versus `6` at 2018 and `2` at 2022; at `tau=2018`, `254` targets are in-band but only `6` have `>=10` post-`tau` parents. The grid was never swept below 2018. | Protocol/statistics failure in D0, not a data failure. The grid was placed inside ChEMBL-37's own curation-lag window, where post-`tau` coverage is near zero for scarce targets. | **REOPEN D0** as Gate D0-R under the per-recipient document-ordered estimand. The 2026-07-31 global-`tau` result stands as a fact about that estimand and is not relaxed. |
| 2026-08-01 | Audit of the Gate D0 dependency-component definition; no model fit | D0 dependency graph | All 193 candidates were joined into one component via shared `document_src_id`. `135,775 / 157,613` corpus rows carry `src_id=1`, the scientific-literature curation channel. | Statistics failure: a curation-channel label was used as a provenance family, forcing a single component and an undefined MDE80. | **STOP** reuse of that dependency graph. D0-R must define components from sequence-identity homology clustering plus document/patent family, and must report the largest component share. |
| 2026-08-01 | Preregistered A2S-MAP core hypothesis `g_U(d)^T z_h` as a source-specific residual reordering direction; not executed | n/a | Not run. Recorded evidence shows the target-conditioned reordering term `E(t,d)` failed in five prior framings, three with a wrong-protein or protein-free control matching the correct-protein arm. | Estimand/identifiability risk, established before execution. | **DEMOTE** A2S-MAP from primary to a conditional route. It may reopen only if the Gate S0-R query-oracle headroom arm shows a target-specific residual direction that a support-only estimator can recover. It may not be rescued by rank, `J`, or capacity. |
| 2026-08-01 | `main.py d0r-roster --out dataset/formal_training/a2s_d0r_roster.v1`; metadata only, seed 1729 | Gate D0-R roster construction | `tests/test_a2s_bir.py::test_support_draws_are_nested_and_distinct` found 15 of 82 recipients carrying five *identical* `k=5` support sets. `nested_draws` keyed draw uniqueness on the ordered tuple, so a support pool of exactly `k` produced five permutations of one set. | Implementation defect in the roster builder. It is the same defect the 2026-07-31 adversarial review flagged in the original D0 ("may repeat a k=5 support set across draws"); the preregistered contract requires five distinct label-blind draws. | **STOP** use of roster `.v1` and any result computed on it. Fixed by keying uniqueness on `frozenset`; rebuilt as `.v2` with 61 recipients and 53 components, still above both D0-R floors. All A2S-BIR results must cite `.v2`. |
| 2026-08-01 | A2S-SDO drift operator, meta-estimated on 206 source targets; no recipient label read | Source targets, document-ordered episodes | Mean later-document minus earlier-document residual is `-0.0135` pKi; the `a2s_bir_nodrift` ablation is indistinguishable from the drift arm. | Mechanism failure of the drift sub-hypothesis. Support and query differ in chemistry, not in mean level, so a scalar shift operator has nothing to correct. | **STOP** the drift term as a load-bearing component; retained only as a reported ablation. Reopen only with a query-conditional (not scalar) shift whose effect survives the order-reversal control. |

## 2026-08-01 Gate D0-R Result And A2S-BIR Model Revision

### Gate D0-R: PASS

**FACT.** `main.py d0r-roster` (`research/a2s_d0r.py`, metadata-only, no affinity
column read) produced the sealed package
`dataset/formal_training/a2s_d0r_roster.v2/` with status `PASS`. The superseded
`.v1` package is retained but must not be cited (see the Failure Ledger entry on
draw uniqueness).

| Quantity | v2 value | D0-R rule |
| --- | ---: | --- |
| Recipients | **61** | `>= 50` |
| Distinct accessions | 61 | — |
| Independent homology components (parasail SW, 40% identity) | **53** | `>= 25` |
| Largest component share | 0.049 | `<= 0.15` |
| Closed source targets | 206 | — |
| Query depth min / q25 / median | 11 / 16 / 23 | min `>= 10` |
| `tau_r` min / median / max | 1981 / 2010 / 2020 | per-recipient |
| Query source-seen fraction | 0.0062 | stratified |
| Query scaffold-cold fraction | 0.8058 | stratified |
| Query/support scaffold overlap | 0.0854 | reported |

**FACT.** All hard overlaps are zero: source/recipient `target_uid`, `accession`,
`document_uid`, `parent`, `assay_id`, and support/query `parent`. Nine recipients
across six components are homology-warm to the source pool and are stratified,
not excluded.

**FACT.** Power at `n = 53` components: `MDE80` is
`0.077 / 0.116 / 0.154 / 0.192` for `sigma_delta = 0.2 / 0.3 / 0.4 / 0.5`.
The material floor `max(MDE80, 0.05)` therefore binds near `0.12` under a
conservative planning SD. **INFERENCE:** this roster is powered for a
`~0.08-0.12` pKi target-macro gain. A gain below that band, even with a
positive bootstrap interval, is reported as under-powered and does not clear
the preregistered PASS rule.

**INFERENCE.** The 2026-07-31 DATA STOP is confirmed as a protocol artefact. The
same sealed corpus yields a roster that clears every D0-R rule once the index
date is per-recipient and dependency is defined by sequence homology rather than
`document_src_id`. Model fitting on this roster is authorized. The global-`tau`
historical estimand and the H0-S identity-collision route remain deferred.

### Measured identifiability of the interaction term

**FACT** (`research/a2s_bir.py`, source targets only, no recipient label read):
after the pooled prior `f0`, the across-target anchor SD is
`tau_b = 0.282` pKi, the within-target residual noise SD is `sigma = 0.997`
pKi, and the across-target scale of a meta-learned global residual code is
`tau_z ~ 0.18-0.19` pKi.

**INFERENCE — the budget calculation.** Resolving one global interaction code
coordinate to its own across-target scale requires roughly
`(sigma / tau_z)^2 ~ 27-31` support labels. The available budget is `k <= 5`.
This is the first quantitative statement in the program of *why* `E(t,d)` has
never been identified: it is not an architecture deficit, it is a sample-size
deficit of roughly six-fold. The identifiability certificate accordingly fires
on `0.000` of episodes at every `k`, and every global-code arm collapses exactly
onto the anchor arm.

**FACT.** The meta-learned support-to-query drift on source targets is
`-0.0135` pKi. **INFERENCE:** the A2S-SDO drift sub-hypothesis is falsified as a
first-order effect; document-ordered support and query differ in *chemistry*,
not in mean level. The `a2s_bir_nodrift` ablation is retained but the drift term
is not load-bearing and no claim may rest on it.

**FACT.** The raw support anchor is actively harmful: `f0_anchor` scores RMSE
`1.6159 / 1.4680 / 1.4394` at `k=1/3/5` against `f0_only` at `1.3724`.
Hierarchical shrinkage repairs it: `f0_anchor_shrunk` scores
`1.3511 / 1.3243 / 1.3132`, beating recipient-only calibration
(`1.5345 / 1.3942 / 1.3523`) at every budget. **INFERENCE:** because
`tau_b << sigma`, the k-shot anchor must be a shrinkage estimate; an unshrunk
support mean is the dominant error source at small `k` and is why several
earlier anchor-style arms looked weak.

### Demoted closed-form baseline: A2S-BIR

**Full name: Budget-constrained Identifiable meta-Residual.** This supersedes
A2S-SDO's drift operator as the primary mechanism and adopts the user's
"interaction-identifiable, budget-constrained meta-residual" framing. Two arms
of the interaction term were implemented and are reported against each other:

```text
yhat_r(q) = f0(q) + b_hat_r + certify(q) * local_r(q)

b_hat_r   : hierarchical (shrunk) anchor from the joint [b, z] posterior
local_r(q): kernel-ridge correction on the k OBSERVED support residuals,
            local_r(q) = k_qs^T (K_ss + lambda I)^-1 (r_S - b_hat_r)
certify(q): per-query coverage gate sigmoid(a * (max_i sim(q, d_i) - c))
```

`K` is a cosine kernel under a **meta-learned diagonal metric initialised at the
raw cosine**, so the learned metric can only improve on the untrained one.
**The meta-learned `lambda` is the identifiability budget**: the realised
effective degrees of freedom `tr(K_ss (K_ss + lambda I)^-1)` are
`0.58 / 1.47 / 2.14` at `k = 1/3/5`, i.e. the model spends roughly 40% of its
nominal label budget. Meta-learned parameters are the diagonal metric, `lambda`,
and the gate `(a, c)`; nothing else adapts at meta-test and no protein feature
enters.

**INFERENCE — why local and not global.** A rank-`m` global code asks the `k`
labels to locate a target in a shared basis, which the budget calculation shows
is impossible at `k<=5`. A support-local correction spends at most `k` degrees
of freedom because it is a convex/kernel combination of the `k` *observed*
residuals, so it is identifiable by construction; what must be meta-learned is
which support compound is informative for a query, and how much local support a
query has. This is a function-class correction discovered by the measured
failure of the global arm, not a capacity increase.

**Prohibited.** The global-code arms (`a2s_bir_global`, `a2s_bir_svd`,
`a2s_bir_random_basis`) are retained only as the measured negative control for
the budget calculation. They may not be revived by increasing `m`, rank, epochs,
or `J`; the reopening condition is a support budget above `~30` labels or a
demonstrated reduction of `sigma`.

## 2026-08-01 A2S-MDK One-Seed Result (Gate K1-R)

### Estimand change

**Primary estimand is now target-specific within-target ranking**, not accuracy:

```text
TS(metric) = metric(operator | recipient's own support)
           - metric(same operator | a different recipient's support)
```

paired per recipient, 5 draws averaged within recipient, bootstrapped 5,000
times over the 53 independent homology components. RMSE is secondary. The
rationale is that a generic ligand-potency improvement is not transfer; only the
part destroyed by substituting another target's support is target-specific.

### Executed runs

`main.py a2s-scao --protocol {ordered,random}`, seed 1729, roster
`a2s_d0r_roster.v2`, RTX 4060, 1,187,523 parameters, 3,000 steps, 4,155 ordered
source episodes, wall 353 s, peak Torch 2,533 MiB. Artifacts:
`reports/active/a2s_mdk_{ordered,random}_seed1729.json`.

### Result: target-specific CI gain, ordered protocol

| arm | k=1 | k=3 | k=5 |
| --- | --- | --- | --- |
| **A2S-MDK** | **+0.0137 [+0.0054, +0.0232]** | **+0.0190 [+0.0070, +0.0321]** | **+0.0249 [+0.0085, +0.0421]** |
| A2S-MDK without contrast | +0.0077 [+0.0004, +0.0156] | +0.0085 [-0.0055, +0.0233] | +0.0091 [-0.0086, +0.0283] |
| pooled ridge fine-tune | +0.0076 [+0.0014, +0.0142] | +0.0161 [+0.0051, +0.0282] | +0.0247 [+0.0092, +0.0420] |
| A2S-SCAO (FiLM/hypernetwork) | -0.0165 [-0.0351, +0.0016] | -0.0186 [-0.0387, +0.0001] | -0.0226 [-0.0453, -0.0021] |

**FACT.** A2S-MDK is positive with a bootstrap lower bound above zero at all
three budgets and monotone in `k`. Spearman and NDCG@10 agree
(`+0.0700 [+0.0269, +0.1180]` and `+0.0265 [+0.0092, +0.0457]` at k=5).
Support-label permutation also degrades it
(Spearman `+0.0359 [+0.0006, +0.0737]`, NDCG `+0.0200 [+0.0040, +0.0385]`),
so the effect depends on the support labels and not only on support chemistry.

**FACT — the one load-bearing innovation.** Removing the counterfactual support
contrast from the meta-objective drops the gain by 45-64% and makes it
non-significant at k=3 and k=5. This is the only component of A2S-MDK that earns
its place.

**FACT — two components are not load-bearing.** The learned budget gate
saturates open (mean 0.970/0.998/1.000; `a2s_mdk_nobudget` is numerically
identical to the full model), and the deep-kernel branch settles at
`alpha = 0.194` while `a2s_mdk_nodeep` matches or beats the full model on
absolute CI (0.5593 vs 0.5584 at k=5). Learned `lambda = 0.0116` on the
trace-normalised kernel, giving effective degrees of freedom 0.99/2.89/4.74 at
k=1/3/5 -- the operator spends essentially its whole label budget.

**FACT — the FiLM/hypernetwork operator is anti-specific.** Under
document-ordered meta-training it is significantly negative at k=5: a different
target's support ranks *better* than the correct one (CI 0.5827 vs 0.5569).

**FACT — the registered differential prediction is NOT confirmed.** A2S-MDK's
target-specific CI gain at k=5 is `+0.0249` under ordered meta-training and
`+0.0317` under random meta-training; pooled fine-tune is `+0.0247` under both.
Only A2S-SCAO shows the predicted collapse (`-0.0226` ordered vs `+0.0139`
random). The prediction therefore holds for amortised parameter-generating
operators, not for the primary model, and may not be cited as support for it.

**FACT — the accuracy criterion fails.** Target-macro RMSE is 1.3028 for
A2S-MDK against 1.3492 for the shrunk-anchor fallback, a gain of 0.046, below
the preregistered material floor `max(MDE80, 0.05)` which binds near 0.12 on
this roster. MAML (1.3465 / CI 0.5534) and ANIL (1.3537 / CI 0.5450) are beaten
on both axes after their divergence was repaired.

**INFERENCE.** A transferable target-specific ranking mechanism exists and is
measurable at k={1,3,5}. A2S-MDK matches the no-meta-learning ridge head at k=5
and roughly doubles it at k=1, where the contrast objective matters most. This
clears the registered ranking criterion and fails the registered accuracy
criterion. Meta-learning is not yet demonstrated to add value beyond a
support-conditioned ridge solve at k=5, so `META_LEARNING_ADDS_VALUE` is not
claimed on one seed.

### Failure Ledger additions

| Date | Scientific hypothesis / command and seed | Data role | Observation | Failure type and alternative explanation | Decision / evidence required to reopen |
| --- | --- | --- | --- | --- | --- |
| 2026-08-01 | A2S-SCAO FiLM/hypernetwork support-conditioned operator; `main.py a2s-scao --protocol ordered`; seed 1729 | 61 recipients / 53 components | Target-specific CI gain `-0.0226 [-0.0453, -0.0021]` at k=5; a different target's support ranks better than the correct one. | Mechanism failure. Query loss is dominated by the target-independent ligand term, so an operator that generates parameters from a set embedding learns to use the support as generic context. | **STOP** parameter-generating amortised operators for this estimand. Reopen only with an objective that makes target-specificity explicit and an architecture whose inner adaptation is the identified estimator. |
| 2026-08-01 | A2S-MDK budget gate and deep-kernel branch; same run | Same | Budget gate saturates open (0.970/0.998/1.000) and is numerically inert; deep branch `alpha=0.194` with `nodeep` matching or beating the full model on absolute CI. | Component-level null. Neither is falsified as an idea, but neither is load-bearing at this roster size and support budget. | **STOP** claiming either as an innovation. Reopen only if a larger roster or a support budget above ~30 labels makes the identifiability gate bind. |
| 2026-08-01 | Registered ordered-vs-random differential prediction | Same | A2S-MDK k=5 TS-CI is `+0.0249` ordered vs `+0.0317` random; pooled fine-tune identical under both. | Prediction not confirmed for the primary model. Non-exchangeability affects amortised operators but not kernel-ridge adaptation. | **STOP** citing the differential prediction as evidence for A2S-MDK. It survives only as a statement about A2S-SCAO. |
| 2026-08-01 | MAML/ANIL baselines, first implementation | Source episodes | MAML returned `nan`, ANIL returned RMSE 384530: the inner loop diverged when regressing raw pKi from 2058-d input. | Baseline implementation defect, found before any comparison was reported. | **CONTINUE** with the repaired baselines (residual target, LayerNorm, zero-init head, inner-gradient clipping, inner LR selected on source episodes only). No comparison against the diverged version may be cited. |

### Unique next action

**Build the Gate D0-R roster through `main.py dataset-run`**, outcome-blind,
from the sealed v4 formal corpus, emitting a new versioned package under
`dataset/formal_training/` without overwriting v4. Required contents: per-recipient
`tau_r`; label-blind nested `S^1 subset S^3 subset S^5`; closed `Q_r`;
Murcko scaffolds (currently absent — only `connectivity_inchikey` exists);
sequence-identity homology components; document/patent-family closure removing
the `169` shared documents and `756` shared parent compounds; the source-seen
query stratum; the frozen `sigma_delta`; and `MDE80`. Do not run Gate S0-R,
fit `c` or `h`, run any baseline, or read a query label until that package is
sealed and D0-R is adjudicated PASS or STOP.

## 2026-08-01 A2S-CMAL Core Objective And Execution Contract

This section supersedes the preceding **Unique next action**. D0-R has since
been sealed and passed; the next scientific action is external five-seed
training of the model specified here, not another local formal run.

### Non-negotiable core

提出并检验一个 **protein-conditioned learnable meta-adaptation operator**：从
abundant source targets 学会如何使用极少 recipient support measurements，
在严格 target-ID unseen 条件下，为每个 query compound 产生 target-specific、
query-dependent affinity ranking correction。最终解必须是从 source episodes
学得的非闭式适配机制；标量校准、插值、相似性检索、ridge/kernel solve 或
固定 posterior update 只能是 baseline。

必须同时保持一个完全不依赖 support/query chemistry leakage 的 support-free
DTA base。base 先训练、后冻结；适配前后的 `use_support=False` 输出必须逐
元素一致。任何只改善绝对 RMSE、却不能在替换 wrong-target support 后失效的
方法，都不能被称为 target-specific transfer。

### Frozen final architecture: A2S-CMAL

`research/a2s_cmal.py` implements:

```text
base:     f_theta(protein, query_ligand)
state:    A_phi(protein, {(support_ligand, measured residual)})
adapted:  f_theta(protein, query_ligand)
          + Delta_phi(protein, state, query_ligand)
```

The support state uses learned self-attention over protein-ligand pair and
measured-residual tokens. Each query cross-attends to that state before a
nonlinear delta head. There is no analytic inner solve, budget gate, deep
kernel, or posterior update. The prior SCAO stop is not silently ignored: its
reopening condition was an explicit target-specificity objective and an
identifiable support intervention. CMAL supplies that intervention and must
still pass the wrong-support tests below; architecture alone is not evidence.

### Load-bearing objective

For every `(recipient protein, correct support, query)` training episode, keep
protein and query fixed and substitute three different-target supports:

1. random different target;
2. ESM-2 protein-hard different target;
3. support-only scaffold+ECFP4 chemistry-matched different target;
4. same-compound wrong-target residual/label swap (constructed online).

Score the correct arm and four counterfactual arms by **frozen-base-anchored
post-adaptation query-ranking gain**, and use InfoNCE requiring the correct
support to win. A wrong arm receives no further reward after becoming worse
than the frozen base. Negative
mining reads neither query label nor query chemistry. The three arms are fused
with the online label-swap arm in one GPU forward. `no-counterfactual` is mandatory; embedding similarity
contrast is a negative control, not the claimed mechanism.

### Sealed model-ready data

- corpus: `dataset/formal_training/chembl37_pki_formal.v4`, 157,613 exact pKi
  measurements, 1,098 targets, 82,646 parent compounds;
- roster: `dataset/formal_training/a2s_d0r_roster.v3`, PASS, 206 source and 63
  recipient targets, 55 independent homology components, 5 draws, k=1/3/5;
- episodes: `dataset/formal_training/a2s_cmal_episodes.v3`,
  `READY_FOR_EXTERNAL_FORMAL_TRAINING`, content SHA-256
  `2df5831bc8a51df93dc54531302327716fcca8900ec43f1aa37f16ed2fb9485a`;
- 30,123 label-blind episodes: 23,127 train, 2,595 validation, 3,456 source
  test, 945 recipient test;
- each episode freezes aligned parent and `measurement_uid` lists from the
  admitted context-main label table; trainer label joins are measurement-exact;
- support/query parent, document, measurement, ordered-time, nested-support, source split
  component, and all three package-level same-target-counterfactual violations: zero.

Strict wording is **unseen target ID**. Nine of 63 recipients are homology-warm;
the 54-target homology-cold subgroup must be reported separately.

### Required formal experiments

The full, source-cited plan is in `experiment.md`. Minimum core matrix:

- five frozen seeds `{1729,1731,1733,1741,1753}` and k=`{1,3,5}`;
- support-free base, calibration, kNN/retrieval, no-meta ridge, A2S-MDK fixed
  posterior, MAML, ANIL, and A2S-CMAL at equal data/representation budget;
- correct support versus random, protein-hard, chemistry-matched wrong support,
  plus label permutation;
- no-counterfactual, each negative type alone, zero/shuffled protein, no
  support attention, embedding-contrast negative control;
- target-macro CI, Spearman, NDCG@10 primary; RMSE, MSE, MAE, Pearson, R2,
  pairwise accuracy and AUPR@7 secondary; component-bootstrap uncertainty;
- all 63 recipient targets and the 54-target homology-cold subgroup separately.

AdaMBind paper reproduction remains a separate PAPER-EXACT track: BindingDB,
KIBA and Davis; random 8:1:1 and CD-HIT-40% cluster 8:1:1 task splits; support
5/40 plus sensitivity 5/10/20/30/40; five replications; MSE, CI, R2, Spearman,
Pearson; published baselines/ablations; KIBA-to-BindingDB transfer. It may not be
presented as the stricter D0-R experiment.

### Formal claim gate

`META_ADAPTATION_MECHANISM_SUPPORTED` requires all of the following across the
five seeds: positive component-bootstrap ranking gain over the support-free
base; advantage over no-meta ridge/fixed posterior and equal-budget meta
baselines; correct support better than every wrong-support arm and permuted
labels; protein shuffle/zero destroys the advantage; no-counterfactual reduces
it; frozen support-free predictions are unchanged. A k=1-only separation must
be written as k=1-only. No new arbitrary effect threshold is introduced.

The one-seed MDK evidence above is now explicitly a pilot. It does not establish
that learned meta-adaptation adds value, and the sentence claiming a mechanism
“exists” must not be cited without the five-seed CMAL tests.

### Compute / GPU contract

Formal recipient training is prohibited on this device. `main.py a2s-cmal
--formal` exits before label access unless `A2S_FORMAL_EXTERNAL=1` is set on the
designated external host. Local `--smoke` reads source labels only and evaluates
source meta-validation.

The former sawtooth defect was per-episode pandas lookup, host transfer and
frequent scalar synchronization. CMAL materializes features, labels, episode
indices and all counterfactual supports on GPU before training, batches episode
sampling, fuses the four support arms, and logs only at five phase intervals.
The measurement-identified v2 250+250-step batch-64 local profile reached 9,342
base and 5,659 adapter episodes/s; 54.5% of telemetry samples were at least 40%
utilization and P90 was 56.0%. Mean 32.8% includes phase boundaries. This verifies the pipeline repair,
not full hardware saturation or a scientific result.

### Current status and next action (updated after source mechanism diagnosis)

- data preprocessing/model-ready package: **COMPLETE**;
- executable implementation and mechanical unit tests: **COMPLETE, EXPERIMENTAL**;
- local source-only CUDA mechanism smoke/profile: **COMPLETE**;
- source validation/holdout mechanism gate: **FAILED**;
- formal five-seed recipient training: **BLOCKED AND NOT RUN**;
- scientific claim: **NOT YET AUTHORIZED**.

The latest component-disjoint diagnostic learned positive correct-vs-wrong
support specificity on one source holdout, but NDCG@10 remained negative and
the validation split deteriorated. This is not a transferable-mechanism result.
The actual current implementation remains `research/a2s_cmal.py`; the legacy
`model/` folder is not the CMAL execution path and must not be promoted or
published. No GitHub commit/push is authorized until a robust source-only
validation plus holdout breakthrough occurs.

Next action: use `reports/active/CMAL_EXTERNAL_AGENT_PROMPT.md` with the three
attachments specified in `reports/active/CMAL_FAILURE_HANDOFF.md`. Diagnose
identifiability and the operator's ranking inductive bias before any further
training. Recipient labels remain sealed.

## 2026-08-01 Superseding A2S Decision After The Balanced V2 Gate

The external review and balanced source-only rerun supersede the preceding next
action. CMAL is frozen. The v2 gate repaired the 97% OOF fold imbalance but did
not obtain a positive lower bound for `Delta_label` at k=1/3/5 or
`Delta_assign` at k=3/5. The decision remains
`NO_GO_INFORMATION_NOT_ADMITTED`.

Current ordered tasks:

1. Complete a label-free same-assay/MMP coverage and component-power census on
   the v2 retained rows.
2. If and only if coverage and MDE are sufficient, preregister one local k=3/5
   G0/G1 assignment-information gate.
3. If the local gate fails with valid synthetic controls, abandon the current
   passive A2S episode construction.
4. If it passes, compare a finite depth-1 SAR grammar with Matsy/MMP/MMS,
   mixed-effects, categorical empirical Bayes, graph-difference, and
   learned-kernel baselines before implementing a meta-operator.
5. Keep locked-source and recipient labels sealed until their separate one-time
   gates are frozen and passed.

There is no positive mechanism breakthrough. No current code may be promoted to
`model/` or uploaded to GitHub.

## 2026-08-01 Final PI Meta-Mechanism Redesign

The complete no-code redesign is frozen in
`reports/active/A2S_FINAL_PI_META_MECHANISM_REDESIGN_2026-08-01.md`.

Five candidate hypotheses were evaluated. The selected highest-potential
scientific direction is `Active Diagnostic Response-Mode Operator`, which
jointly meta-learns diagnostic measurement acquisition and a finite abstaining
query intervention. This is not an admitted model result and changes the task
from passive few-shot adaptation to active k-measurement adaptation.

The next admissible work is Phase 0 only: candidate-pool contract, active-oracle
headroom, label/assignment information admission, power/MDE, and leakage audit.
No implementation may begin until those source-only gates pass.

## 2026-08-02 Passive A2S Meta-Adaptation Supersession

The user has restored the passive A2S-DTA objective and superseded the active
measurement direction. The active objective is again:

> Learn a transferable few-shot meta-adaptation mechanism for unseen-target
> DTA ranking from k={1,3,5} recipient measurements.

TRACE, RIP, KRR, retrieval, uncertainty, and active acquisition are auxiliary
baselines or diagnostics only. Source `locked` and recipient labels remain
sealed.

The PIRS protein-segment interaction-state gate completed on source fit-to-probe
transfer and returned `INTERACTION_STATE_REPRESENTATION_NOT_ADMITTED`. Its
synthetic control passed, but the full-support representation oracle, k-shot
gain, support assignment, protein destructions, and matched representation
controls all failed. R1 is prohibited and nothing is promoted to `model/` or
`script/`.

The active branch is
`research/a2s-conformational-free-energy-state-20260802`. Its biological basis
is conformational selection and population shift: sparse support measurements
may update one or two population logits over physically observed protein states,
which then change state-specific query energy surfaces. The complete formulation
and stop rules are in
`reports/active/A2S_CONFORMATIONAL_FREE_ENERGY_STATE_PROPOSAL_2026-08-02.md`.

The only authorized next action is Gate C0: a label-blind structural coverage
and external apo/holo semantic gate. Do not reuse source probe labels, open
source locked/recipient labels, train an affinity model, implement the learned
support operator, or promote code unless C0 and the subsequent fit-only
representation gate pass.
