# MetaSieve-DTA Experiment History

Last updated: 2026-08-10.

> HISTORICAL EVIDENCE ONLY. This file is not an execution plan.
>
> Current execution authority is, in order: `AGENT_HANDOFF.md`, `task.md`, and
> the active section of `experiment.md`. The canonical completed-experiment
> summary is `report/EXPERIMENTAL_EVIDENCE_LEDGER.md`. Any historical
> instruction that conflicts with those files is superseded.

## Current Project Decision

### F-99 — S7/L2B B4 evidence-integrity and failure-localization audit

The B4 development result remains **FAIL CLOSED** under its preregistered
practical-effect thresholds.  The valid interpretation is narrower than the
original report: B4 is directionally above prevalence, ligand-only,
wrong-protein, motif-shuffle and wrong-ligand controls, with every component
bootstrap lower bound above zero, but only G1 reaches the frozen `0.02` margin.
This supports
`PROTEIN_SIGNAL_IDENTIFIED_BELOW_PREREGISTERED_EFFECT_SIZE`; it does not admit
an exact-residue mechanism to production and does not test affinity direction.

The earlier statement `Representation: THE CAUSE` is **withdrawn as a unique
causal localization**.  The synthetic control excludes gross pipeline
untrainability for a known input-derived target, but it does not distinguish
weak/non-exhaustive MONN interaction labels, atom correspondence error,
pair-coupling limitations, real-task objective/optimization behavior, missing
geometry, or inadequate residue representation.  Frozen ESM2 B5 is therefore
a justified minimal discriminator, not a proven repair.

Four integrity conditions must be closed before B5 is operationally allowed:

1. verify MONN/PDB-to-RDKit atom correspondence by element, degree,
   connectivity and source atom identity rather than heavy-atom count alone;
2. replace index-ordered tie breaking with score-threshold-group tie-aware AP
   and persist per-pair labels/scores so metrics can be recomputed;
3. materialize negative-sampling and control maps plus parameter, activation,
   gradient and checkpoint diagnostics;
4. construct and freeze publication/time closure before opening any
   confirmation cohort.

Held-out B is a scaffold-strict **subset** of held-out A and is retained as a
development robustness analysis, not independent replication. Wrong-protein is
a nuisance corruption control until exact pairwise homology, fold familiarity,
reuse and truncation semantics are fully reported. No threshold is relaxed and
the already inspected development panel is not made untouched again.

The conditional roadmap is frozen at a high level: integrity audit -> one
matched frozen-ESM2 B5 run -> residue/atom/coupling marginal triage if B5 does
not pass -> only then consider label-aware PU/continuous structural teachers or
partner-conditioned geometry. Larger PLMs, generic SSL, typed energy,
affinity, DAVIS, support adaptation and production `z` remain unauthorized.
Any later SSL objective must train protein-ligand coupling rather than only a
single-side embedding. The frozen mathematical operator remains unchanged;
localization is only a candidate biological measurement frontend.

```text
P0_CANONICAL_DTA_V2_PASS
P1A_OPEN_STRUCTURE_PILOT_PASS
P1B_GEOMETRY_GATE_PASS
P1R2B_D0_C_PASS
P1R2B_D1_PASS
P1R2B_E0R2_SYNTHETIC_OBJECTIVE_DESIGN_SOLVER_IDENTIFIED
P1R2B_T_BASIS_R0_RADIAL_BASIS_PARTNER_RECOVERABILITY_IDENTIFIED
S4_AGGREGATE_ESM_ECFP_NOT_PROTEIN_SPECIFIC_SCOPE_CORRECTED
P1R2B_S5_LOCAL_MECHANISM_OBSERVABILITY_REGISTERED_NOT_RUN
EXTERNAL_S5_S9_CLAIMS_NOT_REPRODUCED
SSL_G0_DATA_FEASIBILITY_CONDITIONAL_METADATA_ONLY
REAL_SOURCE_AFFINITY_FROZEN
DAVIS_TRANSFER_FROZEN
P1R2B_T_FROZEN
P2_TO_P4_FROZEN
RECIPIENT_LABEL_READS=0
DAVIS_LABEL_READS=0
```

```text
PAIR_LOCAL_P1B_OBSERVABILITY_NOT_TESTED
POSE_FREE_CLASS_NOT_CLOSED
AFFINITY_ENERGETICS_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

The detailed earlier P1C--X0/L0 verdicts remain below as historical evidence;
they were removed from this current-status block to avoid presenting closed or
superseded branches as active work.

`E-AFF-L0` addresses Claim A (protein affinity location). `E-AFF-X0-B` and `X1`
address Claim B (non-additive protein-by-ligand interaction). These are
different estimands; neither stage replaces or authorizes the other.

The BioLiP2 plus RCSB PDB/CCD route passed P1A and the structural bridge passed
P1B. P1R0 confirmed a legacy readout coordinate defect. P1R1 repaired
permutation semantics and recovered partial correct-protein signal, but failed
the unchanged transfer thresholds. P1R2A isolated a noncollapsed non-additive
interaction component, which remained correct-protein-sensitive but did not
improve the ligand baseline. P1R2B0 found no rescue after PCA32/global pooling;
it did not test pre-pooling local nonlinearity. P1R2B1 amplified compatibility
but not correct-protein affinity direction. F0/F0R did not establish a live-API
corpus. D0-C/D1 subsequently built and governed a release-pinned static corpus.
E0 was authorized and failed its label-free synthetic pre-gate; the source
gate was not run, and typed-interaction T and P2-P4 remain frozen.

The PLINDER MLSB and selected-shard routes below remain failure evidence only.
They are not the mandatory execution path and must not be reopened without new
authorization.

This ledger records failed, rejected, superseded and blocked branches. A
blocked experiment is not a failed model attempt.

## Valid Phase 3 Protocol

- Datasets: DAVIS primary, KIBA replication.
- Seeds: 11, 23, 37.
- Split: whole 40% Smith-Waterman homology clusters assigned to source,
  metaval, or recipient.
- Strict split sizes: DAVIS 220/63/96 and KIBA 139/28/62
  source/metaval/recipient targets; zero clusters crossed a split.
- Evaluation: one fixed `k=5` support set per recipient target, disjoint from
  up to 64 queries; every arm used the same support and queries.
- Representation training: 120 Stage 1 steps and 120 Stage 2 steps per
  architecture. Stage 1 trained support/adaptation plus an auxiliary CSMO with
  protein, ligand, and interaction frozen. Stage 2 also trained interaction.
- Frozen sanity stage: representation frozen and hashed, auxiliary CSMO
  discarded, then a fresh CSMO trained for 600 steps on 1,536 cached source
  episodes.
- ESM outputs, protein projection, and ligand GINE remained frozen. No ESM
  adapter and no full ESM fine-tuning were run.
- Metrics: target-macro CI, Spearman, NDCG@10, and normalized-scale RMSE.
- Intervals: 95% Student-t intervals over the three seed-level target-macro
  values.
- Scalar prediction: midband mean, an engineering readout outside the frozen
  operator.

## Failure Ledger

### F-01: Statistical-only target adaptation collapsed

**Attempt.** The pre-biological model represented a target through a 6-D Morgan
PCA ligand representation and a small set of support-label statistics. Protein
sequence never entered the model.

**Result.** The `k=5` support mean estimated target level with RMSE 0.363 on
DAVIS and 0.389 on KIBA, comparable to or larger than the between-target mean
spread. On metaval, the query-only CSMO view reached objective 0.1963 and beat
all support-dependent views, which were at least 0.2017.

**Reason.** The only target channel was five or fewer noisy affinity labels.
There was no sequence description and no target-ligand interaction channel, so
ignoring support was a rational optimization outcome.

**Action.** Phase 2 added cached ESM representations, ligand graphs, an explicit
interaction module, a Set Transformer, and Q-PMA while leaving CSMO and the band
operator unchanged.

### F-02: `sqrt(0)` poisoned Stage 1 gradients

**Attempt.** Train the new biological statistic through label-spread
coordinates at support sizes `k in {1,3,5}`.

**Result.** The first forward objective was finite, about 0.469, but the first
gradient norm was NaN and the next model state was unusable.

**Reason.** At `k=1`, support variance is exactly zero and the derivative of
`sqrt(var)` is singular.

**Action.** Add a `1e-12` variance floor inside both global and local standard
deviations. This is retained as a load-bearing numerical guard.

### F-03: The first query-conditioning audit used a degenerate support size

**Attempt.** Check that the same support and two queries produce different Q-PMA
adaptations.

**Result.** The audit reported zero change in attention and adaptation.

**Reason.** The selected episode had `k=1`. With one support token, softmax
attention is exactly one for every query, so query-conditioned pooling cannot
change.

**Action.** Require `k>=2` for this structural audit and report `k=1` as a
mathematical degeneracy rather than a model failure.

### F-04: Materialized episode graphs exceeded the local memory budget

**Attempt.** Copy every support and query graph into each episodic batch.

**Result.** A representative KIBA batch required about 5 GB for bond tensors
alone and repeated the same graph encoding many times.

**Reason.** Storage scaled with `batch * (k+1) * atoms^2 * bond_features` even
though drugs and targets were repeated.

**Action.** Batches now carry indices into immutable shared drug and protein
banks. Phase 3 additionally caches outputs of the frozen protein and ligand
towers. Cached and uncached `z` differed by at most `6.7e-16` in the equivalence
check.

### F-05: The default target-ID split was not a strict unseen-target split

**Attempt.** Use a label-independent salted target-ID split.

**Result.** 56.3% of DAVIS recipient targets and 53.4% of KIBA recipient targets
had at least 40% local sequence identity to a training target.

**Reason.** Distinct target identifiers do not imply biological novelty.

**Action.** Phase 3 uses single-linkage 40% Smith-Waterman homology clusters as
the indivisible split unit. No cluster crosses any split boundary.

### F-06: The first full Phase 3 run had unmatched initializations

**Attempt.** Train mean, Set Transformer, and Q-PMA arms with architecture-based
random seeds, then train each frozen-representation CSMO with a different seed.

**Result.** The run completed, but pooling arms differed in their frozen protein
and ligand tower weights, and representation arms differed in downstream head
initialization. Its performance values were discarded and its seed-11 result
file was overwritten.

**Reason.** The comparison changed initialization as well as the mechanism, so
performance differences could not be attributed to biological representation
or pooling.

**Action.** All common frontend weights, training episode draws, and sieve-head
starts are now matched within each seed. Cache attachment verifies a SHA-256 of
the protein and ligand tower parameters.

### F-07: The first no-interaction arm was an unfair inference-only bypass

**Attempt.** Remove interaction only at recipient inference.

**Result.** The smoke arm ran, but it exposed the trained head to a representation
distribution it had never seen.

**Reason.** This measured distribution shift as well as the value of the
interaction module.

**Action.** The valid no-interaction arm trains and freezes its own representation
and fresh CSMO. Only true corruption controls, such as protein/ligand shuffle and
support-label permutation, remain inference-time interventions.

### F-08: The biological representation failed frozen sanity

**Attempt.** Freeze each representation and train only a fresh, matched CSMO.

**Result.** Seed means and 95% intervals were:

| dataset | arm | CI | Spearman | NDCG@10 | RMSE |
|---|---|---:|---:|---:|---:|
| DAVIS | old statistical | 0.6267 [0.5953, 0.6581] | 0.2201 [0.1514, 0.2887] | 0.3545 [0.2353, 0.4737] | 0.2120 [0.2114, 0.2126] |
| DAVIS | biological | 0.5045 [0.3131, 0.6959] | 0.0020 [-0.3311, 0.3351] | 0.1933 [-0.0269, 0.4134] | 0.2119 [0.2118, 0.2121] |
| KIBA | old statistical | 0.5285 [0.5159, 0.5411] | 0.0715 [0.0379, 0.1051] | 0.5477 [0.5206, 0.5748] | 0.1909 [0.1848, 0.1971] |
| KIBA | biological | 0.5158 [0.4536, 0.5780] | 0.0358 [-0.1000, 0.1715] | 0.4777 [0.4386, 0.5168] | 0.1919 [0.1888, 0.1950] |

Biological-minus-old paired contrasts were unfavorable for all ranking metrics
on both datasets. KIBA NDCG was decisively worse by -0.0701
[-0.0915, -0.0486].

**Confirmed proximal reason.** The learned biological state carried much less
query variation. Mean full-`z` pair distance was 0.0231 versus 0.4598 for the
old statistic on DAVIS, and 0.0282 versus 0.4102 on KIBA. The new representation
therefore reached the matched CSMO as an almost constant state.

### F-09: Support labels did not produce useful adaptation

**Attempt.** Permute the five support labels while preserving support compounds,
the label multiset, target, and queries.

**Result.** Positive values would favor real support:

| dataset | CI | Spearman | NDCG@10 | RMSE improvement |
|---|---:|---:|---:|---:|
| DAVIS | -0.0022 | -0.0047 | -0.0011 | approximately 0 |
| KIBA | 0.0004 | 0.0010 | 0.0007 | approximately 0 |

Every interval covered zero. DAVIS point estimates slightly favored permuted
support on all four metrics.

**Confirmed proximal reason.** The model did not use the association between a
support compound and its label. Global mean/spread coordinates survive label
permutation, while the Q-PMA association channel had negligible influence.

### F-10: Query conditioning was structurally present but operationally inert

**Attempt.** Compare Q-PMA to a constant-seed attention intervention under the
same fixed support.

**Result.** Mean pairwise distance in the adaptation coordinates `z[17:26]` was
only 0.00078 on DAVIS and 0.00258 on KIBA. It was exactly zero for mean pooling
and query-independent Set Transformer arms, confirming that the test isolates
query conditioning. Q-PMA performance contrasts were:

| dataset | CI | Spearman | NDCG@10 | RMSE improvement |
|---|---:|---:|---:|---:|
| DAVIS | 0.0004 | -0.0004 | -0.0038 | approximately 0 |
| KIBA | 0.0002 | 0.0003 | -0.0018 | -0.0005 |

**Reason.** Q-PMA changed `z` in the correct structural direction but by too
little to improve the frozen downstream prediction. Mean, Set Transformer, and
Q-PMA intervals substantially overlapped on every metric.

### F-11: Protein information was not used for unseen-target performance

**Attempt.** Replace each target's sequence representation by a deterministic
different target while leaving labels, ligands, and queries unchanged.

**Result.** Positive values would favor the correct protein:

| dataset | CI | Spearman | NDCG@10 | RMSE improvement |
|---|---:|---:|---:|---:|
| DAVIS | 0.0015 | 0.0034 | 0.0015 | approximately 0 |
| KIBA | 0.0005 | 0.0007 | -0.0002 | approximately 0 |

Point estimates were mixed and negligible; every interval covered zero, and the
protein-information stop gate failed on both the primary and replication rules.

**Confirmed proximal reason.** Direct protein coordinates are constant across
queries of one target, so they cannot improve within-target ranking by
themselves. Ranking requires the interaction module to turn protein differences
into ligand-dependent changes. Protein shuffle and the separately trained
no-interaction arm showed that this path had not become useful.

## Root-Cause Assessment for the Final Null

The evidence establishes the following immediate causes:

1. Biological `z` had about 15 to 20 times less full-state query variation than
   the old statistic.
2. Query-conditioned adaptation variation was near zero.
3. Support compound-label pairing, protein sequence, and explicit interaction
   could be changed or removed without reliable degradation.
4. Replacing CSMO with a direct simplex MLP did not consistently rescue
   performance, so the evidence does not isolate CSMO as the primary failure.

The following deeper explanations are plausible but not individually proven by
this experiment:

- The ligand GINE and protein projection were frozen from random initialization.
  ESM features were pretrained, but their task projection was not.
- Only the interaction/support/adaptation modules received 120+120 steps. This
  may be insufficient to learn a target-specific interaction from the frozen
  random ligand tower.
- The regularized band objective does not directly reward CI, Spearman, or
  NDCG, while the reported scalar is an external midband readout. A representation
  can reduce the band objective without learning a useful ranking.
- The 28-D bounded statistic and sparse 2-D CSMO views may compress already weak
  interaction differences. This is a representation-interface hypothesis, not
  a reason to redesign the frozen CSMO without another controlled study.
- Both benchmarks have complete ligand overlap across target splits. The old
  Morgan similarity statistic therefore starts with a strong ligand-geometry
  advantage even though it lacks protein information.

## Reproducibility and Artifacts

- Canonical remote environment: Python 3.12.3, PyTorch 2.5.1+cu124.
- Canonical GPU: NVIDIA GeForce RTX 4080 SUPER.
- The earlier local run used Conda `drug`, Python 3.11.15, PyTorch 2.6.0+cu124,
  and an NVIDIA GeForce RTX 4060 Laptop GPU. It is archived under
  `report/phase3_archive/pre_remote_rerun_20260804`.
- Model dtype: float64; CUDA execution fails closed on CPU.
- Frozen theory SHA-256:
  `3d660448a585662083979c198d42258466cdcca7e0aab197095800cc2d42501e`.
- Raw seed results: `report/phase3_runs/phase3_{DAVIS,KIBA}_s{11,23,37}.json`.
- Aggregate machine-readable result: `report/_archive_sources/2026-08-05_pre_master/phase3_results.json`.
- Representation report: `report/_archive_sources/2026-08-05_pre_master/representation_sanity.md`.
- Adaptation and stop-gate report: `report/_archive_sources/2026-08-05_pre_master/meta_adaptation_audit.md`.

No file under `theory/` was modified, and neither CSMO nor the Band Operator was
redesigned.

## Remote Transition Failures

### F-12: The official Hugging Face endpoint was unreachable remotely

**Attempt.** Download all declared ESM-2 snapshots directly from
`huggingface.co` on the new training host.

**Result.** HTTPS connection attempts to `huggingface.co:443` timed out. The
Hugging Face-compatible mirror was reachable, but the first mirrored download
still delegated large files to Xet CAS and failed with HTTP 401.

**Reason.** The remote network could not reach the official endpoint, and the
mirror did not provide anonymous access to Hugging Face's separate Xet CAS
service.

**Action.** Make the mirror endpoint explicit in the bootstrap command and set
`HF_HUB_DISABLE_XET=1` before importing the client. Ordinary mirrored HTTP
downloads then succeeded. All three remote safetensors match the local
reference SHA-256 values; no weights were uploaded from the workstation.

### F-13: Unfiltered snapshot download fetched redundant weight formats

**Attempt.** Use the default `snapshot_download` file selection.

**Result.** It fetched safetensors, legacy PyTorch binaries, and TensorFlow
weights, growing the cache to 8.2 GB while two redundant 650M files were still
incomplete.

**Reason.** A repository snapshot contains multiple framework formats, while
the offline Transformers loader needs only one model format plus tokenizer and
configuration files.

**Action.** Restrict downloads to `model.safetensors`, `config.json`,
`special_tokens_map.json`, `tokenizer_config.json`, and `vocab.txt`. The
dry-run-first pruning utility removed only undeclared snapshot entries and
unreferenced blobs, reducing the validated cache to 3.1 GB.

### F-14: The first remote Phase 3 attempt lacked scikit-learn

**Attempt.** Run DAVIS seed 11 in the newly provisioned remote environment.

**Result.** Representation staging completed, then construction of the
old-statistical baseline stopped with `ModuleNotFoundError: sklearn`. No seed
JSON had been overwritten.

**Reason.** `scikit-learn` was an implicit dependency of the historical PCA
baseline and was absent from `config/remote/requirements.txt`.

**Action.** Declare and install `scikit-learn>=1.5`, rerun 13 tests, and restart
the seed from the beginning. The complete remote rerun then passed all tests and
theory-hash checks.

### F-15: The uploaded checkpoint was labeled ESM3 but was ESM-2

**Attempt.** Treat the manual `weights/manual_upload/esm3` intake as a newly
available ESM3 representation.

**Result.** Its repository path, configuration, snapshot commit, byte size, and
SHA-256 identify `facebook/esm2_t33_650M_UR50D`. The 2,609,506,392-byte
`model.safetensors` has SHA-256
`a08adabb949fa67ad3c14b509d04fd60368b35007b0095e3358f81200c4f4db0`,
identical to the active ESM-2 650M checkpoint.

**Reason.** A complete ESM-2 Hugging Face cache tree was uploaded into a folder
named `esm3`; the folder name does not determine model architecture.

**Action.** Keep the upload quarantined and run the declared Phase 3 protocol
with the validated ESM-2 cache. No ESM3 result or improvement is claimed.

## Phase 4 and Phase X Recovery Failures

### F-16: Cross-attention plus label-keyed memory did not recover adaptation

**Attempt.** Replace only the biological-to-meta interface with a residue-atom
cross-interaction encoder, label-keyed query-conditioned support memory, and an
auxiliary ranking/association objective. CSMO, the Band operator, and the
frozen theory were unchanged. DAVIS and KIBA used strict unseen-target splits
and seeds 11, 23, and 37.

**Result.** Phase 4 returned `INTERACTION_SIGNAL_NOT_FOUND`. Correct biological
input did not reliably beat support-label permutation, fixed/mean memory,
protein shuffle, or the old statistical state. On DAVIS, the full interaction
path was worse than removing interaction for CI by `-0.0545` with 95% interval
`[-0.1089,-0.0002]`. Ligand shuffle was damaging, while protein shuffle was
neutral. KIBA showed the same null protein and support-label pattern.

**Reason.** The new modules increased latent variation but did not give it
protein-conditioned or label-binding semantics. Softmax memory always assigned
unit mass to five supports, even when none was relevant; raw labels were mixed
without a source prior, residual target, absolute coverage, or coherent
task-specific response function. The interaction fusion also degraded the
useful ligand signal.

**Action.** Keep this implementation isolated under `research/`. Do not promote
it to `model/` or retry it by adding attention blocks, increasing width, or
tuning hyperparameters.

### F-17: Structural support dependence was not useful meta-adaptation

**Attempt.** Run a seed-11 source/metaval boundary audit on the frozen Phase 4
path. Cross four queries with 20 legal support draws per target; intervene on
label binding, support features, support size, wrong-target support, and protein;
measure changes at memory, `z`, simplex, band, and scalar output. Compare a
fixed Morgan-kernel support estimator as a non-neural data-sufficiency control.

**Result.** Both datasets passed order invariance, query-label exclusion, and
nonzero support gradients, but failed protein use, material query-support
interaction, and useful-support gates. Learned label permutation changed the
scalar output by only `0.0015` source-label standard deviations on DAVIS and
`0.0008` on KIBA. Crossed query-support interaction RMS was `0.0154` and
`0.0100`, below the preregistered `0.05` margin. The interaction representation
was worse than ligand-only on both datasets: DAVIS CI difference `-0.0674`
`[-0.0795,-0.0548]`; KIBA `-0.0297 [-0.0572,-0.0054]`.

The fixed Morgan kernel did contain support association signal: label
permutation reduced DAVIS CI from `0.6228` to `0.5037` and KIBA CI from `0.5298`
to `0.4857`. Its relative RMSE benefit was not significant, so this establishes
recoverable ranking information, not a calibrated-affinity gain.

**Reason.** The learned path was sensitive mainly to support count and
wrong-task context, not to correct query-specific compound-label association.
The current protein geometry was low rank, and the interaction module added no
protein-dependent predictive value. The association objective remained near
chance, so longer training is not justified without a separate trainability
test.

**Action.** Return `CHEAP_GATES_FAILED`. Prohibit another full run or capacity
increase. Correct the diagnostic estimands before a future passing result can
authorize training: use at least 12 queries per target, score each support draw
instead of a support ensemble, make feature replacement disjoint from original
supports, gate on legal-support effects, and enforce every leakage assertion.

### F-18: Fixed sequence-physicochemical pair features did not repair Gate 2

**Attempt.** Test the earliest failed Phase Y information boundary with a
source-only, low-capacity sequence-physicochemical x Morgan interaction
fingerprint. A ridge probe and a capacity-matched ligand-only probe were
selected by homology-grouped source CV, then evaluated on strict metaval data
with 20 deterministic cross-component protein derangements. CSMO, Band, the
residual process, recipient targets, and `theory/` were not used or changed.

**Result.** The first descriptor variant used ligand physicochemical channels
instead of the preregistered Morgan interaction and was archived rather than
promoted. The corrected `PhysChem80 x Morgan16` probe also failed. On DAVIS,
pair-minus-ligand CI was `+0.0016 [-0.0225, 0.0258]` and correct-minus-protein-
deranged CI was `+0.0258 [-0.0033, 0.0565]`. On KIBA the contrasts were
`-0.0847 [-0.1159, -0.0533]` and `-0.0217 [-0.0639, 0.0203]`. None meets the
required `+0.03` with a positive DAVIS lower bound and favorable KIBA
direction. DAVIS same-ligand/across-protein CI was above 0.5 by `+0.0729
[0.0377, 0.1047]`, but KIBA was `-0.0332 [-0.1302, 0.0663]` and the primary
target-macro tests failed.

**Reason.** Global and local sequence composition supplies weak protein-varying
features, but without binding-site/contact supervision it does not produce a
ligand-conditioned target response beyond the Morgan ligand baseline. The
narrow crossed-axis result cannot rescue a representation that is neutral or
adverse on the preregistered primary task-level tests.

**Action.** Record `PAIR_GATE_NOT_REPAIRED`. Do not deepen cross-attention,
increase ESM/GNN capacity, implement the residual process, or connect CSMO.
The next authorized boundary is data supervision: construct the capped,
assay-aware ChEMBL37 source pilot and determine whether clean contact/site or
assay-context tasks exist before attempting a supervised residue-atom feature.

### F-19: Combined background launch started only the first pair-probe job

**Attempt.** Start DAVIS and KIBA in one remote `nohup` shell expression.

**Result.** DAVIS completed, but shell backgrounding/precedence left the SSH
session attached and the KIBA command was not started.

**Reason.** The combined remote command did not isolate both background jobs as
independent shell groups.

**Action.** Terminate the attached session, run KIBA explicitly in the
foreground with a dedicated log, then summarize only after both raw JSON files
exist. Both datasets completed successfully and the focused remote suite passed
`6` tests.

## Phase Z Recoverability Failures and Repairs

### F-20: Named controls did not implement their declared mechanisms

**Attempt.** Reuse the Phase 3 control registry and collapse baselines as the
starting point for recoverability experiments.

**Result.** The arm named `uniform_routing` selected the learned stick-breaking
fallback, three collapse baselines referenced stale 14-dimensional coordinates,
and constant `gamma=0.5` made only 6 of 12 declared contexts reachable.

**Reason.** Control names were not tested against constructed mechanisms after
the statistic expanded from 14 to 28 dimensions, and the deployment context
alphabet was not reconciled with the observable data schema.

**Action.** Use a parameter-free fixed uniform gate, remap controls to declared
28-D semantics, reject unknown gate parameterizations, and declare six reachable
contexts. Remote verification passed 13 tests. This is a new deployment, so old
12-context `B` tables/checkpoints are invalid and `B_Z` remains unbuilt.

### F-21: The first hard seal exposed recipient ligand covariates

**Attempt.** Emit source/metaval labels in `permitted.tsv` and retain all-split
row governance without a label column.

**Result.** Labels were absent, but recipient drug IDs and SMILES were still
visible row by row in `governance.tsv`. In-process counters also could not prove
that code had not bypassed the loader to read `dataset/raw`.

**Reason.** The first design interpreted "label seal" too narrowly and retained
more recipient metadata than target/homology governance requires.

**Action.** Replace governance with one target-level row containing only target
ID, sequence, sequence hash, and split. Run Phase Z under unprivileged UID 1000,
make raw tables root-only mode `0700`, and keep sealed manifests root-owned
read-only. Raw-read and manifest-write probes now fail with `PermissionError`;
all four recipient counters are zero. The overbroad v1 seal is not mounted.

### F-22: Sequential ChEMBL screening was operationally non-decisive

**Attempt.** Screen the deterministic 2,500-assay ChEMBL 37 candidate window
sequentially while retaining at most 200 eligible assays.

**Result.** Only about five candidates per minute were screened; the checkpoint
reached 9 accepted assays and 265 rows after roughly 90 decisions. Extrapolation
was several hours before task-sufficiency could even be tested.

**Reason.** Every clearly ineligible assay incurred sequential target,
component, and endpoint activity requests before the cheap count condition was
known.

**Action.** Add a bounded parallel necessary-condition prefilter. Ki, Kd, and
IC50 counts remain separate, and the original standardized connectivity-distinct
acceptance predicate is still rerun. This changes acquisition efficiency only,
not scientific eligibility.

### F-23: The first parallel prefilter lost progress on a disconnect

**Attempt.** Compute all 2,500 assay count screens with 12 worker threads and
write one aggregate file after completion.

**Result.** The official service closed one HTTPS connection after about seven
minutes. `RemoteDisconnected` escaped the retry tuple, and no aggregate
prefilter artifact had been written.

**Reason.** Retry handling covered `URLError` and timeouts but not lower-level
connection resets; aggregation was atomic only at the end rather than resumable
per assay.

**Action.** Retry bounded disconnect/reset exceptions, cache each successful
assay atomically, resume only missing assay IDs, reduce concurrency to six, and
add disconnect plus interruption/resume tests. Remote validation passes 10
focused tests. The partial 9-assay state is not treated as a Gate Z1 result.

### F-24: K-mer-prefiltered homology closure was not a sufficient audit

**Attempt.** Use the repository homology clustering routine and report zero
straddling from the split it constructed.

**Result.** The zero count was partly tautological, and the routine skipped
alignments below a 5-mer containment prefilter. Local 40% identity can exist
despite low global 5-mer containment.

**Reason.** A computational shortcut was promoted into a leakage assertion
without an independent cross-boundary alignment check.

**Action.** For the capped ChEMBL protein set, align every unique protein pair,
union every pair at or above 40% identity before splitting, and independently
report the maximum cross-split identity and count at or above 0.40. Gate Z1 now
fails closed if the exhaustive audit is absent or nonzero.

### F-25: Exhaustive post-acquisition audit used the wrong alignment keyword

**Attempt.** Run all-pairs Smith-Waterman closure after the capped pilot
finished acquiring 141 assays.

**Result.** Acquisition and its atomic checkpoint completed, but audit stopped
before writing a Gate verdict with `TypeError: unexpected keyword argument
'align_max_len'`.

**Reason.** The compatibility wrapper exposes `max_len`; `align_max_len` is the
name used by its caller in the homology-clustering layer, not its public
argument.

**Action.** Change the call to `max_len=None` and rerun the read-only audit from
the completed accepted-assay artifacts. No row, task, split, or model result was
produced by the failed audit.

### F-26: The first capped assay-aware pilot lacked context diversity

**Attempt.** Use a deterministic 2,500-assay ChEMBL 37 census, retaining at
most 200 assays, to establish an assay-aware source/metaval task population.

**Result.** The independently verified first iteration retained 141 assays and
4,463 measurements. It produced 114 legal source tasks and 89 source proteins,
but only 20 source proteins had two or more legal assay contexts. Gate Z1
therefore failed its thresholds of 128 tasks and 32 multi-context proteins.
All unit, censoring, provenance, document-closure, exhaustive homology-closure,
and recipient-seal checks passed.

**Reason.** The evenly spaced assay census found enough distinct proteins but
not enough repeated assay contexts for the same proteins. The missing
information dimension is assay-context diversity, not model capacity or basic
data integrity.

**Action.** Preserve the first pilot by hash and run one preregistered,
target-balanced second-context census for already observed proteins. The repair
must add at least 14 legal source tasks and move at least 12 additional source
proteins to two legal contexts within the 59 remaining accepted-assay slots.
No model or recipient experiment is authorized until the complete Gate Z1
audit is rerun.

### F-27: Gate Z1 passed, but the Gate Z2 cohort lacks independent components

**Attempt.** Run the preregistered target-balanced second-context expansion
within the remaining 59 accepted-assay slots, then promote directly to the k=5
support-value experiment if the corpus became sufficient.

**Result.** The expansion reached 200 assays and 6,910 measurements. Gate Z1
passed with 157 legal source tasks and 47 multi-assay source proteins. The
stricter Gate Z2 census found exactly 32 metaval tasks with at least 19 exact
compounds, but these tasks occupied only 12 independent homology/document
closure components versus the preregistered minimum of 16.

**Reason.** The repair intentionally enriched second assay contexts for the
same 112 proteins. It increased task multiplicity but could not add independent
protein components, so task count and statistical independence diverged.

**Action.** Do not fit the three support estimators yet. Progressively expand
only protein/component diversity using a deterministic protein-balanced source
census, then rerun label-free Gate Z2 eligibility. Do not interpret overlapping
support draws or same-component tasks as independent evidence, and do not start
MRFA, pair-model scaling, or CSMO training before this boundary passes.

### F-28: Reusing the old assay census could not repair target diversity

**Attempt.** First inspect only unprocessed assays from the frozen 5,000-assay
census, then broaden to every unaccepted count-eligible assay while retaining
the original 40% homology, document and split rules.

**Result.** The unprocessed tail contained zero novel legal metaval components.
The complete unaccepted set contained only one, below the six-component
registered target. Both attempts stopped before activity-value acquisition.

**Reason.** The second half of the old census had been deliberately constructed
to find additional contexts for already observed proteins. Reusing it could
increase assay multiplicity but not independent target closure.

**Action.** Build a target-first census from all 12,912 release-pinned ChEMBL37
protein components. Sequence/document/split checks preceded endpoint counts,
and the acquisition plan was frozen before affinity rows were fetched. Six
novel singleton metaval components were registered for full normalization.

### F-29: The component-expansion verifier relaxed a registered survivor rule

**Attempt.** Treat the V1 component expansion as complete after five of the six
registered singleton components survived normalization, because the resulting
37 tasks in 17 components exceeded the older absolute 32/16 cohort floor.

**Result.** `CHEMBL995610` failed the exact standardized-compound requirement,
leaving five accepted components. The producer checked only the older 32/16
floor, and the independent verifier rejected only fewer than four surviving
components. Both therefore reported a pass despite the registered requirement
that six new components survive normalization.

**Reason.** Candidate-plan sufficiency was conflated with post-normalization
survival, then an older global threshold was substituted after attrition. The
37/17 closure arithmetic remains credible, but the V1 repair hypothesis failed
5/6 and cannot authorize estimator fitting.

**Action.** Preserve V1 and its hashes as
`EXPANSION_HYPOTHESIS_FAILED_5_OF_6`. Run one prospective label-blind
continuation in the original salted candidate order, stop after the first
additional singleton component survives, and independently require cumulative
counts of at least six accepted extension components, 38 eligible tasks, and
18 closure components. Freeze the analysis cohort before scoring labels.

### F-30: Gate Z2 freeze v1 was computationally non-viable

**Attempt.** Freeze every task, query, support draw, feature nuisance and
wrong-task donor with the first standalone Gate Z2 implementation.

**Result.** The label-blind process used one CPU core for 25 minutes without
producing a manifest and was terminated by exact PID. Code SHA-256 was
`2ad52c8d479c72ee761a3539734310909520ffa7ff04a746ebd0df7f43a6267a`;
no affinity label or recipient file was read.

**Reason.** The fixed protein-ligand CountSketch interaction feature was
recomputed for the same `(protein, ligand)` pair throughout donor and coverage
construction.

**Action.** Preserve runtime v1 as an engineering failure and memoize the exact
deterministic interaction feature. Verify cache identity by object reuse and
numerical equality before creating a new immutable runtime.

### F-31: Gate Z2 freeze v2 over-constrained wrong-task assay context

**Attempt.** Run the cached label-blind freeze while requiring exact BAO-format
equality for every wrong-task source donor.

**Result.** Freeze completed but retained only 28 tasks in 12 components,
below the registered analysis floor. Runtime v2 recorded code SHA-256
`9a6cc98b536fa0ba315d66b6b80257683359dce3e132a4a84b80fe982284a631`
and `recipient_label_reads=0`; scoring did not start.

**Reason.** Six tasks in four components had no source donor with the exact BAO
format. Exact BAO equality was stricter than the preregistered requirement for
a compatible observable assay context and made the control structurally
unavailable.

**Action.** Preserve runtime v2. Replace opportunistic BAO equality with a
prospectively documented compatibility map and require endpoint, assay type,
numeric-context tolerances, missingness rules and donor-balance audits before
labels are joined.

### F-32: Gate Z2 freeze v3 selected away four independent components

**Attempt.** Use compatible-context source donors and freeze a shared
high-coverage primary population for both non-neural estimators.

**Result.** The verified continuation input contained 38 tasks in 18 closure
components. Five tasks in four components could not produce 20 draws satisfying
the shared high-coverage rule, so the frozen analysis fell to 33 tasks in 14
components. Runtime v3 code SHA-256 was
`32a9a02088ea6e7d4e016b427ab067ed576f0a9a8e088b69e3c3efa57bc5e275`;
all four immutable input hashes matched the previous endpoint,
`recipient_label_reads=0`, and scoring did not start.

**Reason.** Coverage was used as a primary eligibility filter. This changed the
registered 38/18 population and repeated the threshold-substitution failure
that F-29 was intended to close.

**Action.** Preserve runtime v3 with verdict
`GATE_Z2_NOT_RUN_INELIGIBLE:ANALYSIS_COHORT_BELOW_38_18`. Gate Z2 v4 must keep
all 38/18 tasks, all 12 salted queries and the first 20 salted k=5 draws;
coverage is continuous diagnostic information and zero relevance returns the
prior exactly.

### F-33: Pre-score adversarial review found four Gate Z2 estimand defects

**Attempt.** Independently review the v3 scorer before making the separately
stored label file readable.

**Result.** No score was run. The review found that wrong-task source residuals
used an in-sample final prior, estimator B reused an ECFP-only donor without
matching its kernel coverage, relative RMSE was formed after component
aggregation, and nuisance fitting did not separate source-train from
source-calibration components.

**Reason.** The implementation was structurally close to the written protocol
but several nonlinear/provenance details were not encoded as executable
invariants. Each could bias the correct-versus-wrong contrast toward a pass.

**Action.** Freeze the multi-agent repair in
`report/phase_z/07_new_materials_multi_agent_synthesis.md`. Require
component-OOF donor residuals, estimator-specific context/coverage matching,
draw-level paired relative RMSE, and a disjoint source train/calibration split
with regression tests before any label join.

### F-34: Gate Z2 v4 donor failure did not identify a source-data gap

**Attempt.** Run the repaired label-blind v4 freeze on all 38 tasks, 18 closure
components, 12 fixed queries and 20 k=5 draws, with component-OOF donor
residuals, estimator-specific wrong-task matching, draw-level relative RMSE and
a disjoint source train/calibration split.

**Result.** Runtime v4 stopped before label access with
`GATE_Z2_NOT_RUN_INELIGIBLE:WRONG_TASK_MATCH_NOT_IDENTIFIABLE` for estimator B
and task `1529293c22a9c1dfd2fa99cb763198d804d27f6a1bd5bb853635285091545010`.
Code SHA-256 was
`cd6fecb624d8e244ce898ad46d703b37081002d59ce9ab83f169e792e55bac27`;
`recipient_label_reads=0`, DAVIS/KIBA remained unreadable, and scoring did not
start. A label-free reconstruction found ten legal donor tasks in eight source
closure components. The best inspected B candidate for draw 0 exceeded only
`D_mass` (`3.3031078` versus frozen `1.9169567`).

**Reason.** The failure was not a donor-absence certificate. `D_shape` compared
support columns in stored order even though the smoother is invariant to a
joint support/residual permutation; an equivalent column permutation produced
`D_shape=0.6324555` instead of zero in a minimal reproduction. The search also
tested at most one ECFP-optimal plus twenty salted k=5 subsets per donor, while
legal donors contained 18-68 compounds (8,568 to 10,424,128 possible subsets).

**Action.** Preserve runtime v4 as `DONOR_SELECTION_PROTOCOL_NOT_IDENTIFIED`.
Before any label join, run one prospective v5 label-blind feasibility repair:
keep all v4 context rules, population and numeric tolerances fixed; make
`D_shape` permutation invariant; search donor supports deterministically with
an explicit exactness/search audit; and distinguish a feasible witness,
inconclusive search, and proved source-closure failure. Do not download new
source data unless the repaired search first establishes a genuine gap.

### F-35: Gate Z2 v5 serial feasibility audit underused remote CPU

**Attempt.** Run the complete permutation-invariant `38 x 20 x 2` donor
feasibility audit as one process on the replacement remote host.

**Result.** The process remained healthy at one full CPU and about 270 MiB RSS
for 1,893 seconds, but used only one of 128 logical CPUs and produced no partial
or final freeze artifact. It was terminated by its exact PID before label
access. Runtime v5 remains preserved; `recipient_label_reads=0` and scoring
never started.

**Reason.** The v5 repair correctly accumulated every cell before atomic
output, but its task loop was serial. Estimator-native beam and swap searches
are independent across metaval tasks, so the implementation left the remote
hardware idle and could not meet the requested parallel experiment schedule.

**Action.** Keep the v5 search mathematics and registered thresholds unchanged.
Implement a deterministic process-parallel v6 over metaval tasks, aggregate
results in sorted task order, regression-test worker-count invariance, and run
it in a fresh immutable runtime. Do not treat the terminated serial attempt as
a gate result.

### F-36: Gate Z2 v6 masked a definitive context gap in mixed failures

**Attempt.** Start the process-parallel 16-worker freeze after local and remote
tests passed.

**Result.** An independent adversarial code review completed about 86 seconds
after launch and found that mixed failure sets prioritized
`WRONG_TASK_SEARCH_INCOMPLETE` over `ASSAY_CONTEXT_DONOR_MISMATCH`. The main
process and all orphaned workers were terminated by their exact runtime path;
no freeze or score artifact existed, labels remained unreadable, and
`recipient_label_reads=0`.

**Reason.** Internal diagnostics retained per-cell causes, but the top-level
exception selected only one code. The branch order allowed a non-definitive
heuristic-search miss to mask an earlier, proved source context-closure gap.

**Action.** Preserve runtime v6 as a preflight termination. In v7, retain all
observed cause codes and order the primary boundary as context absence, then
exhaustive geometry infeasibility, then nonexact search incompleteness. Add a
mixed-failure regression test before a fresh parallel run.

### F-37: Directory grouping did not satisfy model code consolidation

**Attempt.** Replace the flat model package with `model/biology/` and
`model/core/` subpackages while preserving all implementations and imports.

**Result.** The layout passed 101 local tests and 14 remote focused tests, but it
still contained the same number of implementation files. It improved ownership
but did not meet the requested reduction in code fragmentation.

**Reason.** “Consolidate” was interpreted as directory organization rather
than physical source merging. The user explicitly required fewer files directly
under `model/`, without new subdirectories.

**Action.** Remove the subpackage layout after validation. Merge the six
biological/support modules into `biological.py`, CSMO and Band into
`meta_operator.py`, and five mathematical helpers into `mathematical.py`.
Retain only the four shared boundary files. The final local package has seven
implementation files and passes 103 tests; frozen theory is unchanged.

### F-38: Full-suite verification initially failed at the sandbox temp boundary

**Attempt.** Run the complete consolidated-model test suite with pytest's
default temporary directory, then with a new `C:\tmp` base directory.

**Result.** In both attempts, 73 tests passed and 29 tmp-path-dependent tests
errored before setup with Windows `PermissionError`. No assertion or model
failure occurred.

**Reason.** The managed execution identity could not enumerate the existing
user pytest temp root and could not create the requested `C:\tmp` directory.
This was a filesystem ACL boundary, not a code regression.

**Action.** Create a dedicated writable directory under `D:\MetaSieve` and
rerun with that exact `--basetemp`. The complete suite then passed. A final
lazy-import regression test raised the clean total to 103 tests.

### F-39: Gate Z2 v7 established a source-context gap but not a biological result

**Attempt.** Run the immutable label-blind v7 feasibility audit for every one
of the registered 38 tasks, 20 support draws and two source-tuned estimators
with 16 deterministic workers and corrected mixed-cause precedence.

**Result.** All 1,520 cells were recorded. Of these, 1,286 had a feasible
matched wrong-task donor and 234 were ineligible. Eighty cells, covering both
estimators and every draw for two metaval tasks, had zero compatible source
assay-context donor. A further 154 cells returned
`WRONG_TASK_SEARCH_INCOMPLETE`; because their beam-and-swap search was not
exhaustive, those cells do not prove source-data absence. The primary verdict
was `GATE_Z2_NOT_RUN_INELIGIBLE`, `recipient_label_reads=0`, and no affinity
score was run. The synchronized failure JSON SHA-256 is
`1e93d44de3104bbfafdb005d4c82b2f2b2d9aa7dac49c7d506e258fdd38e970a`.

**Reason.** The registered source corpus lacks compatible observable assay
contexts for two tasks. Separately, the deterministic heuristic did not find a
witness for 154 cells and did not claim global exactness. Pooling these causes
would falsely turn an engineering search limit into a data-gap certificate.

**Action.** Keep every support and recipient label sealed. Acquire only
source-side tasks matching the two missing observable contexts under the
existing provenance and closure rules, and repair search completeness
separately. Do not interpret this result as support failure or architecture
failure.

### F-40: Literature plausibility was not accepted as architecture validation

**Attempt.** Use five independent literature/material reviews to decide whether
the proposed biophysical pair field, residual support adapter and typed state
were promising.

**Result.** The reviews identified transferable mechanisms but produced no
MetaSieve effectiveness evidence. They therefore cannot authorize integration
into `model/`.

**Reason.** A paper precedent establishes plausibility under the paper's data
and estimand; it does not establish that MetaSieve transmits protein-specific
or support-specific information under its frozen theory and closure rules.

**Action.** Start remote reproduction of the official PSICHIC PDBbind v2020
benchmark first, pinned to the authors' code, data, configuration and released
weight. Continue source-only Gate Z2 closure repair in parallel. Promote no new
model code until the registered pair and support information boundaries pass.

### F-41: PSICHIC released-weight reproduction passed, but is not a MetaSieve claim

**Attempt.** Evaluate the authors' released PSICHIC PDBbind v2020 checkpoint on
their temporal test split using the pinned official repository, split,
configuration, degree tensors and checkpoint on the replacement remote GPU.

**Result.** All 363 test rows were scored on CUDA. RMSE was 1.245707, CI
0.760739, Spearman 0.711999, Pearson 0.734565 and MAE 0.954369. The checkpoint
SHA-256 was
`c0d577efff5ccfc6a025d22630a0ff37575420bb94ef3c05dff865bf8ca7811d`;
the prediction CSV SHA-256 was
`2ea7f3860f5052d2c80beef4e8cf4faa0025ff84dc5e5012aaad1edb2d3f55e5`.
No DAVIS/KIBA recipient label was read.

**Reason.** This establishes an executable, interaction-aware published
baseline in the actual remote environment. It does not test MetaSieve support
conditioning and therefore cannot establish the proposed three-module
architecture by itself.

**Action.** Preserve row-level predictions and metrics locally. Reproduce the
paper's Human cold-protein setting next, then use the published baseline only
after the source/support gates authorize a matched MetaSieve arm.

### F-42: First Human cold-protein launch failed before execution

**Attempt.** Create the official cold-protein split and launch training through
one SSH heredoc.

**Result.** Shell delimiter quoting placed the launch statements inside Python,
which raised `SyntaxError`; no split or training process was produced.

**Reason.** The multi-language heredoc boundary was malformed.

**Action.** Reissue the deterministic split as a single Python command, audit
zero protein overlap, and launch the official training command separately.
The corrected run is active as remote PID 21662.

### F-43: Source-context repair plan found too few legal candidates

**Attempt.** Search ChEMBL 37 metadata within the already registered source
closure components for eight Kd binding assays with cellular material, at
least 19 exact standardized measurements, no cross-protein document overlap,
and no recipient access.

**Result.** The remote atomic plan screened 10,477 metadata candidates but
retained only four. Exactly 10,473 failed the Kd count threshold. The frozen
target was eight, so acquisition was not started. Recipient label reads were
zero. The synchronized plan SHA-256 is
`4082fde1178188a91fae26ac5e114721761e8ebf15f1c5fabad3a4d607b1b8cd`.

**Reason.** The existing registered source components do not contain enough
independent high-count Kd cellular-context assays. This is a source population
coverage failure, not evidence that support has no affinity value.

**Action.** Do not lower the registered compound-count or task-count
thresholds. Prepare an identity-only expansion over additional independent
source closure components and require at least eight eligible metadata
candidates before opening any source activity values.

### F-44: Full single-process regression was polluted during collection

**Attempt.** Run the complete remote project suite in one pytest process after
the GPU data-path equivalence check.

**Result.** Collection stopped because `sklearn` and `scipy` had been replaced
by non-package modules before two Phase X/Y tests imported their submodules.
The affected tests passed 6/6 in a clean process. Running every test file in an
independent Python process passed all 17 files and 109 tests.

**Reason.** A test-collection global-module side effect contaminates the
single-process environment. It is not a model, tensor, prediction or training
regression, but the one-process suite cannot currently serve as clean evidence.

**Action.** Retain isolated-file regression as the remote acceptance check and
separately locate the module-polluting collector before restoring a one-process
suite. Do not weaken any scientific gate because of this harness defect.

### F-45: PSICHIC success does not establish MetaSieve support adaptation

**Attempt.** Reproduce the official PSICHIC Human cold-protein classification
benchmark remotely using its published code path and frozen split.

**Result.** The 30-epoch run completed with test AUROC 0.955081, AUPRC
0.957928 and F1 0.869888 over 1,375 rows. This validates a transferable
protein-ligand pair representation baseline under unseen proteins.

**Reason.** PSICHIC does not condition a residual response function on k<=5
measured affinities. Its success therefore cannot answer whether correct
support-label binding improves MetaSieve affinity risk or survives the z/CSMO
handoff.

**Action.** Use the reproduction as the pair-representation reference only.
Continue the registered knowledge-token and governed contact-teacher tests, and
do not promote either into `model/` until protein derangement, ligand-only and
support controls pass.

### F-46: Direct K1 GO-MF supervision did not create target-specific pair value

**Attempt.** Train the same MetaSieve-specific ESM/GINE/atom-residue prior with
and without source-only GO-MF supervision of its existing functional residue
slot. Run three matched seeds on 49 metaval target-assay tasks and use paired
task-block bootstrap.

**Result.** K1 minus baseline CI was +0.0008 with 95% interval
[-0.0197, +0.0189]. RMSE was worse by 0.0306. Correct-protein minus deranged
CI averaged only +0.0067, below the registered +0.03 threshold. All execution
tensors and optimizer states were CUDA resident, so the failure is not a CPU
fallback or loader artifact.

**Reason.** The 89 source proteins and coarse protein-level GO-MF labels do not
provide enough ligand-conditional local supervision. They can describe broad
function without identifying which residue chemistry matters for a query
ligand. A larger KG branch or a paper checkpoint would add capacity without
repairing this missing correspondence.

**Action.** Do not integrate K1 into `model/` and do not run K1+S1. Use the
uploaded ProteinKG25 or a governed UniProt/InterPro corpus only to train the
same sequence student under homology closure, then repeat this exact direct
MetaSieve gate.

### F-47: Residue-site supervision improved the weak pair baseline but not the ligand shortcut

**Attempt.** Supervise the existing atom-residue field with source-only UniProt
binding/active-site residue labels, with no online structure model or embedding
concatenation.

**Result.** S1-site improved RMSE over the weak pair baseline by 0.1280 with a
paired interval excluding zero. It did not beat ligand-only: RMSE was worse by
0.0460 and CI was lower by 0.0094, with both intervals spanning zero.
Correct-protein minus deranged CI averaged +0.0164. The +0.037 seed-17 effect
fell to +0.0078 and +0.0042 in the other seeds.

**Reason.** A residue-only site label is ligand independent. It can regularize
where the model looks, but cannot teach which ligand atom interacts with which
residue or distinguish ligand-specific contact patterns. It is therefore not
the final contact-distillation mechanism.

**Action.** Do not integrate S1-site into `model/`. Acquire governed holo poses
with exact atom-residue contacts and rerun S1 against ligand-only and protein
derangement. Do not label the uploaded intra-protein Contact LMDB as an
atom-residue contact corpus.

### F-48: Task-count labels referred to different eligibility contracts

**Attempt.** Apply the Phase Y minimum of 12 metaval queries before paired
pair-representation scoring, expecting the earlier Gate Z2 count of 38 tasks.

**Result.** All 49 metaval target-assay tasks had at least 12 rows and remained
eligible. The earlier count of 38 additionally required Gate Z2's stricter
exact-compound and support/donor feasibility conditions.

**Reason.** The same word `eligible` had been used for two different estimands.

**Action.** Report 49 tasks for this pair gate and retain 38 only for the Gate
Z2 support-estimator protocol. Preserve both contracts explicitly.

### F-49: Inline PowerShell-to-SSH quoting failed before homology governance

**Attempt.** Launch the ProteinKG25 MMseqs closure pipeline through a nested
inline `bash -lc` command from PowerShell.

**Result.** Local shell interpolation corrupted the remote command before
MMseqs started. No scientific output was produced and no input archive was
modified.

**Reason.** The launcher mixed two shell quoting grammars in one command. This
was an orchestration defect, not a data or model failure.

**Action.** Replaced the inline launcher with the checked-in
`scripts/run_proteinkg25_homology.sh`. The rerun retained 47,736 records after
excluding 189 current-protein homologs and produced zero train/validation
homology-component overlap.

### F-50: Generic ESM cache path starved the GPU during ProteinKG25 scaling

**Attempt.** Reuse the general protein cache extractor for all 47,736 governed
ProteinKG25 sequences with ESM2-150M on CUDA.

**Result.** ESM weights occupied GPU memory, but a ten-second sample showed 0%
SM utilization while the process consumed about 14 CPU cores. The run was
stopped before it produced a cache artifact.

**Reason.** The general extractor copied every full residue window to CPU and
then performed fixed-bank pooling there, forcing a device synchronization and
many small OpenMP pooling operations per batch.

**Action.** Keep frozen ESM forward, mask-aware pooling and 16-bin adaptive
pooling on CUDA; transfer only final FP16 pooled/bank tensors to CPU. A 32-entry
old/new check had identical keys and dtypes, maximum absolute difference
0.00390625 and mean absolute differences below 3e-5. The corrected full run
showed active GPU compute and no CPU fallback.

### F-51: ProteinKG25 scale did not repair K1's ligand-independent semantics

**Attempt.** Govern ProteinKG25 against all current proteins, retain 47,736
nonhomologous proteins across disjoint train/validation homology components,
pretrain the unchanged MetaSieve protein/knowledge interface on 1,850 GO-MF
terms, and repeat the direct K1 gate for seeds 17, 23 and 31.

**Result.** The external teacher learned its task (validation BCE 0.2211 to
0.0990), but scaled K1 did not beat the pair baseline: RMSE delta +0.0070 with
95% interval [-0.1070, +0.1225], and CI delta -0.0076. It was significantly
worse than ligand-only in RMSE by +0.1810 [0.0456, 0.3138]. Correct-protein
minus deranged-protein CI averaged -0.00093.

**Reason.** Data volume was not the earliest remaining boundary. GO-MF is a
ligand-independent protein function target, so even a well-fitted teacher does
not identify which local residue chemistry is relevant to the query ligand.
The transferred interface can retain ontology information while remaining
useless for target-specific protein-ligand discrimination.

**Action.** Reject scaled K1 and do not integrate its checkpoint into `model/`.
Do not run K1+S1. The next pair experiment requires governed holo-complex
atom-residue contacts (PLINDER preferred, BioLiP2/RCSB fallback); broader
assay-context affinity work separately requires the BindingDB context trio and
Papyrus provenance closure.

### F-52: The first split-corrected run used an overbroad metaval estimand

**Attempt.** Replace universal 40% homology closure with the primary
exact-target `fewshot_core_v2` split, retain the old closure as a separate
`novel_family_stress_v1` evaluation, rebuild all governed accepted assays, and
retrain ligand-only, pair baseline and S1-site for three fixed seeds.

**Result.** Core pair baseline improved RMSE over ligand-only by -0.2383 with a
paired 95% interval [-0.3834, -0.0976], but worsened CI by -0.0354
[-0.0596, -0.0114]. Correct-minus-deranged CI was -0.0181. S1-site showed the
same RMSE/ranking conflict. Under novel-family stress, all pair/S1 advantages
and protein-derangement intervals crossed zero.

**Reason.** The first V2 GPU bundle did not encode `episode_role`, so the direct
gate evaluated support-pool and ineligible rows as well as queries. Conflicting
replicate measurements for one task-ligand were also weighted as separate
ligands. The single rolled derangement mapped two stress proteins within their
homology component. These are protocol errors in addition to the observed
calibration/ranking conflict.

**Action.** Preserve the run and mark it superseded. Rebuild a versioned bundle
with episode roles, aggregate each query ligand once, and require exact-target
plus cross-homology derangement before interpreting the split correction. This
correctness rerun does not consume one of the three authorized model repairs.

### F-53: Ranking and residual handoff repairs did not identify pair information

**Attempt.** Correct evaluation to query-only unique task-ligand groups with a
one-to-one cross-homology protein derangement, then test two preregistered
repairs without changing CSMO, Band, theory, training length or hidden width:
(A) add a 1:1 within-task relative-difference MSE; (B) preserve the existing
ligand head and add the existing pair prior as a `tanh`-bounded residual, with
no new parameters.

**Result.** Attempt A improved mean CI over absolute pair by +0.0196, but the
95% interval [-0.0031,+0.0435] crossed zero; it remained -0.0266 below
ligand-only and correct-minus-deranged CI was -0.0150. Its seed-17 positive
protein effect reversed in seeds 23/31. Attempt B reached seed-17 CI 0.5443,
only +0.00987 [-0.03081,+0.05155] over ligand-only, regressed RMSE by +0.03522,
and had unfavorable protein derangement (-0.00695). It failed the cheap gate,
so seeds 23/31 were not run.

**Reason.** The global absolute objective and destructive pair-only handoff
were real secondary defects: ranking loss partly repaired ordering, and the
residual handoff preserved ligand performance. Neither repair made the pair
term depend usefully on the correct protein. The earliest remaining boundary
is the semantic identification of ligand-specific atom-residue relations, not
loss weight, residual scale, capacity or additional epochs.

**Action.** Reject Attempts A/B; do not promote their checkpoints or tune them.
Attempt 3 remains blocked until governed PLINDER holo structures pass the
contact/distance oracle and derangement audit. Consolidated evidence is in
`report/phase_z/dual_protocol_module_failure_diagnosis.md`.

### F-54: Phase Y data validity did not meet the declared Gate 0 contract

**Attempt.** Audit the Phase Y DAVIS/KIBA episodes for strict target and
homology separation, compound/scaffold overlap, assay/document provenance,
source-only transforms, support/query disjointness, hidden target-ID use, and
recipient sealing before constructing a residual meta-process.

**Result.** Target/homology closure, source-only transforms, disjointness and
target-ID exclusion passed. Assay/document provenance failed because the
benchmark episode tables contained no assay or document identifiers. Recipient
sealing also failed because the loader materialized the full raw recipient
dataframe and labels; access counters were declared rather than enforced.
Compound and scaffold overlap were effectively complete and were disclosed,
not hidden. No recipient scoring or tuning was performed.

**Reason.** The legacy benchmark schema could not represent the registered
`target x assay/context` task or prove document closure, and the loader boundary
was weaker than the declared threat model. This is a protocol/provenance
failure, not evidence that recipient labels biased a reported metric.

**Action.** Mark Phase Y Gate 0 failed and prohibit recipient evaluation. Later
governed-source and label-blind exclusion work does not retroactively validate
the Phase Y run. Evidence:
`report/phase_y/PHASE_Y_COMPLETE_SUMMARY.md`.

### F-55: Two non-neural k=5 estimators did not establish general affinity value

**Attempt.** Before implementing the residual process, test two source-tuned
non-neural support estimators, PCA-RBF and Morgan/Tanimoto-power, using correct
k=5 support against label-deranged and source-donor support on DAVIS and KIBA.
The registered primary requirement was greater than 2% relative RMSE benefit,
a positive DAVIS lower interval, and favorable KIBA direction.

**Result.** Against label derangement, DAVIS relative RMSE benefits were only
1.16% for PCA-RBF and 0.61% for Tanimoto-power, below the 2% threshold. KIBA
point estimates were 3.12% and 2.57%, but could not rescue the failed primary
DAVIS requirement. Prior-only and nested k controls were not completed, the
source-donor contrast changed compounds and labels together, and both
estimators shared Morgan-derived ligand information. High-coverage strata and
ranking metrics showed exploratory signal, not general absolute affinity value.

**Reason.** Correct compound-label binding carries ordinal information and may
help under high chemical coverage, but the available k=5 supports did not
produce a statistically identified, dataset-general absolute-affinity benefit
under the declared estimand. The controls were also insufficient for a causal
residual-process promotion.

**Action.** Mark Phase Y Gate 1 failed. Do not build the residual posterior,
handoff state, or CSMO path from this evidence. Evidence:
`report/phase_y/PHASE_Y_COMPLETE_SUMMARY.md`.

### F-56: Phase Y Gates 3-6 were blocked, not failed model experiments

**Attempt.** The registered sequence after support-value and pair-information
validation was: cross-fitted residual-structure testing, a minimal residual
meta-process positive control, bounded handoff to `z`, and unchanged CSMO
evaluation.

**Result.** None of these modules was implemented or trained because Gate 1
and Gate 2 failed. No residuals, posterior, new `z`, deployment map, or CSMO
checkpoint were produced.

**Reason.** The staged protocol intentionally stops at the earliest failed
information boundary. Running downstream modules would confound missing input
information with architecture capacity and violate the registered gate order.

**Action.** Record Gates 3-6 as `BLOCKED_NOT_ATTEMPTED`, not scientific
failures. Do not cite them as evidence against residual processes or CSMO.
Evidence: `report/phase_y/PHASE_Y_COMPLETE_SUMMARY.md`.

### F-57: Corrected E15 confirmed calibration learning without pair information

**Attempt.** Rerun ligand-only, unchanged pair baseline, and S1-site for seeds
17, 23 and 31 on both Core and novel-family stress after fixing the E15
estimand: query-pool rows only, one mean label per task-ligand group, at least
12 unique query ligands, and deterministic one-to-one protein derangement
across exact identity and 40% homology components.

**Result.** On Core, ligand-only/pair/S1 CI was
0.5249/0.4787/0.4937. Pair improved RMSE by -0.2502
[-0.3933,-0.1094] but worsened CI by -0.0461
[-0.0717,-0.0214] and Spearman by -0.1241
[-0.1960,-0.0531]. Correct-minus-deranged CI was -0.0241 for pair and
+0.0046 for S1, with both intervals crossing zero. Under novel-family stress,
ligand-only/pair/S1 CI was 0.5268/0.5162/0.5206; neither RMSE, ranking, nor
protein intervention showed a favorable excluding effect.

**Reason.** The pair route can recover source task/assay calibration, explaining
the Core RMSE gain, but it destroys useful ligand ordering and does not use the
correct protein predictively. S1's ligand-independent residue-site label asks
every ligand for one protein to reproduce the same site target, so it does not
identify ligand-specific atom-residue relations. The lost stress RMSE benefit
rules out a transferable pair response.

**Action.** Return `PAIR_INFORMATION_NOT_IDENTIFIED`. Do not promote baseline,
S1, E16-A, or E16-B checkpoints into `model/`. The complete 10,000-draw paired
task-block result is in `report/phase_z/dual_protocol_v3_results.md`.

### F-58 [HISTORICAL; SUPERSEDED_AS_CURRENT_ACTION]: The third repair lacked selected PLINDER structure shards

> Current execution is P1A BioLiP2 + RCSB; see `AGENT_HANDOFF.md` and F-60.

**Attempt.** Freeze a label-blind PLINDER cohort for the third and final
authorized repair: exact ligand-specific atom-residue contact/distance
supervision with structure, pose, leakage, homology, ligand/scaffold, and
recipient-overlap governance.

**Result.** Metadata screening successfully froze 675 systems, 161 UniProt
targets and 51 pocket components across 24 selected shards. All 608 recipient
target records were resolved, exact recipient targets were excluded, recipient
label reads remained zero, and the machine manifest was hashed. However, none
of the required `systems/*.zip` structure shards is present; metadata alone
cannot verify atom/residue index alignment or extract the contact oracle.

**Reason.** This is a missing-input boundary, not a model failure. Residue-
residue contacts, ligand-independent site annotations, synthetic wrong-protein
affinities, or metadata-only pseudo-contacts are not valid substitutes for
exact holo atom-residue supervision.

**Action.** Record Attempt 3 as `BLOCKED_NOT_ATTEMPTED` and do not consume the
third repair attempt. Download only the 24 fixed shards (about 3.38 GB), verify
their registered GCS MD5 values and pose/index oracle, then run the preregistered
Core seed-17 cheap gate. Manifest SHA-256:
`991ac634d0c652dcbf97c5f252ef8a5e71009173b38e1c14fce814132ba9427d`.
Evidence: `report/phase_z/plinder_shard_selection.md` and
`report/phase_z/plinder_selection_v3/`.

### F-59: Frozen PSICHIC pair latent improved ranking but was protein-insensitive

**Attempt.** Run the independently approved
`R1_FROZEN_PSICHIC_PAIR_LATENT_POSITIVE_CONTROL_V1` on the corrected Core
Episode V3 cohort. Freeze official PSICHIC commit
`cc445aa28d044c6212023705208b1a7704a00622` and checkpoint
`c0d577efff5ccfc6a025d22630a0ff37575420bb94ef3c05dff865bf8ca7811d`.
Compare the same backbone's pre-interaction 200-D ligand latent with its
post-interaction, pre-scalar-head 400-D pair-conditioned latent. Apply
source-only standardization and PCA32, then identical `Linear(32,1)` probes
with 33 trainable parameters and identical initialization. Evaluate query-only
unique task-ligand groups and reuse the correct-pair probe after deterministic
one-to-one protein replacement across exact and 40% homology components.

**Preflight.** The hook contract, row mapping, CUDA path and feature ranks
passed. A 40-source-pair rank smoke plus 32 protein interventions gave centered
rank 39 for both feature families, finite features, and correct-versus-
deranged pair-latent mean L2 4.5032. Full mapping covered 3,841 source and 989
metaval query groups with zero recipient-label reads. The frozen feature cache
SHA-256 is
`d8624941e82374003784bb78de1c19cafff923f8d1eb86e715c84640516f44f3`.

**Result.** At the registered seed-17 cheap gate, ligand/pair/deranged-pair
task-macro RMSE was 1.1521/1.1086/1.1845, CI was
0.51449/0.56450/0.56484, Spearman was 0.02416/0.16990/0.16776 and NDCG@10
was 0.60918/0.64194/0.64008. Pair-minus-ligand CI was a favorable +0.05001
and relative RMSE improved 3.78%, but correct-minus-deranged CI was -0.00033.
The mandatory correct-protein condition failed.

**Reason.** The released frozen representation contains transferable
pair-distribution, ligand or context signal, but its useful ranking variation
is not tied to the correct target protein under the registered intervention.
Large latent movement therefore does not establish biological relation
semantics. Checkpoint strictness is additionally limited: exact PDB-training
overlap was 0 proteins, 22 ligands, 99 scaffolds and 0 pairs, but official 40%
homology and enumerable ESM2 pretraining-population overlap remain `UNKNOWN`.
This result applies only to the specified checkpoint, latent and cohort; it
does not prove protein-ligand relations or exact-contact supervision are
impossible.

**Action.** Return `FROZEN_PAIR_POSITIVE_CONTROL_NOT_IDENTIFIED`. Apply the
registered cheap stop: do not run seeds 23/31, novel-family stress, Round 2,
compression, support/Q-PMA, z, CSMO or Band handoff. Keep the implementation
in `scripts/` and `tests/`; do not promote it into `model/`. Round 3 remains
the independently blocked PLINDER exact-contact hypothesis from F-58. Complete
evidence: `round1_result.txt`, `round2_result.txt`,
`final_experiment_report.txt`, and
`report/phase_z/remote_experiments/round1_frozen_pair_probe_v1/`.

### F-60: Local PLINDER MLSB was the wrong mandatory geometry entry point

**Attempt.** Use the locally available PLINDER MLSB artifact as the mandatory
source for P1 atom-residue contact and distance supervision. The fail-closed
builder required a protein structure, holo ligand coordinates and auditable
atom/residue mappings.

**Result.** The local subset contains 346 receptor `input.pdb` files and ligand
SMILES, but zero holo ligand coordinate files. Synthetic holo sidecar tests,
residue-slot aggregation and CUDA bridge forward/backward checks passed. Real
R0 geometry training and the relation gate were correctly recorded as
`NOT_RUN_FAIL_CLOSED`; no P2-P4 work was authorized.

**Correction.** This result applies only to the downloaded MLSB subset. Full
PLINDER does contain holo systems, but its full footprint is unsuitable as a
mandatory local prerequisite. The result does not falsify the Mechanistic
Interaction Bridge or establish that PLINDER as a project lacks ligand
coordinates.

**Decision.** Replace the mandatory P1 data route with BioLiP2 as the
biological-relevance index and RCSB PDB/CCD as the authoritative coordinate and
chemistry source. Run a governed Pilot-10K before any scale-up. PDBbind is not
mandatory; full PLINDER is deferred to optional robustness validation.

**New gates.** P1A requires at least 10,000 valid regular-ligand holo complexes,
2,000 receptor sequences and 2,000 ligand chemotypes. P1B tests held-out
contact/distance learning and protein/ligand derangements. P1C freezes the
bridge and requires both correct-minus-ligand and correct-minus-deranged CI of
at least +0.03, positive paired lower bounds and favorable Spearman direction
on DAVIS. Only P1C may authorize P2.

**Evidence.** `report/mechanism_refactor/p1_local_structure_audit.json`. This
entry supersedes F-58 only as a statement of the current next action; F-58's
observations and hashes remain valid historical evidence.

### F-61: Open holo P1 passed geometry and failed affinity transfer

**Attempt.** Execute the revised P1A -> P1B -> P1C sequence using BioLiP2 as
the relevance index, RCSB PDB/CCD as coordinate and chemistry truth, frozen
ESM-2, the isolated low-rank bridge, and DAVIS source-to-metaval affinity
transfer. Keep recipient labels sealed and P2-P4 frozen unless P1C passes.

**P1A result.** The corrected v2 pipeline produced 14,906 governed complexes,
14,906 receptor sequences and 2,869 chemotypes after 739 protected-benchmark
homology exclusions. The final sidecar fixed a compact-residue-order error by
globally aligning mmCIF structure residues to the canonical sequence before
slotting. P1A passed; earlier v1 sidecars are superseded.

**P1B result.** On 1,477 exactly controlled structural test complexes, correct
AUPRC was 0.43885 versus 0.05149 deranged-protein and 0.23895
deranged-ligand. Correct expected-distance MAE was 1.97541 A versus 2.66531 A
and 2.13693 A. P1B passed.

**P1C result.** DAVIS target-macro CI was 0.71110 ligand-only, 0.60629 correct
mechanism and 0.60715 deranged mechanism. Correct-minus-ligand was -0.10481
[-0.13493,-0.07613]; correct-minus-deranged was -0.00086
[-0.00587,0.00411]. Correct Spearman 0.22011 was below deranged 0.22161.
P1C failed.

**Decision.** Stop at P1C. Do not scale the corpus, integrate the bridge, or
implement P2-P4. Strong structural geometry learning is not evidence of useful
correct-protein affinity transfer. Complete evidence is in
`report/mechanism_refactor/P1_EXECUTION_REPORT.md`.

### F-62: Readout defect confirmed; MIF repair recovered sub-threshold signal

**P1R0.** Freeze the P1B checkpoint and original P1C transforms, then localize
correct-vs-deranged signal and apply a consistent active-atom permutation. The
relative L2 signal was 0.962 at projected residues, 0.704 at contact logits and
0.766 at distance logits. Orthogonal PCA32 retained 48.781% of scaled delta
energy, rejecting a PCA-only explanation. The unchanged physical pair produced
mean/max Ridge changes 0.00744/0.05528 after atom-coordinate reversal. Verdict:
`P1C_READOUT_NOT_PERMUTATION_INVARIANT`.

**P1R1.** Replace index-wise coordinates with a 288-dimensional MIF: eight
ligand pharmacophore classes by six residue chemistry classes by contact plus
five distance bins. Fit a target-group OOF ligand prior and one target-centered
residual Ridge, reused unchanged for correct, `<40%` deranged, and shuffled arms.

**Result.** MIF permutation error was 5.96e-08 and source/metaval relative L2
was 0.25088/0.22127. CI was 0.71110 ligand-only, 0.71874 correct, 0.69516
deranged and 0.68994 shuffled. Correct-minus-ligand was +0.00764
[-0.00179,0.01685]; correct-minus-deranged was +0.02359
[0.00990,0.03761]. Correct-protein direction became statistically positive,
but both registered +0.03 margins were missed.

**Decision.** `P1R1_FAIL_WITH_PARTIAL_CORRECT_PROTEIN_SIGNAL`. Do not tune the
readout post hoc and do not start P1R2, scale-up, or P2-P4. Evidence:
`report/mechanism_refactor/p1r0_signal_localization_seed17_v1/` and
`report/mechanism_refactor/p1r1_mif_residual_seed17_v1/`.

### F-63: Non-additive MIF is protein-specific but not affinity-incremental

**Registration.** P1R2A froze the P1B checkpoint and all P1R1 affinity
controls. It replaced only the representation with six fixed blocks of 47
chemistry-composition coordinates plus one interaction magnitude, retaining
288 features and 289 Ridge parameters. A complete source `220 x 68` panel
defined ligand/protein main effects and the two-way interaction residual. The
label-free pre-gate required source and metaval relative L2 at least 0.01.

**Audit.** Source candidate variance decomposed into 77.576% ligand main
effect, 20.757% protein main effect, and 1.667% interaction. Interaction
residual correct-vs-deranged relative L2 was 1.35675 source and 1.37465
metaval, so the pre-gate passed before affinity values were accessed.

**Result.** Target-macro CI was 0.71110 ligand-only, 0.70513 correct, 0.68140
deranged, and 0.68502 shuffled. Correct-minus-ligand was -0.00596
[-0.01661,0.00387]; correct-minus-deranged was +0.02373
[0.00828,0.03947]. The non-additive statistic preserves significant
correct-protein contrast but supplies no incremental affinity gain.

**Decision.** `P1R2A_FAIL_INTERACTION_SIGNAL_NOT_AFFINITY_INCREMENTAL`. Stop.
P1R2B typed-interaction retraining requires separate authorization; P2-P4 and
production integration remain frozen. Evidence is in
`report/mechanism_refactor/p1r2a_factorization_seed17_v1/` and
`report/mechanism_refactor/p1r2a_gate_seed17_v1/`.

### F-64: Nonlinear MIF probes do not recover affinity energetics

**Registration.** P1R2B0 returned to pair-local P1R1 MIF and excluded the
complete-grid P1R2A residual. Capacity-matched Ridge, quadratic spline and
two-layer MLP families used fold-local PCA32. Hyperparameters were selected by
5x3 nested source homology-component CV; each frozen family scored metaval once.

**Source OOF.** The 220 source targets formed 26 `<40%`-closed components.
Correct-minus-ligand CI was -0.00263 Ridge, -0.01980 spline, and -0.03026 MLP.
Correct-minus-deranged was -0.00073, +0.00100, and +0.00375. No family showed
source-OOF incremental affinity information.

**Metaval.** Correct-minus-ligand was +0.00441 [-0.00246,0.01130] Ridge,
-0.01019 [-0.02211,0.00421] spline, and -0.00797 [-0.01658,0.00076] MLP.
Correct-minus-deranged was +0.00532, +0.00788, and +0.01662. Nonlinear models
increased wrong-protein penalties without improving the correct ligand baseline.

**Decision.** `CURRENT_MIF_AFFINITY_SEMANTICS_NOT_IDENTIFIED`. The evidence
does not support Ridge functional form as the bottleneck. Retire the PLIP-only
P1R2B definition. `P1R2B_AFFINITY_CALIBRATED_MECHANISTIC_ENERGY` is a future
candidate requiring explicit authorization; it was not implemented. Evidence:
`report/mechanism_refactor/p1r2b0_nonlinear_seed17_v1/`.

### F-65: Reversal-rich pairs retain compatibility, not affinity direction

**P1R2B1 registration.** Source labels froze the top 456 of 2,278 ligand pairs
by ordering entropy before metaval evaluation. Frozen P1R2B0 Ridge/spline/MLP
selections were reproduced. Source uncertainty used homology components and
metaval uncertainty used targets. The unchanged +0.03 dual-contrast rule held.

**B1 result.** The strongest MLP metaval contrast was correct-minus-deranged
+0.03895 [0.01791,0.06010], but correct-minus-ligand was only +0.01915
[-0.00296,0.04096]. Source OOF was +0.01275 [0.00111,0.06707] versus deranged
and -0.00199 [-0.07695,0.00494] versus ligand. Source/metaval selected-pair
label ties were 61.72%/54.67%. Verdict:
`PAIR_COMPATIBILITY_WITHOUT_AFFINITY_DIRECTION`; no E0 or P2 authorization.

**F0 result.** Local ChEMBL 37 manifests report 158 source tasks and
`130/42/10` tasks at 20/32/50 exact compounds, but normalized measurement rows
are not local and DAVIS protected-target 40% exclusion is undocumented.
BindingDB is absent; local Papyrus lacks a license manifest and required
sequence/unit/context fields. `AFFINITY_MECHANISM_PILOT_READY=false`.

**Decision.** Stop after B1/F0. Do not reinterpret B0 as ruling out all
nonlinear biological functions: it rules out only registered low-capacity
readouts of compressed global MIF. E0 and typed-interaction T require new
authorization and the F0 data blockers must be resolved first. Evidence:
`report/mechanism_refactor/p1r2b1_rank_reversal_seed17_v1/` and
`report/mechanism_refactor/p1r2b_f0_source_affinity_census_v1/`.

**Verification.** The complete local suite passed in `drug` with
`185 passed in 55.88s`; both formal artifacts record zero recipient-label
reads and no training authorization.

### F-66: Live ChEMBL API cannot exactly rehydrate the frozen F0 row corpus

**FACT.** F0R fixed the 200 accepted assay paths, byte counts and SHA-256 values
from the prior Phase Z artifact manifest. Five assay files matched exactly.
`CHEMBL1000499` returned 28 current rows but 86,668 bytes/SHA `7782b986...`,
against frozen 86,667 bytes/SHA `33beb2...`. The registered exact-hash Gate
stopped before row-level census. Recipient labels were not read; training was
not authorized.

**INFERENCE.** Live-API rehydration is not sufficiently immutable to replace
the missing accepted row files, even though the service still identifies the
release as ChEMBL37. The manifest-only task counts cannot be promoted to a
training corpus.

**UNTESTED HYPOTHESIS.** A pinned static ChEMBL37 database dump may reconstruct
the old payload. This was not downloaded because F0R authorized only the cheap
exact-API recovery. Do not reopen E0/T/P2 or claim DAVIS homology closure from
this failed run.

**Verification.** F0R report SHA-256 is
`fe4d8cfc28cb3eb2274581b009753520e8388f01360cf1fa6bdac53178dadc05`.
All five retained partial files still match the frozen manifest. The complete
`drug` suite passed `190 passed in 72.48s`.

### F-67: Historical API reproduction is separated from scientific corpus reproducibility

**FACT.** F0R remains `HISTORICAL_REPRODUCIBILITY_FAILURE`; it is closed and
is not rewritten as PASS. The live-API rehydrator is historical-only and is
prohibited from constructing training data. A separate
`P1R2B-D0_RELEASE_PINNED_AFFINITY_CORPUS` stage is authorized, beginning with
the official ChEMBL37 SQLite archive and registered SHA-256
`33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`.

**INFERENCE.** The relevant scientific reproducibility object is an immutable
release plus frozen extraction SQL, deterministic normalization and frozen
governance, ending in hashed canonical rows/tasks/splits. Historical API JSON
serialization is not a valid permanent blocker for that object.

**UNTESTED HYPOTHESIS.** Whole-release, binding-assay, single-protein, exact
Ki/Kd filtering will produce enough independent tasks after DAVIS >=40%
protected homology and document-union closure. D0/D1 must test this before E0;
IC50, BindingDB replication, typed interactions and all model training remain
outside the active D0-C acquisition.

### F-68: Release-pinned ChEMBL37 corpus and closure governance pass

**FACT.** The official 5,764,252,857-byte ChEMBL37 SQLite archive matched
SHA-256 `33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281`.
EnergyPilot.v1 contains 343,562 independently verified Ki/Kd rows. Of 41,619
tasks, 4,092 meet the 20-compound/12-non-tied-comparison contract. D1 retains
3,817 tasks, 697 targets and 253 homology-union-document closure components
after excluding 131 E0-Core targets homologous to protected DAVIS targets.
Five folds contain `1467/588/588/587/587` tasks and have zero homology,
document or union-component straddling. Recipient labels were unread.

**INFERENCE.** D0-C and D1 pass scientific reproducibility and governance. The
largest component contains 38.43% of tasks, so future E0 inference must use
closure components and preserve the fixed fold imbalance; 3,817 tasks are not
3,817 independent units.

**UNTESTED HYPOTHESIS.** A local-before-pooling energy map over frozen P1B
geometry improves correct-protein affinity beyond the OOF ligand prior. No E0
training was authorized or executed. Typed interactions and P2-P4 remain
frozen. Evidence: `report/mechanism_refactor/p1r2b_d0_chembl37_v1/`.

**Verification.** The complete post-D0/D1 `drug` suite passed
`200 passed in 58.88s`.

### F-69: E0 local mechanistic affinity potential is authorized

**FACT.** On 2026-08-06 the user explicitly authorized
`P1R2B-E0_LOCAL_MECHANISTIC_AFFINITY_POTENTIAL` after D0-C/D1 passed. The stage
freezes the P1B frontend and applies a small nonlinear map to local
atom-residue chemistry plus frozen contact/distance geometry before pooling.
Global ligand/protein embeddings and target/assay identifiers are prohibited
MAP inputs. Typed interactions and P2-P4 remain frozen.

**INFERENCE.** This isolates whether identified partner geometry lacks only an
affinity-directed local interpretation. It does not claim a thermodynamic free
energy because the source assays do not share complete physical conditions.

**UNTESTED HYPOTHESIS.** A synthetic trainability control will pass, followed by
both preregistered `+0.03` component-macro source OOF contrasts with positive
closure-bootstrap lower bounds. DAVIS metaval remains inaccessible unless that
source Gate passes.

### F-70: E0 input contract remains feasible after label-blind filtering

**FACT.** An affinity-field-discarding scan found 154,165 governed rows, of
which 152,934 satisfy the existing 128-atom and standard-residue P1B contract.
The remaining label-blind floor is 3,783 tasks, 681 proteins, 245 closure
components and 573 tasks in the smallest fold. This passes the frozen D0/D1
floor. No training or recipient-label access occurred.

**INFERENCE.** Full InChIKey, rather than connectivity key, must index ligand
state caches because stereoisomers can share connectivity. The 93,879 ligand
states should be cached after one frozen GINE pass in variable-length form.

### F-71: E0 stops at the synthetic partner-specificity pre-gate

**FACT.** Exact frozen P1B/ESM local caches were completed for 681 retained
proteins and 93,761 retained ligand states. The registered synthetic holdout
gave CI `0.48454/0.68553/0.64934` for ligand/correct/deranged. The increment
over ligand was `+0.20099`, but correct-minus-deranged was only `+0.03618` and
correct CI was below `0.80`. Permutation error was zero. The synthetic Gate
failed. No real source affinity or DAVIS label was read.

**INFERENCE.** MAP can learn a synthetic ligand-residual component, but the
registered control has not shown sufficiently partner-specific held-out
recovery. This blocks biological interpretation of any real-label failure.

**DECISION.** Fail closed before E0-S. Do not tune and retry, read ChEMBL
affinity, access DAVIS, or authorize typed interactions/P2 without a new
explicit decision.

**Verification.** The complete post-E0 `drug` suite passed
`205 passed in 55.43s`.

### F-72: E0S failure localization is authorized without repair

**FACT.** The E0 synthetic failure remains fail-closed. E0S may only replay
the immutable teacher, frozen P1B geometry, MAP checkpoint, selection and
label-blind manifests. Real source affinity, DAVIS, model updates, new seeds,
hyperparameter retries, typed interactions and repair implementation are
prohibited.

**UNTESTED HYPOTHESIS.** Oracle attainability, T0-to-T3 retention,
derangement semantics, final-checkpoint diagnostics, sample diversity and the
197-row provenance transition can distinguish a mis-specified synthetic
control, insufficient frozen geometry, a MAP realization/optimization defect,
or a valid synthetic-identifiability pass.

### F-73: E0S localizes the synthetic loss to T2 -> T3

**FACT.** The frozen teacher oracle attains holdout correct/deranged CI
`1.00000/0.51250`, partner delta `+0.48750`, and zero ties across 1,520 label
pairs. Its exact `8 x 6 x 5` sufficient statistic reconstructs from frozen P1B
geometry within `2.19e-7`. The MAP checkpoint reaches train correct CI
`0.85313`, but holdout correct/deranged CI remains `0.68553/0.64934`.

**INFERENCE.** The synthetic control and teacher-relevant frozen geometry are
not the failing boundary for this constructed target. The observable loss is
T2 -> T3. Existing final-checkpoint artifacts cannot distinguish hypothesis
class realization from optimization because no training or gradient trace was
persisted.

**CAVEATS.** The lexical holdout covers 8/573 fold-4 tasks and 8/71 fold-4
closure components. One of eight wrong-protein controls has local identity
`0.54545`; this violates `<40%`, but the affected task retains strong oracle
and MAP contrasts and does not explain the aggregate gap. The 197-row corpus
transition is fully reconstructed, and activity-ID sets match exactly.

**DECISION.** `MAP_REALIZATION_OR_OPTIMIZATION_DEFECT`. Keep the original E0
synthetic pre-gate fail-closed. No model repair, real E0-S, typed interaction,
DAVIS or P2-P4 work is authorized by this audit.

**Verification.** The required artifacts are hash-bound by
`p1r2b_e0s_evidence_v1/ARTIFACT_MANIFEST.json`; the complete `drug` suite
passed `207 passed in 54.46s`.

### F-74: Typed-tensor synthetic identifiability repair is authorized

**AUTHORIZATION.** Register `P1R2B-E0R0_TYPED_TENSOR_IDENTIFIABILITY` as a
synthetic-only repair. Compare an analytic centered 240D witness, a learned
240-parameter full energy tensor, a learned 114-parameter CP-rank-6 tensor and
the frozen generic MAP under the original split, optimizer schedule and Gate.

**STOP RULE.** No real affinity, DAVIS, PLIP/typed-interaction training,
production integration, CSMO/Band change or P2-P4 work is authorized. Each
learned tensor head passes or fails independently; no post-hoc model selection
or hyperparameter retry is allowed.

### F-75: Typed tensor is sufficient but not identified by the frozen training protocol

**FACT.** The analytic full tensor passes at correct/deranged CI
`1.00000/0.51250`; raw/residual reconstruction error is `2.00e-8/4.51e-7`.
The analytic CP-rank-6 tensor error is `2.34e-13`. Learned full-240 and
CP-rank-6 reach correct CI `0.69013/0.66776` and partner delta
`+0.05987/+0.06711`; both fail the original Gate.

**INFERENCE.** Biologically typed low-rank coordinates solve the representation
ceiling but not identification under the frozen E0 schedule. Both curves remain
improving with nonzero final gradients, and their tensor cosine to the teacher
is only `0.21025/0.22827`; the run does not support convergence or a learned
typed-tensor success claim.

**DECISION.** Stop. No post-hoc epoch, learning-rate, loss or initializer retry;
no real affinity, PLIP/T, DAVIS or downstream authorization.

**Verification.** E0R0 artifacts are hash-bound by `ARTIFACT_MANIFEST.json`;
the complete `drug` suite passed `210 passed in 53.72s`.

### F-76: Three-agent review authorizes E0R1 objective/design/solver audit

**FACT.** Independent objective-math, implementation and governance audits
agree that E0R0's old rank loss uses residual-only scores against total labels,
while final evaluation adds the ligand baseline. Ridge alpha 10 is not an exact
solver for that objective and rank deficiency alone does not imply holdout
prediction nonidentifiability.

**AUTHORIZATION.** Run synthetic-only E0R1 A/B/C and conditional D. C uses raw
centered float64 Moore-Penrose with no alpha. D uses the corrected point residual
plus residual-difference objective and runs only after C train reconstruction
and unchanged-Gate PASS. A new score-blind derangement map must verify `<40%`
local identity and be hash-frozen before scoring.

**BOUNDARY.** No real affinity, DAVIS, CP/rank-1, PLIP/T, production,
CSMO/Band or P2-P4 work is authorized by this registration or its result.

### F-77: E0R1 confirms objective defect and exact design transport

**FACT.** Residual-only versus total-label ordering conflicts on `370/6080`
pairs. Analytic-teacher old-rank gradient is `2.82e-3`, while point and
residual-difference gradients are approximately `3e-10`. The `640 x 240`
design has rank 225 and condition `1.14e8`, but mean holdout row-space coverage
is `0.999816`.

**FACT.** Unregularized Moore-Penrose reconstructs train residuals at RMSE
`3.18e-8` and attains correct/deranged CI `0.99737/0.61447`. This rules out
typed representation and holdout design transport as the E0R0 primary failure.

**PROCEDURAL RESULT.** Corrected full-batch D attains nominal CI checks but
train RMSE is `7.92e-4`, above the frozen `1e-6` precondition after 526 closure
calls. Its Gate is not evaluated. Verdict:
`NOT_RUN_NUMERICAL_PRECONDITION_FAILED`.

**DECISION.** Preserve A/B/C evidence, stop without numerical retry, and keep
real affinity, DAVIS, PLIP/T, production and P2-P4 frozen.

### F-78: E0R1 artifact and regression verification complete

**FACT.** Every SHA-256 entry in the E0R1 artifact manifest matches its local
artifact. The complete `drug` test suite passes `213 passed in 61.45s`.

**DECISION.** Verification does not change the procedural verdict or authorize
another solve. `NOT_RUN_NUMERICAL_PRECONDITION_FAILED` and all downstream
freezes remain in force.

### F-79: Repository evidence consolidation and failure-code removal

**AUDIT.** The frozen theory is implemented directly for the Band polytope,
law class `K(beta)`, simplex-valued coefficient map, positive ridge, `B(z)`
assembly and CSMO. The biological frontend was only interface-compatible: its
28-dimensional `z` used arbitrary bounded latent projections, while the
P1B-passing contact/distance bridge was not connected to the production
pipeline. No experiment admitted a protein-specific affinity statistic to `z`.

**VERDICT.** `MATHEMATICAL_OPERATOR_IMPLEMENTED_AND_CONTRACT_TESTED;
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z; NO_VALIDATED_END_TO_END_DTA_MODEL`.
This is not a failure of the frozen operator. The earliest unresolved boundary
is validated local geometry to signed transferable affinity semantics.

**CODE CONSOLIDATION.** `model/` now contains only mathematical/operator
primitives, the P1B-validated local encoders, and the geometry bridge. The
failed assembled biological pipeline, support/QPMA/meta-state implementation,
legacy trainer/evaluator, P1C/P1R*/F0R scripts, and their dedicated tests were
deleted. `scripts/` now contains only passed data, structure-geometry, and
release-governance workflows.

**RESEARCH CONSOLIDATION.** The unresolved E0 objective/design/solver chain,
its focused tests, proposal context and required artifacts moved to
`research/e0_identifiability/`. All older research implementations were deleted
after their conclusions were consolidated. This relocation does not authorize
a numerical retry.

**REPORT CONSOLIDATION.** Duplicate master reports, Phase X/Y/Z reports,
superseded smoke artifacts, failed P1C/P1R*/F0R artifacts, the invalid-float32
E0R0 run, and two failed local DAVIS checkpoints were deleted. The initial
verified deletion set contained 54 code/report targets totaling 117.61 MiB,
plus 50 superseded top-level research entries and Python caches. Retained report
evidence is limited to the current split protocol, P1A local audit, final P1B
checkpoint/Gate, and D0-C/D1 PASS artifacts. The theory hash ledger moved to
`theory/FINAL_FROZEN_THEORY/THEORY_HASHES.json`.

**AUTHORITY.** `EVIDENCE_CONSOLIDATION_AND_FAILURE_TRIAGE.md` is the current
evidence map and failure-triage entry point. `README.md`, `task.md`,
`experiment.md`, and `AGENT_HANDOFF.md` were rewritten as concise current-state
documents; this history remains the detailed failure ledger.

**BOUNDARY.** Real E0-S, DAVIS, PLIP/T, production integration, CSMO/Band
changes and P2-P4 remain frozen. Recipient-label reads remain zero.

**VERIFICATION.** The consolidated suite passes `70 passed in 14.76s` in the
`drug` environment. It covers retained production/verified code and the 10
focused E0 research tests. The former 213-test count is historical and included
tests whose only subjects were deleted failed or superseded implementations.

### F-80: Three-agent proposal audit and E0R2 synthetic numerical closure

**AUTHORIZATION.** The user authorized a multi-agent research analysis of the
directional statistical potential proposal and synthetic-only code testing.
Independent mathematics, code-boundary and governance agents agreed that the
proposal must be decomposed. Only E0 numerical closure was executable under the
current label-free boundary; PLIP/type supervision, reference states, real
affinity, few-shot adaptation, DAVIS and production remained frozen.

**CLAIM AUDIT.** The proposal correctly targets geometry-to-signed-affinity as
the unresolved scientific boundary, but several statements were too strong.
Current evidence does not show that P1B has already lost real affinity direction;
the synthetic teacher is recoverable from P1B geometry and chemistry. A
negative log bound/reference ratio is structural log-odds, not free energy, and
its sign convention conflicted with the proposal's positive-favorable wording.
Expected potential under a distance distribution does not preserve the full
uncertainty law. Cross-fitting is leakage-safe residualization, not automatically
a Neyman-orthogonal or causal construction. A few-shot dimension bound `d<=k`
does not establish identification without rank, conditioning and query row-space
coverage.

**PREREGISTRATION.** E0R2 froze the existing 240D design, 32/8 historical split,
E0R1 deranged features, residual/difference target, float64 SVD augmented solve,
`rcond=1e-10`, and the original synthetic CI thresholds. No hyperparameter or
holdout selection was allowed.

**RESULT.** The deterministic solve reached train RMSE `3.188e-8`, maximum
absolute error `2.007e-7`, corrected objective `1.567e-15`, and full-gradient L2
`6.150e-17`. Historical correct/deranged CI was `0.99737/0.61447`, with
correct-minus-ligand `+0.51283` and correct-minus-deranged `+0.38289`. Verdict:
`SYNTHETIC_OBJECTIVE_DESIGN_SOLVER_IDENTIFIED`.

**LIMITATION.** The eight holdout tasks are reused development diagnostics and
the inherited derangement map concentrates wrong-protein reuse. The result
closes synthetic objective/design/solver consistency only. It does not establish
directionality, structural log-odds, affinity energetics, transfer or biological
statistic admission.

**CODE DISPOSITION.** The new contract, preregistration, runner, tests, viewpoint
audit and artifacts remain in `research/e0_identifiability/`. Nothing was moved
to `model/` or normal `scripts/` because a synthetic PASS is below the repository
admission threshold.

**VERIFICATION.** The complete retained suite passes `70 passed in 13.88s` in
the `drug` environment. E0R2 artifact hashes and JSON documents validate, and
`model` exports neither the deleted assembled DTA pipeline nor the research MAP.

### F-81: T-DIR-P0 lightweight structural learnability pilot is negative

**AUTHORIZATION.** The user authorized a lightweight training test after the
multi-agent directional-potential analysis. A research-only pilot was frozen
before PLIP annotation: 24/8/8 train/validation/test complexes, 40 globally
distinct homology groups, zero selected cross-split exact-connectivity or
Murcko-scaffold overlap, PLIP 3.0.1 weak labels, frozen P1B features, and fixed
logistic probes. Real ChEMBL affinity, DAVIS and recipient labels were not read.

**DATA RESULT.** All `40/40` selected structures completed the executed direct
event annotation pipeline, producing 17,556 oracle-near atom-residue candidate
pairs. The primary hydrophobic channel had `53/26/24` positives from `15/7/6`
train/validation/test complexes. Its PLIP event-to-canonical mapping coverage
was `144/147 = 97.96%`.

**MODEL RESULT.** Hydrophobic D0/D1/D2 test AP was
`0.00987/0.03735/0.01996` at prevalence `0.00620`; validation AP was
`0.01906/0.04445/0.01699`. D2 train AP was `0.25496`, showing a large
train-to-held-out collapse. D2 failed the frozen test AP-lift `>=0.10`
condition and did not exceed D0 on validation. Verdict:
`PILOT_LEARNABILITY_SIGNAL_NOT_OBSERVED`.

**PROCEDURAL AUDIT.** A namedtuple traversal defect made all multi-atom
salt/pi/cation-pi mappings uninterpretable. The mapper and a regression test
were fixed after the run, but the used test panel was not rerun. The executed
runner also omitted preregistered shuffle controls and a full ligand chemistry
round-trip audit. These deviations cannot reverse the already negative primary
hydrophobic result, but they prohibit any complete annotation-Gate claim. The
executed and corrected runner hashes are bound in `postrun_audit.json`.

**DECISION.** Keep the pilot code, preregistration, post-run audit and negative
report in `research/e0_identifiability/`. Do not promote code to `model/` or
normal `scripts/`. Full T-DIR, real affinity, DAVIS, production `z`, CSMO/Band
changes and P2-P4 remain frozen. A retry requires a new sealed panel and a
separately registered formal stage; the 8 test complexes are now development
evidence.

**VERIFICATION.** The complete retained suite, including focused namedtuple
mapping regression coverage, passes `70 passed in 13.40s` in `drug`.

### F-82: T-BASIS-R0 fixed radial basis partner recoverability passes

**ROUTE CHANGE.** The deployment contract remains protein sequence plus ligand
2D graph; holo structure is privileged training information only. Sparse PLIP
event taxonomy was downgraded to an auxiliary semantic audit. The primary
hypothesis became a fixed continuous chemogeometric basis, followed only later
by dataset-aware delta-affinity and an identified few-shot section.

**PREREGISTRATION.** T-BASIS-R0 froze a fresh `192/64/64` structure panel,
excluded all 40 T-DIR-P0 records, required 320 distinct homology groups/PDBs/
sequences, P1B-held-out validation/test, complete sequence mapping, and a
one-to-one score-blind wrong-protein map with identity `<0.40` and reuse zero.
No affinity, DAVIS or recipient labels were authorized.

**BASIS.** The privileged teacher is a fixed permutation-invariant
`8 ligand chemistry x 6 residue chemistry x 6 Gaussian radial = 288D` tensor.
The student replaces exact holo distances with fixed radial expectations under
the frozen P1B five-bin distance distribution. Only a shared train-only `6x6`
Ridge radial calibration (`alpha=1e-3`) is learned.

**RESULT.** All 320 complexes completed. On test, train-mean/correct/deranged
standardized MSE was `1.1001/0.5157/0.6875`. Reconstruction gain was
`+0.5312 [0.4433,0.5962]`; wrong-partner degradation was
`+0.1561 [0.1070,0.2007]`. Validation gains were `+0.5324` and `+0.1455`.
All six preregistered conditions passed. Verdict:
`RADIAL_BASIS_PARTNER_RECOVERABILITY_IDENTIFIED`.

**LIMITATION.** The partner contrast changes both residue composition and P1B
distance predictions, so it does not isolate pure pair geometry from protein
marginal chemistry. The result covers two-body radial moments only. It is not
angular/many-body recovery, affinity energetics, universality, transfer,
few-shot identification or biological `z` admission.

**DECISION.** Retain T-BASIS-R0 in research. It authorizes only a separately
registered structure-only angular/many-body privileged-basis study with
marginal-preserving controls. E-AFF, cross-dataset replication, RFSA, DAVIS,
production integration, CSMO/Band changes and P2-P4 remain frozen.

**VERIFICATION.** Generated artifact hashes validate and the complete retained
suite passes `70 passed in 15.15s` in the `drug` environment.

### F-83: E-AFF-P0 finds no population-shared radial affinity direction

**QUESTION.** After T-BASIS-R0 established a recoverable 288D radial
chemogeometric statistic, E-AFF-P0 asked whether one population-shared linear
direction adds within-task Ki/Kd ranking information beyond a closure-OOF
ligand prior. It did not test target-specific coefficients or a richer basis.

**DESIGN.** The score-blind sample contained one task from each of 245 closure
components and 20 ligand states per task (`4,900` rows). The shared direction
was trained only on correct proteins with task-balanced residual differences.
Deranged and marginal-preserving coupling-null arms were evaluation-only and
used the identical direction. The wrong-protein map was one-to-one, had reuse
zero, and maximum exact local identity `0.39394`.

**RESULT.** Component-macro ligand/correct/deranged/null CI was
`0.55225/0.54209/0.54445/0.54210`. Correct-minus-ligand was
`-0.01016 [-0.02069,-0.00018]`, correct-minus-deranged was
`-0.00236 [-0.01301,0.00805]`, and correct-minus-null was
`-0.00001 [-0.00558,0.00530]`. All frozen feasibility conditions failed.
Verdict: `SHARED_DIRECTION_NOT_OBSERVED_H0_DATA_SUPPORTED`.

**INTERPRETATION.** This rejects only the existence of the tested global shared
direction. It does not show that the fixed basis contains no affinity
information, so a score-blind task-local headroom diagnostic was authorized.

**AUDIT.** The independent post-run audit reproduced all component metrics and
bootstrap intervals, verified tensor shape and finiteness, coupling-null
marginals, exact hashes, derangement constraints, and zero DAVIS/recipient
reads.

### F-84: E-AFF-H0A finds task-local headroom without partner specificity

**QUESTION.** H0A tested whether the fixed 288D basis has task-local predictive
headroom even though one shared population direction failed. It was an oracle
headroom diagnostic, not a production model or few-shot section.

**DESIGN.** One deep task was selected score-blind from each of 107 closure
components. Each task used 20 fixed fit ligands and 20 untouched test ligands
(`4,280` rows total). An independent fixed-alpha Ridge direction was fitted on
correct-protein residual differences for each task. Correct, deranged and
coupling-null test arms shared exactly that direction. No wrong-protein example
entered training.

**RESULT.** Component-macro ligand/correct/deranged/null CI was
`0.55404/0.64226/0.63362/0.63817`. Correct-minus-ligand was
`+0.08821 [0.06761,0.10998]`, demonstrating substantial held-out task-local
headroom. Correct-minus-deranged was only
`+0.00864 [0.00338,0.01462]`, below the frozen `+0.03` partner margin.
Correct-minus-null was `+0.00408 [0.00061,0.00747]`. Both Ki and Kd had
positive task-local headroom, but Kd used only 15 components and its
coupling-null contrast was negative. Verdict:
`TASK_LOCAL_RADIAL_HEADROOM_WITHOUT_PARTNER_SPECIFICITY`.

**FAILURE REASON.** The basis is not devoid of affinity-predictive information;
it ranks unseen ligands substantially better after task-local calibration.
However, replacing the protein preserves nearly all of that improvement. The
observed headroom is therefore compatible with ligand chemistry, assay-series
structure, or protein marginal effects and does not identify a
correct-partner-specific affinity section.

**DECISION.** H0-B cross-assay target transport and RFSA are not authorized.
Angular expansion is also not automatically authorized. Any new stage must
first isolate the partner-conditioned component from ligand/series and protein
marginals. Nothing from E-AFF was promoted to `model/` or normal `scripts/`.

**AUDIT.** Independent reconstruction passed all checks: `107` selected tasks,
`4,280` rows, exact 20/20 partitions, one-to-one `<0.40` derangement, coupling-
null marginal error `1.81e-7`, exact task/component/bootstrap metrics, artifact
hashes, and zero DAVIS/recipient reads.

### F-85: E-AFF-H0C shortcut removal does not recover partner affinity

**QUESTION.** H0C tested whether H0A's weak positive partner/coupling contrast
was hidden by a stronger task-local ligand shortcut. It removed both fixed-
tensor marginals and matched the 20-shot information available to a local
ligand nuisance before fitting the interaction residual.

**PREREGISTRATION.** Every H0A task was excluded. Selection used only governed
label-blind rows and ligand structures. One new task per closure component was
required to support 20 support and 20 test ligands with strict within-task
Murcko scaffold separation. P1B, T-BASIS, Ridge `alpha=10`, closure-OOF global
ligand prior and derangement rules remained frozen.

The calibrated tensor contains negative coordinates, so the registered object
was explicitly algebraic rather than probabilistic:
`psi=(phi-phi_null)/total`. The frozen 128D pooled ligand state produced a
five-fold support-cross-fitted nuisance residual. The interaction direction was
trained only on correct-protein residual differences; wrong proteins were
evaluation-only.

**DATA.** The new panel contained 54 tasks from 54 closure components and 2,160
rows. All H0A task overlap and all support/test scaffold overlap were zero.
Derangement was one-to-one with reuse zero and identity `<0.40`.

**RESULT.** Component-macro Global-L/Local-L/correct/deranged CI was
`0.55487/0.59635/0.59244/0.59092`. Local ligand adaptation added `+0.04147`,
but correct-minus-Local-L was `-0.00391 [-0.02040,0.01191]` and
correct-minus-deranged was `+0.00152 [-0.01024,0.01424]`. Neither frozen
condition passed. Verdict: `FIXED_RADIAL_INTERACTION_RESIDUAL_NOT_OBSERVED`.

**FAILURE REASON.** H0A's task-local headroom is reproducible as ligand/series
SAR, but it is not recovered from the pure radial chemistry-distance coupling
after nuisance removal. The tested fixed radial interaction therefore lacks
identified incremental affinity value under the 20-shot scaffold-cold
estimand. This result does not yet say whether real source protein-by-ligand
interaction is too weak/sparse or whether radial coordinates omit required
biology.

**DECISION.** H0-B and RFSA remain unauthorized. Orientation is not authorized
by this failure alone. The next possible stage is a label-light X0 census of
same-document target-by-ligand rectangles, followed only if supported by a
separately registered double-difference Gate.

**AUDIT.** Independent post-run reconstruction passed every check and reproduced
all metrics and intervals exactly. Double-centered marginal error was
`2.29e-17`; DAVIS and recipient reads were zero. The first command session was
orphaned after writing only selection/derangement files. It produced no feature
or metric result and was later removed as disposable runtime residue under F-87.

### F-86: E-AFF-X0 closes ChEMBL crossed interaction as underdetermined

**QUESTION.** X0 asked whether the governed ChEMBL37 Ki/Kd source has enough
independent `2 targets x 2 ligands` rectangles to identify real affinity
interaction before judging the radial basis or opening any value labels.

**LABEL FIREWALL.** Selection started from 152,737 governed activity IDs in
`rows.label_blind.jsonl`. A dedicated SQL projection queried only pinned
ChEMBL37 document, endpoint, assay and target-independent structured context
metadata. It selected no activity value, pAffinity, published value or pChEMBL
field. DAVIS and recipient reads were zero.

**DEPENDENCY CONTRACT.** Exact Ki/Kd endpoints remained separate. Variant
assays were excluded. Rectangles sharing a structured panel were dependent;
panels touching the same D1 homology-document closure component were unioned.
The frozen lower-bound requirement was 245 effective components per endpoint,
derived for an ideal variance test targeting interaction RMS / assay noise
`0.5`, one-sided alpha `0.05`, power `0.80`.

**RESULT.** Ki contained 597 eligible panels, 53,673 cells, 1,794 target pairs
and 1,059,169 nominal rectangles. Dependency closure left only 36 components;
18 were replicate-supported. Kd contained 34 panels, 4,995 cells, 534 target
pairs and 232,875 nominal rectangles, but only 12 components; four were
replicate-supported. Largest-component rectangle fractions were `0.56494` and
`0.76757`. Verdict: `STOP_SOURCE_INTERACTION_UNDERDETERMINED`.

**FAILURE REASON.** The apparent million-scale sample size is
pseudoreplication from shared cells, panels, documents and homologous target
families. The source lacks enough independent crossed units to distinguish
protein-by-ligand interaction variance from assay/panel noise at the registered
sensitivity. This does not assert that biological interaction is absent.

**DECISION.** X1 and X2 are not authorized; no double-difference labels were
read. H0C cannot be reclassified as representation failure, and orientation,
RFSA, DAVIS and production integration remain frozen. Continuing requires a
separately registered label-blind census of a genuinely crossed source or
selectivity corpus; thresholds and closure cannot be weakened post hoc.

**AUDIT.** Independent reconstruction exactly reproduced every panel,
rectangle, dependency and endpoint count, verified existing hashes, found no
affinity fields in SQL/artifacts, and confirmed zero DAVIS/recipient reads.

### F-87: Repository evidence consolidation after H0C/X0

**SCOPE.** The repository was re-audited after the terminal H0C and X0 results.
The cleanup followed an evidence contract: `model/` retains only verified
mathematical/P1B components; `scripts/` retains only passed data, geometry and
governance workflows; research stages retain their preregistration, terminal
result, independent audit and final artifact regardless of whether the result
was positive or negative.

**REMOVED.** Generated Python and pytest caches, two pytest temporary trees, the
E-AFF-P0 launch that stopped after writing only a selection file because of an
import failure, and the first H0C launch that stopped after selection and
derangement were removed. These objects contained no unique feature matrix,
trained weight, metric, verdict or provenance needed to reproduce the final
stages. Their procedural outcomes were already recorded in the corresponding
terminal reports and this ledger.

**RETAINED.** Frozen theory, raw/release data, P1B and D0/D1 PASS evidence, every
final research artifact from E0 through X0, all preregistrations/result reports,
and all independent audit code/tests were retained. Negative scientific results
were not deleted. No production or research model semantics changed.

**STATE.** `AGENT_HANDOFF.md`, the root triage, `task.md`, `experiment.md` and
the research index now agree that H0C is a source negative result, X0 is a data
stop, X1/X2 and all downstream biological integration remain unauthorized, and
there is no automatic next stage.

**VERIFICATION.** The complete `drug`-environment regression passed after the
consolidation: `70 passed in 11.15s`. A final scan found no research imports in
`model/` or normal `scripts/`, no references to removed runtime artifacts, and
no remaining Python/pytest cache directories after post-test cleanup.

### F-88: E-AFF-R0 readout scope audit

**SCOPE.** Registered as a diagnosis of the affinity readout, not of biology and
not as a Gate. No affinity label, DAVIS label or recipient label was read. Inputs
were the repository's own `metrics.concordance`, the published H0C per-task
metrics, and a simulated generative model.

**RESULT.** Within-task concordance changed by `0.0` under per-task prediction
shift, per-task prediction rescale, per-task label shift and per-task label
rescale — exact invariance, not approximate. A simulated predictor holding a
task's affinity level perfectly and nothing else scored exactly `0.5000` at
level-variance shares `0.059/0.200/0.500/0.800/0.941/0.985`, while its
location-sensitive RMSE advantage over a global mean grew from `1.033/1.001` to
`8.832/1.017`. H0C's published per-task contrasts were correct-minus-local
`-0.00391`, deranged-minus-local `-0.00543` and correct-minus-deranged
`+0.00152`, with the geometry term changing `51/54` tasks. Verdict:
`READOUT_BLIND_TO_TASK_LEVEL_AFFINITY_LOCATION|PERFECT_LEVEL_PREDICTOR_SCORES_CHANCE_AT_EVERY_VARIANCE_SHARE`.

**ADDITIONAL DESIGN FACT.** H0C removed the same channel upstream: the geometry
received `y - global_ligand_prior - task_local_ligand_nuisance`, where the
nuisance is fitted on 20 labelled supports of the correct protein's own task and
the resulting `local_score` is added to both the correct and the deranged arm.
Both arms therefore held the correct protein's task level before the contrast.

**DECISION.** No historical verdict is reclassified. P1C, P1R1, P1R2A, P1R2B0,
P1R2B1, E-AFF-P0, H0A and H0C remain valid negative results about within-task
ranking information. Their scope is annotated, not their outcome. R0 does not
show that protein-specific affinity lives in the location channel; it shows that
if it does, the evidence chain could not have detected it. R0 authorizes nothing.

### F-89: E-AFF-X0-FEAS unit feasibility audit

**SCOPE.** Registered as an audit of the X0 estimand before any corpus
acquisition was funded against X0's named continuation. Label-blind: no SQLite
connection, zero affinity fields, zero DAVIS or recipient reads, D0 task-manifest
fields admitted only through an explicit whitelist.

**RESULT.** A rectangle requires two proteins inside one document-keyed panel,
and D1 closure unions every pair of targets sharing a document, so both proteins
of every rectangle already lie in one closure component. Predicted panels
touching more than one closure component `0`; observed `0`. The closure-component
universe of the governed corpus is `245` against a frozen requirement of `245`
per endpoint; only `202` components carry Ki rows and `72` carry Kd rows.
Panel-free rectangle-capable ceilings were Ki `57`, Kd `12`. Recomputing the
closure over the full governed D0 corpus (`37,783` tasks, `4,787` proteins,
`459` components) and over shallower populations raised the best ceiling only to
Ki `97` and Kd `56`; the Ki ceiling *falls* from `97` to `72` as the corpus
grows because added documents merge components. The frozen `245` was
independently re-derived from its stated one-sided chi-square design
(`ratio 1.25`, `alpha 0.05`, `power 0.80`, `df n-1`). Verdict:
`X0_UNIT_REQUIREMENT_UNATTAINABLE_BY_CONSTRUCTION`.

**DECISION.** X0's verdict `STOP_SOURCE_INTERACTION_UNDERDETERMINED`, its
metrics and its artifacts are retained unchanged and are recorded as a
specification-induced stop. Acquiring a more genuinely crossed
source/selectivity corpus is withdrawn as the named continuation under that
unit, because crossing and document-disjointness are produced by opposite kinds
of study. X1 and X2 remain unauthorized.

### F-90: E-AFF-X0-B crossed design re-registration

**SCOPE.** Registered to replace the X0 independence unit and its effective
sample-size calculation only. Unchanged: interaction RMS over assay noise `0.5`,
variance ratio `1.25`, alpha `0.05`, power `0.80`, required effective `n = 245`,
and the `+0.03` affinity Gate margins. Label-blind, running only on the panel
geometry X0 already published.

**RESULT.** Under the model `y[p,t,l] = mu[p] + alpha[p,t] + beta[p,l] +
delta[t,l] + eps`, the double difference cancels every effect additive in target
or in ligand, so independence of `DD` requires cell-disjointness rather than
target/ligand disjointness. Ki packed `11,168` cell-disjoint units across `36`
clusters from `205` distinct target pairs, `224` targets and `19,062` ligands;
Kd packed `1,041` across `12` clusters from `49` pairs, `73` targets and `1,256`
ligands. The conservative target-and-ligand-disjoint comparator gave `705` and
`62`. With `DEFF = 1 + (m_A - 1)rho` and `n_eff = N/DEFF`, breakeven thresholds
were `rho* = 0.0915` (Ki, cap 32) and `rho* = 0.0164` (Kd, cap 125), against
hard cluster bounds `G/245` of `0.1469` and `0.0490`. At `rho = 1` the optimal
design collapses to one unit per cluster and reproduces X0's `36` and `12`
exactly. Verdict:
`X0B_CONDITIONAL_DESIGN_SUPPORTED_KI|X0B_CONDITIONAL_DESIGN_SUPPORTED_KD`.

**DECISION.** This is conditional design support about achievable effective
sample size. It is **not** evidence that protein-by-ligand affinity interaction
exists. X1 remains unauthorized and, when registered, must first estimate `rho`
with uncertainty from exact-assay replicates and abstain without testing if the
upper confidence bound exceeds `rho*`. The protein axis is thin — `205` and `49`
distinct target pairs — so any future positive would generalise over ligands far
better than over proteins.

**LIMITATION.** Packing counts are auditable greedy lower bounds, not maximum
packings. Per-cluster caps are reported as a design option and no data was
subsampled.

### F-91: E-AFF-L0 protein affinity-location Gate, not run

**SCOPE.** First registered stage addressing Claim A: does the correct protein
provide affinity-location information beyond population, ligand-only and
protein-sequence-only baselines? Distinct estimand from X0-B/X1, which address
Claim B. Research-only; nothing entered `model/` or normal `scripts/`.

**PRECONDITIONS PASSED.** The operator/anchor contract froze with
`L0_OPERATOR_AND_ANCHOR_CONTRACT_FROZEN`: all seven frozen anchors lie in the
band polytope, the six logistic anchors satisfy stochastic dominance with
maximum gap `0.0`, band mean intervals rise strictly from `[0.186,0.303]` to
`[0.845,0.889]`, mixture monotonicity violation `0.0`, Hausdorff-`W1` stability
violation `0.0`, mesh unchanged at `M=32`, and all `258` frozen theory files
hash-matched. The location estimand check admitted **Ki only**
(`C1=79`, `C2=218`, `C3=0.630`) and excluded **Kd** (`C1=10 < 30`). R0, X0-FEAS
and X0-B independent audits returned `POSTRUN_AUDIT_PASS`.

**ASSAY NOISE.** Estimated before any arm was scored, from `4,261` replicate
cells and `4,840` degrees of freedom: `sigma_assay = 0.47971`
`[0.47034, 0.48946]` log units, giving `margin_L0 = 0.23985`.

**EXECUTED DESIGN.** `115` tasks from `115` distinct closure components, `20`
ligands each, `2,300` observations, five-fold cross-fitting by closure
component, `355` tasks consumed by P0/H0A/H0C excluded before selection. All
five arms shared folds, estimator, bandwidth rule and a seven-dimensional
bounded input. Component-macro band loss was `0.25662/0.23153/0.23099/0.23107/
0.23108` for `A0/A1/A2/A3/A4`; location error in log units was
`1.24758/1.26684/1.26647/1.26346/1.26351`.

**TERMINAL VERDICT.** `L0_NOT_RUN_NUMERICAL_PRECONDITION_FAILED`. Registered
gate condition 3 compares empirical coverage, registered as containment of the
observed step CDF by the emitted band. On the fixed `33`-point mesh a band must
span nearly the whole unit interval to contain a step, so every arm returned
exactly `0.0` and the statistic carried no information; one of three gate
conditions did not execute as specified. The defect is decidable from the mean
interval widths alone and does not depend on arm performance. The run therefore
fails closed to NOT-RUN rather than reporting a scientific negative.

**SUPPORTING DIAGNOSTIC.** Ligand-only did not beat population-only in location
error, gain `-0.01926`, so no positive control was established and a null
protein result would be uninterpretable regardless. This was not a registered
criterion and is recorded as a diagnostic only. It is consistent with H0C, where
the global ligand prior on the full 128-dimensional state reached only `0.55487`
component-macro CI.

**DECISION.** Claim A remains untested;
`PROTEIN_SPECIFIC_AFFINITY_LOCATION_NOT_YET_TESTED` stands. No historical
verdict is overturned or confirmed. The `115`-task panel is **consumed** and
must not be reused as untouched validation. Rerunning with a repaired coverage
statistic on this panel is prohibited; a corrected Gate requires a new
registration and a fresh panel, and must preregister both an informative
coverage statistic and an explicit ligand-only positive-control precondition.
Nothing was admitted to `z`, `model/` promotion remains blocked, and DAVIS,
recipient labels, X1, X2, angular work, RFSA, theory changes and P2-P4 remain
frozen. DAVIS and recipient reads remain zero.

### F-92: E-AFF-L0R corrected affinity-location Gate, not run

**SCOPE.** The corrected repeat that F-91 prescribed, registered before
selection or scoring, inheriting every frozen element of L0 unchanged and fixing
only the two named defects plus the estimator they implicated.

**REPAIRS.** (1) Mean-interval coverage replaced step-containment coverage and is
now informative, varying from `0.0674` for the population arm to `0.1664` for
ligand-only instead of being identically `0.0`. (2) A registered positive-control
precondition was added: ligand-only must beat population-only on location error
by at least `0.1 * sigma_assay` with a positive bootstrap lower bound. (3) The
fixed-bandwidth kernel conditional CDF, which returns nearly the marginal in
seven dimensions, was replaced by a `k`-nearest-neighbour conditional CDF with
`k = ceil(sqrt(n_train)) = 56`, identical for every arm.

**PANEL.** `195` fresh tasks in `78` closure components, at most three tasks per
component, `20` ligands each, `3,900` observations. All `470` tasks consumed by
P0, H0A, H0C and L0 were excluded; no task was reused. Closure components overlap
L0's because only `78` unconsumed components remained, so L0R is development
evidence, not untouched validation.

**RESULT.** Component-macro band loss `0.24210/0.21298/0.22674/0.22121/0.22119`
and location error `1.12972/1.09552/1.18342/1.15148/1.15165` log units for
`A0/A1/A2/A3/A4`. The positive control returned `+0.03421` log units with 95% CI
`[-0.03304, 0.10793]` against a required `0.04797`. Verdict:
`L0R_NOT_RUN_POSITIVE_CONTROL_ABSENT`.

**DECISION.** No protein verdict was computed or reported. The arm contrasts
`A3-A1`, `A3-A2` and `A3-A4` exist in the artifact but are not evidence about
proteins, because they come from a readout that failed its positive control;
recording them as a protein finding is precisely what the precondition prevents.
`PROTEIN_SPECIFIC_AFFINITY_LOCATION_NOT_YET_TESTED` stands.

**USEFUL MEASUREMENT.** With `sigma_assay = 0.47971` log units, the best
available cross-component location signal in governed ChEMBL37 Ki — ligand
identity itself — is worth roughly `7%` of assay noise and is not separable from
zero across `78` closure components. This bounds what any L0-style location Gate
can achieve on this corpus at this scale, and indicates that a further attempt
must change the design power or the estimand rather than the model.

**BOUNDARIES.** Nothing admitted to `z`; `model/` promotion blocked; DAVIS,
recipient labels, X1, X2, angular and many-body work, RFSA, theory changes and
P2-P4 remain frozen. DAVIS and recipient reads remain zero. No rerun on this
panel is permitted.

### F-93: PKIS/F6I identifiability-resolution intake

**SCOPE.** Integrated the supplied
`MetaSieve_identifiability_resolution_20260808.zip`, standalone final report,
`law_bridge.py` and `component_independent.py` as a new research-only evidence
branch. No production `model/`, normal `scripts/`, frozen theory, DAVIS,
recipient or governed ChEMBL affinity path was modified or read.

**PROVENANCE.** Archive SHA-256
`02d5a3586a90caba6cf6392edc0796d24948b06d49c6a389377d509635d842c1`;
107 entries, 90 files, 1,116,376 uncompressed bytes. All archive paths passed a
path-traversal check. The standalone report, law bridge and F6I implementation
were byte-identical to their archive copies. Exact hashes and the scope decision
are recorded in `research/IDENTIFIABILITY_RESOLUTION_INTAKE.json`.

**CODE VERIFICATION.** The imported PKIS and section-operator suite initially
returned `38 passed, 1 failed`. The failure was an undeclared test-only working-
directory dependency on an unbundled KiSSim CSV, not a model failure. The test
was made self-contained with an equivalent temporary 85-position CSV fixture;
runtime research data loading was unchanged. The suite then returned
`39 passed`. PKIS v1/v2 output hashes match their manifests.

**PROVENANCE LIMITATION.** The historical F0 manifest binds
`ceiling_probe.py` SHA-256 `2c925e...`, while the supplied current script hashes
to `6ca01a...`. Therefore the F0 artifact is preserved as historical evidence
but is not claimed byte-identically rerunnable from the supplied code alone.

**EVIDENCE ACCEPTED.** On consumed kinase activity panels, the supplied branch
supports a component decomposition into a protein-dependent, source-atlas
zero-shot surface and a bounded, protein-independent one-dimensional support
location statistic. The latter is support-multiset permutation invariant and
satisfies `d_adapt = 1 <= k = 5`. PKIS2 component contrasts and the isolated
interaction contrast on Anastassiadis are positive as recorded in the supplied
artifacts.

**UNCHANGED TERMINAL VERDICT.** The registered F6I total verdict remains
`F6I_COMPONENTS_NOT_ADMISSIBLE`. On Anastassiadis the nearest-protein raw MSE is
better by approximately `1.47e-5`; the supplied report may interpret this as an
additive-main-effect mixture, but the preregistered total Gate cannot be
rewritten after observation. PKIS2 and Anastassiadis are consumed development
panels, not fresh external validation.

**THEORY BOUNDARY.** The research `law_bridge.py` passes mass, barycentre,
column-stochasticity, tridiagonal and permutation-invariance tests. It produces
a categorical law on a seven-point mesh followed by mean-preserving Markov
diffusion. The frozen production object instead uses a 33-point mesh, CDF-band
polytope, simplex coefficient map and band assembly. The abstract notation is
compatible but the types are not identical. No production operator equivalence
or theorem-certified biological `z` admission is claimed.

**DECISION.** New formal state:
`KINASE_PANEL_COMPONENT_IDENTIFIABILITY_OBSERVED_IN_DEVELOPMENT;
F6I_TOTAL_GATE_NOT_ADMISSIBLE;
FRESH_ENDPOINT_CONSISTENT_EXTERNAL_ADMISSION_NOT_RUN`. The next scientifically
admissible action is a separately registered, once-only evaluation of the
frozen decomposition on a fresh endpoint-consistent, target/ligand/document-
governed panel. A future PASS could authorize a production-interface
registration only; it would not directly authorize overwriting `model/`, the
frozen theory, CSMO/Band, DAVIS or P2-P4.

### F-94: active-tree consolidation and component-algebra promotion

**AUTHORIZATION.** The user authorized new external-dataset acquisition,
promotion of reusable information into `model/`/`scripts/`, deletion of failed
or ineffective research code after historical consolidation, deep report/data
summarization, and creation of a code-only branch.

**RESEARCH DELETION.** Removed all 325 tracked files under the former
`research/` tree (`86,869,561` bytes) plus ignored interpreter caches. The
deleted tree contained terminal P1C/E0/E-AFF implementations, synthetic retry
artifacts, T-DIR/T-BASIS development artifacts, PKIS mechanism pilots, section
operator pilots, duplicate preregistrations and consumed-panel reports. The
scientific metrics, verdicts, provenance limitations and stop rules were already
recorded in F-1 through F-93 and the evidence-triage document. The exact tree is
recoverable from Git commit `8b7789e`; deletion is a working-tree
consolidation, not destruction of provenance.

**PROMOTED CORE.** Retained only the algebra established by the F6I audit:

```text
prediction(P,L;S) = biological_surface(P,L) + location(S)
```

`location(S)` is one-dimensional, bounded, support-multiset permutation
invariant and independent of the protein representation. It is implemented in
`model/component_statistic.py` but deliberately not exported by `model` and not
connected to biological `z`. `scripts/evaluate_component_statistic.py` provides
a query-label-free JSONL interface and rejects query-label fields. Three new
tests verify the bounds, invariance, protein independence and label firewall.

**NOT PROMOTED.** The PKIS/KLIFS predictor, task-specific section variants,
seven-point categorical-law bridge, PLIP probes, radial affinity heads and all
failed/mixed experimental trainers were not moved into production. The
seven-point bridge is not type-equivalent to the frozen 33-point CDF-band
operator, and F6I's registered total Gate remains not admissible.

**MODEL/SCRIPT AUDIT.** The pre-existing nine `model/` files were all part of
the frozen operator or P1B geometry surface. The pre-existing 36 `scripts/`
files were all part of passed canonical-data, structure-geometry or
release-governance workflows. No additional production file was deleted merely
for having few textual references.

**REPORT CONSOLIDATION.** Added `report/VERIFIED_EVIDENCE_SUMMARY.md`. Retained
only the current split protocol, P1B PASS checkpoint/control artifacts and
D0-C/D1 corpus/governance evidence as active report objects. Historical
research reports now resolve to this ledger and Git history.

**DATA AUTHORIZATION BOUNDARY.** New public data may be downloaded, but affinity
outcome access still requires a preregistration fixing release, endpoint,
closure, controls and Gate. No unspecified dataset was downloaded during this
cleanup, avoiding post-hoc selection of a convenient panel.

**CORE PROBLEM AFTER CONSOLIDATION.** Correct-protein geometry is identified,
and a gauge-separated component interface now exists, but no fresh,
endpoint-consistent source has demonstrated that the biological surface yields
a transferable improvement over both ligand-only and wrong-protein controls.
Consequently affinity energetics and biological-statistic admission to `z`
remain unidentified.


---

## XP1 / XP2 External Crossed-Panel Programme (2026-08-08)

Registered in `research/crossed_panel_identification/PREREG_XP1.md` and
`research/crossed_panel_deployability/PREREG_XP2.md`. Reports in
`report/crossed_panel_identification/XP1_RESEARCH_REPORT.md` and
`report/crossed_panel_deployability/XP2_FINAL_REPORT.md`.

```text
XP1 DECISION: OBJECTIVE_OR_PARAMETERIZATION_FAILURE
XP2 DECISION: CROSSED_INTERACTION_REPRODUCED
              K_LE_5_SECTION_NOT_IDENTIFIED
              PANEL_LOCAL_LOW_RANK_META_LEARNING
              BIOLOGICAL_LANDING_NOT_IDENTIFIED
```

**DATA.** Three external release-pinned panels were acquired. DAVIS, PKIS2 and
Anastassiadis were excluded by governance and never downloaded. No ChEMBL37
affinity value was read at any point; DAVIS and recipient label reads remain 0.

| release | SHA-256 | role |
|---|---|---|
| `metz.xls` (Metz 2011 Nat Chem Biol Table S1) | `81731c4004823bd45fa3898e25d6491d799dfd0e0486fcc8c9c821f9419dd591` | XP1 provenance; XP2 primary source |
| `metz_matrix.csv` | `abe1e3c580478775a352ec5ee78ca565d4c863f0e3e642fdb21d956d8f9d4375` | XP1 primary; XP2 cross-check |
| `aan4368_Table_S2.xlsx` (Klaeger 2017 Science) | `d28b91e62e78e5e011b60da27672875621fef5cdabbea793ac9cce4b98db2c32` | provenance |
| `klaeger_matrix.csv` | `cdf66c7d4e7c1e3a35aeb6995abbfdaf15be80f3e07715524b2bb4449d871010` | cross-platform replication; XP2-F external |
| `KiDatabase.csv` (NIMH PDSP) | `45c9a18ac30f1fad350d1dde186bc1f226c5a75d474ca50f50713852a5637ac6` | independent protein class |

**XP1 POSITIVE RESULTS (retained).** On `BLK-METZ-60` (704 x 82, 34,764 measured
cells) the protein-by-ligand interaction is `59.6%` of affinity variance with
`38%` of the residual reproducible (implied interaction sd `0.442` log units).
The protein-side interaction geometry replicates at `r = 0.8849` across disjoint
compound halves and at `r = 0.5650` across an independent platform with no shared
compounds (label-permutation `p < 5e-4` in both). Under strict KLIFS-group
closure a rank-1..3 support-identified section reached `R2_gamma = +0.160
[+0.109, +0.195]` with derangement specificity `+0.0695 [+0.0524, +0.0874]`,
while every zero-shot protein representation (ESM-2 t30, aligned KLIFS pocket,
pocket physicochemistry, KLIFS conformational state, family, group, homolog kNN)
was indistinguishable from a random-feature null. A synthetic additive panel with
matched margins, noise and truncation gave derangement specificity `+0.00004
[-0.0005, +0.0011]`, falsifying the censoring-artefact explanation.

**XP2-A CORRECTION TO XP1's CENSORING DESCRIPTION.** The journal supplement
encodes 103,118 measured cells, 154,175 left-censored strings at 50 distinct
thresholds (`4.0`-`6.2`), and 405,482 untested blanks; the derived matrix
collapses the latter two to `4.0`. XP1's mask nevertheless admitted 49,457 cells
of which all 49,457 are genuinely measured, matching the supplement at max
`|diff| = 0.0`. XP1's analysis set is therefore correct and no XP1 conclusion
changes; only its single-floor censoring *description*, and the censoring model
used in its destructive control, were approximations. Verdict recorded as
`XP1_EVIDENCE_REPRODUCED`, 18/18 checks.

**XP2 NEGATIVE RESULTS (terminal for this mechanism).** On `BLK-METZ-XP2`
(928 compounds x 147 kinases, 32,849 measured cells, 258 ECFP4-merged
Bemis-Murcko scaffold components, index SHA-256 `7bcb2c05daa4aa5a...`):

1. *Ligand landing succeeded.* Gauge-invariant interaction reconstruction on
   unseen scaffolds reached `R2 = +0.199 [+0.133, +0.261]` (ECFP, `d = 3`)
   against random-feature `+0.025` and mean-loading `+0.024`. The verdict
   `LIGAND_SIDE_DEPLOYMENT_REPRESENTATION_FAILED` does **not** apply.
2. *The identifiable section dimension is exactly `min(k-1, d)`.* Measured
   `0, 1, 2, 3, 3` at `k = 1..5`. At `k = 1` the ridge returns `v = 0` and the
   section arm is identically the additive arm.
3. *At `k <= 5` the section is below the frozen non-negligibility floor.* Under
   protein-group closure only, `R2_gamma` peaks at `+0.0248 [+0.0114, +0.0321]`
   against a registered floor of `0.05`, and `Delta_deploy` over ligand-only
   chemistry never clears zero.
4. *Under simultaneous protein-group and ligand-scaffold closure the section
   loses target specificity altogether.* At the registered primary configuration
   (`k = 5`, `d = 3`): `R2_gamma = +0.0199 [+0.0076, +0.0283]`,
   derangement specificity `+0.00185 [-0.00477, +0.00552]` (CI spans zero), and
   `Delta_deploy = -0.0331 [-0.09108, +0.02326]`. Support design rank was
   `3.00/3` with query coverage `1.000`, so this is not an identifiability
   artefact. XP1's specificity was conditional on ligand reuse.

**GAUGE.** Coordinate-wise loading `R2` in a fold-local gauge was ~0 or negative
in the same runs where the gauge-invariant reconstruction `R2` was clearly
positive. Latent factor coordinates are not stable objects and must never be
assigned hydrogen-bond, hydrophobic, DFG or other biological names. Only fixed
named features and gauge-invariant objects are interpretable.

**INTERFACE.** `research/crossed_panel_deployability/THEORY_INTERFACE_AUDIT.md`
records that the candidate seven-tuple statistic is interface-legal only
conditional on (i) a declared gauge, since `query_coverage` and
`inverse_conditioning` are not `GL(d)`-invariant, (ii) a two-term outer radius
whose second term bounds the unidentified component rather than letting ridge
zero it, and (iii) placement of `validity_flag` and `support_rank` in the finite
context map `kappa` rather than in the continuous sieve coordinates. Abstention
is the existing `p = e_0` simplex vertex and needs no new operator. CSMO, Band,
`K` and the mesh were not modified, and `model/`, `scripts/`, `contracts/` and
`theory/` were not touched at any point in XP1 or XP2.

**NOT ESTABLISHED.** Affinity energetics, biological `z` admission, probability-law
calibration of `K(B(z)F(z))` (never scored in either stage), protein-side
biological landing, and any end-to-end DTA claim.

**DELETION RECORD.** See `report/crossed_panel_deployability/XP2_FINAL_REPORT.md`
section 10 for the disposition of failed implementations. All immutable metrics
are preserved above and in the JSON artifacts under
`report/crossed_panel_identification/` and `report/crossed_panel_deployability/`.

### XP1/XP2 immutable code and artifact hashes (SHA-256, first 32 hex)

| path | sha256[:32] | bytes |
|---|---|---|
| `research/crossed_panel_identification/PREREG_XP1.md` | `51fe525c1171cd3720b2bc606e818808` | 14782 |
| `research/crossed_panel_identification/acquire_kinase_panels.py` | `0f360d02ac0b5f8732b35d575d062e50` | 1889 |
| `research/crossed_panel_identification/acquire_klifs.py` | `2103f052e17838bea0c0c6bda4e55e09` | 1565 |
| `research/crossed_panel_identification/acquire_pdsp.py` | `bffe2d5f0dbcf6d56dca5e1be4ebb6d6` | 858 |
| `research/crossed_panel_identification/build_conformation_features.py` | `f995cf66fa0e620dfc921898da43bb9d` | 4050 |
| `research/crossed_panel_identification/build_protein_features.py` | `f91297f5395217a45a2025d475883bf4` | 3961 |
| `research/crossed_panel_identification/lowrank.py` | `94ee1c86930be9c5d9c6aa0764f553ed` | 3582 |
| `research/crossed_panel_identification/panels.py` | `82b9006fbc5ca3e35f00375933bc1624` | 6843 |
| `research/crossed_panel_identification/pdsp_build.py` | `4fafd1cbe0e57c64b8dbf7977d4c06f3` | 5192 |
| `research/crossed_panel_identification/xp1a_existence.py` | `1e04af414794a2a6039278417b841c10` | 9498 |
| `research/crossed_panel_identification/xp1b_sweep.py` | `cad3c158f3f8142c085236fbb9d0bc12` | 1738 |
| `research/crossed_panel_identification/xp1b_transfer.py` | `f16ae55e939f834a53b42a329dd07944` | 15956 |
| `research/crossed_panel_identification/xp1c_pdsp.py` | `fb7721cadc3076219328e62e4a4a7b5d` | 9312 |
| `research/crossed_panel_identification/xp1d_statistic.py` | `a1ce7033964104db7b40e6f9f0b1a454` | 8375 |
| `research/crossed_panel_identification/xp1e_truncation_control.py` | `baef67903924ad8432ed49283c195a40` | 4373 |
| `research/crossed_panel_deployability/PREREG_XP2.md` | `dd50ec9aedb3df5ac0abeb79a2f3b997` | 11702 |
| `research/crossed_panel_deployability/THEORY_INTERFACE_AUDIT.md` | `4044433219139ade5a8bd759f7625e5c` | 8645 |
| `research/crossed_panel_deployability/acquire_klaeger_structures.py` | `1fd618ea570c71b12949d19f61539565` | 2747 |
| `research/crossed_panel_deployability/xp2_core.py` | `ee873cc7d53f3967d55881bded22aa7f` | 6639 |
| `research/crossed_panel_deployability/xp2_finalize.py` | `470921ed4f6a29911f3690d4adb50531` | 7175 |
| `research/crossed_panel_deployability/xp2_panel.py` | `b9494542315d24c822fec67db9006e28` | 9105 |
| `research/crossed_panel_deployability/xp2a_reproduction_audit.py` | `0c1eaa2d2bc6228144d5f433d55b5312` | 14317 |
| `research/crossed_panel_deployability/xp2b_ligand_landing.py` | `2f26a1e404bb7c7e89cfc3831cd53214` | 7144 |
| `research/crossed_panel_deployability/xp2cd_section.py` | `6789236a9cd9fba9c2eef8002c0bc4c1` | 11639 |
| `research/crossed_panel_deployability/xp2cd_sweep.py` | `c6888040ce7868d5e4cd0f4469d750bc` | 1821 |
| `research/crossed_panel_deployability/xp2e_landing.py` | `305437af1f3d26366a90432b088b9c0c` | 9644 |
| `research/crossed_panel_deployability/xp2f_external.py` | `337df4d60bd9b606e2c8d88575c42c56` | 10460 |
| `report/crossed_panel_identification/XP1_RESEARCH_REPORT.md` | `e8c084bd076ff3d48b96dc181fe9b4d7` | 43629 |
| `report/crossed_panel_identification/xp1a_console.txt` | `60e61e7e96fe220c8d39d3de6d1c1aba` | 2813 |
| `report/crossed_panel_identification/xp1a_existence.json` | `bfb1b8b0d57bc7b57fdf45ac2c09fff2` | 7603 |
| `report/crossed_panel_identification/xp1b_sweep_console.txt` | `95cc9801e7537b1e5b515398e322fe3f` | 49764 |
| `report/crossed_panel_identification/xp1b_sweeps.json` | `51e72e33bba99658d0aa410b2f55ad79` | 151643 |
| `report/crossed_panel_identification/xp1c_console.txt` | `69612a8b1706b3c3e98d1b97da328d29` | 1655 |
| `report/crossed_panel_identification/xp1c_pdsp.json` | `21c66f8a15f0de5a4bf6bef8961b6c2c` | 4544 |
| `report/crossed_panel_identification/xp1d_console.txt` | `730a5592b08b6ded8d578cea797f2ff2` | 1974 |
| `report/crossed_panel_identification/xp1d_statistic.json` | `307a3b4a1a925fa0fbda0856d78276be` | 5689 |
| `report/crossed_panel_identification/xp1e_console.txt` | `83bf7e80a86ab91da08b4f6597e1f618` | 3630 |
| `report/crossed_panel_identification/xp1e_truncation_control.json` | `ab836d6c3d1a6691cedfc22dc46a4568` | 42095 |
| `report/crossed_panel_deployability/DOUBLE_HELD_OUT_RESULT.json` | `70bcc987f600944b8c52e3fc46002f73` | 87671 |
| `report/crossed_panel_deployability/EXTERNAL_REPLICATION_RESULT.json` | `cfc8fda9433d425aec26f9eb145bd9ed` | 3891 |
| `report/crossed_panel_deployability/K5_SECTION_AUDIT.json` | `5d46a245dd6d2ec962e7c1a12fe33cc5` | 33503 |
| `report/crossed_panel_deployability/LIGAND_LANDING_AUDIT.json` | `94a3a621117283f47d006315eb8f380b` | 18749 |
| `report/crossed_panel_deployability/XP1_REPRODUCTION_AUDIT.json` | `8ce175b3271d5764d253a28e116303ad` | 5905 |
| `report/crossed_panel_deployability/XP2E_BIOLOGICAL_LANDING.json` | `00f9a41e64fc0018fbb3ded5b1f7494b` | 4557 |
| `report/crossed_panel_deployability/XP2_FINAL_REPORT.md` | `a15bd43b76f1f73579431665f1609647` | 27990 |
| `report/crossed_panel_deployability/xp2b_console.txt` | `c98eee340505b2119e6b03bb57f2aea4` | 2483 |
| `report/crossed_panel_deployability/xp2cd_console.txt` | `2ec278c09f03b4badbbb891138ed2b3f` | 30778 |
| `report/crossed_panel_deployability/xp2cd_sweeps.json` | `876882af074e3c6af9a6df88428dd73f` | 141591 |
| `report/crossed_panel_deployability/xp2e_console.txt` | `86b88dc53750124752dc5c3ad719c737` | 1045 |
| `report/crossed_panel_deployability/xp2f_console.txt` | `80bd510d103d8438954bb374fe9f5609` | 1325 |

Environment: python 3.11.15, numpy 1.26.4, scipy 1.17.1, pandas 2.3.3, rdkit 2023.09.6, torch 2.6.0+cu124, transformers 4.46.3, scikit-learn 1.9.0; `xlrd` was installed to read the `.xls` supplement. Seeds `{0,1,2,3,4}` throughout; bootstrap seeds are fixed per contrast in source.

**Regression suite: `73 passed` before and after XP2.** `model/`, production `scripts/`, `contracts/` and `theory/` show no modification under `git status` for the whole programme.

Upstream release licences: Metz 2011 and Klaeger 2017 supplements are publisher supplementary data accessed through a public mirror pinned at commit `8ab79cae31c18e49007dcce6dd11f93d2667ab14`; the NIMH PDSP Ki database is a free public NIMH resource; KLIFS is open academic access; PubChem PUG-REST was used for name-to-structure resolution only, with `affinity_values_read = 0`.

### XP1/XP2 implementation removal (2026-08-08)

Recovery commit `3281780` on branch `research/xp1-xp2-crossed-panel` holds the
complete reproducible tree. The following terminal-negative or negligible-effect
implementations were then removed from the active surface. Their conclusions,
metrics and artifact hashes are preserved above and their JSON artifacts are
retained under `report/`.

| removed file | conclusion preserved in | verdict it produced |
|---|---|---|
| `research/crossed_panel_identification/xp1c_pdsp.py` | XP1 report section 6.4, `xp1c_pdsp.json` | PDSP replication below the non-negligibility floor |
| `research/crossed_panel_identification/xp1d_statistic.py` | XP1 report section 6.6, `xp1d_statistic.json` | no protein representation predicts the interaction coordinate across groups |
| `research/crossed_panel_identification/build_conformation_features.py` | XP1 report sections 5 and 6.4 | KLIFS conformational-state arm disqualified on availability (16/82 kinases lack any structure) and circularity |
| `research/crossed_panel_deployability/xp2e_landing.py` | XP2 report section 7, `XP2E_BIOLOGICAL_LANDING.json` | `BIOLOGICAL_LANDING_NOT_IDENTIFIED` |
| `research/crossed_panel_deployability/xp2f_external.py` | XP2 report section 6.3, `EXTERNAL_REPLICATION_RESULT.json` | `EXTERNAL_REPLICATION_FAILED` |

Retained deliberately: `panels.py`, `lowrank.py`, `xp1a_existence.py`,
`xp1b_transfer.py`, `xp1b_sweep.py`, `build_protein_features.py`,
`xp1e_truncation_control.py` (a control that PASSED), the acquisition scripts,
and the whole XP2 core, because `xp2a_reproduction_audit.py` reproduces XP1 by
reading `xp1b_transfer.py` directly and the destructive control is structurally
reusable.

### XP1/XP2 immutable code and artifact hashes (SHA-256, first 32 hex)

| path | sha256[:32] | bytes |
|---|---|---|
| `research/crossed_panel_identification/PREREG_XP1.md` | `51fe525c1171cd3720b2bc606e818808` | 14782 |
| `research/crossed_panel_identification/acquire_kinase_panels.py` | `0f360d02ac0b5f8732b35d575d062e50` | 1889 |
| `research/crossed_panel_identification/acquire_klifs.py` | `2103f052e17838bea0c0c6bda4e55e09` | 1565 |
| `research/crossed_panel_identification/acquire_pdsp.py` | `bffe2d5f0dbcf6d56dca5e1be4ebb6d6` | 858 |
| `research/crossed_panel_identification/build_protein_features.py` | `f91297f5395217a45a2025d475883bf4` | 3961 |
| `research/crossed_panel_identification/lowrank.py` | `94ee1c86930be9c5d9c6aa0764f553ed` | 3582 |
| `research/crossed_panel_identification/panels.py` | `82b9006fbc5ca3e35f00375933bc1624` | 6843 |
| `research/crossed_panel_identification/pdsp_build.py` | `4fafd1cbe0e57c64b8dbf7977d4c06f3` | 5192 |
| `research/crossed_panel_identification/xp1a_existence.py` | `1e04af414794a2a6039278417b841c10` | 9498 |
| `research/crossed_panel_identification/xp1b_sweep.py` | `cad3c158f3f8142c085236fbb9d0bc12` | 1738 |
| `research/crossed_panel_identification/xp1b_transfer.py` | `f16ae55e939f834a53b42a329dd07944` | 15956 |
| `research/crossed_panel_identification/xp1e_truncation_control.py` | `baef67903924ad8432ed49283c195a40` | 4373 |
| `research/crossed_panel_deployability/PREREG_XP2.md` | `dd50ec9aedb3df5ac0abeb79a2f3b997` | 11702 |
| `research/crossed_panel_deployability/THEORY_INTERFACE_AUDIT.md` | `4044433219139ade5a8bd759f7625e5c` | 8645 |
| `research/crossed_panel_deployability/acquire_klaeger_structures.py` | `1fd618ea570c71b12949d19f61539565` | 2747 |
| `research/crossed_panel_deployability/xp2_core.py` | `ee873cc7d53f3967d55881bded22aa7f` | 6639 |
| `research/crossed_panel_deployability/xp2_finalize.py` | `470921ed4f6a29911f3690d4adb50531` | 7175 |
| `research/crossed_panel_deployability/xp2_panel.py` | `b9494542315d24c822fec67db9006e28` | 9105 |
| `research/crossed_panel_deployability/xp2a_reproduction_audit.py` | `0c1eaa2d2bc6228144d5f433d55b5312` | 14317 |
| `research/crossed_panel_deployability/xp2b_ligand_landing.py` | `2f26a1e404bb7c7e89cfc3831cd53214` | 7144 |
| `research/crossed_panel_deployability/xp2cd_section.py` | `6789236a9cd9fba9c2eef8002c0bc4c1` | 11639 |
| `research/crossed_panel_deployability/xp2cd_sweep.py` | `c6888040ce7868d5e4cd0f4469d750bc` | 1821 |
| `report/crossed_panel_identification/XP1_RESEARCH_REPORT.md` | `e8c084bd076ff3d48b96dc181fe9b4d7` | 43629 |
| `report/crossed_panel_identification/xp1a_console.txt` | `60e61e7e96fe220c8d39d3de6d1c1aba` | 2813 |
| `report/crossed_panel_identification/xp1a_existence.json` | `bfb1b8b0d57bc7b57fdf45ac2c09fff2` | 7603 |
| `report/crossed_panel_identification/xp1b_sweep_console.txt` | `95cc9801e7537b1e5b515398e322fe3f` | 49764 |
| `report/crossed_panel_identification/xp1b_sweeps.json` | `51e72e33bba99658d0aa410b2f55ad79` | 151643 |
| `report/crossed_panel_identification/xp1c_console.txt` | `69612a8b1706b3c3e98d1b97da328d29` | 1655 |
| `report/crossed_panel_identification/xp1c_pdsp.json` | `21c66f8a15f0de5a4bf6bef8961b6c2c` | 4544 |
| `report/crossed_panel_identification/xp1d_console.txt` | `730a5592b08b6ded8d578cea797f2ff2` | 1974 |
| `report/crossed_panel_identification/xp1d_statistic.json` | `307a3b4a1a925fa0fbda0856d78276be` | 5689 |
| `report/crossed_panel_identification/xp1e_console.txt` | `83bf7e80a86ab91da08b4f6597e1f618` | 3630 |
| `report/crossed_panel_identification/xp1e_truncation_control.json` | `ab836d6c3d1a6691cedfc22dc46a4568` | 42095 |
| `report/crossed_panel_deployability/DOUBLE_HELD_OUT_RESULT.json` | `70bcc987f600944b8c52e3fc46002f73` | 87671 |
| `report/crossed_panel_deployability/EXTERNAL_REPLICATION_RESULT.json` | `cfc8fda9433d425aec26f9eb145bd9ed` | 3891 |
| `report/crossed_panel_deployability/K5_SECTION_AUDIT.json` | `5d46a245dd6d2ec962e7c1a12fe33cc5` | 33503 |
| `report/crossed_panel_deployability/LIGAND_LANDING_AUDIT.json` | `94a3a621117283f47d006315eb8f380b` | 18749 |
| `report/crossed_panel_deployability/XP1_REPRODUCTION_AUDIT.json` | `8ce175b3271d5764d253a28e116303ad` | 5905 |
| `report/crossed_panel_deployability/XP2E_BIOLOGICAL_LANDING.json` | `00f9a41e64fc0018fbb3ded5b1f7494b` | 4557 |
| `report/crossed_panel_deployability/XP2_FINAL_REPORT.md` | `5d7a5f47f571d41d3cfbe75b914aac8a` | 29674 |
| `report/crossed_panel_deployability/xp2b_console.txt` | `c98eee340505b2119e6b03bb57f2aea4` | 2483 |
| `report/crossed_panel_deployability/xp2cd_console.txt` | `2ec278c09f03b4badbbb891138ed2b3f` | 30778 |
| `report/crossed_panel_deployability/xp2cd_sweeps.json` | `876882af074e3c6af9a6df88428dd73f` | 141591 |
| `report/crossed_panel_deployability/xp2e_console.txt` | `86b88dc53750124752dc5c3ad719c737` | 1045 |
| `report/crossed_panel_deployability/xp2f_console.txt` | `80bd510d103d8438954bb374fe9f5609` | 1325 |

Environment: python 3.11.15, numpy 1.26.4, scipy 1.17.1, pandas 2.3.3, rdkit 2023.09.6, torch 2.6.0+cu124, transformers 4.46.3, scikit-learn 1.9.0; `xlrd` was installed to read the `.xls` supplement. Seeds `{0,1,2,3,4}` throughout; bootstrap seeds are fixed per contrast in source.

**Regression suite: `73 passed` before and after XP2.** `model/`, production `scripts/`, `contracts/` and `theory/` show no modification under `git status` for the whole programme.

Upstream release licences: Metz 2011 and Klaeger 2017 supplements are publisher supplementary data accessed through a public mirror pinned at commit `8ab79cae31c18e49007dcce6dd11f93d2667ab14`; the NIMH PDSP Ki database is a free public NIMH resource; KLIFS is open academic access; PubChem PUG-REST was used for name-to-structure resolution only, with `affinity_values_read = 0`.


## S0-S4 Structural Self-Supervision Programme (2026-08-08)

Registered under `research/ssl_b2_structural_observability/`; artifacts in
`report/ssl_b2_structural_observability/`. Branch `research/ssl-b2-structural`.
Entirely label-free: DAVIS 0, recipient 0, ChEMBL37 affinity 0, any affinity
value 0.

**S1 independent structural test set.** RCSB Search (X-ray, `<= 2.5 A`, bound
non-polymer, one protein entity, released `>= 2024-01-01`) returned 15,003
entries with **exposed overlap 0** against all 10,468 pilot20k PDB ids. 1,476
acquired under CC0-1.0, 1,162 usable complexes, 1,118 with a parsable CCD
ligand. Final block: **1,118 complexes, 621 MMseqs40/80cov protein clusters,
586 Bemis-Murcko scaffolds.** All P1B-exposed ids, including P1B's own val and
test partitions, were treated as exposed. PLINDER and PDBbind recorded as
deliberately not used, with reasons.

**S2 teacher.** Six named channels computed from raw holo coordinates:
directional H-bond, signed electrostatics, hydrophobic burial, aromatic
orientation, steric overlap, pocket burial. Reproducibility audit on 78
complexes: rotation/translation `7.6e-15`, atom permutation `1.6e-14`,
determinism exactly `0.0`, tolerance `1e-9`, **no channel degenerate**. An
implementation defect was found and fixed after the contract was frozen: the
RCSB filter admitted homo-oligomers (one entity, many chains) and the O(nP^2)
neighbour step reached 10.5 GB, so the teacher restricts to protein residues
within 10 A of the ligand. Every channel is defined at `<= 8 A` and bonded
neighbours are `<= 1.8 A`, so this is exactly equivalent; the full invariance
audit was re-run afterwards and channel statistics were unchanged.

**S3 power.** 124 effective independence units per fold; minimum detectable
`R2 = 0.02` at 100% detection under the registered decision rule. The frozen S6
effect floor is `0.02`. **The design is not underpowered**, so the null
contrasts below are real nulls.

**S4 observability, upper bound on the sequence+2D class.** Deviation recorded:
rather than running the P1B checkpoint, the audit asks whether the six channels
are reachable from ESM-2 + ECFP at all. P1B's predicted geometry is a function
of exactly those inputs, so a negative bounds the whole class rather than one
model.

| channel | R2 vs mean | vs random | vs deranged protein |
|---|---|---|---|
| hbond_directional | **+0.268 [+0.166, +0.378]** | **+0.366 [+0.222, +0.505]** | +0.037 [-0.015, +0.084] |
| hydrophobic_burial | **+0.299 [+0.162, +0.454]** | **+0.307 [+0.167, +0.431]** | -0.006 [-0.035, +0.024] |
| steric_overlap | +0.079 [-0.012, +0.137] | +0.058 [-0.005, +0.132] | +0.070 [-0.009, +0.184] |
| pocket_burial | +0.055 [-0.029, +0.138] | +0.052 [-0.029, +0.137] | +0.029 [-0.045, +0.098] |
| aromatic_orientation | +0.033 [-0.013, +0.069] | +0.071 [+0.014, +0.133] | +0.010 [-0.026, +0.041] |
| electrostatic_signed | +0.027 [-0.015, +0.048] | +0.044 [-0.004, +0.103] | -0.006 [-0.031, +0.018] |

Two channels are genuinely observable from deployment inputs and clearly beat
capacity-matched random features. **No channel beats the deranged-protein
control**; every interval spans zero. Substituting a foreign protein's embedding
costs nothing, which localises the predictive information to the ligand side.

**GPU training was NOT authorised.** Three of the four S4 preconditions are met
(teacher reproducible, no channel degenerate, measurable information above
random), but the information is not protein-specific, so a distillation network
would learn ligand chemistry - exactly the population shortcut the programme
forbids. No GPU training was performed at any point in S0-S4; GPU was used only
for frozen ESM-2 inference.

**S4b attribution — the decisive control.** Ligand-only, protein-only and joint
arms on the same cells and split, identical hyperparameter selection:

| channel | LIG-ONLY | PROT-ONLY | BOTH - LIG |
|---|---|---|---|
| hbond_directional | +0.266 [+0.160, +0.391] | +0.009 [-0.026, +0.047] | +0.0015 [-0.0379, +0.0363] |
| hydrophobic_burial | +0.331 [+0.204, +0.485] | -0.017 [-0.044, +0.004] | **-0.0321 [-0.0621, -0.0052]** |
| pocket_burial | +0.077 [+0.027, +0.125] | -0.027 [-0.080, +0.023] | -0.0200 |
| aromatic_orientation | +0.044 [+0.006, +0.078] | +0.000 [-0.010, +0.010] | -0.0137 |
| steric_overlap | +0.043 [+0.009, +0.114] | +0.012 [-0.047, +0.052] | +0.0385 |
| electrostatic_signed | +0.036 [+0.006, +0.054] | +0.001 [-0.013, +0.012] | -0.0101 |

Protein-only is ~0 for all six. The ligand alone explains everything the joint
model explains; adding the protein buys +0.0015 on the best channel (CI spans
zero) and costs -0.032 on hydrophobic burial with a CI excluding zero. The
observable teacher signal is a ligand descriptor.

```text
S-PROGRAMME TERMINAL VERDICT: POSE_FREE_DEPLOYMENT_INPUTS_INSUFFICIENT
```

GPU training was not authorised: three of the four S4 conditions were met, but
the probe reaches R2=0.30 where information exists (so it is not underfitting)
and that information is ligand-side (so distillation would learn chemistry, the
forbidden population shortcut). No affinity value was read; S8 was never entered.
Full account: `report/ssl_b2_structural_observability/S_PROGRAMME_REPORT.md`.

### Post-run interpretation correction and S5 redesign

The S4/S4b numbers above are preserved, but the claim that they form an upper
bound on the complete sequence+2D model class is withdrawn.  The probe used
mean-pooled ESM, ECFP and linear Ridge; it omitted the existing P1B atom-local
GINE states, residue-local ESM states, and atom-by-slot contact/distance tensor.
It therefore establishes only that the tested **aggregate** representation is
ligand-dominated and not protein-specific.  It does not close the pose-free
pair-local route.

Corrected status:

```text
AGGREGATE_ESM_ECFP_PROBE_NOT_PROTEIN_SPECIFIC
AGGREGATE_TEACHER_SIGNAL_IS_LIGAND_DOMINATED
PAIR_LOCAL_P1B_OBSERVABILITY_NOT_TESTED
POSE_FREE_CLASS_NOT_CLOSED
```

The redesign also records two new fail-closed audits.  First, one protein entity
can contain several homo-oligomer chains, while deployment provides one target
sequence; single-chain and interface complexes must be separated.  Second, the
six S2 channels are deterministic structural pseudo-labels rather than physical
ground truth, and their ligand atom chemistry must be rebuilt from a mapped CCD
bond/charge/donor/acceptor/ring contract before mechanistic claims.

`PREREG_S5_LOCAL_MECHANISM_OBSERVABILITY.md` freezes the continuation.  It
requires exact ligand-atom and residue-sequence mapping, a slot-information
ceiling, an actual frozen-P1B observability ladder, a synthetic optimization
control, and only conditionally a small pair-local GPU head.  The existing
1,118-complex S4 block is development-exposed; a new score-blind RCSB block must
be sealed for confirmation.  No affinity label, DAVIS label, biological `z`,
production model change, CSMO/Band change, or P2-P4 authorization follows.

### S0-S4 immutable code and artifact hashes

| path | sha256[:32] | bytes |
|---|---|---|
| `research/ssl_b2_structural_observability/LICENSE_AND_PROVENANCE_AUDIT.md` | `33cb45bfca2a35c321f55ff8a5379bdc` | 4232 |
| `research/ssl_b2_structural_observability/TEACHER_CONTRACT.md` | `5c5763137946630017c0b9c303bd6009` | 4622 |
| `research/ssl_b2_structural_observability/s1_registry_and_exposure.py` | `7081076ed826880d5a932ff8a7925f18` | 6783 |
| `research/ssl_b2_structural_observability/s1b_acquire_independent.py` | `e1b4fdbb2bed0e8b4c3de288830e0680` | 4841 |
| `research/ssl_b2_structural_observability/s2_teacher.py` | `1199bde00e502413636f5c57ea687cb6` | 10002 |
| `research/ssl_b2_structural_observability/s3_power.py` | `3eb30b6854feed90ddb27023ab2f979b` | 3308 |
| `research/ssl_b2_structural_observability/s3s4_observability.py` | `36d725ba832534d1bc5d7a6f2ca3153c` | 13821 |
| `research/ssl_b2_structural_observability/s4b_attribution.py` | `4090ba950b6a843013a14470075f3c67` | 4870 |
| `report/ssl_b2_structural_observability/DATASET_ROLE_REGISTRY.json` | `63ad6710a6909e7c41ff0f4f32ab689c` | 3380 |
| `report/ssl_b2_structural_observability/S4B_ATTRIBUTION.json` | `83198f9931cfaf8f9b89229b404063be` | 6181 |
| `report/ssl_b2_structural_observability/S4_OBSERVABILITY_AUDIT.json` | `4ed05b8c683a861150adffefc0f5a5c5` | 5427 |
| `report/ssl_b2_structural_observability/STRUCTURAL_EXPOSURE_AUDIT.json` | `02c7774f68d45aabb893162a1b0db901` | 676 |
| `report/ssl_b2_structural_observability/STRUCTURAL_POWER_ANALYSIS.json` | `873728e0cc53d4c1d01962af4c0f86f2` | 599 |
| `report/ssl_b2_structural_observability/STRUCTURAL_SPLIT_MANIFEST.json` | `23e76f98a77304abf90232f5b57fe9d2` | 569 |
| `report/ssl_b2_structural_observability/S_PROGRAMME_REPORT.md` | `b8be2b37de541fd1a48a67b290cbb8a9` | 12794 |
| `report/ssl_b2_structural_observability/TEACHER_REPRODUCIBILITY_AUDIT.json` | `3544d05cdb971d3be1a25d6f87e3d0f7` | 1624 |
| `report/ssl_b2_structural_observability/s1b_console.txt` | `9e581b1faac5d63ed0f9fd530a208f45` | 520 |
| `report/ssl_b2_structural_observability/s3_power_console.txt` | `b659122bcad2f785b02c09277013f68f` | 601 |
| `report/ssl_b2_structural_observability/s4_console.txt` | `f279e1fd1b5a5c37798584de2e316e5e` | 5500 |
| `report/ssl_b2_structural_observability/s4b_console.txt` | `b8f3322c038ac001b6d5d7f48141dd82` | 852 |

Independent structural release: RCSB PDB CC0-1.0, 1,476 entries released >= 2024-01-01, acquisition manifest with per-file SHA-256 at `dataset/raw/ssl_b2_independent/acquisition_manifest.json`. gemmi 0.7.5, MMseqs2 repo-pinned. Seeds fixed. GPU used for frozen ESM-2 inference only; no GPU training was performed.

## Global consolidation after S5 redesign (2026-08-08)

The repository was reduced to passed production/data/geometry surfaces and the
single active S5 preregistration.  Terminal-negative research implementations
and duplicate reports were removed from the working tree; exact recovery points
remain Git commits `3281780` (XP1/XP2 recovery), `12a2765` (XP3/XP4/XP5
boundary), and `608decf` (S0-S4 implementation and raw report).

The deleted XP3 census recorded the governing public-data tradeoff: Metz had
928 compounds, 147 proteins but only 8 protein-group closure components;
Klaeger was 93.6% at the measurement floor; PDSP per-report noise was 0.7144 log
units.  The deleted XP4 implementation formed 85 BindingDB panels, 70 protein
clusters and 6,363 cells, but its interaction SD was 0.4058 versus estimated
per-report noise 0.7774 and the bilinear arm had `R2_gamma=-0.00072` relative to
the additive null.  The deleted XP5 pose-free typed basis was worse than the
null (`R2_gamma=-0.00147`) and indistinguishable from foreign-protein control.
These are terminal evidence records, not production code.

S4 numerical evidence is retained in this ledger: aggregate H-bond and
hydrophobic pseudo-labels were predictable by mean-pooled ESM+ECFP, but no
channel beat deranged protein and ligand-only matched the joint probe.  Its
class-wide `POSE_FREE_DEPLOYMENT_INPUTS_INSUFFICIENT` interpretation was
withdrawn because the probe omitted P1B local pair features.  The active status
is `PAIR_LOCAL_P1B_OBSERVABILITY_NOT_TESTED`, governed by the S5 preregistration.

Deleted root proposals (`DRP_MODULE_PROPOSAL.md`,
`IDENTIFICATION_ROADMAP_AND_Z_ADMISSION.md`, and
`SOLUTION_MENU_LITERATURE_INFORMED.md`) were superseded by
`report/CURRENT_RESEARCH_STATUS.md`, `task.md`, and the S5 preregistration.
The unadmitted F6I `component_statistic` implementation, its label-safe wrapper,
and their isolated tests were also removed from `model/`, `scripts/`, and
`tests/`: their algebraic invariants passed, but the component never completed
fresh external biological admission and therefore did not satisfy the stricter
production-retention rule.  It remains recoverable from commit `24a9ae0`.
The post-prune regression contains 70 tests and passes in the `drug` environment;
the three removed tests belonged exclusively to that unadmitted component.

## External SSL report audit and contract hardening (2026-08-09)

The supplied `METASIEVE_SSL_DETAILED_ANALYSIS_REPORT_2026-08-08.md` was checked
against every local branch and all GitHub branch tips. The report describes S5
through S9 metrics, but the repository contains no matching source commits,
input/split manifests, checkpoints, prediction files, bootstrap units, or
result JSON. The report itself states that it did not independently rerun those
experiments. The numerical claims are therefore classified as
`EXTERNAL_CLAIM_NOT_REPRODUCED`; formal status remains `S5 REGISTERED_NOT_RUN`.

The theoretical warning is retained: native-decoy NCE depends on the negative
proposal, and coordinate score matching leaves system-dependent additive
energy constants. This makes ensemble and quantitative difference anchors a
reasonable future hypothesis, not an identified explanation of an unverified
S9 result. MISATO can support MD/QM ensemble observability, while PLINDER can
support pocket/ligand/apo leakage governance; neither alone supplies a complete
binding-free-energy gauge.

Four code defects from the report were independently reproduced and repaired:
the removed `biological.kappa` import in `build_band_operator`, silent
off-domain statistic clamping, zero-partner encoder inputs, and incomplete
deployment semantic hashes. These are contract repairs only and do not change
the scientific verdict or admit a biological statistic to `z`. Detailed triage
is in `report/SSL_DETAILED_REPORT_EVIDENCE_AUDIT.md`.

The accompanying `METASIEVE_SSL_TEST_PLAN_2026-08-08.md` was then reviewed as
a prospective protocol. Its A0/A1/A2 separation, U1/U2/U3 hierarchy,
label-late join and destructive controls were retained. Eight corrections were
registered: missing S5-S9 evidence cannot be assumed; the `k`-dimension law is
not a blanket five-dimensional population-state theorem; wrong-protein margin
and coverage are controls/validity objects rather than ordinary `z`
coordinates; the current 28-coordinate context map cannot assemble a new
five-coordinate sieve; an absolute source loss is empirical calibration rather
than proof of a physical gauge; edge-direction shuffle needs an unambiguous
graph perturbation; U3 cannot pass an unbound-reference test; and a
`k`-dimensional adapter alone is not the frozen theory's conservative
outer-section object.

The only new conditional stage is a metadata-only G0 data-feasibility census.
It separates MISATO bound trajectories, PLINDER linked apo candidates, and
quantitative same-target graphs, and keeps all model training and affinity
labels frozen. See `report/SSL_TEST_PLAN_REVIEW_AND_EXECUTION_BOUNDARY.md` and
`research/ssl_gauge_fixed/PREREG_G0_DATA_FEASIBILITY.md`.

## Experimental evidence consolidation (2026-08-09)

Completed reports, raw retained Gates, current state records and all 94 failure
entries were reconciled into
`report/EXPERIMENTAL_EVIDENCE_LEDGER.md`. The ledger separates verified PASS
evidence, negative model evidence, data/estimand underdetermination, experiments
that were not run, and externally supplied claims that were not reproduced.

No historical verdict, raw Gate artifact, model implementation, training script
or frozen theory file was changed. `README.md`, `project_state.json`,
`report/CURRENT_RESEARCH_STATUS.md` and the compact evidence index now point to
the canonical ledger. This file's top current-status block was refreshed from
the obsolete X0/L0 execution state to the active S5 boundary; the detailed
X0/L0 and earlier records remain preserved below. The consolidated regression
suite passed `75` tests.

---

## F-95: S5 stopped at the Stage 0 environment gate; representation-contract findings recorded (2026-08-09)

**OUTCOME.** `S5_EXECUTION_ENVIRONMENT_UNAVAILABLE` and
`S5_ARTIFACT_RECOVERY_REQUIRED`. No terminal S5 decision-tree verdict is
claimed. Stage status remains `S5_REGISTERED_NOT_RUN`.

**BLOCKER B-1 — no execution environment.** The isolated Linux workspace failed
to start on all seven attempts (`session disk not found: ...sessiondata.vhdx`).
No shell, Python, `gemmi`/`parasail`, `git`, SHA-256 hashing, GPU or `pytest`.
Every numeric quantity in Stages 1-5 is `BLOCKED_NO_COMPUTE`. **The regression
suite was not run**; the last known state is the recorded 75-test pass. All
audits this session were read-only; no process was executed.

**BLOCKER B-2 — the registered pseudo-teacher source is deleted.** `research/`
contains zero `.py` files. Only `__pycache__/s2_teacher.cpython-311.pyc` and
`s3s4_observability.cpython-311.pyc` survive the F-94 consolidation. Recovery is
from commit `608decf` / `8b7789e`. Decompiling the bytecode was rejected: it
yields a look-alike, not the registered implementation, which is precisely what
the R-VERIFY rule forbids.

**BLOCKER B-3 — the S4 exposure registry is absent.** `project_state.json`
asserts 1,118 development-exposed complexes, but no file in the working tree
enumerates them. The preregistration requires the new confirmation block to be
disjoint from *both* exposure registries. One does not exist as data, so the
block cannot be sealed score-blind.

**P1B CONTRACT IS RECONSTRUCTABLE (correction).** Contrary to the risk the
Stage 0 rule guards against, no proxy substitution was needed. Present and
declared: `best.pt` (`checkpoint_sha256 90b0010b...`), geometry supervision
(14,906 pairs, 117 hashed shards), frozen ESM 128-slot bank
(`esm2_t30_150M_UR50D` rev `a695f604...`, hidden 640), ligand bank, MMseqs2
homology split (3,606 groups, PASS), governed complex index and 20,295 raw
mmCIF files. Hash *values* are recorded; hash *verification* is blocked.

**STATIC REPRESENTATION-CONTRACT FINDINGS (all EXACT, derived from frozen source,
verifiable by reading, no execution required).**

1. **M-4.** The geometry slot map `min(127, seq_idx*128//L)` and the ESM slot map
   `(i*128)//L` are the same function; the clamp is never active. The two frozen
   banks join without re-indexing.
2. **M-6.** `contact` is a deterministic function of `distance_bin`:
   `contact <=> d<=6`, `bin<=1 <=> d<6`, differing only on the null set
   `{d=6.0}`. Ladder arms A4 and A5 are **not** information-nested on the label
   side; any gap between them measures two frozen heads, not extra structure.
3. **M-7.** The frozen contract carries **no angular information at any stage**.
   Directional hydrogen bonding and aromatic orientation are `NOT_EVALUABLE` at
   the deployable arm under any model, at any capacity.
4. **M-8 (most consequential).** Directional H-bond was S4's best-recovered
   aggregate channel (`R2 = +0.268 [0.166, 0.378]`) **and** is exactly the
   channel the frozen contract cannot carry to A5. Preregistering it as a core
   S5 channel would have produced structurally guaranteed ligand dominance -- a
   representational artefact indistinguishable in the output from a real
   biological null.
5. **Non-identifiability lemma.** For any additive distance-decaying pair-local
   teacher, `Var[T | frozen pair-local state] > 0` whenever a slot holds >=2
   residues at distinct distances, i.e. whenever `L > 128`. Direction proved;
   **magnitude is the preregistered quantity and remains blocked.**
6. **M-2, M-3, D-3 teacher-fidelity defects.** `_protein_residue_atoms` applies
   no `type_symbol` filter, so deposited protein hydrogens enter the
   min-distance while ligand hydrogens are excluded; residues with no tabulated
   one-letter code become `"X"` and align against the BLOSUM62 `X` row instead
   of failing closed, so the preregistered modified-residue exclusion is not
   implemented in the frozen path; the governed corpus contains cofactors and
   metalloporphyrins (e.g. `101m`/`HEM`) that no candidate channel models.

**H-DILUTION registered.** Slots bin *sequence*, not space, so the slot `min` is
generically attained by the pocket residue (geometry survives comparatively
well) while the ESM slot mean gives that residue weight exactly `1/ceil(L/128)`
(identity diluted linearly in `L`). These are confounded in the base ladder. One
arm `A5c` was added -- explicit 20-dim slot residue-type composition replacing
the mean-ESM slot state, ~20 input dimensions, **zero new parameters, zero new
modules** -- to separate `FROZEN_P1B_LOCAL_STATE_INSUFFICIENT` from a
recoverable pooling defect. Neither branch authorises a larger ESM, ESM
fine-tuning, a deeper GNN, cross-attention, a knowledge graph or fusion.

**AMENDMENT.** `research/ssl_b2_structural_observability/PREREG_S5R_PAIR_LOCAL_AMENDMENT.md`
narrows the core channel set to `{K1 pocket burial coverage, K2 apolar-contact
burial}`, demotes directional H-bond and aromatic orientation to oracle-only
`NOT_EVALUABLE` diagnostics, declares A4/A5 non-nested, and registers the three
teacher-fidelity defects before scoring. No threshold was relaxed, no module
added, no label surface widened.

**NOT CLAIMED.** `MAPPING_OR_SLOT_CONTRACT_FAIL_CLOSED` is explicitly **not**
declared: fail-closed requires a measured retention below the preregistered
requirement, and no measurement exists. Likewise none of
`STRUCTURAL_TEACHER_NOT_REPRODUCIBLE`, `SYNTHETIC_TRAINABILITY_DEFECT`,
`LIGAND_SHORTCUT_ONLY`, `FROZEN_P1B_LOCAL_STATE_INSUFFICIENT` or
`PAIR_LOCAL_STRUCTURAL_MECHANISM_OBSERVED` is available from this session.

**LABEL AUDIT.** `recipient_label_reads = 0`, `davis_label_reads = 0`, zero
ChEMBL/BindingDB affinity value reads, zero PLIP or pose-aware reads. Label
firewall intact. Confirmation block not sealed.

**FROZEN SURFACES.** `theory/FINAL_FROZEN_THEORY/`, `model/`, production
`scripts/`, CSMO, Band, positive ridge, the simplex, the fixed mesh,
`K(B(z)F(z))` and production `z` were read only and not modified. Writes were
confined to `report/ssl_s5/`, one new file under
`research/ssl_b2_structural_observability/`, `history.md` and
`project_state.json`.

**ARTIFACTS.** `report/ssl_s5/` -- `S5_STAGE0_REPOSITORY_AND_EVIDENCE_AUDIT.md`,
`S5_MAPPING_AND_SLOT_RETENTION_AUDIT.md`,
`S5_TEACHER_SPECIFICATION_AND_REPRODUCIBILITY.md`,
`S5_LADDER_AND_SYNTHETIC_CONTROL_EXECUTION_RECORD.md`,
`S5_LABEL_ACCESS_AUDIT.md`, `S5_LITERATURE_EVIDENCE_TABLE.md`,
`S5_FAILURE_LOCALIZATION_GRAPH.md`, `S5_RESULT.json`, `S5_FINAL_REPORT.md`.

**NEXT ACTIONS, IN ORDER.** Restore an environment; run the sequence-length and
slot-multiplicity census over `pilot20k_holo_governed_v1/complexes.jsonl` (pure
text read, minutes, calibrates the whole ceiling question); restore
`s2_teacher.py` and `s3s4_observability.py` from `608decf`/`8b7789e`; recover or
re-register the 1,118-complex exposure registry; verify all declared SHA-256;
run the regression suite; verify M-6 empirically on the shards; only then
execute the Stage 1 oracle ceiling under the S5R amendment.


## F-96: S5R-2 executed; pocket-chemistry enrichment observed, ligand margin does not survive de-proxying (2026-08-09)

**OUTCOME.** Terminal verdict `P1B_PAIR_LOCAL_STRUCTURAL_MECHANISM_OBSERVED` on
the registered 128-slot target, with two bounding qualifications. All three F-95
blockers are discharged. Registered at
`research/ssl_b2_structural_observability/PREREG_S5R2_POCKET_CHEMISTRY_ENRICHMENT.md`
**before any arm was scored**. Full report `report/ssl_s5/S5R2_FINAL_REPORT.md`;
machine-readable `report/ssl_s5/S5R2_RESULT.json`; hashes
`report/ssl_s5/S5R2_ARTIFACT_HASHES.json`.

**B-1 DISCHARGED.** Environment available: `conda run -n drug`, Python 3.11.15,
torch 2.6.0+cu124, CUDA 12.4, RTX 4060 Laptop 8 GB. Regression suite: **75
passed**.

**B-2 DISCHARGED WITHOUT DECOMPILATION.** `s2_teacher.py` and
`s3s4_observability.py` were restored with `git checkout 608decf --` and then
*proved* to be the registered implementations: both compile to bytecode that is
value-identical to the surviving `__pycache__` under CPython 3.11. Frozenset
constant ordering is a marshal artefact of the compiling interpreter and is
therefore compared by value, not by `repr`. This is recovery, not reconstruction.

**B-3 DISCHARGED.** The 1,118-complex S4 exposure registry is recovered from
`608decf` into `report/ssl_b2_structural_observability/` and is enumerable from
`dataset/processed/ssl_b2/teacher_dataset.npz` (`pdbs`, `pclus`, `scaffold`).

**PROVENANCE CORRECTION -- the F-95 audit pinned the wrong banks.** Resolved from
bytes rather than prose, the P1B run manifest's four declared hashes are
`supervision_manifest ad0e4803 -> pilot20k_structure_supervision_v2`,
`records 45907b45 -> pilot20k_homology_split_v2`,
`protein_manifest dc087bda -> pilot20k_esm2_t30_slots128_v1`,
`ligand_bank 823815c2 -> pilot20k_mechanism_ligands_v1.pt`; one hop further,
supervision v2 points to `pilot20k_holo_governed_v2` and
`pilot20k_governance_v2`. F-95 recorded the `_v1` banks. **S5R-2 uses the v2
chain throughout.** Every declared SHA-256 in that chain now verifies
(checkpoint, ligand bank, 117/117 supervision shards, 117/117 ESM shards, all
manifests); the checkpoint's stored `val_loss = 1.9781416454571206` matches the
run manifest exactly.

**T-RULE REGISTERED, AND LOAD-BEARING.** `K1` pocket-burial coverage was demoted
to `TAUTOLOGICAL_REFERENCE` before scoring, because `contact = 1[d <= 6]` **is**
the tensor P1B minimised its loss against. `K1` scores R2 = +0.390. Without
T-RULE this run would have reported training fit as mechanism observability.

**THE REGISTERED COORDINATE.** `e[t] = q[t] - b[t]`, the background-subtracted
pocket residue-type composition, with pocket weight `w(s) = sum_a p_contact[a,s]`
from the frozen P1B posterior and composition taken from the sequence under the
identical slot map `min(127, i*128//L)`. **Zero trainable parameters.** Measured:
gauge `max|sum_t e[t]| = 1.77e-16`, bound `max|e| = 0.566`, effective rank 19,
test sd 0.125. Mean enrichment is positive for G,H,Y,F,C,W and negative for
E,A,L,P,V,K -- the textbook binding-site signature, recovered with no
binding-site label.

**M-6 CONFIRMED EMPIRICALLY.** `contact == (distance_bin <= 1)` verified on
**55,162,405** masked pairs across all 117 v2 shards: **zero violations**.

**THE INFORMATION CEILING, PREVIOUSLY `BLOCKED_NO_COMPUTE`.** 94.0 % of
complexes have `L > 128`; occupied slots hold 2.46 residues on average (p95 5,
max 8); **69.1 % of occupied slots are chemically mixed** (1,303,872 of
1,886,412); zero non-standard residues. Against exact-residue enrichment
recomputed from raw coordinates (n=700, 265 homology groups), the *true slot
target itself* reaches only R2 = 0.516 (APOLAR), 0.468 (NEG), 0.395 (ACCEPTOR),
0.304 (AROMATIC), 0.289 (DONOR), 0.179 (POS). **About half the biological
quantity is destroyed by the slot map before any model is involved.**

**LADDER (test split, 1,449 complexes, 359 homology groups, 198 scaffolds;
widest of two clusterings).** CORE = `APOLAR`, fixed in advance. All five
registered conditions PASS: vs mean +0.4467 [+0.2899]; minus ligand-only
+0.1070 [+0.0240]; minus deranged protein +0.2980 [+0.2289]; minus geometry
shuffle +0.4473 [+0.2907]; minus background +0.2587 [+0.1415]. Destructive
controls collapse exactly: geometry shuffle -0.0006, chemistry shuffle +0.0004,
capacity-matched random -0.0023; the oracle returns +1.0000.

**QUALIFICATION 1 -- the ligand margin does not survive de-proxying.** On the
exact-residue target the CORE margin over ligand-only is **+0.078
[-0.003, +0.152]**, lower bound below zero. Ligand-only alone reaches +0.409.
For `AROMATIC` ligand-only is *better* (-0.105). Only `DONOR` (+0.153 [+0.044])
and `NEG` (+0.211 [+0.085]) retain a separable protein-specific increment --
recorded as a **hypothesis generated by this run, not a result**, since `APOLAR`
was the registered core and is not replaced post hoc.

**QUALIFICATION 2 -- pair specificity is weak.** Wrong protein collapses the arm
(`A5 - AD = +0.356 [+0.279]`), so partner specificity is unambiguous. Wrong
ligand with the correct protein still reaches +0.372, so `A5 - AL = +0.114
[+0.064]`: roughly three quarters of the deployable signal is a protein-only
pocket property that does not depend on which ligand is bound. Per-slot
localisation is weak: +0.085 [-0.004, +0.174].

**H-DILUTION REFUTED.** `A5 - A5e = -0.0075 [-0.0452, +0.0408]`: explicit 20-dim
slot composition (zero parameters) and pocket-weighted frozen ESM (fitted
640-dim readout) are indistinguishable. The ESM slot mean is **not** the binding
constraint. Enlarging ESM, fine-tuning it, deepening the GNN and adding
cross-attention now have a *measured* reason to be refused, not only a policy
reason.

**MODULE PARTICIPATION.** Stage 4 synthetic control PASSES: a
1,726,406-parameter head recovers a known function of the frozen pair state, so
optimisation and the objective are not the defect. The un-normalised synthetic
run left the protein and geometry branches `NOT_CAUSALLY_USED` -- a feature-scale
imbalance (activation variance ligand 0.339 vs protein 0.018), exactly the
*normalisation* defect base prereg section 7 names. The registered repair (frozen
per-branch standardisation; constants, no capacity) fixes it: loss
0.0764 -> 0.0037, a 20x reduction, with **all three branches causally used**. On
the real target un-normalised, all three branches are causally used and ablating
the geometry branch raises loss by 35 %. The repair is not uniformly good: on the
real target its 1e-3 standard-deviation floor amplifies near-constant one-hot
atom features and the ligand branch then reads `NOT_CAUSALLY_USED`. Reported as
measured.

**STAGE 5 DECISIVE -- LEARNING BUYS NOTHING.** Aggregated with exactly A5's
denominator, so that A5 is the special case `g(u) = p_contact * C[j,u]`: the
trained 1.73 M head scores +0.0149 against A5's +0.3943 on CORE (delta
**-0.379**), and is worse on every channel. At the justified budget the
**zero-parameter estimator is the better estimator**, and having no parameters it
has no module that could fail to participate.

**SECTION PROBE (structural dry run, NOT a stage advance).** Adapting
coefficients on named mechanism channels, with the update confined to the
numerically identifiable support row space and coverage computed on the
standardised channel block: at k=5 the support rank is **2.84 of a possible 5**
(median condition number 6.5), and the adaptation gain changes sign with k
(+0.064, +0.002, -0.045, +0.035). Verdict
`SECTION_ADAPTATION_WEAK_AND_NON_MONOTONE_ON_STRUCTURAL_CHANNEL`. The first
version of this probe reported 0 % abstention; that was a defect -- the
intercept's constant column dominates the coverage norm -- and it is recorded
here because it would silently vacate any abstention rule built the same way.

**FAILURE LOCALIZATION.** Rejected: insufficient data; missing partner
specificity; objective or optimization defect. Confirmed: underdetermined
estimand (slot proxy versus exact residues); representation failure (2.46
residues per slot, 69.1 % chemically mixed, ceiling R2 about 0.52); **ligand
shortcut, the dominant defect at the CORE channel**; unidentifiable support
section on the structural channel. Not tested: affinity direction. Full graph in
`report/ssl_s5/S5R2_FAILURE_LOCALIZATION.md`.

**LABEL AUDIT.** `recipient_label_reads = 0`, `davis_label_reads = 0`, zero
ChEMBL/BindingDB affinity value reads, zero PLIP or pose-aware reads. Label
firewall intact.

**FROZEN SURFACES.** `theory/FINAL_FROZEN_THEORY/`, `model/`, `contracts/`,
production `scripts/`, CSMO, Band, the positive ridge, the simplex, the fixed
mesh, `K(B(z)F(z))` and production `z` were read only. Writes were confined to
`research/ssl_b2_structural_observability/`, `report/ssl_s5/`,
`report/ssl_b2_structural_observability/` (git-restored),
`dataset/processed/ssl_s5/`, `history.md` and `project_state.json`.

**WHAT THIS AUTHORISES.** Exactly one separately preregistered source-affinity
increment experiment, whose estimand should be the **exact-residue** quantity and
whose primary reported number should be the increment over **ligand-only**. It
does not admit `e` into production `z`, does not authorise DAVIS, recipient
labels, ChEMBL/BindingDB affinity training, PLIP, pose-aware data, P2-P4, or any
change to the frozen operator. Those remain frozen pending separate
authorisation.


## F-97: S5R-2 terminal verdict withdrawn after independent audit; U0-U3/LSMF not present; S7 registered (2026-08-09)

**THIS ENTRY CORRECTS F-96. F-96 IS NOT DELETED.** Its text stands as the record
of what was originally claimed; the claims listed below are withdrawn from it.
Superseded artifacts are registered with their bytes and hashes in
`report/ssl_s5/S5R2_SUPERSEDED_ARTIFACTS.json`; nothing was deleted.

**CORRECTED VERDICT.**

```text
S5_DATA_OR_PREREGISTRATION_CONTRACT_FAIL_CLOSED
SLOT_PROXY_POCKET_CHEMISTRY_SIGNAL_OBSERVED_AS_DEVELOPMENT_EVIDENCE
EXACT_RESIDUE_PAIR_MECHANISM_NOT_IDENTIFIED
AFFINITY_DIRECTION_NOT_TESTED
S6_NOT_AUTHORIZED
```

`P1B_PAIR_LOCAL_STRUCTURAL_MECHANISM_OBSERVED` is **withdrawn**. Three of its
four load-bearing words were wrong: the scored object is a complex-level
aggregate, not pair-local; the statistic is contact localisation annotated with
residue identity, not a mechanism; and the PASS was unit-dependent.

**D-1, DECISION-CHANGING -- the inference unit was not the registered one.** The
base preregistration and `experiment.md` both say the units are *closure
components*. The implementation bootstrapped homology groups and scaffolds
**separately** and reported the wider interval. Recomputed on the CORE channel
from retained arrays, with no fitting:

| unit | n | cond 1 | cond 3 | A5-AL |
|---|---:|---|---|---|
| homology | 359 | +0.447 [+0.304, +0.578] | +0.298 [+0.229, +0.371] | +0.134 [+0.076, +0.203] |
| scaffold | 198 | +0.447 [+0.290, +0.597] | +0.298 [+0.237, +0.353] | +0.134 [+0.075, +0.216] |
| cells | 546 | +0.447 [+0.310, +0.577] | +0.298 [+0.235, +0.362] | +0.134 [+0.081, +0.196] |
| **closure components** | **54** | **+0.447 [-0.020, +0.497]** | **+0.298 [-0.076, +0.354]** | **+0.134 [-0.161, +0.223]** |

The closure graph has **54** components and the largest holds **1,291 of 1,449 =
89.1 %** of the test split. **Under the registered unit no condition retains a
positive lower bound.**

**D-2, FACTUALLY WRONG -- F-96's module-participation numbers came from a smoke
test.** F-96 reported the un-normalised synthetic run as `0.0760 -> 0.0383` with
the protein and geometry branches `NOT_CAUSALLY_USED`. Those are a
**32-train/16-test/1-epoch smoke test** whose JSON was read before the registered
run overwrote the same filename at 12:06:58. The authoritative artifact
`S5R2_HEAD_SYNTHETIC.json` and the console log `head_syn.log` both record
`train=1200 test=400`, **loss 0.077279 -> 0.002167, all three branches
CAUSALLY_USED**. Consequently the "normalisation defect" never existed at the
registered scale, and the repair applied in response **degraded every measured
outcome**: synthetic 0.00373 vs 0.00217, real 0.00590 vs 0.00301, and it created
the only `NOT_CAUSALLY_USED` verdict in the record via a 1e-3 standard-deviation
floor amplifying near-constant one-hot atom features up to 1000x.

**D-3 -- the head-versus-A5 conclusion is withdrawn, not corrected.** It exists
only in the superseded normalised run, was computed on a complex-level aggregate
the head never optimised, and used calibration fitted on **test** while the
ladder fitted on **train**, so `+0.4467` and `+0.3943` were never comparable. The
better un-normalised head was never compared to A5 at all. "Learning buys
nothing" is not supported.

**D-4 -- wrong-ligand arm was not masked.** `lig_valid` was never intersected, so
67 of 1,449 test rows entered arm `AL` as exact zero vectors. Corrected on 1,382
rows: `A5 - AL = +0.1383 [+0.0827, +0.2037]` versus the published
`+0.1337 [+0.0750, +0.2164]`. Direction unchanged. **`s5r2_ladder.py` is
deliberately left unedited** so its recorded hash still matches what produced
`S5R2_LADDER.json`; the fix is mandatory in S7.

**D-5 -- the exact-residue analysis is post hoc and unreproducible.**
`s5r2_exact_eval.py` (11:53:46) postdates `S5R2_LADDER.json` (11:52:37), so the
de-proxying was designed after the ladder was scored. Neither the 700 `e_exact`
vectors nor the sampled indices were persisted, so the headline
`+0.078 [-0.003, +0.152]` is **not reproducible from retained artifacts**. It
weakened the author's own claim, which is the benign direction, but it is
exploratory evidence.

**D-6/D-7 -- controls looser than stated.** The wrong-ligand policy required only
a different CCD hash and matched atom count, and **2,303 of 14,793 (15.6 %)**
wrong ligands share the correct ligand's Bemis-Murcko scaffold. Partner reuse
reached 18x (protein) and 58x (ligand), unmodelled by any bootstrap. Pairwise
`< 40 %` identity was never verified directly, only inferred from MMseqs2 cluster
membership.

**D-9 -- the evaluation panel is development-exposed.**
`p1b_gate_pilot20k_seed17_v4/gate_report.json` records `split = test`, 1,477 of
1,490 scored, `gate_status = PASS`. The S5R-2 evaluation split **is** the split
the P1B Gate was decided on. It is checkpoint-held-out but
research-development-exposed, and is **not** a confirmation set. Base
preregistration section 3.3 remains **unmet**.

**D-10 -- preregistration chronology is prose, not proof.** `PREREG_S5R2` and all
of `report/ssl_s5/` are untracked; `git log --all` returns nothing for them. The
only evidence is one mutable mtime, and the prereg's last edit (11:26:05)
postdates a census (11:24:05) that reported per-channel retention in which APOLAR
ranks highest. The entire S5R-2 result is therefore `DEVELOPMENT_EVIDENCE_ONLY`.

**D-11 -- wording.** "Proved to be the registered implementations" overstates.
Bytecode equivalence gives **executable-semantic equivalence under CPython
3.11**; `git checkout 608decf` is the actual source-identity evidence and the
report inverted the emphasis.

**WHAT SURVIVES.** The P1B input SHA-256 chain verifies end to end; the v1->v2
provenance correction is right; `contact == (distance_bin <= 1)` holds on
55,162,405 pairs with zero violations; the gauge is exact to 1.77e-16 and the
bound is 0.566; the census (94.0 % `L>128`, 2.46 residues per occupied slot,
69.1 % chemically mixed) reproduces; destructive controls sit at zero
(-0.0006, +0.0004, -0.0023) and the oracle returns exactly 1.0000; the label
firewall is intact; regression passes. **A real, non-artifactual structural
association exists -- as development evidence, at slot resolution, on an exposed
panel, under a non-registered unit.**

**U0-U3 AND LSMF -- NOT PRESENT, THEREFORE NOT REPRODUCED.** An independent
search found **zero** occurrences of `LSMF`, `sparse mechanism field`,
`anchored satellite`, `satellite component` or `u_c(h_i)` anywhere in the tree or
in `git log --all`. Every claimed statistic (549 clusters, 4,387 documents, 336
union components, 64.6 % giant, 71,427 rows, 2,587 tasks, 325 satellite
components) is absent; the apparent numeric hits are coincidental substrings
inside floating-point values, e.g. `71427` inside `0.009371427498929866`. No
preregistration, code, manifest, split, prediction file or label audit exists.
The tokens `U1/U2/U3` **do** occur here, but they denote the apo/holo ensemble
hierarchy of `PREREG_G0_DATA_FEASIBILITY.md`; `U0` does not occur at all. This is
a token collision, not partial corroboration. The only BindingDB material on disk
is the raw XP3/XP4 acquisition zips, which were not opened. Full record:
`report/ssl_s5/U0_U3_AND_LSMF_REPRODUCIBILITY_AUDIT.json`.

Noted without weight: the external summary reports a 64.6 % giant component and
this audit independently measured 89.1 %. Giant-component pseudoreplication under
simultaneous protein and ligand closure is evidently real in this domain, but the
local measurement stands on its own and the external numbers add nothing.

**BOUNDARY ADJUDICATED: A and B jointly, A dominant.** The data *as currently
closed* cannot identify a population protein-ligand interaction -- 54 components
with an 89.1 % giant is not a panel that resolves 0.02 at 95 %, and XP3/XP4
already recorded the public-panel noise floor. But a localised **exact-residue**
mechanism remains a justified hypothesis, because this repository measured that
the 128-slot map destroys about half the biological quantity before any model
runs. **C is rejected** -- that is precisely the withdrawn claim. **D is
rejected** -- the P1B gate's correct-minus-deranged-protein contact AUPRC gap of
+0.353 [+0.341, +0.365] shows protein-specific information does reach the
interface. Full reasoning:
`report/ssl_s5/BOUNDARY_ADJUDICATION_AND_EVIDENCE_COMPARISON.md`.

**PANEL FEASIBILITY -- the decisive practical finding.** The S1b acquisition
registered 15,003 RCSB candidates released >= 2024-01-01, all disjoint from the
10,468 pilot20k-exposed IDs. 1,118 were S4-scored and are exposed; 358 were
downloaded and rejected; **13,885 were never downloaded and never scored**. A
genuinely untouched confirmation panel is therefore available with no new
governance question -- same source, same CC0-1.0 licence, same release rule, same
acquisition script.

**S7 REGISTERED, NOT AUTHORIZED TO EXECUTE.**
`research/ssl_b2_structural_observability/PREREG_S7_EXACT_RESIDUE_LOCAL_MECHANISM.md`.
Exact-residue target with the 128-slot readout demoted to a control arm; panel
drawn score-blind from the 13,885 untouched IDs (seed 20260810, 2,500 records);
**union closure components as the primary inference unit** with a fail-closed
precondition if the largest component exceeds 25 % or fewer than 60 components
exist; wrong-ligand control tightened to require a different Bemis-Murcko
scaffold; all arm masks intersected; per-complex predictions, per-unit SSE and
target vectors persisted before any metric is reported; smoke tests forced to a
`_SMOKE` filename so they can never overwrite a registered artifact. Two **new**
conditions make S7 decisive where S5R-2 was not: `B5 - BL` (pair specificity as a
condition) and `B5 - B4` (the learner must beat the exact-residue zero-parameter
estimator). Thresholds unchanged at 0.02 with LCB > 0.

**S7 may not execute until this preregistration is committed and its SHA-256
recorded in `project_state.json`.** That requirement exists because it is exactly
what S5R-2 lacked.

**LABEL AUDIT.** `recipient_label_reads = 0`, `davis_label_reads = 0`, zero
ChEMBL/BindingDB affinity value reads, zero PLIP or pose-aware reads. Verified
independently: no S5R-2 or audit script opens any protected path; keyword matches
are docstrings only.

**FROZEN SURFACES.** `theory/FINAL_FROZEN_THEORY/`, `model/`, `contracts/`,
production `scripts/`, `weights/`, `config/`, CSMO, Band, positive ridge, the
simplex, the fixed mesh, `K(B(z)F(z))` and production `z` are unmodified --
`git status --porcelain` over all of them is empty. Regression: 75 passed.

**NOT AUTHORIZED.** S6; any affinity read; any training; LSMF training; admission
of any biological coordinate to `z`; any push. Nothing has been committed or
pushed.


## F-98: S7_L2B fails closed at R0 -- the five indexed L2 artifacts are not materialized (2026-08-09)

**TERMINAL VERDICT: `FIVE_ARTIFACTS_NOT_MATERIALIZED`.** Earliest applicable in
the registered hierarchy; R1-R6 were not entered; no GPU work was performed.

**GOVERNING DOCUMENT.** `research/Residue locator/METASIEVE_L2B_FIVE_ARTIFACTS_CONSOLIDATED_REPORT_2026-08-09.md`
(SHA-256 `16e6c62319750b54bcbc054d826fce3712c41473cabfc425fef9cb1c64ad0326`) is
present and was read in full. It is a research synthesis and index. Under the
governing rule its numerical claims become repository evidence only once the
corresponding machine-readable artifacts are located or rebuilt; none were.

**ALL FIVE INDEXED OBJECTS ARE ABSENT** from the working tree and from
`git log --all`. `research/Residue locator/` contains exactly one file -- the
consolidated report itself.

| # | Artifact | Status | Blocks R0 |
|---|---|---|---|
| 1 | `METASIEVE_L2B_NEXT_STAGE_ANALYSIS_2026-08-09.md` | ORIGINAL_UNAVAILABLE | no |
| 2 | `L2A_ORACLE_FACTOR_BUDGET.json` | ORIGINAL_UNAVAILABLE | **yes** |
| 3 | `L2B_INDEPENDENT_DATA_GATE.json` | ORIGINAL_UNAVAILABLE | **yes** |
| 4 | `L2B_RESIDUE_STUDENT_FEASIBILITY.json` | ORIGINAL_UNAVAILABLE | **yes** |
| 5 | `PREREG_L2B_PLM_LOCALIZER.md` | ORIGINAL_UNAVAILABLE | no |

**THE INPUTS ARE ABSENT TOO, SO NO FORMAL REPLACEMENT WAS POSSIBLE.** Searches
for `monn`, `MONN`, `plip`, `PLIP`, `residue_atom`, `anchor`, `satellite` return
**nothing**. There are no residue-atom interaction labels of any kind in this
repository. CD-HIT-2D, the tool the claimed 40 % closure used, is **not
installed** (the project standard is MMseqs2; substituting it would change the
closure definition and therefore the experiment). The single `l2_` filename hit
in the whole tree is the PDB code `5el2` inside PLINDER metadata -- a
coincidence, not evidence.

**No look-alike rebuild was created**, and no substitute corpus was used to
simulate an L2A or L2B result.

**STRUCTURAL BLOCKER THAT DATA ACQUISITION WOULD NOT FIX.** Arm `B4` is defined
as the frozen incumbent exact-residue localizer from the prior low-capacity L2
experiment. **No such checkpoint exists.** The only model checkpoint in the
repository is `report/mechanism_refactor/p1b_pilot20k_seed17_v1/best.pt`, the
P1B contact/distance bridge over a **128-slot** protein index -- a different
estimand at a different protein resolution, and not an exact-residue localizer.
The primary Gate `B5 - B4 >= 0.02 AP` is therefore **unsatisfiable by
construction**. Even a successful MONN acquisition would leave S7_L2B blocked
until an incumbent `B4` is separately built and frozen under its own
preregistration. **S7_L2B is blocked at two independent points and only one of
them is a data problem.**

**THE GOVERNING DOCUMENT'S RUNTIME CLAIM IS FALSE IN THIS ENVIRONMENT.** It
states the runtime "lacks PyTorch, the ESM software and weights, and a GPU".
Measured: **torch 2.6.0+cu124 present, CUDA available, 7.44 GB GPU free,
`transformers` present**. Genuinely absent: the `esm` package and the
`esm2_t33_650M_UR50D` weights (only `esm2_t30_150M_UR50D` is cached -- the model
P1B used). The runtime gap is real but far narrower than stated, and **it is not
the binding constraint**. `PLM_RUNTIME_NOT_FEASIBLE` is **not** claimed: the
650M contract was never attempted, so asserting infeasibility would be a
fabricated verdict.

**EVERY NUMBER IN THE GOVERNING DOCUMENT REMAINS AN EXTERNAL CLAIM**, including
frozen L2 pair AP `0.01786`; the oracle factor table and
`oracle residue - full = 0.23476 [0.22393, 0.24566]`; the ~33x residue-over-atom
ratio; `72,226` same-protein pairs with Jaccard `0.47265` / overlap `0.75766`;
the sequence student `0.08664` vs `0.06798` and its `0.01866` shortfall against
the `0.02` threshold; `4,067`/`701`; `8,646`/`1,328`/`18,871`/`199,015`;
`524`/`157`/`251`/`317`/`2,555`; `2,404`/`1,297,939`/`7,119`; atom-level AP
`0.6767`; and all reported SHA-256, syntax-check and JSON-parse confirmations.
This does not assert they are wrong -- it asserts they are unverifiable here,
which under the evidence-precedence rule is the same as not being evidence.

The withdrawn S5R2, U0-U3 and LSMF claims were **not** reinstated and were not
used anywhere in this step.

**STAGE-NAMING RECONCILIATION.** The instruction requires one stage and forbids
parallel S7/L2B protocols. `PREREG_S7_EXACT_RESIDUE_LOCAL_MECHANISM.md`
(SHA-256 `497ffda5...`) is marked **`SUPERSEDED_BY_STAGE_NAMING`, on hold, not
deleted**. It is a *different* experiment -- exact-residue pocket composition
from RCSB coordinates, versus S7_L2B's residue-atom interaction localization
from PLIP labels. Note for planning: **if the MONN acquisition is not
authorized, the held S7 is the only structural stage this repository can
actually run**, because its data (the 13,885 never-downloaded, never-scored RCSB
candidates) is already present and already governed.

**WHAT R0 WOULD REQUIRE.** (1) authorize and acquire the MONN additional-PDB
data at commit `f2b62ccf...` with a tier/licence/redistribution decision in
`DATASET_ROLE_REGISTRY.json`; (2) rebuild the residue-atom edge corpus under new
names and hashes, **not** reusing the historical `4,067/701` metric identity;
(3) install CD-HIT-2D or formally register a change of closure tool; (4) build
and freeze an incumbent `B4` under its own preregistration; (5) acquire the
`esm` package and 650M weights (~2.5 GB) with revision, licence, URL and weight
SHA-256 recorded; (6) only then write and freeze the single unified S7_L2B
preregistration. Steps 1, 3, 4 and 5 each require explicit authorization; none
was taken.

**GPU AUTHORIZATION CHECKLIST: NOT_AUTHORIZED**, failing at its first item.
No smoke test, no extraction, no training, no scoring, no download, no install.

**LABEL AUDIT.** `recipient_label_reads = 0`, `davis_label_reads = 0`, zero
ChEMBL/BindingDB affinity reads, zero PLIP or pose-aware reads. Firewall intact.

**FROZEN SURFACES.** `theory/`, `model/`, `contracts/`, `scripts/`, `weights/`,
`config/` unmodified. `K(B(z)F(z))`, CSMO, Band, simplex, positive ridge, fixed
mesh and production `z` untouched. The trace-set/minimax structure and the
frozen CDF-band theory were kept separate; no identification between the
trace-set centre/radius and `K(BF)` was made or implied. Nothing committed,
nothing pushed.


## F-99: S7_L2B_R0R -- MONN raw layer reproduced; closure topology fails closed (2026-08-09)

```text
R0R-1 MONN raw provenance reproduction ... PASS
R0R-2 new ligand identity and closure .... FAIL CLOSED
R0R-3..R0R-6 ............................. NOT REACHED
TERMINAL: NEW_CLOSURE_TOPOLOGY_INSUFFICIENT
```

Stopped at the earliest failed boundary. No closure was relaxed to continue, no
`B4` was fabricated, no GPU work was performed, no affinity table was parsed.

**R0R-1 PASS.** MONN cloned and pinned to `f2b62ccf49c18a9502aa0eb0d582c6e0735ef200`
(HEAD verified). Licence recorded verbatim: *"The algorithm and data can be used
only for NON COMMERCIAL purposes."* Declared non-commercial research use; the
clone lives under `dataset/raw/monn/`, gitignored, so **no MONN byte is committed**.

`rebuild_monn_edge_corpus.py --strict-hashes` reproduced every target exactly:
development 12,987 raw / 12,738 mapped / 195,798 binary / 202,766 typed;
additional PDB 1,853 / 1,851 / 9,832 / 9,832; missing atom references 0 in both.
The deterministic outputs are **byte-identical to the supplied
`VERIFIED_RAW_AUDIT.json` hashes** (`9489540b...`, `cbb0ed98...`) and
**bit-identical across two independent runs**.

**ACQUISITION DEFECT FOUND AND REPAIRED -- would have been silent.** All six
source files failed SHA-256 on first clone, each *larger* than the manifest size:
`out7` 58,935,197 vs 53,017,418; `independent_dataset_interaction_dict`
4,719,169 vs 4,234,610; `mol_dict` 59,934,686 vs 59,740,405;
`independent_dataset_mol_dict` 2,267,662 vs 2,259,203. Cause:
`core.autocrlf=true` in the inherited Windows Git configuration, and MONN ships
no `.gitattributes` marking its pickles binary, so Git applied LF->CRLF to
binary pickles. Repaired with `core.autocrlf=false`, `core.eol=lf`,
`git rm --cached -r .`, `git reset --hard`. **Without the supplied manifest
hashes this would have produced a subtly corrupted corpus undetected.** Any
future acquisition of a binary-bearing repository on this machine must disable
autocrlf before checkout. This is a property of the acquisition environment, not
of MONN.

**R0R-2 ligand identities: complete.** All **10,972** required CCD codes resolved,
sanitized, canonicalized and scaffolded. **Zero failures, zero quarantined
complexes.** Frozen policies: RDKit 2023.09.6; `pickle.load(encoding="bytes")` --
`latin1` and the default both fail with *"Bad pickle format: bad endian ID"*, so
this is the only working policy; exact-graph identity =
`sha256(MolToSmiles(isomericSmiles=False, canonical=True))`, collapsing
stereoisomers, the conservative direction; an **empty** Murcko scaffold generates
**no** closure edge, since treating `""` as shared would merge every acyclic
ligand into one artificial blob. Homology used MMseqs2 strictly as a candidate
generator (5,552 pairs) with parasail Smith-Waterman as the authority
(BLOSUM62, gap_open 10, gap_extend 1, identity = matches / local alignment
length, coverage = alignment length / min(len)); at 40 %/80 %, **1,883 accepted
and 3,669 rejected** -- the candidate generator alone would have over-merged
threefold.

**R0R-2 FAIL CLOSED -- the topology.** 14,589 complexes -> **360** union
components, largest **13,595 = 93.19 %**, median component size 1, only 42
components with >=5 complexes, and **18 components straddle development and
additional-PDB**, so the additional-PDB set is not an independent confirmation
cohort even before publication/time closure. A paired whole-component bootstrap
over 360 units where one holds 93 % of the data has effective sample size near
one, while an equal-weight component macro-average is simultaneously dominated by
singletons carrying a handful of edges.

**CAUSE LOCALIZED BY ABLATION -- it is not protein homology.**

| relations | components | largest | fraction |
|---|---:|---:|---:|
| **protein only** (PDB+seq+UniProt+40 % homology) | **1,994** | 453 | **3.11 %** |
| protein + exact ligand graph | 722 | 11,120 | 76.22 % |
| protein + scaffold | 556 | 13,185 | 90.38 % |
| **frozen full rule set** | **360** | **13,595** | **93.19 %** |
| homology only | 13,736 | 22 | 0.15 % |

Protein closure alone partitions beautifully. Adding **exact ligand graph** alone
drives the giant to 76 %. The mechanism is chemical: promiscuous ligands and
cofactors (ATP/ADP/NAD/heme-like) recur across unrelated protein families, so
"same exact ligand graph" transitively bridges the protein universe; shared ring
systems compound it. **A closure that simultaneously closes protein identity and
ligand identity over a promiscuous-cofactor corpus does not partition.** This is
a property of the estimand, not an implementation artifact.

**ONE HISTORICAL NUMBER REPRODUCED.** Development exact sequences = **2,404**,
matching the consolidated report. That single claim graduates to *reproduced*.
`4,067/701`, `8,646`, `524/157` and every AP value remain unverified external
claims and were not used anywhere.

**CORPUS CHARACTERISTICS RECORDED.** All seven PLIP channels clear a 5 %
development-prevalence bar and are evaluable: Hydrogen Bonds 94.0 %, Hydrophobic
86.5 %, Water Bridges 53.5 %, Salt Bridges 33.7 %, pi-Stacking 31.0 %, pi-Cation
10.4 %, Halogen Bonds 7.2 %. The development complete residue x atom matrix is
**403,454,851** cells at a positive rate of **4.85e-4** (~1 in 2,060).

**THE ADJUDICATION THIS REQUIRES (not taken).** Exactly one decision could
unblock R0R-2, and it is a scientific decision about the estimand: define the
inference partition on the **protein side alone** and control ligand leakage by
evaluation design -- held-out ligand graphs plus the `BX` wrong-ligand control --
rather than by closure. Measured projection: development 12,738 complexes /
1,669 components / largest 3.46 % / 458 components >=5; confirmation candidate
(protein-disjoint from development) 710 complexes / 325 components / largest
1.97 % / 348 sequences / 346 ligand graphs / 3,383 positive edges. That candidate
passes **6 of 7** capacity checks and fails `components_ge5_complexes`
(**37** against a required 60). So the alternative is **not a free pass** either.
It would also change the claim: the experiment would generalize over *proteins*,
not over *protein-ligand pairs*, and that must be written into the estimand
rather than discovered afterwards.

**WHAT MUST NOT BE CONCLUDED.** This is outcome **2 -- closure cannot support
inference**. It is not outcome 5. Nothing here says sequence-plus-2D inputs lack
the required information. No model was trained, no AP was computed, no `B4`
exists, no `B5` Gate was attempted, and the previous session's probe showed the
compute path is largely available.

**ARTIFACTS.** `report/s7_l2b_r0r/` -- `R0R_RECONSTRUCTION_REPORT.md`,
`MONN_LOCAL_REPRODUCTION_AUDIT.json`, `NEW_EDGE_CORPUS_MANIFEST.json`,
`NEW_CLOSURE_AND_SPLIT_MANIFEST.json`, `CLOSURE_RELATION_ABLATION.json`,
`PARTITION_FEASIBILITY_PROJECTION.json`, console logs. Code under
`research/s7_l2b_r0r/`. `PUBLICATION_TIME_CLOSURE_AUDIT.json`,
`PREREG_B4_EXACT_RESIDUE_BASELINE.md`, `B4_FROZEN_ARTIFACT_MANIFEST.json`,
`PREREG_S7_L2B_UNIFIED.md` and `GPU_PREFLIGHT_AUDIT.json` were **not** created,
because their steps were never reached and writing them would imply work that
did not happen.

**FROZEN SURFACES.** `theory/`, `model/`, `contracts/`, `scripts/`, `weights/`,
`config/` unmodified. `K(B(z)F(z))`, CSMO, Band, simplex, positive ridge, fixed
mesh and production `z` untouched. No affinity, DAVIS, recipient or few-shot
support label read. Nothing committed, nothing pushed. The pre-existing dirty
worktree was preserved; no unrelated change was reverted or deleted.


## F-100: S7_L2B development gate -- protein signal identified, below preregistered effect size (2026-08-09)

```text
Data / provenance ........... PASS (R0R-1)
Closure construction ........ RESOLVED (protein partition + ligand disjointness filter)
Evaluator contract .......... PASS
Trainability control ........ PASS
Matched baseline B4 ......... TRAINED AND FROZEN
Registered Gate outcome ..... FAIL CLOSED (1 of 5 gates met the effect size)
Escalation rule ............. FIRED -- B5 authorised by the preregistration
B5 frozen ESM2-650M ......... NOT RUN, blocked on weight acquisition
Confirmation cohort ......... SEALED, never opened
```

**PREREGISTRATION COMMITTED BEFORE ANY MODEL EXISTED.**
`research/s7_l2b_r0r/PREREG_S7_L2B_UNIFIED.md`, SHA-256
`2c333f223ae450c566cc62b1a3b276ff59c065c38348005ad9504ac1930b9a92`, commit
`ce186f4`. This is the F-97 lesson applied: chronology is now cryptographic, not
prose. Only the preregistration was committed; nothing was pushed and the dirty
worktree was preserved.

**THE CLOSURE CORRECTION THAT MADE THIS MEASURABLE.** R0R-2 measured that
union-MERGING ligand identity into the inference partition yields a 93.19 % giant
component. Ablation localised the cause to the ligand side (protein closure alone
gives 1,994 components with a 3.11 % largest; adding exact ligand graph alone
drives it to 76.22 %). The correction is to enforce ligand closure as a
**disjointness filter between train and held-out** rather than as a merge
relation. This is *stricter* between the two sets than a merge, while leaving the
partition usable: train 9,758 complexes / 151,065 positives; held-out A 2,415
complexes / 196 components / largest 15.2 % / 36,046 positives; held-out B
(scaffold-strict) 1,881 / 160 / 16.5 % / 29,073. Claim scope stated up front:
generalisation over **proteins**.

**CONTRACT CHECKS.** MONN `atom_name` lists deposited atoms INCLUDING hydrogens
while `mol_dict` is heavy-atom only, and `atom_idx` is a plain identity range,
NOT a mapping into the molecule. The correct mapping is the rank of a slot among
heavy positions. Validated per record: **14,586 of 14,589 pass**, 3 quarantined
on heavy-count mismatch, **zero** positive edges on hydrogens, **zero** residue
indices out of range. Evaluator self-test: AP in float64 with ties broken by a
fixed `(residue_index, atom_slot)` order is bit-identical under a random
permutation of the flattened matrix, returns exactly 1.0 when all rows are
positive, and excludes complexes with no positives. Trainability control: the
identical pipeline recovers a KNOWN function of the frozen inputs at macro-AP
**0.7588** against prevalence 0.0081, so optimisation and objective are **not**
the defect.

**RESULT, held-out A, complete-matrix AP, 196 components.**

| arm | macro-AP |
|---|---:|
| B0 prevalence | 0.00250 |
| BL ligand-only | 0.00450 |
| BP wrong protein | 0.00485 |
| BM motif shuffle | 0.00530 |
| BX wrong ligand | 0.00768 |
| **B4 non-PLM residue** | **0.02295** |

| gate | delta | LCB95 | threshold | |
|---|---:|---:|---:|---|
| G1 B4-B0 | +0.02045 | +0.01719 | 0.02 | PASS |
| G2 B4-BL | +0.01845 | +0.01523 | 0.02 | FAIL |
| G3 B4-BP | +0.01810 | +0.01487 | 0.02 | FAIL |
| G4 B4-BM | +0.01765 | +0.01408 | 0.02 | FAIL |
| G5 B4-BX | +0.01527 | +0.01238 | 0.02 | FAIL |

**Registered outcome: FAIL CLOSED.** Four of five contrasts fall below the
preregistered 0.02 absolute-AP effect size. The threshold was frozen before any
model existed and is **not** relaxed; statistical significance does not override
a preregistered practical-effect requirement. Held-out B reproduces the ordering
(B4 0.02153 vs BL 0.00452, delta +0.01701 [LCB +0.01399]) with **no sign
reversal**, so the result is not an artifact of ligand-scaffold overlap.

**EVERY CONTRAST IS DIRECTIONALLY POSITIVE WITH LCB ABOVE ZERO.** Substituting a
wrong protein collapses B4 from 0.02295 to **0.00485**, essentially back to
ligand-only: almost all of B4's advantage over ligand-only requires the correct
protein. Motif shuffle collapses it to 0.00530, so residue order and position are
load-bearing. Wrong ligand collapses it to 0.00768, well above ligand-only but
far below B4. Ligand-only itself sits at 0.00450 against a prevalence of 0.00250.
The constraint-2 identifiability requirement is met in DIRECTION on every
control; it is the effect SIZE that is not met.

**AUTO-ASSIGNED LABEL WITHDRAWN.** The runner emitted
`S7L2B_LIGAND_ONLY_SHORTCUT` from a decision-logic defect that mapped "G2 did not
pass" straight to the ligand-shortcut verdict. G2 can fail on effect size while
being strongly directional, which is what happened. The label is contradicted by
its own data: B4 is **5.1x** ligand-only and ligand-only barely clears
prevalence. Corrected to
`S7L2B_PROTEIN_SIGNAL_IDENTIFIED_BELOW_PREREGISTERED_EFFECT_SIZE`. A descriptive
label does **not** convert the outcome into a pass; the Gate outcome remains FAIL
CLOSED. The raw run output is retained unedited and the adjudication is in
`S7L2B_DEVELOPMENT_GATE_ADJUDICATED.json`.

**FAILURE LOCALIZATION.** Data: not the cause (six source hashes verified, corpus
bit-identical across runs). Closure: not the cause (196 components, largest
15.2 %, ligand-disjoint from train). Optimization: not the cause (trainability
control at 0.759). Identifiability support: not the cause (all five LCBs above
zero). **Representation: THE CAUSE** -- explicit non-PLM residue features carry
real but weak protein information. Transfer: not tested, confirmation sealed.

**ESCALATION RULE FIRED, AND IT IS EVIDENCE-BACKED.** PREREG section 7 authorises
B5 only if B4-BL fails G2 OR B4 macro-AP stays below 0.10. Both hold. This is
exactly what constraint 9 requires before adding a larger encoder: prior
experiments must explicitly indicate missing information. That evidence now
exists and is measured, and the residue representation is the only component B5
changes. B5 must still clear G1-G5 plus G6 (B5-B4 >= 0.02, LCB > 0).

**B5 NOT RUN -- acquisition blocked.** The frozen `esm2_t33_650M_UR50D` weights
(2,609,621,831 bytes) could not be acquired. The endpoint answers HEAD 200 and
throughput reaches 3.5-6.6 MB/s while flowing, but the transfer drops after
several hundred MB to ~2.1 GB and the client begins a NEW `.incomplete` blob from
zero instead of range-resuming, so progress is not cumulative.
`PLM_RUNTIME_NOT_FEASIBLE` is **NOT** claimed: the GPU contract was never
exercised and the GPU has 7.44 GB free, ample for a frozen 650M forward at window
1000, batch 1. Nothing here is evidence about the PLM hypothesis.

**SELF-CORRECTION RECORDED.** I twice declared the download stalled by polling
blob size with `Get-ChildItem`; on Windows the reported length of an actively
written file lags, and I killed one attempt at 2157 MB of 2610 MB on that false
reading. The correct method is to wait on process exit. This is recorded because
the same mistake would corrupt any future judgement about acquisition
feasibility.

**ARTIFACTS.** `report/s7_l2b_r0r/` -- `S7L2B_DEVELOPMENT_REPORT.md`,
`S7L2B_DEVELOPMENT_GATE.json` (raw), `S7L2B_DEVELOPMENT_GATE_ADJUDICATED.json`,
`GPU_PREFLIGHT_AUDIT.json`, `SPLIT_CONSTRUCTION_CENSUS.json`, plus the R0R-1/2
artifacts. Sealed per-complex predictions for all six arms and the B4 checkpoint
with SHA-256 are under `dataset/processed/s7_l2b_r0r/preds/`. Code under
`research/s7_l2b_r0r/`.

**BOUNDARIES HELD.** No affinity, DAVIS, recipient or few-shot support label was
read. `theory/`, `model/`, `contracts/`, production `scripts/`, `weights/`,
`config/`, CSMO, Band, simplex, positive ridge, mesh and production `z` are
unmodified. All new code is under `research/`. The confirmation cohort was never
scored, and R0R-3 publication/time closure remains unbuilt, so it could not be
opened even on a development pass. Few-shot adaptation remains closed and no
biological statistic is admitted to `z`.


## F-101: I-1/I-2 integrity repairs -- atom correspondence verified; recoverable structure is the RESIDUE MARGINAL, coupling beyond margins is weak (2026-08-09)

Label-side only. No training, no GPU, no affinity read. Held-out A: 2,409
complexes with positives in 196 protein closure components; 2,314 entered the
coupling test.

**I-1 ATOM CORRESPONDENCE -- BLOCKER DISCHARGED.** The mapping under test was
`atom_slot -> rank among non-hydrogen positions in atom_names`, with the RDKit
molecule's own element symbols as the authority. Two earlier necessary
conditions (heavy count equals `mol.GetNumAtoms()`; no positive edge on a
hydrogen) do not establish that the ORDER corresponds, so order was tested
directly at 375,311 positions.

My first two parsers were wrong, and both errors are recorded. A strict
two-letter element parser flagged 166 records, every one of which was the
parser over-matching a single-element name with a positional suffix -- `CL1` is
carbon, `PD` is phosphorus, `SB2` is sulfur. A one-or-two-letter compatibility
rule left 3, of which two were gaps in my symbol list (`OS` osmium, `PR`
praseodymium). The correct, list-free rule is the actual PDB convention: **a
ligand atom name BEGINS with its element symbol**. Under it exactly **one**
record is genuinely incompatible: `5w8v` / CCD `9YP`, position 19, name `CAJ`
aligned to a nitrogen.

Result: **14,585 records admitted; 4 enumerated and quarantined**
(`4xe1`, `5W31_MBO`, `5w8v`, `6GBR_MBO`; sha256 `a64e071d...`). The quarantine
is wired into `s7_dataset.build()`, so the exclusion is part of the data
contract rather than a manual step. Verdict
`ATOM_CORRESPONDENCE_VERIFIED_WITH_ENUMERATED_QUARANTINE`.

**I-2 ANOVA MACHINERY VALIDATED, AND IT IS NOT DOUBLE CENTERING.** The registered
coupling object `G = mu + alpha_i + beta_j + C`, `C = G - Proj_W(additive)`, is
solved as a genuine weighted least-squares ANOVA on the actual mask by weighted
ALS with an identified re-centring step. Self-test both directions: with a
complete mask and uniform weights it reproduces classical double centering to
**2.78e-16**; with non-uniform weights it **differs** by up to **0.423**. Naive
double centering is used only as that oracle and never on real data.

**I-2 THE RECOVERABLE STRUCTURE IS THE RESIDUE MARGINAL.** Macro-AP over
components on the complete residue x heavy-atom matrix:

| object | macro-AP |
|---|---:|
| ligand-only (measured) | 0.00450 |
| Oracle-A true atom marginal | 0.00850 |
| B4 (measured) | 0.02295 |
| **Oracle-R true residue marginal** | **0.21633** |
| **Oracle-RA additive projection** | **0.39075** |
| exact pair (sanity) | 1.00000 |

Residue localisation is worth about **25x** atom propensity. **B4 recovers only
about 6 % of the oracle-marginal ceiling**; headroom to Oracle-RA is **+0.368 AP**
and to Oracle-R alone **+0.193 AP**. The bottleneck is predicting *which residues
bind at all*, not *which residue pairs with which atom*.

This is a fresh measurement on our own corpus, split, evaluator and code. It is
**not** a reproduction of the unverified consolidated-report figures
(0.25262 / 0.02494), which remain external claims; the direction agrees, the
numbers are ours and differ.

**I-2 COUPLING BEYOND MARGINS IS WEAK IN THE LABELS.** Statistic: leading
singular-value share of the marginal-orthogonal residual on the active
submatrix. Null: degree-preserving bipartite rewiring by checkerboard swaps, 20
per complex, 30x(positives) swap attempts, holding every `d_i` and `e_j` exactly
fixed. True mean **0.6124** vs null mean **0.5921**, difference **+0.0203**,
median z **+0.41**, and only **63.1 %** of complexes exceed their own null
(chance 50 %). A real but small excess: the typical complex sits within half a
standard deviation of a degree-matched random matrix. Rewiring was used strictly
as an evaluation control and never as a training negative.

**FALSIFIABLE PREDICTION REGISTERED IN ADVANCE.** Making the coupling term
marginal-orthogonal and having it REPLACE the free pair term remains the right
construction, because it stops the pair term absorbing marginal effects. But this
measurement predicts the coupling head will **not** be the source of a large AP
gain on these labels. A large reported gain from a coupling head on this corpus
should first be suspected of marginal leakage into the coupling term.

**FAILURE LOCALIZATION SHARPENED.** Excluded: label insufficiency for marginals
(Oracle-RA = 0.391, there is much to find); optimization (trainability control
0.759); closure (196 components, largest 15.2 %); atom correspondence (I-1).
Partially confirmed: label insufficiency **for coupling specifically** (median
z 0.41). **Leading explanation: biological representation failure** -- B4
recovers 6 % of the marginal ceiling and the 0.368 AP gap lies in the residue
representation. The prior state record listed representation as one of several
competing explanations; this audit removes "there is nothing to find" from that
list. Section non-identifiability: NOT REACHED, no adaptation attempted.

**WHAT THIS SUPPORTS.** Exactly the registered B5 design as written: change
**only** the residue features, holding the ligand branch, head, rank, sampler,
budget and evaluation mask fixed. It does **not** authorise attention, a geometry
branch, a PLM larger than the registered 650M, an affinity head, or any
additional branch. No affinity, ranking, transfer, few-shot or `z`-admission
claim is made.

**BLOCKERS.** Discharged: `ATOM_CORRESPONDENCE_NOT_FULLY_VERIFIED`. Still open:
tie-aware AP and per-pair prediction materialisation; negative and control
manifest completion; publication/time closure; ESM2-650M weight acquisition.
B5 remains operationally blocked.

**ARTIFACTS.** `report/s7_l2b_r0r/I1_ATOM_CORRESPONDENCE_AUDIT.json`,
`I1_ATOM_QUARANTINE.json`, `I2_COUPLING_IDENTIFIABILITY_AUDIT.json`,
`I2_COUPLING_IDENTIFIABILITY_REPORT.md`, `i2_console.txt`. Code under
`research/s7_l2b_r0r/`. Frozen surfaces unmodified; nothing committed or pushed
beyond the previously committed preregistration.


## F-102: Phase 0 integrity repair complete; B5 PASSES ALL SIX GATES; gain is residue-side and largely generic (2026-08-10)

```text
Phase 0 blockers ............ ALL DISCHARGED
B5 registered Gates ......... 6 of 6 PASS
Failure localization ........ BIOLOGICAL REPRESENTATION, now resolved for residue localisation
Ligand-conditioning ......... WEAK — 92.5% of residue localisation survives a ligand swap
Confirmation cohort ......... SEALED, never opened
```

**PHASE 0 — EVERY BLOCKER DISCHARGED.** Atom correspondence verified at
**375,311** positions; 14,585 admitted, 4 enumerated and quarantined, quarantine
wired into the data contract. Tie-aware AP and sealed per-pair predictions
materialised over **52,062,975** held-out cells per arm as hashed float16.
Negative sampler audited at the contract's own granularity: exactly six per
positive, unique within each positive block, **zero** negatives that are actually
positives (cross-positive repetition 4.8 %, reported, not a violation).
Publication/time closure built and frozen from RCSB. ESM2-650M acquired by
range-resuming curl and **SHA-256 verified** (`c874668852...`, revision
`08e4846e`), then run offline. Phase 0 committed as `139effd` **before** B5 was
scored; the preregistration was already committed as `ce186f4`.

**TWO PHASE-0 CORRECTIONS RECORDED RATHER THAN SMOOTHED OVER.** First, the
earlier determinism comparison was **confounded**: the I-1 quarantine changed
train from 9,758 to 9,757 between the two runs, so the differing checkpoint hash
was expected and was NOT evidence of non-determinism. A same-data test in one
process then gave **bit-identical** state dicts. Second, the ligand-only baseline
is **massively tie-dependent** — optimistic 0.199 versus pessimistic 0.003 —
because a ligand-only model assigns every residue the same score; only the
tie-aware Monte-Carlo expectation is a defensible point estimate for it.

**B5 RESULT — held-out A, 2,409 complexes, 196 protein components, tie-aware
macro-AP over the complete residue x heavy-atom matrix.**

| arm | macro-AP |
|---|---:|
| B0 prevalence | 0.00319 |
| BM5 motif shuffle | 0.00451 |
| BP5 wrong protein | 0.00464 |
| BL ligand-only | 0.00572 |
| BX5 wrong ligand | 0.01968 |
| B4 non-PLM residue | 0.02325 |
| **B5 frozen ESM2-650M** | **0.06960** |

| gate | delta | LCB95 | |
|---|---:|---:|---|
| G1 B5-B0 | +0.06642 | +0.05998 | PASS |
| G2 B5-BL | +0.06388 | +0.05751 | PASS |
| G3 B5-BP | +0.06496 | +0.05849 | PASS |
| G4 B5-BM | +0.06509 | +0.05876 | PASS |
| G5 B5-BX | +0.04992 | +0.04424 | PASS |
| G6 B5-B4 | +0.04635 | +0.04039 | PASS |

Every lower bound clears the frozen 0.02 threshold by at least a factor of two.
**Only the residue features changed** — atom branch, head, rank 32, projected
dimension 128, sampler, optimiser, learning rate, weight decay, epochs, seeds,
split, evaluation mask and tie policy are identical to B4. This is the first Gate
PASS in the S7/L2B programme.

**THE FAILURE IS NOW DEFINITIVELY LOCALISED.** Optimisation had already been
excluded by the trainability control (0.759) and determinism is verified.
Changing only the residue representation tripled pair AP and cleared every Gate,
so the earlier `PROTEIN_SIGNAL_IDENTIFIED_BELOW_PREREGISTERED_EFFECT_SIZE` was
**biological representation failure**, exactly as I-2 predicted from the oracle
budget.

**THE GAIN IS ENTIRELY RESIDUE-SIDE — AND LARGELY GENERIC.** Marginal
decomposition from the sealed predictions, nothing retrained:

| arm | residue-marginal AP | atom-marginal AP |
|---|---:|---:|
| BL | 0.0313 | 0.7246 |
| B4 | 0.0879 | 0.6895 |
| **B5** | **0.2651** | 0.6796 |
| BX5 wrong ligand | 0.2453 | 0.5097 |
| BP5 wrong protein | 0.0434 | 0.6595 |

Residue marginal B5-B4 = **+0.1772 [LCB +0.1601]**. Atom marginal B5-B4 =
**-0.0099 [-0.0218, +0.0044]** — no gain, interval spans zero. ESM2 buys residue
localisation and nothing else; atom propensity is already near its ceiling for a
ligand-only model (BL 0.7246, the highest of any arm), matching I-2's finding
that the true atom marginal is worth only 0.0085 pair AP.

**CAVEAT THAT TRAVELS WITH THE PASS.** Swapping in a wrong ligand leaves
residue-marginal AP at **0.2453 against B5's 0.2651** — about **92.5 %** survives.
B5 predicts a **generic pocket**, not a ligand-conditioned one. The pair-level G5
gap appears in pair scores but may still be explained by additive residue and
atom marginals; it does not identify exact coupling or ligand-specific residue
localisation.
This is precisely the "strong pi alone means generic pocket localisation"
outcome, predicted in advance by I-2's median z of +0.41 against a
degree-preserving rewiring null.

**PUBLICATION/TIME CLOSURE — TIME-FORWARD CONFIRMATION IS INFEASIBLE FROM MONN.**
14,447 entries requested, 14,426 returned, 14,103 with a document key (DOI else
PubMed, never substituting a PDB ID); 323 lack a primary publication identifier
and are quarantined. Development 6,331 documents, additional-PDB 707, with only
**4** shared. But only **2** additional-PDB entries were released on or after the
frozen 2019-01-01 cutoff, and **zero** would qualify at 2024-01-01: MONN was
assembled in 2020 from PDBbind-v2018-era structures and contains no time-forward
holdout by construction. A document-closed confirmation cohort remains
constructible; a time-forward one does not and would need a different source
carrying residue-atom interaction labels.

**BOUNDARY CLASSIFICATION.** DATA/LABEL INSUFFICIENCY: excluded for marginals,
confirmed weak for coupling. BIOLOGICAL REPRESENTATION FAILURE: **resolved for
residue localisation**. OBJECTIVE/OPTIMIZATION FAILURE: excluded.
SUPPORT-SECTION NON-IDENTIFIABILITY: not reached.

**WHAT IS RETAINED AND WHAT IS REFUSED.** Residue localisation is retained as a
**biological statistic candidate**. Exact residue-atom coupling is **not**
claimed. No affinity, ranking, transfer, few-shot or `z`-admission claim is made.
No ChEMBL, BindingDB, DAVIS or recipient label was read.

**ARTIFACTS.** `report/s7_l2b_r0r/` — `P1_B5_REPORT.md`, `P1_B5_GATE.json`,
`P1_MARGINAL_DECOMPOSITION.json`, `P0_SEALED_PREDICTION_MANIFEST.json`,
`P0_DETERMINISM_AND_SAMPLER_CHECK.json`, `P0_ESM2_ACQUISITION_MANIFEST.json`,
`PUBLICATION_TIME_CLOSURE_AUDIT.json`, `I1_*`, `I2_*`, `GPU_PREFLIGHT_AUDIT.json`.
Sealed per-pair predictions and checkpoints under
`dataset/processed/s7_l2b_r0r/`. Code under `research/s7_l2b_r0r/`. Frozen
surfaces unmodified; nothing pushed.

**PHASE 2 PRECONDITION.** Phase 2 may be entered only after a new preregistration
is committed. Its expectation is registered in advance by this evidence: the
residue-first component should improve, and the coupling component is predicted
to be small. A large reported coupling gain on this corpus should first be
suspected of marginal leakage into the coupling term.

### F-103 — Active-tree convergence and pre-B5 archive (2026-08-10)

The repository's active view had accumulated superseded preregistrations,
withdrawn S5/S5R2 implementations, terminated XP/multipanel code, a failed
five-artifact R0 report, and duplicate console logs. Their scientific outcomes
were already consolidated in this ledger and the canonical evidence ledger,
but their presence made obsolete stages look executable.

Seventy-two files were removed from active `research/` and `report/` paths and
packed into `archive/legacy_pre_b5_20260810.zip` before removal. Archive SHA-256:
`a4c916cb09004383c0b55a6c8dd32df70748524874ada5c89a97c5d7ae9012e5` (361,895
bytes). The archive is historical only and must not be used as current evidence
or execution authority.

Retained active surfaces are the passed production primitives, reusable S1-S4
teacher/acquisition tools, the implemented S7/L2B Phase 0/1 path under
`research/s7_l2b_r0r/`, MONN provenance, and Phase 0/1 machine artifacts under
`report/s7_l2b_r0r/`. `task.md`, `experiment.md`, `project_state.json`, and
`report/CURRENT_RESEARCH_STATUS.md` were reduced to the current B5/Phase 2A
boundary. No model, theory, affinity data, DAVIS label, or production `z` was
changed.


## F-103: Phase 2A audit — the LABELS are ligand-conditioned; B5 is not; edge coupling is absent from both (2026-08-10)

```text
Preregistration ............. 4e01401d..., frozen before any metric, NOT committed
Phase 0 contract ............ PASS, 26 artifacts, 7 fail-closed checks
Data identifiability ........ DATA_IDENTIFIABLE
Teacher ligand-conditioned .. YES  (dJ +0.258 [LCB +0.234], rho +0.322)
B5 residue marginal ......... GENERIC  (89% wrong-ligand retention)
B5 coupling ................. real but BELOW the registered 0.01 margin
Teacher edge coupling ....... NOT identified (median z +0.413 vs threshold 2.0)
Label semantics ............. NOT ambiguous
Terminal verdict ............ LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING
```

**THE CORRECTION THIS ENTRY CARRIES.** F-102 reported that a wrong ligand
retained 92.5 % of B5's residue AP and concluded "mostly generic pocket". That
conclusion was about **B5**, and it was being carried forward as if it were also
a statement about the **corpus**. It is not. The F-102 control substituted an
*arbitrary foreign* ligand. Phase 2A replaced it with the scientifically correct
comparison — a real alternative ligand of the same exact construct — against a
noise floor measured from the data itself: two crystals of the same construct
with the *same* ligand.

| component-macro Jaccard of residue masks | value |
|---|---:|
| same construct, SAME ligand, different crystal (noise floor) | 0.6361 |
| same construct, different scaffold-distinct ligand | 0.4165 |
| **T1 dJ, paired over 292 closure components** | **+0.2580 [LCB +0.2344]** |

Registered minimum meaningful effect was 0.05; the observed effect is five times
that. It survives on held-out A alone (+0.1911 [+0.1165], 27 components) and
under the unpaired fallback (+0.2197 [+0.1984]). It is chemistry-associated, not
noise: within construct, mask dissimilarity rises with ligand chemical distance
at Spearman **rho = +0.3221 [LCB +0.2987]**, with a ligand-permutation control
giving median per-construct p = 0.03 and 84.8 % of constructs on the positive
side. 80.4 % of scaffold-distinct pairs change the residue mask meaningfully
(Jaccard <= 0.5 and symmetric difference >= 3 residues).

**Roughly 44 % of the alternative-ligand mask difference is ligand-attributable;
the other ~56 % sits at the replicate noise level.** Any differential objective
must be designed expecting that.

**WHERE B5's LIGAND INFORMATION ACTUALLY LIVES.** Decomposing the sealed logits
on the complete, uniformly weighted mask (orthogonality achieved 1.18e-9 against
a registered 1e-8 tolerance; weighted ALS agreed with double centering to
1.07e-14):

| arm | full | residue marg. | atom marg. | additive | coupling C |
|---|---:|---:|---:|---:|---:|
| **B5** | **0.06975** | 0.04045 | 0.00514 | 0.03983 | **0.01133** |
| B4 | 0.02323 | 0.01619 | 0.00550 | 0.01510 | 0.00637 |
| BX5 wrong ligand | 0.01969 | 0.03595 | 0.00326 | 0.02769 | 0.00346 |
| BP5 wrong protein | 0.00464 | 0.00461 | 0.00492 | 0.00563 | 0.00355 |
| BL ligand-only | 0.00573 | 0.00330 | 0.00573 | 0.00573 | 0.00305 |

A wrong ligand retains **89 %** of B5's residue marginal but only **31 %** of its
coupling term. So B5 *does* use ligand identity — but only through the pair term,
and that term is small.

**BOTH COUPLING CRITERIA FAIL THEIR REGISTERED BAR, AND ARE NOT ROUNDED UP.**

| contrast | delta | LCB95 | margin | |
|---|---:|---:|---:|---|
| B5 coupling - degree-preserving rewiring null | +0.00601 | +0.00461 | 0.01 | FAIL |
| B5 coupling - BX5 coupling | +0.00787 | +0.00620 | 0.01 | FAIL |

Both are clearly separated from zero and both are below the practical margin
that was fixed before the numbers were seen. The teacher's own edge coupling
fails too: median **z = +0.413** with 63.4 % of complexes above their own null,
against a registered threshold of z >= 2.0. This reproduces I-2's +0.41 under a
*stricter* rewiring specification (100xE burn-in, 30xE between samples, 20
independent rewires, **zero** degree-preservation violations) — an independent
confirmation rather than a repetition. Mixing was checked, not assumed: edge
overlap decays 1.000 -> 0.334 -> 0.298 -> 0.292 -> 0.292 at 0/1/5/10/30 swaps per
edge, and successive samples overlap at 0.292, the degree-constrained plateau.

**THE HEADROOM NUMBER THAT REFRAMES THE PROGRAMME.** The well-posed label-fitted
additive ceiling — the AP obtainable by recovering the *true* residue and atom
margins exactly — is **0.3889**. B5's full AP of 0.0698 is **17.9 %** of it; its
additive part 0.0398 is 10.2 %; its residue marginal 0.0404 is 19.8 % of the true
residue-margin ceiling 0.2043. **The bottleneck is the residue marginal, not the
coupling.** A pair-coupling head would optimise a term worth 0.011 while leaving
0.32 of additive AP unclaimed.

The logistic Rasch additive null registered in section 7 is reported but flagged
`rasch_converged: false`: the design is **completely separated** (the matrix is
0.07 % positive, so almost every residue row and atom column has no positive at
all and its coefficient diverges). Its AP is NOT a valid ceiling and is not used
as one. Recorded rather than quietly dropped.

**LABEL SEMANTICS — AUDITED, NOT AMBIGUOUS.** 212,556 typed edges: hydrophobic
68,153; H-bond 57,647; pi-stacking 34,917; salt bridge 25,447; water bridge
17,427; pi-cation 7,754; halogen 1,211. Water-mediated (indirect) edges are
**8.2 %**, below the 20 % threshold, and removing them **strengthens** T1
(dJ 0.258 -> 0.278) rather than reversing it. Metal-mediated edges: 0.0 %.

The registration allowed recording the dense-teacher question as UNRESOLVED *if
no comparator existed*. One did — 2,068 MONN entries already had local mmCIF from
earlier governed stages. Declaring it unresolved would have been a fabrication in
the convenient direction, so amendment 03 required building it. A local
distance teacher was constructed on **1,909** complexes (median mapped sequence
identity 1.000, exhaustive integer-offset scan, ligand copy chosen by a
label-blind rule): **88.1 %** of PLIP positives lie within 5.0 A of a ligand heavy
atom, and only 9.0 % beyond — consistent with the water-bridge fraction. Only 46 %
of geometric neighbours are PLIP positives, which is expected because PLIP
applies chemical *and* geometric criteria; a strict subset is correct behaviour
and is NOT evidence of missing positives. A second interaction-annotation tool
does not exist locally, so that specific comparison remains **UNRESOLVED**.
PU learning and a soft teacher stay unauthorized.

**PHASE 0 CONTRACT.** All seven checks passed fail-closed over 26 hashed
artifacts. The load-bearing one was C3: Phase 1's marginal decomposition indexed
the B5-family memmaps with the B4-family offset table. Both were rebuilt
independently and proved **identical key-for-key and offset-for-offset** (the ESM
filter dropped zero records), so the Phase 1 B5/BX5/BP5 marginals were correctly
aligned. Every sealed prediction hash matched its recorded manifest. Mask:
complete n_res x n_atoms, uniform weights, 52,062,975 cells, 36,237 positives,
density 6.96e-4 — which is what makes classical double centering admissible.

**CENSUS.** 14,585 records, 2,846 exact constructs, 2,916 UniProt ids, 1,994
closure components (largest 3.1 % of records), **1,093** constructs with
scaffold-distinct ligand pairs across **779** components, **292** components with
both pair types, **323,410** scaffold-distinct within-construct pairs, 2,408
replicate pairs. Label-blind power for dJ = 0.05 is >= 0.81 at the paired count
even at the most pessimistic registered sigma of 0.30.

**BOUNDARY CLASSIFICATION.** DATA/LABEL INSUFFICIENCY: **excluded** for
residue-level ligand conditionality; **confirmed** for edge-level coupling.
BIOLOGICAL REPRESENTATION FAILURE: **confirmed and localised** — the labels are
ligand-conditioned, B5's residue marginal is not. OBJECTIVE/OPTIMIZATION FAILURE:
not tested; Phase 2A trains nothing. SUPPORT-SECTION NON-IDENTIFIABILITY: not
reached.

**AUTHORIZED NEXT ACTION AND NOTHING ELSE.** One ligand-conditioned residue
residual head, preregistered as
`research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL.md`
(SHA-256 `ae6d1a01...`): `logit p_r(P,L) = b_r(P) + delta_r(P,L)` with `b_r`
frozen, `delta_r` one low-rank bilinear residual (K <= 8) over existing frozen
states, projected away from constant / pocket-prior / ligand-only directions,
supervised only by the same-protein ligand differential on symmetric-difference
residues, with Gates D1-D5, a replicate-oracle ceiling, and a fail-closed
module-participation audit. **Registered, not authorized. No Phase 2B code
exists.**

**CHRONOLOGY LIMITATION, STATED RATHER THAN CLAIMED AWAY.** Commit authorization
was not granted for this run. The preregistration and its three amendments are
anchored by SHA-256 and embedded in every output artifact, but carry **no git
commit timestamp** — strictly weaker than the guarantee behind `ce186f4` /
`139effd`. The Phase 2A verdict should be read with that attached until the
files are committed.

**ARTIFACTS.** `report/s7_l2b_r0r/` — `PHASE2A_SYNTHESIS.md`,
`PHASE2A_VERDICT.json`, `PHASE2A_INPUT_MANIFEST.json`,
`PHASE2A_DATA_IDENTIFIABILITY_CENSUS.json`, `PHASE2A_CONSTRUCT_GROUPS.json`,
`PHASE2A_TEACHER_CONDITIONALITY.json`, `PHASE2A_MARGINAL_COUPLING_AUDIT.json`,
`PHASE2A_COMPONENT_TABLES.json`, `PHASE2A_LABEL_SEMANTICS.json`,
`PHASE2A_PREREGISTRATION_HASH.json`, `pa3_console.txt`, `pa4_console.txt`.
Code `research/s7_l2b_r0r/pa0..pa5*.py`. Regression 75 passed. No affinity,
DAVIS, KIBA or recipient read. Frozen surfaces unmodified. Nothing committed or
pushed.

## F-104: Phase 2B contract repaired, then stopped fail-closed at its own synthetic control (2026-08-10)

```text
Superseded prereg ........... ae6d1a01..., NEVER EXECUTED, 11 design defects
New prereg R1 ............... 5e6688f6..., committed b9753db BEFORE any code
Contract audit .............. PHASE2B_CONTRACT_PASS, 14/14
Census ...................... matched the registration exactly
Synthetic trainability ...... FAILED, AP_bidir 0.3577 vs required 0.50
Real-label training ......... NOT EXECUTED
Gates R1-R6 ................. NOT SCORED
Terminal verdict ............ PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED
```

**THE PREREGISTRATION WAS REPAIRED BEFORE IT WAS USED, NOT AFTER.** The Phase 2B
preregistration written at the end of Phase 2A contained eleven defects that
would have made its result uninterpretable. It was never executed, so nothing is
withdrawn — a design was replaced before use. It is kept **byte-identical** and
marked `SUPERSEDED_BEFORE_EXECUTION_DESIGN_DEFECT`. Four defects mattered most:
`b_r(P)` was defined as the residue term of a **per-complex** additive
decomposition of B5 pair logits, an object fitted per `(P,L)` and therefore
ligand-dependent, so the same-protein cancellation the entire differential design
rests on would not have held — and it existed only on held-out A, never on the
training split; the projection span `{1, b(P), c(L)*1}` contained two collinear
columns, making the stated projector singular; the primary metric ranked only the
**symmetric-difference** residues, i.e. chose the candidate set using the answer;
and the module-participation audit demanded that detaching the frozen `h_r`
change the result, which is a mathematical no-op that could never be satisfied.

**A SEPARATE INTEGRITY FINDING.** `P1_B5_REPORT.md` no longer matches the hash
recorded for it in the Phase 1 triage (`19c9c205...` -> `dbfe8b92...`). The
change is wording only, in section 3(b), and no number moved. Recorded in
`PHASE1_ARTIFACT_SUPERSESSION.json` rather than reverted or ignored.

**EVERYTHING ELSE PASSED.** `PHASE2B_CONTRACT_PASS` on all fourteen preflight
items: exactly **10,568** trainable parameters with names `{U, V}` and no bias;
`g(L)` atom-permutation invariant at **0.0**; ligand-order swap sign-exact at
**0.0**; the protein-only prior cancels in the same-protein difference at
**2.05e-15** against a 1e-12 tolerance; projection orthogonality **6.19e-15**
against 1e-8; the degenerate-`b` Gram-Schmidt fallback returns rank 1 for a
constant prior and rank 2 otherwise; **zero** train/held-out closure-component
overlap; **zero** held-out ligand-graph overlap with training; **zero** held-out
B scaffold overlap; **zero** affinity-marked paths opened.

Census verified rather than assumed, and it matched the registration: train
**760** constructs / **554** components / **226,765** eligible pairs; held-out A
**174** constructs / **112** components / **46,818** eligible pairs; held-out B
**30,661**. Foreign-pair control coverage **1.000** over all 46,818 held-out
pairs from a 7,546-ligand training pool; within-construct derangement over
**11,123** records with **0** fixed points.

**TWO DEFECTS IN MY OWN PREFLIGHT, FOUND AND FIXED BEFORE THE RUN.** First, the
gradient-reachability probe used `d.sum()` as its objective. The constant
direction is the first column of `Q`, so the projection annihilates exactly that
functional; the gradient came back at ~7e-12 — a "pass" that proved nothing.
Replaced with a fixed generic random linear form. Second, the ESM-availability
check demanded states for every construct in the corpus and failed on 573. Scoped
to what Phase 2B actually requires: **0 missing**. Of the 573, 445 belong to the
**sealed additional-PDB confirmation cohort**, which correctly has no states, and
131 are development constructs whose records were all removed by the ligand-graph
disjointness filter; **0** of them appear in any Phase 2B split.

**THE STAGE STOPPED AT ITS OWN SYNTHETIC CONTROL.** The teacher is a rank-8
projected bilinear differential that lies **exactly** in the candidate's
hypothesis class by construction.

```text
required   AP_bidir >= 0.50
observed   AP_bidir  = 0.3577       chance 0.0376       FAIL
```

The threshold was **not** lowered, nothing was tuned against the synthetic
holdout, and no second seed was tried. Four gauge-invariant diagnostics localise
the shortfall, and none of them is a gate:

| diagnostic | value | reading |
|---|---:|---|
| teacher scored on its own labels | **0.99971** | metric and evaluation code sound |
| in-sample AP, final-epoch sampled pairs | **0.3654** | - |
| held-out AP | **0.3577** | **no generalisation gap**; this is underfitting |
| output-level Pearson r(learned field, teacher field) | **0.717 mean / 0.754 median** | the class IS being fitted |
| parameter movement U / V | **1.808 / 0.426** | the head trained; it did not sit still |

So the failure is **not** the hypothesis class, **not** the evaluation code and
**not** generalisation. It is the registered **optimization budget**: 6 epochs
over at most 8,864 sampled pairs drive the learned field to r ~ 0.75 of the
teacher but not into an exact top-8 ranking at AP >= 0.50.

**AND A SECOND POSSIBILITY THAT MUST BE STATED, NOT ARGUED AWAY.** The 0.50
threshold was set a priori with no scaling curve to calibrate against. Recovering
a field correlation of 0.75 from a random rank-8 teacher is substantial.
Distinguishing "budget too small" from "threshold too strict" is what the next
registration must settle **before** touching real labels — and it cannot be
settled by adjusting either number now, because both were frozen and the
synthetic holdout has been seen.

**NOTHING BIOLOGICAL IS CONCLUDED.** This run does **not** show that the frozen
sequence + 2-D ligand representation lacks a ligand-conditioned residue
correction. `R1`-`R6` were never scored. The Phase 2A finding stands unchanged:
the MONN labels are ligand-conditioned at the residue level (dJ +0.258
[LCB +0.234], chemistry association rho +0.322) and B5's residue marginal is not.

**SOLE AUTHORIZED NEXT ACTION.** Preregister one repair of the Phase 2B
**optimization contract** — not the biology, not the architecture, not the gates.
It must fix the budget and sampler caps from a **measured synthetic scaling
curve**, derive its acceptance threshold from that curve rather than by
assertion, use a **fresh synthetic teacher seed** because 20260905 has been
observed, and leave the architecture, projection, metric, controls and gates
R1-R6 byte-identical. Adding capacity, another PLM, attention, a GNN, a geometry
branch or a typed-interaction branch is **not** an admissible response to a
synthetic-precondition failure and is not proposed.

**DEVICE.** Phase 2B ran on **CPU**, chosen before any result was seen so the
registered bit-exact determinism check is achievable rather than a gamble on
cuBLAS reduction order; the head is 10,568 parameters and every heavy tensor is
frozen. CUDA 12.4 was available and deliberately unused.

**ARTIFACTS.** `report/s7_l2b_r0r/` - `PHASE2B_REPORT.md`, `PHASE2B_GATE.json`,
`PHASE2B_SYNTHETIC_AUDIT.json`, `PHASE2B_INPUT_MANIFEST.json`,
`PHASE2B_CONTROL_MANIFEST.json`, `PHASE2B_TRAINING_TRACE.json`,
`PHASE1_ARTIFACT_SUPERSESSION.json`, `p2b_console.txt`. Code
`research/s7_l2b_r0r/p2b_residue_residual.py`, `p2b_run.py`; audit
`PHASE2B_DESIGN_AUDIT.md`. Tests `tests/test_s7_l2b_phase2b.py`, 25 new.
Regression **100 passed** (75 pre-existing, verified). Commits `0bd1702`,
`b9753db`, `0a8b62e`. No affinity, DAVIS, KIBA or recipient read. Sealed
confirmation cohort not opened. Frozen surfaces unmodified. Nothing pushed.

## F-105: the Phase 2B synthetic control was itself invalid — SYNTHETIC_CONTROL_LOSS_MISALIGNED (2026-08-10)

```text
Prereg S0 ................... 81675578..., frozen before any S0 measurement
S0-A contract ............... PASS
S0-B candidate path ......... PASS
S0-C objective .............. REJECTED
S0-D / S0-E ................. NOT RUN (stage stops at the earliest cause)
Terminal verdict ............ SYNTHETIC_CONTROL_LOSS_MISALIGNED
```

**ONE MEASUREMENT SETTLES IT.** Initialise the student **exactly at the teacher**,
where `AP_bidir = 1.0000` by construction, and run the registered loss, optimizer
and sampler:

| updates | held-out AP | held-out BCE | train AP | train BCE |
|---:|---:|---:|---:|---:|
| 0 | **1.0000** | 0.63637 | 0.99975 | 0.64644 |
| 1 | 0.9732 | 0.63179 | 0.96065 | 0.64233 |
| 10 | 0.5797 | 0.59075 | 0.56338 | 0.60726 |
| 100 | **0.3899** | 0.42286 | 0.38189 | 0.45125 |
| 210 | 0.4963 | 0.38452 | 0.48662 | 0.40955 |

BCE falls monotonically while AP collapses from 1.0 to 0.39. The preregistered
rule — misaligned iff BCE(100) < BCE(0) and AP(100) < AP(0) - 0.05 — fires on
both conditions. **A control whose own answer is destroyed by its own training
procedure cannot adjudicate a student.** The failed Phase 2B precondition
(AP 0.3577 against >= 0.50) therefore carries no information about the candidate:
the quantity being optimised was not the quantity being scored.

**WHY THE RAY AUDIT COULD NOT HAVE SHOWN THIS.** AP is exactly scale-invariant
along the teacher ray, registered in advance as an expectation and confirmed to
machine precision: AP(a=1e-3) = AP(a=1e3) = 0.99999999999999978. Every BCE change
along the ray is metric-free by construction, so the decisive test had to be
**directional**. The ray audit nevertheless produced the number that explains the
mechanism: BCE at the teacher's own scale is **0.63637**, while the minimum along
the ray is **0.34075 at a\* = 20.20**.

**THREE MECHANISMS, SEPARATED AND RANKED.** An addendum diagnostic, run after the
registered rule had already fired and incapable of changing the verdict,
initialised at the ray optimum to remove the radial pressure:

| updates | AP from teacher (a=1) | AP from ray optimum (a\*=20.2) |
|---:|---:|---:|
| 0 | 1.0000 | 1.0000 |
| 10 | 0.5797 | 0.9755 |
| 100 | 0.3899 | 0.9121 |
| 210 | 0.4963 | **0.8913** |

*Primary, SCALE.* AdamW moves each coordinate by roughly the learning rate
irrespective of gradient magnitude: lr 1e-3 against a mean |U\*| entry of 0.0223
is **4.5% of the parameter's own scale per update**. Chasing a 20x scale change
rewrites the direction coordinate-by-coordinate long before the scale is reached.
Remove the radial pressure and the collapse largely disappears.
*Secondary, RESIDUAL MISALIGNMENT.* Even at the ray optimum, BCE still falls
(0.34075 -> 0.33669) while AP degrades (1.000 -> 0.891). The BCE optimum inside
the rank-8 class genuinely is not the teacher — real, but an order of magnitude
smaller.
*Tertiary, KNIFE-EDGE LABELS.* The rank-8/rank-9 gap has median **0.00222**
(quartiles 0.00081 / 0.00222 / 0.00495), so the synthetic top-8 boundary is an
essentially arbitrary tie-break on a continuous field.

One further fact worth stating plainly: the trajectory's endpoint, BCE 0.3845 at
AP 0.4963, is **worse in BCE** than simply rescaling the teacher (0.3408 at
AP 1.0). The registered pipeline was not even minimising its own surrogate well.

**THE PREVIOUS REPORT'S ATTRIBUTION IS DOWNGRADED.** `PHASE2B_REPORT.md` section
5 named "the registered OPTIMIZATION BUDGET" as at fault. That is **not uniquely
established** and the S0 evidence does not support it. The budget was never
tested, because an earlier-applicable cause fired first; and the budget
hypothesis alone cannot explain a pipeline that *discards the answer when handed
it*. The original result
`PHASE2B_NOT_RUN_SYNTHETIC_OR_NUMERICAL_PRECONDITION_FAILED` is preserved
unchanged as historical evidence in `PHASE2B_GATE.json`, `PHASE2B_REPORT.md` and
F-104. It is not rewritten as a PASS, and seed 20260905 was not rerun to rescue
it.

**CONTRACT AND IMPLEMENTATION ARE SOUND.** S0-A materialised and hashed the
complete 48-epoch nested stream (`5d695beb...`): **1,680** optimizer updates,
**166,300** pair presentations, **30,552** unique pairs — three quantities
counted separately, giving 13.5% unique coverage of the 226,765 eligible training
pairs at maximum budget. Antisymmetry error **0.0**; identical-ligand
differential **0.0**; projection orthogonality 5.7e-15; same-seed replay
reproduced the checkpoint hash and predictions bit-identically. Panels frozen
before any metric: train 14,202 by the hash-stratified rule, held-out 46,818
complete. S0-B copied the teacher's U\*, V\* into the production head and
reproduced the teacher to **4.06e-07** relative field error (tolerance 1e-4), AP
agreement **2.2e-16** (tolerance 1e-3) and W = U^T V to **3.46e-08** (tolerance
1e-4), with the tolerances justified in advance from float32 accumulation over
1,280 terms.

**CAUSE LEDGER.** (1) objective/metric misalignment **ESTABLISHED**;
(2) sampled-pair coverage **NOT TESTED**; (3) optimizer updates **NOT TESTED**;
(4) low-rank factorization **NOT TESTED**; (5) implementation/reproducibility
**EXCLUDED**; (6) synthetic generalization **EXCLUDED** — train and held-out AP
track at every checkpoint. Causes 2-4 are neither established nor excluded,
because the preregistration requires stopping at the earliest applicable cause.

**TAIL DIAGNOSTICS at the registered budget from random init**: AP 0.3797,
chance 0.0106, oracle-normalised AP gain 0.3731, Spearman 0.700/0.730, Kendall
tau-b 0.518/0.534, top-8 recall 0.379, top-16 0.522, top-32 0.640, AUPRG 0.0386,
Pearson 0.725/0.752 (one diagnostic among several, never used alone). Spearman
0.73 with top-8 recall 0.38 is the knife-edge label effect in metric form.

**SOLE AUTHORIZED NEXT ACTION.** A separate preregistration for a repaired
synthetic control, written and hashed but **not executed**:
`research/s7_l2b_r0r/PREREG_PHASE2B_S1_REPAIRED_SYNTHETIC_CONTROL.md`
(`4850c7d5ce23db35...`). It repairs the **control only** — the Phase 2B candidate
contract is carried over byte-identical — with R-1 emit the teacher at the
loss-preferred scale; R-2 require the top-8/bottom-8 boundary gap to exceed
0.25 x IQR, excluding and counting knife-edge pairs; and R-3 a mandatory
**alignment certificate**: starting at the emitted teacher, the registered budget
must not cost more than 0.05 AP on every calibration seed, or the control is
still invalid and no student number may be reported as evidence about the
candidate. R-3 is the methodological point S0 forces: a synthetic control must
first prove it is a valid control. The previous one never had to, and that is why
a meaningless number was nearly read as a biological result.

The registered Phase 2B **real-label** loss is deliberately not repaired: S0
raises a genuine concern that the same scale pathology would affect real
training, but combining a control repair with a loss repair is forbidden and must
be its own stage.

**NOTHING BIOLOGICAL IS ESTABLISHED OR EXCLUDED.** No real Phase 2B label decided
anything, no affinity value was read, the sealed verification seed 20260999 was
never touched, and no U/V latent channel is given any biological interpretation.
AP_bidir >= 0.50 was not lowered. Seeds 20260905 (development) and
20260911/12/13 are burned; 20260999 remains sealed.

**ARTIFACTS.** `report/s7_l2b_r0r/` — `PHASE2B_S0_FAILURE_LOCALIZATION_REPORT.md`,
`PHASE2B_S0_VERDICT.json`, `SYNTHETIC_INPUT_AND_STREAM_MANIFEST.json`,
`SYNTHETIC_REPLAY_AND_DETERMINISM_AUDIT.json`,
`SYNTHETIC_CANDIDATE_PATH_WITNESS.json`,
`SYNTHETIC_OBJECTIVE_COMPATIBILITY_AUDIT.json`,
`SYNTHETIC_SCALED_TEACHER_ADDENDUM.json`, plus explicit NOT_RUN records for
`SYNTHETIC_CONTINUOUS_FIELD_WITNESS.json`, `SYNTHETIC_FULL_W_CONVEX_WITNESS.json`,
`SYNTHETIC_EXPOSURE_SCALING_CURVE.json` and
`SYNTHETIC_SEALED_VERIFICATION_GATE.json`. Code `research/s7_l2b_r0r/s0_synth.py`,
`s0_run.py`, `s0_c2_scaled_teacher.py`. Console `s0_console.txt`. Device CPU,
CUDA 12.4 available and deliberately unused. Regression **100 passed**. Nothing
committed: no commit authorization exists for this stage, so the S0 and S1
registrations are anchored by hash only, which is weaker than the guarantee
behind `b9753db` and is recorded rather than claimed away.

## F-106: S2R repaired synthetic trainability; S3R did not identify the real residue direction (2026-08-10)

S0R replay established that the original S0 verdict had been computed on only
2 of 112 components and therefore could not stand as a contract-level result.
The subsequent sequence preserved every historical artifact and isolated the
failure without lowering a Gate:

```text
S0R complete-panel replay ........ contract and panel repaired
S1R factorized pairwise learner .. gauge/scale drift persisted
S2R direct bounded W ............. BINARY_ORDINAL_IDENTIFIABILITY_REPAIRED
S3R real structural transfer ..... REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED
```

S2R removed the non-identifiable U/V factor gauge and unbounded score scale by
training one direct `1280 x 41` matrix on the unit Frobenius sphere. Three fresh
calibration seeds and one sealed seed passed; sealed held-out component-macro
`AP_bidir = 0.662021`. This established optimizer/estimator trainability only.

S3R reused the estimator on real MONN residue-differential labels with frozen
ESM2 residue states, frozen mean-pooled 41-D ligand atom features, 210 fixed
updates and no hyperparameter selection. Primary census: 46,818 pairs and 112
closure components.

| arm | AP_bidir |
|---|---:|
| candidate | 0.035880 |
| chance / zero W | 0.025472 |
| frozen B5 differential | 0.031582 |
| foreign ligand pair | 0.035735 |
| context corruption | 0.032336 |
| trained permuted-label learner | 0.037125 |

R1 candidate-minus-chance was `+0.010408 [LCB +0.006920]`, below the frozen
`+0.05` margin. R2-R5 also failed: `+0.004298`, `+0.000145`, `+0.003544` and
`-0.001245`, respectively, with their registered margins unchanged. The
earliest terminal verdict is `REAL_BINARY_RESIDUE_DIRECTION_NOT_IDENTIFIED`.

This is not an optimization or participation failure. Gradient, movement,
unit norm, score variance, zero-W chance, common masks, identical stream and
bit-exact repeat predictions all passed. The correct scope is narrower: the
current ESM2 residue representation plus a global mean of 41-D ligand atom
features did not identify the real ligand-conditioned residue direction under
closure shift. Phase 2A had already shown that the labels themselves contain
ligand conditionality, so the result localizes the bottleneck to the current
measurement basis/estimand rather than a flat corpus.

Heldout-B was not opened. R6 was superseded before execution because ordinal
ligand differences do not identify an absolute output scale, ligand-feature
origin or directions outside the difference span. Affinity reads remained zero;
`model/`, `scripts/`, `theory/`, CSMO, Band and biological `z` were untouched.

One duplicate `prepare` invocation was correctly rejected by no-clobber after
the valid invocation had already written its manifests. The raw fail-closed
artifact is retained and its chronology is adjudicated in
`PHASE2B_S3R_ORCHESTRATION_ADJUDICATION.json`; it did not control the scientific
verdict. A string-typed `"True"` in the unit-norm subfield is also preserved;
the raw norm is `1.0000000116` and the aggregate participation result is Boolean
PASS.

No active experiment is authorized. The only evidence-aligned future proposal
is a separately registered single-axis ligand-information audit replacing the
global ligand mean with a frozen graph-aware 2D statistic while holding the
protein states, direct-W estimator, loss, split, stream and R1-R5 fixed. Full
repository regression after consolidation: **134 passed** in the `drug`
environment.

## F-107: the ligand mean was a real bottleneck, but graph information does not make the residue direction ligand-specific (2026-08-10)

S3R had failed on the basis `frozen ESM2 residue states x mean-pooled 41-D
ligand atom features` while passing every numerical, participation, firewall
and replay check. The leading hypothesis was that mean pooling destroys the
ligand topology needed for ligand-conditioned residue selection. F-107 tested
that hypothesis on one axis and closed it.

```text
S4R-A label-blind representation audit .. GRAPH_LIGAND_REPRESENTATION_AVAILABLE_AND_INFORMATIVE
S4R single-axis graph-aware transfer .... REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED
```

The audit was registered and committed before any audit code existed, and an
amendment replaced an unverified collision-free claim with unfolded 32-bit
Morgan identifiers plus a measured folding-collision count. It read zero
residue labels. It established that the mean-pooled 41-D basis is collapsed:
pair-difference effective rank `6.183` over 39,435 label-blind heldout-A
ligand-graph pairs, 687 distinct ligand graphs sharing a bit-identical vector,
and `85.2%` of the difference-norm variance explained by heavy-atom-count
difference alone. `m`-xylene and `p`-xylene, identical in atom composition and
different in connectivity, map to the same vector; that case is now a test.

All six audited Morgan candidates cleared every A-gate. The registered
capacity-parsimony rule — smallest `(d, radius)` among the admissible — selected
radius 1, `d = 128`, per-heavy-atom environment counts over a train-only
vocabulary, raising the difference effective rank to `20.93` with `35.5%` of
its energy beyond the baseline's linear span and `163,840` matrix parameters
against the baseline's `52,480`.

S4R then changed only that statistic. Three anchors prove nothing else moved:
the training stream's semantic and file SHA-256 both equal S3R's, the
common-mask SHA-256 equals S3R's, and the `baseline41` arm reproduced the S3R
candidate exactly, `|delta| = 0.0`, with C2 reproducing the S3R R1 interval to
every digit.

| arm | AP_bidir |
|---|---:|
| candidate, graph-aware `d=128` | 0.046856 |
| baseline41, mean-pooled `d=41` | 0.035880 |
| trained permuted-label learner | 0.036293 |
| frozen B5 differential | 0.031582 |
| foreign ligand pair | 0.046212 |
| residue-context corruption | 0.027357 |
| ligand-only / zero-`W` chance | 0.025472 |
| within-construct chemistry shuffle | 0.051322 |

R1 observed `+0.021384 [LCB +0.016064]` against a `+0.05` margin. R2 `+0.015273`,
R3 `+0.000644 [LCB -0.009226]`, R3b `+0.021384`, R4 `+0.019498`, R5 `+0.010563`,
each below its registered margin. The earliest failed boundary is R1, so the
terminal verdict is `REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED`.

Two findings matter more than the verdict. First, the hypothesis was partly
right: the above-chance gain doubled from `+0.010408` to `+0.021384`, the
direct contrast `C1 = +0.010976 [LCB +0.004939]` is above zero, and the
candidate now beats its capacity-matched permuted-label learner by
`+0.010563 [LCB +0.003880]` where S3R had lost by `-0.001245`. The ligand
representation was a genuine, measurable bottleneck and part of the signal it
was hiding is real.

Second, the recovered signal is not ligand-conditioned. Substituting a frozen
foreign ligand pair costs `+0.000644` and the within-construct chemistry
shuffle scores *above* the candidate; the leading singular values of `W` are
nearly flat (`0.209, 0.189, 0.187, 0.177, 0.170`). The learned residue
direction barely depends on which ligand difference is supplied, so what
improved is a construct-level residue-change prior, not ligand-specific residue
selection. The ligand-only arm is exactly chance by construction, since a
residue-constant field lies in `span{1} subset span{Q_P}` and the gauge
annihilates it; R3b is therefore a structural proof rather than an independent
contrast.

Module participation and deterministic replay passed in full: minimum
`|grad W| = 3.2617`, relative movement `1.4139`, unit norm `1.0000000207`,
score variance `1.16e-4`, zero-`W` equal to analytic chance, one shared stream,
and bit-identical repeat predictions.

Governance: heldout-B was neither created nor read, R6 was not opened, affinity
value reads were zero, and `103,116` real structural residue-edge labels were
read across the two created views. Heldout-A had already been consumed by S3R,
so every S4R number is development evidence and none of it is confirmation. No
threshold, seed, margin, budget, capacity or representation was changed after a
metric was read.

The registered stopping rule closes the pose-free ligand representation repair
route, including re-running S4R at `d = 256` or `d = 512`. It authorizes no
attention stack, larger PLM, second protein encoder, parallel branch, pose or
geometry branch, typed channel, affinity supervision, knowledge graph, PU loss
or few-shot adaptation. The open question is no longer representational
richness but whether any pose-free sequence-plus-2D estimand can bind a ligand
substructure to a residue context without geometric correspondence; that needs
a separately governed information stage. Full repository regression:
**159 passed** in the `drug` environment.

## F-108: the estimator does steer on the ligand — it steers somewhere biologically wrong (2026-08-10)

S4R had shown that a foreign ligand pair costs almost nothing. F-108 registered
one mechanism for that and then falsified its own hypothesis.

```text
S5D D1 ligand-steering collapse ....... NOT CONFIRMED, registered mechanism falsified
S5D D2 conditional estimand E1-E3 ..... ALL FAIL
terminal verdict ...................... LIGAND_DIRECTION_COLLAPSE_NOT_CONFIRMED
```

The stage trained nothing, introduced zero parameters, added no representation
and reused the frozen S4R checkpoints byte-for-byte. It did not reopen the
ligand representation route.

The registered claim was that the estimator maps every ligand difference onto
approximately one residue direction per protein, so a foreign ligand yields
nearly the same field. On 131 heldout-A constructs with at least three pairs,
the top principal energy fraction of the mean-centred unit fields was
`rho_graph = 0.4793` against a data-side upper bound `rho_dg = 0.4550`, an
excess of `0.0138` where the registered rule required a median of at least
`0.80` and an excess of at least `0.10`. The median cosine between a pair's
true field and its foreign field was `0.4487` over 46,817 pairs. The estimator
plainly does steer on the ligand, and `rho_base = 0.5758 > rho_graph`
independently confirms that the S4R representation change increased field
diversity exactly as the S4R-A audit predicted.

D2 then aggregated `ap_symdiff_conditional`, an estimator already implemented
and registered in `p2b_residue_residual.pair_metrics` and already aggregated by
the parent Phase 2B runner, which S3R and S4R computed per pair and never
aggregated. Restricting each comparison to the residues that changed makes
pocket membership constant across both classes, so it cancels exactly and
non-parametrically, with no gauge and no tuning. On 40,157 eligible pairs
across 107 closure components, median 7 changed residues and median gain
fraction `0.50`:

| arm | AP_cond |
|---|---:|
| candidate, graph-aware | 0.655030 |
| foreign ligand pair | 0.655470 |
| conditional chance | 0.643744 |
| baseline41, mean-pooled | 0.638830 |
| trained permuted-label learner | 0.628586 |

`E1 = +0.011285 [LCB -0.007749]` against a `+0.05` margin, `E2 = -0.000440
[LCB -0.021814]`, `E3 = +0.026444 [LCB -0.002977]`. Every Gate fails, E1 fails
on its lower bound as well as its margin, and `baseline41` sits *below*
conditional chance.

The joint reading is sharper than S4R alone. Ligand information is not lost
upstream: it reaches the residue field and rotates it by a large angle. It is
not diluted by the metric either: two estimands, one bidirectional and one
pocket-cancelling, agree that the foreign arm ties the candidate to within
`0.0004`. What the ligand determines about the residue direction is simply
unrelated to which residues gained or lost contact. That the whole above-chance
effect vanishes once pocket membership is cancelled is consistent with the S4R
`AP_bidir` gain having been sign-agnostic pocket structure surviving the
two-dimensional `span{1, b^P}` gauge, but S5D registered no test of that and
the observation is not claimed as a result.

A numerical limitation is disclosed rather than hidden: `rho` normalizes before
centring, so a totally collapsed construct floors near `0.87` instead of `1.0`.
The bias is downward in the degenerate limit, it can understate collapse but
cannot manufacture it, and the observed median `0.4793` is nowhere near that
regime. Both properties are asserted as tests.

Governance: heldout-B not created and not read, R6 closed, affinity value reads
zero, one label view opened. **Heldout-A has now been consumed three times, by
S3R, S4R and S5D**, so every number is development evidence and the panel is
weaker as evidence with each look; the registered stopping rule forbids a
fourth estimand variant on it. No threshold, seed, margin, arm or eligibility
rule was changed after a statistic was read.

Both pose-free repair routes are now closed: representation by S4R, estimand by
S5D. The remaining hypothesis is that the missing ingredient is
**correspondence** — which ligand substructure sits against which residue — and
that a pose-free sequence-plus-2D estimand has no channel to supply it. Testing
that is a separately governed information stage about geometry with its own
preregistration, and nothing here authorizes it. Full repository regression:
**174 passed** in the `drug` environment.

## F-109: exact atom-residue correspondence is nearly a function of contact degree (2026-08-10)

S5D had localized the remaining hypothesis to **correspondence** — which ligand
substructure sits against which residue. F-109 tested it audit-only, on a
corpus no MetaSieve stage had ever touched, and closed the route before any
model was trained.

```text
C0 untouched corpus and closure ...... ALL GATES PASS
C1 exact-coupling information ........ FAIL AT C1a
terminal verdict ..................... EXACT_EDGE_COUPLING_NOT_SUPPORTED_BY_TEACHER
C2 geometry-gated router ............. NOT PREREGISTERED, NOT TRAINED
```

An exposure registry unioned every PDB id consumed by P1B, the preflight QC
corpora, MONN/B5/S7/S3R/S4R/S5D and the ssl_b2 independent set: **24,874** ids.
Of 14,169 local raw mmCIF entries, 2,836 were untouched, 2,509 carried a
BioLiP2-relevant ligand, 2,039 systems were admissible and 1,862 survived the
CCD-scaffold rule, over 669 entries, 982 receptor sequences and 340 ligands.
BioLiP2 was annotation-only, PLINDER was not used under its standing licence
audit, affinity reads were zero, and heldout-A was never referenced.

The stage recorded and respected the P1B correction: `contact_prob(i,s)` is a
Bernoulli-like "any residue in slot `s` contacts atom `i`", never additive
contact mass, and multiple residues in one slot may contact the same atom. It
also recorded a structural consequence — P1B is constant across residues within
a slot, so it can gate channels but cannot by itself discriminate inside a slot.

The registered mapping rule `M4` then failed its **own** fail-closed check at
`23/40`. The premise was factually wrong: P1B's sequence comes from BioLiP
column 20, systematically shorter than the mmCIF entity sequence and sometimes
not a prefix. Amendment 01 replaced the rule with the true P1B path — BioLiP
receptor sequence plus parasail alignment, keyed on BioLiP rows — before any
statistic was read, and the check then passed `60/60` on slot assignment. The
C1 execution already running under the rejected mapping was stopped and its
output discarded unread. No Gate, threshold, margin or seed changed.

C0 passed everything: 496 inference components, largest fraction `0.0811`, and
a minimum detectable effect of `0.00453` against a `0.05` requirement, with the
null dispersion estimated from the degree-preserving arm only. The union
closure produced 89 components but exceeded the `0.25` giant-component cap, so
the registered DataSAIL-style two-dimensional fallback was used — the union
giant component was tested, not assumed. The 3-mer prefilter was measured
against exact brute-force alignment: 3,037 true identity edges, 0 missed.

| arm | component-macro within-slot AP |
|---|---:|
| empirical | 0.985611 |
| fixed-degree rewire null | 0.953959 |
| atom shuffle | 0.985611 |
| geometry shuffle | 0.993948 |

`C1a = +0.031652 [LCB +0.029690]` against a `+0.05` margin: FAIL. `C1b` passed
with 162,276 positive units and 100,563 valid 2x2 checkerboards read entirely
from raw coordinates. `C1c` passed but on only 17 cross-entry replicate pairs.

The decisive fact is the ceiling, not the Gate. Empirical within-slot AP is
`0.9856`, so **only `0.0144` of headroom exists above a predictor that ranks a
slot's candidate residues by their contact degree alone**. The `+0.05` margin is
unreachable in principle on this statistic even by an oracle, and the panel was
powered to `0.00453`, so this is an effect-size result and not a detection
failure. The mechanism is chemistry: at the frozen `6.0 A` P1B threshold a slot
holds about three sequence-adjacent, hence spatially adjacent, residues, and if
one is in contact its neighbours usually are too.

Two controls are recorded as degenerate rather than as evidence. Atom shuffle
is an **exact** no-op, `0.985611` to every digit, because the statistic ranks by
a residue-side column sum that row permutation cannot change. Geometry shuffle
scores *higher* than the empirical arm, because breaking sequence adjacency
makes degree-ranking easier — so the real slot partition is the hard case and
"beating the shuffle" would have been the wrong direction of test.

Three routes are now closed by preregistered Gates: representation (S4R),
estimand (S5D) and geometry-gated correspondence (C1). Nothing authorizes
widening the corpus, relaxing the `6.0 A` contract, changing the closure, or
any of the excluded modules. The C0 corpus remains a clean, never-scored asset
with its exposure registry and closure frozen and hashed. Full repository
regression: **193 passed** in the `drug` environment.

## F-110: the crossed ChEMBL37 dependence precondition passes, on a pass that must be read with its bias (2026-08-10)

C1 closed the exact 6 A correspondence route. F-110 turned to the crossed
affinity estimand `DD = y(P1,La) - y(P1,Lb) - y(P2,La) + y(P2,Lb)` and ran the
X1A dependence precondition only. It trained nothing.

```text
X0/X0-B recovery and hashing ..... PASS, 3/3 data files byte-verified
G1 Ki  UCB95(rho) < 0.0915 ....... PASS  3.88e-07
G2 Kd  UCB95(rho) < 0.0164 ....... PASS  1.33e-03
G3 no cluster dominates .......... PASS  Ki 0.0387, Kd 0.2066
G4 effective n >= 245 ............ PASS  Ki 827.0, Kd 604.3
terminal verdict ................. X1_ICC_PRECONDITION_PASSED
X1B .............................. AUTHORIZED, NOT RUN
X2 ............................... NOT AUTHORIZED, NOT TRAINED
```

X0 and X0-B were recovered from `24a9ae0^`; `cells.jsonl`,
`dependency_components.jsonl` and `panels.jsonl` all match their manifest
SHA-256 exactly, while `report.json` does not and is disclosed rather than
repaired. The X0-B design and its statistical unit were not rebuilt. ChEMBL37
was opened only after the preregistration was committed, reading four fields
for 63,859 activity ids already enumerated by the label-blind census; all
carried `standard_relation '='` and a pChEMBL value. BindingDB, DAVIS, KIBA,
PDBbind and recipient reads were zero and no OOF residual was computed.

Two defects were found in this stage's own instruments. First, the registered
ICC estimator was **degenerate**: fitting additive target and ligand effects
within panel forces every panel residual mean to zero, so `var(cluster)` is
identically zero for any dataset — proven on synthetic panels with injected
10/20/30 log-unit offsets. The first execution returned `rho = 0.0000` and
would have PASSED; that result is void, and amendment 01 replaced the
within-panel fit with a global per-endpoint fit while moving no threshold.
Second, G3 and G4 were computed on measurement counts rather than the
registered X0-B cell-disjoint DD unit; both now use the frozen X0-B per-cluster
sizes and reproduce X0-B's capped totals exactly (`Ki sum(min(size,32)) = 827`).

The pass is real against the registered Gates and must be read with three
caveats, all recorded in the artifacts. The additive fit consumes 42% (Ki) and
32% (Kd) of cells as parameters with 12.2% / 14.5% singleton ligands, so `rho`
is credible as a lower bound and the bias direction **favours passing**.
`var(panel)` truncated to zero for both endpoints and `var(cell)` for Ki, so
the nested decomposition is only partly identified. And replicate noise is
99.99998% (Ki) and 99.93% (Kd) of the adjusted variance:

```text
Ki  replicate SD 0.618 log units -> detectable interaction RMS 0.309
Kd  replicate SD 1.860 log units -> detectable interaction RMS 0.930
```

at X0's frozen 0.5 interaction-to-noise ratio. Kd clears the dependence
precondition while being close to unusable for anything smaller than a ten-fold
selectivity swing as an RMS. That is a detectability question, which X1A did
not test and X1B is built to adjudicate through
`I_real^2 = max(0, E[DD^2] - E[v_noise])`.

Nothing beyond X1B is authorized. No trainable component was added, no `q_theta`
was preregistered, the 3D route stays closed, support adaptation was not
implemented, and `model/`, production `scripts/`, `theory/`, CSMO, Band, the
mesh, production `z` and `A(F,z)=K(B(z)F(z))` are unmodified.

## F-111: X1A authorization is withdrawn; direct-DD dependence repair is registered (2026-08-10)

**CORRECTION WITHOUT HISTORICAL REWRITE.** Independent review found that F-110's
amended ICC cannot establish the dependence precondition for X1B. All 310
targets in the label-blind census are confined to one dependency cluster, so
the global target dummies absorb cluster membership. Singleton ligand dummies
also fit their cell means exactly. The remaining nonzero cluster variance is
created after changing weights through cell/panel/cluster averaging and is not
the intended random cluster effect.

The estimand is also mismatched: F-110 measured dependence of signed fitted
residuals, whereas X1B would test `q=DD^2-v_noise`. Positive and negative
interactions may cancel in the former while remaining dependent in the latter.
The final machine artifact used 2,000 bootstrap draws rather than the registered
10,000. Therefore:

```text
historical amended verdict ........ X1_ICC_PRECONDITION_PASSED (retained)
current authorization verdict ..... X1A_ICC_PRECONDITION_NOT_ESTABLISHED
X1B execution ...................... NOT_AUTHORIZED
X2 / GPU training .................. NOT_AUTHORIZED
```

**COMPLETED LABEL-BLIND REPAIR.** The original X0-B deterministic greedy packing
was restored and materialized as individual rectangles. It exactly reproduces
Ki 11,168 rectangles / 36 clusters and Kd 1,041 / 12; frozen caps select 827 and
605. `rectangles.jsonl` SHA-256 is
`22f3e738f4dbc7b53ca9ef23e995e2a398cbca280a9cdde12c546be21500d0a5`.
No affinity value was selected during materialization.

**REGISTERED NEXT ACTION.** `E-AFF-X1A-R_DIRECT_DD_DEPENDENCE` computes DD
directly from four cell means, fits no target/ligand nuisance model, and
estimates dependence of the same `q` statistic used by the interaction test.
Only its PASS may authorize preregistration of X1B. ChEMBL37 X1A scoped label
access is recorded as 63,859 pChEMBL rows; affinity training reads, BindingDB,
DAVIS, KIBA, PDBbind and recipient reads remain zero.

**C1 NUMERICAL CORRECTION.** The previously quoted `0.014389` is
`1-AP_empirical`, the empirical residual to perfection. Maximum possible gain
over the fixed-degree null is `1-0.953959=0.046041`; this is the correct reason
the registered `+0.05` Gate is unreachable. F-109 remains historical and is not
rewritten.

## F-112: X1A-R direct-DD dependence fails; X1B and training stay closed (2026-08-10)

The repair was frozen in commit `4ce54c1` before any direct-DD value was
computed. Exact-assay alignment retained 827 Ki rectangles in 36 dependency
components and 590 Kd rectangles in 12. X1A-R evaluated the intended
cell-interaction-scale statistic `Z=(DD/2)^2-v_D,U` directly, fitting no target
or ligand nuisance model.

```text
Ki  rho_U95 0.120406 > 0.0915; n_eff 200.43 < 245
Kd  rho_U95 0.101078 > 0.0164; n_eff  61.05 < 245
terminal verdict: X1A_R_DEPENDENCE_PRECONDITION_FAILED
```

Cluster dominance itself passed, so the stop is caused by dependence and
effective information, not one oversized component. The run opened 5,986
preselected ChEMBL37 pChEMBL rows, trained no model and used no GPU. BindingDB,
DAVIS, KIBA, PDBbind and recipient reads were zero. The conditionally
preregistered X1B test was not run; X2 was not preregistered and no GPU training
was authorized. A new route requires a separately governed crossed source with
more independent components, not relaxed thresholds or a larger network. Full
repository regression after consolidation: **203 passed** in `drug`.

## F-113: cycle quotient recovers algebraic information but not independent components (2026-08-10)

The proposed interaction quotient was tested label-blind on the frozen X0
cells. Its algebra is valid: the document/context panel graphs contain raw
cycle dimensions 29,677 (Ki) and 3,279 (Kd), far more than the cell-disjoint
rectangle packing. But no exact assay spans multiple targets, so exact-assay
cycle dimension is zero. At panel level, the largest frozen dependency
component contains 48.9% of Ki and 46.1% of Kd cycle dimension; the independent
component counts remain 36 and 12.

```text
CYCLE_QUOTIENT_ALGEBRAICALLY_AVAILABLE_BUT_DEPENDENCY_NOT_REPAIRED
```

Thus cycle projection is retained as the preferred future affinity estimand,
but it does not reverse X1A-R or authorize training on ChEMBL. A BindingDB
curated-article metadata census is preregistered before acquisition. The future
model is deliberately minimal: frozen ESM2/P1B plus one four-coordinate
mechanism response, followed only after confirmation by a closed-form
row-space few-shot section. No parallel attention/GNN/adapter stack is added.
Repository regression after this audit: **207 passed** in `drug` using an
explicit workspace-local pytest base directory because the system temp root
was ACL-inaccessible.

## F-114: BindingDB quotient corpus opens governed development training (2026-08-10)

BindingDB Articles 202608 was joined to the official reaction-set/assay map.
The trusted extractor traversed the monolithic TSV but exposed no numeric value
until the separately frozen endpoint audit. Exact uncensored extraction yielded
24,157 Ki/Kd rows. Ki retained a noise-corrected quotient interaction RMS of
`0.5668 [0.5293, 0.6043]`; Kd was too small for its registered development
sample/rank requirements.

After deterministic removal of inconsistent stereo identities, scaffoldless
ligands and six explicit-H ligands outside the frozen heavy-atom contract, the
Ki training corpus contained 12,457 cells and 320 panels. Strict document,
protein-40%-identity and Murcko-scaffold union closure produced 31 components,
with train/development quotient ranks 6,608/220. The largest component share is
0.8586: optimization is authorized, population inference is not.

## F-115: first open-data affinity training runs; shared radial direction fails (2026-08-10)

The frozen ESM2/P1B/T-BASIS pipeline generated correct, foreign-ligand and
deranged-protein 288D features for all 12,457 cells on CUDA in 117.14 seconds.
One panel-balanced ridge linear response was fitted after train-component CV.

```text
correct RMSE        0.580314
zero RMSE           0.580520
explained fraction  0.000709
correct-zero        +0.000239 [-0.000981, +0.001496] loss reduction
correct-foreign     +0.000870 [-0.001045, +0.002847]
correct-deranged    -0.000817 [-0.003419, +0.001498]
```

Terminal verdict: `CQ_TBASIS_LINEAR_AFFINITY_WITNESS_NOT_OBSERVED`.
This is the first real open-affinity training run in the current programme. It
rejects a single population-shared linear direction on the fixed 288D radial
basis, not open-data training itself and not target-dependent coefficients.
The next scientifically distinct test is a `d<=5` source-learned target
coefficient subspace using dense profiling panels, followed by target-held-out
`k=1/2/3/5` support sections. No production statistic was admitted.
Full repository regression after consolidation: **224 passed** in `drug`.
