# FORT Evidence Ledger

## Established Findings

- A1-v2 stopped at 143/459 exact closures and 32/153 complete units, below the
  frozen 128-unit requirement. Structural transfer is closed.
- Strict dual-cold roles isolate train and development targets, homology
  clusters, scaffolds, connectivities, documents, and assays.
- The k=5 preprocessing audit retains 133 pKi development episodes in 126
  components; pKd retains 51 and has insufficient primary query depth.
- The current one-seed MAML-style pKi pilot did not exceed calibration:
  RMSE/MAE/Spearman were 1.309/1.042/0.072 versus 1.290/1.029/0.084.
- The corrected finite-rank Bayesian Wave 1 run improved calibration slightly
  but failed ligand-only and gradient controls: Bayesian RMSE/MAE/Spearman
  1.693/1.461/-0.043, ligand-only 1.324/1.070/0.098, and gradient
  1.247/1.012/0.073. Flexible-kernel work is blocked.

## Claims Not Established

- Protein-conditioned few-shot reordering gain.
- Query-span, Mamba, graph ligand, or Bayesian posterior benefit.
- Calibrated uncertainty or a valid confirmation result.

## Prohibited Interpretations

- Rows, supports, queries, seeds, folds, and pairs are not independent
  biological units.
- A random target split from external literature is not evidence for strict
  target/homology-cold, scaffold-cold performance.
- More MAML epochs cannot rescue the failed E0 comparison.
- Confirmation and sealed outcomes remain unavailable.

## Current Program Transition

SIMA-DTA is an independent Bayesian target-adaptation program. AdaMBind is an
external reference only: no code, name, task split, or model component is
retained. The primary adapter will be a decoupled Bayesian reordering posterior
with Helmert contrasts, evidence gating, and uncertainty. The current gradient
adapter is retained only as a MAML baseline. The three active modules are the
edge-aware long-context encoder, Bayesian reordering adapter, and
counterfactual query-span episode design.

## Stability Finding

The historical E0 run was unstable by design: sequential target updates had
high variance, train support selection differed from the development roster,
the residual readout started at zero, and the inner loop lacked posterior
covariance or abstention. Those implementation defects were corrected before
Wave 1 and the new CUDA invariants pass, but the corrected Bayesian model still
failed its predictive controls. GPU telemetry from the E0 run was
mean/peak utilization 16.7/44%, mean/peak draw 14.10/17.51 W, and peak memory
5.71 GiB. Windows Task Manager's 3D graph is not a CUDA-compute measurement.

## AnchorDelta History and Failure Reasons

- 2026-07-31: The strict k=5 Bayesian/ICSA route was reviewed and stopped
  below calibration. The governing information bottleneck is that removing
  the intercept and calibration direction leaves at most three independent
  support contrasts; the posterior/gate/kernel cannot create missing
  target-specific ranking information.
- 2026-07-31: AnchorDelta was implemented as a separate testable route. Its
  operator is the exact antisymmetric comparator
  `Delta(p,q,i) = (h(p,q,i) - h(p,i,q)) / 2`, and support labels are used only
  as absolute anchors. The frozen-feature P0 runner uses a fresh TRAIN
  homology-component holdout and target-balanced within-target pairs.
- 2026-07-31: Frozen-feature smoke and precheck runs showed that uniform anchor
  aggregation is invariant to a pure support-label column permutation. This is
  expected algebraically, so the valid label control is wrong-target support
  labels; that control changed absolute RMSE but did not change ranking.
- 2026-07-31: Joint retraining of the interaction encoder and comparator
  improved the single-seed precheck over calibration (RMSE 1.0775 vs 1.2402,
  Spearman 0.2521 vs -0.0026, pairwise 0.5886 vs 0.4872; 23 independent
  components). However, the wrong-protein arm was indistinguishable from the
  correct-protein arm (RMSE 1.0775, Spearman 0.2522, pairwise 0.5887), and all
  correct-minus-wrong-protein component-bootstrap intervals crossed zero.
- Failure reason: the retrained operator learned a transferable ligand/pair
  shortcut or anchor-offset effect, not protein-conditioned affinity change.
  More epochs, a wider head, Mamba/Transformer changes, graph ligands,
  uncertainty heads, or AdaMBind scheduling would not address this diagnosed
  absence of protein-specific signal.
- Decision: `NO_GO_FOR_PROTEIN_CONDITIONED_ANCHORDELTA`. Do not claim a
  protein-conditioned improvement. Any future continuation requires a new,
  independently validated protein-interaction information source or a changed
  data protocol; it must pass wrong-protein, wrong-support, similarity-bin,
  multi-seed, and component-bootstrap gates before architecture expansion.
- 2026-07-31: The TRAIN-only exact-ligand interaction audit found 17,120,772
  document-local rectangles but only 39,166 target-pair/document units. The
  homology-pair bootstrap reversal fraction was 0.385 (pKi 0.358, pKd 0.398),
  while strict same-assay rectangles were zero. pKi unit median absolute
  difference-in-differences was 0.20. The result is
  `AUDIT_RECTANGLES_EXIST_STRICT_ASSAY_COMPARABILITY_LIMITED`, not evidence
  that the current protein representation is identifiable.

## Permanent Firewall State

affinity_training_authorized = false
confirmation_authorized = false
sealed_authorized = false
new_affinity_labels_read_this_round = false
davis_target_confirmation_consumed = false
sealed_test_consumed = false

## Compact Verdict Table

| Date | Verdict |
| --- | --- |
| 2026-07-30 | A1 C1 stopped at frozen coverage gate. |
| 2026-07-30 | Strict k=5 audit: pKi viable for development, pKd secondary. |
| 2026-07-30 | MAML-style E0 did not exceed calibration. |
| 2026-07-30 | Program redefined as an independent Bayesian few-shot model. |
| 2026-07-30 | Corrected Bayesian Wave 1 stopped below ligand-only and gradient controls. |
| 2026-07-31 | AnchorDelta frozen-feature P0 and trainable-encoder precheck failed protein-specificity gate. |

## AdaMBind Literature and Innovation Gate (2026-07-31)

- AdaMBind was reviewed as an external reference before execution. Its new
  component is a query-loss/gradient-similarity task scheduler around a
  MAML-style molecular GNN plus protein 1-D CNN; it is not a new
  protein-ligand interaction operator and cannot create missing target-specific
  ligand reordering.
- MAML, TADAM, Deep Kernel Transfer, and MetaDTA were also reviewed. They
  change adaptation, metric conditioning, covariance transfer, or support
  aggregation, but all require a transferable task signal. The TRAIN-only
  interaction audit and the residual bilinear probe did not establish such a
  signal (protein-free RMSE 1.3550 versus correct-protein RMSE 1.4266;
  component-bootstrap gain -0.0709 [-0.1221, -0.0200]).
- Innovation decision: do not claim a high-performance improvement from
  AdaMBind scheduling or another meta-architecture on the current FORT data.
  Any future architecture expansion requires a new crossed, independently
  validated interaction-information source and correct-versus-wrong-protein
  component-level gates.

## AdaMBind Data Provenance and Split Audit (2026-07-31)

- External data were kept outside the FORT registry under
  `D:\FORT\tmp\adambind-data`; no confirmation or sealed labels were read.
- Training CSVs were downloaded from the pinned AdaMBind GitHub commit
  `01a169a6d62fba0d6c003f47bfba539e55f5b344` and verified against GitHub API
  file sizes/blob IDs. Local SHA-256 fingerprints are:
  `bindingdb-full-data.csv` (42,203 rows, 1,088 targets, 9,862 compounds,
  32,865,741 bytes, `3ebd8dfabd2a20c0dbceba35cc59ba8e6dd44a90798667fef2c9059bab63fbba`),
  `davis-full-data.csv` (30,056 rows, 379 targets, 68 compounds, 25,810,493
  bytes, `dc9331894d5eafa46787632cc0d9754406e5a96eb87980b27d4abe22308a6994e`),
  and `kiba-full-data.csv` (118,254 rows, 229 targets, 2,068 compounds,
  94,281,374 bytes, `7b1e306a2344e38c5d5bbcda6f6112201440bbaa92d5081a4fc054ed83edca24`).
- The paper supplement was downloaded from Figshare DOI
  `10.6084/m9.figshare.30963823.v1`, file `61860844` (`source_data.xlsx`,
  3,981,306 bytes, MD5 `a5e2b3f5d754d169c063f6ecd61b3108`, SHA-256
  `1b73ef01a34578d543070caf0a724c65dcd4d397022c3741864efeaed52cb0ac`, CC
  BY 4.0). The GitHub repository has no declared license; upstream data terms
  must be retained.
- Data integrity findings: BindingDB has two `Infinity` affinity values;
  Davis has 4,284 repeated compound-target rows; KIBA has 597. These are not
  independent observations and the non-finite values require an explicit
  policy before full-data loss/evaluation.
- No official `*_1.txt`, `*_2.txt`, or `*_3.txt` split files were present and
  `cd-hit` was unavailable. `research/adambind_data_audit.py` therefore
  implemented the same target-level random 80/10/10 fallback semantics with
  an independently deterministic NumPy seed 168. Target overlap is zero by
  construction, but exact
  compounds cross splits (BindingDB 554/730/295 train-val/train-test/val-test;
  Davis 68/68/68; KIBA 1,998/1,576/1,506). The per-target `nums=5` support/query
  construction exposes exact pair overlap in Davis (88) and KIBA (18), so it
  must not be interpreted as independent evidence. Full JSON:
  `reports/active/adambind_split_audit_2026-07-31.json`.
- Full provenance report: `reports/active/adambind_data_audit_2026-07-31.md`.

## AdaMBind Execution Ledger and Failure Reasons (2026-07-31)

Each failure below records the command/operation, observed error, root cause,
repair, and post-repair scientific usability.

1. GitHub metadata lookup using `GET /repos/Moohyun-w/AdaMBind/git/trees/01a169...`
   returned HTTP 404. The endpoint requires a tree SHA, not a commit SHA.
   Querying the commit API followed by `/contents/data?ref=<commit>` fixed the
   lookup; file sizes and Git blob IDs then matched. This was read-only and did
   not affect data.
2. An initial raw-download command contained a manually mistyped commit hash
   and returned `404 Not Found`. The URL was corrected to the pinned full hash;
   all three CSV downloads then completed and their SHA-256 values matched the
   API records. The failed request has no scientific use.
3. The unmodified command `D:\anaconda\envs\drug\python.exe create_data.py`
   failed with `FileNotFoundError: /AdaMBind/data//davis-full-data.csv`.
   AdaMBind hard-coded a Linux path. A local run patch added `--root`,
   `--dataset`, and `--max-rows`, and removed the processed filename's trailing
   quote; model architecture and loss code were unchanged.
4. The first patched preprocessing invocation used relative `tmp\\adambind-data`
   paths while the working directory was `tmp\\adambind-source`; PowerShell
   `Resolve-Path` and `--root` failed. Absolute paths from `D:\FORT` fixed the
   invocation. No raw file was changed.
5. After graph construction, loading the generated `.pt` failed with
   `UnpicklingError: ... DataEdgeAttr` because PyTorch 2.6 defaults to
   `weights_only=True`. The trusted local `TestbedDataset` loader was patched
   to use `weights_only=False`; the 2,048-row Davis processed file then loaded
   successfully with graph, target, and label tensors.
6. The first smoke invocation stopped with
   `ValueError: processed rows 2048 != CSV rows 30056`. The smoke directory
   held a full-CSV hardlink but a 2,048-row processed cache. The generated
   smoke CSV was replaced by an independent 2,048-row file and its cache was
   regenerated; the original downloaded CSV was preserved.
7. The attempted CSV repair repeated the same relative-path mistake from the
   source working directory, producing `Resolve-Path`/`Import-Csv` failures.
   The repair was rerun from `D:\FORT` with absolute paths and preprocessing
   completed.
8. The first smoke run reached CUDA but telemetry failed with
   `RuntimeError: Invalid device argument` from
   `torch.cuda.reset_peak_memory_stats(torch.device('cuda:0'))`. Using the
   current CUDA device with the no-argument API fixed telemetry; CUDA training
   then proceeded.
9. The next smoke run failed in AdaMBind `Trainer.train` with
   `UnboundLocalError: local variable 'y'` when `args.noise=0`. The source only
   assigned `y` inside the noise-enabled branch. Initializing `y=data.y` in
   both training paths fixed the deterministic no-noise path; default noise
   behavior was unchanged.

## AdaMBind CUDA Smoke Results (2026-07-31)

- The command pattern was
  `D:\anaconda\envs\drug\python.exe -m research.adambind_smoke` with the
  pinned source, fixed seed 168, official `GAT_GCN`, `Trainer`, and
  `Scheduler`. The runner writes a failure JSON with traceback before raising
  and labels successful outputs as mechanism smoke only.
- Davis: 2,048-row isolated subset, 2 train/1 validation/1 test targets,
  support/query 3/8. CUDA RTX 4060 Laptop GPU; validation RMSE 2.80865, CI
  0.61538, wall 0.985 s, peak allocated 293,755,904 bytes.
- BindingDB: 4,096-row isolated subset (the subset contained no non-finite
  labels), 2/1/1 targets, support/query 1/2. Validation RMSE 1.85753, CI
  1.0, wall 0.742 s, peak allocated 290,421,248 bytes.
- KIBA: 2,048-row isolated subset, 2/1/1 targets, support/query 3/5.
  Validation RMSE 0.47209, CI 0.83333, wall 0.949 s, peak allocated
  288,909,312 bytes.
- All three runs exercised task scoring, task sampling, one support update,
  validation prediction, and the scheduler policy-gradient update. These are
  runtime checks only: the subsets are not strict FORT target/homology/
  scaffold/document/assay-cold benchmarks, and the metrics cannot support a
  claim of superior performance or innovation.

## Verification and Boundary

- New scripts: `research/adambind_data_audit.py` and
  `research/adambind_smoke.py`. Local compatibility edits are confined to the
  downloaded AdaMBind checkout under `tmp/adambind-source`.
- `D:\anaconda\envs\drug\python.exe -m py_compile` passed for the new scripts
  and patched AdaMBind files. Full repository tests passed: `66 passed`.
- The external data and smoke results remain separate from FORT training,
  confirmation, and sealed evidence. They do not overturn the existing
  `NO_GO_FOR_PROTEIN_CONDITIONED_ANCHORDELTA` or interaction-identifiability
  gate, and no AdaMBind-derived performance claim is authorized.

## HTL-DTA Task Migration and Topology Audit (2026-07-31)

- The primary estimand was migrated from strict unseen-target `k=5` adaptation
  to `HTL-DTA`: abundant-to-scarce transfer from data-rich head targets to
  scarce recipient targets, with separate pseudo-tail and natural-tail roles,
  target-macro transfer-gain curves, and an explicit negative-transfer gate.
  The strict task remains a secondary stress test.
- The migration report is
  `reports/active/htl_dta_task_migration_2026-07-31.md`; the control-plane state
  is `HTL_DTA / HTL_1_TARGET_FREQUENCY_TOPOLOGY_AUDIT` in
  `manifests/state.v1.json`.
- A metadata-only audit, `main.py topology-audit --split train`, generated
  `dataset/processed/htl_target_topology.v1.json`. It read only registry
  metadata columns through the TRAIN split filter; affinity and replicate
  labels were not loaded, and no model training ran.
- The audit found 559 pKi and 407 pKd TRAIN targets. Under the candidate rule
  `n_eff < 10`, only 5 pKi targets have at least six unique scaffolds and six
  provenance units, so that threshold cannot meet the primary tail-power floor.
  The unselected `n_eff < 30` rule has 115 pKi upper-bound candidates, but
  document/assay closure, natural-tail provenance, component-aware power, and
  threshold freezing are still required.
- For the provisional abundant-source/scarce-recipient topology (`n_eff >= 100`
  versus `<30`), pKi has 242/193 source/recipient targets and only 12
  same-homology recipients; pKd has 41/256 and only 5 same-homology recipients.
  Most candidate recipients are homology-cold to the source pool, so
  protein-similarity-only transfer is a control, not an assumed mechanism.
- Full software regression after the migration passed (`68 passed`). These
  results are topology evidence only and do not authorize HTL training.

## A2S-DTA Baseline Execution (2026-07-31, SUPERSEDED)

The metrics in this section are retained for audit history only. They were
generated before the global registry-to-feature alignment correction.

- The first target-disjoint A2S source-only baseline ran on CUDA with pKi as
  primary and pKd as secondary. pKi produced 151/137/117 recipient episodes at
  k=1/3/5; pKd produced 119/101/87.
- On pKi, support-compatible source routing gained `+0.034 [-0.004,+0.075]`
  at k=1, then `-0.020 [-0.058,+0.012]` and `-0.013 [-0.044,+0.021]` at k=3/5
  versus recipient calibration. Protein, chemistry, and random routing were
  negative. The scalar adapter is not a universal transfer solution.
- On pKd, support-compatible routing was positive at all budgets
  (`+0.286`, `+0.163`, `+0.132`, all with positive target-bootstrap intervals),
  while random/protein/chemistry controls were mostly negative. pKd remains a
  secondary endpoint and is not pooled with pKi.
- The support-evidence gate reduced some pKi harm but did not create positive
  pKi transfer. Decision: A2S data chain PASS; scalar adapter and current gate
  are diagnostic only. A cross-fitted recipient-conditioned router is the next
  admissible model, with no Mamba/Transformer expansion.

## A2S Alignment Correction And Router Decision (2026-07-31)

- A P0 audit found that the initial A2S runs assigned `source_row` after
  endpoint/split filtering. This local index was incompatible with the global
  `ligand_features.npz` row order. The initial pKi/pKd artifacts are therefore
  quarantined and must not be cited.
- The corrected loader reads the full 343,211-row registry, preserves its
  global row id, verifies feature length and `conn_sha`, and only then filters
  TRAIN rows. The regression suite passed 72 tests.
- The primary corrected ridge gives each source target equal total weight.
  Target-balanced source-support routing is positive on pKi at k=1/3/5:
  `+0.137 [0.070,0.207]`, `+0.112 [0.054,0.174]`, and
  `+0.087 [0.029,0.149]`. Corrected pKd is a separate positive replication
  (`+0.373`, `+0.220`, `+0.242`), never pooled with pKi.
- The target-balanced source-fold cross-fitted router used 726 pseudo-recipient
  episodes but was strongly negative on pKi at k=3/5 (`-1.213` and `-1.315`
  RMSE gain versus recipient calibration). The router route is NO-GO; neural
  architecture training remains blocked pending frozen support draws and a
  genuinely provenance/time-closed natural-tail roster.
- The corrected artifacts now include NDCG@10, benefiting-recipient rate, and
  RMSE-gain AULC. The pKi router AULC is `-1.221`, with mean benefiting rate
  about `0.201`; this confirms the negative-transfer stop is not an artifact of
  a missing aggregate metric.
- The complete A2S-CFRA model specification, literature novelty matrix,
  cross-fitting contract, kill tests, strong-baseline matrix, and final
  judgment are frozen in
  `reports/active/a2s_dta_master_design_2026-07-31.md`.

## A2S Natural-Tail Gate D0 Data Stop (2026-07-31)

- The metadata-only command `main.py natural-tail-audit` read no affinity,
  replicate SD, development, confirmation, or sealed outcomes and performed no
  model training.
- Among 193 candidate pKi recipients, the most optimistic five-draw
  parent/document/assay-closed upper bound is 40 targets; adding distinct
  support scaffolds leaves 34. Both are below the frozen minimum of 50 before
  time, source-family, common-query, or query-scaffold constraints.
- ChEMBL-release temporal ordering leaves at most 22 targets and full temporal
  scaffold/provenance closure leaves 11. A held-out document `src_id` has one
  envelope target and zero fully closed targets, so the strict roster and
  independent component count are both zero and MDE80 is undefined.
- All 193 candidates join into one dependency component through shared source
  family or homology, and the local metadata lack true publication or
  measurement dates. Final decision: `DATA_NOT_READY`; S0, A2S-MAP, MAML,
  AdaMBind, five-seed training, and encoder expansion remain blocked.
- Post-run adversarial review limits the JSON to a feasibility-stop
  certificate, not an admission-grade strict roster: the strict diagnostic
  adds scaffold-cold closure, can repeat support sets across draws, uses
  provisional resource thresholds, and does not recompute emitted overlap.
  The stop remains supported by the separate exhaustive five-distinct-support
  upper bound of 40 and the absence of true time/lineage metadata.

## Ready package separation (2026-08-01)

- Consolidated preprocessing is now `scripts/preprocess.py`; its sealed run
  produced `dataset/processed/a2s_validation_small.v1/`.
- The verified model-facing package is mirrored into
  `dataset/ready/a2s_validation_small.v1/`, indexed by
  `dataset/ready/manifest.json`.
- This is a storage and input-contract update only. Formal historical affinity
  training remains `DATA_NOT_READY`: 12,782 identity collisions remain and
  zero historical affinity values were materialized. No training was run.

## A2S-CMAL Implementation And Failure Ledger (2026-08-01)

| Attempt | Observation | Root cause | Repair / scientific status |
| --- | --- | --- | --- |
| `python -m research.a2s_cmal_data --device cuda`, first run | Failed before creating the output directory with a missing `(target_uid, compound_parent_uid)` lookup for one recipient query compound. | `_first_observation` deduplicated a multi-target table by compound alone. A compound measured against more than one target therefore retained only the first target's row. | Changed the identity key to `(target_uid, compound_parent_uid)` and added a regression test with one compound shared by two targets. No label column was read and no package or scientific result was emitted; the failed attempt is unusable except as preprocessing-key evidence. |
| `pytest -q tests/test_a2s_cmal.py tests/core/test_a2s_cmal_data.py`, first run | Four model-construction tests failed before any forward pass; the five data/objective tests passed. | The helper method `_parameters` shadowed `torch.nn.Module._parameters`, which must remain a parameter dictionary. Calling `.eval()` therefore tried to call that dictionary. | Renamed the helper to `_collect_parameters`; no dataset, checkpoint, metric, or scientific conclusion was produced by this failed run. |
| Intentional local `main.py a2s-cmal --formal` guard test | Exited non-zero before data loading with `formal training is external-only`. | Expected policy gate: formal recipient-label access requires `A2S_FORMAL_EXTERNAL=1` on the designated host. | Guard is working. This is a negative execution test, not a training failure; no recipient label, model step, checkpoint, or metric was produced. |
| First 20+20-step CUDA smoke | Training and source-validation evaluation completed, but metric aggregation emitted `RuntimeWarning: Mean of empty slice`; the 0.9-second training window also yielded only five low-utilization telemetry samples and was too short for a GPU conclusion. | The k=1 label-permutation control is deliberately undefined and represented by `NaN`; direct `np.nanmean` warned when every value in an arm/metric slice was undefined. Telemetry at 0.5-second cadence cannot characterize a sub-second phase. | Replaced direct `nanmean` with an explicit finite-value mean that preserves undefined results without warning. The smoke remains a mechanism/entry-point check only; a longer source-only profile is required before reporting utilization. |

### Successful controls and non-results

- `research.a2s_cmal_data` initially produced the label-blind
  `a2s_cmal_episodes.v1` package. It passed its then-defined identifier audit,
  but is now **superseded and forbidden for formal training** because it froze
  target-parent rather than `measurement_uid`; see the failure entry below.
- The corrected targeted suite passed 9 tests. It verifies support-free
  invariance, protein/support/query dependence, frozen base parameters, the
  post-adaptation ranking contrast, component-safe splits, and label-blind
  counterfactual construction. Unit tests do not establish scientific gain.
- A source-only v2 250+250-step CUDA profile at batch 64 completed without metric
  warnings: base 9,341.9 episodes/s; adapter 5,659.3 episodes/s; telemetry
  mean 32.82%, P90 56.0%, and 54.55% of samples at >=40% utilization. All
  episode tensors were GPU-resident, four counterfactual arms were fused, and
  no pandas lookup or host transfer occurred in the training loop. The report
  is `MECHANISM_SMOKE_ONLY`; no recipient label was read and its metrics must
  not enter a paper table.
- No formal five-seed recipient run has been attempted on this device.
- A documentation audit initially tried to read
  `dataset/ready/DATASET_RECORD_SUMMARY.md` and received `PathNotFound`; the
  existing record is at repository root (`DATASET_RECORD_SUMMARY.md`). The
  corrected path was used. This failure changed no data or code.

### Measurement-identity correction

| Attempt | Observation | Root cause | Repair / scientific status |
| --- | --- | --- | --- |
| Reverse audit of the v1 episode-to-label join before final verification | Among 87,451 episode target-parent pairs, 4,052 mapped to more than one admitted assay context and 587 had no row in `pki_measurements_context_main.parquet`. The trainer instead grouped all exact records for a target-parent. | v1 froze parent identity but not the actual observation. Joining by target-parent could aggregate across assays and years, changing the document-ordered estimand and re-admitting high-noise contexts. | **STOP v1 for formal use.** No recipient label was read during diagnosis. Rebuilt the metadata-only roster as `a2s_d0r_roster.v3` directly from context-main metadata and froze one deterministic `measurement_uid` per support/query item; rebuilt episodes as v2 and changed the trainer to join labels only by that UID. |
| `main.py d0r-roster --out ...v3` and `main.py prepare-cmal-data --output ...v2` | Both completed: roster PASS with 206 sources, 63 recipients and 55 components; episode package has 30,123 episodes and content SHA-256 `d8e9a259f594db293c4e46779cb716b8e52a828a989ae20cdbc5571805877f9b`. | Corrected label-blind construction, not a rescue fitted to outcomes. | All target/accession/document/parent/assay and support/query measurement overlaps are zero. Builder requested no pKi column. v2 is the sole formal CMAL input; no formal training result exists. |

### Counterfactual mechanism diagnosis after the gradient audit

The preceding v2 episode package was subsequently superseded by immutable
`a2s_cmal_episodes.v3`, content SHA-256
`2df5831bc8a51df93dc54531302327716fcca8900ec43f1aa37f16ed2fb9485a`.
v3 keeps the measurement-identified rows and changes only the query-blind
chemical hard-negative rule to maximize support-set Murcko-scaffold Jaccard,
with support-centroid ECFP4 cosine as the tie-break. All leakage counts remain
zero. v1 and v2 remain on disk for provenance but are forbidden as current
CMAL inputs.

| Attempt | Observation | Root cause learned | Repair / scientific status |
| --- | --- | --- | --- |
| Original v3 gradient audit, base 300 + adapter 500, seed 1729 | Counterfactual loss fell from about 1.386 to 1.068 and the train ranking gap reached +0.2229. All adapter modules had nonzero gradients and 10%-32% relative parameter changes. Source validation nevertheless changed by RMSE +0.09846, CI -0.00988, Spearman -0.03085 and NDCG@10 -0.01262; correct support did not beat any wrong-support arm on CI. | The failure was scientific rather than mechanical. The relative contrast could grow without an absolute unseen-target improvement. | **STOP.** No recipient label was read. Added a source mechanism gate and prohibited recipient-label access/checkpoint emission when it fails. |
| Original operator at adapter step 100, with and without counterfactual loss | Both branches already reduced absolute source-validation ranking. Removing the counterfactual term made all three specificity point estimates positive but still left CI/Spearman/NDCG below the frozen base. | The delta head could act as a second query/protein predictor without being algebraically anchored to measured support residuals; late 500-step overfitting and the counterfactual term alone were both insufficient explanations. | Multiplied the learned query-specific scale by the attention-weighted measured residual. Added a regression test requiring exactly zero delta when support labels equal frozen-base support predictions. |
| Residual-anchored operator, adapter step 100 | Validation RMSE improved by 0.0138 and CI became +0.00033 relative to base, but Spearman/NDCG remained negative and protein-hard CI specificity failed. | Removing the query-only bypass repaired one causal defect but did not make the relative contrast safe or transferable. | Retained the measurement anchor; did not call this a breakthrough. |
| Frozen-base-anchored capped-gain InfoNCE, steps 100 and 300 | The objective no longer profited after a wrong arm became worse than base. At step 100 absolute CI/Spearman were only slightly positive, NDCG was negative and all wrong-support specificities were negative; step 300 worsened. | Raw-loss InfoNCE had a destructive-wrong-arm shortcut, but removing it exposed a second chemistry/task-recognition shortcut and weak positive adaptation. | Retained the base-anchored scoring rule because it removes a provable bad incentive; stopped this branch at step 300. |
| Same-compound label-swap counterfactual, step 100 | All validation ranking and four-arm specificity point estimates became positive, but nearly all component-bootstrap intervals crossed zero. A newly evaluated source meta-test changed by CI -0.00036, Spearman -0.00110 and NDCG -0.00020. | A chemistry-only classifier could identify the correct original arm at 51.6% train/54.0% validation versus 25% chance. Label swap removed that shortcut, but the small validation direction did not transfer. | **FAIL, evidence insufficient.** Step 300 overfit both gates. Recipient labels remained sealed. |
| Component-disjoint base/meta-adapter targets with target-balanced sampling, step 100 | Validation CI/Spearman improved +0.00457/+0.00888 but NDCG fell; source holdout CI/Spearman/NDCG all fell. | Base-seen versus unseen residual scale was a real covariate shift (support-residual mean SD +43.5%; query residual SD +19.7%) but not the only cause. | Continued only to the preregistered 300-step diagnostic checkpoint. |
| Component-disjoint branch, step 300 | Source holdout CI/Spearman improved +0.00303/+0.00898 and all four correct-support specificity comparisons were positive for those metrics, but holdout NDCG fell -0.00122 and RMSE worsened 1.6658 to 1.7231. Source validation simultaneously fell CI -0.01011, Spearman -0.03235 and NDCG -0.01551. | The operator learned relative support identity on one split but not a stable beneficial correction. Sign reversal across source splits and top-ranking failure rule out a transferable-mechanism claim. | **STOP; no key positive breakthrough.** Created `reports/active/CMAL_FAILURE_HANDOFF.md` and an English external-agent prompt. No formal run, recipient label access, model-folder promotion, commit, push or GitHub publication occurred. |

## Source Information-Gate Preflight And Probe (2026-08-01)

| Attempt | Observation | Root cause / interpretation | Decision |
| --- | --- | --- | --- |
| Metadata-only source lock | `research/a2s_source_lock.py` created `a2s_source_information_gate_lock_2026-08-01.json`; 559 TRAIN pKi targets became 141 homology+document/assay provenance components with zero cross-role overlap. The largest component contains 380 targets (68.0%). | **FACT:** the provenance closure is internally disjoint but highly concentrated. **INFERENCE:** independent component power is materially smaller than the target count. | Keep as a new locked split; treat the large component as a power limitation. |
| First source-gate execution | OOF base completed and produced `a2s_source_information_gate_oof_2026-08-01.npz`, but episode construction repeatedly rebuilt an 183k-row source-row index for every query. The process was stopped before a result was emitted. | **FACT:** diagnostic implementation performance defect; no scientific result was produced and no locked labels were read. | Fixed by caching the index and moving synthetic feature construction outside the query loop. |
| Source-only G0/G1 probe | Fit/probe labels only; 858 episodes; probe components 10/7/4 at k=1/3/5. Real `Delta_label` ranking-loss means were -0.00123, -0.00358, and -0.00076. Real `Delta_assign` was undefined at k=1, -0.00443 at k=3, and +0.00132 at k=5; all component intervals crossed zero. | **FACT:** the registered information gate failed. **INFERENCE:** the present global source episode distribution does not show stable assignment-specific support-label information in this bounded probe. The result is not an information-theoretic impossibility theorem. | **STOP** before opening locked labels or implementing STOP/SWAP/CSRIO. |
| OOF balance audit | The 380-target provenance component occupies one OOF fold: 176,193 held-out rows versus 5,382 training rows; 97.0% of fit-role held-out rows are in that fold. | **FACT:** the component exclusion contract is satisfied but the source residual geometry is severely unbalanced. **INFERENCE:** the real G0/G1 result is non-confirmatory until a provenance design with adequate fold balance is available. | Do not call the null definitive; require an independent label-free split design or new provenance-rich data. |
| Label-free topology comparison | Homology-only closure has 517 components with maximum size 4. Adding exact document cells gives 163 components with a 347-target (62.1%) giant component; assay cells do not connect targets. | **FACT:** the giant component is driven by the available document provenance, not solely by pipe-token handling. **INFERENCE:** a balanced provenance-closed split likely requires richer source-family/campaign metadata. | Keep the current lock for audit provenance but require a new metadata source before treating a null as confirmatory. |
| Synthetic positive control | The same probe recovered an injected query-dependent support signal with rank-loss gains +0.4181/+0.4219/+0.2619 at k=1/3/5. | **FACT:** the probe has basic signal-detection power after isolating the five label-channel coordinates. | The real null is not explained solely by a disconnected label channel, but k=5 still has only four independent components. |

## Balanced V2 Source Information Gate (2026-08-01)

| Attempt | Observation | Root cause / interpretation | Decision |
| --- | --- | --- | --- |
| First v2 gate launch | Stopped before model fitting because target-filtered source labels still contained rows removed by the v2 quarantine, while the identity table contained retained rows only. | The join direction was inherited from the target-level v1 lock. A row-level lock requires the retained identity table to be the left side of the join. | Reversed the join, added a retained-row hash check, and reran. No locked or recipient labels were requested. |
| Homology-first, row-quarantine v2 lock | Retained 68,782/185,591 source rows. Fit/probe/locked contain 222/110/107 homology components. Cross-role target, homology, document, and assay overlap is zero. OOF fold loads are 7,322/7,322/7,322/7,321 rows. | The v1 97% fold collapse was caused by unioning document-connected targets. Quarantining cross-role provenance rows restores component power but changes the estimand and discards 62.9% of rows. | Accept as a label-free development split; keep the old lock as failure provenance. |
| Corrected nested G0/G1 v2 gate | G0 and G1 were fitted separately at equal capacity, with the label channel masked during both G0 training and evaluation. Delta_label 95% lower bounds were negative at k=1/3/5; Delta_assign lower bounds were negative at k=3/5. Synthetic controls remained strongly positive. | Balanced evidence still does not admit stable incremental information from the current passive support labels. High-data oracle headroom is positive, so this is not evidence that query ranking is intrinsically unimprovable. | `NO_GO_INFORMATION_NOT_ADMITTED`. Do not implement STOP/SWAP, open locked/recipient labels, promote `model/`, commit, or push. The only admissible next study is a label-free same-assay/MMP coverage and power census. |

## A2S-TRACE Stratum Gate And Mechanism (2026-08-01)

| Attempt | Observation | Root cause / interpretation | Decision |
| --- | --- | --- | --- |
| Q1 stratum-resolved information gate (`research/a2s_trace_stratum.py`) | One corpus, one frozen base, one fixed Tanimoto KRR, one residual derangement, 12,246 probe episodes over 110 components. Varying only the support policy and the support-query nearest-Tanimoto stratum: KRR-minus-base CI lower bounds are negative below 0.35 nearest Tanimoto in every policy and +0.023 to +0.048 at nearest Tanimoto >= 0.55. Correct-minus-deranged and correct-minus-norm-matched-wrong-target switch on and off in exactly the same cells. The level channel scored exactly 0.0000 CI in all 45 cells. | **FACT:** support-label information is a property of the support-to-query chemical relation, not of the corpus. The v2 `provenance_disjoint` policy draws queries at mean nearest-Tanimoto 0.19-0.30, i.e. inside the measured null bins, so its global null was a correct measurement of an uninformative stratum; the BindingDB positive was a correct measurement of an informative one. **INFERENCE:** constraint C9 is resolved without overturning either prior result. | `INFORMATION_ADMITTED_IN_A_LOCAL_RELATION_STRATUM`. All later work runs inside the admitted stratum and is judged against fixed Tanimoto KRR, not against the frozen base. Retire `provenance_disjoint` as a training substrate; keep it as the off-stratum control. |
| Transport-scale diagnostic | Multiplying the KRR transport by a constant raises admitted-stratum k=5 CI from 0.5731 to 0.5951 at scale 3. Within-episode SD is 0.96 for the base and 0.46 for the transport, while the base orders at chance. | **FACT:** the frozen base is over-weighted relative to support evidence, and the relative precision of the two channels is a free scalar. **INFERENCE:** any learned gain measured against unscaled KRR would be mostly a scale artefact. | Grant a single meta-learned global scale to the analytic bar (rung R2c) and to the static-mixture baseline. The mechanism claim is measured strictly on top of it. |
| Ranking-surrogate diagnostic | Trained with the convex RankNet logistic loss, the global scale converges to 1.03 and validation CI does not improve; with a bounded smoothed-CI surrogate it converges to ~1.5 and CI improves. | **FACT:** the convex surrogate keeps paying for confidently-wrong pairs, so its optimum sits far below the CI optimum. **INFERENCE:** this is a property of the substrate (chance-level base ordering, large base spread), not a hyperparameter preference. | Use the bounded smoothed-CI surrogate for all ranking work on this substrate; record the mismatch. |
| Q2 TRACE mechanism (`research/a2s_trace.py`) | Learned label-free per-pair transport reliability plus a per-query gate, zero target-specific parameters, bounded by `max_i \|r_i\|`, exactly nesting fixed KRR and Nadaraya-Watson. 19,611 fit / 4,796 inner-validation / 4,029 probe episodes, three seeds, 74-76 components. Global scale gain over fixed KRR: +0.0104 [+0.0042, +0.0169] CI at k=3 and +0.0089 [+0.0027, +0.0149] at k=5. Learned gain over that bar: -0.00015 [-0.00063, +0.00030] at k=3 and -0.00003 [-0.00059, +0.00053] at k=5. Identical null for the low-capacity variant, the protein channel, and against the static kernel mixture; CKA-NNLS loses by 0.017-0.023. Residual-null is a bitwise no-op; correct-minus-deranged is +0.090 and correct-minus-wrong-target +0.060 to +0.083. | **FACT:** all seven admissibility conditions hold and the headline gates M1/M1b fail. **INFERENCE:** per-pair transport reliability adds nothing beyond isotropic Tanimoto in this construction. | `POSITIVELY_CONTROLLED_NULL_LEARNED_TRANSPORT_NOT_ADMITTED`. Do not open `locked` for TRACE; there is nothing to confirm. |
| Synthetic positive control | In an injected world where transport reliability really is pair-dependent, the identical learner and episodes recover +0.0262 [+0.0154, +0.0370] CI at k=3 and +0.0159 [+0.0078, +0.0238] at k=5 over the same bar, i.e. 42% and 27% of the oracle gap. | **FACT:** the pipeline detects an effect roughly 30x larger than the real-data 95% upper bound of +0.0005. | The real null is a measurement, not a failure to look. Report it with its upper bound. |
| Headroom oracles (`research/a2s_trace_headroom.py`) | Hindsight ceilings over fixed KRR on the same probe episodes at k=5: episode-level scale +0.078 to +0.107 CI, episode-level support subset +0.059 to +0.065, per-query support subset +0.353 to +0.366 (absolute CI 0.93-0.95). Mean hindsight-optimal episode scale 2.2-2.5 versus the single learned global scale ~1.5. | **INFERENCE:** the action class is expressive enough to reach near-perfect rankings; what is missing is any label-free predictor of which per-pair choice is right. The largest measured gap with real room is episode-level magnitude, which is the TAMSK claim - now with a ceiling and a bar that already contains the global scale. | Next study is a predictable per-episode transport scale from support-only label-free quantities. No `model/` promotion, no recipient label access, no commit or push occurred. |

## A2S-MODE Pre-Implementation Gates A0-A4 (2026-08-02)

| Attempt | Observation | Root cause / interpretation | Decision |
| --- | --- | --- | --- |
| Redirection to meta-adaptation | The user rejected similarity/reliability weighting as a final contribution and required a learned adaptation state `z_t = A_theta(S_t)` with a small target-specific intervention. Six candidate mechanisms were screened; discrete response-mode selection (A2S-MODE) was selected, with sparse SAR-transformation and continuous state filters rejected on measured grounds. | **INFERENCE:** every mechanism in the record belongs to the class `Delta_q = sum_i w(x_q,x_i,p) r_i`, whose reach Q1 measured to be bounded by support-query chemical distance. A correction of the form `R_theta(z_t, x_q)` contains no support compound and is therefore not distance-limited. | Wrote `reports/active/A2S_MODE_MECHANISM_PROPOSAL_2026-08-02.md`. No mechanism code. Ran Gates A0-A4 first, per programme discipline. |
| Gate A0: per-target headroom by stratum | A per-target ridge head on a 26-dim label-free basis (10 descriptors + 16 Morgan PCs, fit-role statistics only), evaluated on held-out queries, beats the frozen base by +0.0517 [+0.0154, +0.0848] CI at nearest-Tanimoto < 0.20, +0.0774 [+0.0451, +0.1131] at 0.20-0.35, and +0.0847 [+0.0624, +0.1080] pooled. | **FACT:** target-specific ranking structure recoverable by a query-only function exists in every stratum, including the two where fixed KRR, TRACE and all transport operators measure exactly zero. **INFERENCE:** the first positive ranking headroom this programme has measured outside the local-analogue regime. | **PASS.** The A2S-MODE premise is confirmed. The number is an oracle ceiling (it uses the target's abundant labels), not a method. |
| Gate A1: mode sufficiency | k-means over 110 fit-target heads; split-half query selection (never hindsight on the scored queries). The selected mode beats the single global head by +0.0192 [+0.0020, +0.0378] at k=3 / t<0.20, +0.0227 [+0.0045, +0.0410] at k=5 / t>=0.55, and +0.0271/+0.0301 pooled. `M=2` fails; `M` in {3,4,6} pass. | **FACT:** a small discrete set of response modes carries a material part of the A0 headroom. **INFERENCE:** it is not explained by "a better global ligand model" - the global head alone sits at or below the frozen base. | **PASS.** A discrete target state exists and is sufficient. |
| Gate A2: k-shot identifiability | Mode selection from k support residuals reaches accuracy 0.16-0.33 against chance 0.20; every gain lower bound sits below the 0.005 MDE (best: +0.0196 [+0.0036, +0.0373] at k=5 / t>=0.55). The verdict also flips with the nuisance choice `M` (M=3 passes, M in {2,4,6} fail). | **FACT:** the strict conjunction of above-chance accuracy and an MDE-clearing gain does not hold. | **FAIL,** and not interpretable on its own - see Gate A4. |
| Gate A4: synthetic positive control | In a world where each probe target is generated from one dictionary mode plus the measured noise (`sigma`=0.976, level SD=1.898), k-shot selection reaches only 0.32 accuracy at k=5 and every gain interval crosses zero. | **INFERENCE - the decisive localisation.** A k-means dictionary is not separable from k<=5 noisy residuals even when the world is exactly the model. k-means minimises within-cluster variance, which is unrelated to between-mode discriminability `D_k` on a random k-subset. This is precisely the proposal's own argument: `rho_k` is fixed by the data, but `D_k` is an object the outer loop must maximise, and nothing in this gate maximised it. | The A2 failure is estimator-side, not biology-side. It does not refute the discrete-state premise that A0 and A1 confirm. |
| Consequence for the mechanism | The claim narrows to: meta-training a dictionary to maximise worst-case pairwise discriminability at the deployment budget makes a discrete target state k-shot-identifiable, where an unshaped dictionary of matched size, capacity and basis is not. | The unshaped k-means dictionary is the exact nested restriction the shaped model must beat, in the role `a2s_bir_global` played for IDA and `R2b_krr` played for TRACE. Gate A4 becomes the primary sanity rung: shaping must first make the synthetic control recoverable. | `PREMISE_CONFIRMED; K_SHOT_INFERENCE_NOT_YET_IDENTIFIABLE_WITH_AN_UNSHAPED_DICTIONARY`. Implementation authorised with six registered predictions. No `model/` promotion, no recipient label access, no commit or push occurred. |

## A2S-MODE Generalizability And Value Gates G1-G4 (2026-08-02)

| Attempt | Observation | Root cause / interpretation | Decision |
| --- | --- | --- | --- |
| Correction to Gate A0 | A0 stratified the per-target head's gain by **support-query** Tanimoto, an axis a head fitted on ~130 of the target's own rows never sees. The claim "first ranking headroom outside the local-analogue regime" was therefore not established by A0. | **FACT:** mis-stratification. The A0/A1 numbers stand; the interpretation does not. | Withdrew the claim; added a supersession note to `A2S_MODE_GATES_A0_A4_DECISION_2026-08-02.md`. |
| Gate G1: within-target scaffold-disjoint head | A head fitted on one set of Murcko scaffolds and scored on scaffolds it never saw beats the frozen base by +0.0524 [+0.0293, +0.0746] CI over 50 components. Re-stratified by similarity to the head's own training rows: +0.0706 [+0.0425, +0.0986] at >=0.60, +0.0385 [+0.0003, +0.0751] at <0.30 (5 components, underpowered), intervals crossing zero in between. | **FACT:** the object survives scaffold-disjointness, so it is not within-target series memorisation. **INFERENCE:** its reach is longer than k<=5 support transport but not unbounded; "scaffold-cold" and "chemically distant" are different predicates and must stop being conflated. | **PASS,** with the range qualified. |
| Gate G2: intrinsic dimension | Spectrum of 110 source-target heads is nearly flat (top three directions hold 34.7% of variance). Projecting a probe target's own oracle head onto the top-r **source** subspace gives -0.0058 (r=1), -0.0033 (r=2), -0.0007 (r=3), +0.0006 (r=5) versus +0.0524 at full rank 26. | **FACT:** rank-2 projection retains -6% of the full gain. **INFERENCE - decisive:** the directions that matter for a target are essentially orthogonal to the directions along which source targets vary most. Target heads are high-dimensional and idiosyncratic, not draws from a small shared basis. | Refutes A2S-MODE's shared-mode premise (and explains the A2/A4 failure structurally), and refutes the A2S-IDA rank-m code family on this substrate - low rank is not merely unidentifiable but **actively harmful at every label budget**. |
| Gate G3: zero-shot protein prior | A ridge from pooled ESM-2 to head coefficients, fitted on fit targets and applied with no probe label read, scores -0.0185 [-0.0725, +0.0190] CI. | **FACT:** protein sequence does not predict a target's ligand-response head. Consistent with the earlier `TR group not resolvable` result. | **FAIL.** No protein-conditioned zero-shot shortcut exists; adaptation must come from support labels. |
| Gate G4: label learning curve | Empirical-Bayes head from k labels in the top-r source subspace (prior `lambda_j = sigma^2/tau_j^2` from source targets, never tuned on probe), scaffold-disjoint, versus frozen base at full rank: k=1 -0.0067, k=3 +0.0014, k=5 +0.0111 [-0.0028, +0.0245], k=10 +0.0261 [+0.0122, +0.0410], k=20 +0.0432 [+0.0252, +0.0613], k=40 +0.0517 [+0.0332, +0.0735] (the full oracle head is +0.0524). Ranks 1/2/3 are flat-to-negative at every budget. | **FACT:** the knee is at k ~ 10; the best k<=5 cell has a lower bound of -0.0028. The object is fully recovered by 40 labels. **INFERENCE:** the k<=5 deployment budget sits measurably below the knee of this object's learning curve, measured with a prior-regularised estimator that has no free parameters to blame. | **FAIL at k<=5.** `GENERALIZABLE_BUT_NOT_FEW_SHOT_REACHABLE`. |
| Value verdict | At k>=10 a closed-form empirical-Bayes ridge on a 26-dim label-free basis delivers +0.026 to +0.043 CI with no meta-learning of any kind. | **INFERENCE:** under this programme's admissibility rules that is a strong baseline, not a mechanism. The well-posed remaining question is whether meta-learning can move the learning curve left from k~10 to k~5, with G4 as the curve to beat. | Recommended next lever: meta-learn the **basis** so heads become low-rank in the learned representation (the one form of the IDA idea G2 does not refute), with the head spectrum and the rank-2 retained fraction as its own falsifiers. No `model/` promotion, no recipient label access, no commit or push occurred. |

## A2S-RIP Gate R0 And The A2S-HOTSPOT Branch (2026-08-02)

| Attempt | Observation | Root cause / interpretation | Decision |
| --- | --- | --- | --- |
| Literature derivation of A2S-RIP | Few-shot molecular/DTA methods (FS-CAP, MHNfs, PACIA, AdaMBind, ADKF-IFT, APN, CFS-HML) all learn a compact continuous task state and emit a dense correction for every query; none abstains, certifies, or reports a harm rate. Meta-learned conformalisation (Fisch et al. ICML 2021) makes a *threshold* identifiable from few labels; Selective Conformal Risk Control names ranking as future work; conformal selection in chemistry selects compounds to test, not interventions on a ranking. | **INFERENCE:** Gate G2 closes the compact-task-state family on this substrate, so the transferable object must move from the task function to the calibration of the estimator's own uncertainty. The action becomes a sparse, bounded, per-compound certified ranking edit. | Wrote `A2S_RIP_MECHANISM_DERIVATION_2026-08-02.md` with preregistered P1-P7, including a magnitude-matched wholesale control that no prior work runs because none measured the +0.009 CI global-scale artefact. |
| Gate R0a: selection ceiling | Hindsight-selected subsets of the empirical-Bayes head beat the wholesale head by +0.0621 [+0.0542, +0.0705] at k=3 and +0.0748 [+0.0626, +0.0886] at k=5, peaking near 40-50% coverage; +0.0735 versus the frozen base at k=5. The ceiling also beats a magnitude-matched wholesale head by +0.058/+0.067. | **FACT:** applying a noisy head to the right 40% of compounds is worth more than the entire fully-supervised head (+0.052 in G4). The k<=5 problem is *where* to apply the head, not whether it exists. | **PASS.** Ceiling real and not a magnitude artefact. |
| Gate R0b/R0d: is it reachable? | Margin AUC for predicting a correct edit is 0.5536 [0.5273, 0.5779] at k=3 and 0.5553 [0.5289, 0.5796] at k=5; `abs(delta)` alone gives 0.5507/0.5514, so the posterior covariance adds 0.004 AUC and is not load-bearing. The implementable margin rule gains -0.0006 [-0.0170] over base at k=5 and is indistinguishable from random selection at every coverage; against the magnitude-matched control it is +0.0007 [-0.0010, +0.0023] at k=3 and nothing at k=5. | **FACT:** preregistered triggers P4 (random destroys the gain) and P5 (magnitude-matched destroys the gain) both fired. **INFERENCE:** reaching the oracle needs a statistic near AUC 0.75; the two available observables give the same 0.55, so no conformity-score tuning closes a 20-point gap. The limit is the representation, not the estimator or the decision layer. | **A2S-RIP retracted as specified.** The machine verdict string `RIP_CEILING_ADMITTED` refers to R0a only; `decide()` was left unmodified rather than re-fitted after the result. |
| Gate R0c: certification layer | A harm-rate threshold fitted on `fit` transferred to `probe` essentially unchanged (alpha=0.40: 0.399 -> 0.387 at k=5; 0.396 -> 0.397 at k=3), but no coverage met alpha=0.20 or 0.30, and the gain at the certified coverage was +0.0001 [-0.0155, +0.0138]. | **FACT:** Fisch-style cross-task threshold transfer works on this substrate. **INFERENCE:** the certification layer is sound and reusable; there is simply nothing sharp enough to certify. | Retained as a carried-over asset for the successor branch. |
| Competing structural hypothesis, tested immediately | Truncating each target's own head to its top-s **coordinates** versus top-s **source principal directions**: at s=8, +0.0343 [+0.0152, +0.0523] versus +0.0085 [-0.0043, +0.0234]; full head +0.0542. Across 52 targets the top-weighted coordinate takes 20 distinct values of 26. | **INFERENCE - this reinterprets G2 rather than contradicting it.** A set of sparse vectors with *different supports* has an approximately flat covariance spectrum, so G2's flat spectrum is the signature of heterogeneous sparsity, not an absence of structure. Measured effective sparsity is s ~ 8 in a d = 26 generic descriptor basis, and `8*log(26/8) ~ 9-10` is exactly the G4 learning-curve knee at k ~ 10. The theory produces the measured number. | Opened branch `research/a2s-hotspot-sparse-20260802`. |
| New branch theoretical basis | The binding hot-spot principle (Clackson & Wells 1995; Bogan & Thorn 1998) and Free-Wilson additivity (1964) predict that in an interaction-determinant basis a target's head is sparse on a *target-specific* small support. This retrodicts G2 (flat spectrum), G3 (pooled ESM-2 cannot express which residues form a hot spot), G4 (the knee equals `s log(d/s)`), and R0b (a dense estimate is blunt per compound). | **HYPOTHESIS H0:** a representation exists in which the head is 2-3 sparse, moving the knee from k ~ 10 to k ~ 5. Falsifiable: effective `s` must fall below ~6 with top-s truncation retaining >= 60%, and the G4 curve must shift left. | `A2S_HOTSPOT_BRANCH_CHARTER_2026-08-02.md` with gates H1-H5 and a stop rule. H1 already passed in-basis. No `model/` promotion, no recipient label access, no commit or push occurred. |

## A2S PIRS And Conformational-State Successor (2026-08-02)

| Attempt | Observation | Root cause / interpretation | Decision |
| --- | --- | --- | --- |
| PIRS source-only R0 | Synthetic k=1/3/5 control passed. On 52 scaffold-disjoint probe targets in 50 components, the segment full-support oracle CI gain was +0.00460 [-0.00137,+0.01023], low-similarity oracle -0.00077 [-0.00997,+0.00767], and k=5 state gain +0.00025 [-0.00273,+0.00329]. Correct support did not beat deranged or wrong-target residuals, protein-zero/transplant, ligand-only, pooled-protein, or frozen-random coordinates. | **FACT:** the protein-conditioned representation itself has no admitted held-target correction object. The null is upstream of support inference and cannot be repaired by a learned operator. | `INTERACTION_STATE_REPRESENTATION_NOT_ADMITTED`. R1 and promotion prohibited. Artifacts and hashes are in `A2S_INTERACTION_STATE_GATE_DECISION_2026-08-02.md`. |
| Successor biological principle | Conformational selection implies affinity is a population-weighted free-energy difference over physical protein states. Label-blind coverage found UniProt mappings for 232/232 fit targets, AlphaFold structures for 230, experimental PDB structures for 217, and at least two PDB structures for 209. | **HYPOTHESIS:** external physical states can supply a state-specific response representation that coarse sequence segments could not, while k=3/5 updates only one/two population logits. | Opened `research/a2s-conformational-free-energy-state-20260802`. Gate C0 only; no affinity training or source probe reuse. |

## Overall A2S Meta-Adaptation Progress And CFES C0B Stop (2026-08-02)

### Objective and cumulative evidence

The active objective remains unchanged: learn a transferable support-conditioned
adaptation state from abundant source targets that improves ranking for a
strictly unseen target from k={1,3,5} passive measurements. Reliability,
retrieval, uncertainty, calibration, and analytic sparse estimators remain
controls or auxiliary modules, not the final contribution.

The sequential evidence now establishes five boundaries:

1. TRACE admits support-label information only in a local chemical-relation
   stratum; below roughly 0.35 support-query Tanimoto, residual transport and
   learned pair reliability are null.
2. MODE measures genuine scaffold-disjoint target-specific ranking headroom,
   but its raw 26-dimensional response head is high-dimensional and has a
   learning-curve knee near k=10. Low-rank modes and protein-sequence priors do
   not move that object into the k<=5 budget.
3. RIP measures a large hindsight ceiling for sparse ranking intervention, but
   the available margin and posterior-uncertainty observables have AUC only
   about 0.55 and cannot identify the useful edits.
4. Heterogeneous coordinate sparsity explains the flat head spectrum, but PIRS
   shows that coarse protein segments plus ligand features do not create a
   held-target interaction state: its full-support and k-shot oracles,
   support-assignment controls, and protein destructions all fail.
5. CFES tests a genuinely different physical-state hypothesis. Structural
   coverage is broad, but the preregistered external semantic gate below shows
   that the available ligand/pocket representation is not a transferable,
   load-bearing physical interaction object.

These findings do not prove that few-shot unseen-target adaptation is
impossible. They do rule out the tested transport, low-rank mode, uncertainty
selection, coarse interaction-state, and pocket-composition conformational
routes as final mechanisms on the present substrate.

### CFES C0A and C0B execution

- C0A passed label-blind coverage: 232/232 source-fit targets mapped to UniProt,
  231/232 had AlphaFold DB coverage in the final record, 217/232 had an
  experimental PDB structure, 209/232 had at least two PDB structures, and 201
  independent fit components had at least two physical structures.
- C0B was frozen before fitting in
  `reports/active/A2S_CFES_C0B_STRUCTURAL_SEMANTIC_PREREGISTRATION_2026-08-02.md`.
  It loaded only registered raw PLINDER structure/contact columns and official
  train/validation rows. It did not load any affinity column, the outcome-
  exposed processed PLINDER registry, PLINDER test, ChEMBL affinity, source
  probe/locked labels, or recipient labels.
- All 467 A2S source accessions were excluded. After sequence mapping and
  molecule validation, 38,022 train and 561 validation ligand rows remained.
  Purging every train row sharing an accession, exact ligand, Murcko scaffold,
  PDB/system provenance, or registered cluster with validation left 13,094
  train rows. The audit contains 561 rows in 217 independent validation
  clusters; all post-purge overlap counts are zero.
- The CUDA run used an RTX 4060 and completed in 62.30 seconds. The rank-4
  synthetic positive control passed: cross-minus-additive loss gain 1.98281
  with 95% CI [1.74503,2.21817], every fold positive, and ligand/pocket
  shuffling removed more than 100% of the injected effect.
- On real held clusters, the learned rank-16 cross residual beat its additive
  parent by only 0.001719 [0.000402,0.003234]. This small effect was not the
  registered physical semantic object: cross was worse than ligand-only by
  -0.091663 [-0.143368,-0.046753], did not beat the parameter-matched no-cross
  residual (-0.001774 [-0.007113,0.003969]), and did not beat frozen random
  cross features (-0.000170 [-0.001028,0.000736]). Absolute cluster-macro loss
  was 1.13132 for ligand-only, 1.22470 for additive, and 1.22298 for cross.
- One of five seeds (-0.000839) and one of five audit folds (-0.000757) were
  negative. Hydrogen-bond and halogen-bond coordinates were negative, so the
  required contact breadth failed. Pocket shuffle, ligand shuffle, and
  structure transplant removed the tiny increment, but residue randomization
  removed only 20.4%, below the frozen 70% threshold.
- Final machine verdict:
  `CFES_C0B_SEMANTICS_NOT_ADMITTED_STOP_CFES`. C1 affinity representation
  training, C2 support-operator training, and C3 locked-source confirmation are
  prohibited for CFES. No capacity, depth, epoch, or threshold rescue is
  allowed.

### Interpretation, verification, and current boundary

The statistically detectable cross-minus-additive difference is a generic
feature-capacity effect, not evidence of transferable physical pairing: a
protein-free model is substantially stronger and matched/random interaction
controls explain the increment. Therefore conformational population adaptation
has not been shown generalizable or valuable on the operational representation,
despite the general biological plausibility of conformational selection.

The focused CFES/PIRS suite passes 22 tests. Authoritative C0B artifacts are:

- `reports/active/a2s_cfes_semantic_gate_2026-08-02.json`, content SHA-256
  `cbf128a54fb4bd71a228d9045f7459abac325059a07727eefc0f5c34140ac604`;
- `reports/active/a2s_cfes_semantic_gate_records_2026-08-02.parquet`, SHA-256
  `421ae155705231a24fde8149a7cd6a0cf3e0727a9690bde3a4a9c56d12d851ce`;
- `reports/active/a2s_cfes_semantic_gate_weights_2026-08-02.pt`, SHA-256
  `258de719a75552bec29f9a7de0cae8a850ae6ed1f33fb71d73ac24d8822010ca`;
- `research/a2s_cfes_semantic_gate.py` and
  `tests/test_a2s_cfes_semantic_gate.py`.

There is **no major breakthrough** and nothing has been promoted to `model/` or
`script/`. The historical source `probe` outcome was consumed once by PIRS and
may not be reused for model selection. Source `locked` and all recipient labels
remain sealed. Per the registered stop rule and user instruction, the next
research branch must begin from a genuinely different general biological
principle, use source `fit` only until its mechanism is frozen, and preregister
its own representation and k<=5 identifiability gates before implementation.

## 2026-08-02 — Gate T0: is there a transferable adaptation object?

Branch `research/a2s-transfer-object-20260802`. Runner
`research/a2s_transfer_object_gate.py`, tests
`tests/test_a2s_transfer_object_gate.py` (14 passed; `tests/` suite 296 passed).
Source `fit` and `probe` roles only; `locked` and the recipient roster were not
requested. The gate **trains nothing**. Runtime 38 s on the RTX 4060.

**Decision: `TRANSFERABLE_AT_FULL_SUPPORT_BUT_NOT_IDENTIFIABLE_AT_K5`.**

After nine falsified mechanisms, T0 measured the two premises all nine assumed
and none tested.

- **T0A — the headroom is majority measurement context.** A chemistry-free
  document-mean oracle scores +0.0610 [+0.0386, +0.0824], beating the full
  26-dimensional per-target head (+0.0527 [+0.0303, +0.0744]). The own-head gain
  splits into +0.0808 across documents and **+0.0313 [+0.0056, +0.0601]** within
  a document. The offset-free chemical prize is therefore ~60% smaller than the
  figure the programme has quoted since Gate A0, and Q1's "Tanimoto >= 0.55"
  stratum is the same stratum as "same measurement context". On a random
  (non-scaffold-disjoint) split the document oracle reaches +0.1048 of a +0.1137
  headroom, confirming the leakage signature scales as expected.
- **T0B — transfer is real at full support.** Selecting one of 110 source-target
  heads with ~64 recipient labels gives **+0.0266 [+0.0054, +0.0484]** over the
  frozen base and **+0.0346 [+0.0152, +0.0551]** over a pooled head. First
  positive transfer result in the programme. But only 40.7% of heads beat the
  base, the median head scores -0.0257, and only 12.9 of 110 reach half the best.
- **T0C — no label-free shortlist.** Pooled ESM-2 cosine (+0.0049
  [-0.0051, +0.0158]) and library chemotype overlap (-0.0097 [-0.0214, +0.0011])
  are both indistinguishable from a random shortlist.
- **T0D — k<=5 cannot select.** Break-even is **k ~ 20**. At k=5 selection scores
  -0.0104 vs base and -0.0024 vs the pooled head.

**Information account (no free parameter).** Held-out `tau = 0.642`,
`sigma = 1.343`, SNR 0.229, so `bits(k) = (k-1) * 0.1486` = 0.59 bits at k=5,
against `log2(110 / 12.9) = 3.09` bits required. Equating predicts k ~ 21.8; the
measured break-even is k = 20. The synthetic positive control passes (k=5 mode
recovery 18.2% vs 0.9% chance, +0.100 CI); k=1 recovery is exactly chance,
confirming k=1 rank silence is structural.

**Self-correction recorded.** The registered T0D comparator was
`select_k5 - random_selection`, which passed (+0.0155 [+0.0076, +0.0243]). It is
a bad comparator because the library's median member is harmful. The gate was
re-bound to the frozen base and the pooled head — a change made after seeing the
result but **in the failing direction only**, which cannot manufacture a
positive. `tau` was also re-estimated held-out after the first run; the in-sample
value (0.910) overstated the available bits threefold.

**Consequences.** Every future gate must report the same-document contrast, and
the document-mean oracle joins the mandatory control set. The route forward is
not a better estimator or another protein representation but a reduction of the
bits the support labels must supply: the measured requirement is
`M / M_useful <= 1.5` against the current 8.5. The mechanism this licenses,
A2S-FBA, is specified pre-implementation in
`reports/active/A2S_META_ADAPTATION_MECHANISM_DESIGN_2026-08-02.md`, whose first
action is the cheap Gate F1 oracle-ceiling measurement.

Artifacts: `reports/active/a2s_transfer_object_gate_2026-08-02.json` SHA-256
`7f54d1c486c2544d4ce346fc1b9039877cb84826074fd6d53f423896a31c9349`;
`..._records_2026-08-02.parquet` SHA-256
`12c02655f92ce00e4aa94bef01885b6c28073f008adbb2d0a72788bae91974d8`; lock
`6bcf6edc7140cadd2d0a6c0227a120e15401513060bd5bfeccb58e5b30be72ae`.

No promotion to `model/` or `script/`. No major breakthrough.

## 2026-08-02 — Gate T0 revision 2: the revision-1 entry above is CORRECTED

An external review of the Gate T0 entry identified six defects. Every one was
reproduced against the artifacts and is now fixed. The re-run under the corrected
code returns **`NO_TRANSFERABLE_CHEMICAL_HEADROOM_OBJECT_IS_MEASUREMENT_CONTEXT`
with all four gates failing**, superseding the revision-1 verdict.

**Retracted claims.**

1. *"Discrete transfer is real."* The transfer arm was scored on all query pairs
   only. On same-document pairs — the control T0A itself established as mandatory
   — full-support selected-head transfer is **-0.0183 [-0.0435, +0.0054]**, and
   was negative in all six seeds (-0.021 to -0.005). Withdrawn.
2. *T0A admission.* The same-document own-head lower bound is **+0.0046**, below
   the registered 0.005 threshold; across six seeds it never cleared it. The
   offset-free chemical headroom is **unresolved**, not established.
3. *"The information account predicts the break-even with no free parameter."*
   `M_useful` was defined post hoc as heads reaching half a noisy oracle maximum
   (12.8/110, 3.10 bits). Defining it as heads beating base gives 44.8/110 and
   **1.30 bits**. The agreement used the same data on both sides and the k=20
   interval crosses zero. Retracted; retained as an order-of-magnitude heuristic
   with caveats now shipped inside the artifact.
4. *"Trains nothing."* False — closed-form ridge heads are fitted from labels for
   every source target and every recipient support set. No gradient model is
   trained. Corrected.
5. *Support size.* Full support averages **198.1** labels (median 93), not ~64.
6. *Artifact hashes.* The runner hashed the JSON, appended an artifact block
   containing that hash, then rewrote the file, so the recorded digest never
   matched. Revision 1 quoted `7f54d1c4…` for a file hashing to `f05f1fc0…`.

**Reproducibility defect with programme-wide scope.** `build_basis` used
`torch.svd_lowrank`, which draws a random projection and accepts no generator: two
calls in one process differed by 9.13 max-abs, and replaying T0 flipped its
verdict. Replaced with an exact sign-fixed symmetric eigendecomposition of the
fit-role Gram matrix, verified deterministic across processes by a test. **Every
gate built on this basis is affected — A0–A4, G1–G4, R0, HOTSPOT — and their
recorded intervals should be treated as seed-dependent.**

**Firewall correction.** The ledger reserves `probe`: consumed once by PIRS, not
reusable for model selection. T0 evaluates on `probe` and is therefore retained as
**exploratory evidence only**. The proposed Gate F1 was probe-dependent and is
**withdrawn**; successors must use nested `fit`-only development tasks.

**What survives, and it is deterministic.** The document-mean oracle arm does not
depend on the ligand basis at all: **+0.0610 [+0.0386, +0.0824]**, beating the full
per-target chemical head (+0.0519). The new `provenance_audit` measures the split
it was evaluated on: **52/52 targets share documents** between support and query,
**91.1 %** of query rows reuse a support-side document, **88.8 %** reuse a
support-side assay, and 21 targets have every query row from a support-seen
document. Murcko-scaffold disjointness is not document, assay or series
disjointness, and every prior within-target result in this programme inherits this.

**Design consequence.** `A2S_META_ADAPTATION_MECHANISM_DESIGN` revision 2 withdraws
A2S-FBA as a Stage 1 proposal: its premise is retracted, its entropy term rewarded
confident wrong routing, its harmlessness term certified single operators rather
than the deployed mixture, it cannot improve k=1 by construction, and a
method-specific analysis (Modular Meta-Learning, MMAML, CNAPs, VERSA, neural
processes, information-theoretic meta-learning, MetaDTA, AdaMBind) places it inside
existing modular meta-learning. The replacement first action is **Gate D0**:
rebuild the split with simultaneous target/scaffold/document/assay separation,
score on same-document pairs from documents absent from support, and replace the
heuristic bit account with prospective utility plus an empirical information
estimate.

Tests grew to 20 for this gate (302 across `tests/`), adding the four contracts the
first suite lacked: end-to-end determinism, provenance overlap, artifact-hash
integrity, and the same-document estimand.

No promotion to `model/` or `script/`. No major breakthrough. Gate F1 does not run.

## 2026-08-02 — NEA preconditions D0/N0/N1: the programme's terminal measurement

Branch `research/a2s-transfer-object-20260802`. Runner
`research/a2s_nea_preconditions.py` (`main.py a2s-nea-preconditions`), tests
`tests/test_a2s_nea_preconditions.py` (9 passed). Source **`fit` only** — `probe`
was deliberately excluded because PIRS consumed it, so the confirmation roles stay
untouched. `locked` and recipient labels never requested. Deterministic basis, no
gradient model.

**Decision: `NO_CHEMICAL_ADAPTATION_OBJECT_SURVIVES_SEPARATION_STOP_PROGRAMME`.**

**D0 — FAIL.** The same 93 fit targets / 92 components, evaluated under two splits
differing only in what they separate:

| | scaffold-only | separated (scaffold+document+assay) |
|---|---:|---:|
| document / assay overlap | 88.6% / 87.0% | **0.0% / 0.0%** |
| document-mean oracle (chemistry-free) | **+0.0671** [+0.0504,+0.0854] | **+0.0000** [0,0] |
| own head, all pairs | **+0.0610** [+0.0448,+0.0784] | **+0.0044** [-0.0161,+0.0242] |
| own head, same-document pairs | +0.0275 [+0.0125,+0.0435] | **+0.0123** [-0.0041,+0.0281] |

The harness validated itself structurally: on a document-disjoint split the
document oracle can only predict a constant, and it measured exactly zero. The
head's all-pair advantage collapses **93%**. This is not a harder task — base
concordance is essentially unchanged (0.5499→0.5424 all-pair, 0.5037→0.5157
same-document) and the separated split has 5x more within-document pairs
(9,290 vs 1,891). The base performs the same; the head's advantage evaporates.

**Honest bound, not an over-claim.** The same-document point estimate is +0.0123
and positive in both regimes. At SE 0.0082 with 92 components, resolving it to a
lower bound above 0.005 needs **~445 components, ~4.8x the present corpus**. The
terminal statement is therefore *"the effect, if real, is ~+0.012 target-macro CI
and is beyond this corpus's resolution"*, not *"the effect is zero"*.

**N0 — the nuisance quantified.** Across 1,421 documents / 1,403 assays with >=5
rows, per-context offsets explain **68.1% / 68.6% of all residual variance**, with
offset SD **~1.7 pKi** — larger than the target-specific chemical effects being
sought. Log-scale SD is 0.40, so the acting group is the **full affine group**, not
offset-only. This is the quantitative explanation of the entire failure sequence:
any estimator consuming absolute residuals spends its budget on a nuisance roughly
twice the size of the signal.

**N1 — PASS, and it exonerates the mechanism.** Passive k=5 draws contain a
document holding >=2 supports 85.9% of the time (>=3: 54.6%; mean 3.79
within-context pairs); k=1 coverage is exactly 0%, confirming structural k=1
silence. NEA had a deployment path. There was nothing for it to learn.

**Consequences.** The mechanism track stops per the registered stop rule in
`A2S_STAGE1_MECHANISM_SEARCH_V2`; NEA is not implemented and gates N2-N8 do not
run. Nothing is promoted to `model/` or `script/`. The terminal deliverable is the
measurement, with three components: the document oracle beating the chemical head
under conventional splits, the 93% collapse under provenance separation, and the
68% variance attribution to per-context offsets. Reopening requires ~445
provenance-separated components (Papyrus 05.7, BindingDB), not a new architecture.
