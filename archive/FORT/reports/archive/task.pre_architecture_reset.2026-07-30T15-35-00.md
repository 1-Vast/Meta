# FORT Current Research Task

**Last updated:** 2026-07-30

## Current Scientific Verdict

`UBSE-A1-v2` ended at `STOP_C1_A1R_STRUCTURE_COVERAGE_BELOW_FROZEN_MINIMUM`.
Phase A acquired and hash-froze 421 coordinate bodies for 459 instances, but
only 143/459 exact entity-sequence closures and 32/153 complete three-
deposition target units remained. The frozen minimum is 128 complete units.

No threshold, entity-sequence rule, role, or membership may be relaxed. Full
PLIP/ProLIF extraction, C2 coupling, C3 student, C4 affinity, confirmation,
and sealed access are stopped. No deployable `z_int`, structural teacher, or
residue-functional-group mechanism claim exists. P0A, A0W, pose ensembles,
flexible kernels, and A1-S are not rescue paths.

The historic conclusion also stands: current sources did not identify a
transferable strict dual-cold interaction mechanism. Historical work does not
authorize rereading or retraining.

## Active Program

```text
ACTIVE_PROGRAM = SIMA_DTA_FSA_D0
Working name   = SIMA-DTA
Expansion      = Support-Identifiable Mamba Adaptation for Drug-Target Affinity
```

Primary task: `k=5` few-shot adaptation for a previously unseen target, with
target/homology-cold meta-test targets and support-query scaffold-cold query
compounds. Secondary support sizes are `k in {0,1,3,5,10}`; `k=0` is a
reference baseline, not the new primary task.

Hypothesis: a small support set identifies target-specific ligand reordering
beyond target intercept, affinity scale, ligand potency, family identity, and
chemical-neighbour interpolation.

## Permissions And Data Boundary

```text
affinity_training_authorized = false
confirmation_authorized      = false
sealed_authorized            = false
fsa_d0_topology_power_audit  = true
```

FSA-D0 may use only safe metadata, frozen roles, cached non-outcome features,
and TRAIN-only cross-fitted B0 predictions. Do not read new numeric affinity
values, panel_davis target-conditioned confirmation, ChEMBL confirmation, or
sealed outcomes. Analyze pKd and pKi separately unless a conversion contract
is frozen. All CUDA-capable work uses `D:\anaconda\envs\drug\python.exe`;
CPU-only work is labeled. This task does not authorize model training.

## Three Paper-Core Modules

### M1: Long-Context Support-Conditioned Hybrid Adapter

Use frozen or cross-fitted B0 and cached frozen residue-token embeddings. A
lightweight adapter combines bidirectional Mamba propagation, sparse/local/
landmark Transformer retrieval, fixed-boundary support-memory recap, and
support-conditioned FiLM or low-rank adapters. At meta-test, do not update
protein or ligand encoders.

For `S_t = {(d_i,y_ti)}`, a permutation-invariant encoder maps ligand
representation, `B0(t,d_i)`, and `y_ti-B0(t,d_i)` to `c_t in R^q`, initially
`q in {4,8}`. One to three MAML or Meta-SGD steps update only `a_t`, `b_t`,
and `c_t`.

```text
y_hat_cal(t,d)  = a_t + b_t B0(t,d)
y_hat_full(t,d) = y_hat_cal(t,d) + phi_tilde(t,d,c_t)^T U c_t
```

Residualize `phi_tilde` against `[1,B0]` with a frozen ridge rule. The exact
interaction null is `c_t=0`; the primary comparator is intercept-plus-scale
calibration. Required controls are pure Transformer, pure Mamba, hybrid,
protein-free, random task-code, and matched B0 calibration. Retain Mamba only
if it exceeds frozen MDE or is non-inferior with a preregistered memory or
throughput margin. Report cached-adapter and full deployment cost separately.

### M2: Query-Span Support Design

For candidate support compound `d`, compute label-free
`j_td = partial y_hat(t,d) / partial c_t`. Select supports to minimize mean
query-span variance:

```text
mean_q j_tq^T (lambda I + sum_d_in_S j_td j_td^T)^-1 j_tq
```

Require scaffold diversity, support-query scaffold disjointness,
chemical-neighbour caps, frozen label-free ties, and selection frozen before
labels open. Compare uniform random, scaffold-diverse random, k-center,
uncertainty, and query-span selection. The claim is alignment to the
low-dimensional adaptation operator, not novelty of D-optimality.

### M3: Counterfactual Support Identifiability Objective

Each training episode has correct-target support, chemistry-matched wrong-
target support, label-permuted support, and calibration-only adaptation.

```text
L_meta = L_correct
       + lambda_wrong relu(m_wrong - (L_wrong - L_correct))
       + lambda_perm  relu(m_perm - (L_permuted - L_correct))
       + lambda_inc   relu(m_inc - (L_calibration - L_correct))
```

No credit unless correct beats wrong-target, labels beat permutations, full
adaptation beats calibration, and protein-conditioned beats protein-free.
Curriculum statistics are TRAIN-only, stop-gradient, uniformly warmed up, and
bounded by a nonzero uniform floor plus family/target caps.

## Ordered Execution

| Stage | Authorized action | PASS gate | STOP condition |
| --- | --- | --- | --- |
| FSA-D0 | Outcome-safe topology, role, depth, rank, and power audit | Strict k=5 is adequate | Any frozen adequacy floor fails |
| FSA-B0 | Freeze baseline registry and equal budgets | Reproducible role closure | Role, budget, or protected-outcome violation |
| FSA-M1 | Implement M1 after separate authorization | M1 exceeds calibration by frozen MDE | Calibration, protein-free, or random-code matches it |
| FSA-M2 | Add support design | M1+M2 exceeds M1 random support | No query-span gain or neighbour-only gain |
| FSA-M3 | Add identifiability objective | M1+M2+M3 exceeds M1+M2 | Wrong/permuted support matches correct |
| FSA-E0 | One-seed kill test | All controls pass directionally | Any mandatory falsifier fails |
| FSA-E1 | Multi-seed component-held-out development | Paired inference and bootstrap pass | Family, scaffold, document, or provenance dominates |
| FSA-C0 | Separately authorized confirmation | Independent prespecified confirmation passes | No implicit authorization exists |

FSA-D0 determines, per candidate meta-test target and endpoint: support counts
for `k in {1,3,5,10}`, query depth after removal, target/homology closure,
support-query scaffold closure, query-to-training scaffold and neighbour
closure, endpoint/assay/document/provenance, independent components, family
concentration, label-free rank proxy, and TRAIN-only MDE with paired-arm power.
Strict `k=5` must retain adequate targets, components, support/query scaffolds,
query depth, rank, and power before architecture implementation.

No failed gate may be rescued by a larger PLM, more Mamba/attention layers,
larger `q`, more epochs/seeds, extra losses, pose/coordinate inputs, P0A,
A1-derived representations, or support-label leakage.

## Baselines, Metrics, And Inference

Freeze equal budgets for B0 at k=0, intercept-only, intercept-plus-scale,
frozen encoder plus ridge, regularized full-head fine-tuning, standard MAML,
MetaDTA, AdaMBind-compatible, pure Transformer, pure Mamba, and SIMA-DTA.

Metrics: target-macro query RMSE/MAE, within-target Spearman, pairwise ligand-
reordering accuracy, concordance index, full-minus-calibration, correct-minus-
wrong support, query-span-minus-random support, and the k adaptation curve.
Inference is target/homology-component paired with component bootstrap. Rows,
ligands, seeds, and folds are not independent biological units. Freeze MDE and
material thresholds from TRAIN-only episodes before development scoring.

## Planning Cost

Planning estimates, not results: with cached residue tokens, `q<=8`, batch 8,
and 1,024 tokens, peak GPU memory is 4-7 GiB for Transformer-only, 3-6 GiB
for Mamba-only, and 5-8 GiB for hybrid episodes. On the RTX 4060 Laptop GPU,
a 1,000-episode smoke is estimated at 1-3 GPU hours after cache generation.
FSA-D0 is CPU-first: 1-2 engineer-days and less than one GPU hour only if a
Jacobian proxy is admitted. M1/M2/M3 are estimated at 4/2/2 engineer-days
after FSA-D0 PASS.

## Authoritative Artifacts

- `reports/archive/task.pre_sima_dta_transition.2026-07-30T15-12-00.md`
- `reports/archive/history.pre_sima_dta_transition.2026-07-30T15-12-00.md`
- `reports/archive/sima_dta_transition_pre_rewrite_sha256.2026-07-30T15-12-00.txt`
- `reports/active/sima_dta_architecture_blueprint_2026-07-30.md`
- `reports/active/sima_dta_fsa_d0_preregistration_draft_2026-07-30.md`
- `reports/active/sima_dta_novelty_matrix_2026-07-30.md`
- `manifests/sima_dta_fsa_manifest.v1.json`
- `reports/active/ubse_a1v2_falsification_round.v1.json`

## Closed Routes And Reopening Conditions

| Route family | Closed state | Valid reopening condition |
| --- | --- | --- |
| DCST, SISMT, DTIOD | Support or semantic gates failed | New source-resolved interaction object passes frozen gate |
| PCIC, PB-CEC, RDIB, PD-MVR, RBSDD | Provenance, topology, replication, or identification failed | New exact independent lineage and topology |
| OpenMut, MEDIP, Tau | Evidence or teacher admission failed | New exact intervention/source/reliability substrate |
| Pose, P0A, A0W, flexible kernel | No deployable information or wrong role | New admissible source plus preceding gates |
| UBSE G1 and A1-v2 | G1 residual failed; A1 C1 coverage failed | New route only; no A1-v2 relaxation |
| LEXOR routes | Observability or memory control failed | Fixed citable label-memory-free channel |

Detailed historical outcomes, firewalls, and reopening conditions are retained
in `history.md` and the immutable snapshots.
