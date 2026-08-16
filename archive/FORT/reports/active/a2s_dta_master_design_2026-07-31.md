# A2S-DTA Master Design And Falsification Record

Date: 2026-07-31  
Endpoint: pKi primary; pKd independent secondary replication  
Current final judgment: **STOP_FOR_NEGATIVE_TRANSFER**

This document is the authoritative design record for the A2S-DTA objective.
It separates a deployable proposed model from the descriptive support-oracle
control and records the evidence required before any neural architecture is
trained.

## 1. Model name

**A2S-CFRA**: Abundant-to-Scarce Cross-Fitted Recipient-Adaptive routing with
Risk-calibrated Abstention.

## 2. One-sentence core innovation

A2S-CFRA estimates source-specific recipient transfer utility from support-only
evidence using source-fold cross-fitting, then combines only the sparse positive
utility sources and abstains when a calibrated lower confidence bound does not
beat recipient calibration.

## 3. Real problem solved

Public DTA data are target-long-tailed: abundant targets can support stable
source experts while scarce targets provide only `k in {1,3,5}` measured
support labels.  The deployment question is not generic unseen-target
prediction.  It is whether the support labels of a new recipient identify
which already-trained source experts are useful, how much to transfer, and when
to refuse transfer.  The primary task is target-side single-cold: query drugs
may have appeared on source targets; target IDs are strictly disjoint.

## 4. Root causes of current failure

1. Protein similarity is not a sufficient statistic for transfer utility.
2. A scarce support set cannot identify a full recipient model.
3. All-source pooling contaminates the recipient with unrelated target
   mechanisms.
4. The current scalar support-compatible selector uses support labels directly
   and is therefore a descriptive support-oracle control, not a learned
   transfer estimator.
5. Pseudo-recipient query depth and chemistry differ from natural scarce
   recipients, so an untruncated meta-objective can learn the wrong utility.
6. A gate based on an uncalibrated score such as `score > 0` cannot guarantee
   negative-transfer control.

## 5. Falsifiable hypotheses

* H1: Source-fold cross-fitted utility features predict positive recipient-level
  transfer better than protein-only, chemistry-only, random, all-source, and
  scalar support controls.
* H2: A sparse mixture of the top positive-utility source experts improves
  target-macro RMSE and within-target ranking over the strongest no-transfer
  control at at least two of k=1,3,5.
* H3: A support-only calibrated abstention gate reduces negative-transfer rate
  without rejecting most recipients.
* H4: Any retained meta-learning step must beat an equal-budget non-meta linear
  router; otherwise meta-learning is not a contribution.

The current corrected pKi router falsifies H1-H3 for the implemented linear
router: target-balanced router gains are `-1.144`, `-1.213`, and `-1.315` at
k=1,3,5, with only about 20% benefiting recipients.  The hypotheses remain a
design target, not a claim.

## 6. Literature and transfer matrix

| Work/family | A. Established in original work | B. Structural analogy | C. A2S-DTA unverified transfer | D. Failure condition |
| --- | --- | --- | --- | --- |
| MAML/FOMAML | Optimize an initialization for rapid task adaptation | Recipient is a task and support/query episodes are available | Initialization may improve scarce target calibration | Does not identify which donor target to use; can overfit pseudo tasks |
| ANIL | Adapt a small head while freezing representation | Rank-1/low-rank recipient update | Useful equal-budget adaptation control | Frozen representation may not encode cross-target relation |
| MetaDTA/CML | Episodic target-conditioned DTA adaptation | Target-as-task construction | Valid meta-learning baselines | Does not define source-specific donor utility or abstention |
| AdaMBind | Target-as-task support/query and adaptive task scheduling | Scheduler is a task-selection control | Equal-budget scheduler baseline and task curriculum | Scheduling training tasks is not recipient-time donor routing |
| Adaptive task schedulers | Weight tasks by training difficulty or gradient utility | Utility-weighted source episodes | Can supply a meta-objective | Query utility or gradient similarity can leak or optimize wrong task distribution |
| Task-relation graphs | Transfer through learned task relations | Source-recipient relation graph | Graph features can be a router input | Graph similarity may be a chemistry/protein shortcut |
| Sparse MoE/routing | Select a small subset of experts by context | Source experts are DTA target experts | Top-L sparse source mixture | Without support-only utility, it is all-source or nearest-neighbour routing |
| Negative-transfer learning | Downweight harmful domains/tasks | Abstain from harmful source transfer | Calibrated transfer-risk gate | Mera et al. 2025 already uses meta-learned source weighting; this is not novel by itself |
| kNN-DTA | Retrieve chemically similar compounds/targets | Chemistry compatibility control | Strong nearest-neighbour baseline | Apparent gain can be ligand shortcut, not target transfer |
| Mera, Vogt & Bajorath 2025 | Meta-learning controls negative transfer by source-instance weighting | Source utility weighting | Negative-transfer mechanism precedent | Does not establish target-disjoint DTA recipient routing or natural-tail closure |

Novelty is therefore not “first meta-learning in DTA” or “first negative
transfer gate.”  The defensible claim, if it passes, is the combination of
strict target-disjoint recipient-level utility estimation, support-only
source-target routing, and calibrated refusal under endpoint/provenance
closure.

## 7. Innovation module 1: cross-fitted recipient-conditioned utility routing

Each abundant target `h` produces a frozen low-capacity expert `A_h`.  For a
recipient support set, the router uses source reliability, source/recipient
relation, chemical compatibility, and support-only predictive evidence.  The
router is fitted on pseudo-recipients held out from the source fit.  Its
training label is a recipient-level query utility, but natural recipient query
labels are never read during fitting, gating, or early stopping.

This module changes the identifiable quantity from generic target similarity to
`u(r,h) = expected loss reduction from source h for recipient r`.

## 8. Innovation module 2: calibrated negative-transfer abstention

The router predicts both a sparse transfer mixture and a lower confidence
bound for its expected gain.  If the lower bound does not exceed zero against
recipient calibration, transfer weight is set to zero.  The threshold is
calibrated only on held-out pseudo-recipient folds and is frozen before natural
recipient evaluation.  This is a selective prediction rule, not an accuracy
threshold tuned on query labels.

No third innovation is claimed.  Meta-learning is an implementation of module
1 and is retained only if it beats an equal-budget non-meta router.

## 9. Mathematical definition

For source target `h`, fit `A_h(d)` on `D_h` and a target-balanced pooled
baseline `B0(d)` on `D_A`.  For recipient support `S_r`, define source-only
support evidence with leave-one-support-out calibration:

```text
b_{r,h,i}^{(-i)} = mean_{j != i} [ y_{rj} - A_h(d_{rj}) ]
e_{r,h,i} = y_{ri} - A_h(d_{ri}) - b_{r,h,i}^{(-i)}
E_{r,h} = [ RMSE(e), rank_agreement, bias, source_reliability,
             protein_relation, chemistry_overlap, provenance_compatibility ]
```

For pseudo-recipient `r`, the held-out query utility label is

```text
u_{r,h} = L_r(B0 + recipient_calibration) - L_r(A_h + b_r).
```

The linear or ridge utility router is

```text
s_theta(r,h) = theta^T standardize(E_{r,h})
q_{r,h} = softmax(s_theta(r,h) / temperature), h in top-L sources
T_r(d) = sum_h q_{r,h} A_h(d)
g_r = P( u_{r,top-L} > 0 | E_r ) lower-confidence decision
yhat_r(d) = yhat_r^0(d) + g_r [ T_r(d) - yhat_r^0(d) ].
```

`yhat_r^0` is recipient calibration on `B0`, not a query-fitted oracle.  The
meta-objective minimizes squared error of `u` on held-out pseudo-recipient
folds, with target-macro weighting and query truncation/resampling to the
natural scarce query-depth distribution.

## 10. Complete data flow

1. Read only audited public ChEMBL-37 TRAIN rows and preserve global registry
   row IDs.
2. Verify feature-cache row count and `conn_sha` against the full registry.
3. Separate pKi and pKd; never pool them.
4. Deduplicate target/parent/document/assay units and freeze source/recipient
   rosters.
5. Fit target-balanced `B0` and source experts on abundant targets only.
6. Build pseudo-recipient source folds and natural scarce recipient episodes.
7. Fit the utility router and gate on source-fold pseudo-recipient episodes.
8. At natural evaluation, expose only the recipient support labels to evidence
   construction; freeze router, gate, and all nuisance statistics.
9. Report target-macro predictive, ranking, transfer, benefit, abstention, and
   calibration metrics with recipient/component bootstrap intervals.

## 11. Cross-fitting mechanism

Partition abundant sources into three deterministic target folds.  For each
fold, fit `B0` and experts on the other two folds, generate pseudo-recipient
support/query episodes from the held-out fold, and compute query utility only
as a meta-training label.  Fit router parameters on two pseudo-folds and
calibrate the gate on the third, rotating the calibration fold.  Natural
recipient targets are excluded from every router and gate fit.

Source-row IDs and deterministic episode hashes are retained.  No query label
is used for source selection, gate threshold, hyperparameter choice, or early
stopping.  Source quality and source relevance are separate features: quality
is estimated from source-held-out residuals, while relevance is estimated from
recipient support and relation features.

## 12. Meta-learning training and inference

The minimum implementation is closed-form meta-learning: a ridge utility
regressor over cross-fitted pseudo-recipient episodes.  This is deliberately
not MAML or a neural router.  It answers whether meta-learning adds value
before spending GPU budget on gradients.

If and only if the linear router passes pKi natural-tail gates, a one-hidden-
layer utility head may be compared at equal parameter and episode budgets.  A
MAML/FOMAML/ANIL run is a baseline, not the proposed innovation.  At inference,
source experts and router parameters are frozen; only support evidence,
recipient calibration, top-L mixture weights, and the calibrated gate change.

## 13. Negative-transfer abstention

The gate is fit on held-out pseudo-recipient predicted-gain residuals.  It
reports risk-coverage rather than only a point estimate.  A recipient is
abstained when the lower 90% confidence bound of expected transfer gain is
`<= 0`; otherwise the sparse mixture is used.  Report coverage, gain at each
coverage level, negative-transfer rate among accepted recipients, and the
fraction rejected.  Random abstention and no-abstention are required controls.

## 14. Causal source of expected performance

```text
support-only evidence
 -> better source transferability estimate
 -> fewer irrelevant experts in the mixture
 -> more recipient-specific adaptation
 -> lower negative transfer
 -> lower RMSE and better within-target ranking.
```

Expected beneficiaries are scarce recipients whose support residual pattern is
consistent with one or more abundant source experts.  Recipients with no
support-compatible source should abstain and should not be expected to
improve.  RMSE and pairwise ranking should improve first; uncertainty
calibration may improve only if the gate is calibrated. NDCG may not improve
when affinity ties dominate.

## 15. Difference from existing methods

* Pooled, all-source, and target-balanced multi-task models have no
  recipient-specific donor decision.
* Fine-tuning, MAML, FOMAML, and ANIL adapt parameters but do not estimate a
  source-specific transfer utility or refusal rule.
* AdaMBind schedules training tasks; A2S-CFRA chooses donor experts at
  recipient inference under support-only evidence.
* MoE and ordinary routing can select experts, but without query-closed
  cross-fitting their score is not a recipient-level transfer estimator.
* Protein-only and chemistry-only routing are controls for relation shortcuts.
* The scalar support-compatible selector is a support-oracle upper control,
  not a deployable learned router.

## 16. Minimum viable implementation

The existing low-cost implementation is intentionally linear:

```powershell
D:\anaconda\envs\drug\python.exe main.py topology-audit --split train
D:\anaconda\envs\drug\python.exe main.py a2s-baseline --endpoint pKi --out reports/active/a2s_pki_targetbalanced_seed1729.json
D:\anaconda\envs\drug\python.exe main.py a2s-router --endpoint pKi --out reports/active/a2s_router_pki_targetbalanced_seed1729.json
```

It uses 1,034-dimensional frozen ligand features, 242 pKi source targets,
193 recipient candidates, target-balanced ridge fitting, source-fold router
training, and CUDA closed-form solves.  It does not train a new neural
architecture.  The current code's `predicted_utility > 0` gate is explicitly a
diagnostic, not the calibrated lower-bound gate specified above; because the
linear router already fails pKi, the calibrated implementation is not admitted
for post-hoc rescue.

## 17. Cheapest kill tests

| Test | Destroyed information | Expected result | If unchanged |
| --- | --- | --- | --- |
| Random source replacement | Source identity | Gain collapses | Router does not use source relevance |
| Protein relation shuffle | Protein relation | Protein contribution disappears | Protein feature is inert or shortcut |
| Chemistry summary shuffle | Source chemistry summary | Chemistry contribution disappears | Chemistry feature is inert or shortcut |
| Source reliability shuffle | Quality estimate | More negative transfer | Quality feature is not used |
| Support-label permutation | Recipient label/evidence alignment | Transfer gain collapses to calibration | Label leakage or non-support shortcut |
| Support-compound permutation | Support chemistry/evidence alignment | Chemistry-conditioned gain collapses | Router ignores support compounds |
| No-label router | All recipient labels removed | Must match B0 or abstain | Hidden label/query leakage |
| Ligand-only router | Remove target relation and support residuals | Lower than full router | Full router adds no target task signal |
| Source-mean router | Replace source-specific expert with mean | Lower than sparse experts | Selection is not causal |
| Assay/source-only router | Keep only provenance fields | Must not beat full evidence | Provenance shortcut |
| All-source mixture | Set uniform weights | Lower or equal to sparse mixture | Sparse routing is unnecessary |
| Oracle source | Use query labels only as an upper bound | Upper bound above deployable router | No headroom for source selection |
| No abstention | Force `g=1` | Higher coverage, more harm | Gate does not control transfer |
| Random abstention | Match coverage randomly | Worse risk-coverage than calibrated gate | Gate is not informative |
| High-abstention stress | Raise threshold to low coverage | Harm drops but coverage is recorded | Surface gain is invalid |
| Pseudo versus natural | Compare roles under same k | Natural result must be independently positive | Pseudo-only success is not A2S success |

Each test is target-macro and endpoint-separated.  Any unchanged result in the
last column is a mechanism stop, not a reason to add capacity.

## 18. Strong baselines

Required comparison matrix: recipient calibration-only; global ligand-only
B0; pooled PCM; per-recipient ridge and RF; kNN-DTA; pooled plus fine-tuning;
target-balanced multi-task; MAML; FOMAML; ANIL; equal-budget AdaMBind;
all-source averaging; random source; protein-only; chemistry-only;
support-compatible scalar routing; cross-fitted non-meta router; proposed
cross-fitted router; proposed router plus abstention.

Every run records parameters, gradient evaluations, episode count, peak CUDA
memory, wall time, support usage, cross-fitting status, and query-label access.

## 19. Ablation and destructive controls

The two innovation ablations are: (a) remove cross-fitted utility training and
use fixed relation/chemistry routing; (b) remove calibrated abstention and
force transfer.  Additional ablations remove each evidence block, use dense
instead of top-L weights, use row-weighted instead of target-balanced fitting,
and use untruncated instead of scarce-matched pseudo queries.  No ablation may
select its threshold on natural query labels.

## 20. Compute budget

The corrected target-balanced closed-form pKi baseline uses about 38 seconds
wall time and under 50 MiB peak Torch memory on the RTX 4060 Laptop GPU; pKd
uses about 13 seconds.  The cross-fitted pKi router uses about 77 seconds and
117,126 meta rows.  Gradient evaluations are zero.  A neural router is not
authorized until the linear route passes.

## 21. PASS / STOP criteria

PASS requires all of the following: pKi positive gain over the strongest
no-transfer control at at least two k values; paired recipient-bootstrap lower
bound above zero and above the frozen material/MDE floor; superiority to
random, all-source, protein-only, chemistry-only, and scalar controls; a
majority of recipients benefiting; a lower negative-transfer rate after
abstention without near-zero coverage; preserved gain in low chemical
similarity or scaffold-cold strata; at least one independent natural-tail
positive result; and stability across seeds and recipient/component bootstrap.
If meta-learning is retained, it must beat equal-budget non-meta routing.

STOP if pKi fails at k=3/5, routing matches random/all-source, support-label
permutation does not matter, the gate rejects nearly everyone, gains are only
pseudo-tail, simple ridge/RF/kNN matches the proposed model, or any source,
document, assay, provenance, or query leakage is found.

## 22. Innovation audit

The current evidence supports only a benchmark and a descriptive support-oracle
control.  The target-balanced cross-fitted router fails pKi with RMSE-gain
AULC `-1.221`, and its mean benefiting-recipient rate is about `0.201`.  Thus
neither the routing innovation nor abstention innovation is identified.  The
literature audit also prevents claiming novelty for meta-learning, schedulers,
MoE, or negative-transfer weighting in isolation.

## 23. Final judgment

**STOP_FOR_NEGATIVE_TRANSFER**

Unique next action: freeze the current router route and run only the
pre-registered natural-tail document/time/source closure audit; do not train a
neural architecture or add a new routing feature until that audit supplies a
powered recipient roster and a new, source-only linear control is specified.
