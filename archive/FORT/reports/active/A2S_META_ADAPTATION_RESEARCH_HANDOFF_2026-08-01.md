# A2S-DTA Meta-Adaptation Research Handoff and External Research Prompt

Date: 2026-08-01  
Status: research handoff; no model is authorized for implementation  
Scope: abundant-to-scarce drug-target affinity prediction under strict unseen-target evaluation

## How to use this handoff in a chat window

This file is deliberately self-contained because an external research agent may receive only a small attachment packet.

Recommended attachment order:

1. **Required:** this handoff file.
2. **Required:** `research/a2s_cmal.py`, the actual executable CMAL implementation.
3. **Required:** `a2s_cmal_v7_disjoint_source_meta300_seed1729.json`, the latest complete source-only failure result.
4. **If a fourth attachment is possible:** `CMAL_FAILURE_HANDOFF.md`.
5. **If more attachments are possible:** `verify_claims.py`, `CMAL_AUDIT_AND_REPAIR.md`, `a2s_bottleneck_diagnosis_and_minimal_mechanism_2026-08-01.md`, and `a2s_camp_research_proposal_2026-08-01.md`.

If only three files can be attached, use items 1–3. This document contains the essential history, corrected claims, candidate mechanisms, and open questions that would otherwise be lost.

## 1. Non-negotiable scientific objective

The target contribution is a **learned transferable adaptation mechanism**, not a better generic DTA encoder with a few-shot wrapper.

The desired operator must learn across abundant source targets how to use only `k={1,3,5}` measured support affinities from a strictly unseen recipient target to produce target-specific and query-dependent improvements in compound ranking over an identical frozen support-free DTA base.

The final mechanism must not reduce to any of the following:

- intercept or slope calibration;
- label interpolation;
- similarity retrieval;
- ridge or kernel regression;
- a fixed Gaussian-process or Bayesian posterior update;
- simple fine-tuning;
- a larger ligand/protein encoder;
- an ordinary CNP/ANP/set encoder renamed as an adaptation operator;
- a learned functional update that does not establish novelty beyond MetaFun or learned optimizers.

Analytic and conventional methods are mandatory baselines, not forbidden experiments.

The primary scientific endpoint is target-specific compound ranking. CI, Spearman, and NDCG@10 must be reported alongside RMSE/MAE so that calibration gains cannot be misrepresented as ranking adaptation.

## 2. Epistemic rules

Every conclusion must be labeled as one of:

- **FACT:** directly established by code, stored experiment output, or a primary source.
- **INFERENCE:** the most defensible interpretation of available evidence.
- **HYPOTHESIS:** a claim that requires a new experiment.

Do not convert an empirical failure to detect information into an information-theoretic impossibility theorem. Do not call a split untouched after it has been repeatedly inspected. Do not call an operator protein-conditioned merely because a protein embedding is present in its input.

No recipient labels have been authorized for model selection. All immediate diagnosis and admission tests must be source-only.

## 3. Current task and data construction

The current benchmark creates target-level source episodes with a support set `S` and later/publication-ordered query set `Q`. Recipient targets are strictly unseen during source training. Only `k={1,3,5}` early measurements are exposed as support.

Important qualifications:

- **FACT:** recipient early support pools actually contain 6–31 measurements, with median approximately 11; only 1/3/5 are revealed. This is a label-efficiency simulation, not evidence that only k measurements historically existed.
- **FACT:** source support compounds are sampled from an early pool and queries from later documents/years. This is not a verified design-measure-make-test-analyse sequence.
- **FACT:** no project identifier, medicinal chemist/team identifier, compound registration timestamp, test-result timestamp, decision-cycle identifier, or contemporaneous candidate pool is available.
- **FACT:** recipient query/support scaffold overlap is approximately 8.7%; roughly 81% of queries are scaffold-cold; sampled nearest support-query ECFP Tanimoto is approximately 0.223.
- **FACT:** an additional read-only audit found exact support/query `assay_context_uid` overlap for only about 4 of 63 recipient targets across the audited k=5 draws. Source ordered k=5 episodes showed exact-context overlap in about 8.8% of episodes.
- **INFERENCE:** the current task may demand global cross-assay, cross-series extrapolation from chemically distant points rather than a biologically coherent local SAR update.

ChEMBL documents and assays are curated database entities, not medicinal-chemistry project lineage. Publication order must not be described as campaign causality without additional evidence.

## 4. Actual CMAL mechanism

The executable model is `research/a2s_cmal.py`, not the legacy `model/` directory.

For support compound `i`:

\[
\mu_i=f_0(p,x_i),\qquad e_i=y_i-\mu_i.
\]

The adapter encodes support ligand/protein pairs plus residuals, applies support self-attention, attends from each query to support tokens, and predicts

\[
\Delta_q=s_\phi(q,p,S,e)\sum_i a_\phi(q,i,p,S,e)e_i,
\]

\[
\hat y_q=\mu_q+\Delta_q.
\]

Implementation facts:

- **FACT:** the final code computes `delta = delta_scale * measured_residual`.
- **FACT:** gradients eventually reach every adapter module and adapter parameters change materially.
- **FACT:** the first logged adapter step contains a NaN `delta_head` gradient, but later steps are finite. Therefore a permanently disconnected graph is not the primary explanation.
- **FACT:** at k=1 the attention-weighted residual anchor is query-invariant because the one-key attention weight is 1.
- **FACT:** the full k=1 correction can still be query-dependent because `delta_scale` reads the query and residual-conditioned context. Therefore the claim that k=1 query-specific output is impossible in principle is false; only the attention-weight-selection mechanism degenerates.
- **INFERENCE:** the residual-conditioned scale can learn a generic residual-triggered query perturbation without learning the correct support-compound-to-query relationship.

## 5. Latest complete failure result

File: `a2s_cmal_v7_disjoint_source_meta300_seed1729.json`

Run status:

- **FACT:** `formal=false`.
- **FACT:** one seed, `1729`.
- **FACT:** 1,420,290 total parameters; 1,106,049 base-optimizer parameters and 314,241 adapter-optimizer parameters.
- **FACT:** both predefined source gates failed.

### Source validation

| Metric | Frozen base | Adapted | Adapted minus base |
|---|---:|---:|---:|
| RMSE | 1.3743 | 1.5159 | +0.1417, worse |
| CI | 0.5691 | 0.5590 | -0.0101 |
| Spearman | 0.1951 | 0.1628 | -0.0324 |
| NDCG@10 | 0.6898 | 0.6743 | -0.0155 |

### Repeatedly inspected source `meta_test`

| Metric | Frozen base | Adapted | Adapted minus base |
|---|---:|---:|---:|
| RMSE | 1.6658 | 1.7231 | +0.0574, worse |
| CI | 0.5112 | 0.5142 | +0.0030 |
| Spearman | 0.0334 | 0.0423 | +0.0090 |
| NDCG@10 | 0.5976 | 0.5963 | -0.0012 |

Interpretation:

- **FACT:** validation worsens on every main metric.
- **FACT:** the second source split shows tiny CI/Spearman point gains but worsens RMSE and NDCG@10.
- **FACT:** correct-support point estimates often exceed wrong-support arms.
- **INFERENCE:** current specificity is compatible with correct support being less destructive than wrong support; it does not establish beneficial adaptation.
- **FACT:** source validation and the former source `meta_test` have both been repeatedly used during v3–v7 development and can no longer provide untouched confirmation.

## 6. Verified shortcut and control findings

- **FACT:** a four-arm chemistry-only classifier that receives no labels or protein identifies the correct arm at approximately 51.6% on meta-train and 54.0% on meta-validation versus 25% chance.
- **FACT:** the counterfactual loss was repaired so that wrong supports cease receiving reward once they are worse than the frozen base. The old accusation that current InfoNCE still rewards unbounded wrong-arm destruction is no longer accurate.
- **FACT:** in the supplied verification script, the k=1 `label_swap` residual equals the random-arm residual to numerical precision; the arms differ through support chemistry.
- **INFERENCE:** current negative arms can still reveal arm identity through chemistry, provenance, assay context, or residual scale.
- **FACT:** for k>=3, a chemistry-fixed permutation of residual values among the same support compounds is possible and preserves chemistry plus the residual multiset.
- **FACT:** for k=1, no within-support label-assignment permutation exists. A norm-matched sign flip or donor transplant can test residual sensitivity but cannot prove compound-label assignment specificity.

## 7. Corrections to earlier reports

The following claims must not be repeated without qualification:

1. **Invalid overfitting ratio:** comparing one sampled training minibatch MSE with aggregate holdout MSE does not establish an `8.5x` train/holdout overfit factor.
2. **Invalid universal lower bound:** `(sigma/tau)^2 approximately 29` is a useful diagnostic for a particular low-dimensional hierarchical estimator, not a universal sample-complexity lower bound for every learned operator.
3. **Invalid MDK sufficiency claim:** effective degrees of freedom close to k can also indicate weak shrinkage or noise fitting; it does not prove that nearly all useful information has been extracted.
4. **Overstated k=1 impossibility:** the residual anchor is query-invariant, but the total current CMAL delta is not.
5. **Overstated KRR/GP convex-hull claim:** KRR and GP posterior predictions are not universally constrained to the support-label convex hull.
6. **Overstated temporal causality:** forward-versus-reverse asymmetry can be produced by assay drift, publication selection, series replacement, and chemistry progression; it is not proof that later compounds were designed because of the exposed support labels.
7. **Prospective metadata defect:** future query document-year gaps are unavailable when ranking genuinely novel compounds and must not be used as inference-time covariates.

## 8. Committee diagnosis

The ranked diagnosis is:

1. **INFERENCE:** usable incremental support-label information has not been established under the current episode distribution.
2. **INFERENCE:** frozen-base generalization and residual geometry are too weak across unseen source targets; the fixed 1:1 base/adapter target split sacrifices abundant labels.
3. **INFERENCE:** the counterfactual task is vulnerable to arm recognition and residual/provenance shortcuts.
4. **FACT:** evaluation is underpowered for small effects and no longer confirmatory because both source splits were repeatedly inspected.
5. **INFERENCE:** CMAL's unconstrained residual-conditioned point perturbation lacks structural abstention and overcorrects.
6. **FACT:** numerical anomalies exist but do not explain the scientific failure.

Recommendation:

- Freeze CMAL as a failed, mechanically active baseline.
- Retain the data firewall, OOF/frozen-base interface, homology split, GPU episode store, and paired component-level evaluation.
- Replace the fixed 1:1 base/adapter split with leakage-safe component-level cross-fitting of the entire supervised encoder and head.
- Do not implement another high-capacity adapter before the information-admission tests below.

## 9. The actual unknown quantity

Let `Z_Q` denote the frozen base's query-pair ordering errors and let `E_S` denote support residuals correctly assigned to support compounds. The research question is whether

\[
I(Z_Q;E_S\mid p,X_S,X_Q,\mu_S,\mu_Q,M)>0,
\]

where `M` contains only metadata available at inference time.

This conditional mutual-information expression is a conceptual target, not something established by an arbitrary neural probe. A finite-model comparison can only show exploitable incremental predictive information within the tested model classes and statistical power.

At k=1, a support label can plausibly identify an assay/target offset, residual sign, or select a source-learned response motif. At k=3 or k=5, a few within-support contrasts become available. None of these budgets nonparametrically identifies a full unseen-target SAR field.

The most defensible candidate transferable object is therefore not `f_target(x)`. It is a shared rule of the form:

> Given a particular pattern of measured residual evidence, which frozen-base query preferences are justified to change, and when must the operator abstain?

## 10. Three candidate research directions and their current status

### 10.1 Counterfactual Evidence-Gated Discrete Rank-Edit Policy

Status: **HYPOTHESIS; lowest-cost current-data candidate; not paper-ready.**

The originally proposed continuous CSRIO was rejected. For n queries, the adjacent-difference matrix has rank `n-1` on the centered score subspace. Thus

\[
D\Delta=g\odot u,\qquad \mathbf 1^\top\Delta=0
\]

can express an arbitrary centered correction when `u` is unrestricted. Its sparsity penalty is scale-degenerate because `g -> epsilon g` and `u -> u/epsilon` preserve the output while driving `sum(g)` toward zero. A continuous gap edit is also not equivalent to an adjacent swap, and a no-op loss is not a structural guarantee.

A scientifically admissible repair is a truly discrete policy:

\[
a_l\in\{\mathrm{STOP},\mathrm{SWAP}_1,\ldots,\mathrm{SWAP}_{n-1}\},
\qquad l\le B,
\]

with:

- a small, fixed hard edit budget `B`;
- structural null-evidence masking that forces `STOP` when all residual evidence is zero;
- matched wrong/permuted evidence trained to prefer `STOP`, not to become worse than the base;
- query permutation equivariance;
- explicit distractor-insertion, query-subset, and library-size stability tests.

This is a learned support-evidence-to-ranking-action policy, not KRR, posterior inference, or recipient fine-tuning. However, it changes the claim to **fixed-library transductive meta-reranking**. It cannot simultaneously claim independent compound affinity prediction unless candidate-set invariance is established.

### 10.2 Evidence-Activated, Assay-Disentangled SAR Grammar

Status: **HYPOTHESIS; biologically stronger; current data probably inadmissible.**

\[
c=C_\theta\left(p,\{e_i\},\{\tau(x_i,x_j),e_i-e_j\},M_S\right),
\]

\[
\Delta y_q=\operatorname{Apply}_\theta\left(c,\{\tau(x_i,x_q)\},M_q,\mu_q\right).
\]

`tau` must represent explicit transformations: attachment site, matched fragment edit, stereochemistry, physicochemical delta, scaffold relation, and assay compatibility. A valid grammar requires a finite rule vocabulary, activation semantics, bounded composition, auditable execution traces, and an unsupported state. Two generic neural networks called `Controller` and `Apply` are not a grammar.

This direction must beat classical MMP/MMS, mixed-effects, KRR, graph-difference, and learned-kernel baselines. Current low support-query locality and assay discontinuity make it inadmissible unless a predeclared same-assay/MMP-connected source stratum first shows assignment-specific label information.

### 10.3 Selection-Outcome Campaign Adaptation

Status: **HYPOTHESIS; highest long-term biological value; requires a different dataset.**

For verified project histories:

\[
h_l=U_\theta(h_{l-1},x_l,y_l-f_0,\text{assay}_l),
\]

\[
P_\phi(x_q\text{ selected}\mid h_l),\qquad
\Delta y_q=O_\theta(h_l,\tau(H_l,x_q),p).
\]

The selection head models why a compound entered a campaign; the outcome head models potency. This requires project IDs, contemporaneous candidate pools including unselected compounds, registration/test timestamps, series lineage, assay/construct metadata, and preferably quasi-experimental or intervention evidence. Current ChEMBL publication-order episodes cannot identify this model.

## 11. Mandatory pre-implementation gates

### 11.1 Leakage-safe split and base

- Downgrade all previously inspected validation and `meta_test` targets to development only.
- Create a new model-specific locked source split by protein-homology and document/assay provenance components.
- Generate every base residual using full target/component cross-fitting:

\[
\mu_i=f_0^{(-c)}(p_i,x_i),\qquad e_i=y_i-\mu_i.
\]

The entire supervised encoder and head must exclude component `c`; cross-fitting only the final head is insufficient.

### 11.2 Three necessary empirical quantities

Use the same nested architecture and mask only the label channel:

\[
G_0(p,X_S,X_Q,\mu,M),
\]

\[
G_1(p,X_S,X_Q,\mu,E_S,M).
\]

Define:

\[
\Delta_{label}=L(G_0)-L(G_1^{correct}),
\]

\[
\Delta_{assign}=L(G_1^{deranged})-L(G_1^{correct}),
\]

\[
\Delta_{headroom}=L(frozen\ base)-L(high\text{-}data\ source\ oracle).
\]

All three are necessary:

- `Delta_label > 0`: labels provide predictive information beyond chemistry and metadata.
- `Delta_assign > 0`: the correct compound-to-residual assignment matters, not only the residual multiset, assay offset, or added capacity.
- `Delta_headroom > 0`: frozen-base query ranking has meaningful correctable error.

### 11.3 Controls

- k=3/5 chemistry-fixed residual permutation among the same support compounds;
- k=1 norm/assay-matched sign flip or residual transplant, explicitly treated as sensitivity rather than assignment identification;
- residual-null;
- chemistry/norm/assay-matched wrong target;
- residual-noise dose response;
- target shuffle and protein shuffle;
- same-assay/MMP, same-scaffold, and scaffold-cold strata;
- temporal reversal only as an exploratory campaign diagnostic;
- synthetic injected-signal positive controls to show that the probes have adequate power.

Use at least a regularized pairwise linear probe, a capacity-matched DeepSets/CNP probe, and a small MetaFun-style probe. If synthetic positive controls are not detected, a null result is inconclusive.

### 11.4 Admission and stop rules

Proceed only if k=3 and k=5 show, on both development and the new locked split:

- positive component-bootstrap lower bounds for `Delta_label` and `Delta_assign`;
- effects above a practical threshold fixed before locked evaluation through power/MDE analysis;
- replication in at least two probe classes;
- deterioration under increasing residual-label noise;
- non-zero high-data source ranking headroom.

k=1 is exploratory because assignment specificity is not identifiable by within-support permutation.

If information is detected only in a same-assay/MMP stratum, narrow the claim to series-conditioned adaptation. If no information is detected globally or in predeclared biologically valid strata, stop complex meta-adaptation on the current episode construction.

## 12. Required baselines, ablations, and formal success criteria

Mandatory baselines:

- ligand-only DTA;
- frozen support-free DTA base;
- support intercept/slope calibration;
- ridge and KRR;
- closed-form Bayesian/MDK posterior;
- MAML and ANIL;
- CNP/ANP and MetaDTA;
- MetaFun;
- current CMAL;
- equal-capacity generic listwise reranker;
- direct inductive point-score adaptation.

For SAR grammar work, add classical MMP/MMS, hierarchical mixed-effects, graph-difference, and learned-kernel baselines.

Mandatory discrete-rank-policy ablations:

- remove hard `STOP`;
- remove the residual-evidence channel;
- remove counterfactual/no-op training;
- hard edit budget `B=1` versus `B=2`;
- adjacent action versus unrestricted reranking;
- full protein versus no-protein versus adaptation-branch protein shuffle;
- remove chemical relation features;
- correct evidence versus chemistry-fixed residual assignment permutation.

Primary statistics:

- target-macro CI;
- target-macro NDCG@10;
- pairwise proper log loss.

Secondary statistics:

- target-macro Spearman;
- RMSE and MAE;
- action coverage and edits per episode;
- conditional harm rate;
- risk-coverage curve.

The independent resampling unit is the target/homology/provenance component. Average draws and seeds within each component before paired component bootstrap.

A final model claim requires:

- absolute improvement over the identical frozen base, not only correct-versus-wrong separation;
- a locked-split paired lower confidence bound exceeding a preregistered minimum detectable/practically meaningful effect;
- superiority to the strongest analytic and equal-capacity neural baselines;
- loss of the main gain when the claimed innovation is ablated;
- wrong, deranged, and null evidence producing an equivalent no-op rather than deliberate degradation;
- five seeds before the sealed recipient evaluation;
- protein removal/shuffle destroying the adaptation gain before using the phrase `protein-conditioned`;
- query permutation and candidate-set stability before using the phrase `general DTA predictor`.

The final sealed recipient evaluation may be opened once, only after the algorithm, baselines, controls, seeds, metrics, and stop rules are frozen.

## 13. Primary-literature starting points

Use primary sources and check for newer work before making novelty claims:

- Conditional Neural Processes: https://proceedings.mlr.press/v80/garnelo18a.html
- MetaFun: https://proceedings.mlr.press/v119/xu20i.html
- MetaDTA: https://openreview.net/forum?id=yzlif16IASM
- AdaMBind: https://www.nature.com/articles/s41467-026-70554-5
- FS-Mol: https://openreview.net/forum?id=701FtuyLlAd
- PiRank: https://proceedings.neurips.cc/paper/2021/hash/b5200c6107fc3d41d19a2b66835c3974-Abstract.html
- SAM-DTA ligand-only precedent: https://academic.oup.com/bib/article/24/1/bbac533/6955272
- ChEMBL curation and data provenance: https://doi.org/10.1093/nar/gky1075
- SIMPD and real project chronology: https://doi.org/10.1186/s13321-023-00787-9
- Matched molecular series / Matsy: https://doi.org/10.1021/jm500022q
- Activity cliffs: https://doi.org/10.1021/ci3001138
- Public affinity measurement uncertainty: https://doi.org/10.1021/jm300131x

## 14. Unresolved research questions for the external agent

1. Is correct support-label assignment empirically identifiable under the present data distribution, or only within local assay/series strata?
2. Is there a theoretically justified transferable object smaller than a full target function but richer than calibration, similarity, or a fixed linear response?
3. Can a discrete evidence-gated ranking policy be formulated so that abstention and support-label dependence are structural, while avoiding arbitrary candidate-set dependence?
4. Does a rank-edit action class have meaningful oracle headroom at k=3/5, and can its advantage over unrestricted listwise scoring be made load-bearing?
5. Can EASG be defined as a genuine finite grammar with auditable rules rather than a generic mixture-of-experts?
6. What public or realistically obtainable dataset can support same-project, same-assay, time-resolved campaign claims?
7. How can protein context be made biologically meaningful: pooled whole-sequence embedding, binding-pocket residues, construct/variant/cofactor state, or another representation?
8. What formal partial-identifiability, abstention, or regret statement would make the method more than an application-specific decoder?
9. Which recent papers already cover support-conditioned reranking, selective meta-learning, evidence-gated policies, or target-conditioned SAR transformation?
10. Under what exact evidence should the entire current A2S-DTA paradigm be abandoned?

---

# Copy-Paste Prompt for an External Research Agent

```text
You are an independent research committee reviewing a failed abundant-to-scarce drug-target affinity (A2S-DTA) meta-learning program. Your task is not to repair the existing model by default. Your task is to determine whether a scientifically identifiable, transferable, learned adaptation mechanism exists under the current data and k-shot constraints, and, if it may exist, to define the strongest defensible form of that mechanism.

You may receive only a partial attachment packet because this is a chat-window review. The intended priority is:

1. A2S_META_ADAPTATION_RESEARCH_HANDOFF_2026-08-01.md — the self-contained research record and corrected committee synthesis.
2. research/a2s_cmal.py — the actual executable CMAL implementation.
3. a2s_cmal_v7_disjoint_source_meta300_seed1729.json — the latest complete source-only failure result.
4. CMAL_FAILURE_HANDOFF.md — optional extended chronology and data contract.
5. Optional audits and verification scripts if the interface permits more attachments.

Begin by stating exactly which attachments you can inspect. Do not silently invent the content of missing files. Treat numerical facts copied into the handoff as reported facts, but independently challenge every interpretation and novelty claim.

SCIENTIFIC OBJECTIVE

The desired contribution is a learned transferable adaptation mechanism. It must learn across abundant source targets how to use only k={1,3,5} support affinity measurements from a strictly unseen recipient target to improve target-specific, query-dependent compound ranking over an identical frozen support-free DTA base.

Calibration, interpolation, similarity retrieval, ridge/KRR, fixed GP or Bayesian posterior updates, simple fine-tuning, larger encoders, and generic support-conditioned set prediction are baselines only. Do not present a method as novel merely because it is neural, non-closed-form, attention-based, recurrent, listwise, or called an operator.

No recipient labels may be used for diagnosis, architecture selection, hyperparameter selection, or stopping decisions. The immediate work must remain source-only.

EPISTEMIC STANDARD

Label every substantive statement as:

- FACT: directly established by an attachment, executable code, stored result, or primary source.
- INFERENCE: an evidence-based interpretation.
- HYPOTHESIS: a claim requiring a new test.

Do not turn a finite-model null result into an information-theoretic impossibility theorem. Do not call repeatedly inspected source splits untouched. Do not infer protein conditioning from the mere presence of a protein vector.

CURRENT VERIFIED RESULT

The current CMAL is mechanically active but scientifically unsuccessful. It predicts

delta_q = learned_query_scale(q,p,S,e) * sum_i attention(q,i,p,S,e) * support_residual_i.

At k=1, the attention-weighted residual is query-invariant, but the total delta can remain query-dependent through the residual-conditioned query scale.

The latest one-seed, non-formal source-only run reports:

- Source validation: RMSE 1.3743 -> 1.5159; CI 0.5691 -> 0.5590; Spearman 0.1951 -> 0.1628; NDCG@10 0.6898 -> 0.6743.
- Repeatedly inspected source meta_test: RMSE 1.6658 -> 1.7231; CI 0.5112 -> 0.5142; Spearman 0.0334 -> 0.0423; NDCG@10 0.5976 -> 0.5963.
- Both predefined source gates failed.
- Correct support often beats wrong-support arms by point estimate, but correct support does not produce consistent absolute gain over the frozen base.
- Gradients reach the adapter and parameters change, so a dead graph is not the main explanation.

Data diagnostics include approximately 0.223 nearest support-query ECFP Tanimoto, about 81% scaffold-cold queries, approximately 8.7% scaffold overlap, and very limited exact support/query assay-context overlap. Publication/document order is not verified medicinal-chemistry campaign lineage.

Chemistry-only arm recognition is above chance. k>=3 permits chemistry-fixed residual-assignment permutation among the same support compounds. k=1 has no true within-support assignment permutation; sign-flip or norm-matched donor controls test residual sensitivity only.

PRIOR COMMITTEE VERDICT

CMAL should be frozen as a failed baseline. The primary unknown is whether correctly assigned support residuals provide exploitable incremental information about frozen-base query-ranking errors after conditioning on protein, support/query chemistry, base predictions, and inference-available metadata.

The committee considered three directions:

1. A discrete evidence-gated rank-edit policy with actions STOP or a small number of adjacent SWAP operations.
2. An evidence-activated, assay-disentangled SAR transformation grammar.
3. A selection-outcome campaign model for genuine project-level chronological data.

The original continuous sparse-rank formulation was rejected because its adjacent-difference parameterization can express any centered score correction, its soft sparsity gate is scale-degenerate when update magnitudes are unbounded, a continuous gap edit is not an adjacent swap, and a no-op loss is not a structural no-op guarantee.

The repaired discrete policy is only a conditional research direction. It changes the claim to fixed-library transductive meta-reranking unless it passes candidate-set stability tests. EASG likely requires same-project, same-assay, matched-series data. The campaign model requires project IDs, decision cycles, contemporaneous candidate pools, and selection/outcome separation; current ChEMBL publication-order episodes cannot identify it.

YOUR REQUIRED ANALYSIS

1. Independently audit the current scientific failure. Rank the likely causes among task misdefinition, absent or weak support-label information, source/recipient distribution shift, weak frozen-base geometry, shortcut learning, objective mismatch, operator misspecification, numerical optimization, and evaluation contamination. Explain what evidence discriminates these causes.

2. Formally analyze identifiability at k=1,3,5. Define the smallest plausible transferable object. State what can and cannot be learned from each support budget without pretending that k labels identify an arbitrary target-specific function.

3. Audit the corrected information-admission test:

   Delta_label = L(G0) - L(G1_correct)
   Delta_assign = L(G1_deranged) - L(G1_correct)
   Delta_headroom = L(frozen_base) - L(high_data_source_oracle)

   G0 and G1 must share the same nested architecture, with only the label channel masked. Residuals must come from component-level out-of-fold predictions of the complete supervised encoder and head. Assess false-positive and false-negative risks, proper scoring rules, synthetic positive controls, power, multiple biological strata, and the limits of interpreting a probe loss gap as conditional information.

4. Deeply review the discrete STOP/SWAP policy. Determine whether it is genuinely meta-learning or merely support-conditioned learning-to-rank. Examine structural null-evidence abstention, hard action budgets, action supervision, permutation equivariance, distractor insertion, query subset stability, library-size shift, transductive dependence, score calibration, oracle action-class headroom, and novelty relative to CNP/ANP, MetaDTA, MetaFun, learned optimizers, PiRank, and selective prediction.

5. Deeply review EASG. Provide a precise finite rule vocabulary, activation semantics, composition rules, capacity bound, unsupported state, and execution trace if it can be made real. Compare it with classical matched molecular pairs/series, Matsy, hierarchical mixed-effects, graph-difference prediction, activity-cliff methods, relation networks, and learned kernels. Reject it if it remains two generic networks under a new name.

6. Deeply review the campaign direction. Specify the minimum project-level data needed to distinguish compound-selection policy from biochemical outcome. Address adaptive sampling, missing-not-at-random data, publication bias, time-varying confounding, unavailable counterfactual candidates, temporal reversal, and prospective metadata leakage. State whether any public dataset can support the required claim.

7. Search current primary literature, including work published through 2026, for the closest prior art in few-shot DTA, conditional neural processes, functional meta-learning, learned update rules, selective/abstaining meta-learning, support-conditioned ranking, molecular relation networks, matched-series transfer, and real medicinal-chemistry temporal validation. Provide direct links and distinguish exact precedent from analogy.

8. Propose three to five candidate mechanisms only if they survive the prior-art and identifiability audit. For every candidate provide:

   - one-sentence contribution;
   - mathematical definition;
   - exact learned object;
   - source meta-training procedure;
   - meta-test procedure;
   - why k=1/3/5 is or is not feasible;
   - why it is not calibration, kernel regression, posterior inference, retrieval, CNP/MetaFun, or ordinary listwise prediction;
   - structural leakage and shortcut controls;
   - maximum scientific risk;
   - a decisive source-only falsification experiment.

9. Act as a hostile Nature Machine Intelligence, NeurIPS, and Bioinformatics reviewer. Reject renamed components, application-only novelty, candidate-set leakage, unsupported causal language, selective cherry-picking, and claims based only on correct>wong separation. State which venue, if any, each surviving direction could plausibly target and what evidence is missing.

10. End with a decision tree, not a vague roadmap:

   - conditions under which the current A2S-DTA episode construction must be abandoned;
   - conditions under which the claim must be narrowed to same-assay or series-conditioned adaptation;
   - conditions under which a discrete rank policy is admitted;
   - conditions under which a different dataset is mandatory;
   - exact evidence required before any model implementation;
   - exact evidence required before opening sealed recipient labels.

MANDATORY BASELINES AND CONTROLS

Include ligand-only DTA, frozen base, calibration, ridge, KRR, a closed-form Bayesian/MDK posterior, MAML, ANIL, current CMAL, CNP/ANP/MetaDTA, MetaFun, an equal-capacity generic listwise reranker, and a direct inductive point-score adapter. Add classical MMP/MMS and graph-difference baselines for SAR grammar claims.

Controls must include chemistry-fixed label assignment permutation for k>=3, k=1 norm/sign sensitivity controls, residual-null, chemistry/norm/assay-matched wrong support, label-noise dose response, target/protein shuffle, temporal reversal, support/query relation strata, query permutation, distractor insertion, query subset, and library-size shift.

Primary metrics are target-macro CI, target-macro NDCG@10, and pairwise proper log loss. Secondary metrics are Spearman, RMSE, MAE, action coverage, edits per episode, conditional harm, and risk-coverage. The target/homology/provenance component is the statistical unit.

Do not write code or propose a broad hyperparameter sweep. Do not force a winner. If no candidate is currently paper-worthy, state that clearly. The most valuable result may be a rigorous admission test or a justified rejection of the current paradigm.

OUTPUT FORMAT

Part 1 — Independent factual audit and corrected failure diagnosis.
Part 2 — Identifiability analysis for k=1/3/5 and definition of the transferable object.
Part 3 — Primary-literature and novelty audit.
Part 4 — Three to five surviving mechanisms, or a clear statement that fewer survive.
Part 5 — Hostile reviewer verdict for each mechanism.
Part 6 — Minimal source-only experiment matrix with pass/fail criteria.
Part 7 — Final decision tree and recommended research direction.
```
