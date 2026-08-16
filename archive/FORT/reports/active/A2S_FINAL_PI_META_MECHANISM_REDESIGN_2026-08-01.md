# A2S-DTA Final PI Review: Learned Transferable Meta-Adaptation

Date: 2026-08-01  
Status: final no-code mechanism redesign  
Evidence status: source-only development; locked-source and recipient labels remain sealed

## Executive decision

**FACT:** The balanced v2 source gate repaired the former 97% OOF-fold
concentration, retained 222/110/107 fit/probe/locked homology components, and
still found no positive lower bound for real `Delta_label` at k=1/3/5 or
`Delta_assign` at k=3/5. Synthetic label-channel controls passed, while a
high-data target oracle showed positive query-ranking headroom.

**INFERENCE:** The dominant failure is not lack of neural capacity. The present
passive episode asks k<=5 arbitrary, chemically distant and frequently
cross-assay measurements to identify a continuous target-specific residual
field. The observation design and adaptation object are mismatched.

**HYPOTHESIS:** The highest-potential direction is an **Active Diagnostic
Response-Mode Operator (ADRO)**: abundant source tasks jointly teach which
compounds are diagnostically worth measuring and how their outcomes select a
small, auditable, query-dependent intervention with structural abstention.

This is not an established mechanism result. It changes the formulation from
passive few-shot adaptation to **active k-measurement adaptation under the same
label budget**. The current ChEMBL data can support only a retrospective
finite-pool test, not a claim about historical medicinal-chemistry decisions.

---

## Part 1 - Final scientific problem

Let target t have protein p_t, a frozen support-free DTA predictor f_0, a legal
unlabelled measurement pool C_t, k measured support items S_t, and an untouched
query set Q_t. Source and recipient target/homology roles are disjoint.

The final problem is not simply

\[
\min_f \sum_{(p,x,y)} \ell(f(p,x),y),
\]

but

\[
\theta^*=\arg\min_\theta
\sum_{T_i\in\mathcal T_{source}}
L_{Q_i}\left(f_0+A_\theta(C_i,S_i,p_i),Q_i\right),
\]

subject to

\[
|S_i|=k\in\{1,3,5\},\qquad
A_\theta(C,\varnothing,p)=0,
\]

and no recipient-specific SGD or analytic solve.

The scientific claim must require all of the following:

1. The same frozen f_0 retains support-free DTA prediction capability.
2. Correct recipient evidence produces an absolute improvement over f_0 in
   target-macro ranking and affinity error, not merely less harm than wrong
   support.
3. The improvement depends on correctly assigned labels, target evidence and
   query identity.
4. The operator learned across abundant source targets how to acquire and
   interpret evidence; it is not calibration, interpolation, retrieval,
   kernel/GP inference, a closed-form posterior, or fine-tuning.
5. Unsupported or contradictory evidence produces an exact no-op rather than
   an unconstrained correction.

The transferable object is therefore not an arbitrary function f_t(x). It is:

> a source-learned rule for acquiring a few diagnostic observations and
> converting them into one of a small number of bounded query interventions,
> including the decision to do nothing.

---

## Part 2 - Why the current route failed

### 2.0 Task, data and split facts

**FACT:** The current executable CMAL path is `main.py a2s-cmal ->
research/a2s_cmal.py`. The repository `model/` directory is a legacy adaptation
stack and is not the implementation that produced the CMAL failure results.

**FACT:** The frozen formal CMAL episode package contains 30,123 label-blind
episodes built from 206 declared source targets and 63 recipient targets across
55 homology components. Target/accession/document/parent/assay role overlaps
are zero, and support/query measurement identity is frozen. Nine recipient
targets are homology-warm; therefore the full cohort supports a strict unseen
target-ID claim, while only the 54-target subgroup supports a homology-cold
claim. Publication/document ordering is not a verified DMTA chronology.

**FACT:** The balanced v2 information gate is a separate source-only diagnostic
over the ChEMBL dual-cold TRAIN registry, not the formal 206/63 recipient split.
It starts from 559 source pKi targets, assigns homology components before row
quarantine, and leaves its locked-source role plus every formal recipient label
sealed. Its negative result is evidence about support-label information under
that diagnostic construction; it is not a recipient performance estimate.

### 2.1 Root-cause audit

| Failure axis | Evidence | Scientific consequence |
| --- | --- | --- |
| Data/task construction | nearest support-query Tanimoto about 0.223; about 81% scaffold-cold queries; scaffold overlap about 8.7%; exact assay overlap rare; publication order is not a verified campaign | Passive supports often do not instantiate a relation that can be transferred to queries |
| Split/evaluation | old validation and `meta_test` were repeatedly inspected; v1 provenance closure put 97% of OOF rows in one fold | Old results are development evidence; v1 null was non-confirmatory |
| Balanced information gate | v2 OOF folds are approximately 20% each; real label/assignment intervals still cross zero; synthetic controls are strongly positive | The global passive episode has not admitted usable assignment-specific information |
| Identifiability | after centering, rank(H_k Phi_S)<=k-1 | k=1 gives zero contrasts; k=3 at most two; k=5 at most four, before assay nuisance |
| Adaptation object | CMAL estimates an unrestricted continuous query residual field | The object has more degrees of freedom than support evidence can identify |
| Inductive bias | CMAL can recognize arms/chemistry and make residual-triggered query perturbations; correct support can remain worse than f_0 | Specificity can be shortcut recognition rather than beneficial adaptation |
| Objective | MSE plus broad pairwise ranking is not aligned with top-tail NDCG, but RMSE/CI/Spearman/NDCG all failed together | Loss mismatch contributes but is not the primary cause |
| Optimization | adapter gradients and material parameter changes are present | The main failure is scientific, not a disconnected graph or too few epochs |

### 2.2 What information can actually transfer at k<=5

For

\[
r_t(x)=b_t+\phi(x)^\top z_t+\epsilon,
\qquad
H_ke_S=H_k\Phi_Sz_t+H_k\epsilon,
\]

the support can identify only a low-dimensional projection in the row space of
H_k Phi_S. It cannot identify a global unseen-target SAR surface.

| Information source | Defensible role | Identifiability boundary |
| --- | --- | --- |
| Protein | prior routing among source-learned response hypotheses; base DTA input | Previous correct-protein probes failed; any adaptation claim needs protein shuffle/removal evidence |
| Ligand chemistry | defines diagnostic candidates, transformations and query applicability | Chemistry without labels is task geometry, not target adaptation; similarity alone is retrieval/kernel |
| Support labels | OOF residual sign, magnitude and correctly assigned contrasts | k=1 has no assignment test; k=3/5 can identify only a few contrasts |
| Residual patterns | evidence relative to the same frozen f_0 | Must be full encoder/head component-OOF; in-sample residuals create shortcuts |
| Task latent | at most a low-cardinality state supported by diagnostic measurements | A continuous high-dimensional latent is not credible at k<=5 |
| Response modes | finite source-learned interventions or molecular rule states | Identifiable only if selected supports separate the modes above assay noise |
| Uncertainty | controls whether an intervention is safe | Uncertainty or abstention alone is not adaptation |
| Assay/context | nuisance control and possibly mode routing | Future query metadata is unavailable; publication/assay IDs cannot be used as campaign causality |

### 2.3 Why previous adaptation definitions were insufficient

**Bayesian posterior / MDK / ridge / KRR.** These estimate offsets,
low-dimensional coefficients or similarity-weighted residuals. They are useful
and mandatory baselines, but their task update is fixed or analytically solved.
Historical MDK gains largely match ordinary ridge, so they do not establish a
learned transferable adaptation law.

**MAML / ANIL / MetaDTA / AdaMBind.** These are direct prior art for target-as-
task DTA. A gradient update can fit k labels, but it does not solve the
information problem: with noisy, chemically distant supports, it estimates an
underidentified parameter direction. AdaMBind changes task scheduling around a
MAML learner, not the information contained in a fixed support set.

**CMAL contrastive adaptation.** CMAL is genuinely learned and mechanically
active, but it learns the wrong object. Its continuous residual correction can
use residual scale or arm identity without learning a valid support-query
relation. Counterfactual losses repaired specific shortcuts, yet absolute
ranking remained unstable and the balanced v2 label/assignment gate failed.

**Final diagnosis:** the central problem is not merely that the model did not
learn. The passive task asks it to learn an object that the available support
usually does not identify.

---

## Part 3 - Five candidate meta-learning mechanisms

All candidates preserve

\[
\hat y_q=f_0(p_t,x_q)+\Delta_q,
\qquad \Delta_q=0\text{ under null/unsupported evidence}.
\]

Everything in this part is a **HYPOTHESIS** unless explicitly labelled as
established prior art or a reviewer inference. These are candidate mechanisms,
not positive experimental results.

### Candidate 1 - Active Diagnostic Response-Mode Operator (ADRO)

**HYPOTHESIS - learned object.**

**One-sentence contribution.** Meta-learn which compounds should consume the
k-label budget and map their outcomes to a finite, sparse, query-dependent
response intervention.

**Learned object and source supervision.** Source episodes supervise an
acquisition policy pi_theta, evidence-state transition U_theta, finite
intervention dictionary {R_m}_{m=1}^M, and NULL gate. Query outcomes supervise
whether the selected measurements distinguish response modes and whether the
chosen intervention improves the identical f_0.

\[
x_j=\pi_\theta(h_{j-1},p_t,C_t\setminus x_{<j}),
\]

\[
r_j=y_j-f_0(p_t,x_j),\qquad
h_j=U_\theta(h_{j-1},\phi(p_t,x_j),r_j),
\]

\[
z_k=g_\theta(h_k)\in\{\mathrm{NULL},1,\ldots,M\},
\qquad
\Delta_q=\sum_{m=1}^M z_{k,m}R_m(p_t,x_q).
\]

M is small and frozen before locked evaluation; z is hard/sparse at test time;
offset-only R_m is not an innovative mode.

**Why k can suffice.** The measurements are chosen to discriminate a finite
mode vocabulary, not to cover the entire SAR surface. k=1 can support only a
binary/NULL diagnostic and cannot establish assignment specificity. k=3/5 can
test concordance and distinguish a few modes.

**Why it is different.** There is no recipient SGD, learned initialization,
analytic posterior, support-similarity interpolation, or donor retrieval.
Unlike passive CNP/MetaFun/meta-LTR, the load-bearing hypothesis is joint
experimental design plus a hard finite intervention and structural NULL.
Generic models can express parts of it; novelty exists only if these constraints
are empirically load-bearing.

**Minimum validation.** On identical source candidate pools compare random,
diversity, uncertainty, D-optimal and active KRR/GP acquisition. Require active
oracle headroom first. Then require learned policy > all acquisition baselines,
positive selected-support Delta_label and k=3/5 Delta_assign, absolute NDCG/CI
and RMSE improvement over f_0, and degradation under wrong support, label
permutation, protein shuffle and candidate-pool controls. Kill the direction if
random policy is equivalent or active-selected supports still fail information
admission.

**Largest risk.** The task contract changes, and target-specific candidate
pools can leak task identity. ChEMBL permits only retrospective finite-pool
simulation, not a historical campaign-policy claim.

**INFERENCE - PI verdict:** highest scientific potential; not yet admitted.

### Candidate 2 - Evidence-Activated SAR Program Induction (EASPI)

**HYPOTHESIS - learned object.**

**One-sentence contribution.** Compile same-assay matched-pair measurements into
a finite molecular transformation program and execute only exactly supported
rules.

For a support transformation e=(x_e^-,x_e^+),

\[
d_e=(y_e^+-y_e^-)-[f_0(p,x_e^+)-f_0(p,x_e^-)],
\]

\[
z_e=G_\theta(p,\tau_e,c_e,d_e)
\in\{\mathrm{UP},\mathrm{DOWN},\mathrm{NULL},\mathrm{UNSUPPORTED}\},
\]

\[
\Delta_q=\sum_{\tau\in P_t}a_\tau(z_\tau,p)
\mathbf 1[\tau\text{ exactly applies to }x_q].
\]

Source tasks supervise transformation reliability, protein/assay routing and
held-out matched-pair order. The vocabulary, attachment context, frequency
floor and depth-1 executor are source-frozen.

**Why k can suffice.** Two measurements form one contrast; k=3/5 can support a
few rules and consistency checks. k=1 must abstain. The claim is local
series-conditioned adaptation, not global scaffold-cold DTA.

**Why it is different.** Exact discrete transformation semantics replace
similarity weighting and unrestricted set-to-function mapping. It is not novel
if it reduces to Matsy/MMP statistics, copy-support-sign, mixed effects,
categorical empirical Bayes, graph-difference or a learned kernel.

**Minimum validation.** First perform a label-free same-assay/MMP coverage and
component-power census. If feasible, require a local k=3/5 assignment gate,
protein/assay routing ablations, rule/context/sign shuffles, exact no-op on
unsupported queries, and superiority to all classical rule baselines. Kill on
insufficient coverage or baseline equivalence.

**Largest risk.** Strict closure may leave too few independent components, and
the available assay IDs may not represent comparable protocols.

**INFERENCE - PI verdict:** most biologically concrete direction on current data, but only
conditionally admissible.

### Candidate 3 - Campaign-State Transition Operator (CSTO)

**HYPOTHESIS - learned object.**

**One-sentence contribution.** Learn a finite medicinal-chemistry state
transition induced by each measurement, assay change and selection decision,
instead of treating support as an unordered set.

\[
h_j=U_\theta(h_{j-1},\phi(x_j),y_j-f_0(p,x_j),a_j,b_j,d_j),
\]

\[
\Delta_q=\sum_m g_m(h_k)R_m(p,x_q),
\]

where a_j, b_j and d_j encode assay context, batch/context and the actual
selection action. Source projects supervise next-round ranking, future outcomes
and state transition consistency.

**Why k can suffice.** Each observation updates a finite campaign state rather
than estimating a function. This is plausible only if state transitions repeat
across projects.

**Why it is different.** It is not ordinary recurrent CNP if time, assay and
decision semantics are load-bearing. It has no recipient parameter update,
retrieval or analytic solve. If removing chronology/decisions has no effect, it
collapses to a standard sequence/set learner.

**Minimum validation.** Project-disjoint prospective splits; contemporaneous
candidate pools; time reversal, selection shuffle and assay shuffle controls;
comparison with chronological CNP/RNN and no-decision/no-assay variants.

**Largest risk.** Current ChEMBL lacks verified project IDs, decisions,
candidate pools and outcome timestamps. RetroDMTA supports retrospective
prioritization, but not automatic causal separation of historical selection and
biochemical response.

**INFERENCE - PI verdict:** high biological value, mandatory new data, paused now.

### Candidate 4 - Selective Discrete Rank-Edit Automaton (SDREA)

**HYPOTHESIS - learned object.**

**One-sentence contribution.** Restrict adaptation to STOP or at most B=1/2
discrete edits of the frozen ranking.

\[
a_\theta(S,p,Q)\in\{\mathrm{STOP}\}\cup\mathcal E_B(Q),
\qquad
\hat\pi_Q=E_{a_\theta}(\pi_{f_0}).
\]

Source query labels supervise legal-action regret or a permutation-invariant
bounded-action oracle.

**Why k can suffice.** The operator chooses among a few actions rather than a
continuous field; weak evidence returns STOP.

**Why it is different.** It is operationally distinct from calibration,
posterior inference and fine-tuning, but not a new general method family: it is
a selective support-conditioned meta-learning-to-rank policy. Adjacent edits
also imply fixed-library transductive dependence.

**Minimum validation.** Before implementation, require passive k=3/5
information admission and B=1/2 oracle headroom. Test query permutation,
distractor insertion, subsets, library-size shift and equal-capacity generic
meta-LTR. Kill if constraints are not load-bearing.

**Largest risk.** Candidate-set shortcuts and weak novelty.

**INFERENCE - PI verdict:** rejected as the paper's core; retain only as an action-class
diagnostic.

### Candidate 5 - Uncertainty-Gated Response-Mode Composition (UGRMC)

**HYPOTHESIS - learned object.**

**One-sentence contribution.** Learn whether current evidence is reliable enough
to execute a sparse composition of finite response interventions.

\[
e=E_\theta(p,S),\qquad (c,w)=G_\theta(e),
\]

\[
\Delta_q=c(e)\sum_{m=1}^M w_m(e)R_m(p,x_q),
\qquad c\in\{0,1\}.
\]

Source query outcomes supervise mode-specific intervention regret and selective
risk. Null, contradictory or unsupported evidence must give c=0 exactly.

**Why k can suffice.** The task is evidence concordance and selection among a
few modes, not recovery of a task function. k=1 permits only a simple
mode/NULL decision; k=3/5 can check agreement.

**Why it is different.** Uncertainty controls a target-specific intervention,
not only an interval, temperature or prediction rejection. Nevertheless, a CNP
or mixture-of-experts can represent the mapping; mode sparsity, risk supervision
and hard NULL must be load-bearing or novelty disappears.

**Minimum validation.** Compare at fixed coverage with support-free uncertainty,
selective CNP, conformal/selective prediction and chemistry/protein-only gates.
Require correct labels > shuffled/wrong labels at the same coverage, monotonic
label-noise response, and absolute improvement over f_0. Kill if performance is
explained only by abstaining more.

**Largest risk.** It becomes selective prediction or an amortized categorical
posterior under a new name.

**INFERENCE - PI verdict:** useful safety module for ADRO, insufficient as an independent
core innovation.

### Candidate ranking

| Candidate | Scientific clarity | k-shot identifiability | Current data feasibility | Novelty potential | Decision |
| --- | --- | --- | --- | --- | --- |
| ADRO | high | high if diagnostics separate modes | retrospective only | highest | select |
| EASPI | high/local | high for matched contrasts | unknown until coverage census | medium-high | conditional |
| CSTO | high | plausible finite state | no | high with new data | pause |
| SDREA | medium | bounded action | global gate failed | low | reject as core |
| UGRMC | medium | finite modes | global gate failed | low-medium | merge into ADRO |

### Prior-art boundary

**FACT:**

The following are exact or near-exact precedents and must be included as strong
baselines rather than described as absent literature:

- [AdaMBind](https://doi.org/10.1038/s41467-026-70554-5) and
  [MetaDTA](https://openreview.net/forum?id=yzlif16IASM): exact few-shot,
  target-as-task DTA precedents.
- [FS-CAP](https://pmc.ncbi.nlm.nih.gov/articles/PMC11267577/),
  [CNP](https://proceedings.mlr.press/v80/garnelo18a.html), and
  [MetaFun](https://proceedings.mlr.press/v119/xu20i.html): context-conditioned
  prediction and learned functional-update precedents.
- [Meta-learning to rank for sparsely supervised queries](https://doi.org/10.1145/3698876):
  exact precedent for support-conditioned meta-ranking.
- [Matsy](https://doi.org/10.1021/jm500022q): direct matched-series transfer
  precedent for any finite SAR grammar claim.
- [RetroDMTA](https://pubs.rsc.org/en/content/articlelanding/2026/dd/d5dd00387c):
  retrospective medicinal-chemistry prioritization data, but not complete
  counterfactual candidate or causal selection-policy evidence.

### Why each candidate is not a renamed standard method

| Candidate | MAML/ANIL/MetaDTA/AdaMBind boundary | Bayesian/GP/kernel boundary | Retrieval boundary | Fine-tuning boundary |
| --- | --- | --- | --- | --- |
| ADRO | jointly learns measurement acquisition and a hard finite intervention; no inner gradient | no analytic posterior/function solve; active GP is a baseline | cannot access donor outcomes or choose source episodes | recipient parameters remain frozen |
| EASPI | compiles an exact finite molecular program, not a gradient-adapted predictor or generic context regressor | no similarity-weighted continuous prediction; categorical Bayes is a baseline | an applicable source-frozen rule is executed on an unseen transformation context rather than copying a neighbour label | no recipient parameter update |
| CSTO | updates an explicit campaign state using action/assay semantics rather than model parameters or task scheduling | no fixed likelihood/covariance update | does not retrieve historical campaigns | transition law is frozen at meta-test |
| SDREA | emits a bounded action rather than adapting predictor weights, but remains close to meta-LTR | no continuous posterior or kernel smoother | no donor lookup | no SGD; novelty still fails if generic meta-LTR matches it |
| UGRMC | chooses whether/which intervention is safe rather than adapting a CNP/MAML head | no analytic posterior; an amortized-posterior interpretation is a major rejection risk | no source outcome lookup | all model weights remain fixed |

### Candidate-specific validation contracts

| Candidate | Baselines | Required ablations | Negative controls | Leakage checks | Pass criterion | Kill criterion |
| --- | --- | --- | --- | --- | --- | --- |
| ADRO | random/diversity/D-optimal/uncertainty/active KRR-GP acquisition; frozen base; MAML/ANIL/MetaDTA/CNP/MetaFun | no policy, no state update, no dictionary, no protein routing, no NULL gate | wrong support, shuffled labels, matched donor, k>=3 assignment derangement, k=1 sign/null/norm, label-noise dose | candidate-pool arm classifier, query-label firewall, target/homology/document/assay closure, pool permutation/size | learned acquisition and response operator beat strongest equal-budget baselines; absolute ranking and affinity improve; Delta_label/Delta_assign pass | no active-oracle headroom, random-equivalent policy, failed assignment gate, or non-load-bearing constraints |
| EASPI | copy-support-sign, Matsy/MMP/MMS, mixed effects, categorical empirical Bayes, graph-difference, learned kernel | no protein routing, no assay routing, unrestricted executor, no UNSUPPORTED state | rule/sign/context/transformation shuffle and wrong-assay pair | exact same-assay provenance, source-frozen vocabulary, target/homology split, query applicability audit | local k=3/5 gate passes and finite grammar beats all classical rule baselines on new components | insufficient coverage/power, unstable rule sign, baseline equivalence, or k=1 adaptation claim |
| CSTO | chronological CNP/RNN, frozen base, no-state and pooled temporal models | remove decision, assay, batch, chronology or selection head one at a time | time reversal, selection shuffle, assay shuffle, outcome-null | project-disjoint split, contemporaneous candidate pool, timestamp and propensity audit | prospective future-round ranking/affinity improves and outcome-driven state remains load-bearing | publication order only, MNAR unaddressed, chronology ablations inert, or cross-project failure |
| SDREA | frozen ranking, equal-capacity generic meta-LTR, PiRank-style and direct point-score adapters | remove STOP, edit budget or action restriction | wrong/permuted/null support and action-label shuffle | query permutation, distractor, subset, library-size and tie-breaking audit | positive B=1/2 oracle headroom and stable absolute gain beyond generic meta-LTR | passive information gate fails, candidate-set instability, or constraints do not matter |
| UGRMC | fixed uncertainty threshold, support-free uncertainty, selective CNP, conformal/selective prediction | no mode, no hard gate, no risk supervision, no query dependence | shuffled/wrong/null support and noise-dose series | fixed-coverage comparison, target/protein shuffle, OOD calibration on locked components | at equal coverage, correct evidence improves ranking/affinity and mode/gate ablations reduce it | gain comes only from lower coverage, label assignment is ignored, or uncertainty does not predict harm |

---

## Part 4 - Selected mechanism: ADRO

### One-sentence paper contribution

> Learn from abundant source targets how to select a budgeted set of diagnostic
> compounds and translate their measured outcomes into a finite, abstaining,
> query-dependent response intervention for a strictly unseen target.

### Core modules

The load-bearing mechanism consists of: a frozen support-free DTA interface, a
legal candidate-pool contract, sequential diagnostic acquisition policy,
assignment-preserving evidence encoder, recurrent low-dimensional state update,
finite response-mode dictionary, structural NULL/harm gate, and bounded query
intervention executor. The evaluator and counterfactual objective establish
attribution but are not claimed as additional model innovations.

### Why this is the highest-potential hypothesis

1. It directly addresses the v2 failure: passive support did not admit stable
   label-assignment information, so ADRO changes information acquisition rather
   than increasing estimator capacity.
2. Its adaptation object is finite and auditable; k measurements discriminate
   modes instead of estimating an arbitrary function.
3. The contribution is experimentally falsifiable: acquisition advantage,
   mode separability, assignment specificity and downstream ranking can each
   fail independently.
4. It preserves the base DTA predictor through an exact support-free/no-op path.
5. It reflects an actual drug-discovery decision: which compound should be
   measured next to resolve uncertainty about a target-specific response.

### Core mathematical operator

For source task T_i=(C_i,Q_i,p_i), policy roll-out produces

\[
S_i^\pi=\{(x_j,y_j)\}_{j=1}^{k},
\qquad x_j\sim\pi_\theta(h_{j-1},p_i,C_i\setminus x_{<j}).
\]

The state update and intervention are

\[
h_j=U_\theta(h_{j-1},\phi(p_i,x_j),y_j-f_0(p_i,x_j)),
\]

\[
\hat y_q=f_0(p_i,x_q)+
g_{safe}(h_k)\sum_{m=1}^M g_m(h_k)R_m(p_i,x_q).
\]

The source meta-objective is

\[
\theta^*=\arg\min_\theta\sum_i
\Big[
L_{rank}(\hat y_{Q_i},y_{Q_i})
+\lambda_a L_{aff}(\hat y_{Q_i},y_{Q_i})
+\lambda_{cf}L_{assignment}
+\lambda_h L_{harm/abstain}
\Big],
\]

with a hard k-label budget. Exact coefficients and thresholds must be
preregistered from source development; they are not a broad tuning dimension.

### Why existing methods cannot substitute for the claim

- MAML/ANIL/AdaMBind optimize how parameters or tasks are updated after labels;
  they do not jointly learn which recipient measurements identify a finite
  response intervention.
- KRR/GP/Bayesian active learning selects points to reduce model uncertainty and
  performs analytic function inference; they are strong acquisition baselines,
  not the same learned intervention object.
- CNP/MetaFun can approximate the mapping, but do not by themselves establish
  a hard diagnostic budget, finite response modes, exact NULL behavior or
  acquisition causality. If a capacity-matched generic CNP matches ADRO and the
  constraints can be removed without loss, the novelty claim fails.
- Retrieval selects historical analogues; ADRO may not query source labels or
  copy donor outcomes at meta-test.
- A stronger backbone improves f_0, not the learned law for adapting a new task.

### Non-negotiable caveat

**FACT:** ADRO has not passed an active information gate. It is the highest-
potential hypothesis, not the current winning model. If active selection does
not make label assignment informative, the correct decision is to abandon this
task/data contract.

---

## Part 5 - Module-level implementation plan, no code

### Phase 0 - Data and information admission

1. **Candidate-pool contract.** Define a label-free C_t and disjoint Q_t for
   every source target; selected compounds leave the query evaluation. Freeze
   target/homology/document/assay/parent/scaffold audits and candidate-pool
   sizes before labels.
2. **Retrospective-claim boundary.** Current ChEMBL pools test budgeted
   diagnostic selection only. Do not claim campaign policy or causal DMTA.
3. **Active oracle gate.** Before a neural policy, measure whether any legal
   k-step selection can improve mode separability and ranking over random,
   diversity and D-optimal selection.
4. **Power/MDE.** Use target/homology/provenance components as independent
   units and freeze the practical effect threshold before the locked role.

### Required new modules

1. **Frozen Base Interface**
   - exposes f_0(p,x) and fixed ligand/protein representations;
   - guarantees bit-identical support-free predictions before/after meta-training;
   - generates full encoder/head component-OOF residuals for source episodes.

2. **Candidate-Pool and Task Contract**
   - represents legal unmeasured compounds and masks selected/query items;
   - prevents query labels, future metadata and candidate-pool identity leakage;
   - supports permutation and library-size tests.

3. **Diagnostic Acquisition Policy**
   - selects one compound at each of k steps;
   - is supervised through downstream source-query utility and mode
     discriminability;
   - cannot retrieve source labels or use recipient gradients.

4. **Support Evidence Encoder**
   - encodes measured OOF residual, compound identity/relation, protein and
     admissible assay context;
   - preserves compound-label assignment and is permutation equivariant.

5. **Recurrent Adaptation-State Updater**
   - updates a bounded low-dimensional evidence state after each measurement;
   - exposes state trajectories for null, sign-flip and label-noise audits.

6. **Finite Response-Mode Dictionary**
   - contains a small fixed M of bounded query-dependent intervention
     primitives;
   - excludes pure calibration-only modes from the innovation claim;
   - freezes M and capacity before locked evaluation.

7. **Structural Abstention Gate**
   - forces exact Delta=0 for empty, residual-null, unsupported or insufficient
     evidence;
   - learns intervention harm, not merely predictive uncertainty.

8. **Query Intervention Executor**
   - composes sparse modes with f_0;
   - produces inductive per-compound scores when claimed, or explicitly narrows
     the claim to fixed-library transductive ranking if candidate-set tests fail.

9. **Meta-Ranking and Counterfactual Objective**
   - jointly scores absolute improvement over f_0 and assignment specificity;
   - uses chemistry-fixed derangement for k>=3 and sign/null/norm controls for
     k=1;
   - includes a harm cost so wrong evidence cannot be rewarded for destructive
     predictions.

10. **Locked Component Evaluator**
    - reports target-macro CI, NDCG@10 and pairwise proper log loss as primary;
    - reports Spearman, RMSE, MAE, acquisition coverage, abstention and
      conditional harm as secondary;
    - averages draws/seeds within target before paired component bootstrap.

### Minimum experiment ladder

1. Frozen support-free base, ligand-only base and high-data target oracle.
2. Random, diversity, D-optimal, uncertainty and active KRR/GP acquisition.
3. Calibration, ridge, KRR and closed-form Bayesian/MDK posterior.
4. Fine-tuning, MAML, ANIL, MetaDTA, AdaMBind-style scheduler, CNP/ANP,
   MetaFun and current CMAL.
5. ADRO without acquisition, without state update, without mode dictionary,
   without protein routing and without abstention.
6. Wrong support, shuffled labels, chemistry-matched donor, k>=3
   chemistry-fixed assignment permutation, k=1 sign/null/norm controls,
   protein/target shuffle and label-noise dose response.
7. Candidate-pool classifier, query permutation, distractor insertion, query
   subset and library-size shift.

### Success criteria

ADRO is admitted only if all are true on source development and then once on a
separately frozen locked-source role:

1. Active-selected support has positive Delta_label and k=3/5 Delta_assign
   component-CI lower bounds.
2. Correct evidence produces absolute positive target-macro CI/NDCG and lower
   RMSE/MAE relative to the identical f_0.
3. ADRO beats the strongest acquisition and adaptation baseline under the same
   candidate pool, k and representation budget.
4. Removal of acquisition, state update, finite modes or abstention materially
   reduces the effect.
5. Label-noise produces monotonic degradation; wrong/permute/null controls do
   not improve predictions.
6. Protein removal/shuffle weakens the claimed protein-conditioned component;
   otherwise the claim is narrowed and cannot call the adaptation protein-
   conditioned.
7. Candidate-set and leakage audits pass, and results are stable across the
   preregistered seeds and independent components.

### Stop rules

- No active-oracle headroom: stop before policy implementation.
- Learned acquisition equals random/diversity/D-optimal/active GP: reject the
  acquisition contribution.
- Active supports fail Delta_label or k=3/5 Delta_assign: abandon the active A2S
  episode construction.
- Correct support beats wrong support but not f_0: scientific failure.
- Generic CNP/meta-LTR matches ADRO and structural constraints are not
  load-bearing: novelty failure.
- Only ranking or only calibration improves: narrow the claim; the requested
  combined ranking/affinity objective is not achieved.
- Candidate-pool leakage or future metadata dependence: invalidate the run.
- Locked-source failure: do not open recipient labels.

## Final PI decision

**HYPOTHESIS:** ADRO is the best research direction because it attacks the
information bottleneck directly and makes the adaptation object finite. EASPI
is the correct low-cost current-data fallback if a label-free MMP coverage
census passes. CSTO requires a different project-level dataset. SDREA is not
novel enough, and UGRMC belongs inside ADRO as a safety mechanism.

**FACT:** No mechanism is presently admitted, no new model code should be
written, and no current `model/` implementation is ready for publication. The
next valid result may still be a source-only STOP.
