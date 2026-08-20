# CIIP-2 Research Report: Toward a Deployable Protein-Conditioned Interaction Potential

Date: 2026-08-20
Author role: independent research scientist (successor cycle)
Status: Phase 0-2 complete (audit, diagnostics design, candidate comparison);
Phase 3 preregistration issued separately; no production code touched.
Scope guards: no training was run for this report; all numbers below are either
(a) quoted from frozen CIIP-1A artifacts, or (b) computed read-only from those
artifacts (scripts inline in the session log, no new state written to stages).

---

## 1. Formal problem definition

### 1.1 Data contract (frozen, verified)

- Panel: Duong-Ly kinase interaction panel, 97 constructs x 183 ligands,
  raw single-concentration % inhibition (never relabeled).
- Admitted single-point WT/variant pairs: 65; ESM-covered subset: 49
  (12 parents have >= 2 single-mutant pairs; 6 parents have >= 2 covered
  TRAIN pairs: ABL1 4, KIT 5, EGFR 3, FGFR3 2, PDGFRA 2, RET 4).
- Endpoint target (frozen in CIIP-1A):
  d_p(L) = y(P_v,L) - y(P_w,L);  c_p(L) = d_p(L) - mean_L' d_p(L').
  Centering removes the ligand-invariant mutation-wide shift (mutation main
  effect). What remains is the ligand-dependent interaction component plus
  noise. % inhibition at one concentration is a functional competition
  readout; it is NOT Ki/Kd/pK/DeltaDeltaG and is never interpreted as such.
- Current split: pair-level 60/20/20 stratified per parent (39/13/13), so the
  same parent appears in train and test. CIIP-1B (parent-disjoint) has not run.

### 1.2 Deployment object

A single scalar field s_theta(P, L) defined for any construct sequence P and
ligand L, with every scientific contrast derived from its finite differences:

    protein contrast      g(P_w->P_v, L) = s(P_v,L) - s(P_w,L)
    ligand contrast       g(P, A->B)     = s(P,B)   - s(P,A)
    double contrast       D_hat          = cross-difference of the above

Antisymmetry, identity-zero, and cycle consistency hold by construction.
Deployment input: full construct sequence + ligand. Mutation coordinates,
MSAs, structures, and target IDs are never deployment inputs.

### 1.3 Endpoint semantics separation (domain A hygiene)

| concept | measured here? | proxy |
|---|---|---|
| mutation effect on stability/folding | no | - |
| mutation effect on kinase function | partially (one dose, competition readout) | y itself |
| catalytic activity change | no | - |
| ATP-competitive vs allosteric mode | not separable from this panel | tracer assay context only |
| direct binding affinity (Ki/Kd) | NO | single-dose % inhibition only |
| ligand-dependent interaction change | YES (the estimand) | c_p(L) |

Claims in this programme never cross these boundaries.

### 1.4 Structural decomposition under study

    y(P,L) = mu + m_P(P) + m_L(L) + I(P,L) + a_L(y) + eps

where a_L is the ligand-specific assay transform (single-dose link) and eps
the assay noise. Within-pair ligand-centered contrasts cancel m_P and m_L
exactly and leave Delta I_p(L) + Delta eps. Feature-side confounding (any
regressor correlated with parent identity) is the residual identification
threat and is handled by design (Section 6).

---

## 2. What CIIP-1A established and did not establish (Phase 0 audit)

### 2.1 Frozen results (verified against artifacts)

Seven-arm control experiment (prereg 39d02166..., result 55b0de12...,
adjudication ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED), 49 covered pairs,
pair-level split, 9 test pairs / 6 parents, single seed:

| arm | R2 | Spearman | dead-zone sign | nonconstant | slope med |
|---|---:|---:|---:|---:|---:|
| oracle_local_esm_correct | +0.0075 | +0.331 | 0.702 | 9/9 | 0.466 |
| family_preserving_shuffle | -0.0036 | -0.007 | 0.498 | 9/9 | 0.907 |
| random_local_window | +0.1291 | +0.320 | 0.705 | 9/9 | 1.076 |
| ligand_only | +0.0000 | nan | 0.494 | 0/9 | nan |
| ligand_invariant_shift | +0.0000 | nan | 0.494 | 0/9 | nan |
| random_protein | +0.0004 | +0.001 | 0.493 | 9/9 | 2.085 |
| free_pairwise | +0.0265 | +0.125 | 0.562 | 9/9 | -0.009 |

Contextual propagation audit (read-only, cdd6e0a8...): mean L2 delta norms
mutation site 4.011, radius-6 window 1.161, non-site context 0.057, full
sequence 0.074; mutation erasure (X-substitution) exactly zeroes residue
deltas for 49/49 pairs.

### 2.2 Proven

1. ESM-2 residue states sense the mutation (49/49 pairs; site delta 0.531 vs
   random-window 0.027). Representation-level mutation sensitivity is real.
2. The signal is localized (site >> local window >> distant context) but
   propagates: a distant window is not mutation-null.
3. The oracle mutation-centered window feature has no detectable
   centered-response increment over matched controls under the frozen
   low-capacity potential and pair-level split (correct - random = -0.1217,
   parent bootstrap [-0.4569, +0.0327]).
4. KLIFS pocket one-hot collapses structurally (nonzero on 3/9 test pairs);
   ESM restores nonconstant outputs (9/9) but not accuracy.

### 2.3 NOT proven (and now measured by this audit)

1. Whether ANY deployable sequence feature carries protein-conditioned
   interaction information (no arm used a deployment-legal input that beat
   the shared ligand pattern; see 2.4).
2. Parent-disjoint behavior (split not parent-disjoint).
3. Binding affinity relevance (endpoint is functional, single-dose).
4. Absence of mechanism (n=9 test pairs; bootstrap CI half-width ~0.2-0.4;
   the study is a coarse filter, not proof of absence).

### 2.4 New read-only findings of this audit (computed 2026-08-20)

(a) The ligand_only arm is structurally constant. Because s = alpha(P)^T
psi(L) and protein features were zeroed, Delta_P s(L) = 0 identically
(0/9 nonconstant). The arm therefore understates the ligand-only baseline
rather than measuring it.

(b) The true shared-ligand-pattern baseline. A train-pair mean centered
profile m(L) achieves cell R2 = 0.1313 on the same 9 covered test pairs
(median profile: 0.1004). This matches the random-window arm (0.1291).
Mechanistic reading: the random-window arm succeeds through a stable,
low-variance Delta alpha that keys a rank-8 ligand-response basis; what it
generalizes is the panel-shared resistance pattern, not distributed mutation
information. The mutation-site window is idiosyncratic per pair and dilutes
that shared pattern (hence 0.0075 < 0.1291). The verdict is unchanged, but
the interpretation sharpens: no arm demonstrated protein-conditioned value
beyond the shared ligand pattern.

(c) Hierarchical variance decomposition of the covered interaction signal:
    parent-shared component 134.8%^2, mutation-specific residual 89.7%^2
    (40% of total, noise-inclusive; no replicates exist to separate them).

(d) Same-parent sibling structure is real: predicting a pair's centered
    profile from its siblings (leave-one-sibling-out) gives per-pair median
    R2 = 0.293 (IQR [-0.168, 0.557]; sibling profile correlation mean 0.386
    vs random-pair 0.024). Under the CIIP-1A pair-level split this ceiling
    was legally reachable (siblings in train); no arm approached it.

(e) Family-level transfer of the shared component fails trivially: a
    family-mean profile (Manning groups; here dominated by TK) predicts
    held-out parents at cell R2 = -0.021. So the parent-shared interaction
    structure does not transfer by category priors; whether it transfers by
    sequence features is exactly the CIIP-1B deployment question.

(f) Assay structure: 23.0% of panel cells fall outside [0,100]; the WT panel
    is heavily ceiling-loaded (per-ligand WT mean median 90.9%; 99/183
    ligands with WT mean > 90; 0 ligands < 10; per-ligand WT sd across
    parents median 17.3). Single-dose contrasts against near-ceiling WT
    activity are one-sided-compressed: interaction sensitivity is
    concentrated in the ~84 mid-zone ligands. This is an assay-link effect
    (a_L above), not a biological interaction property, and must be modeled
    or stratified, not averaged over.

(g) Power: with 9 test pairs, R2 differences below ~0.3 are not resolvable
    per split (bootstrap CI half-widths from the frozen run: 0.2-0.4). Any
    successor claim must aggregate over seeds and splits and predeclare
    graded interpretations.

### 2.5 The actual bottleneck, restated

The bottleneck is not mutation sensitivity (proven present) and not
nonconstant output (proven achievable). It is that the current single
bilinear term entangles two different components with different transfer
behavior: (i) a panel-shared ligand-response pattern (ligand-side,
parent-free, R2~0.13) and (ii) a protein-conditioned deviation (parent-keyed
structure worth ~+0.16 more for seen parents, currently unpredicted). A
deployable, transferable interaction learner must factor these explicitly,
then test whether the deviation part is predictable from sequence under
parent-disjoint evaluation. The failure of the oracle arm is explained by
this entanglement, and no larger encoder fixes it.

---

## 3. Literature map (Phase 2A/2B) and what each field contributes

Inspected anchors are marked [I]; items marked [K] are well-established
results used as background and not load-bearing for any claim below.

### 3.1 Bioinformatics / protein modeling (domain A)

- ESM-1v / ESM-2 zero-shot mutation effects [I]: masked-marginal likelihood
  ratios predict variant effects without supervision (Meier et al., NSMB
  2021; Lin et al., Science 2023). Relevance: the per-residue "surprise"
  channel log p(a_i | sequence) is a DEPLOYABLE per-sequence feature that is
  mutation-sensitive without mutation coordinates.
- ProteinGym [I] (Notin et al. 2022/2024): substitution-effect benchmark;
  shows PLM scores are strong but endpoint-specific; functional-assay
  readouts are not affinity readouts.
- SaProt [I] (Su et al., ICLR 2024): structure-aware vocabulary (Foldseek
  tokens) improves mutation-effect and binding-related tasks; deferred to a
  later arm (adds a structure pipeline dependency).
- MSA Transformer / coevolution [K]: conservation and coevolution carry
  interaction-relevant constraints; MSA availability for a deployable
  construct is not guaranteed; kept as an optional representation arm.
- PremPLI [I] (Wang et al. 2022): physics-based mutation effects on
  protein-ligand binding; requires structures and predicts DeltaDeltaG-like
  affinity changes - a different estimand, kept separate.
- Platinum [I] (Witte et al., NAR 2015): ~1k measured mutation-ligand
  affinity changes; future Ki/Kd bridge candidate, not usable for the
  functional endpoint here.
- MdrDB [I] (2023, Commun Chem): mutation-drug-resistance benchmark; pooled
  endpoints across assays; not pooled with ours (endpoint-mixing ban).
- eSIG-Net [I] (Nat Methods 2026): interaction language model decoding
  single-mutation effects on PPI networks; uses mutation-site-focused
  encoding plus joint task learning; supports the privileged-teacher idea
  (oracle-site information accelerates learning) but does not address our
  deployment constraint; cited as the closest mutation-aware interaction LM.
- Cold-start DTA evaluation practice [I/K]: entity-disjoint splits are the
  accepted standard (Davis/BindingDB lines, CS-DTA critiques); current
  pair-level split does not meet it; CIIP-1B does.
- AdaMBind [I] (Nat Commun 2026): meta-learning, task-adaptive few-shot DTA;
  evidence that support-conditioned adaptation is the current frontier for
  small-sample DTA; informs candidate C5 (no test-time gradients there
  either).

### 3.2 Vision / multimodal (domain B)

- Conditional computation and query-dependent routing [K] (MoE, Perceiver
  cross-attention): the mechanism that lets one deployable trunk consult
  different input regions per query (here: per ligand).
- Multiple-instance learning and weakly-supervised localization [I] (Ilse et
 . 2018 attention-MIL and successors): bag-level labels, instance-level
  attention; the exact formal analogy for sequence-level interaction labels
  with residue-level mechanisms. Attention weights are routing statistics,
  never contact probabilities (binding constraint).
- Privileged information / LUPI [I] (Vapnik & Vashist 2009; Lopez-Paz et al.
  2016 unifying distillation and LUPI): teacher uses oracle inputs at
  training only; student is deployable; the distillation gap measures
  recoverability of privileged information.
- Disentangled and counterfactual representations [K]: encourage factorized
  codes; here the factorization is enforced structurally (main effects vs
  interaction) rather than by regularizers.
- IRM [I] (Arjovsky et al. 2019) and its critiques (Rosenfeld et al. 2020):
  invariance learning needs many environments; with 12-21 parents the
  environment set is too small for IRM-style training - REJECTED as a
  training method, retained as an evaluation lens (per-parent stability).
- Group DRO [I] (Sagawa et al. 2020): worst-group objectives; usable as an
  optional training variant for per-parent robustness once >10 training
  parents exist.
- Gradient-free test-time adaptation [K] (BN-statistic adaptation, feature
  realignment): forbidden-adjacent in spirit; we adopt NONE (no test-time
  adaptation of any kind, gradient or not, beyond deterministic forward).

### 3.3 NLP / meta-learning (domain C)

- Conditional/Attentive Neural Processes [I] (Garnelo et al. 2018; Kim et
  al. 2019): amortized task adaptation from a support set with no test-time
  optimization; the correct vehicle for few-shot DTA support conditioning
  (candidate C5).
- Set Transformer [I] (Lee et al. 2019): permutation-invariant attention
  pooling for support sets.
- Prototypical / FEAT / Matching networks [I] (Snell 2017; Ye 2020; Vinyals
  2016): task embeddings, task-adapted metrics; informs how a parent-panel
  support set is summarized into a task token.
- Hypernetwork task conditioning [K] (HyperFormer): weights generated from
  task context; heavier than needed here; deferred.
- Task-invariant representations under shift [K]: motivates reporting
  per-parent error stability rather than pooled aggregates only.

### 3.4 System identification / causal inference / SciML (domain D)

- Functional ANOVA and orthogonal decompositions [I] (Stone 1994; Hooker
  2004-2007; recent exact fANOVA for categorical inputs): main effects and
  interactions are identifiable under orthogonality constraints; with
  correlated features the decomposition is convention-dependent - hence
  enforcing it architecturally (centering layers), not post hoc.
- Double/debiased ML and Neyman orthogonality [I] (Chernozhukov et al. 2018;
  Newey & Robins 2018): cross-fitted nuisances make interaction estimation
  robust to first-stage errors (second-order sensitivity). Basis of the
  cross-fitted ligand-pattern residualization (candidate C3).
- CATE meta-learners [I] (Kuenzel et al. 2019; Nie & Wager 2020 R-learner):
  mutation-as-treatment framing; our ligand-only nuisance is exactly the
  outcome model; the centered objective is the treatment-contrast loss.
- Hierarchical Bayes shrinkage / partial pooling [K]: the correct analysis
  lens for parent-level effects (our variance decomposition implements its
  frequentist limit); a full HB panel model is kept as an analysis ceiling
  tool, not a deployment model.
- Tensor factorization / panel factor models [I]: shared-pattern structure
  extraction; again analysis tooling (we used profile-mean and sibling
  ceilings), not deployable sequence functions.
- Integrability / cycle consistency [K]: requiring contrasts to come from
  one scalar field is the discrete analogue of conservative-vector-field
  constraints; already enforced by the potential parameterization.

---

## 4. Candidate mechanisms (Phase 2C): four classes plus two rejected

Notation: h_i(P) in R^d are frozen ESM-2 residue states of construct P;
z_L is the ligand feature (ECFP4 2048 now; drop-in replaceable); Lambda_p
is pair p's ligand panel; Delta_P s(L) = s(P_v,L) - s(P_w,L).

### C1. Ligand-Conditioned Residue Router (LCRR) - deployable mechanism

Math:
    q(L)  = MLP_phi(z_L)                          ligand query
    a_i   = softmax_i( (W_k h_i)^T q(L) / sqrt(d) )  residue routing
    alpha(P,L) = sum_i a_i(L,P) W_v h_i  in R^r   interaction coordinate
    beta(L)    = MLP_psi(z_L)  in R^r             ligand interaction basis
    s(P,L)     = alpha(P,L)^T beta(L)
Prediction of the centered contrast:
    c_hat_p(L) = Delta_P s(L) - mean_L' Delta_P s(L')

Why it is not the failed bilinear form: in CIIP-1A s = alpha(P)^T psi(L)
has a ligand-INdependent alpha, so Delta_P s(L) = (Delta alpha)^T psi(L)
is a single separable bilinear term - it cannot express "this mutation
changes the response of ligand class A but not class B" beyond one shared
inner product. With LCRR, alpha depends on L (conditional computation), so
Delta_P s(L) can be nonseparable in (mutation, ligand).

1. Math definition: above. Rank r in {4,8,16}; heads in {1,4}; frozen ESM.
2. Inputs: construct sequence (any construct, WT or variant), ligand
   feature. Output: scalar s(P,L).
3. Transferable carrier: the routing function (which residue states a
   ligand consults) and the shared interaction basis beta. The carrier is
   ligand-conditioned pooling weights, not entity embeddings; no target ID.
4. Main-effect removal: within-pair, per-panel centering of predictions
   (mutation main effect cancels); ligand main effect never enters the
   centered contrast by construction; explicit b_L head exists only for the
   raw-endpoint auxiliary loss.
5. Mutation-free deployment: native - the router sees only the sequence.
6. Few-shot support/query: support ligands enter as conditioning tokens
   via C5; the router itself needs no support.
7. Zero/few-shot DTA: zero-shot = s(P,L) directly; few-shot = C5 residual
   head on top of s.
8. Training objective: weighted centered MSE
   sum_p sum_L in Lambda_p w_L (c_p(L) - c_hat_p(L))^2 + lambda_wd ||theta||^2
   with assay-gain weights w_L (Section 5.3); optional teacher term (C4).
9. Shortcuts: (a) router collapse to conserved motifs -> family-keyed
   shared-pattern prediction (detected by parent-disjoint eval and
   family-prior baseline); (b) ligand-similarity leakage (detected by
   ligand-permutation and ligand-only nuisance); (c) memorizing parents
   (detected by parent-disjoint split); (d) attention maps over-read as
   mechanism (forbidden; routing stats reported as stats only).
10. Negative controls: within-pair ligand-label permutation (train-time and
    eval-time), random-protein swap, mutation-erasure invariance (X-out),
    same-parent wrong-mutation assignment, ligand-invariant-shift arm.
11. Data needed: current frozen contract (49 covered pairs) suffices for
    smoke and CIIP-1A-scope claims; CIIP-1B needs the parent-disjoint split
    (12 multi-pair parents; ~5-6 held out).
12. Compute: model 0.2-2 M params; 49x183 ~ 9k cells; minutes/seed on GPU;
    full multi-seed multi-split grid < 2 h GPU; < 2 GB memory.
13. Compatibility: drops into the existing UnifiedPotential as a new s
    parameterization; contrasts g(.,.) and D_hat unchanged; production
    model/ untouched (successor stage dir only).
14. Minimal runnable version: single-head, r=8, frozen ESM, mean-pooled
    ablation arm + router arm, 1 seed, pair-level split, control arms
    wired.
15. If it fails: (a) router == mean-pool everywhere -> conditional routing
    unnecessary at this scale (a real finding: sequence length of useful
    routing evidence is too short); (b) parent-disjoint R2 not above the
    ligand-pattern prior -> sequence-level interaction transfer NOT
    SUPPORTED on this endpoint/panel (the central question answered
    negatively at this power); (c) succeeds only with oracle site input ->
    interaction evidence exists but is not deployable (privileged-only).

### C2. Orthogonal Interaction Decomposition (OID) - constraint layer

Math: decompose the raw endpoint model
    y_hat(P,L) = mu + b_P(P) + b_L(L) + s(P,L)
with architectural orthogonality enforced in the forward pass at loss time:
    s_bar(P,L) = s(P,L) - mean_{L in Lambda_p} s(P,L)          (panel axis)
    s_tilde(P,L) = s_bar(P,L) - mean_{P in B} s_bar(P,L)
                 + mean_{P in B, L in Lambda} s_bar             (protein axis)
s_tilde is the only term that touches centered targets. Batch-level protein
axis uses only training constructs (leakage guard). This is fANOVA-style
orthogonalization implemented as mean-subtraction layers (differentiable,
no closed-form parameter solve), so interaction cannot reabsorb main
effects BY CONSTRUCTION rather than by regularization.
1-3: same inputs as C1; carrier = s_tilde only.
4. Main-effect removal: structural (above).
5-7: unchanged from C1 (it is a layer, not a model).
8. Objective: same weighted centered MSE on s_tilde contrasts; raw-endpoint
   auxiliary loss trains b_P, b_L, mu (never mixed into interaction claims).
9. Shortcut: batch-composition dependence of the protein-axis centering -
   mitigated by fixed construct-grouped batches, reported as a caveat.
10. Controls: verify b_P/b_L absorb >99% of raw-endpoint trainable variance
    before s_tilde unlocks (grad-gate diagnostics, as in CIIP-1A audit).
11-12: zero extra data/compute.
13. Compatible with any s parameterization including C1.
14. Minimal version: centering layers + variance-accounting diagnostics.
15. If it fails: interaction R2 drops when orthogonality is enforced ->
    prior nonzero results were main-effect contamination (a decisive
    negative-control outcome).

### C3. Cross-Fitted Orthogonal Interaction Estimation (CFOIE) - estimation layer

Math (R-learner analogy, mutation-as-treatment):
    Stage A (nuisance): ligand-pattern model m_hat(L), gradient-trained,
       K-fold cross-fitted OVER PARENTS (fit K-1 folds, predict held-out).
    Stage B: residual r_p(L) = c_p(L) - m_hat^{(-p)}(L); interaction loss
       sum w_L (r - c_hat)^2.
Neyman orthogonality: first-order insensitivity to m_hat errors, so the
"must beat the ligand pattern" requirement is embedded in estimation, not
only in evaluation. No closed-form solver anywhere; both stages trained by
SGD.
1-3: unchanged; carrier = residual interaction only.
4. Main-effect removal: ligand-pattern (shared) component explicitly
   removed before fitting; protein main effect already canceled in target.
5-7: unchanged.
8. Objective: Stage-B weighted centered MSE.
9. Shortcut: if m_hat is too weak, s relearns the shared pattern (monitor
   via ||c_hat - m_hat|| attribution); if too strong (overfit folds), real
   interaction gets residualized away (monitored by foldwise residuals).
10. Controls: fold-concordance of m_hat; permutation of pair identity
    during Stage A (must not change Stage-B ceiling).
11-12: negligible extra compute (Stage A is a small MLP on 183 ligands).
13. Compatible; recommended default objective.
14. Minimal version: K=3 folds by parent, m_hat = 2-layer MLP on ECFP.
15. If it fails: s cannot beat cross-fitted residual zero -> the shared
    pattern is all there is at this power (a clean NOT_SUPPORTED verdict
    for protein-conditioned transfer).

### C4. Privileged Mutation-Site Teacher Distillation (PMSTD) - optional auxiliary

Math: teacher T = C1 with an augmented input (oracle site index appended to
residue states, or mutation-window delta features); trained on train pairs;
student S = deployable C1; joint loss
    L = L_label(S) + lambda_d ||Delta_P s_S - stopgrad(Delta_P s_T)||^2.
Deployment uses S only. Distillation gap (teacher-student divergence on
held-out parents) measures how much site-privileged information is
recoverable from sequence - itself a scientific readout (LUPI framing).
1-3: carrier = routing/basis of S; T is training-only.
4. Main effects: unchanged (contrasts).
5. Mutation-free deployment: yes (S never sees coordinates).
6-7. few-shot/DTA: unchanged.
8. Objective: label + distillation terms.
9. Shortcut: student mimics teacher's shortcuts -> mitigate by distilling
   only the centered contrast, not raw scores; report teacher's own control
   compliance.
10. Controls: teacher with wrong-site privilege (must degrade); student
    with lambda_d = 0 (= C1).
11. Data: same contract; teacher arms are diagnostics.
12. Compute: 2x training cost.
13. Compatible.
14. Minimal version: teacher = C1 + site-index channel.
15. If it fails: gap stays large -> site information not recoverable from
    sequence features (bounded deployability); if teacher itself fails
    controls -> oracle-site interaction signal absent on this endpoint.

### C5. Support-Conditioned Neural Process (SCNP) - few-shot bridge (Phase 6)

Math: for a parent pair p with support S = {(L_j, c_j)}_{j<=k}:
    t_p   = SetTransformer({ (z_Lj, c_j) })                task token
    h_P   = protein pair encoding (from C1 trunk, frozen or fine-tuned)
    c_hat(L_q | S) = ANP cross-attention(h_P, t_p, z_Lq)
    + Gaussian uncertainty head -> abstention when sigma > tau.
No test-time gradients; k=0 reduces to the s_theta prior. Used only after
deployment diagnostics pass; evaluated in the P-line few-shot protocol
(k in {0,1,2,3,5,10,20,40}) with frozen selection rules.
Shortcut: support leakage of query labels (forbidden by construction);
uncertainty miscalibration (reported via calibration curves).
If it fails: amortized adaptation adds nothing over s_theta prior ->
interaction field already captures transferable structure or none exists.

### Rejected or deferred (with reasons)

- IRM-style invariance training: needs many environments; <= 21 parents.
  REJECT as training method.
- Tensor / hierarchical-Bayes panel factor models as the deployed s_theta:
  not sequence-conditioned functions (cannot score unseen constructs);
  retained as analysis ceilings only.
- Free pairwise transformation heads as production: violate potential
  integrability (already diagnostic-only).
- Structure-aware encoders (SaProt, AlphaFold tokens): deferred to a
  labeled optional arm after sequence-only mechanisms are resolved.
- Any test-time adaptation (gradient or statistic-based): forbidden.
- Ridge/kernel/closed-form solvers anywhere in the deployed path:
  forbidden; all components SGD-trained.
- Pooling more panels (Davis Kd + Duong-Ly): endpoint mixing forbidden.

---

## 5. Feasibility, risks, and the negative-control matrix

### 5.1 Data feasibility (census, verified against frozen artifacts)

| requirement | status |
|---|---|
| 49 covered pairs for CIIP-1A-scope | available |
| parent-disjoint CIIP-1B split | formable: 12 parents with >=2 pairs |
| same-parent wrong-mutation test | train-feasible: 6 parents with >=2 covered train pairs (20 pairs with siblings) |
| within-pair ligand permutation | feasible: ~183 ligands/pair, many permutations |
| mutation-erasure predictive test | feasible: X-substitution verified exact (49/49); requires ESM forward on erased sequences (cache + ~67 forwards, minutes) |
| family-preserving shuffle | feasible via Manning groups; weak within TK (17/21 parents are TK) - reported as caveat |
| assay-gain weights | feasible from WT rows of train parents only |
| scaffold-cold ligand strata | feasible: Bemis-Murcko scaffolds on 183 ligands |
| replicates for noise floor | NOT available - noise separated only via structure assumptions (flagged) |

### 5.2 Identifiability and shortcut risks (ranked)

1. Shared-pattern entanglement (highest): any arm can reach R2~0.13 by
   keying the panel-shared ligand pattern. Mitigation: CFOIE residual
   objective + explicit ligand-prior baseline in every table.
2. Parent memorization: pair-level split allows parent-keyed prediction
   (ceiling +0.16 over shared pattern). Mitigation: CIIP-1B as the primary
   deployment claim; pair-level results labeled as such.
3. Assay-link compression: 54% of ligands near WT ceiling; unweighted fits
   average over regimes where interaction is unobservable. Mitigation:
   assay-gain weights + saturation strata reporting.
4. Small-sample overfitting: 49 pairs, 9k cells, effective independent
   units = parents. Mitigation: capacity caps, frozen features, multi-seed
   multi-split, hierarchical bootstrap by parent.
5. Selection leakage: any tuning on test labels. Mitigation: grouped
   checkpoint selection on val only; frozen metrics code; prereg before
   training.
6. Attention over-interpretation: routing weights are statistics; never
   contact probabilities (hard reporting rule).

### 5.3 Assay-gain weighting (new, deployable-compatible)

Per-ligand sensitivity from the single-dose Hill link: with
y = 100 C/(C+IC50), |dy/d log10 IC50| is maximized at y=50 and vanishes at
the ceiling/floor. Empirical weight (train parents' WT rows only):
    w_L proportional to ybar_W(L) (100 - ybar_W(L)), normalized.
This downweights the 99 ceiling-dominated ligands and upweights the ~84
mid-zone ligands where interaction is observable. Weights are frozen before
training; stratified (unweighted) metrics are always reported alongside.

### 5.4 Negative-control matrix

| control | destroys what | expected if mechanism real | verdict if violated |
|---|---|---|---|
| train-time ligand-label permutation (within pair) | ligand-conditioning of c | R2 -> ~0 | any surviving signal is assay artifact -> STOP |
| eval-time permutation distribution | significance of ordering | observed R2 > 95th pct of null | treat as null result |
| random-protein swap (pair partner) | protein-conditioning | R2 -> 0 | entity-keyed artifact |
| mutation-erasure (X-out both constructs) | mutation information | deployable model degrades toward WT-prior | if unchanged, mutation not used |
| same-parent wrong-mutation assignment | mutation specificity | correct assignment beats wrong | only parent-context learned |
| family-preserving profile prior | trivial transfer | model beats -0.02 baseline | sequence adds nothing |
| ligand-prior baseline (train mean profile) | shared pattern | model beats 0.13 (CIIP-1A) / prior (CIIP-1B) | no protein-conditioned value |
| free-pairwise ceiling | integrability cost | potential within ~1 SE of free head | potential class too restrictive |
| teacher wrong-site privilege | privileged validity | teacher degrades | teacher shortcut |

### 5.5 Metrics contract (all reported, no NaN laundering)

centered MSE; raw-endpoint MSE (auxiliary only, never mixed); R2; explained
variance; Spearman; Pearson; dead-zone sign accuracy (|c|>10); slope and
scale-recovery (regression of truth on prediction; diagnostic only);
variance recovery var(pred)/var(true); nonconstant coverage; rank-evaluable
denominator; per-parent table; per-mutation-class (gatekeeper /
activation-loop / juxtamembrane / other, mapped by canonical position);
scaffold novelty strata; assay saturation strata (WT-mean bands);
calibration/uncertainty (C5 only); parent-level bootstrap CIs;
correct-vs-control paired bootstrap; multi-split stability. Any NaN rank
metric stays undefined.

---

## 6. Recommended mainline (unique)

**OLR-Potential: Orthogonalized Ligand-Routed interaction potential** =
C1 (LCRR deployable field) + C2 (orthogonal decomposition layer) + C3
(cross-fitted residual objective) with assay-gain weights; C4 (site-teacher
distillation) as a pre-registered optional ablation; C5 (SCNP) reserved for
the Phase-6 few-shot bridge.

Why this one and only this one:
- Trainable: <2 M params on frozen features, 9k cells, minutes/seed.
- Deployable: full sequence + ligand only; no coordinates, no test-time
  adaptation, no closed-form solvers anywhere.
- Unified: one scalar s_theta(P,L); every contrast is a finite difference
  (same object serves CIIP, BindingDB bridge later, few-shot prior).
- Falsifiable: each component is a single toggle (router vs mean-pool; OID
  on/off; residual vs raw target; gain weights on/off; teacher on/off), so
  any claim reduces to a univariate experiment.
- Addresses the audited failure directly: the entangled bilinear term is
  replaced by an explicit factorization (shared pattern via nuisance +
  protein-conditioned deviation via ligand-conditioned routing), and the
  assay-link effect is modeled instead of averaged.
- Connects to few-shot DTA: the same s_theta is the k=0 prior; C5 adds
  support conditioning without gradients.

Expected signature if the biology cooperates: parent-disjoint residual R2
in (0, 0.29] with controls destroyed and per-parent stability. Expected
signature if it does not: all deployable arms collapse to the ligand prior
under parent-disjoint evaluation - itself the answer to the deployment
question at this power.

### Graded claim ladder (what each phase can and cannot claim)

R1 representation: deployable features sense mutations. (Already SUPPORTED
for ESM residue states; extend: surprise channel.)
R2 identification: mutation information changes ligand-conditioned
prediction beyond controls (pair-level scope first).
R3 deployment: parent-disjoint transfer above ligand prior with controls
destroyed (THE deployable claim).
R4 bridge: few-shot DTA gains via s_theta prior (Phase 6, separate line).
R5 mechanism: only R3 + cross-endpoint replication would support
interaction-mechanism language; never asserted from this panel alone.

---

## 7. Minimal experiment sequence and stop rules

Phase 3 (done at report issue): successor preregistration + structural
tests (no training).
Phase 4: single-seed smoke, CIIP-1A scope (pair-level), full control set;
gate: pipeline runs end-to-end, metrics finite, ligand_only arm implements
the true profile prior (fixing the audited design gap), permutation control
executes.
Phase 5: 5 seeds x 3 frozen splits, CIIP-1A scope, then CIIP-1B
parent-disjoint (primary), matched controls, parent bootstrap.
Phase 6 (conditional on R3 support): BindingDB Native Ki bridge /
Tanimoto incremental correction / SCNP few-shot.

STOP (any one triggers formal termination of the mechanism line):
1. all deployable representations indistinguishable between
   mutation-erased and correct inputs;
2. ligand permutation does not destroy performance (assay artifact);
3. same-parent wrong-mutation indistinguishable from correct;
4. deployable router requires mutation coordinates to succeed;
5. parent-disjoint split cannot be formed (data side);
6. all positive effect driven by a single parent (leave-one-parent-out
   sign flip);
7. endpoint effective signal insufficient (shared-pattern residual variance
   below noise scale estimated in simulation);
8. candidate works only on oracle subsets;
9. candidate needs closed-form solves or test-time optimization;
10. candidate cannot unify under one s_theta(P,L).

Final verdicts are issued per claim ladder level as SUPPORTED /
NOT_SUPPORTED / UNRESOLVED, with UNRESOLVED explicitly meaning "power or
data insufficient at this panel scale", never silently upgraded.

---

## 8. Immediate next actions

1. Freeze the Phase-3 preregistration for OLR-Potential (issued alongside
   this report in tools/research/stageCIIP2_olr_potential_20260820/).
2. Structural tests (no training) for: centering correctness, contrast
   antisymmetry/cycle-zero, no-coordinate input contract, permutation
   wiring, leakage guards (train-only gain weights, fold-fitted nuisance).
3. Phase-4 single-seed smoke upon green tests.

---
---

## 9. Phase 3-4 execution log and methodological findings (2026-08-20)

The successor stage (tools/research/stageCIIP2_olr_potential_20260820) was
implemented and structurally qualified before any real-data training.
Findings that materially changed the protocol (each frozen in a dated
addendum BEFORE real-data runs):

1. AM-1 (erasure degeneracy): for ANY potential-form model, full
   mutation-erasure is degenerate - erased WT and erased variant states are
   identical, so the contrast is identically zero. The informative control
   is site-erasure of the variant construct (context preserved, site
   residue masked), implemented as C-erased-site.
2. AM-2 (OID identity theorem): adding any separable main-effect component
   to s leaves the panel-centered contrast EXACTLY unchanged (proposition
   asserted numerically in test T3). Orthogonal decomposition therefore
   cannot change contrast predictions - it is an invariance, not a modeling
   choice. A3-oid is retained only as an equivalence/determinism check.
   This is a useful negative result: several "orthogonal interaction"
   designs in the literature constrain contrasts that are already invariant.
3. AM-3 (panel width): per-pair ligand panels have width 179-183; profile
   priors must be built on a padded panel before averaging.
4. AM-4 (router skip): pure attention routing optimizes too slowly at this
   scale; alpha = mean-pool skip + routed deviation makes the router a
   strict superset of the bilinear arm (the ladder stays falsifiable).
5. AM-5 (budget): lr 1e-3/300 epochs undertrains severely (train MSE still
   falling at cap); lr 3e-3 + cosine, 900 epochs converges on the synthetic
   instrument (train MSE 384 -> 12).
6. ADD-2 (evaluation definition): scoring on cross-fitted residual targets
   makes the analytic ligand-prior ill-posed (it predicts exactly the
   removed component; observed R2 -3.0), inflating any Delta-R2. Primary
   evaluation and checkpoint selection moved to the FULL centered target
   (the deployable quantity); residual R2 demoted to secondary diagnostic.
   Deployable predictions for nuisance-trained arms compose the ligand-only
   nuisance with the model contrast.
7. Instrument redesign: the first planted field keyed on raw mean ESM
   states, but kinase mean states are so similar cross-parent that the
   prior absorbed nearly everything (between-parent spread ~16 vs the
   measured 134.8). The plant now keys on parent-DEVIATION fields scaled to
   the measured between-parent variance, keeping it linear in sequence
   features (learnable in principle).

These corrections were all discovered on synthetic or structural grounds;
no real-data result existed when they were frozen. The instrument
qualification verdict under the final protocol is recorded in
QUALIFICATION.json and gates Phase 4/5 interpretation.

### 9.1 Final instrument verdict (frozen consequence)

INSTRUMENT_UNDERPOWERED, robust across estimator family (A1 bilinear
0.105, A2 router -0.327, A5 composed 0.560-0.566 total R2; only +0.0156..
+0.0296 above the analytic ligand-pattern prior 0.5361), rank {2,8}, and
selection rule (val early-stop optimal; last-epoch checkpoint -0.379).
A planted parent-conditioned field linear in real sequence features, scaled
to the measured between-parent variance (134.8%^2 = half the interaction
variance), is learnable on train (train R2 0.67) but transfers at ~5% of
its recoverable magnitude at n=39 train pairs / 13 val pairs.

Scientific reading: the binding constraint is not representation, not
architecture, and not optimization-on-train - it is the number of
independent mutation conditions. The panel's 49 pairs over 18 parents
cannot support honest selection of a parent-conditioned interaction
function at the pre-registered 0.25 Delta-R2 standard. Every real-data
deployment claim from this panel is therefore UNRESOLVED (power), and the
field needs either multi-panel pooled same-endpoint data or a fundamentally
larger mutation-conditioned panel before the deployability question can be
adjudicated. This is the single most important output of the CIIP-2 cycle.

### 9.2 Real-data Phase 4 outcome (single seed, split S1)

Smoke gates: pipeline PASS, A0-sanity PASS (0.1313), nonconstant 9/9 PASS,
permutation-destroyed FAIL (C-perm 0.1818 >= best correct arm A4 0.1732,
A5 0.1388). No deployable arm separates from the permutation control; the
largest model increment over the ligand-pattern prior is +0.042 R2 (A4),
noise-level at this power. Phase 5 (multi-seed parent-disjoint) is NOT
authorized by the frozen chain and could only add noise-level measurements.

Programme verdict (final, this cycle):
- R1 representation: SUPPORTED. R2 identification: NOT SUPPORTED at this
  power. R3 deployment: UNRESOLVED (power; instrument-planted transferable
  field unrecoverable). R4 few-shot DTA bridge: BLOCKED. R5 binding
  interpretation: NOT CLAIMED.
- Overall: UNRESOLVED (power) - the panel cannot adjudicate deployable
  protein-conditioned interaction learning in either direction at the
  pre-registered standard. Successor requires >=100 mutation conditions
  across >=30 parents on one endpoint, or a Ki/Kd DeltaDeltaG corpus.



