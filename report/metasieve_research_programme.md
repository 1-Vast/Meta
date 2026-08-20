# Cold-Target DTA: Programme Review, Literature Assessment, and Staged Research Design

**Scope.** Independent review of the stated objective, the current protocol, and the synthetic qualification programme; literature verification against primary sources; and a finite staged research design. Theory and design only — no repository inspection, no code, no implementation patches.

**Epistemic marking.** Statements marked **[V]** are verified against a primary or official source with a URL. Statements marked **[I]** are my interpretation or inference. Statements marked **[U]** are unverified claims taken from the brief that I could not check and that should be treated as provisional.

---

## 1. Executive judgment: is the present objective overconstrained?

**Yes — but not because any individual constraint is unreasonable. The protocol is overconstrained because it asks one experimental surface to settle four logically distinct questions at once, two of which are in direct tension.**

The four questions currently entangled:

| # | Question | Nature |
|---|---|---|
| Q-perf | Can we reach MSE ≤ 1.00 pK² at k = 0…5 on cold targets? | Performance |
| Q-mech | Is there a *transferable protein-conditioned interaction* signal at all? | Mechanism / attribution |
| Q-repr | Which protein representation carries it? | Representation |
| Q-train | Which training mechanism extracts it? | Optimization |

The structural error is this. **Q-perf is a variance-explained question; Q-mech is a residual-after-main-effects question. They are measured on different quantities and cannot share a surface.**

Concretely: the brief simultaneously requires (a) strict removal of protein and ligand main effects, and (b) MSE ≤ 1.00 pK². But most of the explainable variance in heterogeneous DTA labels *is* main effects — target-level potency level, assay level, ligand promiscuity. Your own evidence says so: at k = 0 much of the error is target/assay calibration error **[U]**. Once you subtract the main effects, the residual interaction variance is small, and "MSE ≤ 1.00 pK² on the residual" is a far more severe requirement than the headline number implies. Conversely, if the main effects are retained, MSE ≤ 1.00 pK² becomes reachable — but every unit of the improvement is attributable to calibration, which you have correctly ruled inadmissible as a mechanism claim.

**Therefore: the performance target and the attribution requirement must never be evaluated on the same surface.** This is the single most important structural correction in this report. Everything downstream follows from it.

A second overconstraint is subtler and is a genuine protocol defect rather than a philosophical one:

> **CD-HIT40 on full-length sequence does not produce pocket-cold splits.** For kinases in particular, the ATP site is strongly conserved across the family; two kinases below 40% *global* identity can have near-identical KLIFS 85-residue pockets. KLIFS exists precisely because a consistent 85-residue binding-site alignment spans the whole family **[V]** (https://pubs.acs.org/doi/10.1021/jm400378w, https://academic.oup.com/nar/article/44/D1/D365/2502606).

So the current split is simultaneously *too strict* (it destroys the homology information a real project would have) and *not strict enough* (it does not guarantee pocket novelty, which is what actually matters for an interaction claim). The split unit should be **pocket similarity over the aligned binding-site positions**, not global sequence identity. **[I]**

**Classification of the current protocol**, against the four options offered in the brief: it is **(4) a mixture of several distinct deployment scenarios that should be separated**, with the double-cold arm being **(3) an overconstrained mechanism-identification problem** when it is asked to also carry a performance target.

Maximal strictness is not maximal scientific value. Strictness buys protection against a specific confounder. Stacking every strictness simultaneously means each confounder is defended against, but the residual effect becomes too small to detect with the number of independent clusters available, and a null result becomes uninterpretable — you cannot distinguish "no transferable interaction exists" from "the design had no power." A test that cannot distinguish those two outcomes is not a strong test; it is an expensive one.

---

## 2. The current protocol versus published cold-target / few-shot DTA protocols

The most important finding from the literature survey is not about any single method. It is this:

> **The three leading few-shot bioactivity models — MetaDTA, FS-CAP, and ActFound — are all protein-blind. None of them uses the target sequence at all. They define the state of the art in few-shot bioactivity prediction without any protein-conditioned mechanism.**

- **MetaDTA** (ICLR 2022 workshop; Attentive Neural Processes) states explicitly that its formulation *does not directly use the target protein information*; the target-specific affinity function is inferred from the support pairs alone **[V]** (https://openreview.net/pdf?id=yzlif16IASM, https://iclr.cc/virtual/2022/8615).
- **FS-CAP** is explicitly *target-free / ligand-based*. Context compounds and their activities are encoded via a multiplicative featurization; a separate encoder handles the query; there is no protein input. Trained on BindingDB with test targets held out of training **[V]** (https://pubs.acs.org/doi/10.1021/acs.jcim.4c00485, https://arxiv.org/abs/2311.16328).
- **ActFound** (Nature Machine Intelligence 2024) uses Siamese pairwise *within-assay* relative learning plus MAML-style meta-learning across ~35k–50k assays from ChEMBL and BindingDB. Its stated limitation is that it "cannot consider the meta-data for each assay, such as the sequence of the protein target" **[V]** (https://www.nature.com/articles/s42256-024-00876-w, https://www.biorxiv.org/content/10.1101/2023.10.30.564861v1.full.pdf).

Two consequences follow.

**(a) Your ligand-only control arm is not a weak baseline. It is the published state of the art.** Any protein-conditioned method must beat a *well-tuned* FS-CAP/ActFound-class ligand-only model, not a naive Tanimoto-weighted average. If your protein-conditioned model is compared against a weak ligand-only arm, a positive margin is uninformative.

**(b) There is a genuine open niche.** Nobody has demonstrated that protein conditioning adds transferable signal on top of a strong support-conditioned ligand model. That is precisely your Core Task 1. The field's silence here is evidence that it is hard, not that it is unexplored trivia.

Other reference points:

| Work | Task | Protein given? | Support/query same target? | Support size | Test-time gradients? | Protein split | Ligand split | Cold surface | Assay identity | External pretraining | Protein attribution tested? | Explicable by ligand similarity or calibration? | Truly zero-shot? |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **AdaMBind** (Nat Commun 2026) https://www.nature.com/articles/s41467-026-70554-5 | DTA regression | Yes (GNN drug + protein encoder) | Yes (task = one protein) | few-shot | Yes (MAML inner loop) | "random" and "novel" task split **[U]** | not scaffold-controlled **[V]** as reported | protein-cold | no | no | no | **plausibly yes** | No — requires support |
| **MetaDTA** https://openreview.net/pdf?id=yzlif16IASM | DTA regression | **No** | Yes | few-shot | No (ANP) | target-held-out | none | protein-cold | no | no | n/a | by construction ligand-only | No |
| **FS-CAP** https://pubs.acs.org/doi/10.1021/acs.jcim.4c00485 | Continuous activity | **No** | Yes (same assay) | 1–16 | No | assay held out | none | assay-cold | implicit | no | n/a | by construction ligand-only | No |
| **ActFound** https://www.nature.com/articles/s42256-024-00876-w | Bioactivity regression | **No** | Yes (same assay) | ~16 typical | Yes (MAML + kNN-assay FT) | assay held out; ChEMBL↔BindingDB overlap filtered by correlation **[V]** | scaffold generalization reported | assay-cold, scaffold | yes (assay is the task) | no | n/a | by construction ligand-only | No |
| **ZeroBind** (Nat Commun 2023) https://www.nature.com/articles/s41467-023-43597-1 | **Binary DTI**, AUROC/AUPRC | Yes (AlphaFold protein graph + subgraph IB pocket) | task = protein | 0 / few | Yes | unseen proteins | unseen drugs | both | no | yes (structures) | pocket subgraphs shown, not control-tested | not excluded | Yes, for classification |
| **PSICHIC** (Nat Mach Intell 2024) https://www.nature.com/articles/s42256-024-00847-1 | Affinity + functional effect | Yes, residue tokens, sequence-only | n/a | n/a | No | reported splits **[U]** | — | — | no | no | interpretable fingerprints, not control-tested | not excluded | n/a |
| **FS-Mol** (NeurIPS D&B 2021) https://github.com/microsoft/FS-Mol | Binary activity (task = protein) | **No** | Yes | 8–128 | varies | disjoint task sets | random/scaffold | assay-cold | yes | varies | n/a | ligand-only by construction | No |
| **PCM** (van Westen, Cortés-Ciriano et al.) https://pubs.rsc.org/en/content/articlelanding/2011/md/c0md00165a | Affinity regression | Yes (target descriptors) | n/a | n/a | No | typically **leave-one-target-out within a family** | usually random | protein-cold, within-family | varies | no | rarely | often yes | Yes, within-family |

Three protocol-level observations:

1. **PCM's classical validation is leave-one-target-out within a target family, not CD-HIT40 cross-family cold.** PCM's own literature is explicit that full extrapolation to novel compounds *on* novel targets is not achieved **[V]** (https://academic.oup.com/ib/article/6/11/1023/5199186). Your protocol is substantially harder than the entire PCM tradition's standard evaluation. That is defensible as a stress test; it is not a fair comparison surface.

2. **The DTI/DTA field has a documented shortcut problem, and the community's response has been split redesign, not attribution controls.** Hidden ligand bias in DUD-E and the Human dataset is well established (TransformerCPI's label-reversal experiments; DrugBAN's cold-pair split motivation) **[V]** (https://arxiv.org/pdf/2208.02194, https://academic.oup.com/bib/article/26/5/bbaf491/8260789). An OOD assessment paper shows that association by drug alone or target alone carries predictive power and therefore constitutes leakage under interaction-level splitting **[V]** (https://www.biorxiv.org/content/10.1101/2022.04.20.488898.full.pdf). **Your control matrix is more rigorous than anything in the published DTA literature.** That is a real asset and should be preserved — but it is an *attribution* instrument, and attribution instruments belong on an attribution surface, not on the performance leaderboard.

3. **Simple baselines are stubbornly strong, which raises the bar for any claimed innovation.** On FS-Mol, Random Forest on fingerprints outperforms single-task GNNs and the Molecule Attention Transformer **[V]** (https://arxiv.org/pdf/2305.09481). On activity cliffs, descriptor-based ML beat deep learning across 30 targets in the MoleculeACE benchmark **[V]** (https://pubs.acs.org/doi/10.1021/acs.jcim.2c01073). Any protein-conditioned architecture must clear these, not just clear other deep models.

---

## 3. Decomposition of deployment scenarios

The brief's four scenarios are well drawn. My assignment of roles:

| Scenario | Description | Role | Rationale |
|---|---|---|---|
| **A — cold-target lead optimization** | unseen target; 5–40 support compounds, possibly same series; sequence/family/predicted structure available; need within-target SAR and ranking | **Primary performance benchmark** | This is the modal real deployment. It is where a genuine advance has commercial and scientific value, and where the data density actually exists. Metric emphasis: per-target Spearman, activity-cliff sign, centered error. |
| **B — cold-target virtual screening** | unseen target; k ≈ 0–2; large library; need enrichment and uncertainty | **Secondary performance benchmark** | Metric emphasis: enrichment at top 1%/5%, calibration, abstention quality. MSE is close to the wrong metric here. |
| **C — strict double-cold generalization** | unseen protein component *and* unseen scaffold; no support; main effects removed | **Attribution benchmark only. No performance target.** | This is the mechanism instrument. Its output is a binary admission decision plus an effect size, not an MSE. |
| **D — mutation-specific selectivity** | WT/mutant constructs; same ligands across variants; dense panels | **Biological positive control and the strongest available mechanism surface** | Same-ligand WT vs mutant is the cleanest possible protein-conditioned contrast: the ligand is held exactly fixed, so any signal *must* be protein-conditioned. Promote this from "validation" to "primary mechanism evidence." |

**The most consequential recommendation in this section: Scenario D should be the mechanism surface, not Scenario C.**

Reasoning **[I]**. In Scenario C, the correct-protein arm must beat controls on a residual whose size is unknown, using a protein representation whose relevance is unknown, on a covariate shift that is maximal. Three unknowns, one measurement. In Scenario D, the ligand is held identical, the two proteins differ by one residue, and the estimand — signed change on the same ligand — is *definitionally* not obtainable from ligand similarity, target-level calibration, or ligand-only retrieval, because those are all constant across the pair. **Scenario D isolates the claim architecturally rather than statistically.** That is a much stronger instrument and it uses data that exists.

---

## 4. Revised definition of Core Task 1

Current CT1 asks one question that bundles existence, representation, and regime. Split it:

**CT1-a — Pipeline capability (synthetic, near-closed).**
*Claim:* the gradient-trained learner and its full pipeline can recover a known, identifiable, feature-conditioned interaction on protein-cold, ligand-cold, and double-cold surfaces, under the same optimizer, batching, and checkpoint policy used in production.
*Status:* the train-only oracle already recovers the planted truth on all three cold surfaces **[U]** — so recoverability is established and this is now an *optimizer/parameterization* question, not a science question. One bounded diagnostic remains (Section 5).
*What it proves:* implementation capability only. Never biology.

**CT1-b — Protein-conditioned contrast on real data with support permitted (k ≥ 1).**
*Claim:* on real, legally matched data, the correct protein representation beats every matched protein control at fixed support size.
*Surface:* Scenario D (WT/mutant same-ligand pairs), primary; Scenario C at k ≥ 1, secondary.
*What it proves:* that local protein information carries ligand-specific transferable signal. It does **not** prove zero-shot capability.

**CT1-c — Zero-shot protein-conditioned interaction (k = 0).**
*Claim:* the same, with no support labels.
*What it proves:* the strongest version. Required only for Scenario B/C claims.

**Ordering — and here I disagree with the brief.** Rule G currently says: do not proceed to few-shot adaptation until zero-shot protein-conditioned interaction has passed its own attribution gate. **I recommend inverting this.** Run CT1-b before CT1-c.

Justification **[I]**:
- Attribution validity comes from the *control matrix*, not from k = 0. The shuffled / family-preserving-shuffled / matched-wrong / residue-permuted controls are equally valid at k = 3 as at k = 0, because they hold everything constant except the protein identity fed to the model.
- CT1-b has strictly more statistical power. The support set absorbs the target-level offset, which is the dominant nuisance term; the protein representation then only has to explain the *shape* residual. Removing the largest nuisance raises the signal-to-noise of the very quantity you are trying to detect.
- CT1-b is a valid falsifier of CT1-c. If protein information cannot beat matched-wrong protein *even with support*, it will not at k = 0. So CT1-b is a cheap, high-power, one-directional screen.
- Conversely, failing CT1-c first tells you almost nothing, because it is confounded with the calibration problem you already know dominates at k = 0.

**What failure of CT1-c should and should not block:**

| Failure | Blocks |
|---|---|
| CT1-a fails (learner cannot recover a known identifiable interaction) | Closes the *current learner family as parameterized*. Does not block representation work with a different learner. |
| CT1-b fails (correct protein indistinguishable from matched-wrong on real WT/mutant pairs) | Rejects the *proposed protein representation family*. Does not block ligand-only few-shot performance work. |
| CT1-c fails but CT1-b passes | Blocks **zero-shot protein-mechanism claims only**. Few-shot protein-conditioned claims remain live. |
| CT1-b and CT1-c both fail across representations | Blocks all protein-interaction-innovation claims. Performance track continues as an explicitly ligand-conditioned few-shot contribution. |

**Nothing in this list blocks all few-shot modelling.** Given that the published SOTA is protein-blind, holding practical few-shot work hostage to a protein-mechanism gate risks spending the whole programme on a mechanism the field has not yet demonstrated anyone can extract, while the deliverable that would actually be publishable and useful sits unbuilt.

---

## 5. Finite stopping rule for synthetic qualification

### 5.1 Diagnosis of what the synthetic programme has already established

Taking the brief's account at face value **[U]**: a train-only oracle recovers the planted, span-restricted truth on protein-cold, ligand-cold, and double-cold surfaces, while the gradient learner fits training interactions near-perfectly and fails on cold proteins, with 33–43% of learned protein-map energy in directions unconstrained by the training protein features.

**That combination is diagnostic, and it does not indicate a scientific problem. It indicates an optimizer/parameterization pathology.** The reasoning **[I]**:

- Identifiability is established (the oracle recovers it), so it is not non-identifiability.
- Train loss is near-zero, so it is not model misspecification or insufficient capacity.
- Cold error is large while the identifiable component is recoverable, so it is not insufficient data.
- The residual explanation is that the learner's parameters carry energy in directions the training data cannot constrain, and the optimizer put it there.

A quantitative check that should be run before anything else: with a nominal 32-dimensional protein feature space at train rank 28, an isotropic random initialization deposits roughly (32−28)/32 ≈ 12.5% of its energy in the exact null space. Observing 33–43% is well above that. So either the "identifiable span" is being thresholded to include weakly-excited directions, or null-space energy is *accumulating during training*.

**The single measurement that discriminates these:** log ‖P_null · A‖²_F / ‖A‖²_F as a function of training step, and separately report energy by singular-value decile of the train protein feature matrix rather than as a binary in/out-of-span split. A flat trajectory implicates initialization; a rising trajectory implicates the optimizer or the parameterization.

Two mechanisms produce a rising trajectory, and neither is addressed by span-initialization plus mild L2:

1. **Adaptive optimizers do not preserve non-axis-aligned subspaces.** For an unfactored bilinear map, the gradient with respect to the protein index lies in the row space of the training protein features by construction, so plain SGD from a zero initialization can never leave that span. Adam's diagonal preconditioner rescales coordinates and therefore maps an in-span gradient to a generally out-of-span update. This leaks on every step, with zero train-loss penalty.
2. **Factorization leaks.** If the map is a product `A = L·R` with the protein index on `R`, the update contribution `ΔL·R` has protein-side row space equal to that of `R`, not of the training features. Any null-space content in `R` is amplified rather than starved.

**Consequence for the frozen successor:** span-initialization plus mild L2 targets only the initialization hypothesis. If mechanism 1 or 2 is operative, that experiment will fail *even though the learner family is fine*. **A failure of the span-init + L2 experiment must not, on its own, close the bilinear family.** That would be a Type-II error against the leading hypothesis.

### 5.2 Answers to the eight questions in RQ3

**1. Minimum positive control.** One experiment: a truth of the learner's own functional form, restricted to the identifiable training span, recovered by the *production* training path — same optimizer, same batching, same checkpoint policy, same seed discipline — evaluated on all three cold surfaces against the oracle ceiling and the shuffled-protein floor. Nothing more is required to establish pipeline capability.

**2. Should the truth be matched to the learner's functional form?** **Yes, for this purpose.** The synthetic instrument's job is to test the *pipeline*, holding model class correct. A mismatched truth confounds optimization failure with misspecification, and you then cannot attribute a failure. Testing misspecification is a legitimate but separate question, and it should be asked on real data, where the true functional form is the thing you actually care about. This is also the guard against the self-fulfilling-benchmark risk: because the truth is matched *by design and declared as such*, a PASS is explicitly and only a statement about implementation capability, and cannot be mistaken for evidence about biology.

**3. Defensible number of further synthetic successors: one.** And it should be governed by a hard admission rule: *a new synthetic stage may be created only if it tests a pre-registered discriminating measurement whose outcome changes the terminal decision.* Redesigning the planted truth is not a discriminating measurement. Changing the parameterization or the optimizer to test a named mechanism is.

**4. Is identifiable-subspace parameterization (A = V_train·G) a valid diagnostic?** **Valid as a diagnostic; invalid as a deliverable.** It is legitimate under your constraints — V_train derives from training *inputs* only, uses no labels and no truth weights, and G remains gradient-trained. It answers exactly one question: *is the bilinear learner family capable in principle, with the identifiable constraint made explicit?* It does make the task easier, and that is the point: it converts an ambiguous failure into an attributable one.

But it does **not** transfer. Real protein feature matrices — KLIFS one-hot over 85 positions, or local ESM windows — are effectively full-rank with a decaying spectrum. There is no exact null space to project out; the pathology becomes ill-conditioning rather than unidentifiability. **The durable fix is spectral regularization or a trust region on the protein map, not a hard span projection.** Report the diagnostic as such and carry forward the conditioning control, not the projection.

**5. Should the model pass every endpoint / missingness / censoring level?** **No.** Those levels test *robustness*, which is only meaningful conditional on *existence*. Establish existence at the cleanest level; test robustness on real data, where the censoring and missingness structures are the actual ones you face rather than a synthetic approximation of them. Requiring a full robustness ladder before existence is established is the principal mechanism by which synthetic programmes become unbounded.

**6. Threshold calibration.** Replace heuristic thresholds (0.30 Spearman, 0.70 sign accuracy, +0.05, +0.03) with power-calibrated ones. Procedure:
- Simulate the **null** by running the shuffled-protein arm ≥ 200 times with different seeds. The 95th percentile of the null distribution of each statistic is the floor.
- Run the **oracle** to obtain the achievable ceiling on the same surface.
- Pre-register the gate as a fraction of the oracle–null gap (e.g. the learner must recover ≥ 50% of the gap), and pre-register the **minimum detectable effect** given the number of independent clusters.
- Report the power of the design explicitly. If the design cannot detect the effect size you would consider scientifically meaningful, a null result is uninterpretable and the design must be fixed before it is run.

This converts arbitrary numbers into statements with operating characteristics, and it makes "we failed" and "we had no power" distinguishable — which the current thresholds do not.

**7. Failures that close the learner family.** The bilinear family is closed if and only if: with a *matched-form* truth, verified identifiable, span-restricted parameterization made explicit, plain-SGD-from-zero and Adam arms both run, and checkpoint selection verified not to monitor a leaking surface, the gradient learner still fails to recover cold-surface performance meaningfully above the null floor. Any failure that has not eliminated the optimizer and parameterization mechanisms above does not close the family.

**8. When to move to real biology.** **Now, in parallel — not after.** The synthetic programme is *incapable in principle* of producing evidence about biology; it can only falsify the pipeline. Serializing biology behind it converts a bounded engineering question into the critical path for the entire programme. Running the real matched-pair control benchmark in parallel is not a governance violation; it is the only way to bound total elapsed time, and the real benchmark is small and cheap (Section 16).

### 5.3 The stopping rule, stated

> **Synthetic qualification terminates after at most one further stage.** That stage is a named diagnostic successor with a frozen preregistration containing: (i) the null-space-energy trajectory measurement; (ii) arms for plain-SGD-from-zero, Adam, and the explicit span reparameterization; (iii) power-calibrated thresholds derived from a ≥200-run shuffled-protein null and an oracle ceiling. At its conclusion a terminal PASS/FAIL on the bilinear learner family is recorded. **No synthetic stage may be created whose purpose is to redesign the planted truth.** If the diagnostic passes, synthetic work is closed and does not resume. If it fails, the family is closed and the programme moves to a different learner family — also without resuming synthetic redesign.

---

## 6. Ranked protein–ligand representation directions

Assessment criteria as requested. Ranked by expected decisive information per unit cost, given that the mechanism claim is the bottleneck.

**Rank 1 — Aligned-position *difference* representation for the WT/mutant estimand.**
Feed only the KLIFS-aligned per-position difference between the two proteins, and predict a signed ΔΔ.
- Biological plausibility: high. Data requirement: low (panels exist). Pose-free: yes. **Shortcut resistance: maximal — identically zero for identical proteins, so it cannot encode protein identity at all.** Cold-parent transfer: directly testable by splitting on parent. Zero-shot: yes for the ΔΔ estimand. Few-shot: yes. Interpretability: high (per-position attributions are meaningful). Cost: very low. Risk of repeating past failures: lowest of any candidate, because target-ID memorization is architecturally excluded rather than statistically controlled.
- This is the highest-value representation direction and it is currently not in your ranked list at all.

**Rank 2 — KLIFS-aligned residue tokens + z-scale physicochemical features, with fragment-level ligand tokens and a low-rank position×fragment interaction operator.**
- Plausibility high; the 85-position alignment is a curated, family-wide standard **[V]** (https://academic.oup.com/nar/article/44/D1/D365/2502606). Pose-free. Interpretable. Cheap. Constrained to kinases — which is acceptable: *make the mechanism claim where the alignment exists, then test transfer.*
- Shortcut risk: moderate. 85 one-hot positions can still identify a target. Mitigate with a low-rank bottleneck on the protein path and mandatory residue-permutation controls.

**Rank 3 — Local ESM at aligned positions (never pooled) + explicit physicochemical features.**
- Adds evolutionary/contextual information that one-hot lacks. But note the sharpest caveat in this report: **PLM likelihood-based variant effect prediction is ligand-agnostic by construction.** ProteinGym-style zero-shot scores are log-likelihood ratios over sequence alone **[V]** (https://arxiv.org/pdf/2507.02624, https://arxiv.org/abs/1909.09157 is unrelated — correct reference: https://proteingym.org and Meier et al. NeurIPS 2021). A model whose protein pathway is dominated by PLM likelihood will produce *ligand-invariant* mutation effects — which is exactly the protein-main-effect shortcut you have already ruled inadmissible.
- Shortcut risk: **high.** ESM embeddings encode identity strongly. Requires the full control matrix.

**Rank 4 — Global-plus-local decomposition with an explicitly gated global path.**
- Rather than removing main effects, *model them in a separate, deliberately capacity-limited path that receives gradient only from the absolute-level loss*, and route the interaction claim through the local path. This preserves the performance benefit of calibration while making its contribution measurable and separable. Strongly recommended as an architectural principle regardless of which local representation wins.

**Rank 5 — Interaction-only interface featurization (CORDIAL-style).**
- The PNAS 2025 CORDIAL work argues that generalization failure to novel targets stems from bias toward structure-specific correlations, and addresses it by restricting the model to interface physicochemistry **[V]** (https://www.pnas.org/doi/10.1073/pnas.2508998122). Conceptually the most shortcut-resistant family. Cost: requires a pose or a credible pocket proxy, which the brief says you do not have. Defer, but note as the direction to take if pose becomes available.

**Rank 6 — Sequence-only residue-token cross-attention with physicochemical constraints (PSICHIC-class).**
- PSICHIC keeps residue and atom tokens until an interaction fingerprint weighted by intermolecular-force scores, from sequence alone, and reports matching or exceeding structure-based methods on affinity **[V]** (https://www.nature.com/articles/s42256-024-00847-1). This validates the architectural family you want. Open question is entirely whether it survives strict pocket-cold + scaffold-cold with the control matrix — which nobody has tested.

**Rank 7 — MSA / conservation / co-evolution priors.**
- Useful as *position weighting* (which of the 85 positions matter), not as a signal source. Ligand-agnostic. Cheap. Include as a prior, not as a representation.

**Rank 8 — Predicted monomer structure / predicted pockets, no complex pose.**
- Moderate cost, modest expected return without a pose. ZeroBind's pocket-subgraph approach is the reference point **[V]** (https://www.nature.com/articles/s41467-023-43597-1) but it is a *binary classification* method — AUROC/AUPRC, not pK regression. Do not import its performance claims into a regression setting.

**Rank 9 — Retrieval-augmented protein-family context.**
- Likely to *help performance and hurt attribution*. If used, it must be a declared, ablated component, and any gain it produces must be reported as retrieval, not interaction.

**Rank 10 (negative comparator only) — Global mean-pooled ESM.**
- Retain solely as the control that should fail. Your observation that it behaves as a target identifier **[U]** is consistent with the design: mean-pooling over a full-length sequence is nearly invariant to a single substitution.

**Two cross-cutting cautions.**
- *Attention weights are not mechanistic explanations.* Interpretability claims from attention maps must be validated against held-out structural ground truth or not made.
- *Ensembling across conformational/pocket states* is worth trying only after a single-state model has passed attribution; otherwise it adds capacity that can absorb shortcuts.

---

## 7. Ranked training innovations

Constraint recap: the core innovation must live in training; must be gradient-trained and normally optimized; no closed-form adaptation; no deployment-time query labels or query gradients.

**The published landscape constrains what can count as novel:**
- Within-assay/within-target *ligand-pair* relative learning is **taken** by ActFound **[V]**.
- Task = protein with adaptive task sampling weighted by query loss and support/query gradient similarity is **taken** by AdaMBind **[V]** (https://www.researchgate.net/publication/401811999).
- Support-set aggregation into a task encoding is **taken** by FS-CAP and MetaDTA **[V]**.

**Rank 1 — Antisymmetric within-ligand protein-difference supervision (ΔΔ across protein pairs at fixed ligand).**
- This is the open slot. All published relative learning is *within-target across ligands*. Nobody trains on *within-ligand across proteins* with hard structural constraints.
- Impose exchange-antisymmetry f(P₁,P₂,L) = −f(P₂,P₁,L) and identity-zero f(P,P,L) = 0 **architecturally**, not as a penalty. Then the mechanism is unable to express a ligand-only or target-identity solution: both vanish under the constraint.
- This is simultaneously a training innovation and an attribution guarantee, which is unusually efficient.
- Directly served by Scenario D data.

**Rank 2 — Level/shape gradient disentanglement with enforced gradient ownership.**
- Absolute-level loss updates only the calibration path; relative/interaction losses update only the interaction path; measure and report the gradient conflict between them (cosine similarity per step, per parameter group). PCGrad supplies the standard measurement and mitigation apparatus **[V]** (https://arxiv.org/pdf/2001.06782), with the caveat from CAGrad that PCGrad lacks a clear per-step objective and its random ordering matters **[V]** (https://arxiv.org/pdf/2110.14048).
- Use conflict as *instrumentation first*: the fraction of interaction-path gradient explained by the level loss is a direct, quantitative version of the "is this just calibration?" question.

**Rank 3 — Censoring-aware pairwise likelihood.**
- Essential, not optional, for the kinase panels (Section 8). Single-concentration percent-inhibition data is censored at both ends. A Tobit-style or interval likelihood on *differences* is needed so that pairs where both members are saturated contribute correctly (i.e., contribute almost nothing) instead of contributing a spurious zero.

**Rank 4 — Per-dataset heads on a shared trunk with mandatory external-data ablation.**
- Not novel, but it is the correct handling of label-semantic incompatibility and it makes the external-data contribution reportable.

**Rank 5 — Assay-balanced episodic sampling.** Cheap, reduces assay confounding, ablatable.

**Rank 6 — Uncertainty-aware abstention.** High practical value in Scenario B; weak as a mechanism claim.

**Rank 7 — Activity-cliff curriculum.** Plausible given the documented failure of models on cliffs **[V]** (https://pubs.acs.org/doi/10.1021/acs.jcim.2c01073), but it optimizes a ligand-side property and will not support a protein-mechanism claim.

**Demoted — MAML-style inner-loop adaptation as the *core* innovation.** See Section 8.

**Excluded by your own constraints** — ridge/closed-form adaptation, test-time query-label optimization.

---

## 8. Decision on AdaMBind / MAML

**Decision: use as a strong baseline; adopt the adaptive task-sampling component only as an ablatable module; reject the inner loop as the core innovation; reject entirely for k = 0.**

Grounds:

1. **Zero-shot incompatibility.** AdaMBind's task *is* a protein with its ligands and affinities; the meta-test procedure adapts on a support set. With k = 0 there is no support set and the inner loop is undefined. It cannot be the mechanism for the zero-shot arm of your objective.

2. **The inner loop is probably not where the benefit lives.** Raghu et al. showed via layer-freezing and representational analysis that MAML's effectiveness is dominated by *feature reuse*, not rapid learning; ANIL, which removes the inner loop for everything but the head, matches MAML on few-shot classification and RL **[V]** (https://arxiv.org/abs/1909.09157, https://openreview.net/forum?id=rkgMkCEtPB). If that carries over, the correct minimal form of AdaMBind is ANIL-like: a well-trained shared representation plus a light task-specific head. That is much cheaper and much easier to attribute.

3. **The comparison in the source paper does not isolate meta-learning.** AdaMBind's baselines were fine-tuned for 5 steps on the same support sets **[V]** (https://www.researchgate.net/publication/401811999), which controls for *support-label access* but not for fine-tuning budget or for the adaptive sampling component separately. To separate meta-learning from extra support labels, fine-tuning budget, and task sampling, you need a four-way factorial: {pretrain-only, pretrain + equal-budget fine-tune, MAML, ANIL} × {uniform sampling, adaptive sampling}. Without that factorial, any adopted gain is unattributable — which is exactly the failure mode the brief is trying to escape.

4. **Adaptive task sampling is the genuinely portable idea.** Weighting tasks by query loss and support/query gradient consistency is a sensible, cheap, architecture-agnostic component. Take it; ablate it; do not build the identity of the contribution on it, since it is already published.

---

## 9. Dataset roles and label semantics

**Governing rule: never merge incompatible endpoints into one regression target.** Shared trunk, separate heads, separate metrics, separate claims.

| Dataset | Endpoint — what it *actually* measures | Units | Role | Hard cautions |
|---|---|---|---|---|
| **BindingDB Ki** | equilibrium inhibition constant | pKi | **Primary DTA performance surface** | Heterogeneous provenance; overlaps ChEMBL heavily — ActFound had to filter test assays by cross-source correlation **[V]**. Deduplicate across sources before any split. |
| **Davis** | dissociation constant, 72 inhibitors × 442 kinases **[V]** (https://www.nature.com/articles/nrd3647) | pKd | Secondary performance surface | Heavily censored at the assay ceiling; a large mass of values sits at the detection floor. Uncensored MSE on Davis is partly a measure of how well you predict a censoring artifact. Report censoring-aware metrics or restrict to the dynamic range. |
| **KIBA** | composite score integrating Ki/Kd/IC50 | **arbitrary** | Ranking-only surface | **Not an affinity in physical units.** Never report "MSE in pK²" on KIBA. CI/Spearman only. |
| **Saifudeen kinase panel** (Nat Biotechnol 2026, https://doi.org/10.1038/s41587-026-03090-8) | **percent inhibition of catalytic activity at a single 1 μM screening concentration**, HotSpot radiometric filter-binding assay; 92 inhibitors × 758 kinases = 409 WT + 349 variants (311 mutants + 38 fusions) **[V]** (https://www.reactionbiology.com/insights/blog/mapping-the-kinome/, https://oncodaily.com/oncolibrary/lung-oncology/kirhub) | % inhibition, bounded [0,100] | **Primary mechanism-qualification surface (Scenario D)** | **Not Ki, Kd, pK, or DTA affinity — never label it as such.** Doubly censored (floor and ceiling). 38 fusions are not single-point mutants and must be excluded from the WT/single-mutant estimand. Reported measurement counts differ between secondary sources (~65,000 vs ~290,000); verify from the paper itself. |
| **Anastassiadis 2011** | percent *remaining* catalytic activity; 178 inhibitors × 300 kinases, functional assay **[V]** (https://www.nature.com/articles/nbt.2017, https://pubmed.ncbi.nlm.nih.gov/22037377/) | % activity | **Preferred cross-study replication for Saifudeen** | Shares the assay platform lineage (Reaction Biology / HotSpot; Haiching Ma is an author on both **[V]**), making the endpoints far more compatible than a cross-platform comparison. This is a better replication surface than Duong-Ly for endpoint reasons. |
| **Duong-Ly / Anastassiadis 2016 mutant profiling** (https://pubmed.ncbi.nlm.nih.gov/26776524/) | mutant kinase inhibitor profiling | % activity | Cross-study replication only | As the brief states: WT reference values originate from another study, so within-study WT–mutant matching is not available. Correct as written. |
| **PKIS / PKIS2** | percent inhibition, published probe sets | % inhibition | Auxiliary replication | Same endpoint family as above; same censoring issues. |
| **Platinum** (https://academic.oup.com/nar/article/43/D1/D387/2439527) | experimentally measured **ΔΔG of ligand binding** upon mutation, with the explicit curation requirement that WT and mutant be measured *by the same group under the same conditions*; ~1,000 mutations over 250 protein–ligand complexes **[V]** (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4935856/) | kcal/mol | **Gold-standard ΔΔ affinity surface, small** | Small but it is the only properly matched, physically united ΔΔ affinity resource. Its curation rule is exactly your "legally matched pair" criterion — adopt their standard verbatim. An alchemical FEP benchmark on it covers 134 mutations / 17 proteins / 27 ligands **[V]** (https://pubs.acs.org/doi/10.1021/acscentsci.8b00717), which gives a realistic sense of scale. |
| **MdrDB / mutation-resistance resources** | mixed, often structure-inferred | mixed | Auxiliary; verify provenance per entry | Frequently lacks experimental structures; DDMuffin had to predict them computationally **[V]** (https://www.biorxiv.org/content/10.1101/2025.07.19.665665v1.full). |
| **HIV genotype–phenotype resistance** | fold-change in IC50 vs reference | fold-change | Strong second mechanism surface, different protein family | Dense, many variants × many drugs. Excellent independent replication of the Scenario D claim outside kinases. |
| **ProteinGym / DMS** | **ligand-agnostic** organismal or biochemical fitness | assay-specific | **Negative control — not training signal** | See below. |
| **ChEMBL assay collections** | heterogeneous | mixed | Shared-trunk pretraining corpus, per-assay heads | Confidence scores must be used; FS-Mol restricted to single-protein-target assays with ≥32 measurements **[V]** (https://github.com/microsoft/FS-Mol). |
| **FS-Mol** | binary activity, task = protein | — | Standard few-shot comparison surface | 125 of 157 test tasks are kinases **[V]** (https://arxiv.org/pdf/2305.09481) — so FS-Mol "generality" is largely kinase generality. |

### 9.1 Two dataset-level controls that are currently missing and are high-value

**(a) ProteinGym as a ligand-invariance negative control.**
PLM-based variant effect scores are ligand-agnostic. Therefore: compute the correlation between your model's predicted WT→mutant change and a ligand-agnostic DMS fitness / PLM likelihood score. **If your predicted ΔΔ correlates with ligand-agnostic fitness as strongly as it does with the ligand-specific measured ΔΔ, you have not learned a ligand-conditioned interaction — you have learned "this mutation is disruptive."** That is a protein main effect wearing an interaction costume. This control is cheap and I believe it is the most likely way your programme would produce a false positive on real data.

**(b) The ATP-competition confounder for functional kinase panels.**
Percent inhibition at a fixed inhibitor concentration in a catalytic-activity assay depends, for ATP-competitive inhibitors, on the ratio [ATP]/Km,ATP via the Cheng–Prusoff relation. **Activating mutations change Km,ATP.** So a mutation can shift apparent percent inhibition *with no change whatsoever in inhibitor affinity*. A model can therefore "predict ΔΔ" by predicting activation status — again a ligand-agnostic protein main effect.
**Diagnostic:** a ΔΔ driven by this confounder is approximately *ligand-invariant across ATP-competitive inhibitors* and should be near-absent for allosteric/type-III inhibitors. So: (i) add a **ligand-invariant mutation-effect control arm** that predicts a single per-mutation scalar shared across all ligands, and require the full model to beat it; (ii) stratify results by inhibitor binding mode where annotation permits. If the correct-protein model does not beat the ligand-invariant-shift model, the apparent interaction signal is this confounder. **[I]**

These two controls belong in the frozen control matrix alongside shuffled and matched-wrong protein.

---

## 10. Evaluation hierarchy

Five levels. Each is specified with split unit, support sizes, allowed information, primary metrics, uncertainty method, required baselines, the claim it supports, and the claim it cannot support.

### L1 — Practical target-cold few-shot (Scenario A) — *primary performance benchmark*
- **Split unit:** protein target, clustered by **pocket similarity** (aligned binding-site identity), not global CD-HIT40.
- **Support:** k = 5, 10, 20, 40. Chemical relatedness between support and query **permitted and reported**.
- **Allowed:** sequence, family, MSA, predicted monomer structure, assay metadata, public pretrained models.
- **Primary metrics:** per-target Spearman; centered RMSE; activity-cliff sign accuracy; per-target Pearson. Secondary: uncentered RMSE, CI.
- **Uncertainty:** cluster bootstrap over protein components.
- **Required baselines:** tuned FS-CAP-class ligand-only; ActFound-class relative learner; Random Forest on ECFP per target; nearest-neighbour Tanimoto.
- **Supports:** a practical performance claim.
- **Cannot support:** any protein-mechanism claim.

### L2 — Target-cold zero-shot (Scenario B)
- **Split unit:** as L1. **Support:** k = 0.
- **Primary metrics:** enrichment at top 1% and 5%; global and per-target Spearman; **target-level calibration error reported separately from centered error**.
- **Required baselines:** ligand-promiscuity prior; global-mean predictor; family-mean predictor.
- **Supports:** screening utility.
- **Cannot support:** interaction attribution — a family-mean predictor can score well here.

### L3 — Scaffold-novel stratification (a *stratification of L1/L2*, not a separate run)
- Report every L1/L2 metric split by Bemis–Murcko scaffold novelty relative to training. **If the gain vanishes on scaffold-novel queries, classify it as retrieval/memorization.** This is the brief's own stopping rule G, and it is correct.

### L4 — Strict double-cold stress test (Scenario C) — *attribution only, no performance target*
- **Split unit:** pocket-cold × scaffold-cold. **Support:** k = 0 and k = 3, both reported.
- **Allowed:** deliberately restricted; declare exactly what is withheld.
- **Primary metrics:** effect size versus the power-calibrated null floor and oracle ceiling; **not** MSE against an absolute target.
- **Uncertainty:** component-level cluster bootstrap; report the number of independent clusters and the minimum detectable effect.
- **Supports:** existence of transferable protein-conditioned interaction under maximal shift.
- **Cannot support:** any performance claim. Do not report an MSE number from L4 as a headline.

### L5 — Mutation/selectivity biological positive control (Scenario D) — *primary mechanism surface*
- **Split unit:** **parent protein**, plus ligand scaffold for the double-cold variant. Never split within a parent.
- **Estimand:** signed WT→single-mutant change on the **same ligand**, restricted to legally matched pairs with compatible construct and substrate context, both members within the assay's dynamic range.
- **Support:** k = 0 for the ΔΔ head (the ligand is held fixed, so no support is needed for the contrast).
- **Primary metrics:** sign accuracy on pairs exceeding a pre-registered effect threshold; Spearman on signed ΔΔ; magnitude calibration. Report the number of *usable* pairs after censoring exclusion prominently — I expect this to be the binding constraint.
- **Uncertainty:** cluster bootstrap over parent proteins.
- **Required control arms:** full matrix (Section 11), including the two new arms from §9.1.
- **Supports:** the protein-conditioned interaction claim, on real data, with the ligand held exactly constant.
- **Cannot support:** a pK-affinity claim, because the endpoint is percent inhibition. Replicate on Platinum for the affinity-unit version.

### Metric discipline across all levels
- **Never report a single pooled metric across heterogeneous targets.** Report per-target distributions and their medians with cluster CIs.
- **Always report centered and uncentered error together.** The gap between them *is* the calibration component, and reporting it makes the "calibration shortcut" question quantitative instead of rhetorical.
- Report MSE/RMSE, CI, Spearman, Pearson, activity-cliff sign accuracy, enrichment (L2 only), centered error, and target-level calibration error.
- Apply the brief's rule G verbatim: if only MSE improves while CI/Spearman or centered error worsens, reject the candidate as a calibration shortcut.

---

## 11. Mechanism-attribution control matrix

Every mechanism claim (L4, L5) must run all arms on **matched rows, matched initializations, matched batch order, matched checkpoint policy, matched capacity**.

| Arm | What it holds constant | What it removes | What its failure proves |
|---|---|---|---|
| **Correct protein** | — | — | reference arm |
| **Protein-blind (ligand-only)** | ligand, support | all protein information | that protein information adds anything at all |
| **Additive-only** | main effects | the interaction term | that the signal is interaction, not additive |
| **Shuffled protein** | marginal protein distribution | protein–ligand correspondence | crude correspondence |
| **Family-preserving shuffled protein** | family identity | within-family specificity | that the signal is finer than family |
| **Similarity-matched wrong protein** | representation-space distance | correct identity | that the signal is not a smooth similarity gradient |
| **Capacity-matched random protein** | parameter count | information content | that the gain is not capacity |
| **Residue-position permutation** | composition | positional structure | that positional alignment matters |
| **No-interaction-head** | encoders | the interaction operator | that the operator is doing the work |
| **Ligand-invariant mutation shift** *(new, §9.1b)* | per-mutation scalar | ligand specificity | **that the effect is not the ATP-competition / activation confounder** |
| **Ligand-agnostic fitness predictor** *(new, §9.1a)* | PLM/DMS score | ligand conditioning | **that the effect is not "this mutation is disruptive"** |
| **Diagnostic oracle** | — | — | ceiling for effect-size calibration |

**Credit rule (from the brief, endorsed and strengthened):** a correct-protein gain counts only if it comes from *improving the correct branch*, measured against the fixed null floor. It cannot be credited if the margin arises from degrading corrupted branches. Operationally: pre-register the correct arm's absolute performance target relative to the shuffled-protein null distribution, and check it *before* computing any between-arm margin.

---

## 12. Realistic performance-development track

Runs in parallel with Track B. Deliverable: a defensible few-shot cold-target DTA advance, honestly labelled.

- **P1.** Reproduce and tune strong protein-blind baselines (FS-CAP-class, ActFound-class, per-target RF) on BindingDB-Ki with pocket-clustered target-cold splits and full scaffold-novelty stratification. **This defines the bar and is non-negotiable prerequisite work.**
- **P2.** Shared trunk with per-dataset heads; explicit external-data ablation; assay-balanced episodic sampling.
- **P3.** Add censoring-aware pairwise loss; evaluate on Davis where censoring is severe.
- **P4.** Add the level/shape split with gradient-ownership enforcement; report conflict statistics. Even with a protein path that turns out to be uninformative, this is a legitimate contribution: it makes the calibration/interaction decomposition measurable.
- **P5.** Add the protein path **only after** Track B admits a representation.
- **Promotion criterion at each step:** improvement over the previous step on L1 per-target Spearman *and* centered RMSE, with cluster CIs, and no degradation in scaffold-novel strata.

**Labelling discipline:** if the advance comes from P1–P4, it is a *few-shot bioactivity* contribution, not a protein-interaction contribution, and must be written up as such.

---

## 13. Strict scientific-validation track

- **B0.** Terminal synthetic diagnostic (Section 5.3). Bounded, one stage, power-calibrated, then closed either way.
- **B1.** Construct the Scenario D matched-pair benchmark: Saifudeen legal WT/single-mutant same-ligand pairs, both members in dynamic range, fusions excluded, split by parent protein. **Report the usable-pair count and the power analysis before running any model.**
- **B2.** Run the full control matrix on B1 with the simplest adequate models first — including non-deep baselines. The question at this stage is *does any method detect the effect*, not *does my method win*.
- **B3.** Representation comparison on B1 under identical head and budget: aligned-position difference; KLIFS one-hot + z-scales; local ESM; local ESM + physchem; global pooled ESM (negative); shuffled/random/permuted controls.
- **B4.** Replication: Anastassiadis (same assay lineage) → Platinum (affinity units, small) → HIV genotype–phenotype (different family).
- **B5.** Only after B2–B4 admit a representation: L4 double-cold attribution on BindingDB.
- **B6.** Relative-training innovation (antisymmetric protein-difference supervision) evaluated against endpoint-only and ligand-difference-only supervision, with gradient ownership and conflict reported.

---

## 14. Staged roadmap with promotion and stop criteria

| Stage | Action | Promote if | Stop / redirect if |
|---|---|---|---|
| **S0** | Null-space trajectory measurement on existing artifacts | trajectory identified as flat or rising | — (measurement only; ~1 day) |
| **S1** | Terminal synthetic diagnostic: SGD-from-zero, Adam, span-reparameterized arms; power-calibrated thresholds | learner recovers ≥ pre-registered fraction of oracle–null gap on double-cold | if all arms fail → **close bilinear family**, move to a different interaction operator; **do not** redesign the truth |
| **S2** *(parallel with S1)* | Build B1 benchmark + power analysis | usable pairs ≥ minimum detectable-effect requirement | if usable pairs are too few after censoring exclusion → Scenario D on Saifudeen is underpowered; pivot primary mechanism surface to HIV genotype–phenotype or Platinum |
| **S3** | B2 control matrix with simple models | correct protein beats **every** control arm including the two new ones | if correct protein ≈ matched-wrong or ≈ ligand-invariant-shift → **reject the protein representation family**; protein-mechanism claims closed; performance track continues |
| **S4** | B3 representation comparison | one local representation beats all controls and beats global-pooled | if only global-pooled works → the signal is target identity; reject |
| **S5** | B4 replication | effect reproduces in ≥ 1 independent surface | if it does not replicate → single-study artifact; withdraw |
| **S6** | B6 relative-training innovation | antisymmetric protein-difference supervision beats endpoint-only and ligand-difference-only on B1, with no gradient-ownership violation | if not → keep as ablation, not as the core innovation |
| **S7** | L4 double-cold attribution on BindingDB | effect exceeds power-calibrated floor with cluster CIs | if not → the claim is scoped to Scenario D, stated honestly |
| **S8** | Merge into performance track (P5) | L1 improves with the protein path and does not degrade scaffold-novel strata | if L1 gains vanish on scaffold-novel → retrieval, not interaction |

**Parallelism rule:** the performance track (P1–P4) runs continuously and is never blocked by Track B. The only gate is *labelling*: no protein-interaction claim may be made without the corresponding Track B admission.

---

## 15. Risks: leakage, assay confounding, target-ID memorization, scaffold recall

| Risk | Mechanism | Detection | Mitigation |
|---|---|---|---|
| **Source overlap leakage** | ChEMBL and BindingDB share primary literature and patents; ActFound had to filter test assays by cross-source correlation **[V]** | cross-source assay correlation screen | deduplicate before splitting; report the number removed |
| **Pocket leakage under CD-HIT40** | global identity < 40% does not imply pocket dissimilarity, especially in kinases | compute aligned-pocket identity between train and test targets and report the distribution | **split by pocket similarity, not global identity** |
| **Target-ID memorization** | protein encoder learns an identity code; performance is retrieval | matched-wrong and family-preserving-shuffle arms; low-rank bottleneck ablation | difference-only protein representation (Rank 1); capacity bottleneck on the protein path |
| **Scaffold recall** | query scaffold seen in training against a different target | Bemis–Murcko novelty stratification at every level | report all metrics stratified; apply rule G |
| **Assay confounding — endpoint mixing** | IC50/Ki/Kd/%inhibition merged into one regression target | label-semantics audit | per-dataset heads; never merge |
| **Assay confounding — ATP competition** | Km,ATP changes with activating mutations shift apparent %inhibition without affinity change | ligand-invariance of predicted ΔΔ; stratify by inhibitor binding mode | ligand-invariant-shift control arm (§9.1b) |
| **Ligand-agnostic fitness shortcut** | model predicts "disruptive mutation" rather than ligand-specific ΔΔ | correlation of predictions with DMS/PLM fitness scores | ligand-agnostic fitness control arm (§9.1a) |
| **Censoring artifact** | Davis floor mass, panel floor/ceiling saturation | fraction of pairs at bounds; censored vs uncensored metric gap | censoring-aware likelihood; dynamic-range restriction; report usable-pair counts |
| **Calibration shortcut** | MSE improves via level, CI/Spearman flat or worse | centered vs uncentered gap; target-level calibration error | reject per rule G |
| **Power masquerading as a null** | too few independent clusters to detect a real effect | pre-registered minimum detectable effect | power-calibrated thresholds (§5.2.6); never report an underpowered null as a negative finding |
| **Underpowered PASS via degraded controls** | margin created by corrupting negative arms | absolute check of the correct arm against the null floor before any margin | credit rule (§11) |

---

## 16. Data, compute and experimental requirements

**Data.** All required resources are public. Saifudeen supplementary data / KiRHub portal; Anastassiadis 2011 supplementary; Platinum; KLIFS (85-position alignments, weekly-updated **[V]**); BindingDB; Davis; KIBA; ChEMBL; FS-Mol; ProteinGym; Stanford HIVdb. No new experiments required for S0–S7.

**The binding constraint is not data volume — it is *usable pair count* after legal matching and censoring exclusion.** From the panel dimensions (92 inhibitors × 311 point mutants, minus fusions, minus construct/substrate mismatches, minus pairs where either member is saturated at the floor or ceiling), the usable set could plausibly fall by an order of magnitude from the nominal cell count. **Compute this number before designing anything else.** If it is small, the entire mechanism programme is power-limited and the design must change — that is a decision-relevant fact and it is obtainable in hours.

**Compute.** Modest. Nothing in S0–S6 requires large-scale training: model sizes < 10M parameters, KLIFS features are 85 positions, panel data is order 10⁴–10⁵ rows. Single-GPU days, not weeks. ESM embedding extraction is a one-time cost. S7 on BindingDB is larger but still single-node. **GPU capacity is not the bottleneck and should not be treated as a success metric.**

**Human effort.** The dominant cost is curation: legal-pair construction, construct/substrate compatibility, censoring annotation, pocket-similarity clustering. Budget accordingly; this is where errors of the kind already found in the synthetic programme are most likely to recur.

---

## 17. The single highest-value next action

**Construct the Scenario D legal-pair benchmark and run the power analysis — before running any model, and before waiting for the synthetic diagnostic.**

Concretely: take the Saifudeen panel; restrict to WT/single-point-mutant pairs on the same ligand with compatible construct and substrate context; exclude fusions; exclude pairs where either member is saturated at the assay floor or ceiling; cluster by parent protein; and report (i) the number of usable pairs, (ii) the number of independent parent clusters, (iii) the distribution of effect sizes, and (iv) the minimum detectable effect under a parent-level cluster bootstrap.

Why this and not the synthetic diagnostic **[I]**:

- The synthetic question has become an *engineering* question with a known answer shape. You already know a recoverable answer exists (the oracle finds it) and you already know the learner fails cold. The remaining work discriminates between named optimizer/parameterization mechanisms. It is worth one bounded stage, but its outcome cannot change what is true about biology.
- The Scenario D power number is a **genuine unknown that gates the entire programme** and is obtainable in hours rather than weeks. If the usable-pair count is too small, then no amount of synthetic qualification, representation engineering, or training innovation can produce a credible mechanism claim on that surface, and you need to know that *now* so you can pivot to HIV genotype–phenotype or Platinum.
- It is also the cheapest possible falsification of the whole protein-conditioned programme. Once the benchmark exists, running the control matrix with simple, non-deep models answers "does *any* method detect a ligand-specific mutation effect that survives matched-wrong-protein and ligand-invariant-shift controls?" If the answer is no for everything, the programme's central hypothesis is unsupported by the best available data, and that conclusion is worth far more than another synthetic stage.

Run S1 in parallel. It is bounded and independent.

---

## 18. Concise final architecture hypothesis

Offered conditionally: this is the architecture *supported by the analysis above*, and it should be built only as Track B admits its components.

**Core idea: support fixes the level; protein fixes the shape.** This single principle is what the existing evidence actually points to — k = 0 error is dominated by calibration; k ≥ 2 gains come from ligand similarity; both of those are *level* phenomena. The protein-conditioned claim lives entirely in the shape term, so it should be architecturally isolated there.

- **Ligand encoder:** fragment/atom tokens (GNN or fragment-decomposed fingerprint). Not pooled before the interaction operator.
- **Protein encoder:** KLIFS 85 aligned positions; per-position features = [z-scale physicochemical, local ESM window, conservation weight]. **Never globally pooled before the interaction operator.**
- **Interaction operator:** low-rank position × fragment bilinear map or sparse (entmax) cross-attention, bottlenecked to rank r ≪ 85 so the protein path cannot carry an identity code.
- **Two output paths with enforced gradient ownership:**
  - *Level head* — receives gradient only from the absolute-endpoint loss; deliberately capacity-limited; consumes support labels when k ≥ 1.
  - *Interaction head* — receives gradient only from the relative losses; consumes no support labels; this is the zero-shot path.
- **Dual relative supervision:**
  - within-target ligand-difference (established to work — ActFound **[V]**);
  - **within-ligand protein-difference, exchange-antisymmetric and identity-zero by construction** (the novel component, and simultaneously an attribution guarantee, since ligand-only and identity solutions vanish under the constraint).
- **Likelihood:** censoring-aware (interval/Tobit) on both endpoints and differences.
- **Heads:** one per dataset, preserving label semantics; shared trunk; external-data contribution ablated and reported.

**Why this configuration and not the current bilinear one.** It retains the bilinear interaction operator — which the synthetic work suggests is capable once the identifiability/conditioning pathology is fixed — but removes the three failure modes the evidence has repeatedly exposed: global pooling (which makes the protein a target identifier), unconstrained capacity in the protein path (which permits memorization), and an undifferentiated loss (which lets calibration masquerade as interaction). The antisymmetric protein-difference supervision is the only component that is simultaneously novel, gradient-trained, and structurally incapable of expressing the shortcuts you have spent this programme ruling out.

**Predicted failure mode, stated in advance:** the most likely way this architecture produces a false positive is via the ligand-agnostic fitness / ATP-competition confounders in §9.1. The two new control arms exist specifically to catch it. If it survives those on B1 and replicates on Anastassiadis and Platinum, the mechanism claim is credible.

---

## Appendix: primary sources consulted

**Methods**
- AdaMBind — https://www.nature.com/articles/s41467-026-70554-5 · code https://github.com/Moohyun-w/AdaMBind
- MetaDTA — https://openreview.net/pdf?id=yzlif16IASM · https://iclr.cc/virtual/2022/8615
- FS-CAP — https://pubs.acs.org/doi/10.1021/acs.jcim.4c00485 · https://arxiv.org/abs/2311.16328
- ActFound — https://www.nature.com/articles/s42256-024-00876-w · https://www.biorxiv.org/content/10.1101/2023.10.30.564861v1.full.pdf · https://github.com/BFeng14/ActFound
- ZeroBind — https://www.nature.com/articles/s41467-023-43597-1
- PSICHIC — https://www.nature.com/articles/s42256-024-00847-1
- CORDIAL — https://www.pnas.org/doi/10.1073/pnas.2508998122
- DrugBAN (cold-pair split, hidden ligand bias) — https://arxiv.org/pdf/2208.02194
- PCM foundations — https://pubs.rsc.org/en/content/articlelanding/2011/md/c0md00165a · https://academic.oup.com/ib/article/6/11/1023/5199186

**Benchmarks and evaluation**
- FS-Mol — https://github.com/microsoft/FS-Mol · https://datasets-benchmarks-proceedings.neurips.cc/paper/2021/file/8d3bba7425e7c98c50f52ca1b52d3735-Paper-round2.pdf
- MoleculeACE (activity cliffs) — https://pubs.acs.org/doi/10.1021/acs.jcim.2c01073
- DTI OOD assessment / leakage — https://www.biorxiv.org/content/10.1101/2022.04.20.488898.full.pdf
- DTI benchmark survey — https://academic.oup.com/bib/article/26/5/bbaf491/8260789

**Data resources**
- Saifudeen et al., Nat Biotechnol 2026 — https://doi.org/10.1038/s41587-026-03090-8 · summaries: https://www.reactionbiology.com/insights/blog/mapping-the-kinome/ · https://oncodaily.com/oncolibrary/lung-oncology/kirhub
- Anastassiadis et al. 2011 — https://www.nature.com/articles/nbt.2017 · https://pubmed.ncbi.nlm.nih.gov/22037377/
- Davis et al. 2011 / selectivity commentary — https://www.nature.com/articles/nrd3647
- Duong-Ly et al. 2016 (mutant kinase profiling) — https://pubmed.ncbi.nlm.nih.gov/26776524/
- KLIFS — https://pubs.acs.org/doi/10.1021/jm400378w · https://academic.oup.com/nar/article/44/D1/D365/2502606
- Platinum — https://academic.oup.com/nar/article/43/D1/D387/2439527 · mCSM-lig https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4935856/
- FEP ΔΔG benchmark on Platinum — https://pubs.acs.org/doi/10.1021/acscentsci.8b00717
- DDMuffin (ΔΔG transfer learning, MdrDB) — https://www.biorxiv.org/content/10.1101/2025.07.19.665665v1.full

**Learning theory / optimization**
- ANIL / rapid learning vs feature reuse — https://arxiv.org/abs/1909.09157 · https://openreview.net/forum?id=rkgMkCEtPB
- PCGrad — https://arxiv.org/pdf/2001.06782
- CAGrad (critique of PCGrad) — https://arxiv.org/pdf/2110.14048
- ProteinGym / PLM zero-shot variant effect — https://arxiv.org/pdf/2507.02624 · https://arxiv.org/pdf/2405.06729
