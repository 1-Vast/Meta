# MetaSieve core model research summary

Date: 2026-08-11
Scope: replacement of the failed core modules of the few-shot DTA model.
Status: derivation and implementation complete; **no new training was executed
this cycle** (see §6).

```text
CENTERING_GAUGE_COLLAPSE_DERIVED
TBASIS_PARTNER_LOSS_LOCALIZED_TO_COMPOSITION_MARGINALIZATION
CENTERED_SECTION_NECESSARY_BUT_PROVABLY_INSUFFICIENT
CROSSED_CONTRAST_MECHANISM_SPECIFIED_AND_IMPLEMENTED
NO_NEW_TRAINING_EVIDENCE_THIS_CYCLE
BIOLOGICAL_REPRESENTATION_REPAIR_REQUIRED
```

---

## 1. Current failure mechanism

Two results, one derived and one localized, replace the previous informal
diagnosis ("the section is mostly calibration").

### 1.1 The centring-gauge collapse

Write the pair coordinate as an affine function of a ligand embedding,

$$m(P,L)=A_P\,\varphi(L)+b_P,\qquad A_P\in\mathbb R^{d\times H}.$$

Apply the R2 minimal repair (ligand-only population $\mu_L$, explicit intercept
$b_t=\frac1k\mathbf 1^\top r$ with $r=y_S-\mu_L(L_S)$, centred positive ridge):

$$\hat y_q=\mu_L(L_q)+b_t+(m_q-\bar m_S)^\top M_c^\top\big(M_cM_c^\top+\lambda I\big)^{-1}(r-b_t\mathbf 1).$$

Substituting the affine form:

$$M_c=\big[A_P(\varphi_i-\bar\varphi_S)\big],\qquad m_q-\bar m_S=A_P(\varphi_q-\bar\varphi_S),$$

so **every** occurrence of the protein is through the Gram metric
$G_P=A_P^\top A_P$; $r$ and $b_t$ contain no protein term at all. Therefore:

| condition on $A_P$ | consequence |
|---|---|
| $b_P$ arbitrary | annihilated exactly by centring |
| $A_P$ orthogonal (shared $G$) | correct/correct $\equiv$ wrong/wrong, **identically**; the correct-vs-wrong-protein contrast is exactly zero |
| $A_P=c_P A$ | protein reduces to one scalar shrinkage knob, vanishing as $\lambda\to0$ |

This is decision-relevant and, as far as this cycle can establish, was not
previously stated: **the recommended minimal repair removes the protein's only
additive channel.** Making the section calibration-orthogonal is necessary to
stop mis-attributing the offset, but it *increases* wrong/wrong invariance
unless the pair representation makes the within-task ligand **metric**
protein-specific. The V1 2×2 factorial (`V0 wrong/wrong 1.765` vs
`correct 1.800`; `V1-B 3.866` vs `3.890`) is the empirical signature of exactly
this regime.

Executable witnesses: `tests/test_v2_crossed_contrast_section.py::
test_centered_section_is_blind_to_protein_when_the_metric_is_shared` and
`::test_anisotropic_metric_restores_protein_dependence`.

### 1.2 Where the 288-D T-BASIS destroys partner identity

The frozen readout (`research/e0_identifiability/run_tbasis_radial.py`,
`aggregate_basis` / `slot_composition`) is

$$F[a,r,k]=\frac1A\sum_i\sum_s c_i[a]\;\mathrm{comp}_P[s,r]\;\mathrm{radial}[i,s,k],$$

with 8 atom channels, 6 residue chemistry classes, 6 radial centres,
$\mathrm{MECHANISM\_RESIDUE\_SLOTS}=128$.

Three structural facts follow from the code, not from a fit:

1. `slot_composition` bins residues by **relative sequence position**
   (`index * 128 // len(sequence)`) and row-normalizes. It is a
   length-normalized, order-free residue-class composition profile of the
   **whole chain**. Pocket residues are scattered in sequence, so every slot
   averages pocket with non-pocket residues; the pocket is never isolated.
2. Each row of $\mathrm{comp}_P$ is a probability simplex, so
   $\sum_r F[a,r,k]$ is **exactly** independent of $\mathrm{comp}_P$. The
   residue-chemistry axis can only redistribute mass that the marginal already
   fixes.
3. If $\mathrm{comp}_P\approx f$ (protein-independent background composition),
   then $F\approx f\otimes D(P,L)$ with
   $D[a,k]=\frac1A\sum_i c_i[a]\sum_s\mathrm{radial}[i,s,k]$ — a
   slot-marginalized contact-mass profile. Two consequences: five of the six
   residue dimensions carry no independent information (rank-1 along $r$), and
   the protein survives only as an **aggregate contact-degree statistic**.

Consequence 3 is precisely the shape that produces a near-additive design.
It explains, without any new experiment, the three recorded observations:
additive explained fraction `0.9807`, fixed-ligand partner dispersion `0.0513`,
and the A1 selectivity probe in which T-BASIS (`0.926`) was *worse* than the
zero predictor (`0.585`) and the rewired coupling null (`0.635`).

The R2 correction stands: the frontend does condition distance logits on ESM
residue states. The loss is not "the frontend never sees residues" — it is
**the aggregation marginalizes them against a near-constant weight**.

`research/meta_fewshot/v2_tbasis_partner_localization.py` measures all four
quantities label-free (composition dispersion, residue-mode rank-1 fraction,
reconstruction of the 288-D vector from its 48-D residue marginal, partner
shift direction) and emits one of four localization verdicts.

---

## 2. Modules retained

| Module | Status | Reason |
|---|---|---|
| Frozen law operator $\mathsf A(F,z)=K(B(z)F(z))$ | unchanged, disconnected | no admitted biological $z$; §10 |
| Source-learned function family → support-identifiable section | retained as the core principle | unchanged; §10 |
| Protein-as-task episodic learning | retained | matches the estimand and the AdaMBind protocol level |
| Positive closed-form differentiable ridge section | retained | support-only, stable, primal/dual exact |
| Support budget $k\le5$, $d\le5$ | retained | frozen capacity ceiling (F20/CP-3) |
| Support-conditioned adaptation improves real BindingDB MSE | retained as fact | `1.916` vs `8.711` (d=0) on main-v0 |
| Explicit intercept + centred section | retained, and **re-scoped** | necessary for honest attribution; §1.1 shows it is not a partner-specificity mechanism |
| Frozen ESM residue states, GINE atom states, P1B distance/contact bridge | retained as inputs | the failure is localized downstream of them |
| Full control battery (zero / foreign / permuted / wrong-protein / wrong-wrong) | retained and extended | §7 |

---

## 3. Modules rejected or replaced

| Module | Decision | Evidence |
|---|---|---|
| `slot_composition` + `aggregate_basis` (288-D T-BASIS readout) | **replaced** | §1.2; near-additivity and the negative A1 probe are consequences of the aggregation, not of the encoder |
| Uncentred section | replaced | calibration confound (R2-E0) |
| "Centred section alone fixes partner specificity" | **rejected as a mechanism** | §1.1 derivation |
| Free full-rank pair prior $\mu(P,L)$ MLP | rejected | V1 pair-d0 `3.806` worse than ligand-d0 `3.084` |
| Larger pair MLP / larger $d$ / wider readout | rejected | V1-A/V1-B both worse than v0 in absolute MSE |
| RFMS reserved block | remains closed | $\Xi(c_{\text{correct}}-c_{\text{wrong}})\neq0$ unproven |
| Q-PMA, unrestricted MAML, second foundation model | remain closed | the failure is not support-attention capacity |
| Random wrong-protein negative affinity labels | forbidden | fabricates non-binders and optimizes against the control Gate |
| CSMO integration | remains closed | no admitted biological $z$ |

---

## 4. Literature-supported reasoning

The failure has an exact name in the meta-learning literature. In
[Meta-Learning without Memorization](https://arxiv.org/abs/1912.03820)
(Yin, Tucker, Zhou, Levine, Finn, ICLR 2020) the meta-learner outputs one model
that solves every task while ignoring the support set, because the task
distribution is not mutually exclusive. MetaSieve exhibits the *complementary*
degeneracy: the support set is used, but only for its mean. Their remedy —
regularize so that adaptation must carry information — motivates the design
principle here, though their information-theoretic regularizer is not adopted
directly because it does not distinguish *which* information adaptation carries.

For separating assay calibration from within-task structure, the strongest
precedent is
[MBP / ChEMBL-Dock](https://academic.oup.com/bib/article/25/1/bbad451/7469349),
which classifies relative rankings **within** a bioassay precisely because
inter-assay offsets are systematic and intra-assay comparisons are not. This is
the same principle as the explicit intercept, and it establishes that the
within-task contrast is a legitimate supervision target rather than a trick.

For cross-target supervision, direct regression on the **affinity difference
between two targets for the same compound** is reported to beat deducing
selectivity from two independent predictors
([A1/A2A selectivity](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7222572/);
see also the [kinase profiling review](https://pmc.ncbi.nlm.nih.gov/articles/PMC11419124/)).
That is the first-order cross-protein difference. The mechanism proposed in §9
is its second-order form.

Nothing in the searched literature applies a *crossed* (second-order) measured
difference through a support-conditioned closed-form section. The individual
ingredients are all prior art and must be cited as such: closed-form adaptation
([R2-D2](https://openreview.net/forum?id=HyxnZh0ct7),
[ALPaCA](https://arxiv.org/abs/1807.08912)); bilevel deep-kernel adaptation for
molecular few-shot
([ADKF-IFT](https://arxiv.org/abs/2205.02708), ICLR 2023); protein-task
episodic DTA ([AdaMBind](https://www.nature.com/articles/s41467-026-70554-5));
interaction-aware physicochemical pair representations
([PSICHIC](https://www.nature.com/articles/s42256-024-00847-1)); side-information
factorization ([Macau](https://arxiv.org/abs/1509.04610)) and PCM cross-terms.
The novelty claim in §9 is therefore restricted to the training-level
construction, and is stated as a hypothesis to be tested, not as priority.

---

## 5. Existing-data exploitation findings

* **Where partner information still plausibly exists.** ESM residue states are
  consumed *before* the loss (they condition the bridge distance logits), and
  the `deranged/natural` feature shift of `2.878` shows the statistic does
  respond to protein substitution. What is missing is an *affinity direction*,
  not partner sensitivity. The loss is at the aggregation (§1.2), which is a
  repairable stage, not at ESM.
* **KLIFS pockets.** The XP2 package already carries 147 length-85 KLIFS pocket
  strings, while `esm2_t30_kinase_pocket85.npz` covers only 82 kinases. A
  pocket-restricted residue set is therefore available for the kinase panels
  without new data acquisition — but the bank must be regenerated and hashed
  before any coverage claim is made.
* **Measured crossed supervision exists.** The V1 source census found 1,002
  within-panel ligand groups and 1,820 same-panel/same-ligand partner groups
  crossing CD-HIT-40 families. Metz XP2 (32,849 measured cells, 31,775
  interaction df) and PDSP core (10,701 / 8,313) are dense crossed panels.
* **None of it is fresh.** Every panel closes to one dependency component
  (Metz, PDSP) or has 91.2% of cells in the largest of eight (BindingDB); 117
  of 129 BDB proteins and 825 of 2,845 ligands overlap main-v0 exactly. These
  are **development** supplies. `v2_rectangle_census.py` enforces the
  distinction by reporting dependency components and largest-component share
  alongside the rectangle count.
* **Consumed splits.** main-v0 test, meta-validation and A1 are consumed. Every
  script written this cycle drops `meta_test` rows *before* reading `pK`.

---

## 6. Experiments actually executed

**None.** This must be stated plainly.

The Cowork Linux sandbox failed to start for the entire session
(`session disk not found`), so no Python could be executed: no training, no
numeric diagnostic, no test run. Every number quoted in this document is taken
from an earlier, already-recorded MetaSieve artifact and is cited as such. No
number in this document is new, and no verdict here rests on a new measurement.

What was produced instead is the complete, runnable specification of the
mandatory experiments:

| Artifact | Purpose | State |
|---|---|---|
| `research/meta_fewshot/v2_nested_section_ablation.py` | the mandated 5-arm nested comparison with full control battery | implemented, **not run** |
| `research/crossed_interaction/v2_crossed_transfer_ladder.py` | does a crossed interaction term exist **and transfer** (C0–C5) | implemented, **not run** |
| `research/meta_fewshot/v2_tbasis_partner_localization.py` | label-free localization of the §1.2 loss | implemented, **not run** |
| `research/meta_fewshot/v2_rectangle_census.py` | label-free supply census for crossed supervision | implemented, **not run** |
| `research/meta_fewshot/v2_crossed_contrast_section.py` | the candidate mechanism (§9) | implemented, **not trained** |
| `tests/test_v2_crossed_contrast_section.py` | algebraic guarantees and the §1.1 witness | implemented, **not run** |

Execution order is fixed and must not be reordered:

```text
v2_tbasis_partner_localization      (label-free, cheap, no Gate)
v2_rectangle_census                 (label-free, decides whether X-CON has supply)
v2_crossed_transfer_ladder          (does any crossed term transfer at all?)
v2_nested_section_ablation          (the mandated intercept-vs-section comparison)
    only then: train v2_crossed_contrast_section
```

The transfer ladder is deliberately placed **before** the nested ablation.
Its arm C4 uses no protein features whatsoever and estimates the target factor
from the same $k$ support labels. If C4 does not beat the intercept baseline C1
and the additive ligand prior C3 under protein-family-cold and scaffold-cold
evaluation, then there is no transferable crossed structure for any biological
encoder to represent, and the correct action is to stop rather than to repair
the frontend.

---

## 7. Metrics and controls

Primary quantities, all reported at **target macro and dependency-component /
CD-HIT-cluster macro**, with one-sided 95% bootstrap lower bounds:

* `MSE(centred section) − MSE(pair intercept)` — the registered primary;
* within-target Pearson, Spearman and concordance index — absolute MSE
  improvement alone is explicitly insufficient;
* the full support battery: correct vs zero, foreign, permuted;
* the protein battery: correct vs wrong query, wrong support, and the absolute
  Gate `MSE_ww − MSE_cc > 0`;
* the 2×2 synergy $\Delta_{\rm BioMeta}=(M_{cf}-M_{cc})-(M_{wf}-M_{wc})$ as a
  mechanism indicator only, never as a substitute for the absolute Gate.

Two controls are **added** this cycle because the existing battery cannot
separate family memorization from partner identity:

1. **Within-family wrong protein.** Every wrong-protein donor used so far came
   from a *different* CD-HIT-40 cluster. A model that only learns a family-level
   metric passes that control. The new required donor is a different sequence
   from the *same* CD-HIT-40 cluster.
2. **Calibration-leakage check.** Shifting every support label by a constant
   must move the prediction by exactly that constant
   (`calibration_leakage`). Any deviation is calibration re-entering the
   section.

`PAIR_POP_INTERCEPT` is a permanent standing baseline and stays in every report
regardless of outcome.

### Adversarial accounting

Every candidate gain must be assigned to one of these before it is called
biology:

| Adversarial explanation | Discriminating control |
|---|---|
| target intercept | `PAIR_POP_INTERCEPT` arm; centred section by construction cannot absorb a constant |
| ligand memorization | scaffold-cold query relative to support **and** to training |
| family memorization | new within-family wrong-protein donor (above) |
| assay / document effect | document-blocked rectangles; component macro |
| scaffold overlap | Murcko-scaffold closure in the split and in episode drawing |
| support-label distribution | permuted-support arm; foreign-support arm |
| coordinate self-consistency (the V1 finding) | absolute `MSE_ww − MSE_cc` Gate, not the 2×2 interaction |
| dependency-driven significance | component-macro bootstrap; largest-component share reported |
| selection on consumed data | `meta_test` physically dropped before any label read |

---

## 8. Failed attempts and why they failed

Three repair attempts are on the record, each rejected for a distinct reason.

**Attempt 1 — more readout capacity (V1-A / V1-B).** A shared full-rank pair
prior plus a residual pair adapter, with measured within-panel and partner
difference losses. Rejected: pair-d0 `3.806` was worse than ligand-d0 `3.084`,
and V1-A/V1-B absolute MSE (`4.204` / `3.890`) were worse than v0 (`1.800`).
V1-B did stabilize the permuted contrast, but the wrong/wrong Gate still failed
(`0.066`, LCB `−0.123`). *Diagnosis:* the auxiliary contrast loss was applied to
the **prior**, which is exactly the additive channel that the crossed estimand
removes; it could not reach the section.

**Attempt 2 — the reserved-fiber section (RFMS).** Rejected before training:
non-constant $c_0$ is not a partner-specificity certificate, since the exposure
$\Xi=b-a(A^\top A+\lambda I)^{-1}A^\top C$ may annihilate the partner
difference; $d_{\rm support}\le k<d$ is incompatible with the standing
$d\le5,k=5$; and its wrong-protein hinge would have fabricated non-binders and
optimized directly against the final control Gate.

**Attempt 3 — the centred section as the partner fix.** Rejected **by
derivation** this cycle (§1.1), not by a training run: centring annihilates
$b_P$ and leaves only $G_P$, so on a representation whose $A_P$ is nearly
protein-independent it makes wrong/wrong invariance exact rather than
approximate. It remains mandatory as the attribution fix and as the registered
baseline arm; it is no longer a candidate mechanism for partner specificity.

Per the iteration rule, attempt 4 must change the information source and the
estimand together, not enlarge the network. §9 does both: the information source
becomes pocket-resolved residue states entering a protein-conditioned **metric**,
and the estimand becomes the measured crossed difference.

---

## 9. Proposed architecture and training mechanism

The target shape required by the brief,

$$\text{explicit calibration channel}+\text{interaction-specific biological function family}+\text{support-identifiable centred section},$$

is instantiated as follows. Implementation:
`research/meta_fewshot/v2_crossed_contrast_section.py`.

### 9.1 Calibration channel (retained, unchanged)

$b_t=\frac1k\mathbf 1^\top r$ is explicit; the section sees only $r-b_t\mathbf 1$
and centred coordinates. A constant support shift moves the prediction by
exactly that constant. Label-dependent degrees of freedom per task:
$1+\operatorname{rank}(M_c)\le k$.

### 9.2 Interaction-specific biological function family

Replace the pooled 288-D readout by a protein-conditioned **map**, not a
protein-conditioned vector:

$$m(P,L)=W(P)^\top\varphi(L),\qquad W(P)=W_0+\sum_{j=1}^{J}g_j(P)\,W_j,\qquad \|g(P)\|_2\le1 .$$

$g(P)$ is computed from **pocket-restricted** residue states (KLIFS pocket for
kinases; gated top-$m$ residue selection otherwise), never from a global
length-normalized composition. $J$ is small (default 8): the protein selects a
member of a small learned family of interaction bases; it does not receive a
free embedding, and no wider MLP is introduced.

This is the minimum change that defeats §1.1: the centred kernel becomes
$\langle\varphi_q-\bar\varphi_S,\,G_P(\varphi_i-\bar\varphi_S)\rangle$ with
$G_P=W(P)W(P)^\top$ genuinely protein-dependent. $G_P$ is invariant under the
latent gauge $W(P)\mapsto W(P)O$, so it is a legitimate reported quantity.

**The encoder alone is explicitly not the innovation.** An unsupervised $W(P)$
can still collapse to a shared metric; §9.3 is what forbids it.

### 9.3 The training-level mechanism — X-CON (crossed contrast)

For a measured rectangle $\{P_1,P_2\}\times\{L_1,L_2\}$ with all four cells
measured under the same document/panel, and legitimate $k$-shot support
episodes $S_1,S_2$ that exclude both rectangle ligands and their scaffolds:

$$\mathcal L_\times=\Big(\big[\hat y(P_1,L_1)-\hat y(P_1,L_2)\big]-\big[\hat y(P_2,L_1)-\hat y(P_2,L_2)\big]-\Delta_\times y\Big)^2,$$

where $\hat y$ are the model's **own few-shot predictions** through the centred
section and $\Delta_\times y=y_{11}-y_{12}-y_{21}+y_{22}$.

The total objective is $\mathcal L_{\rm episodic}+w\,\mathcal L_\times$: the
absolute term is required because $\Delta_\times$ identifies no offset, and the
crossed term is required because the absolute term is satisfiable by
calibration alone.

Why this and not another auxiliary loss:

* $\mu_L$ cancels identically — it appears with the same value in both protein
  brackets. **The ligand prior cannot reduce this loss.**
* $b_t$ cancels identically within each protein. **Target-level calibration
  cannot reduce this loss.** (Test:
  `test_ligand_prior_and_task_calibration_cannot_produce_a_crossed_prediction`
  — with the section switched off the predicted crossed difference is exactly
  zero.)
* Every protein and ligand main effect is annihilated on the label side too.
* What survives is $G_{P_1}-G_{P_2}$. If the model sets $G_P$ constant,
  $\mathcal L_\times$ is bounded below by the dispersion of $\Delta_\times y$
  and cannot be reduced.
* It uses **only measured cells**. No unmeasured pair is treated as a
  non-binder.

This is a constraint on the existing predictor, not additional capacity, and it
is applied to the object that actually carries the failure — which is precisely
where V1-B went wrong (§8).

### 9.4 Why both correct pairing and correct identity become necessary

* **Correct support (ligand, label) pairing.** $\Delta_\times$ is a difference
  of *within-task* ligand contrasts. The section's within-task slope is
  identified only from correctly paired support rows; permuting support labels
  preserves $b_t$ exactly and destroys the slope. A model trained on X-CON must
  therefore separate correct from permuted support, which the incumbent barely
  does (v0 permuted contrast `0.130` target macro).
* **Correct protein identity.** After centring, the protein enters only through
  $G_P$; X-CON supervises exactly $G_{P_1}-G_{P_2}$ against a measured,
  main-effect-free quantity. A model that succeeds at X-CON provably has
  non-constant $G_P$, hence correct/correct $\neq$ wrong/wrong by construction —
  which is the Gate the project has failed four times.

---

## 10. Exact theory-to-model mapping

The frozen operator is unchanged and remains disconnected.

| Frozen object | Model realization | Status |
|---|---|---|
| Family $\mathcal F$ | source-learned function family: $\varphi$, $\{W_j\}$, $\mu_L$ | trained on source components, then frozen and hashed |
| Archive | meta-training tasks (one protein = one task) | unchanged |
| Auxiliary $c_b$ (fiber restriction) | the protein, entering as the family selector $g(P)$ and hence the metric $G_P$ | **this is the change**: $c_b$ now restricts the *function family*, which is its defined role (CI-A1), instead of contributing an additive offset |
| Support $S_b$, $k\le5$ | the $k$ measured support pairs | unchanged |
| Section cut $S_\varepsilon(\tilde y\mid c)$ | explicit intercept + centred positive ridge | unchanged |
| Support capacity: at most $k$ continuous dimensions (F20/CP-3) | $1+\operatorname{rank}(M_c)\le k$, enforced in code | preserved |
| Gauge invariance (IB-4/CR-5): no reported quantity may depend on latent coordinates | reported metric is $G_P=W(P)W(P)^\top$, invariant under $W\mapsto WO$ | preserved and unit-tested |
| Permutation symmetry of the support (IB-1/AP-1) | intercept and centred Gram are both permutation-symmetric | preserved |
| Affine equivariance (IB-5/AP-4) | explicit intercept absorbs translation; the centred kernel is equivariant to scale | preserved |
| Useful-iff (CI-A3): auxiliary information helps iff it changes the joint window | X-CON is precisely a measurement of whether the protein changes the joint window, and the zero-gradient property under constant $G_P$ is its operational form | **this is the new correspondence** |
| $\mathsf A(F,z)=K(B(z)F(z))$, radius, coverage flags | **not connected** | no admitted bounded $z$; coverage remains an engineering surrogate, not the frozen certificate |

Two boundaries restated: the ridge $\lambda$ is not the frozen simplex $\mu$;
kernel power is not a frozen certificate. A point predictor is not the law
operator, and upstream success does not inherit the law theorem.

---

## 11. Remaining scientific risks

1. **Signal-to-noise of $\Delta_\times$.** A second difference of four measured
   cells has up to four times the variance of one cell. Mitigation: restrict to
   same-document rectangles where errors are correlated; aggregate to
   panel-level statistics before inference. This is the single largest risk to
   X-CON and it is not yet measured.
2. **Rectangle dependence.** Rectangles sharing a protein, ligand or document
   are not independent, and the local panels already close to giant components.
   All inference must be at component level; `v2_rectangle_census.py` reports
   the largest-component share and refuses the confirmation label.
3. **Family-level rather than protein-level metric.** $g(P)$ may learn a family
   selector, which the existing cross-family wrong-protein control cannot
   detect. The new within-family donor control (§7) is required and has never
   been run.
4. **The prerequisite may simply fail.** If arm C4 of the transfer ladder does
   not beat C1/C3 under double-cold evaluation, no encoder repair is warranted
   and the correct terminal verdict changes.
5. **No fresh confirmation supply exists.** Even a clean development PASS cannot
   be confirmed with current data. A prospective cohort isolated by protein
   family, scaffold and document, with ≥30 eligible targets and adequate
   independent components, must be built before any partner-specificity claim.
6. **Pocket definition coverage.** KLIFS pockets cover kinases; the gated
   residue-selection fallback for non-kinase targets is unvalidated and could
   silently reintroduce a global-composition statistic.
7. **Absolute performance remains poor.** main-v0 `R2 = −1.244`, Pearson
   `0.097`. Nothing proposed here is yet a competitive DTA claim, and no
   mechanism result should be reported as one.

---

## 12. Terminal verdict

The centred section is retained as the mandatory attribution fix and as the
registered baseline arm, but it is proved insufficient as a partner-specificity
mechanism (§1.1). The partner information loss in the current pair
representation is localized to a specific, repairable aggregation stage (§1.2).
A training-level mechanism that makes correct pairing and correct identity both
necessary is derived, implemented and unit-specified (§9). No new training was
executed this cycle (§6), so no candidate can be advanced on evidence, and the
prerequisite transfer question is unanswered.

```text
BIOLOGICAL_REPRESENTATION_REPAIR_REQUIRED
```

Production `model/` and `scripts/` are unmodified. CSMO, Band and
$\mathsf A(F,z)=K(B(z)F(z))$ are unchanged and disconnected. No Gate was opened
and no production migration is authorized.
