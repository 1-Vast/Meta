# Core Task 1 — Verifiable Conditional Interaction for Cold-Target DTA
## A Complete Research Program with Explicit References

- Date: 2026-08-20
- Kind: research-program / design document (theory + design; **no training launched by this file**)
- Location: `report/research_ideas/`
- Author role: independent senior research scientist
- Status relative to governance: Core Task 1 is currently **UNRESOLVED on current local assets**
  (frozen in `report/CORE_TASK1_UNRESOLVED_TERMINAL_20260817.md`). This document does NOT
  overturn that verdict. It is a proposed **completion program** consistent with it: UNRESOLVED
  is a legal terminal state, and this program defines the smallest set of experiments that could
  move the verdict to SOLVED or to FALSIFIED-AS-TESTED under pre-registered controls.
- Supersedes/refocuses: `report/research_ideas/CORE_TASK1_INTERACTION_PROOF_PROGRAM_20260820.md`
  (earlier draft; this file is the complete version with an explicit reference section).
- Companion design inputs: `report/research_ideas/ciip/CIIP2_RESEARCH_REPORT_20260820.md`
  (OLR-Potential, one candidate mainline), `report/research_ideas/CIIP_SUCCESSOR_STAGE_RESEARCH_PLAN_20260820.md`
  (identifiability prequel).

Inline citations use `[n]`; the full bibliography is in Section 9, grouped by theme, with each
reference's role stated. Repository artifacts cited as grounding are listed separately in Section 10
(they are frozen evidence, not external literature).

---

## 0. Executive summary

**中文摘要。** 核心任务一 = 证明一个可部署模型"真的利用了蛋白—配体条件交互"来做冷靶点 DTA。
本报告把它重定义为一个**可证伪的五柱完成判据集**（CT1-1…CT1-5），并给出 **5 条可行方案**
（模型创新点），每条都满足：(i) 交互项结构上显式、可分离；(ii) 在冷靶点切分下 load-bearing
（消融它会按预注册幅度掉点）；(iii) 由反事实/消融/负控可验证，而非靠 attention 热图；
(iv) 仅用可部署信息（序列 + 合法先验 + 配体结构 + 少量 support 标签，无突变坐标）。
配套一套与架构无关的**反事实交互验证套件 CIVS**（7 项），把"模型用了交互"变成可证伪命题。
推荐主线：**方案 B（可部署配体条件残差路由）+ 方案 D（差分原生训练目标）**；方案 A 作最小可识别
参照、方案 C 作 do()-干预可验证层、方案 E 作 few-shot 桥。所有参考文献在第 9 节逐条列明并说明用途；报告同时与本地历史对齐，第 2b 节逐条列出已被证伪/关闭的分支并从全部方案中排除，作为 CIVS 负控复验。

**English.** Core Task 1 is reframed from "build a better DTA model" to "produce falsifiable
evidence that a protein–ligand conditional interaction term is (a) deployable, (b) load-bearing
on cold targets, and (c) not reducible to protein main effect, ligand main effect, family key,
assay batch, or random-context shortcut." Five literature-grounded, deployable model-innovation
schemes are specified, plus one architecture-agnostic verification suite (CIVS) that is the
methodological core. Completion requires all five evidence pillars to pass under pre-registered,
bootstrap-stable tests; falsification at any pillar is a legitimate terminal outcome. The program is integrated with the local history: Section 2b enumerates every already-falsified / closed branch, excludes it from all schemes, and re-tests each as a CIVS negative control.

---

## 1. Core Task 1 as a five-pillar completion criteria set

A model M completes Core Task 1 only if **all** five pillars hold, each pre-registered,
keyed-rng, parent/family-cluster bootstrapped (2000 draws), with leave-one-family-out sign
stability. Bootstrap means are never used as point estimates; Spearman is undefined (never 0)
for constant predictions.

- **CT1-1 Deployability.** M uses only protein sequence / legal protein priors (ESM residue
  states, KLIFS pocket, licensed conservation) + ligand structure + few support labels. No
  mutation coordinate, no co-complex structure at inference, no test-time label leakage.
- **CT1-2 Cold-target transfer.** On unseen proteins/families, unseen chemotypes, and at least
  one independent dataset, the interaction-driven prediction is stable on MSE/RMSE and CI/
  Spearman/Pearson vs strong baselines {Tanimoto transport, ligand-only floor, additive
  main-effect model, frozen-ESM linear probe}. Both k=0 and k∈{1,2,3,5}.
- **CT1-3 Interaction load-bearing.** Ablating/zeroing the explicit interaction pathway
  degrades cold-target performance by a pre-registered margin (parent/family-cluster bootstrap
  lo2.5 > 0). The model must NEED the interaction term, not merely contain it.
- **CT1-4 Interaction-specificity.** The signal survives matched counterfactual and negative
  controls: same-parent wrong-mutation, family-preserving shuffle, random local window,
  protein-invariant shift, ligand-only floor, assay/family-key permutation. The signal must be
  DESTROYED by interaction-scrambling controls and shift predictions in the LABEL direction
  under ligand-swap / protein-edit counterfactuals.
- **CT1-5 Metric superiority.** Stable gain on MSE/RMSE and CI/Spearman over the baseline set,
  not attributable to a re-fit level term.

Completion = CT1-1 ∧ CT1-2 ∧ CT1-3 ∧ CT1-4 ∧ CT1-5.

### Why the level wall makes Core Task 1 the ONLY path

The frozen boundary (repo artifact `BOUNDARY_20260817_NIGHT.md`) establishes that on the
governed BindingDB-Ki double-cold protocol, "the level term is assay history," and the protocol
removes the transferable part: legal protein probes explain ≤ ~25.9% of level variance; the
pocket prior failed its identifiability gate; LM conditioning gave ranking gains but no MSE
movement. Consequence: on a cold target you CANNOT recover the level by memorizing the protein.
The only transferable object is the **conditional interaction** — how binding changes as a joint
function of protein state and ligand. Therefore proving the model uses conditional interaction is
not a side quest; it is the mechanism by which the k=0 level wall can be broken. Any scheme that
still leans on a memorizable level term will hit the same wall.

---

## 2. Why prior work could not complete Core Task 1

- **CIIP-1A** tested an ORACLE mutation-coordinate ESM window on an entangled centered estimand
  with a miscalibrated random-window negative control and adjudicated
  `ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED`. It proved the estimand was near-unidentifiable, not that
  interaction is absent. It never tested a deployable representation or a cold-target split.
- **The centered estimand trap.** c = parent-level component (achievable ceiling ~Spearman 0.58
  WITHOUT reading the mutation) + mutation-specific residual (~2/3 energy, idiosyncratic, no
  cross-mutation sharing) + panel-shared ligand pattern (~R² 0.13). A single centered target
  mixes all three, so no arm could attribute its score to conditional interaction.
- **CIIP-2 (OLR-Potential)** proposes the right decomposition but is unexecuted and is ONE point
  in design space; it has no within-family parent-profile ceiling arm and no fit-time
  counterfactual measurement.
- **Cold-target DTA literature** (Section 3a) overwhelmingly claims interaction via attention
  heatmaps, binding-site-overlap case studies, or cold-split benchmark metrics — none of which
  is a load-bearing or counterfactual proof. This is exactly the gap CIVS fills.

### 2b. Local-history integration — already-falsified branches are excluded, not re-run

A hard rule of this program: no scheme may re-open a branch the local history has already
closed at the pre-registered standard. Each closed branch below is (i) excluded from all five
schemes by construction, and (ii) re-tested inside CIVS as a negative control, so the program
actively confirms it stays closed rather than silently assuming it. This is what makes the
program cumulative with the local record instead of repetitive of it.

| # | Already-closed branch | Local verdict and source | How this program handles it |
|---|---|---|---|
| 1 | Oracle mutation-coordinate local ESM window | ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED; correct-vs-random R² diff −0.1217, parent bootstrap [−0.4569, +0.0327] (CONTROL_REPORT.md; CIIP-1A) | Scheme B uses full-sequence deployable ESM + ligand routing, no coordinate. The oracle window is re-run as a CIVS-6 negative architecture control. |
| 2 | Mutation-specific interaction V1 at 38 rows | unidentifiable above the noise envelope, θ −0.406 [−0.704, −0.073] (CURRENT_MODEL_EVIDENCE.md; CIIP-1A) | No scheme claims mutation-specific recovery from 38 rows. Attribution is field/interaction-level; the estimand is re-derived under a power table (CT1-0) before any training. |
| 3 | CIIP-1B | forbidden, never authorized (task.md; EVIDENCE_LEDGER.md) | Not proposed; no scheme opens it. |
| 4 | BindingDB Bridge / CIIP potential into BindingDB | forbidden, never authorized (task.md) | Not proposed; no scheme transfers potential into BindingDB. |
| 5 | P0 production mechanism integration | frozen, not authorized (task.md) | Research-only; production model/ and scripts/ untouched. |
| 6 | Pocket prior / Manning family-category prior | FAILED identifiability gate, R² −0.021 (CIIP-2 diagnostics) | Pocket is at most an optional, non-load-bearing field in Scheme A. No scheme makes pocket/family-category the load-bearing path. CIVS-5 family-preserving shuffle + family-key permutation must not reproduce the signal. |
| 7 | Centered single-target estimand | entangled: parent-level ceiling Spearman 0.579 WITHOUT mutation + idiosyncratic residual ~2/3 energy + shared-ligand R² 0.13 (CIIP_SUCCESSOR §4) | Scheme D trains only on pairwise differences (main effects cancel); Schemes A/B use a cross-fitted residual (DML [25]). A single centered target is never re-fit for attribution. |
| 8 | Pure main-effect / additive cold-target level | level wall: level term = assay history; protein probes ≤25.9% level variance; LM conditioning gives no MSE movement (BOUNDARY; CURRENT_MODEL_EVIDENCE.md) | No scheme leans on a memorizable level. The nuisance/level head is separate and audited apart from the interaction head. |
| 9 | Ligand-global / shared-ligand pattern as interaction | low power: LOO R² 0.060; shared-ligand train-mean profile R² 0.1291 (new-session diagnostics) | Ligand-only is a floor / negative control (CIVS-1, CT1-4), never the claim. |
| 10 | Attention heatmaps / binding-site overlap as proof | forbidden evidence (Jain & Wallace [12]; repo constraint) | CIVS reports them DESCRIPTIVE only; proof is ablation + counterfactual direction. |
| 11 | Cross-platform residual sharing with Davis | CLOSED: pairs involving Davis are negative (CORE_TASK1_UNRESOLVED_TERMINAL) | No scheme imports Davis residuals; cold splits use governed BindingDB-Ki double-cold + DTI-DG [10]. |
| 12 | Context-only prediction | NOT evaluated as a mechanism (propagation measured, no predictor fit) | Scheme E actually FITS a support-conditioned predictor and verifies it by support manipulation, converting the unevaluated gap into a tested claim. |
| 13 | Censored single-platform panels | severe censoring: Davis 71.2%, Metz 60.4%, Klaeger 93.5% at detection floors (TERMINAL) | CT1-0 censored-exclusion audit before any label use; functional % inhibition never relabelled Ki/Kd/pK. |
| 14 | Duong-Ly parent-profile as mutation-specific signal | parent-level-dominated: parent-profile LOPO ceiling Spearman 0.579 WITHOUT reading the mutation | Schemes claim interaction from ligand×residue modulation, not parent-profile; parent-profile is a ceiling reference only. |
| 15 | MSA / coevolution features | externally blocked (no UniRef snapshot locally) | Deployability (CT1-1) uses sequence + legal priors only. |
| 16 | meta_test | sealed, never read | Never read; final evaluation only on frozen governed splits. |
| 17 | W1 training | never trained, paused (TERMINAL) | Not re-opened; out of scope. |

**Consequence.** The five schemes occupy ONLY the not-yet-ruled-out region: deployable
sequence+ligand conditional interaction, tested by load-bearing ablation and counterfactuals.
Every closed branch above is re-tested as a negative control so the program confirms it stays
closed. Nothing in Schemes A–E re-litigates a branch the history has already adjudicated.

---

## 3. Literature synthesis and applicability boundaries

Notation: [PR] peer-reviewed, [PP] preprint/theory. "Fits current assets" = usable with the
ESM residue cache + KLIFS + ECFP + Duong-Ly + BindingDB-Ki, no new external data/structures/MSA.

### 3a. How cold-target DTA papers claim interaction (and their proof gap)

| # | paper | mechanism | how interaction is claimed | proof gap CIVS fills |
|---|---|---|---|---|
| [1] | ZeroBind, Nat Commun 2023 [PR] | protein-specific zero-shot DTI via ligand-subgraph ↔ binding-site-subgraph matching | subgraph-motif transfer to unseen proteins; binding-site overlap | binding-site library often structure-derived; interaction is a matching heuristic, not load-bearing tested |
| [2] | Task-conditioned DTI, NeurIPS 2022 | conditions on task | task embedding | conditioning ≠ proof of interaction usage |
| [3] | CoAff-DTI, J Biomed Inform [PR] | PLM + affinity-guided fine-grained interaction | affinity-guided attention | attention heatmaps, no counterfactual |
| [4] | Meta-learning task-adaptive DTA, Nat Commun | meta-learning task adaptation | meta-generalization | task adaptation ≠ interaction attribution |
| [5] | Meta-learning inductive logistic matrix completion (kinase inhibitors) | low-rank interaction transfer | matrix-completion transfer | closest few-shot kinase analogue to Scheme E |
| [6] | Dual-modality binding-site-informed DTA, npj Digit Med 2025 [PR] | fuses binding-site info | site-informed fusion | binding site often structure-derived |
| [7] | Contrastive kinase-inhibitor activity & selectivity, Nat Commun 2025 [PR] | contrastive pretraining | selectivity gap | contrastive signal, no counterfactual load-bearing test |
| [8] | EGA-DTA, target-conditional gating | target-conditional gating | gating = conditional interaction | closest architecture analogue to Scheme B; gating ablation not standardized |
| [9] | Dr Kinase: drug-resistance hotspots | resistance hotspot prediction | hotspot enrichment | biology prior for the mutation-specific flank |
| [10] | TDC / DTI-DG benchmark | domain-generalization splits | cold-domain protocol | provides the split discipline for CT1-2 |
| [11] | CS-DTA (entity-disjoint cold-start) | entity-disjoint eval | cold-target protocol | already internalized by the repo |

**Synthesis of the gap.** Every paper above establishes TRANSFER (a cold-split metric) but none
establishes that the interaction term is LOAD-BEARING or INTERACTION-SPECIFIC. Attention is not
explanation [12], and shortcut learning is the default failure mode [13]. The missing artifact is
a standardized counterfactual + ablation battery — CIVS (Section 6) — which is the core
methodological innovation of this program and is applicable to ANY scheme.

### 3b. Cross-domain mechanism transfer (raw material for the schemes)

| source field | mechanism | canonical source | transfers to |
|---|---|---|---|
| CV / VQA | FiLM conditioning (per-feature γ,β from a conditioning input) | Perez et al., AAAI 2018 [14] | Scheme B: ligand generates γ,β to modulate protein residue states |
| CV / generative | HyperNetworks (one net generates another's weights) | Ha, Dai, Le, ICLR 2017 [15] | Scheme B stronger variant: ligand generates the readout weights |
| CV / VQA | Modular co-attention (residue ↔ atom) | Yu et al., CVPR 2019 (MCAN) [16] | Scheme B/C co-attention readout |
| NLP | ESIM cross-sentence attention | Chen et al., ACL 2017 [17] | residue–ligand cross-encoding |
| NLP | Counterfactually-augmented data | Kaushik, Wallace, Lipton, ACL 2020 [18] | CIVS counterfactual construction + Scheme D negatives |
| NLP | Hypothesis-only / annotation-artifact baselines | Poliak EMNLP 2018; Gururangan NAACL 2018 [19] | the ligand-only and protein-invariant floors |
| RecSys | FM / Field-aware FM (explicit 2nd-order interaction) | Rendle ICDM 2010; Juan RecSys 2016 [20] | Scheme A: exact variance attribution of the interaction term |
| ML / interpretability | Concept Bottleneck Models (intervenability by construction) | Koh et al., ICML 2020 [21] | Scheme C: do()-intervention on the interaction concept |
| ML / explanation | Counterfactual GNN explanations; causal-spurious decoupling for OOD | [22] | counterfactual test design; nuisance decoupling |
| Meta-learning | MAML/ANIL (feature reuse); CNP/NP (support-conditioned prediction) | Finn ICML 2017; Raghu ICLR 2020; Garnelo ICML 2018 [23] | Scheme E: support-conditioned interaction readout |
| Bioactivity ML | ActFound pairwise meta-learning | Feng et al., Nat Mach Intell [24] | Scheme D: difference-native objective |
| Causal estimation | DML / R-learner / X-learner / CFR (nuisance residualization) | Chernozhukov 2018; Nie & Wager 2021; Shalit ICML 2017 [25] | the cross-fitted residual objective in Schemes A/B |
| Genetics | Marginal epistasis test (detect interaction without joint features) | Crawford et al., PLoS Genet 2017 [26] | detection-first fallback if estimation is underpowered |

### 3c. Out of scope (requirements we do not have)

Structure-required ΔΔG methods (PremPLI-class, need co-complex coordinates); MSA-required methods
(DeepSequence/EVE-class, UniRef snapshot absent locally); large external assay corpora for
ActFound-style pretraining. These are boundary references only.

---

## 4. The design principle (the innovation thesis)

An interaction claim is scientific only if the interaction is:

1. **Structurally explicit** — a named, separable term/pathway I(P,L), not an emergent property
   of a black-box concatenation. This gives an ablation handle.
2. **Load-bearing under the cold split** — removing I(P,L) measurably degrades cold-target
   performance beyond the additive main-effect model. This rules out "present but unused."
3. **Counterfactually verifiable** — under ligand-swap and protein-edit interventions the
   prediction moves in the LABEL direction; under interaction-scrambling and nuisance
   transformations it does not. This rules out family-key / assay-batch / random-context shortcuts.
4. **Deployable** — computable from sequence + legal priors + ligand alone.

Every scheme instantiates all four. CIVS (Section 6) operationalizes them identically across
schemes, so the schemes are directly comparable and Core Task 1 has one uniform standard of proof.

---

## 5. Five feasible schemes (model innovation points)

Each scheme lists: mechanism; why it completes Core Task 1; CIVS hooks; literature grounding;
innovation point; risk; feasibility with current assets. All are <2M params on frozen features,
trainable in minutes/seed on CPU-class budgets, matching repo governance (no closed-form, no
ridge, end-to-end gradient training, keyed rng).

### Scheme A — FAIM: Field-Aware Interaction Machine (minimal identifiable reference)

- **Mechanism.** Decompose protein into fields {ESM-global, KLIFS-pocket, optional family/
  conservation} and ligand into fields {ECFP, scaffold, pharmacophore}. Interaction is an
  explicit field-aware bilinear I(P,L)=Σ_{f,g} v_f(P)^T W_{fg} v_g(L). Main effects m_P, m_L are
  separate and orthogonalized (OID-style centering). Target = cross-fitted residual (DML/
  R-learner [25]) so I explains only the nuisance-removed part.
- **Why it completes CT1.** The interaction term has EXACT variance attribution — drop W and you
  have the additive model; the cold-target gap is the interaction contribution by construction.
  It is the smallest model for which CT1-3/CT1-4 are cleanly measurable, so it is the mandatory
  identifiability reference against which B–E are judged.
- **CIVS hooks.** Ablation = zero W_{fg}; counterfactual ligand-swap changes within-protein
  ranking in the label direction; field-pair shuffle destroys it; ligand-only / protein-invariant
  floors are the negative controls [19].
- **Grounding.** FM/FFM [20]; DML/R-learner [25].
- **Innovation point.** The first DTA model with per-field, exactly-attributable interaction
  terms and a pre-registered load-bearing gate.
- **Risk.** Limited expressivity for nonlinear interaction; may underfit. Acceptable — its job is
  attribution, not maximal accuracy. **Feasibility: high** (existing features, tiny param count).

### Scheme B — LCR: Ligand-Conditioned Residue routing (deployable FiLM/hypernetwork)

- **Mechanism.** Frozen ESM-2 RESIDUE states h_i (full sequence, deployable, NO coordinate). A
  small ligand-conditioned generator emits per-residue/per-region modulation (γ_i(L),β_i(L))
  (FiLM [14]) or gates g_i(L) (target-conditional gating [8]); modulated states
  h'_i=γ_i⊙h_i+β_i are pooled for the interaction readout; a parallel additive path bypasses
  modulation. This is the deployable cousin of CIIP-1A's oracle window: the ligand, not an
  annotation, decides where the protein representation is read.
- **Why it completes CT1.** CT1-1 (sequence+ligand only); CT1-2 (routes on cold families because
  routing depends on ligand × residue state, not a memorized target); CT1-3 (γ→1,β→0 collapses to
  additive; the cold-target gap is the interaction); CT1-4 (counterfactual ligand swap shifts
  within-protein ranking; ligand-conditioned routing on a cold protein still shifts predictions).
- **CIVS hooks.** Modulation ablation; ligand-swap counterfactual; routing-weight vs known
  binding sites is DESCRIPTIVE ONLY (never proof, per Jain & Wallace [12]); shuffled-generator
  negative control.
- **Grounding.** FiLM [14]; HyperNetworks [15]; MCAN/ESIM co-attention [16][17]; EGA-DTA gating
  [8]; CIIP-2 LCRR.
- **Innovation point.** Deployable ligand-routed residue field — the interaction is conditional
  modulation of the protein representation by the ligand, requiring no mutation coordinate; it
  directly repairs CIIP-1A's failure mode (oracle dependence) while keeping residue-level richness.
- **Risk.** Routing may collapse to a global (family-key) shift. CIVS nuisance-invariance +
  family-preserving shuffle detect this; a per-region routing-entropy regularizer can be
  pre-registered as a mitigation. **Feasibility: high** (ESM residue cache exists; generator tiny).

### Scheme C — ICB: Interaction Concept Bottleneck (verifiability by construction)

- **Mechanism.** Force prediction through a low-dimensional interaction concept z=I(P,L)∈R^d
  (d≈8–16). z is trained to predict the within-protein ligand contrast (centered);
  ŷ=b_P+b_L+h(z). Concepts may be weakly supervised (KLIFS pocket class × ligand class
  posteriors) but this is optional.
- **Why it completes CT1.** A concept bottleneck is INTERVENABLE by construction (Koh et al.
  [21]): you can run do(z→z') and check that the predicted affinity moves as the concept
  semantics demand. This is the most direct answer to "does the model use the interaction?" —
  you literally intervene on the interaction and watch the output. CT1-3 by ablating z; CT1-4 by
  do()-counterfactuals.
- **CIVS hooks.** do(z) intervention direction test; z-ablation cold-target collapse; concept
  interpretability; ligand-swap changes z in the label direction.
- **Grounding.** Concept Bottleneck Models [21]; counterfactual explanations [22].
- **Innovation point.** The first DTA model whose interaction usage is provable by do()-
  intervention rather than post-hoc attribution.
- **Risk.** Bottleneck too tight → underfit; weak-supervision concept labels noisy. Mitigate with
  a pre-registered capacity sweep on val only. **Feasibility: medium-high.**

### Scheme D — DNT: Difference-Native Training (identifiability-by-design; the training-module innovation)

- **Mechanism.** NEVER train on absolute affinity. Train only on pairwise differences: within-
  protein ligand-ranking, within-ligand protein-ranking, and cross contrasts (ActFound-style
  pairwise meta-learning [24]). Main effects cancel by construction, so ANY predictive success is
  definitionally interaction-driven — there is no main-effect pathway to leak into. For deployment,
  a separate, separately-audited nuisance/level head restores absolute scale and is reported apart
  from the interaction head.
- **Why it completes CT1.** Identifiability BY DESIGN — the strongest answer to CT1-4, because the
  training objective makes a main-effect shortcut impossible rather than merely controlled.
  Within-panel contrasts are the estimand shown robust to assay variability and isomorphic to
  ActFound's deployed design [24].
- **CIVS hooks.** Identity-zero antisymmetry (bitwise); support/query label isolation; cold-target
  within-protein Spearman as headline; counterfactual negative pairs (same-parent wrong-mutation)
  as decisive negatives [18].
- **Grounding.** ActFound [24]; pairwise/contrast training; the repo's own pairwise arms.
- **Innovation point.** A training module (satisfies the "at least one innovation in the training
  module" rule) that makes interaction the ONLY learnable signal; pairs with any of A/B/C as scorer.
- **Risk.** Absolute MSE needs the separate level head (which hits the level wall — acknowledged;
  the interaction head is judged on ranking/contrast, the level head audited separately). Pairwise
  sparsity on cold targets; mitigate with k-shot support (couples to Scheme E).
  **Feasibility: high** (matches existing pairwise/contrast infrastructure).

### Scheme E — SCNPI: Support-Conditioned Neural Process Interaction (few-shot bridge)

- **Mechanism.** For a cold target, a k-shot support set S={(L_i,y_i)} is encoded by a
  permutation-invariant set encoder (DeepSets/NP [23]) into a context c(S); the query prediction
  for L_q conditions the interaction readout on c(S). The set encoder is the transfer mechanism:
  it carries the target's conditional-interaction fingerprint from a few labelled ligands to new
  ones.
- **Why it completes CT1.** CT1-2 across the k∈{0,1,2,3,5} curve; CT1-3 by support-ablation
  (empty support → k=0 prior); CT1-4 by support-label isolation (shuffle support labels → k-shot
  gain destroyed), support permutation invariance, query-label isolation, and wrong-protein-
  support counterfactual (must destroy predictions).
- **CIVS hooks.** k-shot improvement curve; support-label isolation test; permutation invariance
  to ~1e-6; wrong-support counterfactual.
- **Grounding.** CNP/NP [23]; MAML/ANIL [23]; meta-learning kinase matrix completion [5]; the
  repo's P-line CNP re-adjudication.
- **Innovation point.** An explicit, manipulable few-shot transfer of the conditional interaction;
  the interaction is read out conditioned on support, so its usage is verifiable by support
  manipulation.
- **Risk.** Sparse/noisy support on cold targets; NP variance. Mitigate with the deterministic
  DeepSets variant already re-adjudicated in the repo. **Feasibility: high** (P-line infra exists).

---

## 6. CIVS — Counterfactual Interaction Verification Suite (the proof standard)

CIVS is architecture-agnostic and is run identically for every scheme. It is the methodological
core of this program: it turns "the model uses interaction" into a falsifiable claim. All tests
are pre-registered, keyed-rng, parent/family-cluster bootstrapped (2000 draws), with
leave-one-family-out sign stability; bootstrap means never used as point estimates; Spearman
undefined (never 0) for constant predictions. Counterfactual construction follows [18]; nuisance
decoupling follows [25][22]; negative baselines follow [19].

1. **Cold-target transfer** (CT1-2): double-cold + DTI-DG-style splits [10][11]; MSE/RMSE, CI,
   Spearman, Pearson vs {Tanimoto transport, ligand-only, additive main-effect, frozen-ESM
   linear}; k∈{0,1,2,3,5}.
2. **Load-bearing ablation** (CT1-3): zero the interaction pathway; cold-target gap with
   bootstrap lo2.5 > 0 across families.
3. **Counterfactual ligand-swap** (CT1-4): fix a cold protein, swap ligand features at inference;
   within-protein ranking must move in the label direction (sign accuracy / Spearman of
   counterfactual shift vs true shift above a frozen threshold).
4. **Counterfactual protein-edit** (CT1-4): fix a ligand, edit the protein (deployable: ESM on
   the edited sequence / pocket perturbation); within-ligand ranking shifts in the label direction.
5. **Nuisance-invariance** (CT1-4): predictions invariant (within tolerance) to assay/family-key/
   batch permutations and to random-context (non-informative) protein perturbations; random-window
   / shuffled-pocket controls must NOT reproduce the interaction signal.
6. **Interaction-scrambling negative architecture** (CT1-4): same architecture with the
   interaction pathway shuffled/randomized at train time must collapse to baseline — rules out the
   architecture accidentally encoding main effects.
7. **Identifiability gate** (CT1-4/CT1-5): the interaction effect survives family-cluster
   bootstrap and leave-one-family-out; a pocket/family-key identifiability probe (as in the frozen
   boundary) must not absorb it.

A scheme completes Core Task 1 only if CIVS 1–7 pass at pre-registered thresholds. Attention
heatmaps / binding-site overlap may be reported as DESCRIPTIVE supplements only, never as evidence
for any pillar (Jain & Wallace [12]).

---

## 7. Recommended main line, staging, and promotion gates

### Recommended main line

**Scheme B (deployable ligand-conditioned residue routing) as the mechanism + Scheme D
(difference-native objective) as the training module**, with **Scheme A** as the mandatory
minimal identifiable reference and **Scheme C** as the verifiability/do()-intervention layer on
the same backbone; **Scheme E** reserved for the few-shot (k≥1) bridge. This combination:
satisfies "at least one innovation in the training module" (D) and one in the architecture (B);
is fully deployable (no coordinate); is falsifiable component-by-component (each is a single
toggle); and directly repairs the two audited failure modes (oracle dependence; entangled centered
estimand). It is compatible with, and a strict generalization of, CIIP-2's OLR-Potential
(B≈LCRR, D≈CFOIE residual). If governance runs CIIP-2 first, this program's CIVS battery and
Schemes A/C must be prepended by addendum.

### Staged execution (all stages preregistered + SHA-256; CPU smoke → single seed → multi-seed
### only after structure + negatives pass)

- **CT1-0 Governance & data audit** (read-only): freeze inputs (ESM residue cache, KLIFS, ECFP,
  Duong-Ly, BindingDB-Ki governed splits); coverage / parent-overlap / ligand-overlap /
  assay-semantics / censoring audit; leakage audit; re-derive the CIIP_SUCCESSOR §4 decomposition
  under preregistration; power table + MDE per CIVS test; freeze thresholds by dated addendum
  BEFORE any training. meta_test never read; test labels only for final evaluation of frozen arms.
- **CT1-1 Scheme A (FAIM)**: smallest identifiable interaction; run full CIVS. Establishes the
  attribution reference and validates the CIVS pipeline end-to-end.
- **CT1-2 Scheme B (LCR)**: deployable routing; run full CIVS; compare to A.
- **CT1-3 Scheme D objective on B**: difference-native training; CIVS; compare to the centered-
  residual objective (a univariate toggle).
- **CT1-4 Scheme C (ICB) do()-intervention** on the B backbone: the verifiability layer; CIVS +
  concept-intervention tests.
- **CT1-5 Scheme E (SCNPI)**: k-shot bridge; CIVS k-curve + support tests.

### Promotion gates and stop rules

A scheme advances only if CIVS 1–7 pass at frozen thresholds on a single seed AND reproduce on
≥3 seeds with family-cluster bootstrap lo2.5 > 0 and LOFO sign stability. No promotion on a
single seed/parent/few pairs; no metric shopping; UNRESOLVED is a legal terminal state. Any scheme
that collapses to the ligand-only / additive / family-key floor under CIVS-6, or whose
counterfactual direction test (CIVS-3/4) fails, is falsified-as-tested and closed with a boundary
entry. If ALL schemes fail CIVS-2 (cold-target transfer), Core Task 1 is recorded as falsified
under current legal data and deployable constraints and the boundary document is updated — a
legitimate terminal outcome.

---

## 8. Success criteria, deliverables, and constraints

### Success (Core Task 1 complete)

A promoted scheme passes CIVS 1–7 with: cold-target ranking/contrast gain over the strong baseline
set (bootstrap-stable, LOFO-stable); a significant load-bearing ablation gap; counterfactual
direction tests above frozen thresholds; nuisance-invariance within tolerance; and metric
superiority not reducible to a re-fit level term. Completion is reported per-pillar (SUPPORTED /
NOT SUPPORTED / UNRESOLVED), never as a single number.

### Deliverables per stage

PREREGISTRATION.md + SHA-256; threshold addendum (frozen pre-training); data-audit JSON/MD;
RESULT.json (machine-readable all metrics); REPORT.md (civil verdict table + authorization block);
commands.jsonl; structure + data-contract tests; SHA256SUMS; FAILURES.md; append-only sync to
history.md / task.md / EVIDENCE_LEDGER.md.

### Production constraints (binding)

No change to model/ or production scripts/ during research; no oracle mutation-coordinate ESM in
any deployment path; no CIIP potential into BindingDB, no CIIP-1B, no BindingDB Bridge start
without its own assay-semantics qualification; no biological-mechanism claim from correct-vs-
random-window differences; no success claim from single seed/parent/few pairs; no larger
backbone/budget to mask non-identifiability; functional % inhibition never relabelled Ki/Kd/pK;
context-propagation magnitude never cited as predictive value; attention heatmaps never cited as
proof of interaction usage; failure never reported as biological absence of protein-conditioned
signal.

---

## 9. External references (bibliography, by theme)

Every inline citation `[n]` above resolves here. Each entry states the venue/year/identifier and
its role in this report. Where an author list is not independently verified, the entry cites
title + venue + year + URL (all independently checkable); no authorship is asserted beyond what
was verified.

### 9.1 Cold-target / zero-shot DTA (the target-domain literature)

- **[1] ZeroBind.** "ZeroBind: a protein-specific zero-shot predictor with subgraph matching for
  drug-target interactions." Nature Communications, 2023.
  https://www.nature.com/articles/s41467-023-43597-1
  Role: the flagship zero-shot cold-target DTI method; motivates Scheme B's ligand-routed
  subgraph/residue matching, and exemplifies the "transfer by matching heuristic, not load-bearing
  proof" gap CIVS fills. (SJTU Pan Xiaoyong group; verified via SJTU SEIEE news.)
- **[2] Task-conditioned DTI.** NeurIPS 2022 paper on task-conditioned drug–target interaction.
  https://dev.neurips.cc/virtual/2022/57454
  Role: shows conditioning on task ≠ evidence of interaction usage (CT1-4 motivation).
- **[3] CoAff-DTI.** "CoAff-DTI: Fine-grained drug–target interaction prediction using pre-trained
  language models and affinity-guided mechanisms." Journal of Biomedical Informatics, 2026.
  https://www.sciencedirect.com/science/article/abs/pii/S1532046426001000
  Role: PLM + affinity-guided fine-grained interaction; its interaction evidence is attention
  heatmaps (proof gap). Authors verified as Peng, Liu, et al. via Semantic Scholar.
- **[4] Meta-learning task-adaptive DTA.** "A meta learning and task adaptive approach for drug
  target affinity prediction." Nature Communications.
  https://www.nature.com/articles/s41467-026-70554-5
  Role: meta-generalization ≠ interaction attribution; informs Scheme E and CT1-2 protocol.
- **[5] Meta-learning inductive logistic matrix completion for kinase inhibitors.**
  https://ouci.dntb.gov.ua/en/works/4bggWQNl/
  Role: closest few-shot kinase analogue to Scheme E (low-rank interaction transfer).
- **[6] Dual-modality binding-site-informed DTA.** npj Digital Medicine, 2025.
  https://www.nature.com/articles/s41746-025-01464-x
  Role: binding-site-informed fusion; flags that binding site is often structure-derived, which
  our deployability constraint (CT1-1) forbids.
- **[7] Contrastive kinase-inhibitor activity & selectivity.** "Enhancing kinase-inhibitor activity
  and selectivity prediction through contrastive learning." Nature Communications, 2025.
  https://doi.org/10.1038/s41467-025-65869-8
  Role: contrastive signal for kinase selectivity; motivates Scheme D's contrast objective; its
  interaction claim lacks a counterfactual load-bearing test. (Lead author verified as Li Yu-Quan.)
- **[8] EGA-DTA.** "EGA-DTA: An Energetic-Geometric Augmented Graph Neural Network With
  Target-Conditional Gating for DTA Prediction."
  https://visualize.jove.com/42599860
  Role: closest architecture analogue to Scheme B (target-conditional gating); its gating ablation
  is not standardized — a gap CIVS standardizes.
- **[9] Dr Kinase.** Drug-resistance hotspot prediction resource.
  https://digitalcommons.library.tmc.edu/uthshis_docs/678/
  Role: biology prior for the mutation-specific flank (kinase resistance hotspots), used only as
  descriptive context, never as a feature at inference.
- **[10] TDC / DTI-DG benchmark.** "Therapeutics Data Commons / Drug-Target Interaction Domain
  Generalization." arXiv:2102.09548. https://arxiv.org/abs/2102.09548
  Role: provides the cold-domain split discipline for CT1-2.
- **[11] CS-DTA (entity-disjoint cold-start).** PMC13161074.
  Role: entity-disjoint evaluation protocol already internalized by the repo.

### 9.2 Explanation faithfulness and shortcut learning (why CIVS exists)

- **[12] Attention is not Explanation.** Jain, S. & Wallace, B.C. NAACL 2019. arXiv:1902.10186.
  https://aclanthology.org/N19-1357/
  Role: the core reason attention heatmaps / binding-site overlap are DESCRIPTIVE only, never
  evidence of interaction usage.
- **[13] Shortcut learning in deep neural networks.** Geirhos, R. et al. Nature Machine
  Intelligence, 2020. https://www.nature.com/articles/s42256-020-00257-z
  Role: motivates the negative-control / counterfactual design of CIVS-5/6.

### 9.3 Cross-domain mechanisms (CV / NLP / RecSys / meta-learning / causal)

- **[14] FiLM.** "FiLM: Visual Reasoning with a General Conditioning Layer." Perez, E. et al.
  AAAI 2018. https://ojs.aaai.org/index.php/AAAI/article/download/11671/11530
  Role: the conditioning primitive (γ,β) used by Scheme B.
- **[15] HyperNetworks.** Ha, D., Dai, A., Le, Q.V. ICLR 2017. https://arxiv.org/abs/1609.09106
  Role: stronger Scheme B variant where the ligand generates readout weights.
- **[16] MCAN.** "Deep Modular Co-Attention Networks for Visual Question Answering." Yu, Z. et al.
  CVPR 2019. https://openaccess.thecvf.com/content_CVPR_2019/html/Yu_Deep_Modular_Co-Attention_Networks_for_Visual_Question_Answering_CVPR_2019_paper.html
  Role: modular co-attention (residue ↔ atom) for Schemes B/C.
- **[17] ESIM.** "Enhanced LSTM for Natural Language Inference." Chen, Q. et al. ACL 2017.
  https://aclanthology.org/P17-1152/
  Role: cross-encoding attention pattern for residue–ligand readout.
- **[18] Counterfactually-Augmented Data.** "Learning the Difference that Makes a Difference with
  Counterfactually-Augmented Data." Kaushik, D., Wallace, B.C., Lipton, Z.C. ACL 2020.
  arXiv:1909.12434. https://aclanthology.org/2020.acl-main.711/
  Role: the blueprint for constructing CIVS counterfactuals and Scheme D decisive negatives.
- **[19] Hypothesis-only / annotation-artifact baselines.** Poliak, A. et al. EMNLP 2018
  (https://aclanthology.org/D18-1024/); Gururangan, S. et al. NAACL 2018
  (https://aclanthology.org/N18-2017/).
  Role: justify the ligand-only and protein-invariant floors as negative controls.
- **[20] FM / Field-aware FM.** Rendle, S. "Factorization Machines." ICDM 2010
  (https://researchr.org/publication/Rendle10/bibtex); Juan, Y. et al. "Field-aware Factorization
  Machines for CTR Prediction." RecSys 2016. https://dl.acm.org/doi/10.1145/2959100.2959134
  Role: explicit 2nd-order interaction with exact variance attribution → Scheme A.
- **[21] Concept Bottleneck Models.** Koh, P.W. et al. ICML 2020.
  https://proceedings.mlr.press/v119/koh20a.html
  Role: intervenability-by-construction → Scheme C's do()-intervention proof.
- **[22] Counterfactual explanations & causal-spurious decoupling.** Counterfactual GNN
  explanations: https://axi.lims.ac.uk/paper/2410.15165 ; causal–spurious feature decoupling for
  OOD: https://www.sciencedirect.com/science/article/abs/pii/S0957417426022633
  Role: counterfactual test design and nuisance decoupling for CIVS and Schemes A/B.
- **[23] Meta-learning & conditional processes.** MAML: Finn, C. et al. ICML 2017; ANIL: Raghu, A.
  et al. ICLR 2020 (https://openreview.net/forum?id=rkgMkCEtPB); CNP/NP: Garnelo, M. et al.
  "Conditional Neural Processes." ICML 2018 (http://proceedings.mlr.press/v80/garnelo18a.html).
  Role: support-conditioned prediction and feature reuse → Scheme E.
- **[24] ActFound.** Feng, B. et al. "ActFound: an evidence-based foundation model for drug-target
  affinity prediction via pairwise meta-learning." Nature Machine Intelligence. DOI
  10.1038/s42256-022-00581-6. GitHub: https://github.com/BFeng14/ActFound
  Role: the deployed pairwise-meta-learning design that Scheme D's difference-native objective is
  isomorphic to.
- **[25] Causal double/debiased estimation.** DML: Chernozhukov, V. et al. Econometrics Journal,
  2018; R-learner/X-learner: Nie, X. & Wager, S. Biometrika 108(2), 2021
  (https://academic.oup.com/biomet/article/108/2/299/5911092); CFR: Shalit, U. et al. ICML 2017
  (http://proceedings.mlr.press/v70/shalit17a.html).
  Role: nuisance residualization / cross-fitted residual objective in Schemes A/B and CIVS-7.
- **[26] Marginal epistasis test.** Crawford, L. et al. PLoS Genetics, 2017.
  https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006869
  Role: detection-first fallback if joint-feature estimation is underpowered.

---

## 10. Repository artifacts cited as grounding (frozen evidence, not external literature)

These are internal frozen records. They are read-only inputs; none is modified by this document.

- `report/CORE_TASK1_UNRESOLVED_TERMINAL_20260817.md` (+ `.json`) — Core Task 1 terminal
  verdict = UNRESOLVED on current local assets; the state this program proposes to move.
- `report/BOUNDARY_20260817_NIGHT.md` — frozen cold-target boundary: level wall, assay-history
  level term, ≤25.9% level-variance by legal probes, pocket-prior identifiability failure.
- `tools/research/stageCIIP_potential_bridge/CONTROL_REPORT.md` — `ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED`.
- `tools/research/stageCIIP_context_propagation_20260820/CONTEXT_PROPAGATION_REPORT.md` —
  representation-level context propagation; context-only prediction not evaluated.
- `report/research_ideas/ciip/CIIP2_RESEARCH_REPORT_20260820.md` — OLR-Potential candidate.
- `report/research_ideas/CIIP_SUCCESSOR_STAGE_RESEARCH_PLAN_20260820.md` — identifiability prequel
  and the four-gate CT1 contract this report consolidates.
- `report/research_ideas/CORE_TASK1_INTERACTION_PROOF_PROGRAM_20260820.md` — earlier draft this
  complete report supersedes.
- `report/CURRENT_MODEL_EVIDENCE.md`, `report/EVIDENCE_LEDGER.md` — current baseline and
  chronological claim provenance.

---

## 11. What this report changes relative to the prior documents

1. Refocuses the entire programme on Core Task 1 as the SOLE objective and turns it into a
   five-pillar, falsifiable completion criteria set.
2. Adds the missing artifact in both the repo and the cold-target DTA literature: CIVS, a
   standardized counterfactual + load-bearing + nuisance-invariance battery that makes "the model
   uses interaction" provable rather than asserted.
3. Supplies five concrete, deployable, literature-grounded model-innovation schemes (A–E) with a
   recommended combination, instead of a single point design.
4. Explicitly couples Core Task 1 to the frozen level wall: conditional interaction is identified
   as the ONLY transferable object on a cold target, so proving its use is the mechanism for
   breaking the wall.
5. Keeps every prior safeguard: preregistration + SHA-256, sealed meta_test, no test-label
   leakage, no closed-form, single-seed-then-multi-seed, UNRESOLVED as a legal terminal state.
6. Integrates the local history: Section 2b enumerates every already-falsified / closed branch (oracle-coordinate ESM, CIIP-1B, BindingDB bridge, P0 production, pocket/family prior, centered estimand, main-effect level, ligand-global pattern, Davis residuals, censored panels, MSA, meta_test, W1) and excludes each from all schemes, re-testing it as a CIVS negative control.
