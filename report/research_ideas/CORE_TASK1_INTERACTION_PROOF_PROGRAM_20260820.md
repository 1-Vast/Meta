# CORE TASK 1 COMPLETION PROGRAM
## Verifiable Conditional Interaction for Cold-Target DTA

Date: 2026-08-20
Role: independent senior research scientist (design only; no production change, no training launched by this document)
Supersedes/refocuses: CIIP_SUCCESSOR_STAGE_RESEARCH_PLAN_20260820.md (kept as the identifiability prequel)
Companion: report/research_ideas/ciip/CIIP2_RESEARCH_REPORT_20260820.md (OLR-Potential, one candidate)
Status: DESIGN. Execution requires freezing a preregistration per Section 7 and recording SHA-256.

---

## 0. Executive summary

**中文摘要**：核心任务一 = 证明模型"真的利用了蛋白—配体条件交互"，用于冷靶点 DTA。
本报告把它重定义为一个**可证伪的完成判据集**（5 条证据柱），并给出 **5 条可行方案**
（模型创新点），每条都满足：(i) 交互项结构上显式、(ii) 在冷靶点切分下 load-bearing、
(iii) 由反事实/消融/负控可验证、(iv) 仅用可部署信息。配套一套与架构无关的
**反事实交互验证套件 CIVS**，这是"模型真的用了交互"这一断言变成可证伪命题的关键。
推荐主线：**可部署配体条件残差路由（方案 B）+ 差分原生训练目标（方案 D）**，以
字段感知双线性（方案 A）作最小可识别参照、交互概念瓶颈（方案 C）作可干预性验证、
support-条件神经过程（方案 E）作 few-shot 桥。所有方案共享 CIVS；任一方案的"交互
贡献"都由消融落差 + 反事实方向性 + 负控摧毁三项联合裁决，而非 attention 热图。

Core Task 1 is reframed from "build a better DTA model" to "produce
falsifiable evidence that the protein-ligand conditional interaction term is
(a) deployable, (b) load-bearing on cold targets, and (c) not reducible to
protein main effect, ligand main effect, family key, assay batch, or
random-context shortcut." Prior work could not complete this because it
either tested an oracle representation (CIIP-1A), tested an entangled
estimand (the centered contrast dominated by the parent-level component), or
claimed interaction via attention/benchmark metrics without a load-bearing
or counterfactual test. This program closes that gap.

---

## 1. Core Task 1, redefined as a completion criteria set

A model M completes Core Task 1 only if ALL five pillars hold, each with a
pre-registered, bootstrap-stable test:

- **CT1-1 Deployability**: M uses only protein sequence / legal protein
  priors (ESM residue states, KLIFS pocket, conservation if licensed) +
  ligand structure + few support labels. No mutation coordinate, no
  co-complex structure at inference, no test-time label leakage.
- **CT1-2 Cold-target transfer**: on unseen proteins/families, unseen
  chemotypes, and at least one independent dataset, M's interaction-driven
  prediction is stable (MSE/RMSE, CI, Spearman, Pearson vs strong baselines:
  Tanimoto transport, ligand-only floor, additive main-effect model,
  frozen-ESM linear probe). Zero-shot (k=0) AND few-shot (k in {1,2,3,5}).
- **CT1-3 Interaction load-bearing**: ablating/zeroing the explicit
  interaction pathway degrades cold-target performance by a pre-registered
  margin (parent/family-cluster bootstrap lo2.5 > 0). The model must NEED
  the interaction term, not just contain it.
- **CT1-4 Interaction-specificity (counterfactual)**: the interaction signal
  survives matched counterfactual and negative controls — same-parent
  wrong-mutation, family-preserving shuffle, random local window,
  protein-invariant shift, ligand-only floor, assay/family-key permutation.
  The signal must be DESTROYED by interaction-scrambling controls and must
  shift predictions in the label direction under ligand-swap / protein-edit
  counterfactuals.
- **CT1-5 Metric superiority**: stable gain on MSE/RMSE and CI/Spearman over
  the strong baseline set, not attributable to a re-fit level term.

Completion = CT1-1 ∧ CT1-2 ∧ CT1-3 ∧ CT1-4 ∧ CT1-5. Falsification at any
pillar is a legitimate terminal outcome and must be recorded with the exact
estimand, split, and controls.

### Why the level wall makes Core Task 1 the ONLY path

Frozen boundary (report/BOUNDARY_20260817_NIGHT.md): on the governed
BindingDB-Ki double-cold protocol the "level term is assay history," and the
protocol removes the transferable part. Legal protein probes explain at most
~25.9% of level variance; the pocket prior failed its identifiability gate;
LM conditioning gave ranking gains but no MSE movement. Consequence: on a
cold target you CANNOT recover the level by memorizing the protein. The only
transferable object is the CONDITIONAL interaction — how binding changes as a
joint function of protein state and ligand. Therefore Core Task 1 (prove the
model uses conditional interaction) is not a side quest; it is the mechanism
by which the k=0 level wall can be broken. Any scheme that still leans on a
memorizable level term will hit the same wall.

---

## 2. Why prior work could not complete Core Task 1

- **CIIP-1A** tested an ORACLE mutation-coordinate ESM window on an entangled
  centered estimand with a miscalibrated random-window negative control; it
  adjudicated ORACLE_LOCAL_SIGNAL_NOT_SUPPORTED. It proved the estimand was
  near-unidentifiable, not that interaction is absent. It never tested a
  deployable representation or a cold-target split.
- **The centered estimand trap** (CIIP_SUCCESSOR diagnostics): c = parent-
  level component (achievable ceiling ~Spearman 0.58 WITHOUT reading the
  mutation) + mutation-specific residual (~2/3 energy, idiosyncratic, no
  cross-mutation sharing) + panel-shared ligand pattern (~R2 0.13). A single
  centered target mixes all three, so no arm could attribute its score to
  conditional interaction.
- **CIIP-2 (OLR-Potential)** proposes the right decomposition (LCRR + OID +
  CFOIE + assay-gain + optional PMSTD + reserved SCNP) but is unexecuted and
  is ONE point in design space; it has no proposition-A lane, no
  within-family parent-profile ceiling arm, and its C-erased is an
  evaluation-time control rather than a fit-time B-measurement.
- **Cold-target DTA literature** (Section 3a) overwhelmingly claims
  interaction via attention heatmaps, binding-site overlap case studies, or
  benchmark metrics on a cold split — none of which is a load-bearing or
  counterfactual proof. This is the exact gap this program fills.

---

## 3. Literature synthesis and applicability boundaries

Notation: [PR] peer-reviewed, [PP] preprint/theory. "Fits current assets" =
usable with ESM residue cache + KLIFS + ECFP + Duong-Ly + BindingDB-Ki, no
new external data/structures/MSA.

### 3a. How cold-target DTA papers claim interaction (and their proof gaps)

| paper | mechanism | how interaction is claimed | proof gap this program fills |
|---|---|---|---|
| ZeroBind, Nat Commun 2023 [PR] https://www.nature.com/articles/s41467-023-43597-1 | protein-specific zero-shot DTI via ligand-subgraph ↔ binding-site-subgraph matching | subgraph-motif transfer to unseen proteins; binding-site overlap | needs a structure-derived binding-site library; not purely sequence-deployable; interaction is by matching heuristic, not load-bearing test |
| Task-conditioned DTI, NeurIPS 2022 https://dev.neurips.cc/virtual/2022/57454 | conditions on task | task embedding | conditioning ≠ proof of interaction usage |
| CoAff-DTI, Expert Syst Appl [PR] https://www.sciencedirect.com/science/article/abs/pii/S1532046426001000 | PLM + affinity-guided fine-grained interaction | affinity-guided attention | attention heatmaps, no counterfactual |
| Meta-learning & task-adaptive DTA, Nat Commun https://www.nature.com/articles/s41467-026-70554-5 | meta-learning task adaptation | meta-generalization | task adaptation ≠ interaction attribution |
| Meta-learning inductive logistic matrix completion for kinase inhibitors https://ouci.dntb.gov.ua/en/works/4bggWQNl/ | matrix completion + meta-learning | low-rank interaction transfer | closest few-shot kinase analogue to Scheme E |
| Dual-modality binding-site-informed DTA, npj Digit Med 2025 https://www.nature.com/articles/s41746-025-01464-x | fuses binding-site info | site-informed fusion | binding site often structure-derived |
| Enhancing kinase-inhibitor activity & selectivity via contrastive learning, Nat Commun 2025 [PR] https://doi.org/10.1038/s41467-025-65869-8 | contrastive pretraining | selectivity gap | contrastive signal, no counterfactual load-bearing test |
| EGA-DTA: target-conditional gating for DTA https://visualize.jove.com/42599860 | target-conditional gating | gating = conditional interaction | closest architecture analogue to Scheme B; gating ablation not standardized |
| Dr Kinase: drug-resistance hotspots https://digitalcommons.library.tmc.edu/uthshis_docs/678/ | resistance hotspot prediction | hotspot enrichment | biology prior for mutation-specific flank |
| TDC / DTI-DG benchmark, arXiv 2102.09548 [PR-adjacent] https://arxiv.org/abs/2102.09548 | domain-generalization splits | cold-domain protocol | provides the split discipline for CT1-2 |
| CS-DTA (entity-disjoint cold-start), PMC13161074 | entity-disjoint eval | cold-target protocol | already internalized by the repo |

**Synthesis of the gap**: every cold-target DTA paper above establishes
TRANSFER (a cold-split metric) but none establishes that the interaction term
is LOAD-BEARING or INTERACTION-SPECIFIC. Attention is not explanation
(Jain & Wallace, NAACL 2019, https://aclanthology.org/N19-1357/), and
shortcut learning is the default failure mode (Geirhos et al., Nat Mach
Intell 2020, https://www.nature.com/articles/s42256-020-00257-z). The
missing artifact is a standardized counterfactual + ablation battery — that
is CIVS (Section 6), and it is the core methodological innovation of this
program, applicable to ANY scheme.

### 3b. Cross-domain mechanism transfer (the raw material for the schemes)

| source field | mechanism | canonical source | transfers to |
|---|---|---|---|
| CV / VQA | FiLM conditioning (per-feature γ,β from a conditioning input) | Perez et al., AAAI 2018, https://ojs.aaai.org/index.php/AAAI/article/download/11671/11530 | Scheme B: ligand generates γ,β to modulate protein residue states |
| CV / generative | HyperNetworks (one net generates another's weights) | Ha, Dai, Le, ICLR 2017, https://arxiv.org/abs/1609.09106 | Scheme B stronger variant: ligand generates the interaction readout weights |
| CV / VQA | Modular co-attention (residue ↔ atom cross attention) | Yu et al., CVPR 2019 (MCAN), https://openaccess.thecvf.com/content_CVPR_2019/html/Yu_Deep_Modular_Co-Attention_Networks_for_Visual_Question_Answering_CVPR_2019_paper.html | Scheme B/C co-attention readout |
| NLP | ESIM cross-sentence attention | Chen et al., ACL 2017 (TACL/ACL), https://aclanthology.org/P17-1152/ | residue-ligand cross-encoding |
| NLP | Counterfactually-augmented data | Kaushik et al., ACL 2020, https://aclanthology.org/2020.acl-main.711/ | CIVS counterfactual construction + Scheme D negatives |
| NLP | Hypothesis-only / annotation-artifact baselines | Poliak et al. EMNLP 2018 https://aclanthology.org/D18-1024/ ; Gururangan et al. NAACL 2018 https://aclanthology.org/N18-2017/ | the ligand-only and protein-invariant floors |
| RecSys | Factorization Machines / Field-aware FM (explicit 2nd-order interaction) | Rendle ICDM 2010; Juan et al. RecSys 2016 https://dl.acm.org/doi/10.1145/2959100.2959134 | Scheme A: exact variance attribution of the interaction term |
| ML / interpretability | Concept Bottleneck Models (intervenability by construction) | Koh et al., ICML 2020, https://proceedings.mlr.press/v119/koh20a.html | Scheme C: do()-intervention on the interaction concept |
| ML / explanation | Counterfactual GNN explanations; causal-spurious decoupling for OOD | https://axi.lims.ac.uk/paper/2410.15165 ; https://www.sciencedirect.com/science/article/abs/pii/S0957417426022633 | counterfactual test design; nuisance decoupling |
| Meta-learning | MAML / ANIL (feature reuse) ; CNP/NP (support-conditioned prediction) | Finn ICML 2017; Raghu ICLR 2020 https://openreview.net/forum?id=rkgMkCEtPB ; Garnelo ICML 2018 http://proceedings.mlr.press/v80/garnelo18a.html | Scheme E: support-conditioned interaction readout |
| Bioactivity ML | ActFound pairwise meta-learning | Feng et al., Nat Mach Intell 2024, https://github.com/BFeng14/ActFound | Scheme D: difference-native objective |
| Causal estimation | DML / R-learner / X-learner / CFR (nuisance residualization) | Chernozhukov 2018 https://academic.oup.com/ectj/article/21/1/C1/5056401 ; Nie & Wager 2021 https://academic.oup.com/biomet/article/108/2/299/5911092 ; Shalit ICML 2017 http://proceedings.mlr.press/v70/shalit17a.html | the cross-fitted residual objective in Schemes A/B |
| Genetics | Marginal epistasis test (detect interaction without joint features) | Crawford et al., PLoS Genet 2017, https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006869 | detection-first fallback if estimation is underpowered |

### 3c. Out of scope (requirements we do not have)

Structure-required ΔΔG methods (PremPLI-class, need co-complex coordinates);
MSA-required methods (DeepSequence/EVE-class, UniRef snapshot absent
locally); large external assay corpora for ActFound-style pretraining. These
are boundary references only.

---

## 4. The design principle (the actual innovation thesis)

An interaction claim is scientific only if the interaction is:

1. **Structurally explicit** — a named, separable term/pathway I(P,L), not an
   emergent property of a black-box concatenation. This gives an ablation
   handle.
2. **Load-bearing under the cold split** — removing I(P,L) measurably degrades
   cold-target performance beyond the additive main-effect model. This rules
   out "the interaction is present but unused."
3. **Counterfactually verifiable** — under ligand-swap and protein-edit
   interventions the prediction moves in the LABEL direction; under
   interaction-scrambling and nuisance transformations it does not. This
   rules out family-key / assay-batch / random-context shortcuts.
4. **Deployable** — computable from sequence + legal priors + ligand alone.

Every scheme below instantiates all four. The shared CIVS battery (Section 6)
operationalizes them identically across schemes, so the schemes are directly
comparable and Core Task 1 has one uniform standard of proof.

---

## 5. Five feasible schemes (model innovation points)

Each scheme lists: mechanism; why it completes Core Task 1; CIVS hooks;
literature grounding; innovation point; risk; feasibility with current
assets. All are <2M params on frozen features and trainable in minutes/seed
on CPU-class budgets, matching repo governance (no closed-form, no ridge,
end-to-end gradient training, keyed rng).

### Scheme A — FAIM: Field-Aware Interaction Machine (minimal identifiable reference)

- **Mechanism**: decompose protein into fields {ESM-global, KLIFS-pocket,
  optional family/conservation} and ligand into fields {ECFP, scaffold,
  pharmacophore}. Interaction is an explicit field-aware bilinear
  I(P,L) = Σ_{f,g} v_f(P)^T W_{fg} v_g(L). Main effects m_P, m_L are separate
  and orthogonalized (OID-style centering). Target = cross-fitted residual
  (DML/R-learner) so I explains only the nuisance-removed part.
- **Why it completes CT1**: the interaction term has EXACT variance
  attribution — drop W and you have the additive model; the cold-target gap
  is the interaction contribution by construction. It is the smallest model
  for which CT1-3/CT1-4 are cleanly measurable, so it is the mandatory
  identifiability reference against which B–E are judged.
- **CIVS hooks**: ablation = zero W_{fg}; counterfactual ligand-swap changes
  within-protein ranking in label direction; field-pair shuffle destroys it;
  ligand-only/protein-invariant floors are the negative controls.
- **Grounding**: FM/FFM/DeepFM; DML/R-learner.
- **Innovation point**: the first DTA model with per-field, exactly-
  attributable interaction terms and a pre-registered load-bearing gate.
- **Risk**: limited expressivity for nonlinear interaction; may underfit.
  That is acceptable — its job is attribution, not maximal accuracy.
- **Feasibility**: high; existing features; tiny param count.

### Scheme B — LCR: Ligand-Conditioned Residue routing (deployable FiLM/hypernetwork)

- **Mechanism**: frozen ESM-2 RESIDUE states h_i (full sequence, deployable,
  NO coordinate). A small ligand-conditioned generator emits per-residue or
  per-region modulation (γ_i(L), β_i(L)) (FiLM) or gates g_i(L) (target-
  conditional gating); modulated states h'_i = γ_i⊙h_i + β_i are pooled for
  the interaction readout; a parallel additive path bypasses modulation.
  This is the deployable cousin of CIIP-1A's oracle window: the ligand, not
  an annotation, decides where the protein representation is read.
- **Why it completes CT1**: CT1-1 (sequence+ligand only); CT1-2 (routes on
  cold families because routing depends on ligand × residue state, not
  memorized target); CT1-3 (γ→1, β→0 collapses to additive; the cold-target
  gap is the interaction); CT1-4 (counterfactual ligand swap shifts
  within-protein ranking; ligand-conditioned routing on a cold protein still
  shifts predictions).
- **CIVS hooks**: modulation ablation; ligand-swap counterfactual; routing-
  weight vs known binding sites is DESCRIPTIVE ONLY (never proof, per
  Jain & Wallace); shuffled-generator negative control.
- **Grounding**: FiLM; HyperNetworks; MCAN/ESIM co-attention; EGA-DTA target-
  conditional gating; CIIP-2 LCRR.
- **Innovation point**: deployable ligand-routed residue field — the
  interaction is conditional modulation of the protein representation by the
  ligand, requiring no mutation coordinate; directly addresses the CIIP-1A
  failure mode (oracle dependence) while keeping the residue-level richness.
- **Risk**: routing may collapse to a global (family-key) shift. CIVS
  nuisance-invariance + family-preserving shuffle detect this; a per-region
  routing-entropy regularizer can be pre-registered as a mitigation.
- **Feasibility**: high; ESM residue cache exists; generator is tiny.

### Scheme C — ICB: Interaction Concept Bottleneck (verifiability by construction)

- **Mechanism**: force prediction through a low-dimensional interaction
  concept z = I(P,L) ∈ R^d (d ≈ 8–16). z is trained to predict the within-
  protein ligand contrast (centered); ŷ = b_P + b_L + h(z). Concepts may be
  weakly supervised (KLIFS pocket class × ligand class posteriors) but this
  is optional.
- **Why it completes CT1**: a concept bottleneck is INTERVENABLE by
  construction (Koh et al.): you can run do(z → z') and check the predicted
  affinity moves as the concept semantics demand. This is the most direct
  possible answer to "does the model use the interaction?" — you literally
  intervene on the interaction and watch the output. CT1-3 by ablating z;
  CT1-4 by do()-counterfactuals.
- **CIVS hooks**: do(z) intervention direction test; z-ablation cold-target
  collapse; concept interpretability; ligand-swap changes z in the label
  direction.
- **Grounding**: Concept Bottleneck Models (Koh ICML 2020); counterfactual
  explanations.
- **Innovation point**: the first DTA model whose interaction usage is
  provable by do()-intervention rather than post-hoc attribution.
- **Risk**: bottleneck too tight → underfit; weak-supervision concept labels
  noisy. Mitigate with a pre-registered capacity sweep on val only.
- **Feasibility**: medium-high.

### Scheme D — DNT: Difference-Native Training (identifiability-by-design; the training-module innovation)

- **Mechanism**: NEVER train on absolute affinity. Train only on pairwise
  differences: within-protein ligand-ranking, within-ligand protein-ranking,
  and cross contrasts (ActFound-style pairwise meta-learning). Main effects
  cancel by construction, so ANY predictive success is definitionally
  interaction-driven — there is no main-effect pathway to leak into. For
  deployment, a separate, separately-audited nuisance/level head restores the
  absolute scale and is reported apart from the interaction head.
- **Why it completes CT1**: this is identifiability BY DESIGN — the strongest
  possible answer to CT1-4, because the training objective makes a
  main-effect shortcut impossible rather than merely controlled. Within-panel
  contrasts are exactly the estimand shown robust to assay variability
  (Nelen 2025) and isomorphic to ActFound's deployed design.
- **CIVS hooks**: identity-zero antisymmetry (bitwise); support/query label
  isolation; cold-target within-protein Spearman as headline; counterfactual
  negative pairs (same-parent wrong-mutation) as decisive negatives.
- **Grounding**: ActFound (Nat Mach Intell 2024); KronRLS/SimBoost pairwise
  DTA; within-panel contrast (Nelen J Cheminform 2025); the repo's own
  pairwise arms (U/V, arm 7).
- **Innovation point**: a training module (satisfies the "at least one
  innovation in the training module" rule) that makes interaction the ONLY
  learnable signal; paired with any of A/B/C as the scorer.
- **Risk**: absolute MSE needs the separate level head (which hits the level
  wall — acknowledged; the interaction head is judged on ranking/contrast,
  and the level head is audited separately). Pairwise sparsity on cold
  targets; mitigate with k-shot support (couples to Scheme E).
- **Feasibility**: high; matches existing pairwise/contrast infrastructure.

### Scheme E — SCNPI: Support-Conditioned Neural Process Interaction (few-shot bridge)

- **Mechanism**: for a cold target, a k-shot support set S = {(L_i, y_i)} is
  encoded by a permutation-invariant set encoder (DeepSets/NP) into a context
  c(S); the query prediction for L_q conditions the interaction readout on
  c(S). The set encoder is the transfer mechanism: it carries the target's
  conditional-interaction fingerprint from a few labelled ligands to new ones.
- **Why it completes CT1**: CT1-2 across the k∈{0,1,2,3,5} curve; CT1-3 by
  support-ablation (empty support → k=0 prior); CT1-4 by support-label
  isolation (shuffle support labels → k-shot gain destroyed), support
  permutation invariance, query-label isolation, and wrong-protein-support
  counterfactual (must destroy predictions).
- **CIVS hooks**: k-shot improvement curve; support-label isolation test;
  permutation invariance to ~1e-6; wrong-support counterfactual.
- **Grounding**: CNP/NP (Garnelo ICML 2018); MAML/ANIL; meta-learning kinase
  matrix completion; the repo's P-line CNP re-adjudication (arm 6, AD2).
- **Innovation point**: an explicit, manipulable few-shot transfer of the
  conditional interaction; the interaction is read out conditioned on
  support, so its usage is verifiable by support manipulation.
- **Risk**: sparse/noisy support on cold targets; NP variance. Mitigate with
  the deterministic DeepSets variant already re-adjudicated in the repo.
- **Feasibility**: high; P-line infrastructure exists.

---

## 6. CIVS — Counterfactual Interaction Verification Suite (the proof standard)

CIVS is architecture-agnostic and is run identically for every scheme. It is
the methodological core of this program: it turns "the model uses
interaction" into a falsifiable claim. All tests are pre-registered,
keyed-rng, parent/family-cluster bootstrapped (2000 draws), with
leave-one-family-out sign stability; bootstrap means never used as point
estimates; Spearman undefined (never 0) for constant predictions.

1. **Cold-target transfer** (CT1-2): double-cold + DTI-DG-style splits;
   MSE/RMSE, CI, Spearman, Pearson vs {Tanimoto transport, ligand-only,
   additive main-effect, frozen-ESM linear}; k∈{0,1,2,3,5}.
2. **Load-bearing ablation** (CT1-3): zero the interaction pathway;
   cold-target gap with bootstrap lo2.5 > 0 across families.
3. **Counterfactual ligand-swap** (CT1-4): fix a cold protein, swap ligand
   features at inference; within-protein ranking must move in the label
   direction (sign accuracy / Spearman of counterfactual shift vs true shift
   above a frozen threshold).
4. **Counterfactual protein-edit** (CT1-4): fix a ligand, edit the protein
   (deployable: ESM on the edited sequence / pocket perturbation); within-
   ligand ranking shifts in the label direction.
5. **Nuisance-invariance** (CT1-4): predictions invariant (within tolerance)
   to assay/family-key/batch permutations and to random-context (non-
   informative) protein perturbations; random-window / shuffled-pocket
   controls must NOT reproduce the interaction signal.
6. **Interaction-scrambling negative architecture** (CT1-4): same
   architecture with the interaction pathway shuffled/randomized at train
   time must collapse to baseline — rules out the architecture accidentally
   encoding main effects.
7. **Identifiability gate** (CT1-4/CT1-5): the interaction effect survives
   family-cluster bootstrap and leave-one-family-out; a pocket/family-key
   identifiability probe (as in the frozen boundary) must not absorb it.

A scheme "completes Core Task 1" only if CIVS 1–7 pass at pre-registered
thresholds. Attention heatmaps / binding-site overlap may be reported as
DESCRIPTIVE supplements only, never as evidence for any pillar.

---

## 7. Recommended main line, staging, and promotion gates

### Recommended main line

**Scheme B (deployable ligand-conditioned residue routing) as the mechanism +
Scheme D (difference-native objective) as the training module**, with
**Scheme A** as the mandatory minimal identifiable reference and **Scheme C**
as the verifiability/do()-intervention layer on the same backbone; **Scheme E**
reserved for the few-shot (k≥1) bridge. This combination: satisfies "at least
one innovation in the training module" (D) and one in the architecture (B);
is fully deployable (no coordinate); is falsifiable component-by-component
(each is a single toggle); and directly repairs the two audited failure modes
(oracle dependence; entangled centered estimand). It is compatible with, and
a strict generalization of, CIIP-2's OLR-Potential (B≈LCRR, D≈CFOIE
residual); if governance runs CIIP-2 first, this program's CIVS battery and
Schemes A/C must be prepended by addendum (see CIIP_SUCCESSOR §9.5).

### Staged execution (all stages preregistered + SHA-256; CPU smoke → single
### seed → multi-seed only after structure + negatives pass)

- **CT1-0 Governance & data audit** (read-only): freeze inputs (ESM residue
  cache, KLIFS, ECFP, Duong-Ly, BindingDB-Ki governed splits); coverage /
  parent-overlap / ligand-overlap / assay-semantics / censoring audit;
  leakage audit; re-derive the CIIP_SUCCESSOR §4 decomposition under
  preregistration; power table + MDE per CIVS test; freeze thresholds by
  dated addendum BEFORE any training. meta_test never read; test labels only
  for final evaluation of frozen arms.
- **CT1-1 Scheme A (FAIM)**: smallest identifiable interaction; run full
  CIVS. Establishes the attribution reference and validates the CIVS
  pipeline end-to-end.
- **CT1-2 Scheme B (LCR)**: deployable routing; run full CIVS; compare to A.
- **CT1-3 Scheme D objective on B**: difference-native training; CIVS;
  compare to the centered-residual objective (a univariate toggle).
- **CT1-4 Scheme C (ICB) do()-intervention** on the B backbone: the
  verifiability layer; CIVS + concept-intervention tests.
- **CT1-5 Scheme E (SCNPI)**: k-shot bridge; CIVS k-curve + support tests.
- **Promotion gates**: a scheme advances only if CIVS 1–7 pass at frozen
  thresholds on single seed AND reproduce on ≥3 seeds with family-cluster
  bootstrap lo2.5 > 0 and LOFO sign stability. No promotion on a single
  seed/parent/few pairs; no metric shopping; UNRESOLVED is a legal terminal
  state.

### Stop rules

Any scheme that collapses to the ligand-only / additive / family-key floor
under CIVS-6, or whose counterfactual direction test (CIVS-3/4) fails, is
falsified-as-tested and closed with a boundary entry. If ALL schemes fail
CIVS-2 (cold-target transfer), Core Task 1 is recorded as falsified under
current legal data and deployable constraints, and the boundary document is
updated — a legitimate terminal outcome per the programme goal.

---

## 8. Success criteria, deliverables, and constraints

### Success (Core Task 1 complete)

A promoted scheme passes CIVS 1–7 with: cold-target ranking/contrast gain
over the strong baseline set (bootstrap-stable, LOFO-stable); a significant
load-bearing ablation gap; counterfactual direction tests above frozen
thresholds; nuisance-invariance within tolerance; and metric superiority not
reducible to a re-fit level term. Completion is reported per-pillar
(SUPPORTED / NOT SUPPORTED / UNRESOLVED), never as a single number.

### Deliverables per stage

PREREGISTRATION.md + SHA-256; threshold addendum (frozen pre-training);
data-audit JSON/MD; RESULT.json (machine-readable all metrics); REPORT.md
(civil verdict table + authorization block); commands.jsonl; structure +
data-contract tests; SHA256SUMS; FAILURES.md; append-only sync to
history.md / task.md / EVIDENCE_LEDGER.md.

### Production constraints (binding)

No change to model/ or production scripts/ during research; no oracle
mutation-coordinate ESM in any deployment path; no CIIP potential into
BindingDB, no CIIP-1B, no BindingDB Bridge start without its own assay-
semantics qualification; no biological-mechanism claim from correct-vs-
random-window differences; no success claim from single seed/parent/few
pairs; no larger backbone/budget to mask non-identifiability; functional %
inhibition never relabelled Ki/Kd/pK; context-propagation magnitude never
cited as predictive value; attention heatmaps never cited as proof of
interaction usage; failure never reported as biological absence of
protein-conditioned signal.

---

## 9. What this report changes relative to the prior documents

1. Refocuses the entire programme on Core Task 1 as the SOLE objective and
   turns it into a five-pillar, falsifiable completion criteria set.
2. Adds the missing artifact in both the repo and the cold-target DTA
   literature: CIVS, a standardized counterfactual + load-bearing +
   nuisance-invariance battery that makes "the model uses interaction"
   provable rather than asserted.
3. Supplies five concrete, deployable, literature-grounded model-innovation
   schemes (A–E) with a recommended combination, instead of a single point
   design.
4. Explicitly couples Core Task 1 to the frozen level wall: conditional
   interaction is identified as the ONLY transferable object on a cold
   target, so proving its use is the mechanism for breaking the wall.
5. Keeps every prior safeguard: preregistration + SHA-256, sealed meta_test,
   no test-label leakage, no closed-form, single-seed-then-multi-seed,
   UNRESOLVED as a legal terminal state.
