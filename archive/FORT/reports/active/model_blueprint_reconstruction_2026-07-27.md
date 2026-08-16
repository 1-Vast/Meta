# Model Blueprint Reconstruction — Audit, Literature Synthesis, and Selected Model (2026-07-27)

**Type:** first-principles reconstruction of the model-side blueprint in `task.md`.
**Method:** full direct read of `task.md` (1,293 lines) and `history.md` Parts VIII–XII; two codebase/report
exploration agents; four literature campaigns (dual-cold DTA, protein/ligand foundation models,
matrix-completion/hierarchical statistics, LLM roles); one Oracle stress-test with one documented deviation.
**No experiment was run, no model trained, no label read, no API called.** `sealed_test_consumed=false`;
`confirmation_labels_read=true` (pre-existing; this record read nothing). `lexor/` untouched.

---

## 1. Audit of the current model blueprint

The following defects were identified in the 2026-07-27 `task.md` by direct reading. The reconstruction in
`task.md` Part 7 resolves them; historical statements are preserved or marked superseded, never silently
edited.

### 1.1 The `k=0` contract contradicts every zero-shot model in the ledger

Part 1.3 states `k=0` must **exactly return `b(d)`**. Parts 8.6 and 9.2 then define zero-shot models whose
`k=0` prediction is `b(d) + g(t,d)` with nonzero `g`. Part 8.6's bridge sentence ("with zero support labels
the `k=0` contract holds vacuously") is a reinterpretation, not a resolution. The contract was written for
the BM0-class support-posterior model (whose only target-dependent term is the posterior, with posterior
mean zero at `k=0`) and was silently carried into a ledger whose active models have a direct zero-shot
interaction term. **Resolution (adopted):** the contract is revised explicitly — the zero-shot model is
`b(d) + φ(d)ᵀw̄ + γ(aᵀk_t)(cᵀφ(d) − μ_ref)`; the support-posterior layer is dormant and keeps `Δ(0)=0`;
the ligand-side correction direction is centered under the training ligand reference measure so the
correction carries **zero within-target mean by construction** (it cannot act as a target intercept for any
target, seen or unseen). See §6 of this report for the Oracle critique and why centering answers it.

### 1.2 The §2.7 "600-component floor" was derived from the estimand §2.8 retired

§2.7's binding floor (≥600 components; 618 at 0.03; 1,391 at 0.02) is computed from the paired dispersion of
the **containment-mean contrast** (`inside/signal` per target, averaged). §2.8 then registers that **no
containment mean anywhere in the program is a usable effect size** (unbounded ratio; selection on a noisy
denominator). A floor derived from a retired estimand has no scientific standing. **Resolution (adopted):**
the containment-mean estimand is retired from active use; power is restated per-estimand (§9 of this
report): predictive paired ranking contrast ≈ **423** multi-family components (from the program's own
registered arm-heterogeneity figure MDE80 0.0614 at 101 components, `reports/active/panel_power_pd1.json`,
scaled by the √n rule); mechanism gate on `γ̂` via sign/signed-rank ≈ **70–155** components for median shifts
with per-component positivity 0.65–0.60. The multi-family requirement stands independently of power.

### 1.3 §2.5 cites Stage D magnitudes that §8.10/§8.11 forbid citing, and misreads its own theorem

§2.5 quotes `true − global = −0.232, LCB −0.083` as if it were an effect size. §8.10 registers that this
number flips sign under a one-target perturbation (mean → +0.565, driven by a 56× outlier) and that Stage D
magnitudes **must not be cited as effect sizes**. §2.5 also reads the dirty-IMC literature (Chiang–Hsieh–
Dhillon 2015; Yang & Ma arXiv:2605.17189) as proving that inexact row side-information "destroys the
cold-start sample-complexity advantage regardless of how `B` is parameterised." The same papers, and the
ledger's own §8.2, say the opposite: error degrades as **O(ε)** under bounded misspecification (optimal,
graceful), DirtyIMC's `XMYᵀ + N` absorbs misspecification in `N`, and Yang & Ma's penalized IMC/MC
interpolation exists precisely because the advantage degrades rather than vanishes. **Resolution (adopted):**
the pooled-ESM-2 closure stands on the *empirical* record (six parameterisations failed on one coordinate;
`esm − shared_global` sign p = 0.0013; `esm − own derangement` p = 0.82), the magnitudes are removed, and
the theorem is cited correctly — shared-plus-specific is the theoretically endorsed response to inexact side
information, not a prohibited one.

### 1.4 §2.2 overstates the null

"B0 ligand-only potency is the only repeatedly confirmed predictive signal" is true at the family-disjoint
tier (+0.700, RECRO L0) but erases the registered single-source protein-conditioned positives: PB6
(CFRI − random-protein **+0.0511 [+0.0129, +0.0905]** on Metz), KirHub A1-strict (true ESM − B0 **+0.0290
[+0.0083, +0.0497]**, beats shuffle/random, does not beat group centroid), TR-0 G1 (group − ligand
**+0.0400 [+0.0181, +0.0616]**). External corroboration: Mattsson & Walters (bioRxiv 2026.06.29.735309)
find ligand-only baselines at r = 0.66 (FEP+ 4) / 0.36 (OpenFE) on current benchmarks, collapsing to
r = 0.14 in the strictest novelty tier. **Resolution (adopted):** §2.2 is scoped — B0 is the only signal
confirmed *provenance-family-disjoint*; small protein-conditioned zero-shot signals exist on single-source
substrates, are not resolvable beyond coarse taxonomy, and have never been certified provenance-disjoint.

### 1.5 Innovation-budget accounting contradicts itself across Parts 8/9

Part 8's header says the PARC P1+P2 pair "spends the `rules.md` §3 model-side budget in full"; §8.11 says
"P2 was never built, so the model-side innovation budget is **not** consumed"; Part 9 then claims the budget
"is available." **Resolution (adopted):** a single accounting rule — an innovation is spent only when its
module is **built and run**; registered-but-never-built modules (P2, Mamba ladder rungs, K-LBP v1) spend
nothing. The active blueprint spends exactly 2 (§10 of this report).

### 1.6 Blocked architecture plans remain on the active path

Part 4 (Mamba admissibility + novelty boundary) and Part 5 (E0→F4 ladder with F2's five-rung Mamba
architecture comparison) remain titled "active execution sequence" although every supervision gate that
could unlock them has failed terminally (F0-P `PAPYRUS_F0_RAW_PROVENANCE_INSUFFICIENT_STOP`; the three
pretraining-anchor candidates all stopped; LEXOR L0 `LEXOR_L0_CORPUS_FRAME_INSUFFICIENT_STOP`, L1 locked).
E0's outcome is an engineering fact (`MAMBA_E0_PASS_COMPONENTS_AUTHORIZED_ONLY`) that supplies no
supervision. **Resolution (adopted):** both plans are moved to a marked *Superseded execution plans* part;
their outcomes are preserved as evidence; the Mamba-in-Transformer ladder and F1–F4 are no longer the
active path. This is not a claim that Mamba is useless; it is a claim that an architecture comparison
without a supervision gate is not an executable plan.

### 1.7 Duplicate and incompatible task definitions

The ledger carries at least four task definitions: within-target ranking (Spearman, Parts 5/8), containment
fraction (Part 8), `γ̂` measurement (Parts 8/9), absolute RMSE with safety gates (throughout). The primary
estimand is never chosen. **Resolution (adopted):** one primary estimand — target-component macro Spearman,
paired full-model vs nested null, component bootstrap — with RMSE (per-target query-centered, both ways
reported), calibration, and selective risk as secondary. Containment retired (§1.2).

### 1.8 Modules without a load-bearing role still hold contract status

Part 1.3's support-posterior contract governs a layer whose every route failed or fell below floor
(`BM0_FAIL_STOP`, `BM1_RR_FAIL_STOP`, `PANEL_GATE_PC_FAIL_STOP`; the ORRC v2 derivation that positive
target-conditioned precision cannot flip signs; `RB_DR_QMAPD_ORACLE_INFORMATION_FAIL_STOP` measuring the
k=4→16 channel at +0.0154 against a 0.0452 floor). **Resolution (adopted):** the layer is *dormant*, not
deleted — the contract is retained for it (`Δ(0)=0`, support labels never enter an encoder) and its
reactivation requires a support-information audit that beats the RB-DR-QMAPD floor. Zero-shot and few-shot
are separated as distinct models (an explicitly permitted conclusion of this reconstruction).

### 1.9 LEXOR power thresholds vs model-side power requirements

LEXOR L0/L4 gate at ≥30 provenance-independent components, MDE80 ≤ 0.03. The model-side predictive contrast
needs ≈423. §9.10 already flags that LEXOR could pass L0–L4 at real cost and still be unable to resolve the
model at L5. **Resolution (adopted, cross-track note only — `lexor/` and task.md Part 5 are unchanged):** a
LEXOR L4 pass at 30 components does **not** unlock the model's predictive gate; LEXOR's power thresholds
must be re-derived from paired arm-contrast dispersion before any LEXOR output is consumed by the model
track. This is registered as a precondition in the model blueprint, not as an edit to the LEXOR design.

### 1.10 Verification note on one citation

A background literature agent claimed the target-mirroring reference (bioRxiv 2026.06.29.735309) might be a
pMHC paper. Verified directly: the DOI resolves to Mattsson & Walters, *Identifying and Addressing
Systematic Data Leakage in Protein-Ligand Affinity Benchmarks* (bioRxiv, posted 2026-06-30, CC-BY 4.0;
NTAB benchmark, github.com/bamattsson/ntab). The citation in the LEXOR literature section is **correct**;
the agent's counter-claim is rejected.

---

## 2. Historical diagnosis — what failed, what is unresolved, what is underpowered

Classification follows the five categories registered in `history.md`'s reading guide. Full per-report
classification was produced by the digest agent and is consistent with this table except where noted.

| route | verdict | class | one-line reason |
| --- | --- | --- | --- |
| Physical pose/pocket/docking (FIRE, Gate-P/P2A, BridgeFIRE) | program CLOSED | falsified (well-powered) | native-pose upper bound below ligand-marginal; physics degrades descriptors |
| Pair-compatibility contrastive pretraining | `MULTISEED_FAIL` | mechanism does not transfer | retrieval AUC real (0.635>0.512); downstream null on 3 seeds |
| Assay-offset de-noising (Module A) | `DENOISE_FAIL` | mechanism does not transfer | 24% measurement de-noise, macro-ρ −0.051 cold-target (3/3 seeds) |
| GO-hierarchy protein pretraining | `BIOLOGICAL_PRIORS_NOT_YET_LOAD_BEARING` | does not transfer | real > shuffled, lost to random init (3-seed) |
| BM0 / BM1-RR / PANEL-PC posterior adapters | `*_FAIL_STOP` | mechanism absent at floor | correct support ≈ wrong/permuted support; RMSE-only gains |
| Transformer+Bayes+meta (TBM) / HIER prior | `*_FAIL_REVIEW` / `PROTEIN_CONDITIONED_PRIOR_NOT_LOAD_BEARING` | representation-specific | point gain +0.02 with protein controls destroyed |
| HQ-GBMA Stage D / PARC M0 | `*_FAIL_STOP` → **underpowered** | underpowered + representation-specific | 77 components vs 618–1,391 (containment estimand); pooled-ESM closure stands on sign tests |
| SCGD / QACO (Reinecke) | `SCGD_FAIL_STOP` / `QACO_FAIL_STOP` | mechanism absent on substrate | wrong-target support ≥ correct; protein-free ≥ protein |
| SI0 residual-covariance | `SI0_COVARIANCE_ALIGNMENT_FAIL_STOP` | underpowered/below-floor | +0.0271 < +0.03; median query-support Tanimoto 0.185 |
| MIF-NK MNI-0 | `MIF_MNI0_FAIL_STOP` | representation-specific | joint ≡ ligand-only; target destruction below threshold |
| TR-0 (taxonomic resolution) | `TR0_PREMISE_FAIL_STOP` | falsified premise | own-group-cold retains 77% of the group gain |
| RECRO R0 → L0 | `RECRO_SIGNAL_EXPLAINED_BY_PROVENANCE` | confounded | +0.334 → +0.090 (LCB<0) family-disjoint; 79.75% exact duplicates |
| SAFSA / MMP-X / TCOPA / Papyrus F0-P | four `*_STOP` | data-blocked | supervision structure absent after firewalls |
| KirHub U1 double-difference | `SIGNAL_PRESENT_CONFIRMATORY_GATE_BLOCKED` | data-blocked (real signal) | 34 gene units; no replicates |
| OpenBind U2 local support | `OPENBIND_K4_SIGNAL_ORACLE_UNCERTAIN_STOP` | representation-specific (real local SAR) | single protein; oracle LCB crosses zero |
| AMOB O0 | `OPEN_DATA_INSUFFICIENT_FOR_AMOB` | confounded/unlicensed | +0.434 without assay/doc IDs; superseded by RECRO L0 |
| CROSSDOC C3 | `CROSSDOC_C3_FAIL_PUBLIC_OVERLAP` | underpowered (signal credible) | +0.4946 on 13 units < 30 gate |
| CAPIT/ASPIRE C2 pocket | gate failed by 0.0014 | underpowered (ligand-warm oracle) | pocket beat ESM +0.0556 with all destruction controls |
| PANEL-EVIDENCE P0 / Mamba E0 | `*_PASS` | engineering only | no scientific credit, unlocks nothing |

**What truly failed (falsified mechanisms):** physical-structure load-bearingness; support-posterior/
positive-precision adaptation as a reordering channel; the group→target taxonomic resolution premise;
provenance-confounded "cross-document target signal".

**What remains unresolved (not falsified):** pocket restriction vs random positions (G3 directionally
positive, unresolved at 77); structure-token (3Di/SaProt) coordinate (never run — dependency-blocked);
mechanism-knowledge coordinate (K-LBP v1's card — registered, never audited); cross-family transfer of the
shared basis (untestable on current substrates).

**What is underpowered (do not read as refutation):** HQ-GBMA Stage D, PARC M0, KirHub A1 (effect 0.0290
vs 0.030 floor), CROSSDOC (13 units), RECRO R0 (MDE80 0.195 > 0.10), Reinecke (MDE80 0.0668).

**One correction to the digest agent:** the agent's RECRO entry inverts the corrective order. The
authoritative sequence (history.md Part XI) is R0 first (+0.334 cross-document) and **L0 second**, with L0's
provenance-family attribution (+0.0901, LCB<0; correct−wrong −0.0744) **overturning** R0's reading. This
report follows history.md.

---

## 3. Literature synthesis (four campaigns)

### 3.1 Strict dual-cold DTA — what the field has and has not established

Twelve+ primary sources (DeepDTA, GraphDTA lineage; UdanDTI drug-bias trap; DebiasedDTA; Q-BAFNet;
Ricci-GraphDTA; Co-Diffusion; CS-DTA; HGRL-DTA; MixingDTA; EBD-DTI; DTIAM; ProSmith; TransDTAP;
Golts et al. benchmark curation). **Not one paper reports a replicated dual-cold result with destruction
controls** (target shuffle, random protein, matched wrong target, ligand-only baseline). All evaluated
benchmarks are kinase-dominated (Davis/KIBA/Metz) with no binding-profile firewall. The strongest external
diagnosis converges with this program's: UdanDTI's drug-bias trap (dual-branch models collapse to the
ligand branch — the same signature as CFRI's ligand-shuffle costing nothing); Mattsson & Walters' target
mirroring and ligand-only r = 0.66 on FEP+ 4. **What the literature does not establish:** that any
architecture beats a ligand-only baseline under simultaneous accession/homology/profile/scaffold/neighbour/
provenance isolation; that reported dual-cold gains (CI ~0.80 on Davis, −42% MSE) are anything but
uncontrolled leakage. **Consequence for novelty:** a model with a full destruction battery under all eight
isolation axes is differentiated by evaluation rigor alone; the novelty claim must nevertheless rest on the
mechanism, not the evaluation.

### 3.2 Protein/ligand foundation models — what is information-disjoint from pooled ESM-2

Six representation families are genuinely disjoint from the falsified coordinate: (1) Foldseek-3Di
structure tokens (SaProt; SaBAN-DTI SOTA on LIT-PCBA AUROC 68.2 vs DrugCLIP 57.2; GenSPARC robust to AF2
structures); (2) ESM-IF1 pure-geometry encoder (residue-level, no sequence information); (3) ProSST
2048-state structure tokens (ProteinGym SOTA); (4) explicit 3D pocket embeddings (Uni-Mol pocket); (5)
structure-based residue graphs (3DProtDTA: cold-kinase-family MSE 0.305, CI 0.880); (6) joint co-folding
(AlphaRank/AlphaDTA — heavy, leakage-suspect via PDBbind memorization; rejected). AlphaFold DB makes
families (1)-(5) universally computable for an unseen target. Ligand side: MolCLR/iMolCLR (subgraph-removal
invariances; scaffold-split evaluated) and Uni-Mol 3D conformations resist scaffold memorization better than
SMILES LMs; CS-DTA's ablation found the pretrained encoders, not the interaction module, drive cold-start
robustness. **Consequence:** a structure-pocket coordinate is a legitimate, dependency-gated reopening
(PARC M0 never ran its 3Di arm; the closure was of *pooled ESM-2*, not of structure).

### 3.3 Matrix completion and hierarchical statistics — the correct theoretical frame

Jain & Dhillon 2013 (IMC, exact features); Chiang–Hsieh–Dhillon 2015 (DirtyIMC `XMYᵀ + N`, graceful
degradation); **Yang & Ma arXiv:2605.17189 (inexact side information degrades error as O(ε) optimally;
penalized IMC/MC interpolation)**; Jalali et al. 2010 (dirty shared-plus-specific multitask); Gelman & Hill
(partial pooling: new group → population mean with full uncertainty); Garnelo 2018/Kim 2019 (neural
processes underfit under task shift — matches the in-house posterior/meta failures); Amemiya & Fuller 1984
(errors-in-variables attenuation — the program's V_t-subtraction discipline is exactly this correction);
Peters et al. 2016 / Rothenhäusler et al. 2021 (invariant prediction/anchor regression — environments here
are assay/provenance, not the target axis, so invariance shrinks toward potency-only: rejected as a
predictor); Angelopoulos et al. 2023 (prediction-powered inference — relevant to certifying compiled
covariates against a gold subset); conformal prediction (exchangeability violated under target shift;
component-level split-conformal with caveats). **Consequence:** the selected model's shape —
shared low-rank + small feature-driven correction with an explicit interpolation gate — is precisely the
shared-plus-specific/dirty form the theory endorses for inexact side information. §2.5's reading is
corrected (§1.3).

### 3.4 LLM roles and knowledge-conditioned DTA

Ranked by evidence: **(1) fixed-schema deterministic knowledge encoding** — KGE_NFM (knowledge graph + NFM,
Nat. Commun. 2021, evaluated under protein cold-start), Harrison 2026 (community heuristics compiled to
fixed features; F1 0.38→0.59 on understudied families; PAINS false positives −58%), Yao 2025
(knowledge-based regularization), WideDTA (protein domains/motifs match full-sequence performance —
condensed structured target features suffice); **(2) frozen encoder embeddings** (TwinBooster assay-text
integration); **(3) offline evidence compiler** (RAG without parameter memorization); (4) teacher
distillation (contamination cascades); **(5) direct prediction — REJECTED** (Busch et al. 2026 blinding
study: closed-book LLM chemistry is memorization-dominated; Guo & Ding 2026: LLM SAR baselines win 1.9% of
156 comparisons; ChemPro/FGBench: degradation on deep chemical reasoning); (6) retrieval controller (poor
verification). **Consequence:** the Mechanism Card role (offline, fixed-schema, no affinity numerics,
contamination-probed) is the only admissible LLM participation; a closed-book probe is decisive and can
kill the route; a deterministic knowledge-graph/database encoding is the mandatory non-LLM baseline arm.

---

## 4. Candidate hypotheses and decision matrix

Six candidates were specified in full (new information, mathematical object, transfer arguments, parameter
count, data, expected effect, leakage channel, cheapest falsification, relation to failed routes). Four
serious ones follow; H5/H6 were rejected at specification.

**H1 — Mechanism-Card rank-1 correction (retained K-LBP core).** New information: time-frozen compiled
mechanism attributes (site chemistry, catalysis class, plasticity, conservation). Object: scalar coordinate
`aᵀk_t`. Target transfer: cards compilable for any target with literature (Task B) or from sequence alone
(de-identified arm, Task A). Ligand transfer: `φ(d)` is a descriptor function. Pair transfer: rank-1 gated
product with nested null. ~30–40 target-side parameters. Expected effect: +0.02–0.05 (if real). Leakage:
taxonomy collapse, study depth, hidden SAR, benchmark memorization — each gated. Cheapest falsification:
S1 card⊥taxonomy (no affinity model). Relation: replaces the falsified pooled-ESM coordinate.

**H2 — Structure-pocket coordinate (deterministic).** New information: predicted local structure
(Foldseek-3Di / ESM-IF1 pocket field from AlphaFold). Object: pocket structure-composition vector `s_t`
in the same rank-1 form. Target transfer: AFDB universal. Leakage: none beyond structure prediction itself.
Cheapest falsification: pocket-set shuffle + composition-vs-position (replicating the p=0.86 finding) +
card-vs-ESM redundancy, all affinity-free. Risks: within-kinome-only evidence; HonestAffinity's external
finding that pocket priors hurt strict tiers; dependency cost (two downloads). Never run in-house.

**H3 — Protein-free shared basis + hierarchical pooling only (the honest null route).** New information:
none (sequence-similarity/taxonomy only). Object: shared `B` + pooled component deviations. This is the
`γ ≡ 0` nested null with better shrinkage — the fallback if every coordinate audit fails. No novelty claim;
it is the comparison, not the candidate.

**H4 — Environment-anchored invariant residual.** Rejected as predictor: environments are assay/provenance,
not the target axis; RECRO L0 shows the target-specific part collapses family-disjoint, so an invariant
model provably shrinks to potency-only. Would formalize the negative, not escape it.

**H5 — LLM teacher distilled interaction model.** Rejected: contamination (Busch et al.), no positive DTA
evidence, dominated by H1's deterministic-schema version.

**H6 — Joint co-folding / contrastive complex model.** Rejected: physical-structure program closed
(well-powered); pair-compatibility pretraining `MULTISEED_FAIL`; MIF-NK ligand-prior collapse; PDBbind
memorization channel; heaviest compute.

**H7 — Mutation-anchored double-difference.** Rejected for the main route: the transferable variable
(per-target mutant response) is not available for an arbitrary unseen target at inference; confirmatory
gate blocked at 34 gene units.

### Decision matrix (scores 0–5; weighted equally; higher is better)

| criterion | H1 card | H2 structure | H3 null-route | H4 invariant |
| --- | :-: | :-: | :-: | :-: |
| information novelty | 4 | 4 | 1 | 2 |
| biological plausibility | 4 | 4 | 3 | 2 |
| statistical identifiability | 4 | 4 | 5 | 2 |
| dual-cold compatibility | 5 | 5 | 5 | 2 |
| ligand-only-collapse risk (5 = low) | 4 | 4 | 5 | 3 |
| target-identity-shortcut risk (5 = low) | 4 | 5 | 5 | 4 |
| benchmark-memory risk (5 = low) | 3 | 5 | 5 | 4 |
| parameter efficiency | 5 | 5 | 5 | 3 |
| data availability now | 3 | 2 | 5 | 3 |
| falsifiability | 5 | 5 | 5 | 3 |
| repository compatibility | 4 | 3 | 5 | 3 |
| computational feasibility | 3 | 3 | 5 | 4 |
| novelty vs published work | 4 | 3 | 1 | 2 |
| **total** | **52** | **51** | **52** | **37** |

**Selection (smallest defensible):** H1 and H2 are not stacked. The blueprint separates the **estimator**
(fixed: shared-global completion + rank-1 scalar-gated correction) from the **coordinate** `k_t` (audited;
exactly one enters the model). H1 is the primary coordinate for Task B, H2 the primary for Task A and the
deterministic competitor; the cheapest audits (R1, affinity-free) decide which coordinates survive to the
R4 head-to-head `γ̂` comparison; H3 is the nested null and the fallback verdict. This satisfies: credible
information advantage (compiled mechanism / predicted pocket), clear mathematical role (one scalar
coordinate, one direction), observable destruction tests (C1–C8 + pocket shuffles), plausible data path
(Metz train audits now; structure dependency gated; card API budget gated), nested fallback (`γ=0` →
shared-global, bit-exact), publishable novelty (§12).

---

## 5. The selected model (full specification)

Registered in `task.md` Part 7 as the active blueprint (supersedes K-LBP v1). Summary here.

**Tasks (reported separately, never mixed).**
Task A (closed-world): inputs = target sequence + deterministic derived features (predicted-structure arm
dependency-gated) + ligand structure. Task B (time-frozen knowledge-augmented): inputs += pre-cutoff
external mechanistic evidence compiled to the Mechanism Card; no exact query-pair affinity evidence.

**Prediction equation.**
`ŷ(t,d) = b(d) + φ(d)ᵀw̄ + γ·(aᵀk_t)·(cᵀφ(d) − μ_ref)`, `‖a‖=‖c‖=1`, `γ≥0`,
`μ_ref = E_train[cᵀφ(d)]` (ligand-reference centering → zero within-target-mean correction by
construction). `γ=0` recovers the shared-global model bit-exactly (nested null); the absolute null is B0.

**Components.**
`b(d)`: frozen validated B0 MLP (1034→128→1). `φ(d)`: shared ligand interaction basis — reference
instantiation is the PD-M-validated fixed 64-d design; a learned MLP variant (1034→128→r, r≤32) is a
registered same-budget arm decided at R3. `w̄`: shared-global coefficient fitted by the exact observed-edge
projection with nuclear-norm regularization (ORRC/PD-M machinery; convex, KKT-audited). `k_t`: ONE audited
coordinate — Mechanism Card (Task B; fixed schema `d_k ≤ 32`; target-intrinsic fields only; no affinity,
no family/group label, no source counts; missing/contradiction/abstain states; version-frozen compiler;
every numeric recomputed in code) or structure-pocket vector (Task A candidate; Foldseek-3Di/ESM-IF1 pocket
composition; dependency-gated) or sequence-only de-identified card (Task A candidate); pocket-conservation
vector admissible as a third arm under the same dependency authorization. `γ,a,c`: rank-1 scalar-gated
correction (~30–40 effective target-side parameters).

**Training.** Stage 1: cross-fitted `b(d)`; `(φ, w̄)` by regularized completion on the projected residual.
Stage 2: `(γ, a, c)` on the residual-after-shared with **within-target pairwise ranking loss** (the
program's one validated reusable objective; scale-invariant; correct when the unseen target's baseline is
unknowable) **+ `λ_γ|γ|` + unit-norm constraints**, cross-fitted on homology components; errors-in-variables
correction (subtract empirical `V_t`) whenever `γ` is estimated from cross-fitted coefficients. **No
counterfactual loss** (destruction is evaluation-only — see §8); no other auxiliary loss.

**Inference.** Zero-shot: one forward pass; no target-specific fitting; no support labels anywhere in the
active blueprint. Uncertainty: component-level split-conformal (Module B lineage: coverage@90 ≈ 0.86,
monotone selective risk) with the exchangeability caveat stated; card abstention states propagate to
selective prediction. Few-shot layer: dormant (`Δ(0)=0` contract retained); reactivation requires a
support-information audit beating the RB-DR-QMAPD floor.

**Parameter accounting.** Ligand-side: B0 298K (validated) + shared basis (fixed-φ: spectral; learned-φ arm
≈140K). Target-side: `w̄` (r) + `a` (d_k) + `c` (r) + `γ` ≈ **40–100 total amortized over all targets**
(≈0.4–1.0 per training target) — against 27,072 (≈304/target) for the retired ProteinGrassmann and 630 for
the retired §8.2 MLP.

## 6. Why it can generalize dual-cold (target / ligand / pair separately)

**Target transfer.** The correction sees the target only through `k_t`, which is defined for any target
(card: compiled from sequence and/or pre-cutoff literature; structure arm: computed from sequence via
predicted structure). No target identity, family label, accession, or binding profile enters. The shared
basis contributes target-generic interaction structure that PD-M showed transfers across held homology
components (3/8 directions feature-explainable and component-transferable after Holm correction).
**Ligand transfer.** `b(d)` and `φ(d)` are functions of descriptors, not identity; the scaffold/parent/
Tanimoto firewalls make the evaluation ligands genuinely novel; the ligand basis is shared across all
targets and trained on ~10²–10³ ligands, far from the target-side capacity wall. **Pair transfer.** The
pair term is a rank-1 product of a target scalar and a ligand direction — the smallest deviation from the
protein-free null that can still reorder a target's ligands; `γ̂` with a component bootstrap is itself the
preregistered measurement of whether the coordinate carries anything. The honest boundary: cross-family
transfer of the shared basis is an **untested hypothesis** (all identified substrates are kinase-only);
the blueprint claims a kinase-scale mechanism audit plus the exact conditions under which cross-family
claims become testable, nothing more.

## 7. Simplicity audit

Every retained module and its indispensable contribution:

| module | indispensable contribution no other module provides |
| --- | --- |
| `b(d)` ligand-only base | the validated potency prior; the absolute null; RMSE anchor |
| `φ(d)ᵀw̄` shared-global | the only empirically validated transferable interaction structure (PA4, PD-M); the nested null |
| `γ(aᵀk_t)(cᵀφ(d)−μ_ref)` | the sole target-specific channel; one scalar coordinate; the object under audit |
| exact observed-edge projection | provable main-effect removal; prevents base absorption (I0/CFRI failure mode); convex audited fit |
| within-target ranking loss | the only validated reusable objective; endpoint-aligned; baseline-unknowable-correct |
| `V_t` errors-in-variables subtraction | prevents manufactured effects from treating estimated coefficients as noise-free |
| component split-conformal | validated reliability layer (Module B); uncertainty without extra innovations |
| destruction family (evaluation-only) | the falsification apparatus; no training-time role |

Nothing else is in the model. No Transformer/Mamba/GNN stack, no cross-attention, no MoE, no docking or
pose channel, no support posterior, no auxiliary self-supervised loss, no free LLM embeddings.

## 8. Leakage and contamination audit

**Eight isolation axes** (unchanged + §2.6): target accession; homology/pocket family; binding-profile
correlation; ligand parent; Bemis–Murcko scaffold; high-similarity neighbourhood; assay/document
provenance; component dependence.

**Card-specific (evaluation-only family; training uses no corruption):**
C1 named vs de-identified vs evidence-redacted card · C2 closed-book affinity probe (train cells only;
cell-level power; a recovery kills the route) · C3 wrong-target card (exposure-matched) · C4 random card at
matched `d_k` and field marginals (**pipeline-free**) · C5 field/section ablation · C6 card⊥taxonomy
(TR-0 precedent; any field predictable from family/group is struck) · C7 popularity-only arm must be null
(**pipeline-free**; study-depth is the registered fourth channel) · C8 disjoint corruption families rule —
retained for any future training-time corruption use.

**Mandatory LLM contamination analysis (answers).** Target-name memory → C1. Exact affinity memorization →
C2 (Busch et al. methodology: blinding + progressive masking). Hidden SAR retrieval → partially mitigated
only; registered honestly: a field can be a true mechanism claim AND a compressed SAR summary; the residual
risk is why card compilation for evaluation targets prefers de-identified-sequence cards and why the
route's primary claims wait on a powered substrate. Benchmark-paper memorization → C2 + time-frozen cutoff
+ the program's note that a failed probe does not prove absence of contamination. Target popularity /
study depth → schema ban on confidence/source-count features + C7. Retrieval of chemically equivalent
query compounds → ligand-side firewalls are independent of the card; C3 wrong-target card breaks the
chemotype-target link. Pretraining knowledge predating the cutoff but containing evaluation labels → the
closed-book probe is run on train cells where labels are known; version freeze + prompt hash per card.
Closed-source version drift → compiler version-frozen and hashed; a model/prompt change voids every card.
Training vs evaluation corruptions → independently seeded, disjoint families, declared before any use.

**Oracle deviation record.** The Oracle recommended retaining a directional counterfactual margin penalty
in training. **Rejected, documented:** artifacts shared between true and corrupted cards are invariant to
the penalty (it cannot detect what it cannot vary); the systematic confounds the penalty targets (study
depth, popularity) are covered by pipeline-free controls (C4, C7, C2); re-adding the penalty reintroduces
the C8 self-manufacture hazard for zero marginal coverage. All other Oracle recommendations were adopted:
rank-2 sensitivity at R4 (non-gating); ligand-reference centering of the correction (stronger than the
Oracle's query-intercept fix — adopted instead of it); structure-first audit sequencing with the head-to-head
at R4; conservation arm registered as admissible under the structure dependency; cross-family scoped to a
hypothesis.

## 9. Data and power plan (LEXOR-independent)

**Usable now (no `lexor/` edit):** Metz train cells (12,560 cells, 101 components; interaction identified,
PA2 0.663) — R1/R3/R4; ChEMBL-37 train raw per-record extract (assay/doc IDs; provenance families via the
RECRO method) — R1 study-depth and provenance audits; Reinecke train anchors — secondary audit substrate;
KLIFS 85-residue pockets — structure-arm kinase instantiation; AlphaFold DB + Foldseek/ESM-IF1 or SaProt —
structure arm (two user-authorized downloads, dependency gate); UniProt/InterPro/PDB annotations —
deterministic no-API proxy card fields for the R1 dry run before any API budget; PLINDER raw tables —
ligand-warm diagnostics only.

**Blocked for prediction:** every identified substrate is kinase-only; the only multi-family graph with
depth (ChEMBL-37 dual-cold) fails the identifiability floor (P0_CYCLE_A 0.356 < 0.5); no open substrate
meets ≈423 multi-family components × ≥40 scaffold-diverse query ligands (`NO_OPEN_POWERED_INDEPENDENT_PANEL`;
BindingDB-native 38 targets ≥40; Davis 12 queries/target; SPD median 14). **Unlocks:** a future LEXOR L5
substrate **only after** its power thresholds are re-derived per §1.9; a prospective factorial panel; or a
newly identified multi-family panel meeting the same bar. The model status is therefore **audit-open,
prediction-blocked** — and the audits can each kill the route.

**Power restatement (per-estimand, planning values; each stage re-derives its own threshold before
running):** predictive paired ranking contrast ≈ **423** multi-family components (MDE80 0.0614 @ 101 →
101×(0.0614/0.03)²); mechanism gate on `γ̂` (sign/signed-rank) ≈ **70–155** components for positivity
0.65–0.60 — Metz's 77–101 sits inside this band, hypothesis-generating only (within-kinome); containment
estimand retired (§1.2).

## 10. Cheapest-first execution ladder

| stage | content | primary statistic | pass | fail-stop | a pass does NOT prove |
| --- | --- | --- | --- | --- | --- |
| R0 | historical + mathematical audit (this report) | — | registered | — | — |
| R1 | external-information audit, affinity-free: card⊥taxonomy (C6), card-vs-ESM redundancy, study-depth (C7), deterministic-proxy card dry run; structure arm: pocket⊥taxonomy, structure-vs-ESM, composition-vs-position | redundancy/orthogonality effect sizes vs frozen thresholds | coordinate(s) survive to R4 | `R1_COORDINATE_IS_TAXONOMY_STOP`, `R1_COORDINATE_IS_POPULARITY_STOP`, `R1_COORDINATE_REDUNDANT_WITH_ESM_STOP` | that the coordinate predicts affinity |
| R2 | contamination audit (card route only; explicit API authorization required): closed-book probe on train cells (C2), named vs de-identified vs redacted (C1) | closed-book recovery rate; named−deidentified gap | card eligible for R4 | `R2_CLOSED_BOOK_RECOVERY_STOP`, `R2_NAME_MEMORY_STOP` | absence of contamination (a failed probe is not proof) |
| R3 | synthetic estimator sensitivity under **empirical `V_t`** | `γ̂` recovery of known `γ*` within CI across heteroscedastic regimes | estimator certified | `R3_ESTIMATOR_INSENSITIVE_NO_DECISION` | that real `γ` is nonzero |
| R4 | one-seed train-only mechanism pilot on Metz train: `γ̂` head-to-head across surviving coordinates; rank-2 sensitivity (non-gating); C3–C5 | `γ̂` + component sign/signed-rank | coordinate-specific `γ̂` LCB>0 with all controls | `R4_COORDINATE_NOT_LOAD_BEARING_STOP` | predictive accuracy; cross-family transfer |
| R5 | powered predictive pilot | paired component-macro Spearman contrast vs nested null and B0 | LCB>0 at frozen floor | `R5_PREDICTIVE_GATE_FAIL_STOP` | — |
| R6 | multiseed reproduction | seed×arm interaction | reproduces R5 | `R6_NOT_REPRODUCIBLE_STOP` | — |
| R7 | independent confirmation | preregistered confirmatory contrast | confirm | `R7_CONFIRMATION_FAIL_STOP` | — |

R5–R7 are **blocked** (§9). R1 and R3 are authorized to plan; R2 additionally requires the explicit API
budget authorization; R4 requires R1+R3 passes and, for the card arm, R2.

## 11. Rejected alternative blueprints (with reasons)

Transformer/Mamba architecture ladder (F2) — no supervision gate; architecture size is not an innovation.
CFRI-class cross-attention interaction head — failed Gates Z/PB; 282K target-side parameters against the
capacity wall. HQ-GBMA Grassmann map — 27,072 parameters on ~89 targets; underpowered gate; over-parameterized
by two orders of magnitude. BM0/BM1/PC posterior adapters — positive-rescaling ablation; cannot flip signs;
below floor. SCGD/QACO — failed on Reinecke with correct-support losses. Joint co-folding/contrastive —
physical closure + `MULTISEED_FAIL` + memorization channel. Learned KG-embedding target features —
unauditable, popularity-confounded; deterministic schema encoding is the retained form of that idea.
Unrestricted RAG / free LLM embeddings — unauditable, version-drift and contamination channels.
LLM-direct prediction — memorization-dominated (Busch et al.); wins 1.9% of comparisons (Guo & Ding).
Teacher distillation — contamination cascades into the student. Invariant/anchor-residual — environments
are provenance, not targets; provably shrinks to potency-only (RECRO L0). Mutation double-difference —
transferable variable unavailable at dual-cold inference; data-blocked.

## 12. Novelty statement (restrained)

If R1–R4 pass, the defensible claim is: **the first strict dual-cold DTA formulation in which (i) the
target-specific zero-shot correction is constrained to a single scalar-gated rank-1 direction along an
audited, fixed-schema mechanism coordinate, (ii) the protein-free shared-global completion is the exact
nested null, and (iii) contamination/taxonomy/popularity destruction is a first-class preregistered gate
family.** This is a claim about the *exact combination* — transferable information + low-complexity
mechanism + strict dual-cold + falsification design. It is explicitly **not**: first LLM in DTA (precedent
exists), first multimodal DTA, first Bayesian DTA, or first mechanism-aware model (KGE_NFM, WideDTA,
Harrison precede it). No claim is made before the audits; a null `γ̂` with a certified estimator (R3) and a
complete destruction battery is itself a publishable negative.

## 13. Claims not yet earned

The card carries affinity-relevant information beyond taxonomy, ESM, and popularity (R1 decides). The LLM
compiler is uncontaminated for this panel (R2 decides; a pass is not proof of absence). The estimator
recovers `γ` under realistic noise (R3 decides). `γ̂ > 0` on Metz train (R4 decides). Anything at all about
multi-family or cross-family prediction (R5+, blocked). The shared basis transfers across families
(untested hypothesis). Few-shot support adds anything (dormant; below floor at last measurement).
