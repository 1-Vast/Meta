# K-LBP — Knowledge-Conditioned Low-Rank Biological Prior: critical evaluation and blueprint revision

**Date:** 2026-07-27. **Type:** blueprint revision + critical analysis. **No experiment run, no model
trained, no data read, no LLM/API call.** `sealed_test_consumed=false`;
`panel_development_labels_read=false`; `confirmation_labels_read=true` (pre-existing).

**Recommendation: `K_LBP_REGISTERED_AUDIT_ONLY__IMPLEMENTATION_DEFERRED_TO_LEXOR`.**
K-LBP is accepted as the registered candidate replacement for the failed protein-coordinate strategy. Its
audit stages (S0–S2) are authorised because they can *kill* the route cheaply. Its model stages (S3–S4)
are **not** authorised, because §2.7 is unchanged and no available substrate can resolve the gate.

---

## 1. The statement this revision is required to make first

> **A better biological coordinate cannot solve insufficient statistical power.**

`gamma` is a scalar contrast between two nested models. Its resolution is set by the number of independent
components and the per-component dispersion — **not** by the quality of `k_t`. Metz supplies 77 evaluable
components; the ranking contrast needs ~423 and the containment estimand needs 618–1,391 (§2.7). A perfect
Mechanism Card would still yield a `gamma_hat` whose interval spans zero. K-LBP changes the
*representation*; it does not change the *identifiability*.

### 1.1 The three-way distinction the project must maintain

| claim | status | basis |
| --- | --- | --- |
| **Biological signal is absent** | **NOT established anywhere in this program** | no adequately powered test has ever been run |
| **A specific biological coordinate is not identifiable on a specific substrate** | **established, narrowly** | pooled ESM-2 on Metz: `esm − shared_global` median −0.117, sign p = 0.0013, Wilcoxon p = 0.0003; `esm − esm_wrong_target` p = 0.82 (indistinguishable from its own derangement) |
| **The data are insufficient** | **established and binding** | §2.7: 77 components vs 423–1,391 required |

**Previous failures did not disprove biological priors.** `PARC_M0_COORDINATE_NOT_LOAD_BEARING_STOP` and
`HQGBMA_STAGE_D_FAIL_STOP` were underpowered by 18×–181× in units (§8.11) and are re-read as
`UNDERPOWERED`. The only *positively identified* negative is about pooled ESM-2 specifically, and even
that is a statement about one coordinate on one kinase panel.

---

## 2. Critical analysis of the proposal

### 2.1 Do Mechanism Cards carry genuinely new information relative to sequence embeddings?

**Unknown, and — importantly — testable without touching affinity.**

*For* [THY]: a card asserts curated, source-bound propositions that are not linear functionals of
sequence and that ESM-2's masked-token objective was never asked to expose — cofactor identity,
documented conformational states, presence of a covalent-targetable nucleophile in the site, reported
cryptic pockets, induced-fit propensity. ESM-2 encodes evolutionary sequence statistics; these are
*experimental annotations about mechanism*. They are different in kind, not merely a different projection.

*Against* [PRJ]: two measured facts push the other way.
1. On the identified panel, **aligned residue position added nothing over pocket amino-acid composition**
   (`parc − pocket_composition`, p = 0.86). Many card fields are positional/structural claims about the
   site; that is weak but real counter-evidence that site geometry is where the recoverable signal sits.
2. **TR-0 is the precedent to fear.** The coarse-kinase effect is real (`group − ligand_only` +0.0400
   [+0.0181, +0.0616]) but **not resolvable to the KLIFS group** (`group − group_cold` +0.0091
   [−0.0076, +0.0256]): a generic all-other-kinases centroid retains ~77% of the gain. A Mechanism Card
   whose fields are largely determined by kinase group would reproduce that null exactly, and would look
   like a success on every control except the own-group-cold falsification.

**Consequence for the plan:** the first substantive audit is not an affinity experiment. It is
**card ⊥ taxonomy**: measure how much of the card is predictable from family/group labels and from pooled
ESM. If a card is a relabeling of taxonomy, the route dies for the cost of compiling 111 cards, with no
affinity model and no power requirement.

### 2.2 Memorization risk — three distinct channels, and the standard controls catch only two

| channel | what leaks | caught by |
| --- | --- | --- |
| **Target-name memory** | the model recalls facts keyed to "CDK2" rather than reasoning from evidence | **named vs de-identified-sequence card** |
| **Benchmark memorization** | the model has seen Metz (Nat. Chem. Biol. 2011) and can reproduce its values | **closed-book affinity probe** |
| **Hidden SAR retrieval** | the card encodes *which chemotypes bind this target* without ever stating a number or a compound | **neither, reliably** |

The third is the dangerous one and it must be stated plainly. A card field such as "prefers small
hydrophobic gatekeeper-adjacent substituents" is simultaneously a legitimate mechanism claim and a
compressed SAR summary. It would **pass** target-shuffle, matched-wrong-target and random-feature
controls, because it genuinely is target-specific information — it just is not *transferable biology*, and
it would evaporate on a target the model has never read about, which is exactly the dual-cold case.

Partial mitigations, in decreasing strength: (i) run K-LBP on a substrate whose documents the compiler has
not seen — i.e. **a LEXOR-recovered corpus, not Metz**; (ii) de-identified-sequence arm as the primary
card; (iii) field-level ablation to locate which fields carry the effect; (iv) a schema that forbids
chemotype-valued fields outright. **Residual risk cannot be driven to zero on a public 2011 panel.** This
is, on its own, a sufficient reason to defer implementation to LEXOR rather than to run K-LBP on Metz.

### 2.3 A fourth leakage channel not in the proposal: study depth

`confidence`, `n_sources` and `contradiction` are **popularity statistics**. Well-studied targets have more
literature, higher card confidence, *and* more ChEMBL data, deeper panels and different chemotype
coverage. If confidence enters `k_t` as a feature, the model can learn "well-studied ⇒ different SAR",
which is a provenance/curation artifact of exactly the class §2.2 forbids.

**Registered rule:** evidence-confidence fields may be used **only** as observation weights or as an
abstention trigger. They must never enter `a^T k_t`. `n_sources` must never be a feature in any form.

### 2.4 Does K-LBP solve a model problem or only change the representation?

**Only the representation.** It is a strictly better-shaped representation than pooled ESM or PARC —
correct degrees of freedom, nested null, no 3Di dependency, field-level ablatable — but it does not touch
the binding constraint. §2.7 applies to K-LBP verbatim.

### 2.5 Do current datasets have enough power? No.

| endpoint | measured dispersion | components needed | available |
| --- | --- | ---: | ---: |
| target-macro Spearman contrast (predictive) | Metz arm-heterogeneity MDE80 0.0614 at 101 components, floor 0.03 | **≈ 423** | 101 |
| containment (train-only mechanism) | MAD-sigma 0.266 | **618** (0.03) / **1,391** (0.02) | 77 |

### 2.6 Recommendation on timing

**Audit now (S0–S2); implement never on Metz; defer S3–S4 to LEXOR.** This is the "only audited" +
"postponed until LEXOR creates a suitable substrate" option, in that order. The audits are worth running
now precisely because they are *asymmetric*: each can terminate the route cheaply, and none can authorise
it.

---

## 3. Mathematical correction to the proposed form

The proposal is

```text
Phi(t,d) = Phi_0(d) + gamma (a^T k_t)(b^T h_d) v
```

Three corrections are needed before this is estimable.

**(i) `v` is redundant for a scalar prediction.** With `g(t,d) = w^T Phi(t,d)`, the correction collapses to
`gamma (a^T k_t)(b^T h_d)(w^T v)`, so `v` is absorbed into the scale. `v` becomes meaningful only if the
*basis itself* must move — which matters only for a future support posterior over the coefficient. Since
few-shot is deferred (`Delta_info = +0.0154 [−0.0155, +0.0464]`, below the usable floor), **drop `v` now**
and reinstate it only if few-shot is reauthorized.

**(ii) `h_d` must be the low-dimensional ligand coordinate, not the raw descriptor.** Degrees of freedom,
against ~89 training targets:

| choice of `h_d` | parameters `a`+`b`(+`v`)+`gamma` | per target |
| --- | ---: | ---: |
| raw 1034-d Morgan+descriptor | 24 + 1034 + 6 + 1 = **1065** | 12.0 |
| 64-d ORRC ambient ligand PCA | 24 + 64 + 6 + 1 = 95 | 1.07 |
| **`r`-d shared-basis coordinates (recommended)** | 24 + 6 + 1 = **31** | **0.35** |

Putting `b` on the raw descriptor reintroduces on the *ligand* side exactly the over-parameterization that
sank the protein side (27,072 params / 89 targets ≈ 304 per target). Pin `b` to the already-validated
shared basis.

**(iii) The scale is non-identifiable, so the penalty must be applied to a normalized parameterization.**
`gamma`, `a`, `b` trade off multiplicatively, so `lambda_gamma |gamma|` is meaningless unless
`||a|| = ||b|| = 1` and `gamma >= 0` are enforced. Otherwise the optimizer shrinks `gamma` and inflates
`a` at zero cost, and the reported `gamma_hat` is uninterpretable.

**Corrected minimal form (30 effective parameters):**

```text
g(t,d) = phi(d)^T w_bar  +  gamma * (a^T k_t) * (c^T phi(d))
         ------------------  ---------------------------------
         shared-global        rank-1 biological correction
         (unchanged)          ||a|| = ||c|| = 1,  gamma >= 0
gamma = 0  ==>  exactly the current shared-global model
```

### 3.1 A hazard in the counterfactual objective

`L_counterfactual` forces true cards to beat corrupted cards. If the **same corruption draw** is used for
training and for the destruction control, the control is guaranteed to pass and proves nothing —
the objective would manufacture its own gate result. **Registered rule:** corruption draws used in the
training objective and those used in evaluation must come from disjoint, separately seeded families, and
the evaluation corruption must be held out and declared in the preregistration before training.

---

## 4. Mechanism Card schema (fixed, frozen, family-agnostic)

Target: `d_k <= 32` after encoding. Every field is a *target-intrinsic mechanism* property. No affinity,
no compound, no chemotype, no SAR, no benchmark value, no family/group label.

```text
A. Functional state                     encoding
   enzyme_or_receptor_class             one-hot(6)   kinase|protease|GPCR|NR|channel|other
   cofactor_requirement                 one-hot(5)   ATP|NAD/FAD|SAM|metal|none
   catalytic_mechanism_class            one-hot(4)
   documented_conformational_states     ordinal(0-3)

B. Site chemistry (qualitative only)
   site_polarity_class                  ordinal(0-2)
   hbond_donor_availability             ordinal(0-2)
   hbond_acceptor_availability          ordinal(0-2)
   buried_hydrophobic_subpocket         binary
   charged_residue_in_site              one-hot(3)   acidic|basic|none
   metal_coordination_in_site           binary
   covalent_targetable_nucleophile      binary

C. Plasticity
   induced_fit_reported                 binary
   cryptic_pocket_reported              binary
   local_disorder_near_site             ordinal(0-2)

D. Site conservation
   site_sequence_conservation           ordinal(0-2)
   allosteric_site_documented           binary

E. Evidence metadata  -- NEVER a feature; weights and abstention only
   per_field_confidence                 {high|med|low|absent}
   contradiction_flag                   binary
   source_ids                           list of PMID / UniProt / PDB accessions
   compiler_version, prompt_hash, temperature, seed, retrieval_date
```

**Hard schema rules.** (1) Section E never enters `a^T k_t`. (2) Any field found in the S1 audit to be
predictable from family/group labels above a frozen threshold is **struck from the card**, not
down-weighted. (3) Absent fields are encoded as an explicit missing indicator and abstain; they are never
imputed. (4) Every non-missing field must carry ≥1 resolvable source id, and a card with an unresolvable
source is rejected, not repaired. (5) The compiler is version-frozen and hashed into the manifest; a
prompt or model change voids every card.

---

## 5. Controls

The five mandated, plus three that the analysis above shows are necessary.

| # | control | detects | required outcome |
| --- | --- | --- | --- |
| C1 | named-target card vs **de-identified-sequence card** | target-name memory | de-identified arm must not be materially worse; if named ≫ de-identified, the gain is recall |
| C2 | **closed-book affinity probe** | benchmark memorization | model must not recover panel affinities without being shown data |
| C3 | wrong-target card (exposure-matched derangement) | biological relevance | true ≫ wrong, LCB > 0 |
| C4 | random card at matched `d_k` and matched field marginals | information content | true ≫ random, LCB > 0 |
| C5 | evidence-ablation (leave-one-field-out, leave-one-section-out) | which fields carry the effect | effect must localise to ≥1 section, not be diffuse-and-unremovable |
| **C6** | **card ⊥ taxonomy audit** (mutual information with family/group; TR-0 precedent) | card is a relabeled taxonomy | card must retain information after conditioning on group |
| **C7** | **study-depth control** — card vs `n_sources`/confidence-only features | popularity shortcut | popularity-only arm must be null |
| **C8** | **held-out corruption family** (§3.1) | self-fulfilling counterfactual objective | evaluation corruptions disjoint from training corruptions |

Plus, unchanged and mandatory: provenance-family-disjoint sensitivity, the §2.6
binding-profile-correlation firewall, and the §2.8 estimator rules (ratio of sums, no denominator
selection, sign/signed-rank primary inference, denominator CV reported).

---

## 6. Staged plan, with what each stage can and cannot conclude

| stage | what it does | cost | can kill? | can authorise? |
| --- | --- | --- | --- | --- |
| **S0 leakage audit** | C2 closed-book probe on train cells; C1 named vs de-identified | 1 bounded API budget; ~12,560 train cells of power (cell-level, **not** component-limited) | **yes** | no |
| **S1 informativeness audit** | C6 card ⊥ taxonomy; card vs pooled ESM redundancy; C7 popularity control. **No affinity model.** | 111 cards, no training | **yes** | no |
| **S2 synthetic sensitivity** | estimator recovery of a known `gamma*` using the **empirical `V_t`**, not `sigma^2 I` | ~1 GPU-hour | yes (if insensitive → no decision possible) | no |
| **S3 train-only mechanism gate** | `gamma_hat` with sign/rank inference | — | — | **BLOCKED by §2.7** |
| **S4 predictive gate** | paired arm contrast | — | — | **BLOCKED by data gate** |

S2 explicitly repairs the disclosed defect in PARC M0's G0 control, which used `V_t = sigma^2 I` and
therefore demonstrated sensitivity only under well-conditioned noise.

### 6.1 Stop conditions (any one is terminal)

1. C2: the compiler recovers panel affinities closed-book at material accuracy → **cards are contaminated
   at the root on this substrate**; `K_LBP_S0_CONTAMINATED_STOP`.
2. C1: named card materially beats de-identified card → the signal is recall, not biology;
   `K_LBP_S0_NAME_MEMORY_STOP`.
3. C6: the card is predictable from family/group above the frozen threshold → TR-0 reproduction;
   `K_LBP_S1_CARD_IS_TAXONOMY_STOP`.
4. C7: the popularity-only arm matches the full card → `K_LBP_S1_STUDY_DEPTH_SHORTCUT_STOP`.
5. S2: the estimator cannot recover a known `gamma*` under empirical `V_t` →
   `K_LBP_S2_ESTIMATOR_INSENSITIVE_NO_DECISION`.
6. On an eventually powered substrate: `gamma_hat` interval spans zero, or any of C3–C5 fails →
   `K_LBP_S3_NOT_LOAD_BEARING_STOP`.

---

## 7. Relationship with LEXOR

K-LBP and LEXOR are **complementary and share one firewall discipline**, and K-LBP is downstream.

* **Shared rule:** no LLM output is ever a label, a numeric, or a hidden state. LEXOR transcribes verbatim
  measurement strings; K-LBP emits fixed-schema qualitative mechanism attributes. Both recompute every
  number deterministically in code. Both are version-frozen, prompt-hashed and source-bound.
* **Shared control:** the closed-book contamination probe is the same instrument in both programs
  (LEXOR §6.6.2, K-LBP C2). An S0 result is reusable by LEXOR and vice versa.
* **Dependency:** LEXOR supplies the *substrate*; K-LBP supplies the *coordinate*. K-LBP S3/S4 unlock only
  on `LEXOR_L5_PROVENANCE_DISJOINT_TARGET_SIGNAL_IDENTIFIED` **and** a §2.7-compliant component count, or
  on a prospective factorial panel meeting the same bar.
* **Why this ordering also reduces contamination:** a LEXOR-recovered corpus is the only realistic route to
  targets whose measurement documents the compiler has not memorized — which is the *only* real mitigation
  for the hidden-SAR-retrieval channel (§2.2).
* **A power discrepancy that must be resolved first.** LEXOR L4's gate is ≥30 provenance-independent
  components. §2.7 requires ≥600 for a containment gate and ≈423 for a paired ranking contrast. These
  differ by an order of magnitude. Part of the gap is estimand, but the Metz record already documents that
  MDE computed from retraining noise is optimistic by ~1.8× in spread once a *different arm* adds
  per-component heterogeneity. **L4's threshold should be re-derived from paired arm-contrast dispersion
  before any extraction budget is spent**, otherwise LEXOR could pass L0–L4 at real cost and still be
  unable to resolve K-LBP at L5.

---

## 8. What this revision removes from the blueprint

* PARC's MLP correction `Delta = W_2 sigma(W_1 u_t)` (630 params, ~7.1 per target) — **retired** as
  over-parameterized; replaced by the 30-parameter rank-1 gated form.
* `v` in the proposed correction — **removed** as redundant for a scalar prediction.
* `b` acting on the raw 1034-d ligand descriptor — **replaced** by `c` on the `r`-dimensional shared basis.
* Any use of evidence confidence, `n_sources` or contradiction counts as predictive features —
  **prohibited**.

No new module beyond the single rank-1 biological correction is added. Nothing in Part 3's closed-route
list is reopened: K-LBP uses no pose, complex, docking score, contact map, interaction field, atom–residue
network, cross-attention, MoE, GNN, transformer, taxonomy label, assay id or provenance id.
