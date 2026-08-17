# Stage T preregistration — true-MMP transformation space

Frozen **before any census statistic, any transformation-key label analysis and
any trained-arm evaluation metric was computed.** Only two things were known
when this file was written: (a) the pre-aggregation provenance artifacts exist
and hash-match the corpus manifest, and (b) the counts of repeated-measure
groups in `meta_train` (233 groups with >1 raw row: 136 cross-panel, 97
within-panel). Those are provenance counts, not labels and not results; they
are recorded here so the T0 design is grounded rather than speculative.

Nothing in this file may be changed after the first T1 admission statistic or
the first T2 evaluation metric is read. A post-hoc threshold change voids the
stage.

## 0. Estimand and why this is not Stage S

Stage S rejected `phi(L) + global alpha(P) -> FiLM potential` over
**whole-molecule** pairs: correct and shuffled proteins produced the same
apparent gain (+0.0065 vs +0.0085 Pearson), the protein channel was learned as a
target-identity key (+0.109 in-distribution, +0.0065 out of component, identical
under permutation), and only the ligand-only signed-SAR field retained signal.
`method_ladder/CLOSURE_MAP.md` family 8 records that **pairs were never
matched-molecular-pair identified**, so the MMP constraint of that family has
never been instantiated.

This stage tests a **different estimand**. Write

    delta_y(t, tau) = mu_tau + delta(t, tau) + noise

for a directed transformation `tau` applied in target `t`. `mu_tau` is the
generic chemical effect; `delta(t, tau)` is the target-specific response. The
decisive quantity is the **crossed double difference**

    D(tau, t1, t2) = delta_y(t1, tau) - delta_y(t2, tau)

which removes the target-level affinity offset **and** `mu_tau` exactly. The
question is whether protein biology predicts `D` on unseen protein components
better than shuffled and matched-wrong proteins. A model that scores on `D`
cannot be scoring on target level, on generic medicinal chemistry, or on a
target key — all three cancel algebraically.

This is an admissibility test, not a performance claim. It is the final one for
the protein-conditioned SAR-latent-space direction.

## 1. Governance

* Governed BindingDB-Ki corpus `main_v0`, double-cold split `v1`, mounted for
  all model and episode work through the **physically isolated split view**
  (`bindingdb_ki_double_cold_v1_views`). The `meta_test` label artifact is not
  present on that surface.
* `meta_val` is used for nothing: not training, not selection, not thresholds,
  not early stopping, not a census. A parsed-AST test fails the suite if its
  name appears as a string constant in any stage module.
* Frozen `scripts/internal_validation.py` partition of `meta_train`
  (`PARTITION_SEED = 20260818`): 227 fit / 31 internal-validation components.
  Training reads fit only; reporting reads internal-validation only.
* Stable SHA-256 seeds throughout (`scripts.qpsmp_data.stable_seed`). Python
  `hash()` is forbidden and tested for by AST parse.
* Labels never touch feature construction, split construction, MMP definition,
  transformation-key granularity or hyperparameter selection. Every rule below
  that could in principle peek is stated as a **label-blind** rule.

### Disclosed isolation exception for the T0 provenance step

T0 must read pre-aggregation measurements, which live in a single all-label
artifact (`dataset/processed/crossed_interaction/bindingdb_202608/
exact_labels.jsonl.gz`, SHA-256 `6a44a151…c823dcc`, which matches the corpus
manifest's `labels_sha256`; and `metadata_projection.jsonl.gz`, matching
`projection_sha256`).

Procedure: the allow-list of `source_row_id`s is derived **from the physically
isolated split view, `meta_train` only**; the label artifact is then streamed
and a row is retained **in the same pass** only if its id is allow-listed.
Non-allow-listed values are never bound to any retained structure. The result is
written once to a `meta_train`-only provenance cache and **every downstream
module reads only that cache.**

This is **logical exclusion after parsing**, which is strictly weaker than the
model path's physical isolation, and it is disclosed as such in
`T0_RELIABILITY.json`. `meta_val` and `meta_test` rows are both excluded by the
allow-list. No `meta_val` or `meta_test` value enters any statistic in this
stage.

## 2. Stage T0 — measurement-reliability audit

**Purpose: measure the reliability of the supervision. This is not, and will not
be reported as, a universal benchmark MSE floor.**

Aggregation rule actually used by the corpus (from `manifest.json.cleaning`,
quoted verbatim): `within_panel_aggregation: "median"`, then
`cross_panel_pair_aggregation: "equal-panel median"`, endpoint
`"exact positive uncensored Ki"`, transform `pKi = 9 - log10(Ki[nM])`,
`panel_id_role: "BindingDB assay proxy"`.

Three provenance levels, from `metadata_projection` + `exact_labels`:

| level | definition | interpretation |
|---|---|---|
| L1 | rows sharing `(panel_id, assay protocol_sha256)` | technical / same-protocol repeat |
| L2 | rows sharing `panel_id`, differing protocol | within document+endpoint+target, different assay |
| L3 | rows of one `(target, ligand)` across `panel_id`s | between assay / between document |

Reported: counts and coverage at each level; robust dispersion (median absolute
deviation and the interquartile range of centred residuals, alongside the sd) of
raw `pK` within each level; and the implied variance components.

**Uncertainty of a difference label under the actual aggregation rule.** For a
pair `(a, b)` of ligands in one target, `delta_y = y_b - y_a` where each `y` is
the aggregated median. The reported quantity is `sigma^2_Delta` decomposed into
the part that cancels when both cells share a panel and the part that does not.

**Frozen, label-blind confidence strata** (defined by provenance only, so no
label value can move a pair between strata):

* `S1 same_panel_single` — the two cells share a panel and each has
  `panel_count == 1`;
* `S2 same_panel_multi` — the two cells share a panel, at least one has
  `panel_count > 1`;
* `S3 cross_panel` — no shared panel. **Weak/noise stratum. Never pooled into
  the primary bank.**

Uncertainty weights, if used, are `w = 1 / sigma^2_Delta(stratum)` with the
`sigma^2` taken from T0 and **fixed before T2**.

**If the repeated-measure subset is too small or too selected to identify a
component, T0 must say the quantity is not identifiable and no number is
invented.** T0 must also report the selection bias of the repeated-measure
subset (which compounds get measured twice are not a random sample) and the
admission-time filtering the corpus already applied
(`admission.conflicting_ligands = 333`, `cross_panel_pairs = 618`).

## 3. Stage T1 — true-MMP census

Transformations are built with **RDKit's supported MMPA machinery**
(`rdkit.Chem.rdMMPA.FragmentMol`, Hussain–Rea), single cut
(`minCuts=1, maxCuts=1`), default bond pattern, `resultsAsMols=False`,
isomeric SMILES throughout. No ad-hoc SMILES string manipulation anywhere.

Label-blind construction rules, all frozen here:

* **core / R assignment**: of the two fragments, the one with more heavy atoms
  is the core; ties broken by canonical SMILES sort;
* **pair admissibility**: same target, same Ki endpoint, identical core SMILES
  (including `[*:1]` and stereochemistry), different R;
* **attachment context**: the core atom bearing `[*:1]`, recorded as
  `(element, aromatic, in_ring, degree, formal_charge)`;
* **stereochemistry**: isomeric SMILES retained on core and both R groups; a
  transformation whose only difference is stereochemical is flagged
  `stereo_edit = True` rather than collapsed;
* **charge**: `formal_charge(R_b) - formal_charge(R_a)` recorded on every
  transformation;
* **canonical direction**: the pair is ordered by canonical SMILES sort of
  `(R_a, R_b)`; `delta_y` is taken in that direction. This is a function of
  structure only, so no label can choose the direction. The inverse
  transformation maps to `-delta_y` by construction, and a test asserts it;
* **deduplication**: a repeated `(target, core, R_a, R_b)` observation is
  reduced to one row by **corpus order of the lower cell index** — never by
  label value.

**Two transformation-key granularities, both fixed now, neither chosen on
results:**

1. **exact key** — `(core-attachment context, R_a isomeric SMILES,
   R_b isomeric SMILES)`;
2. **coarse key** — `(reduced attachment context (element, aromatic),
   R_a canonical SMILES with stereochemistry removed, R_b likewise)`.

**Primary bank requires an overlapping identical governed `panel_id`.**
Cross-panel MMPs are built and reported as the S3 weak stratum and are never
silently pooled.

Bipartite evidence graph `G = (targets, transformation keys, observations)`.
Reported: total observations; targets and protein components; key counts at both
granularities; target-degree and component-degree per key; keys with >= 3 targets
and >= 3 components; connected components of `G`; same-panel vs cross-panel;
scaffold and transformation overlap between fit and internal populations;
activity-cliff and stereo-edit counts; T0 uncertainty per observation; and
degree concentration (top-key share) so that a handful of transformations cannot
carry the stage.

### Frozen T1 admission thresholds — verbatim

1. at least **2,000** same-panel fit observations;
2. at least **50** fit targets;
3. at least **30** transformation keys each spanning **>= 3 targets and
   >= 3 components**;
4. at least **300** internal observations;
5. at least **10** internal protein components.

Thresholds 1–3 are evaluated on the **exact** key. Threshold 3 is additionally
reported on the coarse key; a coarse-only pass does **not** admit T2, it is
recorded as a coverage observation.

**If any threshold fails, stop immediately.** Write a negative `REPORT.md`. Do
not train T2 or any neural MMP model.

## 4. Deployment-coverage audit

On the frozen, label-blind nested `k = {1,2,3,5}` episode banks
(`QPSMPData.fixed_nested_episode_banks`, the same construction the governed
protocol uses), estimate

    C_k = P(at least one support-query pair forms a valid MMP)

reported separately for the exact and coarse keys, and stratified by `k`,
protein component and ligand novelty. Query labels are not read.

If `C_k` is low, the report must state that MMP may be usable as a **training
signal** only, and cannot serve as a universal reference-based inference
mechanism.

## 5. Stage T2 — crossed-double-difference identifiability

Runs **only if T1 passes.**

`delta_y(t, tau)` is aggregated robustly (median) within compatible panels,
carrying the T0 stratum. Double differences `D(tau, t1, t2)` are formed **only
within one transformation key** — never across incompatible keys.

Discriminator: small, ordinary gradient-trained. Inputs are
**frozen protein embeddings / residue summaries** and **structured
transformation descriptors** derived from core, removed fragment, added
fragment and attachment context. Explicitly forbidden as model input: target ID
embedding, document, assay, panel, target index, component ID. Forbidden
methods: ridge, analytic solver, pseudoinverse, closed-form estimator,
test-time label gradient.

**Structural guarantee.** The protein-pair prediction is a difference of a
per-protein response:

    D_hat(tau, p1, p2) = R(tau, p1) - R(tau, p2)

so identity (`p1 = p2` gives exactly 0), antisymmetry and protein-cycle
consistency hold for every parameter setting. Tests assert this before training.

Two evaluation surfaces:

1. **protein-component cold** with repeated transformation keys;
2. **protein-component cold + transformation-disjoint** (or
   transformation-family-cold), run only if T1 proves the coverage for it.

Splits are made before training, never place related protein components on both
sides, and keep every instance of a transformation key (or declared family)
together.

### Matched arms — same fixed budget

| arm | inputs |
|---|---|
| **A** | zero / generic-transformation baseline |
| **B** | transformation descriptor, no protein biology |
| **C** | correct protein + transformation |
| **D** | stable shuffled-protein control |
| **E** | similarity-matched wrong-protein control |
| **F** | within-transformation shuffled double-difference labels |

**No "correct beats wrong" counterfactual loss is trained.** Stage S measured
that such a loss passes by making the wrong branch explode (hard-wrong MSE
12.62 against its own correct 4.09). Wrong proteins are **evaluation controls,
not optimization targets.**

Matched-wrong proteins: a different CD-HIT40 component; selected on
`meta_train`-only frozen protein similarity; as similar as feasible in PLM
space; no shared assay/document programme; **only the protein input is
replaced.**

## 6. Statistics

Reported: MSE and MAE of `D`; Pearson and Spearman; sign accuracy and
CI/concordance where meaningful; per-transformation and per-component results;
same-panel uncertainty strata (S1/S2/S3); transformation-frequency strata;
chemical and protein novelty strata.

**Rows sharing a target or a transformation are correlated, so raw-row
bootstrap is forbidden.** The preregistered estimator is a **two-way
(multiway) cluster bootstrap** over protein components and transformation keys:
each draw samples components with replacement and transformation keys with
replacement independently, and a row's multiplicity is the product of its
component draw count and its key draw count. Both arms of a contrast are
re-scored on the identical draw. 2,000 draws, seed 20260820. **Effective
independent units** are reported as `min(#components, #transformation keys)`
alongside every interval.

## 7. Frozen T2 gate — verbatim

The route passes only if **all** hold:

1. correct protein minus shuffled protein Pearson **>= +0.05**;
2. its clustered-bootstrap 95% lower bound is **> 0**;
3. correct protein minus matched-wrong protein Pearson **>= +0.05**;
4. its 95% lower bound is **> 0**;
5. correct protein improves both error and ranking/sign metrics;
6. label shuffle destroys the effect;
7. the result is not confined to one transformation, one component or
   high-frequency edits;
8. the correct protein's output change is **aligned with truth, not merely
   large** (the Stage P / Stage S failure mode: measured as the correlation
   between the protein-induced shift and the truth residual);
9. no target-key shortcut is reproduced on a fit-unsampled bank;
10. the protein/transformation-cold result does not reverse the
    protein-only-cold result.

## 8. Stop rule

**If T2 fails, protein-conditioned SAR latent space is formally stopped under
the current BindingDB protocol.** It is not rescued by DrugBAN, PSICHIC-style
attention, Cartesian tensors, conformers, MSA, more capacity, more seeds,
meta-learning or threshold changes.

If T2 passes, nothing is integrated into production. The output is a proposal
for the next isolated stage: transformation/edit tokens querying residue-level
protein regions with **no pooled-protein bypass**. Cartesian ligand geometry
stays unauthorized until a residual audit shows a specific stereochemical or
conformational failure.

## 9. What each stage does and does not establish

* **T0** measures supervision reliability. It is **not** a universal MSE floor.
* **T1** measures whether the transformation graph is identifiable at all.
* **T2** measures protein x transformation interaction.
* None of the three, alone or together, is a zero-shot or few-shot DTA
  performance claim.
* One seed may reject a hypothesis; it cannot establish final performance.

## 10. Verification required before training

Deterministic MMP decomposition; inverse transformation and sign consistency;
attachment and stereochemistry preservation; no cross-target or cross-panel
contamination; the physical meta-test seal; stable banks across
`PYTHONHASHSEED` values; antisymmetry, identity and protein-cycle consistency;
no dead trainable parameters; no evaluation-label path into inputs or splitting.
