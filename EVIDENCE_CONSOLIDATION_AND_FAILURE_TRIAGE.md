# MetaSieve Evidence Consolidation And Failure Triage

Updated: 2026-08-08

## 0-PRE-2. Structural self-supervision programme (2026-08-08)

```text
POSE_FREE_DEPLOYMENT_INPUTS_INSUFFICIENT
```

The only untested information rung was genuine 3D structural mechanism. S1 formed
an independent zero-exposure RCSB CC0 test set (1,118 complexes, 621 protein
clusters, 586 scaffolds). S2 froze a six-channel coordinate teacher, reproducible
at machine precision with no degenerate channel. S3 showed the design detects
`R2 >= 0.02` at 100% power. S4 found two channels observable from sequence+2D and
beating random features, but **no channel beats the deranged-protein control**.
S4b localised the entire observable signal to the **ligand**: protein-only
`R2 ~ 0` on all six channels, and adding the protein to the ligand costs
`-0.032` on hydrophobic burial with a CI excluding zero. GPU training was not
authorised because the distillation target is a ligand descriptor. See
`report/ssl_b2_structural_observability/S_PROGRAMME_REPORT.md`.

## 0-PRE. Terminal outcome of the external programme (2026-08-08)

```text
PUBLIC_DATA_INSUFFICIENT_FOR_IDENTIFICATION
DEPLOYMENT_INPUTS_INSUFFICIENT
```

XP3 quantified the governing trade-off: single-lab profiling panels have real
reproducible interaction but only ~8 protein-group independence units (kinase
taxonomy cap); literature aggregations have 70-85 units but their interaction
sits below the measurement noise floor. XP4 confirmed the latter directly
(`sigma=0.777`, `gamma` sd `0.406`, chemistry ceiling `R2=-0.539`). XP5 closed
the pose-free typed rung on the panel that does have signal (derangement
specificity exactly `+0.00000`). Full account:
`report/FINAL_IDENTIFIABILITY_REPORT.md`.

## 0. XP1 / XP2 External Crossed-Panel Programme (2026-08-08)

```text
XP1: OBJECTIVE_OR_PARAMETERIZATION_FAILURE
XP2: CROSSED_INTERACTION_REPRODUCED
     K_LE_5_SECTION_NOT_IDENTIFIED
     PANEL_LOCAL_LOW_RANK_META_LEARNING
     BIOLOGICAL_LANDING_NOT_IDENTIFIED
     EXTERNAL_REPLICATION_FAILED
```

XP1 established on external release-pinned panels that a protein-by-ligand
interaction exists, is large (`59.6%` of panel variance), low-rank, reproducible
across disjoint compound halves (`r = 0.885`) and across an independent platform
(`r = 0.565`), and that a support-identified section transfers to unseen kinase
groups at `k = 16` while no zero-shot protein representation does.

XP2 then tested whether that section is deployable and **closed the mechanism**:

- XP2-A reproduced XP1 from immutable artifacts (18/18 checks) and corrected
  XP1's single-floor censoring *description* without changing its analysis set.
- XP2-B **passed**: the ligand loading transfers to unseen Bemis-Murcko
  scaffolds (`R2 = +0.199 [+0.133, +0.261]` vs random `+0.025`), so `u(L)` is
  not a lookup table.
- XP2-C: the identifiable section dimension is exactly `min(k-1, d)` — measured
  `0,1,2,3,3` at `k = 1..5` — so a frozen `k <= 5` caps it at four and gives
  **zero** at `k = 1`; magnitude never exceeds `R2_gamma = 0.025` against a
  frozen `0.05` floor.
- XP2-D: under simultaneous protein-group and ligand-scaffold closure the
  derangement, permutation, zero-adaptation and random-correction controls all
  have intervals containing zero, and **random ligand features reproduce the
  entire remaining gain**. XP1's specificity was conditional on ligand reuse.
- XP2-E: every zero-shot protein representation lies within `0.002` of a random
  protein embedding.
- XP2-F: direction transfer to Klaeger kinobeads fails
  (`Delta_interaction = +0.00346 [−0.00234, +0.00863]`).

Nothing was promoted. `model/`, production `scripts/`, `contracts/` and
`theory/` were not modified; CSMO, Band, `K` and the mesh are untouched; the
regression suite is `73 passed` before and after. Full record in `history.md`,
`report/crossed_panel_identification/XP1_RESEARCH_REPORT.md` and
`report/crossed_panel_deployability/XP2_FINAL_REPORT.md`.

The probability-law operator `K(B(z)F(z))` was **never scored** in either stage.
Interaction `R2` is not law calibration.

## 1. Current Scientific Verdict

```text
KINASE_PANEL_COMPONENT_IDENTIFIABILITY_OBSERVED_IN_DEVELOPMENT
F6I_TOTAL_GATE_NOT_ADMISSIBLE
FRESH_ENDPOINT_CONSISTENT_EXTERNAL_ADMISSION_NOT_RUN
WITHIN_TASK_RANKING_DIRECTION_NOT_IDENTIFIED
PROTEIN_SPECIFIC_AFFINITY_LOCATION_NOT_YET_TESTED
CROSSED_INTERACTION_EXISTENCE_NOT_YET_TESTED
AFFINITY_ENERGETICS_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
NO_VALIDATED_END_TO_END_DTA_MODEL
```

Retained from earlier consolidations, unchanged:

```text
MATHEMATICAL_OPERATOR_IMPLEMENTED_AND_CONTRACT_TESTED
GEOMETRY_IDENTIFIED
PAIR_COMPATIBILITY_IDENTIFIED
FIXED_RADIAL_BASIS_PARTNER_RECOVERABILITY_IDENTIFIED
TASK_LOCAL_RADIAL_AFFINITY_HEADROOM_OBSERVED
CORRECT_PARTNER_AFFINITY_SECTION_NOT_IDENTIFIED
FIXED_RADIAL_INTERACTION_RESIDUAL_NOT_OBSERVED
CHEMBL_CROSSED_SOURCE_INTERACTION_UNDERDETERMINED
PKIS_ZERO_SHOT_SURFACE_PLUS_ONE_DIMENSIONAL_LOCATION_OBSERVED
```

### 1.1 Two separate estimands

**Claim A** — does the correct protein provide affinity-*location* information
beyond population, ligand-only and protein-sequence-only baselines? Addressed by
`E-AFF-L0`.

**Claim B** — does a non-additive protein-by-ligand affinity *interaction*
exist? Addressed by `E-AFF-X0-B` and then `X1`.

**L0 and X1 are different estimands.** Neither replaces, implies nor authorizes
the other.

### 1.2 Scope of the historical affinity negatives

The P1C/P1R\*/E-AFF-P0/H0A/H0C evaluations primarily used within-task
concordance or another rank-based metric. Those experiments remain valid and
every verdict and metric is retained unchanged. Their negative conclusions are
limited to **within-task ranking information**, because within-task concordance
is invariant to task-wise prediction shifts and positive rescaling. They did not
directly test whether the correct protein provides transferable
affinity-location information. No historical FAIL, CLOSED, NOT-RUN or negative
result is reclassified as PASS.

The frozen theory and the bioinformatics code are **interface-compatible but
not deeply integrated**. The implementation realizes the mathematical output
chain

```text
z -> F(z) in simplex -> B(z)F(z) -> K(beta)
```

but the biological coordinates supplied to `z` have not passed the required
protein-specific affinity gate. The former 28-dimensional state used arbitrary
bounded projections of protein, ligand, pair and support latents. P1B's validated
contact/distance bridge was deliberately not connected to that state. Therefore
the mathematics constrains the output object, while biology has not yet determined
an identified statistic coordinate system.

## 2. Theory-To-Code Evidence Map

| Frozen object | Retained code | Evidence status | Triage |
|---|---|---|---|
| valid CDF-band polytope | `model/bands.py` | contract tests pass | retain |
| `K(beta)` and Hausdorff-W1 | `model/mathematical.py` | exact/finite-grid tests pass | retain |
| simplex coefficient map `F(z)` | `model/meta_operator.py` | simplex invariants and frozen-HN embedding tested | retain |
| positive ridge and band loss | `model/mathematical.py` | theory-aligned unit tests pass | retain |
| `B(z)` assembly | `model/meta_operator.py` | validity and manifest checks pass | retain |
| observable bounded `z(S,Q,gamma)` | no admitted production implementation | previous arbitrary latent state failed biological gates | remove failed implementation |
| protein/ligand local states | `model/encoders.py` | used by P1B geometry PASS | retain |
| contact/distance geometry | `model/mechanism.py` | P1B correct-partner Gate PASS | retain |
| fixed privileged radial basis | `history.md` | T-BASIS-R0 partner-recoverability PASS | historical development evidence |
| population-shared basis-to-affinity calibration | `history.md` | E-AFF-P0 negative | failure record only |
| task-local basis-to-affinity calibration | `history.md` | headroom positive; partner margin failed | mixed record only |
| centered radial interaction residual | `history.md` | no increment over support-matched ligand nuisance | failure record only |
| crossed affinity interaction data | `history.md` | X0 effective components insufficient before label read | source-data stop |
| X0 independence unit and its 245 requirement | `history.md` | X0-FEAS: requirement exceeds the unit's own ceiling | specification-induced stop |
| cell-disjoint crossed design | `history.md` | X0-B: conditional design support pending `rho` | historical design record |
| within-task rank readout | `history.md` | R0: exactly invariant to task-level affinity location | scope annotation only |
| ordered-anchor location interface | `history.md` | L0 not run, precondition failed | failure record only |
| PKIS/KLIFS zero-shot biological surface | `history.md` | component-level development signal on consumed panels | consumed development evidence |
| protein-independent support location `tau` | `model/component_statistic.py` | bounded one-dimensional permutation-invariant algebra | internal, not exported |
| seven-point categorical law bridge | `history.md` | local invariants passed; not type-equivalent to frozen production CDF bands | removed |

The theory proves properties conditional on a legal statistic and coefficient
map. It does not prove that ESM/GINE features, QPMA, CSMO views, or any proposed
biological coordinate contain affinity information. Those are empirical claims
and remain gated separately.

## 3. Consolidated Evidence Ledger

### Passed and retained

| Stage | Established fact | Retained evidence |
|---|---|---|
| P0 | canonical DTA/data contracts | `scripts/data*.py`, sealed-data tests |
| P1A | governed open holo corpus | structure acquisition/index/governance scripts |
| P1B | correct-partner contact/distance geometry | P1B checkpoint and final Gate artifacts |
| D0-C | release-pinned ChEMBL37 Ki/Kd corpus | static-release extractor and D0 report |
| D1 | homology/document closure | governance scripts and frozen splits |
| E0R1-C | synthetic design transports the teacher | research-only Moore-Penrose witness |

### Failed and removed after consolidation

| Stage | Terminal evidence | Historical conclusion |
|---|---|---|
| legacy biological frontend / Phase 3-4 | support/protein interaction not identified | no validated assembled model |
| P1C | correct and deranged affinity CI overlap | geometry is not affinity semantics |
| P1R0/R1 | old readout invalid; invariant MIF still below Gate | readout repair insufficient |
| P1R2A | interaction variance `1.67%`; affinity increment negative | non-additive MIF not affinity-incremental |
| P1R2B0 | Ridge/spline/MLP did not rescue MIF | global readout capacity not the bottleneck |
| P1R2B1 | wrong-protein penalty without positive ligand increment | compatibility is not affinity direction |
| F0/F0R | live API could not reproduce immutable rows | historical API route closed |
| E0 synthetic | correct-partner recovery below frozen Gate | real affinity training never authorized |
| E0R0 learned heads | analytic realization passed; learned heads failed | representation sufficiency separated from training |

### Synthetic numerical boundary closed in research

E0R1 confirmed an objective-semantics defect and obtained an exact pinv witness.
E0R2 then used a separately registered float64 augmented least-squares solve of
the corrected point-residual plus residual-difference objective. It reached
train RMSE `3.188e-8`, full-gradient L2 `6.150e-17`, and historical
correct/deranged CI `0.99737/0.61447`. Verdict:
`SYNTHETIC_OBJECTIVE_DESIGN_SOLVER_IDENTIFIED`.

This closes only the synthetic realization/solver boundary. Directional/type
identifiability, reference-state structural log-odds, real affinity calibration,
few-shot mechanism adaptation, production admission and P2-P4 remain untested
or frozen.

### Structural typed-interaction feasibility did not replicate

T-DIR-P0 used PLIP weak labels and frozen P1B features on 40 label-blind selected
complexes. Direct hydrophobic/H-bond mappings were operational, but the primary
hydrophobic D2 arm overfit: train/validation/test AP was
`0.25496/0.01699/0.01996`. It failed the frozen held-out feasibility criteria.
The stronger low-dimensional D1 descriptive result cannot be promoted because
the pilot is small, arms are not capacity matched, and required shuffle and
chemistry-roundtrip controls were omitted. Group-based salt/pi channels also
had a post-run mapping defect and are uninterpretable.

The strongest allowed conclusion is
`PILOT_LEARNABILITY_SIGNAL_NOT_OBSERVED`. This is not evidence that typed
interactions are impossible; it is evidence that the tested high-dimensional
frozen local-state probe is not ready for a formal Gate.

### Fixed radial basis recovery passes on a fresh panel

T-BASIS-R0 replaced sparse named events with a fixed analytic 288D two-body
chemogeometric radial teacher. A small shared radial calibration recovered the
teacher from frozen sequence+2D P1B distance distributions. On 64 held-out test
complexes, reconstruction gain was `0.5312 [0.4433,0.5962]`; replacing the
protein increased error by `0.1561 [0.1070,0.2007]`. The test panel excluded
all T-DIR-P0 records and was held out from P1B training.

This admits the two-body radial basis as a research-stage structural object,
not as biological `z`. Partner degradation still mixes protein composition and
geometry. Angular/many-body recoverability and affinity direction remain
separate untested boundaries.

### Affinity falsification separates headroom from partner specificity

E-AFF-P0 fitted one shared 288D direction using task-balanced residual
differences and a closure-OOF ligand prior. Across 245 closure components,
correct-minus-ligand was `-0.01016 [-0.02069,-0.00018]` and
correct-minus-deranged was `-0.00236 [-0.01301,0.00805]`. A universal radial
affinity direction was therefore not observed.

E-AFF-H0A then fitted independent task-local directions on 20 ligands and
tested each on 20 untouched ligands in 107 closure components. The resulting
correct-minus-ligand gain was `+0.08821 [0.06761,0.10998]`, proving that the
fixed tensor carries substantial task/series-local ranking information.
However, correct-minus-deranged was only
`+0.00864 [0.00338,0.01462]`, below the preregistered `+0.03` partner margin.
The same head retained nearly all of its value after replacing the protein.

The current data therefore support neither "the radial basis has no affinity
information" nor "a target-specific radial section is identified." They
support a narrower conclusion: task-local affinity headroom exists, while its
correct-protein-specific component remains insufficient. H0-B and RFSA are not
authorized, and adding orientation is not automatically justified.

### Shortcut removal does not recover radial interaction affinity

E-AFF-H0C excluded every H0A task and selected 54 new tasks/components with
strict 20/20 Murcko scaffold-disjoint support/test partitions. A frozen 128D
ligand-state nuisance used the same 20 support labels. The interaction head then
received only cross-fitted nuisance residual differences and the fixed
double-centered tensor `psi=(phi-phi_null)/total`.

Global-L, Local-L, correct and deranged CI was
`0.55487/0.59635/0.59244/0.59092`. Thus task-local ligand adaptation added
`+0.04147`, while the correct radial residual changed Local-L by
`-0.00391 [-0.02040,0.01191]`. Correct-minus-deranged was only
`+0.00152 [-0.01024,0.01424]`. Verdict:
`FIXED_RADIAL_INTERACTION_RESIDUAL_NOT_OBSERVED`.

This localizes H0A's headroom primarily to ligand/series adaptation and closes
the proposed fixed-basis shortcut-removal rescue. It still does not distinguish
missing biological coordinates from insufficient source interaction labels.
That distinction requires a real target-by-ligand double-difference census and
Gate; it cannot be supplied by RFSA, derangement training, or the frozen
mathematical operator.

### Crossed ChEMBL evidence is underdetermined after dependency closure

E-AFF-X0 queried only label-blind metadata for the 152,737 governed source
activity IDs. It formed target-independent structured panel contexts and
counted `2 targets x 2 ligands` rectangles without selecting any affinity value
field.

Ki had 1,059,169 nominal rectangles across 597 panels, but shared panel and D1
homology-document closure reduced these to 36 effective components. Kd had
232,875 nominal rectangles across 34 panels but only 12 effective components.
Replicate-supported counts were 18 and 4. Both are far below the frozen lower-
bound requirement of 245 independent units. The largest dependency component
alone contains 56.49% of Ki and 76.76% of Kd rectangles.

Verdict: `STOP_SOURCE_INTERACTION_UNDERDETERMINED`. This is a data/estimand
stop, not evidence that real affinity interaction is absent and not evidence
that radial or angular biology is sufficient. X1 and X2 are prohibited because
opening labels cannot repair missing independent crossed support.

### The X0 stop was specification-induced

E-AFF-X0-FEAS audited the X0 estimand rather than the source, label-blind. A
rectangle requires two proteins inside one document-keyed panel, and D1 closure
unions every pair of targets sharing a document, so both proteins of every
rectangle always already lie in one closure component. Predicted panels touching
more than one closure component: `0`; observed: `0`. Crossing can therefore
never create an X0 unit, and effective components are bounded by the
closure-component universe.

That universe is `245` for the entire governed corpus, against a frozen
requirement of `245` **per endpoint**. Only `202` components carry Ki rows and
`72` carry Kd rows. Recomputing the closure over the full governed D0 corpus and
over shallower task populations raises the best ceiling only to Ki `97` and Kd
`56`. The frozen `245` was independently re-derived from its stated one-sided
chi-square design and is arithmetically correct; the defect is that an IID
sample-size formula was bound to a unit whose universe is `245` in total.

Verdict: `X0_UNIT_REQUIREMENT_UNATTAINABLE_BY_CONSTRUCTION`. The X0 verdict,
metrics and artifacts are retained unchanged; this is a scope annotation, not a
reclassification. It also means acquiring a more genuinely crossed
source/selectivity corpus does not repair the census under that unit, since
crossing and document-disjointness are produced by opposite kinds of study.

### X0-B is conditional design support only

E-AFF-X0-B re-registered the independence unit as the **cell-disjoint
rectangle** and the effective sample size as a design-effect adjusted count,
leaving the effect size `0.5`, variance ratio `1.25`, alpha `0.05`, power
`0.80`, the `245` requirement and the `+0.03` affinity Gate margins unchanged.

Ki packs `11,168` units across `36` clusters from `205` distinct target pairs;
Kd packs `1,041` across `12` clusters from `49` pairs. With `n_eff = N/DEFF` and
`DEFF = 1 + (m_A - 1)rho`, the design meets `245` only if intra-cluster
correlation is at most `rho* = 0.0915` (Ki) or `0.0164` (Kd); the hard bound
`n_eff <= G/rho` makes cluster count, not rectangle count, the binding resource.
At `rho = 1` the optimal design collapses to one unit per cluster and reproduces
X0's `36` and `12` exactly, so X0 was this model's total-correlation corner.

This is a statement about achievable effective sample size. It is **not**
evidence that protein-by-ligand affinity interaction exists. X1 remains
unauthorized until `rho` is estimated with uncertainty and its upper confidence
bound is gated against `rho*`.

### The affinity readout was scoped to ranking only

E-AFF-R0 diagnosed the instrument rather than the biology, reading no affinity
label. Against the repository's own `metrics.concordance` on the real H0C
partition, within-task concordance changed by `0.0` under per-task prediction
shift, prediction rescale, label shift and label rescale — exactly invariant, not
approximately. A simulated predictor holding a task's affinity level perfectly
and nothing else scores exactly `0.5000` at every variance share tested, up to
`0.985`, while its advantage under a location-sensitive error grows from
`1.033/1.001` to `8.832/1.017`.

H0C additionally removed that channel upstream: the geometry was shown
`y - global_ligand_prior - task_local_ligand_nuisance`, where the nuisance is
fitted on 20 labelled supports of the correct protein's own task and the result
is added to **both** the correct and the deranged arm. Both arms therefore held
the correct protein's task level before the contrast was taken.

R0 does **not** show that protein-specific affinity lives in the location
channel. It shows that if it does, nothing in the evidence chain could have
detected it. All prior verdicts stand as statements about within-task ranking.

## 4. Failure Localization

```text
sequence + ligand graph
        |
        v
local encoders ------------------------------- PASS
        |
        v
correct-partner contact/distance geometry ---- PASS (P1B)
        |
        +--> fixed radial basis recovery ------- PASS (T-BASIS-R0)
        |
        +--> angular/many-body basis ----------- NOT TESTED
        |
        v
protein-specific affinity direction ---------- NOT IDENTIFIED
        |
        +--> shared radial direction ---------- NOT OBSERVED (E-AFF-P0)
        |
        +--> task-local radial headroom -------- OBSERVED (E-AFF-H0A)
        |
        +--> correct-partner local section ----- BELOW GATE (E-AFF-H0A)
        |
        +--> centered radial interaction ------- NOT OBSERVED (E-AFF-H0C)
        |
        +--> real crossed source support -------- UNDERDETERMINED (E-AFF-X0)
        |         (specification-induced; E-AFF-X0-FEAS)
        |
        +--> crossed design under cell-disjoint - CONDITIONAL ON rho (E-AFF-X0-B)
        |
        +--> protein affinity LOCATION ---------- NOT YET TESTED (E-AFF-L0 registered)
        |         every result above is rank-scoped (E-AFF-R0)
        |
        v
biological statistic admitted to z ----------- BLOCKED
        |
        v
A(F,z) end-to-end biological claim ------------ NOT ESTABLISHED
```

The earliest unresolved boundary is not the frozen operator. It is the map from
validated local geometry to a signed, transferable affinity statistic. No change
to CSMO, Band, simplex, positive ridge, `K`, or mesh `h` is justified by current
evidence.

## 5. Repository Disposition

### Production/verified surface

- `model/`: frozen mathematical primitives plus P1B local encoders/bridge only.
- `scripts/`: successful data, structure-geometry and release-governance workflows.
- `contracts/`: active data, ligand and P1B mechanism contracts.
- `report/`: current split protocol and immutable PASS evidence only.

### Research surface

- `research/README.md` is now only the boundary for future preregistered work.
- Terminal synthetic, negative, mixed and consumed-development implementations
  were removed after consolidation into `history.md`.
- The exact deleted tree remains recoverable from Git commit `8b7789e`.

### Canonical evidence index

| Evidence | Canonical location | Retention reason |
|---|---|---|
| Frozen mathematics | `theory/FINAL_FROZEN_THEORY/` | authoritative theorem and scope contract |
| P1B geometry PASS | `report/mechanism_refactor/p1b_*` | accepted checkpoint and partner-control evidence |
| D0-C/D1 PASS | `report/mechanism_refactor/p1r2b_d0_chembl37_v1/` | immutable corpus and closure governance |
| E0 through X0 research | `history.md`, Git `8b7789e` | consolidated terminal result and recoverable deleted artifacts |
| PKIS/F6I component research | `history.md`, Git `8b7789e` | component-level development evidence and exact deleted sources |
| Consolidated report | `report/VERIFIED_EVIDENCE_SUMMARY.md` | retained PASS evidence and current boundary |
| Current scientific summary | this file | single theory-to-evidence interpretation |
| Full historical ledger | `history.md` | negative results, supersession and deletion record |

### Deleted failed/superseded surfaces

The following categories are removed only after their conclusions were copied to
this file and `history.md`:

- failed P1C/P1R*/F0R scripts and their dedicated tests;
- the failed assembled biological pipeline and legacy trainer/evaluator;
- superseded smoke, invalid-numeric and failed-gate artifact directories;
- old Phase X/Y/Z reports, generated summaries and duplicate master reports;
- obsolete checkpoints associated with failed local biological runs.
- interpreter/test caches and incomplete launch artifacts that contain no
  feature, metric, model or unique provenance result.

Raw datasets, source releases, governed splits, P1B PASS checkpoints, D0/D1 PASS
evidence and the frozen theory are not deleted.

## 6. Failure-Triage Decision Rules

1. Do not restore a deleted experimental implementation merely to reproduce a
   failed metric; use the numerical record in `history.md`.
2. Code enters `model/` or normal `scripts/` only after its registered Gate passes.
3. A synthetic-only PASS cannot authorize real affinity or production integration.
4. A correct-vs-deranged gain is insufficient unless correct also improves over
   ligand-only under the registered estimand.
5. Mathematical interface compatibility is not evidence of biological depth.
6. H0C localized task-local gain to ligand/series adaptation, and X0 then found
   the current ChEMBL crossed estimand underdetermined. No ChEMBL X1/X2 label
   look is authorized. Continuing requires a newly registered label-blind census
   of a more genuinely crossed source corpus; H0-B/RFSA and angular expansion
   remain frozen.

7. A rank-based negative is evidence about ranking only. A location-based
   negative would be evidence about location only. Neither generalizes to the
   other estimand, and neither is evidence about interaction.
8. X0-B's conditional design support is not permission to read interaction
   values. X1 requires a separate registration whose first step estimates `rho`
   and abstains if its upper bound exceeds `rho*`.
9. A component-level result on consumed PKIS/Anastassiadis panels does not
   override a registered total `NOT_ADMISSIBLE` verdict and does not identify
   Ki/Kd affinity energetics.
10. The research seven-point categorical-law bridge shares the abstract
    `z -> F -> B F -> K` notation but is not the frozen 33-point CDF-band
    operator. Interface resemblance is not type equivalence or production
    admission.

## 7. Current Stop Boundary

```text
AUTHORIZED: repository consolidation and E0R1 evidence preservation
COMPLETE READOUT AUDIT: E-AFF-R0 rank metric blind to task-level location
COMPLETE SPECIFICATION AUDIT: E-AFF-X0-FEAS requirement above its unit ceiling
COMPLETE DESIGN RE-REGISTRATION: E-AFF-X0-B conditional on rho
NOT RUN: E-AFF-L0 Claim A location Gate, numerical precondition failed
NOT RUN: E-AFF-L0R corrected Gate, positive control absent
COMPLETE IN RESEARCH: E0 synthetic objective/design/solver closure
COMPLETE NEGATIVE PILOT: T-DIR-P0 annotation/learnability feasibility
COMPLETE STRUCTURAL PASS: T-BASIS-R0 fixed two-body radial recoverability
COMPLETE SOURCE NEGATIVE: E-AFF-P0 shared direction not observed
COMPLETE MIXED DIAGNOSTIC: E-AFF-H0A headroom without partner specificity
COMPLETE SOURCE NEGATIVE: E-AFF-H0C centered interaction residual not observed
COMPLETE DATA STOP: E-AFF-X0 ChEMBL crossed source underdetermined
COMPLETE RESEARCH INTAKE: PKIS/F6I component-level development evidence
REGISTERED TOTAL F6I VERDICT: NOT ADMISSIBLE
RESEARCH LAW BRIDGE: INVARIANTS PASS, PRODUCTION EQUIVALENCE NOT ESTABLISHED
NOT RUN: fresh endpoint-consistent external admission
NOT AUTHORIZED: E-AFF-X1/X2 value access
NOT AUTHORIZED/FROZEN: H0-B, angular basis, replication, RFSA, DAVIS, production, P2-P4
NEXT AUTOMATIC STAGE: none
```

Two registered continuations now exist, addressing different estimands.

For **Claim A**, `E-AFF-L0` is registered and not run: a five-arm
location Gate under the frozen operator, with held-out proteins, a mandatory
protein-sequence-only control, no shared support nuisance, an interval-score
primary metric and a `0.5 * sigma_assay` margin required simultaneously against
ligand-only, sequence-only and deranged.

For **Claim B**, `X1` requires a separate registration whose first step
estimates intra-cluster correlation and abstains if its upper confidence bound
exceeds X0-B's `rho*`. Acquiring a new crossed corpus is no longer the named
continuation, because X0-FEAS showed the previous unit made that action
self-defeating.

No existing research file is an authorization to proceed automatically, and
neither continuation authorizes the other.
