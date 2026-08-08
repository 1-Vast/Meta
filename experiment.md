# MetaSieve Experiment Registry

Updated: 2026-08-08

## Terminal Evidence

| Stage | Status | Scientific meaning |
|---|---|---|
| P0 | PASS | canonical data and sealing contracts work |
| P1A | PASS | open holo structure corpus is governable |
| P1B | PASS | geometry depends strongly on the correct protein |
| P1C | FAIL | frozen geometry readout adds no affinity direction |
| P1R1 | FAIL | invariant chemistry MIF is protein-sensitive but subthreshold |
| P1R2A | FAIL | non-additive MIF does not improve ligand prior |
| P1R2B0 | CLOSED | nonlinear global readouts do not rescue MIF |
| P1R2B1 | CLOSED | compatibility signal is mainly wrong-protein penalty |
| F0R | HISTORICAL FAILURE | live API is not immutable provenance |
| D0-C | PASS | ChEMBL37 static release corpus is reproducible |
| D1 | PASS | homology/document closure has zero cross-fold leakage |
| E0 | FAIL-CLOSED | synthetic partner-identification Gate failed |
| E0R0 | MIXED | analytic tensor passes; learned heads fail |
| E0R1 | PROCEDURAL NOT-RUN | pinv passes; deterministic solve misses RMSE precondition |
| E0R2 | SYNTHETIC PASS | deterministic augmented solve closes corrected objective/design/solver |
| T-DIR-P0 | PILOT NEGATIVE | PLIP direct-event pipeline works, but frozen local-state D2 does not replicate held-out learnability |
| T-BASIS-R0 | STRUCTURAL PASS | fixed two-body radial basis is recoverable and correct-partner dependent |
| E-AFF-P0 | SOURCE NEGATIVE | no population-shared 288D radial residual-affinity direction was observed |
| E-AFF-H0A | MIXED/STOP | task-local held-out affinity headroom exists, but correct-partner specificity is below the frozen margin |
| E-AFF-H0C | SOURCE NEGATIVE | support-matched ligand nuisance succeeds, but centered radial interaction adds no affinity or partner value |
| E-AFF-X0 | DATA STOP | nominal ChEMBL rectangles collapse to too few independent panel/closure components; X1 is not authorized |
| E-AFF-R0 | READOUT SCOPE AUDIT | within-task concordance is exactly invariant to task-level affinity location, so the historical affinity negatives are scoped to ranking information |
| E-AFF-X0-FEAS | SPECIFICATION AUDIT | the frozen 245-component requirement exceeds the structural ceiling of X0's own unit at every governed population |
| E-AFF-X0-B | CONDITIONAL DESIGN | cell-disjoint units reach the unchanged 245 requirement only if intra-cluster correlation is small; not evidence of interaction |
| E-AFF-L0 | NOT RUN | Claim A Gate failed a registered precondition; a gate condition used a degenerate coverage statistic, so Claim A is untested |
| E-AFF-L0R | NOT RUN | corrected Gate; both defects fixed, but the registered ligand-only positive control failed, so no protein contrast is interpretable |
| PKIS mechanism v2 | MIXED DEVELOPMENT | product-kernel interaction transfers on PKIS2 but does not establish a fresh endpoint-consistent affinity result |
| Section F6I | COMPONENT EVIDENCE / TOTAL NOT ADMISSIBLE | protein-dependent zero-shot surface and one-dimensional support location separate the named gauge; registered total Gate remains not admissible |
| Research law bridge | INVARIANTS PASS / RESEARCH ONLY | seven-point categorical-law bridge passes local invariants but is not type-equivalent to the frozen production CDF-band operator |
| XP1-A | EXISTENCE IDENTIFIED | on a complete crossed panel the protein-by-ligand interaction is 59.6% of affinity variance, 40.3% of it reproducible, rank-saturating by 8; its protein-side geometry replicates across disjoint compound halves (`r=0.885`) and across an independent platform (`r=0.565`) |
| XP1-B | MIXED / DECISIVE | under strict kinase-group closure a rank-1..3 section identified from `k` support labels reaches `R2_gamma=+0.160 [+0.109,+0.195]` and is target-specific (`+0.0695 [+0.0524,+0.0874]` vs deranged support); every zero-shot protein representation, ESM-2 included, reaches `~0` |
| XP1-C | DIRECTIONAL REPLICATION | PDSP `Ki` (GPCR/transporter, radioligand assay, MMseqs2 40% closure) reproduces the pattern at `A4 R2_gamma=+0.024 [+0.007,+0.052]`, below the preregistered `0.05` negligibility floor |
| XP1-D | REPRESENTATION LOCALISED | the interaction coordinate is predictable from the aligned KLIFS pocket at `R2=0.52` leave-one-protein-out and at `-0.087` under group closure, indistinguishable from a random-feature null; learning curve flat over n=20..65 |
| XP1-E | DESTRUCTIVE CONTROL PASSED | a synthetic additive panel with the same margins, noise and 40% truncation yields `Delta_specific=+0.00004 [-0.0005,+0.0011]`; left-censoring cannot manufacture the XP1-B result |
| XP2-A | REPRODUCED | `XP1_EVIDENCE_REPRODUCED`, 18/18 checks; all five release checksums match; XP1's measured cells match the journal supplement at max `|diff|=0.0`; XP1's single-floor censoring *description* corrected (50 distinct thresholds plus untested blanks), analysis set unchanged |
| XP2-B | PASS | ligand loading transfers to unseen Bemis-Murcko scaffolds: gauge-invariant reconstruction `R2=+0.199 [+0.133,+0.261]` (ECFP, `d=3`) vs random-feature `+0.025` and mean-loading `+0.024`; `u(L)` is not a lookup table |
| XP2-C | FAIL (below floor) | identifiable section dimension is exactly `min(k-1,d)`, measured `0,1,2,3,3` at `k=1..5`; at `k=1` ridge returns `v=0` and the arm is identically additive; `R2_gamma` peaks at `+0.0248 [+0.0114,+0.0321]` against the frozen `0.05` floor; `Delta_deploy` over ligand-only never clears zero |
| XP2-D | FAIL (specificity lost) | under simultaneous protein-group and ligand-scaffold closure at the primary `k=5, d=3`: `R2_gamma=+0.0199 [+0.0076,+0.0283]`, derangement specificity `+0.00185 [-0.00477,+0.00552]` (CI spans zero), `Delta_deploy=-0.0331`; design rank `3.00/3`, coverage `1.000`, so not an identifiability artefact. XP1's specificity was conditional on ligand reuse |
| XP3 | CENSUS | public data offers either low noise with ~8 protein-group units (single-lab profiling) or 70-85 independent units with noise exceeding the interaction (literature aggregation); no accessible source offers both |
| XP4 | DATA-LIMITED NULL | 85 independent BindingDB panels, 70 protein clusters: per-report `sigma=0.777`, `gamma` sd `0.406` below the `0.650` noise alone transmits; within-panel chemistry-neighbour ceiling `R2=-0.539` |
| XP5 | FAIL | ten fixed named physicochemical complementarity channels, computed not predicted: `R2_gamma=-0.0015`, derangement specificity exactly `+0.00000 [-0.00069,+0.00066]`; shuffling the pairing improves the fit |
| S1 | INDEPENDENT SET FORMED | RCSB CC0, released >=2024, 15,003 candidates with **exposed overlap 0** against all 10,468 P1B-exposed ids; 1,118 usable complexes, 621 protein clusters, 586 scaffolds |
| S2 | TEACHER PASS | six named 3D channels from raw holo coordinates; rot/trans `7.6e-15`, permutation `1.6e-14`, determinism `0.0`, no channel degenerate |
| S3 | WELL POWERED | 124 effective units per fold, minimum detectable `R2=0.02` at 100% detection; S6 floor frozen at `0.02` before opening the test set |
| S4 | MIXED / DECISIVE | H-bond `R2=+0.268 [+0.166,+0.378]` and hydrophobic burial `+0.299 [+0.162,+0.454]` beat random features, but **no channel beats deranged protein** |
| S4b | LIGAND-SIDE | protein-only `R2 ~ 0` for all six; ligand-only matches the joint model; adding the protein costs `-0.032` on hydrophobic burial with CI excluding zero |
| S5 | NOT AUTHORISED | the distillation target is a ligand descriptor; training it would learn the forbidden population shortcut |
| XP2-G | CONDITIONAL | interface legality requires a declared gauge, a two-term outer radius, and placement of the discrete coordinates in `kappa`; abstention is the existing `p=e_0` vertex; CSMO/Band/`K`/mesh untouched |

All detailed metrics and former report conclusions are preserved in `history.md`
and `EVIDENCE_CONSOLIDATION_AND_FAILURE_TRIAGE.md`.

The PKIS/F6I preregistrations, code and compact artifacts were audited, committed
in `8b7789e`, consolidated into `history.md`, and then removed from the active
working tree. PKIS2 and Anastassiadis are consumed development panels. The F6I
registered total verdict remains `F6I_COMPONENTS_NOT_ADMISSIBLE`; this update is
not a production PASS.

## Estimand Separation

`E-AFF-L0` addresses **Claim A**: does the correct protein provide
affinity-location information beyond population, ligand-only and
protein-sequence-only baselines? `E-AFF-X0-B` and `X1` address **Claim B**: does
a non-additive protein-by-ligand affinity interaction exist?

These are different estimands. Neither stage replaces, implies nor authorizes
the other.

The P1C/P1R\*/E-AFF-P0/H0A/H0C evaluations used within-task concordance or
another rank-based metric. Those experiments remain valid; their negative
conclusions are limited to within-task ranking information, because within-task
concordance is invariant to task-wise prediction shifts and positive rescaling.
No historical verdict is reclassified.

Repository retention is evidence-based rather than outcome-based: terminal
negative stages retain their preregistration, result, independent audit and
final artifact, while caches and launches that produced no feature/metric/model
are not experiments and are removed after being logged in `history.md`.

## Current Registered Research

E0R2, T-BASIS-R0, E-AFF-P0 and E-AFF-H0A are complete. The retained research
package records the synthetic E0 evidence, directional-potential claim audit,
negative sparse-event pilot, fixed-radial structural PASS, and governed
ChEMBL37 source-affinity diagnostics. Only the registered E-AFF stages read
source Ki/Kd labels; DAVIS and recipient labels remain prohibited.

E0R1 key evidence:

- ordering conflict: `370/6080 = 6.0855%`;
- design rank: `225/240`, mean holdout coverage `0.999816`;
- pinv train RMSE: `3.18e-8`;
- pinv correct/deranged CI: `0.99737/0.61447`;
- deterministic solver train RMSE: `7.92e-4 > 1e-6`;
- formal verdict: `NOT_RUN_NUMERICAL_PRECONDITION_FAILED`.

E0R2 closure:

- train RMSE: `3.188e-8`;
- corrected objective: `1.567e-15`;
- full-gradient L2: `6.150e-17`;
- historical correct/deranged CI: `0.99737/0.61447`;
- verdict: `SYNTHETIC_OBJECTIVE_DESIGN_SOLVER_IDENTIFIED`.

T-DIR-P0 structure-only pilot:

- selected complexes: `24/8/8`, 40 distinct homology groups;
- annotation completion: `40/40`; oracle-near pairs: `17,556`;
- hydrophobic D0/D1/D2 test AP: `0.00987/0.03735/0.01996`;
- hydrophobic test prevalence: `0.00620`;
- D2 train/validation/test AP: `0.25496/0.01699/0.01996`;
- verdict: `PILOT_LEARNABILITY_SIGNAL_NOT_OBSERVED`;
- post-run limitation: group-channel mapping defect and omitted shuffle/
  chemistry-roundtrip controls; no same-panel retry.

T-BASIS-R0 fixed radial basis:

- fresh split: `192/64/64`, 320 distinct homology groups;
- T-DIR-P0 overlap: `0`; derangement max identity: `0.32632`; reuse: `0`;
- fixed basis: `8 x 6 x 6 = 288D`;
- test mean/correct/deranged MSE: `1.1001/0.5157/0.6875`;
- test reconstruction gain: `0.5312 [0.4433,0.5962]`;
- test partner gain: `0.1561 [0.1070,0.2007]`;
- verdict: `RADIAL_BASIS_PARTNER_RECOVERABILITY_IDENTIFIED`;
- scope: structural two-body radial moments only; no affinity reads.

E-AFF-P0 population-shared source feasibility:

- sample: one score-blind task from each of 245 closure components, 4,900 rows;
- component-macro ligand/correct/deranged/null CI:
  `0.55225/0.54209/0.54445/0.54210`;
- correct-minus-ligand: `-0.01016 [-0.02069,-0.00018]`;
- correct-minus-deranged: `-0.00236 [-0.01301,0.00805]`;
- correct-minus-coupling-null: `-0.00001 [-0.00558,0.00530]`;
- verdict: `SHARED_DIRECTION_NOT_OBSERVED_H0_DATA_SUPPORTED`.

E-AFF-H0A task-local headroom diagnostic:

- sample: 107 tasks from 107 closure components; 20 fit and 20 untouched test
  ligand states per task, 4,280 rows;
- component-macro ligand/correct/deranged/null CI:
  `0.55404/0.64226/0.63362/0.63817`;
- correct-minus-ligand: `+0.08821 [0.06761,0.10998]`;
- correct-minus-deranged: `+0.00864 [0.00338,0.01462]`;
- correct-minus-coupling-null: `+0.00408 [0.00061,0.00747]`;
- verdict: `TASK_LOCAL_RADIAL_HEADROOM_WITHOUT_PARTNER_SPECIFICITY`;
- decision: H0-B/RFSA not authorized because the frozen `+0.03` partner margin
  failed despite strong task-local headroom.

E-AFF-H0C fixed radial interaction residual:

- fresh tasks: 54 tasks/components, excluding all H0A tasks;
- split: 20 support/20 test per task with zero Murcko scaffold overlap;
- global/local/correct/deranged CI: `0.55487/0.59635/0.59244/0.59092`;
- correct-minus-Local-L: `-0.00391 [-0.02040,0.01191]`;
- correct-minus-deranged: `+0.00152 [-0.01024,0.01424]`;
- verdict: `FIXED_RADIAL_INTERACTION_RESIDUAL_NOT_OBSERVED`;
- decision at H0C completion: shortcut-removal rescue closed and X0 was the
  next discriminator; X0 has since completed with a data-stop verdict.

E-AFF-X0 crossed source census:

- label firewall: `152,737` governed IDs; zero affinity value fields selected;
- Ki: 597 panels, 1,059,169 nominal rectangles, 36 effective components,
  18 replicate-supported components;
- Kd: 34 panels, 232,875 nominal rectangles, 12 effective components,
  4 replicate-supported components;
- frozen lower-bound requirement: 245 independent components per endpoint;
- largest dependency rectangle fraction: Ki `0.56494`, Kd `0.76757`;
- verdict: `STOP_SOURCE_INTERACTION_UNDERDETERMINED`;
- decision: X1/X2 not authorized; no affinity double differences were read.

E-AFF-R0 readout scope audit:

- readout under test: within-task concordance, macro-averaged over closure
  components; theory-controlled functional: Hausdorff-`W1` between law classes;
- invariance deviation under per-task prediction shift, prediction rescale,
  label shift and label rescale: `0.0`, `0.0`, `0.0`, `0.0`;
- perfect task-level predictor concordance at level variance shares
  `0.059/0.200/0.500/0.800/0.941/0.985`: `0.5000` at every share;
- same predictor's location-sensitive RMSE advantage over a global mean across
  the same sweep: `1.033/1.001` to `8.832/1.017`;
- H0C published per-task contrasts: correct-minus-local `-0.00391`,
  deranged-minus-local `-0.00543`, correct-minus-deranged `+0.00152`, with the
  geometry term changing `51/54` tasks;
- verdict: `READOUT_BLIND_TO_TASK_LEVEL_AFFINITY_LOCATION|PERFECT_LEVEL_PREDICTOR_SCORES_CHANCE_AT_EVERY_VARIANCE_SHARE`;
- affinity labels read: none.

E-AFF-X0-FEAS unit feasibility audit:

- closure-component universe: `245`; frozen requirement: `245` per endpoint;
- E0-Core components with endpoint data: Ki `202`, Kd `72`;
- panel-free rectangle-capable components: Ki `57`, Kd `12`;
- best ceiling over all governed populations: Ki `97`, Kd `56`;
- panels touching more than one closure component: `0`, as predicted;
- frozen `245` independently re-derived from its stated chi-square design;
- verdict: `X0_UNIT_REQUIREMENT_UNATTAINABLE_BY_CONSTRUCTION`;
- consequence: the X0 stop is specification-induced; its verdict is retained.

E-AFF-X0-B crossed design re-registration:

- unchanged: effect size `0.5`, variance ratio `1.25`, alpha `0.05`, power
  `0.80`, required effective `n = 245`, affinity Gate margins `+0.03/+0.03`;
- changed: unit becomes the cell-disjoint rectangle; effective `n` becomes
  design-effect adjusted;
- Ki: `11,168` units, `36` clusters, `205` target pairs, `224` targets,
  `19,062` ligands, `rho* = 0.0915`, hard bound `G/245 = 0.1469`;
- Kd: `1,041` units, `12` clusters, `49` target pairs, `73` targets,
  `1,256` ligands, `rho* = 0.0164`, hard bound `G/245 = 0.0490`;
- at `rho = 1` the optimal design reduces to `36` and `12`, exactly X0's counts;
- verdict: `X0B_CONDITIONAL_DESIGN_SUPPORTED_KI|X0B_CONDITIONAL_DESIGN_SUPPORTED_KD`;
- scope: conditional design support only; X1 remains unauthorized until `rho` is
  estimated with uncertainty and gated against `rho*`.

## Frozen Boundaries

X1/X2, H0-B, angular/many-body extension, reference-state potentials, RFSA and
support adaptation, DAVIS and recipient labels, production integration,
CSMO/Band/theory changes, and P2-P4 remain frozen or not authorized. Recipient
and DAVIS label reads remain zero.

An `E-AFF-L0` pass would establish at most
`PROTEIN_SPECIFIC_AFFINITY_LOCATION_IDENTIFIED_IN_SOURCE`. It would not
authorize biological `z` admission, `model/` promotion, DAVIS or recipient
evaluation, X1, angular or many-body development, RFSA, theory modification,
P2-P4, an end-to-end DTA claim, or a physical binding-free-energy
interpretation. Production or biological-`z` admission still requires a
separately registered independent-source replication and a sealed novel-target
transfer Gate.
