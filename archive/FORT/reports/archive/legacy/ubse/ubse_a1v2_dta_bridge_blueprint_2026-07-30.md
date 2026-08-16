# UBSE-A1-v2 finite DTA bridge blueprint

Date: 2026-07-30  
Status: binding design and stage-budget record, not an event-extraction,
student-training, affinity-loading, confirmation-scoring, or sealed-test
authorization  
Current decision: `NO_DECISION_A1V2_STAGE1_EXTRACTOR_TOPOLOGY_PENDING`

## 1. Final task and role of A1

The final task remains zero-support continuous affinity prediction for a
jointly unseen target and unseen ligand scaffold under target-homology,
chemical-neighbour, assay, document, provenance, confirmation, and sealed
test closure.

A1-v2 is not the final model and cannot succeed by establishing an
interesting structural label alone. It is a candidate training-time
interaction supervision source. Its contribution exists only if the full
chain is demonstrated:

```text
reliable typed residue-FG events
  -> non-marginal pair-specific coupling
  -> deployment-side student predicts that coupling
  -> the frozen coupling coordinate is used by an affinity model
  -> strict dual-cold DTA improves beyond matched baselines and its MDE
```

## 2. Current stage in the required ten-field format

1. Current stage: Stage 1, construction and topology.
2. Unique hypothesis: the frozen A1 roles contain enough correctly resolved,
   independently auditable residue-by-functional-group-by-event topology to
   support a powered coupling test.
3. Allowed data: frozen safe metadata, locators, package/version metadata,
   and later only coordinate bodies admitted by a separate hash/download
   contract. Forbidden data: every affinity value, A1-C event/contact value,
   confirmation score, sealed outcome, and any un-audited model output.
4. Strongest matched null: a pair-conditioned rank-one event rate with the
   same deployment inputs and comparable parameter budget.
5. Destructive controls: wrong target, matched wrong ligand, FG
   derangement, residue/position destruction, event-channel permutation,
   fixed-margin transport, and fixed-mass/dustbin controls.
6. Independent units: A1-R has 153 target units but only 124 frozen
   PDB/PubMed/physical-ligand components; worst-case Kish effective size is
   98.7722. A1-C has 512 metadata-isolated primary candidates, but its
   realized legal checkerboard count is still unknown.
7. Numeric result so far: 1,035/1,035 locators and 997/997 strict HEAD
   availability pass; no event value has been read.
8. Decision: `NO_DECISION_A1V2_STAGE1_EXTRACTOR_TOPOLOGY_PENDING`.
9. Direct DTA effect: none yet.
10. Distance to trainable affinity model: the DTA bridge is now specified,
    but Stages 1-3 must pass before affinity can be loaded.

## 3. Literature and novelty boundary

DOI metadata was independently checked through Crossref on 2026-07-30.

| Work | Original task | Transferable mechanism | Boundary here |
|---|---|---|---|
| PLIP, DOI `10.1093/nar/gkv315` | deterministic complex interaction profiling | typed event extraction | one extractor is not independent biological replication |
| ProLIF, DOI `10.1186/s13321-021-00548-6` | interaction fingerprints | independent implementation vocabulary | agreement may still reflect shared geometric definitions |
| Arpeggio, DOI `10.1016/j.jmb.2016.12.004` | atom-level interaction calculation | possible second extractor/manual audit aid | event taxonomies require a frozen crosswalk |
| LINKER, DOI `10.1021/acs.jcim.6c00527` | sequence/SMILES prediction of residue-FG interactions | closest student and required strong baseline | typed maps, PLIP labels, and two-stage DTA are already covered |
| PSICHIC, DOI `10.1038/s42256-024-00847-1` | sequence/ligand physicochemical interaction fingerprints | deployable pair representation | public training membership is not adequate for the present confirmation firewall |
| ATOMICA, DOI `10.1101/2025.04.02.646906` | self-supervised interface representations from complexes | interface-as-object pretraining | preprint; structure-dependent and not strict sequence/2D dual-cold affinity evidence |
| interaction-aware LLM perspective, DOI `10.1038/s42004-025-01883-7` | compositional biological interaction models | motivates pair rather than entity representation | a perspective does not supply an identifiable observation |
| CleanSplit, DOI `10.1038/s42256-025-01124-5` | data-bias-resistant affinity evaluation | ligand/pocket/template-aware closure | split hygiene is necessary but not a model innovation |
| DynamicBind, DOI `10.1038/s41467-024-45461-2` | ligand-specific complex prediction | deployment-side structural upper bound | training/template membership must be audited; pose prediction is not DTA gain |
| NeuralPLexer, DOI `10.1038/s42256-024-00792-z` | state-specific complex prediction | cofold upper bound | unavailable or contaminated structure cannot enter the primary arm |
| ConPLex, DOI `10.1073/pnas.2220778120` | contrastive DTI from protein-language space | matched direct pair baseline | pair contrast can still learn target/scaffold shortcuts |
| AdaMBind, DOI `10.1038/s41467-026-70554-5` | target-wise 5/40-shot adaptation | later secondary diagnostic | it changes the estimand and is forbidden for primary `k=0` |
| UCE, DOI `10.1038/s41586-026-10689-z` | universal cell embeddings | tokenized object/atlas organization | masking and an atlas add no observation channel |

The earlier A1 correction cited `10.1016/j.cels.2019.08.002` as MONN, but
that DOI is a microbial ecological-network paper. Crossref identifies MONN
as `10.1016/j.cels.2020.03.002`. The cited PLANET DOI
`10.1021/acs.jcim.1c01475` returns no Crossref record and is not used here.

The surviving novelty claim is deliberately narrow:

> A source- and membership-closed typed-event teacher, certified across
> depositions and an independent extractor/manual chemistry protocol, whose
> deployment student is constrained to predict fixed-margin residue-FG
> coupling beyond its own pair-conditioned rank-one null, and whose final
> affinity contribution is an exact-null linear residual with no
> target-dependent bypass.

Neither typed interaction labels, PLIP distillation, OT, predicted monomers,
masked-event learning, a two-stage pipeline, nor a universal interaction
space is novel by itself.

## 4. Finite four-stage program

### Stage 1: construction and topology

Purpose: prove that events can be defined, independently audited, split, and
powered. Freeze before coordinate-body access:

- PLIP package/commit and dependencies;
- ProLIF/Arpeggio or a manual chemistry sample and adjudication protocol;
- seven event channels and exact cross-extractor mapping;
- ligand functional-group SMARTS, overlap priority, graph symmetry classes,
  charge/protomer rules, and impossible-cell masks;
- model, assembly, altloc, occupancy, chain, water, metal, covalent-ligand,
  missing-residue, and receptor-construct rules;
- component-bootstrap seed 1729, 2,000 replicates, MDE, and event-specific
  power calculation;
- coordinate URL/body hashes and write-once extraction ledgers.

Hard stop: locator/extractor disagreement, zero legal checkerboard rank,
insufficient independent units, or fewer than the preregistered powered
event types. No student is trained.

### Stage 2: teacher coupling

Purpose: prove that the observed object contains pair-specific coupling.
For target `t`, ligand `l`, residue `i`, functional group `g`, and event
type `k`, let

```text
Y[t,l,i,g,k] in {0,1}
M[t,l,i,g,k] in {0,1}
```

be the observed event and legal-cell mask. Compare the full tensor to the
pair-conditioned separable rate

```text
lambda_rank1[t,l,i,g,k] =
    c[t,l,k] * r[t,l,i,k] * u[t,l,g,k].
```

The binding certificate is a legal within-complex 2-by-2 cycle:

```text
Omega =
  log lambda[i,g] + log lambda[j,h]
  - log lambda[i,h] - log lambda[j,g].
```

Every rank-one row/column model has `Omega = 0`. With structural masks, the
identifiable dimension is the rank of the legal checkerboard-cycle design,
not `(R-1)(G-1)` by assertion.

Required controls are the rank-one null, correct pair, within-family wrong
target, size/FG/scaffold-matched wrong ligand, FG derangement, residue and
position destruction, event permutation, fixed-margin balanced transport,
fixed-mass partial transport, dustbin-only, LINKER-form, and additive
matched-capacity models.

Hard stop: the correct coupling does not beat every matched marginal and
destruction control with a positive component-bootstrap lower bound. No
affinity is loaded.

### Stage 3: deployable student

Purpose: predict the certified coupling from inputs available at strict
dual-cold inference.

Primary inputs:

```text
protein sequence tokens       X_t: [L, 320] frozen ESM2 residue features
optional audited monomer      S_t: [L, d_s] local invariant geometry
ligand atom graph             X_l: [N, d_a]
functional-group assignment   A_l: [N, G]
residue hidden state          H_t: [L, 128]
FG hidden state               F_l: [G, 128]
event types                   K = 7
```

The primary arm is no-P0A. A future clean P0A-v2 may supply a proposal score
`p_t: [L,1]` only as an auxiliary arm after its separate NO-GO findings are
closed. It may not be the only route to a residue.

The student first learns a pair-conditioned rank-one null:

```text
c[t,l,k]       scalar nonnegative burden
r[t,l,:,k]     softmax over residues
u[t,l,:,k]     softmax over functional groups
Pi0[k]         = MaskedIndependenceProject(r[k], u[k], M[k])
```

Both null and full model see the same deployment inputs. The null is trained
and then frozen. The full model adds only a rank-32 residue-FG compatibility
kernel. A balanced Sinkhorn/IPF projection uses the frozen `r` and `u` as
exact real-real marginals:

```text
score[i,g,k] = <U_k H_t[i], V_k F_l[g]> / sqrt(32)
Pi1[k] = FixedMarginProject(exp(score[k]), r[k], u[k], M[k]).
```

On full support, `Pi0 = r outer u`; on masked support it is the unique frozen
independence projection on the same legal cells used by `Pi1`. Zero
compatibility must recover `Pi0` cell by cell. If a dustbin is needed,
its masses, scores, masks, and real-real normalization are fixed before
events and the same exact recovery test applies.

The fixed-size coupling-only representation is

```text
z_int[k] =
  sum_(i,g) (Pi1[i,g,k] - Pi0[i,g,k])
            * (W_h H_t[i] elementwise_mul W_f F_l[g])

z_int in R^(7*32) = R^224.
```

It is exactly zero at the rank-one null and contains no explicit row/column
marginal coordinate. This does not mathematically exclude entity identity
from learned compatibility features; target-only, ligand-only, family, and
scaffold destruction remain binding empirical gates.

Training losses on A1-S fit are:

```text
L_cell   = masked event likelihood for c * Pi1
L_cycle  = logistic loss on preregistered legal checkerboard directions
L_int    = L_cell + lambda_cycle * L_cycle
```

Matched wrong/destruction arms are evaluation controls, not arbitrary
negative examples that could teach dataset identity. Frozen A1-S development
selects the one preregistered setting. A1-C is read once only after the
student, thresholds, and controls are frozen.

Gradient path:

- frozen PLM and optional monomer source initially receive no gradients;
- the rank-one null is trained first and frozen;
- coupling training updates only ligand/residue adapters and compatibility;
- P0A receives no A1 gradient;
- no affinity gradient exists in Stage 3.

Hard stop: full-minus-rank-one coupling, cycle direction, and semantic
destruction gates fail, or target-only/ligand-only/position-free/matched
direct baselines match the student. No affinity is loaded.

### Stage 4: minimal DTA prototype

Purpose: test whether the certified deployment coordinate improves actual
strict dual-cold affinity prediction.

Freeze the Stage-3 student. Use the existing component-cross-fitted
ligand-only B0 and separate pKi/pKd heads:

```text
y_hat[t,l,e] = B0[l,e] + theta[e]^T z_int[t,l]
theta[e] in R^224.
```

The exact nested null is `theta = 0`, which recovers B0. The affinity head
has no raw target, protein embedding, pair-burden, row/column marginal, source
identifier, assay identifier beyond the allowed endpoint head, or alternate
target-dependent path. This prevents the interaction representation from
being an ignorable prompt.

Minimal Stage-4 loss is the same endpoint-specific affinity loss used by
the matched B0 plus regularization on `theta`. There is no artificial
`L_use` reward: architectural exclusion of the bypass and destructive
evaluation establish use. The structure loss does not receive affinity
gradients in the first prototype.

Required equal-budget arms:

- strongest reproducible B0;
- target-only and ligand-only;
- additive two-tower;
- pair-conditioned rank-one/pair-burden only;
- LINKER-form and monomer-augmented LINKER;
- matched-capacity direct interaction model without A1 supervision;
- frozen A1 coupling model;
- wrong-target, wrong-ligand, FG-shuffled, event-shuffled, and
  position-shuffled A1 coordinates;
- zeroed interaction coordinate and deleted-branch ablation.

Primary evaluation uses the existing strict dual-cold development protocol,
component-macro target Spearman, RMSE, and a ranking/calibration auxiliary.
At preregistration, the required gain is at least the larger of the
then-current empirical MDE and the existing 0.0586 macro-Spearman MDE, with a
positive component-bootstrap 95% lower bound, all seeds positive, no family,
scaffold, source, or assay dominance, and RMSE no worse than 1.02 times the
strongest matched baseline. Destruction must remove the preregistered
fraction of the gain.

Only a Stage-4 development PASS can request confirmation scoring. Sealed test
access remains a later, separate authorization.

## 5. Leakage firewall

- A1-R events never enter student fitting.
- A1-S fit/development and A1-C must remain closed on target/homology, PDB,
  PubMed, connectivity, scaffold, ECFP4, template, pocket, and inherited
  model membership. Current P0A fails the inherited-membership condition and
  is excluded from the primary arm.
- A1-C event values are unavailable for architecture, threshold, epoch,
  event-type, extractor, FG, or control selection.
- ChEMBL affinity is unavailable until Stages 1-3 pass and a Stage-4
  preregistration is frozen.
- Development, confirmation, and sealed affinity roles never construct the
  event representation or B0 cross-fits.
- No test-time affinity adaptation, transductive normalization, query-batch
  fitting, or structure retrieval using protected entities is allowed.

## 6. Compute budgets and checkpoints

- Stage 1: CPU parsing, graph construction, extractor execution, and
  component/power audit; no GPU justification.
- Stage 2: CPU sparse contingency/cycle analysis; GPU only if a registered
  tensor calculation is measurably faster.
- Stage 3 smoke: one maximum-shape AMP forward/backward plus synthetic
  exact-null tests, then a one-seed bounded pilot. Full three-seed training
  is forbidden until independent scientific and implementation reviews pass.
- Stage 4: smallest linear frozen-student prototype first; no joint
  fine-tuning or larger backbone until the frozen prototype passes.

Every long run requires atomic per-seed/epoch checkpoints, verified resume,
persistent stdout/stderr and heartbeat, utilization/memory/power/temperature
telemetry, a wall-time budget, and fail-fast finite/progress gates.

## 7. P0A and backup disposition

P0A-v2 currently has a separate prelaunch NO-GO. Because the primary student
is no-P0A, that 5.5-hour retraining is deferred and does not block Stages 1-3.
It may return only after its implementation findings are fixed and A1-direct
establishes a powered coupling object that justifies the proposal arm.

No new public-data backup is active. Only if A1 stops for unreliability,
zero coupling rank, insufficient topology, deployment-student failure,
deployment unavailability, or lack of novelty may one backup be proposed.
The existing prospective, provenance-separated cycle-closing acquisition
design is the current fallback information source; it is not authorized by
this blueprint.

## 8. Binding next action

The only next action is to preregister the Stage-1 extractor, FG, assembly,
component-inference, power, and coordinate-body contracts. Coordinate bodies
remain locked until that document and its tests pass. P0A-v2 CUDA training,
A1 event extraction, student training, affinity loading, confirmation, and
sealed access remain unauthorized.
