# Preregistration — S7/L2B Phase 2A

## Audit-only attribution of B5: marginal, coupling, and label-side ligand conditionality

Stage identifier: `S7_L2B_PHASE2A_MARGINAL_COUPLING_AND_LIGAND_CONDITIONALITY_AUDIT`

Written: 2026-08-10.
Repository commit at registration: `623602e76b7d4f445af069014782278163183d59`.
Governing prior registration: `research/s7_l2b_r0r/PREREG_S7_L2B_UNIFIED.md`
(SHA-256 `2c333f223ae450c566cc62b1a3b276ff59c065c38348005ad9504ac1930b9a92`,
commit `ce186f4`).

This document is frozen before any Phase 2A metric is computed. Every threshold,
estimand, inference unit, seed, tolerance and verdict rule below is fixed here.
No quantity defined here may be altered after any result is observed.

### Chronology limitation, stated in advance

Commit authorization was not granted for this run. This registration is
therefore anchored by its SHA-256 and by the fact that every Phase 2A output
artifact records that hash, but it is **not** anchored by a git commit
timestamp. This is a strictly weaker chronological guarantee than the one used
for `ce186f4`/`139effd`, and it is recorded here rather than claimed away. The
Phase 2A verdict must be read with that limitation attached until the file is
committed.

---

## 1. Scientific question

B5 (frozen ESM2-650M residue features) passed all six registered structural
development Gates, but 92.5% of its residue AP survived substitution of an
arbitrary foreign ligand. The pair score may therefore be fully explained by

```text
generic protein-pocket residue propensity  +  ligand-atom propensity
```

with no ligand-conditioned residue-atom coupling. Phase 2A determines, without
training anything, **where the ligand-conditioning is absent**: in the labels,
in the observable inputs, or in the B5 realization.

Two causes must be separated before any repair is designed:

1. the MONN teacher itself may contain little within-protein across-ligand
   variation (a data/label property);
2. the teacher may contain such variation that B5 fails to recover (a model
   property).

## 2. Scope and prohibitions

Phase 2A is audit-only. It:

- trains and fine-tunes nothing;
- reads no ChEMBL, BindingDB, DAVIS, KIBA, recipient or other affinity value;
- downloads no corpus;
- adds no PLM, attention stack, geometry branch, typed-interaction branch,
  affinity head, PU loss, knowledge graph or parallel module;
- modifies nothing under `model/`, `scripts/`, `theory/`, CSMO, Band, mesh or
  the frozen operator `A(F,z) = K(B(z)F(z))`;
- lowers no existing Gate;
- selects no sample, channel, mask or threshold after examining results.

All Phase 2A code lives under `research/s7_l2b_r0r/`. All Phase 2A outputs live
under `report/s7_l2b_r0r/`. No historical immutable artifact is rewritten.

## 3. Frozen definitions

| Symbol | Definition |
|---|---|
| exact construct | `seq_key` = SHA-256 of the UniProt sequence used for residue indexing. Residue indices are comparable **only** within one `seq_key`. |
| protein | `uniprot_id`. |
| ligand identity | `graph_key` = SHA-256 of the canonical non-isomeric SMILES of the heavy-atom molecule. |
| scaffold | Bemis–Murcko scaffold SMILES from RDKit. |
| replicate pair | two records, same `seq_key`, same `graph_key`, different `pdb_id`. |
| alternative-ligand pair | two records, same `seq_key`, different `graph_key`. |
| scaffold-distinct pair | an alternative-ligand pair whose two Murcko scaffolds are both non-empty and different. |
| inference unit | protein closure component (`protein_components`, unchanged from Phase 1). Atom–residue rows are never inference units. |
| residue mask | `R(P,L) = { r : exists a with Y_ra = 1 }` for that record. |
| `G` | sealed per-pair float16 prediction matrix of an arm, complete `n_res x n_atoms`. |
| `Y` | binary teacher matrix, complete `n_res x n_atoms`. |

Scope split, fixed here:

- **label-side (teacher) analysis** uses the **full admitted corpus**, because it
  is a property of the labels and involves no model. Held-out-A-restricted
  values are additionally reported.
- **B5-side analysis** uses **held-out A only**, the split on which predictions
  were sealed.

## 4. Phase 0 — contract and artifact audit (executed before any Phase 2A metric)

### 4.1 Artifacts hashed

Raw corpus (`monn_development_edge_corpus.jsonl.gz`,
`monn_additional_pdb_edge_corpus.jsonl.gz`, `raw_corpus_summary.json`);
`I1_ATOM_QUARANTINE.json`; `r0r2_closure/homology_alignment_stats.json`;
all six `sealed_preds/heldoutA_*.f16.dat` plus `heldoutA_index.json` and
`ap_tables.json`; all four `sealed_preds_b5/heldoutA_*.f16.dat` plus
`ap_tables_b5.json` and `B5_checkpoint.pt`; the ESM2 residue cache and index;
`P0_SEALED_PREDICTION_MANIFEST.json`; `P1_B5_GATE.json`;
`P1_MARGINAL_DECOMPOSITION.json`; `PUBLICATION_TIME_CLOSURE_AUDIT.json`;
`I1_ATOM_CORRESPONDENCE_AUDIT.json`; `I2_COUPLING_IDENTIFIABILITY_AUDIT.json`;
`PREREG_S7_L2B_UNIFIED.md`; this file; the MONN `mol_dict` pickles.

### 4.2 Fail-closed checks

| id | check |
|---|---|
| C1 | every listed artifact exists; SHA-256 recorded. |
| C2 | recomputed prediction-file hashes equal those recorded in `P0_SEALED_PREDICTION_MANIFEST.json` (`arm_sha256`) and `P1_B5_GATE.json` (`per_pair_prediction_sha256`). |
| C3 | the held-out A offset table is rebuilt twice: once as `p0_seal_predictions` built it, once as `p1_run_b5` built it (with the ESM availability filter). Both must equal each other **and** the sealed `heldoutA_index.json`, key for key and offset for offset. |
| C4 | for every arm, file size in bytes equals `2 * total_cells`. |
| C5 | all ten arms are addressed by one index, hence identical rows and identical masks. |
| C6 | metadata required by the census is present on every admitted record: `seq_key`, `uniprot_id`, `pdb_id`, `ligand_ccd`, `graph_key`, `scaffold`, `cohort`, `n_res`, `n_atoms`, `edges`, `positive_typed_edges`; and publication/time fields are available in `PUBLICATION_TIME_CLOSURE_AUDIT.json`. |
| C7 | the set of label fields read is enumerated; the set of affinity sources opened is empty. |

`C3` is load-bearing: `p1_marginal_decomposition.py` indexed the B5-family
memmaps with the B4-family index. If the two offset tables ever differed, every
Phase 1 marginal number for B5/BX5/BP5 would be misaligned. This registration
requires that to be proved, not assumed.

Any failed check terminates the stage with
`PHASE2A_CONTRACT_OR_ARTIFACT_FAIL_CLOSED` and no further metric is computed.

## 5. Phase 1 — data identifiability census

Reported: number of proteins and exact constructs; distribution of
ligand/scaffold depth per construct; publication and release-year distribution;
closure-component count; largest-component fraction; independent multi-ligand
group depth per component; residue-mask availability; exact-ligand,
connectivity-key and Murcko-scaffold overlap between paired records; whether
each comparison is within the same exact construct; and label-blind power.

**Label-blind power.** The preregistered minimum meaningful effect for teacher
conditionality is `dJ_min = 0.05` Jaccard units. Power is computed from the
component structure only, over a preregistered assumed-SD grid
`sigma in {0.10, 0.15, 0.20, 0.25, 0.30}`, two-sided `alpha = 0.05`, reporting
achieved power at the observed component count and the component count required
for 80% power. No observed label variance enters the power statement.

**Sufficiency thresholds (frozen).** The corpus is data-identifiable iff

- `D1`: at least **30** distinct closure components contain at least one
  within-construct scaffold-distinct alternative-ligand pair; **and**
- `D2`: at least **100** within-construct scaffold-distinct alternative-ligand
  pairs exist in total; **and**
- `D3`: every record entering any pair carries a non-empty residue mask.

Failure of `D1`, `D2` or `D3` terminates the stage with
`PHASE2A_DATA_NOT_IDENTIFIABLE`.

## 6. Phase 2 — teacher ligand-conditionality ceiling

The controlling insight fixed here: **the noise floor is measurable from the
data**. Two crystal structures of the *same* construct with the *same* ligand
differ only by experimental and crystallographic variation. That replicate
Jaccard is the correct comparator for alternative-ligand Jaccard. An arbitrary
foreign ligand is a corruption control, not a biological negative, and is
reported separately and labelled as such.

**T1 — primary estimand.**

```text
dJ = mean Jaccard(replicate pairs) - mean Jaccard(alternative-ligand pairs)
```

Component-macro: average within closure component first, then over components.
Paired component bootstrap over components contributing both pair types
(10,000 resamples, seed `20260819`); the unpaired two-sample component bootstrap
is reported as a fallback when the paired set is smaller than 10 components.
`T1` PASSES iff `dJ >= 0.05` with one-sided 95% LCB `> 0`.

**T2 — cross-ligand residue AP retention.** For each ordered within-construct
pair, rank the residues of the target record by a source mask and score residue
AP against the target mask. Source masks: the record's own mask (degenerate,
reported for completeness), a replicate mask, an alternative-ligand mask, the
leave-one-out protein-marginal frequency mask over all other ligands of the
construct, and prevalence. The decisive contrast is
**alternative-ligand AP vs replicate AP**.

**T3** positive-residue gain and loss rates `|R_a \ R_b| / |R_a|` and
`|R_b \ R_a| / |R_a|`.

**T4** within-construct between-ligand residue variance
`V = sum_r ybar_r (1 - ybar_r) / sum_r ybar_r`, compared against the same
quantity computed on replicate-only constructs.

**T5** scaffold-distance sensitivity: Spearman `rho` between mask dissimilarity
`1 - Jaccard` and ligand chemical distance `1 - Tanimoto` on RDKit Morgan
fingerprints (radius 2, 2048 bits, frozen), computed within construct and
summarized by component-macro.

**T6 — ligand-permutation control.** Within each construct, permute the
assignment of ligand chemical descriptors to residue masks and recompute `T5`.
200 permutations, seed `20260820`. Reported as a permutation p-value and z.
`T6` PASSES iff the observed `rho` exceeds the permutation null with
`p <= 0.05` in the direction of positive association.

**T7** per-construct and component-macro summaries, and the fraction of pairs
with a *biologically meaningful* residue change, frozen as
`Jaccard <= 0.5 AND |R_a symmetric-difference R_b| >= 3 residues`.

**Teacher verdict.** The teacher is declared ligand-conditioned iff **`T1` and
`T6` both pass**. This conjunction is the conservative reading: `T1` establishes
that variation exceeds the replicate noise floor, `T6` establishes that the
variation is associated with ligand chemistry rather than unstructured. If `T1`
passes and `T6` fails, the verdict is still `TEACHER_GENERIC_POCKET_ONLY` with
qualifier `variation_present_but_not_chemistry_associated`, because a
ligand-conditioned head cannot be trained on unpredictable variation.

Weak teacher variation is a property of the labels and must never be reported as
model failure.

## 7. Phase 3 — marginal and coupling decomposition

On the actual observed pair mask with the frozen pair weights. For held-out A
the mask is the complete `n_res x n_atoms` matrix with uniform weights; this is
verified in Phase 0, and only under that verified condition is classical double
centering admissible. The general weighted-ALS solver
(`weighted_additive_fit`, already contract-tested) is used regardless, with
double centering retained solely as a numerical self-test.

**Teacher additive null (label-fitted oracle, never a deployable arm).**

```text
logit P(Y_ra = 1) = mu + alpha_r + beta_a
```

fitted per complex by Newton/IRLS with ridge `1e-6` for separation control.

**B5 additive projection (label-free, deployable decomposition).**

```text
G_add = argmin_A || W^(1/2) (G - A) ||^2,   A in span{1, alpha_r, beta_a}
C     = G - G_add
```

**Orthogonality tolerance (frozen).**

```text
|| X^T W C ||_F / (1 + || C ||_F) <= 1e-8
```

Violation is a fail-closed contract error, not a reported number.

Pair AP (tie-aware, component-macro) is reported for full `G`, `G_add` and `C`
for arms `B5`, `B4`, `BX5`, `BP5`, `BL`, `B0`. Deployable prediction
decompositions, label-fitted oracle ceilings and evaluation nulls are reported
in three separately labelled blocks and are never mixed.

## 8. Phase 4 — matched attribution controls

Under identical masks: (1) full B5; (2) residue marginal; (3) atom marginal;
(4) weighted additive prediction; (5) marginal-orthogonal coupling residual;
(6) wrong ligand; (7) wrong protein; (8) within-complex pair shuffle where
valid; (9) degree-preserving bipartite rewiring.

**Degree-preserving rewiring specification (frozen).**

- legal swap: pick positive edges `(r1,a1)`, `(r2,a2)` with `r1 != r2`,
  `a1 != a2`, `A[r1,a2] = 0`, `A[r2,a1] = 0`; set `A[r1,a1] = A[r2,a2] = 0` and
  `A[r1,a2] = A[r2,a1] = 1`. Every row and column degree is exactly preserved.
- burn-in: `100 * E` successful swaps discarded before the first sample, where
  `E` is the number of positive edges.
- swaps between successive samples: `30 * E` successful swaps.
- independent rewires: **20** per complex; seed `20260817`.
- non-switchable graphs: a complex is `non_switchable` if it has fewer than 4
  positives, fewer than 3 active residues, fewer than 3 active atoms, or fewer
  than 50% swap success within `8x` the target attempt budget. Such complexes
  are excluded from the null **with their count reported**, and every point
  estimate is additionally recomputed on the switchable subset so the exclusion
  cannot move a verdict silently.
- exact degree preservation is asserted on every sample; a violation is a
  fail-closed contract error.
- mixing diagnostics: edge-overlap fraction against the original edge set at
  `{0, 1, 5, 10, 30} x E` swaps (decay curve), and the overlap between
  successive independent samples.

Rewiring is an attribution null only. It is never a non-binder label and never a
training negative.

Inference: component-macro tie-aware AP, effect sizes, and component bootstrap
intervals (10,000 resamples, seed `20260818`). Confidence intervals are never
derived by treating pair rows as independent.

## 9. Phase 5 — label-semantics audit

Using only already available evidence:

- interaction-type census from `positive_typed_edges` and `positive_event_edges`;
- the water-mediated fraction. Water bridges are **indirect** contacts; the
  primary teacher decomposition is therefore recomputed with them removed as a
  **preregistered sensitivity**, fixed here rather than chosen later;
- coverage: positives per complex, per residue, per atom, and density;
- missing atom/residue mappings from the I-1 quarantine, hydrogen-landing edges,
  out-of-range residue indices;
- consistency with a dense distance/contact teacher **if one exists for these
  complexes**;
- disagreement with a second frozen interaction tool **if one exists**;
- dependence on interaction type and distance threshold.

If the required comparator does not exist for this corpus, the question is
recorded as `UNRESOLVED` with the reason. Missing-positive contamination is
never inferred from low AP. PU learning or a soft teacher may be authorized only
if ambiguity is positively demonstrated.

`LABEL_SEMANTICS_AMBIGUOUS` requires positive demonstration, defined here as:
at least **20%** of positive edges are indirect/water-mediated **and** the `T1`
conclusion reverses when they are removed; **or** a second frozen teacher exists
and disagrees on at least **20%** of edges.

## 10. Terminal verdict — frozen precedence

Evaluated strictly in this order; the first satisfied rule is the verdict.

```text
1. any Phase 0 check fails
     -> PHASE2A_CONTRACT_OR_ARTIFACT_FAIL_CLOSED
2. D1, D2 or D3 fails
     -> PHASE2A_DATA_NOT_IDENTIFIABLE
3. T1 fails OR T6 fails
     -> TEACHER_GENERIC_POCKET_ONLY
4. label-semantics ambiguity positively demonstrated (section 9)
     -> LABEL_SEMANTICS_AMBIGUOUS
5. BC true
     -> EDGE_COUPLING_ALREADY_IDENTIFIED
6. BC false and TC true
     -> EDGE_COUPLING_PRESENT_B5_ABSENT
7. BC false and TC false
     -> LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING
```

where

- **TC** (teacher edge coupling present): the teacher's marginal-orthogonal
  coupling statistic exceeds its degree-preserving rewiring null with median
  `z >= 2.0` **and** at least 60% of switchable complexes above their own null.
- **BC** (B5 edge coupling identified): the component-macro tie-aware pair AP of
  B5's coupling residual `C` exceeds **both** (a) its degree-preserving rewiring
  null and (b) the wrong-ligand arm `BX5` coupling residual, each by
  `>= 0.01` with one-sided 95% LCB `> 0`.

These verdicts are mutually exclusive and exhaustive by construction. The
verdict is chosen by the rules above, never by which reading is most
interesting.

## 11. Mandatory verdict-to-action mapping

| verdict | only permitted next action |
|---|---|
| `PHASE2A_CONTRACT_OR_ARTIFACT_FAIL_CLOSED` | repair only the missing evidence contract; run nothing else. |
| `PHASE2A_DATA_NOT_IDENTIFIABLE` | preregister a metadata-only census for a new multi-ligand structural corpus; do not train. |
| `TEACHER_GENERIC_POCKET_ONLY` | close the exact-coupling claim on the current MONN teacher; preregister a dense same-protein multi-ligand structural corpus; do not repair B5 on the same labels. |
| `LABEL_SEMANTICS_AMBIGUOUS` | preregister one dense continuous-coordinate or audited soft-teacher reconstruction; do not change the learner. |
| `LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING` | preregister one ligand-conditioned residue residual head. |
| `EDGE_COUPLING_PRESENT_B5_ABSENT` | preregister the T0–T3 input-observability ladder first; do not train a pair head. |
| `EDGE_COUPLING_ALREADY_IDENTIFIED` | do not repair; preregister one sealed independent structural confirmation. |

No Phase 2A verdict authorizes affinity labels, DAVIS, KIBA, recipient data,
few-shot adaptation, admission of any statistic into `z`, or any change to
CSMO, Band, mesh, or `A(F,z) = K(B(z)F(z))`.

## 12. Registered expectations

Registered in advance so that a surprising result cannot be re-narrated after
the fact:

1. The teacher is expected to show **some** ligand conditionality, because
   different ligands occupy different sub-pockets, but the margin over the
   replicate noise floor is expected to be **modest**.
2. B5's coupling residual is expected to be **weak**, because I-2 measured
   coupling beyond degree-preserving margins at median `z = +0.41` and 92.5% of
   residue AP survived foreign-ligand substitution.
3. A **large** reported B5 coupling gain on this corpus should first be
   suspected of marginal leakage into the coupling term, and the orthogonality
   tolerance in section 7 is the check that must be inspected first.

## 13. Deliverables

`PHASE2A_PREREGISTRATION_HASH.json`, `PHASE2A_INPUT_MANIFEST.json`,
`PHASE2A_DATA_IDENTIFIABILITY_CENSUS.json`,
`PHASE2A_TEACHER_CONDITIONALITY.json`,
`PHASE2A_MARGINAL_COUPLING_AUDIT.json`,
`PHASE2A_ATTRIBUTION_AND_REWIRING.json`,
`PHASE2A_LABEL_SEMANTICS.json`, `PHASE2A_VERDICT.json`,
`PHASE2A_SYNTHESIS.md`, and — only if the verdict authorizes it — exactly one
conditional next-stage preregistration.

Every artifact records this file's SHA-256, the repository commit, all seeds,
all software versions, the label fields read, the affinity read count, the
inference unit, and the numerical tolerances actually achieved.
