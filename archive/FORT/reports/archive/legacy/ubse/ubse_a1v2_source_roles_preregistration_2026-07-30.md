# UBSE-A1-v2 source-role and membership preregistration

Date: 2026-07-30  
Status: frozen metadata construction protocol; no coordinate/event access  
Decision boundary: this protocol may freeze candidate roles only. It cannot
authorize coordinate-body download, event extraction, model training,
affinity access, confirmation scoring, or sealed access.

## 1. Purpose

A1-v1 cannot identify residue-functional-group coupling and its former
88-panel audit is no longer fresh. A1-v2 therefore requires three roles to be
fixed before any new coordinate or event value is read:

1. `A1-R`: cross-deposition and cross-extractor/manual-review reliability;
2. `A1-S`: source fit and development roles derived from the historical A0
   source after chemical-neighbour remediation;
3. `A1-C`: a new confirmation candidate role whose binary, typed, and
   residue-FG labels have never been read.

The source-role construction is label blind. A metadata pass does not imply
that A1-C is confirmation-clean: inherited supervised-model membership,
pocket/template similarity, predicted-structure membership, locator
resolution, coordinate completeness, and event topology are separate hard
gates.

## 2. Frozen inputs and firewall

Expected SHA-256 values:

- `closed_registry.parquet`:
  `7905e4edf88073f564baa4b2d4fb50d496432bc4e15e97cccbfa0766b1b0638d`;
- `ubse_a0c_3d_event_sources_v2.parquet`:
  `adc72f142e515c47ea18d20d7af08f6a434a30202a1483dcb362115062a068d5`;
- `BioLiP.txt.gz`:
  `c92229bbc8c55c3bd84a9813c3e278ba62f4cfa44e6315cc98d9bf63ed64b6ec`;
- accepted current P0A result:
  `d4be1fe9fa1a87faa5db8f390587f742a629d634ef35d9364e89cb92163a7a61`.

The only columns that the source-role preflight may project from the closed
registry are:

```text
target_key
sequence
accession
pdb_id
receptor_chain
ligand_ccd
ligand_chain
ligand_serial
pubmed
conn
scaffold
heavy_atoms
```

It must not project either binding-residue column, `affinity_presence`, any
coordinate body, any PLINDER interaction field, any typed event, any P0A
prediction, or any affinity/outcome field. It must record
`binding_residue_fields_loaded=false`, `event_fields_loaded=false`,
`coordinate_files_read=false`, `affinity_fields_loaded=false`, and
`sealed_test_consumed=false`.

## 3. Frozen similarity definitions

Protein closure uses upper-cased exact sequences, unique 4-mers, and

\[
\operatorname{containment}(a,b)
=\frac{|K_4(a)\cap K_4(b)|}
{\min(|K_4(a)|,|K_4(b)|)}.
\]

`containment >= 0.4` is a conflict. Equality is not admissible.

Ligand closure uses RDKit `2023.09.6`, Morgan radius 2, 2,048 bits,
`useChirality=False`, and `useFeatures=False`. Full exact Tanimoto is used;
ANN or sampled comparisons are forbidden. `Tanimoto >= 0.5` is a conflict,
including equality.

All PDB identifiers are lower-cased. All other identifiers are stripped and
compared as exact strings. Empty identifiers are ineligible.

Selection seed is `1729`. Hash order uses SHA-256 over the ASCII serialization

```text
<role salt>|1729|<ordered identity fields>
```

and then the ordered identity fields as a collision tie-break. The role salts
are `UBSE-A1R-V2` and `UBSE-A1C-V2`.

## 4. A1-R reliability role

### 4.1 A0 closure

Start from the safe registry projection. Exclude:

- every target whose sequence has containment at least `0.4` to any A0
  target;
- every exact A0 PDB, PubMed, connectivity, or scaffold;
- every ligand whose ECFP4 Tanimoto to any A0 ligand is at least `0.5`.

The expected intermediate counts before A1-R topology are:

- A0: 957 targets, 2,833 PDBs, 871 PubMeds, 2,668 connectivities, and
  897 scaffolds;
- direct A0 target-neighbour closure: 10,849 registry targets excluded;
- exact-resource-closed pool: 39,057 rows and 28,841 targets;
- after the A0 ECFP4 gate: 33,373 rows and 24,670 targets.

### 4.2 Unit definition

An A1-R target is eligible only when a three-instance unit can be formed:

- two `correct` instances share exact
  `(target_key, conn, ligand_ccd)`;
- the two correct instances have different PDBs and different PubMeds;
- one `wrong` instance has the same target and a different connectivity;
- the wrong PDB and PubMed are each different from both correct instances;
- all three rows survive every A0 closure above.

For a target with several legal triples, choose one by:

1. same-scaffold wrong ligand before different-scaffold wrong ligand;
2. smaller absolute correct-versus-wrong heavy-atom difference;
3. role-salted SHA-256 order of the complete three-instance identity;
4. lexical identity order.

The expected topology is 222 eligible correct units over 155 targets.
Order the 155 chosen target triples by the A1-R target hash and freeze the
first 153 as primary A1-R targets. The two unselected targets are not
post-result reserves.

A1-R events must never enter A1-S fitting. Locator and extractor completion
must retain at least 128 of the 153 frozen targets. Failure gives
`STOP_UBSE_A1R_TYPED_EVENT_TEACHER_UNRELIABLE`; targets may not be replaced
after any coordinate or event value is read.

Agreement from one extractor across PDB/PubMed depositions is called
cross-deposition repeatability, not independent reliability.

## 5. A1-C confirmation-candidate role

### 5.1 Reference closure

Use the union of all A0 instances and all 459 selected A1-R instances as the
reference set. Exclude every registry row conflicting with that union on:

- target sequence containment;
- exact target, PDB, PubMed, connectivity, or scaffold;
- ligand ECFP4 Tanimoto.

Require a parseable single-component ligand and a nonempty locator identity.
The preflight must report the heavy-atom range rather than tune it after
selection.

### 5.2 Conflict-free packing

Collapse exact duplicate coordinate identities. Sort remaining rows by:

1. increasing maximum exact-resource degree over target, PDB, PubMed,
   connectivity, and scaffold;
2. increasing target row degree;
3. increasing ligand target degree;
4. A1-C role-salted SHA-256;
5. lexical coordinate identity.

Greedily accept one instance only if it conflicts with neither the reference
set nor an already accepted instance. Accepted instances must have:

- unique target, PDB, PubMed, connectivity, and scaffold;
- pairwise protein containment below `0.4`;
- pairwise ligand ECFP4 Tanimoto below `0.5`.

Freeze the first 512 accepted instances as A1-C primary and the next 64 as
metadata-only reserves. A reserve may replace a primary only for a
preregistered HEAD-availability or byte-safe locator-resolution failure, in
frozen order, and only before any coordinate body or event value is read.
Event absence, checkerboard absence, performance, family, chemistry, or
model output may never trigger replacement.

A1-C is a one-complex-per-target confirmation candidate role. It is not a
two-ligand target panel. The required residue-FG checkerboard is within one
complex and event type.

## 6. Historical A1-S role remediation

The existing A0 roles are historical and immutable. Before A1-S event
extraction, construct a new metadata-only A1-S role view with precedence:

```text
fit > development > legacy_audit
```

Keep A0 fit fixed. Remove an entire later-role panel if any of its ligands
has ECFP4 Tanimoto at least `0.5` to any earlier retained role. Exact
target/homology, PDB, PubMed, connectivity, and scaffold closure remains
mandatory.

The preregistration-time audit found:

- validation: 7/140 violating ligands in 5 panels, maximum `0.6984127`;
- legacy audit: 12/197 violating ligands in 7 panels, maximum `0.7767857`,
  before applying development-to-legacy-audit precedence.

The final retained panel counts and all cross-role maximum similarities must
be emitted by the preflight. The historical G1 audit remains development
evidence only; it can never become the new paper-level confirmation role.

## 7. Inherited supervised-model membership

The accepted current P0A weights read binding-residue labels from 62,849
rows over 38,781 targets. A1-R and A1-C candidates are drawn from that legal
training complement. Therefore the current P0A is not clean for A1-C.

Before any primary A1-C evaluation, both are mandatory:

1. retrain `P0A-v2` after A1-C is frozen, excluding the full A1-C target
   homology, PDB, PubMed, connectivity, scaffold, and ECFP4-neighbour
   closure from contact supervision;
2. retain a no-P0A A1-v2 arm under the same non-P0A inputs and evaluation.

The current P0A may be reported only as architecture/proposal feasibility or
as a clearly transductive non-primary control. Any other inherited supervised
model must have zero exact/homology target membership overlap with A1-C.

Predicted monomer, cofold, pocket, and template sources remain pending. Their
training cutoff, template use, exact held PDB membership, target sequence
membership, and pocket/structure neighbourhood must be audited before they
enter a primary A1-C arm.

## 8. Locator and extractor gates

After metadata roles are frozen, extend the A0C byte scanner. It may decode
only BioLiP columns `0,1,4,5,6,18,19,20` and must byte-skip affinity columns
`13-16`. It must recover BioLiP column 20 as `mmcif_auth_seq_id` by a strict
unique join; filename serial and auth sequence ID remain distinct.

Every retained instance must have:

- one unique raw join;
- a nonempty scalar `mmcif_auth_seq_id`;
- a unique mmCIF coordinate-instance identity;
- one frozen official RCSB URL.

Before coordinate access, separately freeze:

- primary PLIP package/commit and all dependency hashes;
- a genuinely independent second extractor, or a frozen manual chemistry
  sample and adjudication protocol;
- functional-group SMARTS, graph-symmetry collapse, altloc/occupancy, model,
  water, metal, covalent-ligand, missing-residue, chain, and assembly rules;
- selected URL list and download/body hash contract.

## 9. Gates and decisions

`SR-0 Identity/firewall`

- all frozen hashes and versions match;
- only safe columns were projected;
- no coordinate, binding-residue, event, affinity, confirmation, or sealed
  value was read.

`SR-1 A1-R topology`

- exactly 153 primary targets and 459 instances;
- expected 155-target candidate pool is reproduced;
- each target satisfies the two-correct-plus-wrong definition;
- at least 128 targets later complete locator and extractor resolution.

`SR-2 A1-C topology`

- exactly 512 primary and at least 64 ordered reserves;
- primary target/PDB/PubMed/connectivity/scaffold counts are each 512;
- all exact, target-containment, and ECFP4 conflict counts are zero.

`SR-3 A1-S remediation`

- all retained cross-role exact and similarity conflicts are zero;
- panel removal follows the frozen precedence only.

`SR-4 Inherited membership`

- every primary supervised component has zero A1-C target/homology
  membership;
- current P0A must fail this gate until P0A-v2 is retrained.

`SR-5 Locator/extractor readiness`

- all strict joins and locator uniqueness checks pass;
- extractor/manual-review, functional-group, and coordinate contracts are
  frozen before body access.

`SR-6 Freshness and power`

- at least 400 A1-C primary targets remain coordinate complete;
- after one-time frozen extraction, at least 200 independent A1-C targets
  contain a legal within-complex residue-FG checkerboard;
- at least four event types each cover 50 A1-C targets.

The old threshold of 40 confirmation checkerboards is withdrawn as
underpowered. With 512 independent units, a directional-accuracy truth of
0.60 has only about 64% power to place a lower bound above 0.55; all
directional results must therefore report the realized independent
checkerboard count and a power-qualified interval.

If `SR-0` through metadata `SR-3` pass while current P0A membership fails, the
only allowed decision is:

`FREEZE_A1V2_METADATA_ROLES_RETRAIN_P0A_AND_COMPLETE_MEMBERSHIP_GATES`

Any topology/firewall failure gives:

`STOP_UBSE_A1V2_SOURCE_ROLES_OR_CONFIRMATION_INADEQUATE`

No metadata result can directly return an event-extraction or model-training
authorization.

