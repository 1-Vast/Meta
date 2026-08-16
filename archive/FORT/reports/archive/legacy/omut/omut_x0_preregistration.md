# OpenMut `OMUT-X0` preregistration

**Frozen before execution:** 2026-07-28.
**Predecessor:** `OMUT-D1`, verdict `OMUT_D1_TOPOLOGY_ADEQUATE`.
**Authority:** `task.md`, Route A, `OMUT-X0`.

**Implementation amendment A1 (2026-07-28, after one execution was manually
terminated before producing any result):** the first implementation followed ChEMBL
`next` links in memory and emitted no progress until completion. After more than 13
minutes in an unobservable network call chain, the process was terminated; no X0
result or projection cache existed and no outcome value had been requested. The
replacement stores each completed, label-free ChEMBL page in a release-keyed cache,
verifies the projected schema and content hash on reuse, and logs progress. This
changes transport observability and recoverability only. Sources, fields, candidate
set, gates, thresholds, and verdict rules below are unchanged.

**Implementation amendment A2 (2026-07-28, after the first cached execution failed
before producing a result):** the run completed all 143 deferred WT queries and then
encountered ChEMBL accession `B1LRJ1`. UniProt returns this record as inactive and
deleted, with no sequence. The implementation had incorrectly applied the
`X0_UNIPROT_COMPLETE` fatal requirement for the 17 frozen D1 accessions to every
additional accession discovered in the full ChEMBL census. The replacement keeps
the requirement fatal for all D1 accessions, explicitly records unresolved
additional accessions, and excludes their components from exact-sequence primary
topology. No source, field, candidate, gate, threshold, or verdict rule is changed.
No result was produced and no numeric outcome was requested or materialized by the
failed execution.

## 1. Question and boundary

`OMUT-X0` asks whether the label-free topology reported by D1 survives an
evidence-bound reconstruction:

> Are the candidate WT-to-single-substitution components real, sequence-resolved
> constructs with endpoint, assay locator, document provenance, and enough shared
> ligands to justify opening `OMUT-I0`?

This is still a no-outcome stage. It may read sequence, construct, ligand identity,
endpoint *presence*, assay identity/context, document identity, institution, source,
and source-row locators. It must not materialize or retain a numeric affinity,
relation, censor, unit, temperature, pH, observed difference, empirical variance,
reliability estimate, or empirical MDE. Endpoint presence is produced as a Boolean
inside the blind projector; the numeric cell is never parsed or returned.

An X0 pass is evidence availability, not interaction, mechanism, or prediction
evidence. It cannot authorize a neural model. It unlocks I0 only.

DAVIS-Complete remains excluded under the immutable `OMUT-F0`
`preserve_confirmation` policy. PLATINUM remains excluded because its data rights are
unresolved. MdrDB is a source index only and is not an independent outcome lineage.

## 2. Frozen sources

1. The SHA-256-verified local BindingDB 202607 Articles archive used by D1.
2. UniProtKB REST records for exactly the D1 `k=4` accessions. The response release
   header, accession, entry version, sequence length, sequence SHA-256, protein name,
   organism, and similarity-family text are recorded. Full sequences are used in
   memory but are not written to the result.
3. The local KLIFS 2026-07-22 kinase-information snapshot. A UniProt accession match
   supplies the frozen kinase group and family. Non-KLIFS proteins use the first
   UniProt `SIMILARITY` family statement; absence is reported as `unclassified`.
4. The ChEMBL web API release exposed by the D0-frozen service. All 119,801
   variant-linked activities are scanned through a server-side `only=` projector.

No external annotation, paper interpretation, LLM retrieval, or manually selected
gold subset may override a sequence or provenance failure.

## 3. BindingDB evidence reconstruction

### 3.1 Candidate set

The candidate set is exactly the 62 `(accession, mutation-token)` pairs that cleared
`k=4` in D1. X0 may disposition them but may not add candidates by relaxing D1's
regex, lowering `k`, or treating multi-substitution names as single substitutions.

### 3.2 Blind projected fields

The projector retains only:

- ligand InChIKey, BindingDB monomer ID, target name, curation source;
- article DOI, BindingDB entry DOI, PMID, institution, and one-based source row;
- number of target chains;
- for each chain, the BindingDB sequence and Swiss-Prot/TrEMBL accession;
- Boolean presence for `Ki`, `Kd`, `IC50`, and `EC50`.

Numeric endpoint cells are tested for empty/non-empty status only. Their contents are
never parsed, copied, hashed, compared, logged, or written. All relation, censor,
unit, pH, temperature, kinetic, and other outcome fields are forbidden.

### 3.3 Unique mutation-chain and construct validation

For token `W p M` assigned by D1 to accession `A`:

1. UniProt sequence `A[p]` must equal `W`.
2. The target row must contain exactly one D1 token.
3. Exactly one chain carrying accession `A` must be compatible.
4. With an explicit `[start-end]` range, chain length must equal the range length,
   and the chain sequence must equal the corresponding UniProt slice at every
   position except the one token position, where it must contain `M`.
5. Without a range, chain length must equal the canonical UniProt length and satisfy
   the same one-difference condition.
6. Slash-separated or repeated ranges are tried independently. More than one valid
   chain/range resolution is ambiguous and is rejected.
7. The exact WT construct is reconstructed by changing only `M` back to `W`. Ordered
   accessions and SHA-256 hashes of all chains define the construct key. A WT
   observation must have zero D1 tokens and match this exact reconstructed key.

This deliberately rejects numbering conventions, isoforms, complexes, truncations,
or additional substitutions that cannot be resolved without guessing. In
particular, a mutation present on another chain is not credited to chain 1.

### 3.4 Components

Components are keyed by accession, token, exact WT construct key, and endpoint.
`Ki` and `Kd` are separate. `IC50` and `EC50` are registry-only sensitivity
endpoints and cannot clear X0. A shared ligand requires an evidence-bound mutant row
and exact-WT row with the same InChIKey and endpoint.

For every primary component the output records shared-ligand count, distinct
documents, documents that contain both sides of at least one shared-ligand pair,
source-row locators, and whether a source-native assay identifier exists. BindingDB
has no source-native assay identifier in the frozen schema, so this limitation is
reported rather than replaced with an invented key.

## 4. Family and concentration audit

Candidate and verified components are counted by broad family and accession. KLIFS
group/family is authoritative where present. For other accessions, the normalized
UniProt family statement is retained without manually merging unrelated families.

The output reports the largest accession and broad-family share. Row, ligand, seed,
and paper expansion do not increase the independent biological component count.

## 5. ChEMBL census and deferred WT comparison

### 5.1 Full variant census

The activity projection is:

```text
activity_id, assay_chembl_id, assay_variant_accession,
assay_variant_mutation, target_chembl_id, molecule_chembl_id,
document_chembl_id, standard_type
```

with `assay_variant_mutation__isnull=false`, pages of 1,000, and no row cap below
the API-declared total. ChEMBL's unavoidable `type` alias is tolerated present and
dropped unread exactly as in D1. Completion requires fetched rows to equal the
declared total and every page to satisfy the frozen projection.

Only `Ki` and `Kd` groups with at least four distinct variant-side molecules can
possibly clear X0. For each such `(accession, mutation, target, endpoint)` group,
the deferred WT query requests non-variant activities for the same target, endpoint,
and only those molecule IDs, in deterministic chunks. This restriction is lossless
for the frozen `k=4` question.

Assay records are fetched for candidate variant and WT observations through the
non-affinity assay resource. The nested `variant_sequence` and categorical assay
context are retained. Variant constructs are checked against UniProt by the same
single-difference rule. A WT assay with no explicit sequence is not called an exact
construct merely because its target identifier is canonical.

Document identifiers are passed through the `OMUT-F0` provenance function. A ChEMBL
document and a BindingDB DOI/PMID naming the same primary paper are one provenance
unit, not two database replications.

Network or API-filter failure is fail-closed. A partial census cannot be extrapolated.

## 6. Gates and verdicts

All gates below are frozen before execution:

| gate | condition |
| --- | --- |
| `X0_D1_BOUND` | D1 result and preregistration bindings pass and its verdict is `OMUT_D1_TOPOLOGY_ADEQUATE` |
| `X0_NO_OUTCOME_MATERIALIZED` | projector schema and result contain no numeric affinity, relation, censor, unit, temperature, pH, observed-difference, empirical-variance, reliability, or empirical-MDE field |
| `X0_UNIPROT_COMPLETE` | every D1 `k=4` accession has a versioned UniProt record and sequence |
| `X0_ALL_D1_CANDIDATES_DISPOSITIONED` | every one of the 62 D1 candidates has exactly one explicit verification status |
| `X0_BINDINGDB_CONSTRUCT_EXACT` | every primary BindingDB component uses a uniquely verified mutant sequence and an exact reconstructed WT construct |
| `X0_FAMILY_MAPPED` | every D1 accession has a KLIFS, UniProt, or explicit `unclassified` family disposition |
| `X0_CHEMBL_VARIANT_CENSUS_COMPLETE` | fetched variant activities equal the live API-declared total |
| `X0_CHEMBL_WT_COMPARISON_COMPLETE` | every ChEMBL `Ki`/`Kd` variant group capable of `k=4` has a complete, filtered WT comparison or an explicit zero result |
| `X0_PROVENANCE_COLLAPSED` | every document count is produced through the F0 provenance rule and cross-database identifiers are not counted twice |
| `X0_PRIMARY_TOPOLOGY_ADEQUATE` | at least 25 exact `Ki`/`Kd` components have at least four shared ligands, at least six broad families are represented, and no single accession contributes more than 50% |
| `X0_ASSAY_EVIDENCE_AVAILABLE` | every primary observation has a source-native assay locator or is excluded from the primary topology |

Verdicts:

- all gates pass: `OMUT_X0_EVIDENCE_REGISTRY_ADEQUATE`;
- execution and firewall gates pass but primary topology or assay evidence fails:
  `OMUT_X0_EVIDENCE_INADEQUATE_STOP`;
- any execution, binding, census, projector, or provenance gate fails:
  `OMUT_X0_INCOMPLETE_STOP`.

The adequacy thresholds are not relaxed after seeing results. Failure is not rescued
with a larger model, embeddings, extra epochs, more seeds, imputation, paper count,
database duplication, or synthetic labels.

## 7. Deliverables

- `research/omut_x0.py`
- `tests/test_omut_x0.py`
- `reports/active/omut_x0.json`
- `reports/active/omut_x0_decision.md`
