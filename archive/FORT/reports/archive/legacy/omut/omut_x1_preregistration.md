# OpenMut `OMUT-X1` preregistration

**Frozen before formal execution:** 2026-07-28.
**Predecessor:** `OMUT-X0`, verdict
`OMUT_X0_EVIDENCE_INADEQUATE_STOP`.
**Stage boundary:** label-free evidence recovery only.

## 1. Question

X0 found 225 ChEMBL `Ki`/`Kd` WT-versus-single-mutant candidates with at
least four shared ligands, but no WT assay exposed an explicit sequence in
`variant_sequence`. X1 asks a narrower source-recovery question:

> Can source-native ChEMBL assay descriptions and assay groups identify exact
> WT/mutant construct pairs under the same documented assay context, without
> reading an affinity value?

ChEMBL's deposition guidance defines an assay record as one assay instance,
asks descriptions to include isoform/mutation details, and provides
`ASSAY_GROUP` for assays considered comparable by the depositor. The ChEMBL
variant-curation paper also warns against assuming that a canonical target
means a full-length unmodified WT construct. X1 therefore uses descriptions
only when they state an explicit construct span or explicitly state
full-length protein.

## 2. Discovery contamination and confirmation boundary

A pre-preregistration discovery probe inspected one component:
`P15056`, `V600E`, `CHEMBL5145`, `Kd`. It found paired descriptions in one
document for BRAF `S429` to `E741`, with and without V600E. That accession is
permanently excluded from every X1 primary count, adequacy threshold, and
family threshold. It remains a positive software fixture only.

No affinity, relation, censoring, unit, observed difference, empirical
variance, or empirical MDE value was requested or read in the discovery
probe.

## 3. Frozen sources and fields

1. X0's release-keyed, projected ChEMBL 37 variant census.
2. ChEMBL activity records projected to:
   `activity_id`, `assay_chembl_id`, `assay_variant_accession`,
   `assay_variant_mutation`, `target_chembl_id`, `molecule_chembl_id`,
   `document_chembl_id`, and `standard_type`.
3. ChEMBL assay records projected to:
   `assay_chembl_id`, `target_chembl_id`, `document_chembl_id`,
   `variant_sequence`, `description`, `assay_group`, `assay_type`,
   `confidence_score`, and `src_assay_id`.
4. ChEMBL document records projected to:
   `document_chembl_id`, `doi`, and `pubmed_id`.
5. Versioned UniProt records and the frozen local KLIFS snapshot used by X0.

Numeric activity fields, relations, units, comments, validity flags,
censoring, pH, temperature, and endpoint values are forbidden. Free-text
assay descriptions are retained as evidence text but are not parsed for any
numeric outcome.

All network responses are stored in a release-keyed projected cache with
request parameters, allowed fields, and response SHA-256. Cache reuse fails
closed on a release, schema, request, or hash mismatch.

## 4. Candidate set

The candidate set is every X0 ChEMBL `(accession, single-substitution token,
target_chembl_id, endpoint)` group with at least four distinct variant-side
molecules. Endpoints remain separate and only `Ki` and `Kd` are primary.

Every candidate receives one explicit disposition. The threshold cannot be
lowered and candidates cannot be added from description-only mutation mining.

## 5. Exact description pairing

For a variant and WT activity on the same molecule:

1. target, endpoint, and nonempty primary document identifier must match;
2. the variant assay sequence must be consistent with exactly the named
   substitution against the versioned UniProt sequence;
3. the WT activity must have no ChEMBL variant mutation and its assay must
   have no variant sequence;
4. both descriptions must resolve to the same single construct span:
   - an explicit residue/amino-acid range in both descriptions; or
   - explicit `full length`/`full-length` wording in both descriptions;
5. the mutation position must lie inside the span and the UniProt reference
   residue must equal the token's WT residue;
6. after removing only the exact mutation token, the words `wild type`,
   `wild-type`, `WT`, `mutant`, `mutation`, and `variant`, punctuation and
   repeated whitespace, the remaining description signatures must be
   nonempty and exactly equal; alternatively, both nonempty source-native
   `assay_group` values must be equal;
7. assay type must match.

There is no fuzzy matching, learned text similarity, manually selected
exception, DOI-neighbor imputation, target-ID inference, or source-row
proximity rule. Descriptions containing `unknown origin` cannot establish an
exact construct unless an explicit span or explicit full-length statement is
also present.

A recovered shared ligand must have at least one qualifying assay pair. A
component is keyed by accession, mutation token, target, endpoint, and
construct span. To prevent construct expansion from inflating independent
biological sample count, only the span with the largest recovered ligand
count is retained per `(accession, token, target, endpoint)`, with lexical
span order as the frozen tie-breaker.

## 6. Family and independence accounting

KLIFS group is authoritative for mapped kinases. Other proteins use the first
normalized UniProt similarity-family statement, or `unclassified`.

Component, accession, broad-family, document, and construct-span counts are
reported separately. Rows, ligands, assay pairs, papers, databases, endpoints,
and construct spans do not create additional independent base proteins.

The discovery accession `P15056` is excluded before all primary counts.

## 7. Gates

| gate | condition |
| --- | --- |
| `X1_X0_BOUND` | X0 result is present, firewall-safe, and has the frozen inadequate verdict |
| `X1_NO_OUTCOME_MATERIALIZED` | no forbidden activity or empirical-result field occurs in requests, caches, or result |
| `X1_CENSUS_BOUND` | release and all 119,801 projected variant rows match X0 |
| `X1_WT_QUERIES_COMPLETE` | every k-capable target/endpoint/molecule query is complete |
| `X1_ASSAY_DOCUMENT_FETCH_COMPLETE` | every requested assay and document locator is resolved |
| `X1_ALL_CANDIDATES_DISPOSITIONED` | every frozen X0 k-capable group has exactly one status |
| `X1_EXACT_PAIR_RULE_ENFORCED` | every recovered ligand has exact document, endpoint, span, context, and sequence evidence |
| `X1_DISCOVERY_ACCESSION_EXCLUDED` | `P15056` contributes zero primary components |
| `X1_FAMILY_MAPPED` | every primary accession has KLIFS, UniProt, or explicit unclassified disposition |
| `X1_PRIMARY_TOPOLOGY_ADEQUATE` | at least 25 non-BRAF components, at least six broad families, and largest accession share at most 50% |

Verdicts:

- all gates pass: `OMUT_X1_DESCRIPTION_REGISTRY_ADEQUATE`;
- execution and firewall gates pass but primary topology fails:
  `OMUT_X1_DESCRIPTION_REGISTRY_INADEQUATE_STOP`;
- any execution, binding, cache, or firewall gate fails:
  `OMUT_X1_INCOMPLETE_STOP`.

An adequate verdict unlocks `OMUT-I0` outcome and empirical-power
preregistration only. It does not authorize a predictor, representation, or
mechanism claim.

## 8. Deliverables

- `research/omut_x1.py`
- `tests/test_omut_x1.py`
- `reports/active/omut_x1.json`
- `reports/active/omut_x1_decision.md`
