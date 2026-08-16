# OpenMut `OMUT-X2` preregistration

**Frozen before execution:** 2026-07-28.
**Predecessor:** `OMUT-X1`, verdict
`OMUT_X1_DESCRIPTION_REGISTRY_INADEQUATE_STOP`.
**Stage boundary:** label-free source-accessibility topology only.

## 1. Question

X1 recovered four exact non-discovery components and found 25 additional
components with one to three exact paired ligands. X2 asks:

> After requiring all non-construct evidence to match exactly, are enough
> unresolved construct records linked to openly accessible primary full text
> or supplements to make deterministic, quote-bound construct recovery worth
> executing?

X2 is an optimistic source-reachability gate. A pass authorizes full-text
evidence extraction, not affinity reading, power analysis, representation
fitting, or prediction.

## 2. Frozen candidate and discovery boundary

The candidate set is exactly the 279 X1 ChEMBL groups. `P15056` remains a
discovery-only accession and contributes zero primary count.

For each shared molecule, an unresolved pair must already match on:

- molecule, target, endpoint, and nonempty primary document;
- exact single-substitution variant sequence against versioned UniProt;
- WT activity with no variant mutation and WT assay with no variant sequence;
- assay type;
- exact normalized context signature after removing only the named mutation,
  WT/mutant markers, and explicit construct-span/full-length phrases; or a
  shared nonempty source-native assay group.

Pairs with conflicting explicit spans, different descriptions after frozen
normalization, different documents, other mutations, missing locators, or
unresolved reference sequence are excluded.

There is no fuzzy text similarity, language-model decision, assay-ID
adjacency, target-name imputation, DOI-neighbor propagation, or manual
exception.

## 3. Source availability

ChEMBL document records are projected to:

```text
document_chembl_id, doi, pubmed_id, title, journal, year, patent_id
```

For journal articles with DOI or PMID, the Europe PMC REST API is queried with
an exact identifier. The result records only article identity, PMCID,
`isOpenAccess`, `inEPMC`, and supplement availability. A document is
full-text-accessible only when Europe PMC reports both open access and an EPMC
full-text record. Search-result snippets and abstracts do not count as full
text.

Patent documents are reported separately and do not count as accessible in
this gate. Publisher pages, shadow libraries, browser-rendered snippets, and
manually supplied copies are excluded from the formal X2 count.

All ChEMBL and Europe PMC responses are cached by source version or exact
request, projected schema, and SHA-256.

## 4. Recoverable upper bound

For each `(accession, mutation, target, endpoint)`:

1. retain the X1 exact paired ligands;
2. add distinct unresolved ligands whose qualifying primary document is
   full-text-accessible;
3. count the component as source-recoverable when the union has at least four
   ligands.

This is an upper bound: full-text availability does not prove that the
construct is stated. Multiple documents, assays, ligands, spans, or endpoints
do not create additional independent mutations.

Families use the same KLIFS/UniProt mapping as X1. `P15056` is removed before
all primary calculations.

## 5. Gates

| gate | condition |
| --- | --- |
| `X2_X1_BOUND` | X1 frozen result, preregistration, firewall, and inadequate verdict bind |
| `X2_NO_OUTCOME_MATERIALIZED` | no numeric activity, relation, censor, unit, empirical variation, or observed difference is requested or retained |
| `X2_ALL_CANDIDATES_DISPOSITIONED` | all 279 X1 groups receive one exact/near/unmatched status |
| `X2_NEAR_PAIR_RULE_ENFORCED` | every near pair satisfies all frozen non-construct checks and has no conflicting explicit span |
| `X2_DOCUMENT_METADATA_COMPLETE` | every near-pair ChEMBL document locator is resolved |
| `X2_EPMC_QUERIES_COMPLETE` | every near-pair journal identifier has a complete exact Europe PMC disposition |
| `X2_DISCOVERY_ACCESSION_EXCLUDED` | `P15056` contributes zero primary source-recoverable components |
| `X2_FAMILY_MAPPED` | every primary source-recoverable accession has a family disposition |
| `X2_SOURCE_TOPOLOGY_ADEQUATE` | at least 25 non-BRAF source-recoverable components, at least six broad families, and largest accession share at most 50% |

Verdicts:

- all gates pass: `OMUT_X2_FULLTEXT_RECOVERY_FEASIBLE`;
- execution/firewall gates pass but source topology fails:
  `OMUT_X2_FULLTEXT_RECOVERY_INSUFFICIENT_STOP`;
- any execution, cache, projector, or binding gate fails:
  `OMUT_X2_INCOMPLETE_STOP`.

## 6. Deliverables

- `research/omut_x2.py`
- `tests/test_omut_x2.py`
- `reports/active/omut_x2.json`
- `reports/active/omut_x2_decision.md`
