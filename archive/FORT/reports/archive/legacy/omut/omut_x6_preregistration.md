# OpenMut `OMUT-X6` preregistration

**Frozen before execution:** 2026-07-28.
**Predecessor:** `OMUT-X5`,
`OMUT_X5_SUPPLEMENT_CONSTRUCT_REGISTRY_INADEQUATE_STOP`.
**Stage boundary:** label-free ChEMBL reagent-locator inventory.

## 1. Question

The ChEMBL variant curation guidance identifies exact sequences and supplier
catalog numbers as preferred construct evidence. X0-X5 exhausted explicit
assay descriptions and accessible papers for the current candidates. X6
asks:

> Do WT assay descriptions for the frozen near-exact candidates contain
> enough exact official reagent locators to justify a source-native
> supplier-registry verification stage?

X6 does not query a supplier, follow a URL, fetch a product page, read an
abstract/body/supplement, or materialize any activity outcome.

## 2. Frozen source and candidates

1. The X5-bound X2 near-pair candidate reconstruction.
2. The corresponding ChEMBL 37 assay records and descriptions.
3. The existing `P15056` discovery exclusion.

Only WT-side assay descriptions belonging to a frozen near pair are scanned.

## 3. Exact locator patterns

The following case-insensitive labels may introduce a locator:

```text
catalog, catalogue, cat., product, product no., item, item no., SKU
```

The locator token must contain at least one digit, be 3-40 characters, and
contain only letters, digits, period, underscore, slash, or hyphen.

An explicit HTTPS URL is also a locator. Bare assay-platform names, company
names, gene names, mutation identifiers, database accessions, paper
citations, phone numbers, concentrations, and unlabelled number strings are
not exact locators.

Recognized supplier/context names are projected only to make a locator
actionable:

```text
BPS Bioscience, Carna Biosciences, DiscoverX, Eurofins, Invitrogen,
Millipore, PerkinElmer, Promega, Reaction Biology, SignalChem,
Thermo Fisher, Upstate
```

A labeled token without a supplier or URL is retained as
`unresolved_supplier`, never counted as actionable.

## 4. Candidate and topology rules

A candidate is `actionable_locator_k4` only when:

- it is non-BRAF and reference-resolved;
- it has at least four X2 near-pair ligands;
- at least one WT assay used by those ligand pairs contains an exact locator;
- that locator has an HTTPS URL or a recognized supplier context.

One locator may support multiple ligands in the same candidate, but each
candidate is counted once. No construct is inferred and no X1/X5 component
is added at X6.

## 5. Gates

| gate | condition |
| --- | --- |
| `X6_X5_BOUND` | X5 result hash, preregistration, firewall, and inadequate verdict bind |
| `X6_CANDIDATES_BOUND` | X2 candidate reconstruction and ChEMBL assay fetch are complete |
| `X6_WT_ONLY` | every scanned assay is a WT-side assay in a frozen near pair |
| `X6_LOCATOR_RULE_ENFORCED` | every actionable locator has a frozen label/URL form and supplier/URL context |
| `X6_DISCOVERY_ACCESSION_EXCLUDED` | `P15056` contributes zero actionable candidates |
| `X6_NO_OUTCOME_MATERIALIZED` | no activity outcome or source body is requested or retained |
| `X6_EXTERNAL_VERIFICATION_TOPOLOGY` | at least 25 actionable `k>=4` candidates across at least six broad families and no accession supplies more than 50% |

Verdicts:

- all gates pass: `OMUT_X6_REAGENT_VERIFICATION_FEASIBLE`;
- execution/firewall gates pass but locator topology fails:
  `OMUT_X6_REAGENT_VERIFICATION_INSUFFICIENT_STOP`;
- any binding, fetch, locator, or firewall gate fails:
  `OMUT_X6_INCOMPLETE_STOP`.

## 6. Deliverables

- `research/omut_x6.py`
- `tests/test_omut_x6.py`
- `reports/active/omut_x6.json`
- `reports/active/omut_x6_decision.md`
