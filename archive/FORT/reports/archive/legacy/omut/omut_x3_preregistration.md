# OpenMut `OMUT-X3` preregistration

**Frozen before execution:** 2026-07-28.
**Predecessor:** `OMUT-X2`, verdict
`OMUT_X2_FULLTEXT_RECOVERY_INSUFFICIENT_STOP`.
**Stage boundary:** label-free licensed-open-version discovery only.

## 1. Question

X2 found a theoretical 110-component near-exact topology but only nine EPMC
open-full-text documents, yielding 16 source-recoverable components. X3 asks:

> Do exact-DOI OpenAlex records expose enough additional accepted or
> published open versions with explicit licenses to clear the frozen 25
> component / six family source-accessibility gate?

X3 does not fetch or parse a PDF, XML body, supplement, abstract, snippet, or
activity outcome.

## 2. Frozen sources

1. X2's 87 near-pair ChEMBL documents and frozen exact/near pair rules.
2. X2's Europe PMC source dispositions.
3. OpenAlex Works API queried by exact normalized DOI.

Each OpenAlex response is cached by DOI, response SHA-256, and a projected
location schema.

## 3. Admissible open location

An OpenAlex location counts only when all conditions hold:

- the work DOI exactly equals the requested DOI;
- `is_oa` is true;
- version is `acceptedVersion` or `publishedVersion`;
- a HTTPS PDF or landing-page URL is present;
- license is one of:
  `cc0`, `cc-by`, `cc-by-sa`, `cc-by-nc`, `cc-by-nc-sa`,
  `cc-by-nd`, or `cc-by-nc-nd`;
- the host is not ResearchGate, Academia.edu, Sci-Hub, LibGen, or a search
  result/cache service.

Submitted/preprint versions, `bronze` links without a license, unknown
licenses, closed locations, DOI resolver links without an open artifact, and
manually discovered URLs are excluded.

An X2 EPMC-accessible document remains accessible. Otherwise a document
becomes accessible only through at least one admissible OpenAlex location.

## 4. Topology

The exact X2 component calculation is repeated with the union of EPMC and
admissible OpenAlex documents. `P15056` remains discovery-only and is removed
before all primary counts.

## 5. Gates

| gate | condition |
| --- | --- |
| `X3_X2_BOUND` | X2 result, preregistration, firewall, and insufficient verdict bind |
| `X3_NO_OUTCOME_MATERIALIZED` | no activity outcome or document body is requested or retained |
| `X3_ALL_DOIS_DISPOSITIONED` | every near-pair DOI has one exact OpenAlex result, not-found result, or explicit no-DOI status |
| `X3_LICENSE_RULE_ENFORCED` | every new accessible location satisfies DOI, OA, version, URL, license, and host rules |
| `X3_DISCOVERY_ACCESSION_EXCLUDED` | `P15056` contributes zero primary source-recoverable components |
| `X3_FAMILY_MAPPED` | every primary source-recoverable accession has a family disposition |
| `X3_SOURCE_TOPOLOGY_ADEQUATE` | at least 25 non-BRAF source-recoverable components, at least six broad families, and largest accession share at most 50% |

Verdicts:

- all gates pass: `OMUT_X3_LICENSED_FULLTEXT_RECOVERY_FEASIBLE`;
- execution/firewall gates pass but topology fails:
  `OMUT_X3_LICENSED_FULLTEXT_RECOVERY_INSUFFICIENT_STOP`;
- any execution, cache, DOI, or license gate fails:
  `OMUT_X3_INCOMPLETE_STOP`.

## 6. Deliverables

- `research/omut_x3.py`
- `tests/test_omut_x3.py`
- `reports/active/omut_x3.json`
- `reports/active/omut_x3_decision.md`
