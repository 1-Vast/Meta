# OpenMut `OMUT-X3C` preregistration

**Frozen before execution:** 2026-07-28.
**Predecessors:** `OMUT-X2`,
`OMUT_X2_FULLTEXT_RECOVERY_INSUFFICIENT_STOP`; `OMUT-X3`, transport-stopped
without a result.
**Stage boundary:** label-free Crossref license/link discovery only.

## 1. Question

X2 found a theoretical 110-component near-exact topology but only nine EPMC
open-full-text documents, yielding 16 optimistic source-recoverable
components. X3C asks:

> Do exact-DOI Crossref records expose enough additional explicitly
> licensed version-of-record or accepted-manuscript full-text links to clear
> the frozen 25 component / six family source-accessibility gate?

X3C does not fetch or parse a linked PDF, XML, HTML body, supplement,
abstract, snippet, or activity outcome.

## 2. Frozen sources and query

1. X2's 87 near-pair ChEMBL documents and frozen exact/near-pair rules.
2. X2's Europe PMC source dispositions.
3. Crossref REST API `GET /works/{doi}` queried by exact normalized DOI,
   without an invented identity.

Only the response DOI, `license`, and `link` fields are projected. The full
response is not retained. Each projection is cached with the requested DOI,
HTTP-response SHA-256, and canonical projection SHA-256.

## 3. Admissible Crossref location

A Crossref record contributes a new accessible document only when all
conditions hold:

- the returned DOI exactly equals the requested normalized DOI;
- at least one license entry has a recognized Creative Commons URL:
  `CC0`, `CC-BY`, `CC-BY-SA`, `CC-BY-NC`, `CC-BY-NC-SA`,
  `CC-BY-ND`, or `CC-BY-NC-ND`;
- that license applies to the version of record (`vor`) or accepted
  manuscript (`am`), not only text-and-data mining (`tdm`);
- if the license has a start date, it is no later than the run time;
- at least one HTTPS full-text link has content version `vor` or `am`, or
  omits content version; explicit `tdm` content versions are rejected;
- when a link declares `vor` or `am`, an active admissible license must apply
  to that same version;
- that link declares PDF, XML, or HTML content, and its intended application
  is `text-mining` or is absent;
- the host is not a DOI resolver, ResearchGate, Academia.edu, Sci-Hub,
  LibGen, a search/cache service, or a URL-shortening service.

License and link metadata are work-level Crossref deposits and do not prove
that the license is attached to a particular link. X3C therefore reports an
optimistic licensed-link upper bound, never a recovered construct. Unknown
licenses, publisher links without an admissible license, licenses without a
full-text link, delayed licenses not yet active, `similarity-checking` links,
and manually discovered URLs are excluded.

An X2 EPMC-accessible document remains accessible. Otherwise a document
becomes accessible only through an admissible Crossref license-plus-link
record.

## 4. Topology

The exact X2 component calculation is repeated with the union of EPMC and
admissible Crossref documents. `P15056` remains discovery-only and is removed
before all primary counts.

## 5. Gates

| gate | condition |
| --- | --- |
| `X3C_X2_BOUND` | X2 result, preregistration, firewall, and insufficient verdict bind |
| `X3C_NO_OUTCOME_MATERIALIZED` | no activity outcome, abstract, or document body is requested or retained |
| `X3C_ALL_DOIS_DISPOSITIONED` | every near-pair document has one exact Crossref result, not-found result, or explicit no-DOI status |
| `X3C_LICENSE_LINK_RULE_ENFORCED` | every new Crossref-accessible record satisfies the frozen DOI, license, version, date, URL, content, application, and host rules |
| `X3C_DISCOVERY_ACCESSION_EXCLUDED` | `P15056` contributes zero primary source-recoverable components |
| `X3C_FAMILY_MAPPED` | every primary source-recoverable accession has a family disposition |
| `X3C_SOURCE_TOPOLOGY_ADEQUATE` | at least 25 non-BRAF source-recoverable components, at least six broad families, and largest accession share at most 50% |

Verdicts:

- all gates pass: `OMUT_X3C_LICENSED_LINK_RECOVERY_FEASIBLE`;
- execution/firewall gates pass but topology fails:
  `OMUT_X3C_LICENSED_LINK_RECOVERY_INSUFFICIENT_STOP`;
- any execution, cache, DOI, or admissibility gate fails:
  `OMUT_X3C_INCOMPLETE_STOP`.

## 6. Deliverables

- `research/omut_x3c.py`
- `tests/test_omut_x3c.py`
- `reports/active/omut_x3c.json`
- `reports/active/omut_x3c_decision.md`
