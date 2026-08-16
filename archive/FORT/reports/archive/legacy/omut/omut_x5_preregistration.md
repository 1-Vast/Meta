# OpenMut `OMUT-X5` preregistration

**Frozen before execution:** 2026-07-28.
**Predecessor:** `OMUT-X4`,
`OMUT_X4_CONSTRUCT_REGISTRY_INADEQUATE_STOP`.
**Stage boundary:** official Europe PMC supplementary-archive construct
recovery.

## 1. Question

X4 verified nine EPMC main-text documents but recovered no new exact
candidate-level construct relation. All nine records declare supplementary
files. X5 asks:

> Do their official Europe PMC supplementary archives contain explicit
> common WT-construct evidence for enough candidates to clear the unchanged
> 25 component / six family / 50% concentration gate?

X5 does not add or retain affinity outcomes.

## 2. Frozen source set and transport

The only requests are:

```text
https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles
```

for the nine X4 transport-complete EPMC records. Europe PMC documents that the
endpoint returns available supplementary files as a ZIP archive. No
publisher search, DOI landing page, Crossref expansion, manually substituted
URL, or non-X4 document may enter.

For each response, record the PMCID, URL, final host, HTTP status, media type,
byte length, response SHA-256, ZIP validity, and a member inventory with
member name, uncompressed size, suffix, and SHA-256. Reject path traversal,
encrypted members, nested archives, executables, and members over 100 MiB.
Raw archives and members are not retained after projection.

## 3. Admissible text-bearing members

Only `.pdf`, `.xml`, `.html`, `.htm`, and `.txt` members are parsed. Tables in
CSV/Excel form, images, chemical structure files, and unknown binaries are
inventoried but not read. This deliberately prevents activity matrices from
entering through a table parser.

- XML/HTML/TXT: parse text blocks with section provenance.
- PDF: extract page-numbered text with `pdfplumber`; if a candidate is
  accepted from a PDF, render the cited page with Poppler and visually verify
  the construct statement before the verdict is frozen.

## 4. Outcome firewall

The X4 forbidden text patterns and construct vocabulary are unchanged.
Stored fragments may contain mutation identifiers and construct spans, but
may not contain `Ki`, `Kd`, `IC50`, `EC50`, `pKi`, `pKd`, potency, affinity,
inhibition percentages, or concentration units. No raw member is serialized
to a report or cache.

## 5. Candidate-level construct evidence

The exact X4 evidence forms are reused:

1. the named mutant was generated from the named WT construct/template;
2. all named mutants were generated from one named WT construct/template and
   the candidate mutation token appears in the same construct section/member;
3. WT and mutant share one explicit full-length state or residue span and
   expression construct.

Evidence is keyed by document, member, accession, mutation token, ChEMBL
target, and endpoint. One accepted mutation never licenses another mutation
automatically. X1 exact description pairs remain accepted.

## 6. Gates

| gate | condition |
| --- | --- |
| `X5_X4_BOUND` | X4 result hash, preregistration, firewall, and inadequate verdict bind |
| `X5_SOURCE_SET_FROZEN` | requests are exactly the nine X4-complete EPMC records |
| `X5_ARCHIVES_DISPOSITIONED` | every request has a valid archive or explicit transport/archive failure |
| `X5_MEMBER_POLICY_ENFORCED` | every parsed member has an allowed suffix and passes size/path/encryption/nesting rules |
| `X5_NO_OUTCOME_MATERIALIZED` | retained projections contain no forbidden outcome pattern |
| `X5_EXACT_CONSTRUCT_RULE_ENFORCED` | every accepted candidate/member has one frozen evidence form |
| `X5_DISCOVERY_ACCESSION_EXCLUDED` | `P15056` contributes zero primary components |
| `X5_FAMILY_MAPPED` | every primary component has a nonempty family disposition |
| `X5_SOURCE_TOPOLOGY_ADEQUATE` | at least 25 non-BRAF components, at least six broad families, and largest accession share at most 50% |

Verdicts:

- all gates pass: `OMUT_X5_SUPPLEMENT_CONSTRUCT_REGISTRY_ADEQUATE`;
- execution/firewall gates pass but topology fails:
  `OMUT_X5_SUPPLEMENT_CONSTRUCT_REGISTRY_INADEQUATE_STOP`;
- any binding, archive-disposition, cache, member-policy, or firewall gate
  fails: `OMUT_X5_INCOMPLETE_STOP`.

## 7. Deliverables

- `research/omut_x5.py`
- `tests/test_omut_x5.py`
- `reports/active/omut_x5.json`
- `reports/active/omut_x5_decision.md`
