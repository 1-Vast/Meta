# OpenMut `OMUT-X4` preregistration

**Frozen before execution:** 2026-07-28.
**Predecessor:** `OMUT-X3C`,
`OMUT_X3C_LICENSED_LINK_RECOVERY_FEASIBLE`.
**Stage boundary:** source transport and outcome-free construct-statement
recovery.

## 1. Question

X3C found an optimistic 32-component topology from 11 licensed/open document
records. X4 asks:

> Do those frozen source bodies explicitly establish a common WT construct
> for enough candidate WT/mutant assay pairs to retain at least 25
> non-BRAF components across six broad families, with no accession supplying
> more than 50%?

X4 does not add affinity outcomes. It may read source bodies only to project
construct provenance and assay-context statements.

## 2. Frozen sources

The source set is exactly the 11 X3C-accessible documents:

- the nine Europe PMC `fullTextXML` records already accepted by X2;
- the two exact Crossref version-of-record XML links accepted by X3C.

No search result, manually substituted URL, abstract endpoint, DOI landing
page, shadow library, social-sharing site, or newly discovered document may
enter X4.

Each request records requested URL, final URL, HTTP status, media type,
byte-length, response SHA-256, and a canonical projection SHA-256. Raw source
bodies are not retained after projection.

## 3. Transport

A source body is transport-complete only when:

- the response is HTTP 200 after HTTPS redirects;
- the final host is Europe PMC or the exact host accepted by X3C;
- the payload is nonempty XML or HTML with no authentication, rate-limit, or
  error envelope;
- the response parses as a document.

Authorization failures and placeholder/error documents are explicit failed
dispositions, not inaccessible-evidence negatives.

## 4. Outcome firewall and projection

Only section paths and text fragments relevant to construct identity may be
retained. A fragment must be in Methods, Experimental, Materials,
Supplementary Methods, a construct table/caption, or an equivalent
source-preparation section and contain:

- a candidate mutation token or a statement explicitly applying to all
  mutants; and
- at least one construct term: `construct`, `cDNA`, `clone`, `plasmid`,
  `vector`, `full-length`, an explicit residue span, `site-directed
  mutagenesis`, `wild type`, `WT`, `template`, `expressed`, or `transfected`.

Fragments containing `Ki`, `Kd`, `IC50`, `EC50`, `pKi`, `pKd`, potency,
affinity, inhibition percentages, or concentration units are rejected from
the projection. The stored artifact may contain no numeric affinity value.

## 5. Candidate-level construct evidence

A near-exact candidate/document pair is accepted only if the projected source
states one of:

1. the named mutant was generated from the named WT construct/template;
2. all named mutants were generated from one named WT construct/template and
   the candidate mutation token appears in the same construct section;
3. the WT and mutant are both assigned the same explicit full-length state or
   residue span and expression construct.

Target identity is bound by exact ChEMBL document DOI plus accession/target
metadata. A generic statement that mutations were made, canonical UniProt
identity, an unspecified commercial protein, or a mutation token outside a
construct section is insufficient. Evidence is candidate-specific; one
accepted mutation does not automatically license every mutation in the same
paper.

Existing X1 exact description pairs remain accepted without source-body
recovery. X2 near pairs enter only through an accepted candidate/document
disposition.

## 6. Topology and gates

| gate | condition |
| --- | --- |
| `X4_X3C_BOUND` | X3C result, preregistration, firewall, and feasible verdict bind |
| `X4_SOURCE_SET_FROZEN` | requests are exactly the 11 X3C-accessible records |
| `X4_TRANSPORT_DISPOSITIONED` | every frozen source has a complete or explicit failed transport disposition |
| `X4_NO_OUTCOME_MATERIALIZED` | projections contain no forbidden outcome field or numeric affinity pattern |
| `X4_EXACT_CONSTRUCT_RULE_ENFORCED` | every accepted near pair carries one of the three frozen candidate-level evidence forms |
| `X4_DISCOVERY_ACCESSION_EXCLUDED` | `P15056` contributes zero primary components |
| `X4_FAMILY_MAPPED` | every primary component has a nonempty family disposition |
| `X4_SOURCE_TOPOLOGY_ADEQUATE` | at least 25 non-BRAF components, at least six broad families, and largest accession share at most 50% |

Verdicts:

- all gates pass: `OMUT_X4_CONSTRUCT_REGISTRY_ADEQUATE`;
- execution/firewall gates pass but topology fails:
  `OMUT_X4_CONSTRUCT_REGISTRY_INADEQUATE_STOP`;
- any binding, transport-disposition, cache, or firewall gate fails:
  `OMUT_X4_INCOMPLETE_STOP`.

## 7. Deliverables

- `research/omut_x4.py`
- `tests/test_omut_x4.py`
- `reports/active/omut_x4.json`
- `reports/active/omut_x4_decision.md`
