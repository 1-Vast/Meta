# PANEL-EVIDENCE P0 local decision

## Verdict

`PANEL_EVIDENCE_P0_LOCAL_EVIDENCE_HARNESS_PASS__LLM_UNTESTED`

This no-API audit validates only deterministic table evidence linkage on a known local source. It
does not test LLM semantic extraction, discover a new document, establish an independent
provenance family, or authorize any model training.

## Frozen source and gold

* Extracted matrix: `41589_2023_1459_MOESM3_ESM.xlsx`, sheet `Kinobeads Drugmatrix - all`.
* Layout gold: `41589_2023_1459_MOESM4_ESM.xlsx`, sheet `CATDS target`. It is a differently formatted table in the same
  paper and therefore checks extraction/layout fidelity, not biological replication.
* Compound evidence: `41589_2023_1459_MOESM2_ESM.xlsx` structure tables, with RDKit canonical parent connectivity.
* Target candidates: frozen local reviewed-UniProt cache; no network or API call was made.

## Results

| measure | value |
| --- | ---: |
| accepted evidence-bound cells | 9346 |
| exact tuple precision vs long-table gold | 1.000000 |
| exact tuple recall vs raw-source-eligible long-table gold | 1.000000 |
| raw-value mismatches | 0 |
| long-table positive values rejected because source cell was nonpositive | 20 |
| source-coordinate roundtrip errors | 0 |
| blank/unreported cells abstained | 365641 |
| nonpositive cells rejected | 22 |

## Admission boundary

All emitted records carry source file, sheet, row/column coordinates, raw value/unit/endpoint,
compound evidence, target-candidate evidence, and one declared provenance family. They remain blocked
from F0 because this source is already known, has one campaign/provenance family, and does not expose
the construct/independent-replication facts required by the strict dual-cold contract.

The next admissible step is a blind, API-free fixture suite containing ambiguous headers, missing
units, censoring, and unresolved entities. API use remains unauthorized until a separately
preregistered semantic-extraction test passes those checks.
