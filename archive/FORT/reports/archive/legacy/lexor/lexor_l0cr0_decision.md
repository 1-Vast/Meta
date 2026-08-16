# LEXOR L0C-R0 Calibration Decision

## Verdict

`L0CR0_MECHANICAL_OBSERVABILITY_PASS__SOURCE_NONCOUNTABLE`

This is a local calibration on an already development-seen source. It is not a
new discovery source, a blind fixture, an L0 pass, or an L1 authorization.

## Mechanical Observability

| measure | value |
| --- | ---: |
| accepted evidence-bound cells | 9346 |
| source-coordinate roundtrip errors | 0 |
| first ledger SHA-256 | `3e0e96ece56dd5f715f1f88d7a151ff51d0121194279129f3d85f6c2250e712a` |
| second ledger SHA-256 | `3e0e96ece56dd5f715f1f88d7a151ff51d0121194279129f3d85f6c2250e712a` |
| deterministic replay | True |

All three local source hashes, registered table locators, required headers, and
coordinate roundtrips passed before the ledger was emitted.

## Source Eligibility

`countable_for_l0 = False`. Every emitted record
has `construct_status = not_reported_in_local_table`; the source table also
supplies only a gene label, while the frozen local cache supplies an accession
candidate rather than source-verified accession evidence. Query depth after all
firewalls is therefore intentionally `null`. The source also represents exactly
one provenance family and cannot satisfy L0's independent-environment portfolio.

## Boundaries

* historical exact supplementary-file URLs remain incomplete, so the manifest
  is calibration-only and cannot certify a new L0C acquisition;
* the local reviewed-UniProt cache supplies an accession candidate only and
  does not repair the missing construct evidence;
* no external network, LLM/API call, model training, confirmation-label read,
  or sealed-test read occurred;
* no L0, L1, API, or training authorization was granted.
