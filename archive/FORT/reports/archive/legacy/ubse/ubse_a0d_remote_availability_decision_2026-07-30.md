# UBSE-A0D Corrected Remote-Availability Decision

**Date:** 2026-07-30  
**Decision:** `REQUEST_UBSE_A1V2_SOURCE_AND_TOPOLOGY_PREREGISTRATION`

## Outcome

All 2,833 unique official RCSB mmCIF URLs in the accepted A0C v2 manifest
returned HTTP 200 with `application/gzip`. Overall, fit, validation, and audit
coverage are each 100%, so the corrected A0 remote-availability gate passes.

This removes the network-addressability WAIT only. It does not authorize
coordinate parsing, validate the ligand instance inside a coordinate file,
extract typed events, establish a functional-group checkerboard, or unlock
affinity.

## Frozen input

Manifest:

`dataset/public/biolip2/processed/ubse_a0c_3d_event_sources_v2.parquet`

SHA-256:

`adc72f142e515c47ea18d20d7af08f6a434a30202a1483dcb362115062a068d5`

The input contains 3,467 complex rows over 2,833 role-disjoint PDB entries.

## Availability

| Role | Checked | Available | Coverage |
|---|---:|---:|---:|
| fit | 2,496 | 2,496 | 100% |
| validation | 140 | 140 | 100% |
| audit | 197 | 197 | 100% |
| total | 2,833 | 2,833 | 100% |

Every response was:

- request method `HEAD`;
- status 200;
- content type `application/gzip`;
- final origin `https://files.rcsb.org`;
- zero response-body bytes.

Attempt histogram:

- first attempt: 2,832;
- second attempt: 1;
- third attempt: 0.

There were no unavailable URLs or failure reasons. Execution wall time was
76.686 seconds with at most 16 concurrent connections.

## Ledger integrity

| Check | Result |
|---|---:|
| rows | 2,833 |
| unique PDB IDs | 2,833 |
| unique requested URLs | 2,833 |
| duplicate PDB/URL rows | 0 |
| cross-role PDB overlap | 0 |
| non-HEAD methods | 0 |
| response-body bytes | 0 |
| successful final-origin violations | 0 |

## Gates

All seven frozen gates pass:

- A0D-1 manifest identity;
- A0D-2 exact URL/role topology;
- A0D-3 RCSB origin firewall;
- A0D-4 overall availability;
- A0D-5 per-role availability;
- A0D-6 HEAD-only response firewall;
- A0D-7 ledger integrity.

## Firewall

The audit loaded no:

- coordinate bytes;
- parsed coordinates;
- event labels;
- binding-residue labels;
- affinity fields or values;
- development/confirmation outcomes;
- sealed outcomes.

## Next boundary

A1-v1 remains blocked by the coupling-identifiability correction. The next
permitted action is an A1-v2 source/topology preregistration that freezes:

- a fresh untouched confirmation role;
- ligand-neighbour, pocket, template, and model-membership closure;
- a precise coordinate-instance and extractor contract;
- a true residue-functional-group checkerboard power gate;
- pair-conditioned rank-one, fixed-margin/fixed-mass, and dustbin nulls.

Only after that preregistration may coordinate files be downloaded and
parsed.

## Authoritative artifacts

- `reports/active/ubse_a0d_remote_availability_preregistration_2026-07-30.md`
- `reports/active/ubse_a0d_remote_availability.json`
- `reports/active/ubse_a0d_remote_availability_ledger.parquet`
- `research/ubse_a0d_remote_availability.py`
- `tests/test_ubse_a0d_remote_availability.py`

Result SHA-256:

`e7cab03a2ae0457e16101dd505a60fa3d88bf47c625d20b082b18708dc640667`

Ledger SHA-256:

`05d15da5ea73301783bf86cd424bc2542915704383b24b73465d24c8ce50b872`
