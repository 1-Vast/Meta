# UBSE-A0D Corrected Remote-Availability Preregistration

**Frozen:** 2026-07-30 before full URL enumeration or execution  
**Scope:** HTTP metadata only; no coordinate body or scientific label

## Question

Does the accepted A0C v2 manifest point to enough currently available
official RCSB mmCIF files to justify an A1-v2 source/topology
preregistration?

This gate tests source availability only. It does not parse coordinates,
extract events, validate ligand instances inside coordinates, establish
teacher reliability, or authorize affinity.

## Frozen input

`dataset/public/biolip2/processed/ubse_a0c_3d_event_sources_v2.parquet`

SHA-256:

`adc72f142e515c47ea18d20d7af08f6a434a30202a1483dcb362115062a068d5`

Expected unique PDB URLs:

| Role | Unique PDB URLs |
|---|---:|
| fit | 2,496 |
| validation | 140 |
| audit | 197 |
| total | 2,833 |

Cross-role PDB/URL overlap must be zero.

## URL and transport contract

Every URL must match exactly:

```text
https://files.rcsb.org/download/{lowercase_pdb_id}.cif.gz
```

No alternative host, query, fragment, credential, or user-provided redirect
target is accepted. Redirects may be followed only when the final host is
still `files.rcsb.org`.

For each unique URL:

1. send HTTP `HEAD` with a fixed FORT user agent;
2. use at most 16 concurrent connections;
3. use a 20-second per-attempt timeout;
4. retry connection/timeout errors, HTTP 429, and HTTP 5xx at most twice;
5. do not retry stable HTTP 4xx other than 429;
6. never issue a GET or read a response body.

A URL is available only when the final response is HTTP 200 and its content
type is `application/gzip`, `application/x-gzip`, or
`application/octet-stream`.

The ledger records role, PDB ID, requested/final URL, attempts, status,
content type, availability, and a normalized error reason. It contains no
coordinate bytes.

## Frozen outputs

- `reports/active/ubse_a0d_remote_availability_ledger.parquet`
- `reports/active/ubse_a0d_remote_availability.json`

Neither output may overwrite an existing artifact.

## Gates

All identity and closure gates are binding:

1. **A0D-1 manifest identity:** the input hash matches.
2. **A0D-2 URL topology:** exactly 2,833 unique PDB URLs with exact role
   counts 2,496/140/197 and zero cross-role overlap.
3. **A0D-3 origin firewall:** every requested and successful final URL obeys
   the official RCSB HTTPS origin/path contract.
4. **A0D-4 overall availability:** at least 95% of all 2,833 URLs are
   available.
5. **A0D-5 role availability:** fit, validation, and audit availability are
   each at least 95%.
6. **A0D-6 response firewall:** all requests are HEAD-only; zero response
   body bytes, coordinate parses, event labels, binding-residue labels,
   affinity fields/values, development outcomes, confirmation outcomes, or
   sealed outcomes are loaded.
7. **A0D-7 ledger integrity:** one row per frozen URL, no duplicate PDB/URL,
   all response/status fields finite or explicitly missing, and the ledger
   hash is recorded.

## Decisions

Pass:

`REQUEST_UBSE_A1V2_SOURCE_AND_TOPOLOGY_PREREGISTRATION`

If the gate cannot execute because every representative URL fails with a
local socket, DNS, TLS, or policy error:

`WAIT_UBSE_A0D_NETWORK_RUNTIME`

If execution is possible but real HTTP/content-type availability fails the
overall or any role threshold:

`STOP_UBSE_A0D_REMOTE_COORDINATE_COVERAGE_INADEQUATE`

No decision unlocks coordinate parsing until the corrected A1-v2
identifiability, fresh-confirmation, membership, extractor, and functional-
group checkerboard contracts are frozen.

## Preliminary environment probe

Before this full preregistration, three representative official URLs
(`5yic`, `7e2u`, and `7sly`) returned HTTP 200 with
`application/gzip` under an approved read-only network check. Those three
responses are environment evidence only and are not counted toward the
2,833-URL gate.
