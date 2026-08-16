# UBSE-A1-v2 HEAD-only remote-availability preregistration

Date: 2026-07-30  
Status: frozen before any A1-v2 network request  
Scope: HTTP metadata only; no coordinate body or scientific label

## 1. Question and boundary

Are enough of the 997 frozen official RCSB mmCIF URLs available to retain at
least 128 complete A1-R reliability units and exactly 512 A1-C primary
candidates after the already frozen reserve rule?

This audit may record only HEAD response metadata. A pass does not download or
parse coordinates, establish event coverage, complete `SR-5`, pass `SR-6`,
admit P0A, or authorize affinity, confirmation scoring, Stage-2, or sealed
access.

## 2. Frozen inputs and implementation

| Input | SHA-256 |
|---|---|
| selected 997-URL manifest | `284fdb770b5539476facc8a3e873dac676d1df5f836060735eac6fb1979adbd7` |
| A1-R locator manifest | `da788694eaf0f072752a8b7dfbaac6c9ddcd9bf57fcf9709364a1dc7902f0083` |
| A1-C locator manifest | `629bcd2e4ad1164b4a7d0dd330687735b5da94c15c3a0d650c475b0887b96491` |
| locator result | `cb49ebf701234c1ef41133ea6d9238fc47dabdf2437d0503e60c98d1f990f2e2` |

Frozen implementation:

- `research/ubse_a1v2_remote_availability.py`:
  `0eebfee49528964c2b52d51a2ac5fdfe444f7a904247206eaa112d497ee0a936`;
- `tests/test_ubse_a1v2_remote_availability.py`:
  `d3fdf2b8e90c7bba3147f0313129ecd4d5dd225e474da1b123cb6d8689703eaa`;
- synthetic verification before freeze: `5 passed`.

Any implementation or threshold change requires an amendment before another
network execution.

## 3. Frozen URL topology

| Role | Unique PDB/URL rows |
|---|---:|
| `a1r` | 421 |
| `a1c_primary` | 512 |
| `a1c_reserve` | 64 |
| Total | 997 |

Every URL must equal:

```text
https://files.rcsb.org/download/{lowercase_pdb_id}.cif.gz
```

PDB identities are globally unique across all three roles.

## 4. Transport firewall

For every URL:

1. issue HTTP `HEAD` only;
2. use user agent `FORT-UBSE-A1V2-H0/1.0`;
3. use at most 16 concurrent connections;
4. use a 20-second per-attempt timeout;
5. retry connection/timeout errors, HTTP 429, and HTTP 5xx at most twice;
6. do not retry stable non-429 4xx responses;
7. follow a redirect only if the final URL still satisfies the exact official
   RCSB origin and path contract;
8. never issue `GET` or read a response body.

A URL is available only for HTTP 200, a valid final origin, and content type
`application/gzip`, `application/x-gzip`, or
`application/octet-stream`.

The ledger records role, PDB ID, requested/final URL, request method,
attempts, response presence/status/content type, origin validity,
availability, and normalized failure reason. `response_body_bytes` must be
zero in every row.

## 5. Frozen reserve rule

Available frozen A1-C primaries retain their original `role_rank`. For each
unavailable primary in increasing `role_rank`, promote the earliest available
reserve in frozen reserve order. Locator or HEAD failure is the only legal
promotion reason. Event coverage, checkerboard availability, chemistry,
family, model output, or performance can never trigger replacement.

The effective primary manifest records the frozen role/rank, effective
primary rank, replaced primary rank, and selection reason.

## 6. Outputs

- `reports/active/ubse_a1v2_remote_availability_ledger.parquet`
- `dataset/public/biolip2/processed/ubse_a1v2_a1c_head_complete_primary.parquet`
- `reports/active/ubse_a1v2_remote_availability.json`

All outputs are create-once. A complete zero-response network-runtime failure
must raise `WAIT_UBSE_A1V2_HEAD_NETWORK_RUNTIME` before creating outputs so
the identical frozen execution can be retried in a network-enabled runtime.

## 7. Gates

1. `H0-1 Frozen inputs`: all input hashes match.
2. `H0-2 URL topology`: exact role counts 421/512/64 and 997 globally unique
   PDB/URL identities.
3. `H0-3 Origin firewall`: every requested and available final URL obeys the
   official RCSB contract.
4. `H0-4 HEAD response firewall`: all methods are HEAD and every
   `response_body_bytes` value is zero.
5. `H0-5 A1-R availability`: at least 128 targets have all three frozen
   locator instances available.
6. `H0-6 A1-C primary after reserves`: exactly 512 unique effective primary
   targets remain under the frozen replacement order.
7. `H0-7 Ledger integrity`: all 997 frozen request identities occur exactly
   once and attempts stay within the frozen maximum.

Pass:

```text
FREEZE_A1V2_HEAD_AVAILABILITY_KEEP_COORDINATE_BODIES_LOCKED
```

Complete zero-response runtime failure:

```text
WAIT_UBSE_A1V2_HEAD_NETWORK_RUNTIME
```

Source/topology failure:

```text
STOP_UBSE_A1V2_REMOTE_AVAILABILITY_INADEQUATE
```
