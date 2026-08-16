# UBSE-A1-v2 strict HEAD-firewall verification preregistration

Date: 2026-07-30  
Status: frozen before H0V network execution  
Scope: transport correction only; coordinate bodies and labels remain locked

## 1. Reason for verification

The original H0 execution returned 997/997 HTTP 200
`application/gzip` responses on the first attempt, and every final URL equaled
its requested URL. No redirect was observed. However, its client permitted
redirect following and its `response_body_bytes=0` field was a protocol
constant rather than the transport library's downloaded-byte counter.

H0V preserves the original artifacts and repeats the exact URL set with a
stricter transport:

- `follow_redirects=False`;
- streaming HEAD;
- response stream never iterated;
- `response.num_bytes_downloaded` recorded for every request.

## 2. Frozen inputs

| Input | SHA-256 |
|---|---|
| selected URL manifest | `284fdb770b5539476facc8a3e873dac676d1df5f836060735eac6fb1979adbd7` |
| initial H0 ledger | `fa288a0f60374c0f4e81e5fedc38e2c9215675e515eea8acbb6e6336759cfea3` |
| initial H0 result | `64448c919073aba6c06a6240785d97bae211d0edd768fa2cadc4eeb33aae9c7b` |

Frozen implementation:

- `research/ubse_a1v2_head_firewall_verification.py`:
  `272293c9173e3646679c006cc85f4808044be834cabc22a9468557f57bd756ea`;
- `tests/test_ubse_a1v2_head_firewall_verification.py`:
  `36ee8b0f90cffb0851b9ca7acd9e0b30b7211d56d303020994616b66005d4b34`;
- synthetic verification: `4 passed`.

## 3. Frozen transport

- exact 997 official RCSB URLs;
- HTTP method `HEAD`;
- `follow_redirects=False`;
- streaming response with no body iteration;
- maximum 16 concurrent connections;
- 20-second timeout;
- exactly two permitted retries for connection/timeout, 429, or 5xx;
- user agent `FORT-UBSE-A1V2-H0V/1.0`;
- accepted content types remain `application/gzip`,
  `application/x-gzip`, and `application/octet-stream`.

A redirect response is rejected without a second request. A row is available
only when status is 200, content type and exact URL origin/path pass, final URL
equals requested URL, and actual downloaded bytes equal zero.

## 4. Outputs and decisions

- `reports/active/ubse_a1v2_head_firewall_verification_ledger.parquet`
- `reports/active/ubse_a1v2_head_firewall_verification.json`

Both are write-once.

Pass requires all 997 identities to reproduce availability with zero redirect
following, zero body-stream iteration, and zero actual downloaded bytes:

```text
FREEZE_A1V2_H0V_STRICT_HEAD_FIREWALL_KEEP_BODIES_LOCKED
```

Complete absence of HTTP responses gives:

```text
WAIT_UBSE_A1V2_H0V_NETWORK_RUNTIME
```

Any transient availability divergence is recorded as WAIT rather than used
to alter roles:

```text
WAIT_UBSE_A1V2_H0V_AVAILABILITY_DIVERGENCE
```

No H0V outcome can unlock a coordinate GET, extractor, event, affinity,
confirmation score, Stage-2 model, or sealed test.
