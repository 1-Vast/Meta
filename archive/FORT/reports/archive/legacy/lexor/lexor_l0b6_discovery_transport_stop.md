# LEXOR L0B6 Discovery Transport Stop

Date: 2026-07-27

## Outcome

`LEXOR_L0B6_DISCOVERY_TRANSPORT_STOP`

One preregistered metadata request returned HTTP 400 before any candidate
metadata was accepted, reduced, or written. The runner stopped immediately.
Because it wrote no per-request progress artifact before stopping, the failing
endpoint and request position were not identified. No discovery artifact,
candidate inventory, L0 audit, LLM/API call, raw measurement read, model
training, or FORT-label access occurred.

## Interpretation

The exact L0B2 OpenAlex query shape returned HTTP 200 in L0B5 only under the
transport probe, but L0B6's uninstrumented 400 cannot be attributed to that
OpenAlex request: it may instead have occurred at any later request, including
Zenodo. L0B6 remains frozen and will not be rerun; a protocol-only,
per-endpoint diagnostic is required before any later discovery plan can be
justified.

## Firewall State

* raw measurement files read: `False`
* external metadata accepted: `False`
* LLM API called: `False`
* model trained: `False`
* FORT labels read: `False`
* sealed test consumed: `False`
