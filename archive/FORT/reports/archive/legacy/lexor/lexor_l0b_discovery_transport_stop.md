# LEXOR L0B Discovery Transport Stop

Date: 2026-07-27

## Outcome

`LEXOR_L0B_DISCOVERY_TRANSPORT_STOP`

One preregistered metadata request returned HTTP 400 before any candidate
metadata was accepted, reduced, or written. The runner stopped immediately, as
registered. Because it wrote no per-request progress artifact before stopping,
the failing endpoint and request position were not identified. No discovery
artifact, candidate inventory, L0 audit, LLM/API call, raw measurement read,
model training, or FORT-label access occurred.

## Interpretation

This is an endpoint-parameter compatibility failure, not evidence about the
literature corpus. The invalid plan is retained unchanged at
`manifests/lexor_l0b_discovery_plan.v1.json`; it will not be rerun. A separate
L0B1 preregistration is required before any corrected request is made.
