# LEXOR L0B1 Discovery Transport Stop

Date: 2026-07-27

## Outcome

`LEXOR_L0B1_DISCOVERY_TRANSPORT_STOP`

One preregistered metadata request returned HTTP 400 before any candidate
metadata was accepted, reduced, or written. The runner stopped immediately.
Because it wrote no per-request progress artifact before stopping, the failing
endpoint and request position were not identified. No discovery artifact,
candidate inventory, L0 audit, LLM/API call, raw measurement read, model
training, or FORT-label access occurred.

## Interpretation

This is a second endpoint-parameter compatibility failure, not evidence about
the literature corpus. The frozen L0B1 plan remains unchanged at
`manifests/lexor_l0b1_discovery_plan.v1.json` and will not be rerun. Any
subsequent discovery attempt must use a separately preregistered plan after the
transport fault has been identified without accepting candidate metadata.

## Firewall State

* raw measurement files read: `False`
* external metadata accepted: `False`
* LLM API called: `False`
* model trained: `False`
* FORT labels read: `False`
* sealed test consumed: `False`
