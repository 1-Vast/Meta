# LEXOR L0B2 Discovery Transport Stop

Date: 2026-07-27

## Outcome

`LEXOR_L0B2_DISCOVERY_TRANSPORT_STOP`

One preregistered metadata request returned HTTP 400 before any candidate
metadata was accepted, reduced, or written. The runner stopped immediately.
Because it wrote no per-request progress artifact before stopping, the failing
endpoint and request position were not identified. No discovery artifact,
candidate inventory, L0 audit, LLM/API call, raw measurement read, model
training, or FORT-label access occurred.

## Interpretation

L0B2 removed all nested and optional OpenAlex projection fields, retaining only
`id,doi,title`, but the uninstrumented stop does not establish that the HTTP
400 came from OpenAlex. It therefore cannot reject or confirm an explanation
about OpenAlex projection fields. The runner also did not retain the endpoint
error body, so it identifies neither the failing source nor its incompatible
parameter.

The frozen L0B2 plan remains unchanged at
`manifests/lexor_l0b2_discovery_plan.v1.json` and will not be rerun. Any
subsequent transport diagnosis must be independently preregistered, avoid
candidate acceptance, and retain only a bounded protocol-error description.

## Firewall State

* raw measurement files read: `False`
* external metadata accepted: `False`
* LLM API called: `False`
* model trained: `False`
* FORT labels read: `False`
* sealed test consumed: `False`
