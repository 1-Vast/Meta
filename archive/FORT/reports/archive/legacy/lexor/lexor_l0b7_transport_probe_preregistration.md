# LEXOR L0B7 Header Transport Probe Preregistration

Date: 2026-07-27

## Scope and firewall

L0B5 returned HTTP 200 for the exact failed discovery URL under the
transport-probe User-Agent, while L0B6 returned HTTP 400 under the discovery
User-Agent. This independent diagnostic holds the URL fixed and uses only the
registered discovery User-Agent. Successful response bodies are not read,
parsed, retained, or serialized. It is not discovery and cannot produce a
candidate inventory, L0 audit, or L1 authorization.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `candidate_metadata_accepted=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`

## Frozen protocol

Plan: `manifests/lexor_l0b7_transport_probe_plan.v1.json`

<!-- LEXOR_L0B_TRANSPORT_PROBE_PLAN_SHA256: b39c10e780ccf36c1acdd2efe31d1eb46b20afa05487195f347e6dad5e030237 -->

The runner may issue exactly one unauthenticated HTTPS GET request. The result
may distinguish a header-specific incompatibility from a response-body handling
or endpoint-state issue only; it cannot establish source eligibility, query
depth, or any LEXOR gate.
