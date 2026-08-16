# LEXOR L0B5 Exact-Query Transport Probe Preregistration

Date: 2026-07-27

## Scope and firewall

L0B3 and L0B4 established that the L0B2 parameter names, minimal projection,
and page size return HTTP 200 for a synthetic non-candidate query. This
independent diagnostic tests the exact failed L0B2 first-query shape. Its
single actual search phrase is explicitly registered in the plan through
`allow_exact_search=true`; successful response bodies are not read, parsed,
retained, or serialized. It is not discovery and cannot produce a candidate
inventory, L0 audit, or L1 authorization.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `candidate_metadata_accepted=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`

## Frozen protocol

Plan: `manifests/lexor_l0b5_transport_probe_plan.v1.json`

<!-- LEXOR_L0B_TRANSPORT_PROBE_PLAN_SHA256: ed4fcd154a3f314c30e6cac5bccc8fdf4cec9e63416ee1dfec98eebe43045a31 -->

The runner may issue exactly one unauthenticated HTTPS GET request. The result
may distinguish keyword-specific incompatibility from a transient failure only;
it cannot establish source eligibility, query depth, or any LEXOR gate.
