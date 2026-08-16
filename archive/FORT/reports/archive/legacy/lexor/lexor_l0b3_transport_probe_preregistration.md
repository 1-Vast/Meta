# LEXOR L0B3 Transport Probe Preregistration

Date: 2026-07-27

## Scope and firewall

This is a protocol-compatibility diagnostic following the frozen zero-data
L0B2 stop. It is not literature discovery and cannot produce a candidate
inventory, L0 audit, or L1 authorization. Each request carries the same
synthetic non-candidate search token. Successful response bodies are not read,
parsed, retained, or serialized; HTTP error responses may contribute only a
bounded parameter-error description.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `candidate_metadata_accepted=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`

## Frozen protocol

Plan: `manifests/lexor_l0b3_transport_probe_plan.v1.json`

<!-- LEXOR_L0B_TRANSPORT_PROBE_PLAN_SHA256: a4d2717b0525cc80bec10ea3ba764ad827b284b6bf19b355efe90496e5d652e7 -->

The runner may issue exactly five unauthenticated HTTPS GET requests. It adds
the OpenAlex parameters incrementally: base search, `per-page`, filter, sort,
and minimal `select`. The result may identify a compatible query shape only;
it must not make a claim about source eligibility or query depth.
