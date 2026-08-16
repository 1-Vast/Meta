# LEXOR L0B4 Page-Size Transport Probe Preregistration

Date: 2026-07-27

## Scope and firewall

L0B3 established that the L0B2 parameter names and minimal projection are
compatible at `per-page=1`. This independent diagnostic tests the sole
uncovered difference in the failed discovery request: `per-page=200`. It is
not discovery and cannot produce a candidate inventory, L0 audit, or L1
authorization. Successful bodies are not read, parsed, retained, or
serialized; an HTTP error may contribute only a bounded protocol description.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `candidate_metadata_accepted=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`

## Frozen protocol

Plan: `manifests/lexor_l0b4_transport_probe_plan.v1.json`

<!-- LEXOR_L0B_TRANSPORT_PROBE_PLAN_SHA256: 5fc19f52db0683409a3dce2a0d2381c7fd73485f1c10b6cd08bcc61b7b3737d9 -->

The runner may issue exactly two unauthenticated HTTPS GET requests. The first
tests page size alone; the second uses the complete L0B2 OpenAlex parameter
shape with that same page size. The result may identify a compatible page size
only and may not establish source eligibility or query depth.
