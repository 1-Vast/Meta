# LEXOR L0B8 Zenodo Transport Probe Preregistration

Date: 2026-07-27

## Scope and firewall

The attribution correction establishes that prior discovery HTTP 400 stops
cannot be assigned to OpenAlex. This independent protocol-only probe tests all
three exact Zenodo request shapes from the discovery plan under the same
metadata-only User-Agent. Successful response bodies are not read, parsed,
retained, or serialized. It is not discovery and cannot produce a candidate
inventory, L0 audit, or L1 authorization.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `candidate_metadata_accepted=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`

## Frozen protocol

Plan: `manifests/lexor_l0b8_transport_probe_plan.v1.json`

<!-- LEXOR_L0B_TRANSPORT_PROBE_PLAN_SHA256: 15c4b51d5d29a16a931b47682a947d5ae72ab8ab7f5af8d6337e0c25f1512568 -->

The runner may issue exactly three unauthenticated HTTPS GET requests. The
result may identify a Zenodo protocol incompatibility only; it cannot establish
source eligibility, query depth, or any LEXOR gate.
