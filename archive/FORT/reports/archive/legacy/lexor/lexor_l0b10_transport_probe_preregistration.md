# LEXOR L0B10 Zenodo Anonymous Page-Size Probe Preregistration

Date: 2026-07-27

## Scope and firewall

L0B9 established from Zenodo's protocol response that anonymous requests may
not use `size=100` and must use at most 25. This independent diagnostic
confirms `size=25` with a synthetic non-candidate query. Successful bodies are
not read, parsed, retained, or serialized. It is not discovery and cannot
produce a candidate inventory, L0 audit, or L1 authorization.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `candidate_metadata_accepted=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`

## Frozen protocol

Plan: `manifests/lexor_l0b10_transport_probe_plan.v1.json`

<!-- LEXOR_L0B_TRANSPORT_PROBE_PLAN_SHA256: 0233896dbffe887ba2cd908865c48a3b9e79603c1530c4eeca0c710035c283aa -->

The runner may issue exactly one unauthenticated HTTPS GET request. The result
may confirm the anonymous page-size limit only; it cannot establish source
eligibility, query depth, or any LEXOR gate.
