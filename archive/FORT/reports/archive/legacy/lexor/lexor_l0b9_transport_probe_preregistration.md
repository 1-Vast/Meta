# LEXOR L0B9 Zenodo Parameter-Isolation Probe Preregistration

Date: 2026-07-27

## Scope and firewall

L0B8 localized the prior uninstrumented discovery failure to Zenodo: all three
registered Zenodo request shapes returned HTTP 400. This independent probe
uses a synthetic non-candidate query and adds the three Zenodo parameters
incrementally (`q`, `size`, `sort=mostrecent`). Successful bodies are not read,
parsed, retained, or serialized. It is not discovery and cannot produce a
candidate inventory, L0 audit, or L1 authorization.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `candidate_metadata_accepted=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`

## Frozen protocol

Plan: `manifests/lexor_l0b9_transport_probe_plan.v1.json`

<!-- LEXOR_L0B_TRANSPORT_PROBE_PLAN_SHA256: 80d3773c36197d467c0debd35ccbb54f64e8a3d8f19384e744a9b12718418815 -->

The runner may issue exactly three unauthenticated HTTPS GET requests. The
result may identify one incompatible Zenodo parameter only; it cannot establish
source eligibility, query depth, or any LEXOR gate.
