# LEXOR L0B Discovery Preregistration: Public Metadata Expansion

Date: 2026-07-27

## Scope and firewall

This is a new discovery phase after the local-only L0 stop. It is not an
amendment to, or continuation of, the frozen local inventory. It may query only
the public metadata endpoints and fields frozen in
`manifests/lexor_l0b_discovery_plan.v1.json`. It must not read an article,
supplement, raw measurement table, affinity value, chemical structure file,
FORT development/confirmation/sealed label, `.env` credential, or model API.

Public metadata HTTP requests are allowed only for the bounded plan below. They
are not LLM/API extraction calls. The later L0B audit itself must be local and
must not perform a network request.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`

## Frozen discovery protocol

Plan: `manifests/lexor_l0b_discovery_plan.v1.json`

<!-- LEXOR_L0B_DISCOVERY_PLAN_SHA256: c6942a50995be3d49d11271a099f00fc8e94592e469024a63e312d76a2a12aa2 -->

The runner may make at most seven unauthenticated HTTPS GET requests, one for
each predeclared query. It records request URLs without credentials, response
status, retrieval time, and SHA-256 of the reduced metadata artifact. It stores
only the approved metadata fields. Responses are deduplicated by DOI and then
source identifier.

No candidate is made eligible from a title, aggregate record count, or inferred
chemical diversity. An explicit post-firewall scaffold-diverse query-ligand
count must already be present in the allowed metadata fields. Missing remains
missing.

## Transition to L0B

After discovery, the resulting metadata artifact is frozen. A separate
`lexor_l0b_preregistration.md` must bind its candidate inventory SHA-256 before
the local `research/lexor_l0.py` audit runs. The L0B audit uses the unchanged
provenance-family, open-license, and MDE80 gates. A failing L0B does not permit
expanding queries after seeing the result; another expansion requires a new
discovery plan and preregistration.
