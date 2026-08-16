# LEXOR L0B6 Discovery Preregistration: Post-Probe Fresh Discovery Run

Date: 2026-07-27

## Scope and firewall

L0B1 and L0B2 are retained as zero-data HTTP 400 stops and will not be rerun.
L0B3-L0B5 then established, without accepting candidate metadata, that the
OpenAlex parameter shape and the exact L0B2 first query return HTTP 200. L0B6
is therefore a fresh, independent, no-retry discovery run, not an amendment to
any frozen failed plan.

It may read only the seven registered public bibliographic/repository metadata
responses. It may not read article full text, supplementary files,
measurements, structures, FORT labels, `.env`, or call a model API.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`

## Frozen discovery protocol

Plan: `manifests/lexor_l0b6_discovery_plan.v1.json`

<!-- LEXOR_L0B_DISCOVERY_PLAN_SHA256: 04c59f86a55bde1c0a7fbd65a4492240471c20256e608d9093c9b0c4466b2fef -->

The runner may make exactly seven unauthenticated HTTPS GET requests under this
plan. It records only the approved metadata fields and may not infer a
scaffold-diverse query count. Any transport error stops L0B6 without a partial
inventory.

## Transition

If discovery completes, a separate L0B6 audit preregistration must freeze the
candidate-inventory hash before a local L0 audit runs. No discovery artifact
may be added after that audit starts.
