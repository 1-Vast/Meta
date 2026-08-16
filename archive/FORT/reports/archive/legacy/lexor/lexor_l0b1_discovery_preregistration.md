# LEXOR L0B1 Discovery Preregistration: Corrected Public Metadata Transport

Date: 2026-07-27

## Scope and firewall

This is an independent replacement for the zero-data L0B transport stop. It
uses the same bounded metadata discovery question and all the same no-full-text,
no-measurement, no-`.env`, and no-model-API restrictions. It changes only the
OpenAlex `select` parameter, removing the two unconfirmed optional fields that
caused the HTTP 400 before any metadata was observed.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`

## Frozen discovery protocol

Plan: `manifests/lexor_l0b1_discovery_plan.v1.json`

<!-- LEXOR_L0B_DISCOVERY_PLAN_SHA256: 086b8ef392478bf9687524a34761787789340cfed8d6e7b56e435dba937d9bee -->

The runner may make exactly seven unauthenticated HTTPS GET requests under this
plan. It records only the approved metadata fields and may not infer a
scaffold-diverse query count. Any transport error stops this L0B1 run without a
partial inventory.

## Transition to L0B1

If discovery completes, a separate L0B1 audit preregistration will freeze the
candidate inventory hash before the local audit runs. No discovery artifact may
be added after that audit starts.
