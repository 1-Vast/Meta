# LEXOR L0B2 Discovery Preregistration: Minimal OpenAlex Projection

Date: 2026-07-27

## Scope and firewall

This is an independent replacement for the zero-data L0B1 transport stop. It
uses the same bounded metadata-discovery question, candidate queries, and
seven-request budget as L0B1. It may not read article full text, supplementary
files, measurement values, structures, FORT labels, `.env`, or call a model
API.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`

## Registered correction

The prior stops preserved HTTP status but not the OpenAlex error body. They do
not establish which parameter caused HTTP 400. L0B2 therefore makes the
smallest bounded projection change that removes all nested and optional
OpenAlex fields: its OpenAlex `select` value is exactly `id,doi,title`.
It retains the endpoint, filter, ordering, queries, Zenodo source, accepted
license list, eligibility rule, and seven-request budget. The runner cannot
infer a query-depth count, and no title or metadata aggregate may be used as a
substitute.

## Frozen discovery protocol

Plan: `manifests/lexor_l0b2_discovery_plan.v1.json`

<!-- LEXOR_L0B_DISCOVERY_PLAN_SHA256: a94a0be96ee6e87690123abd3ff4095780a01ff0c1e0f1399333ad0096d15376 -->

The runner may make exactly seven unauthenticated HTTPS GET requests under this
plan. It records only approved bibliographic/repository fields and may not
infer a scaffold-diverse query count. Any transport error stops this L0B2 run
without a partial inventory.

## Transition

If discovery completes, a separate L0B2 audit preregistration must freeze the
candidate-inventory hash before a local L0 audit runs. No discovery artifact
may be added after that audit starts.
