# LEXOR L0B11 Discovery Preregistration: Verified Anonymous Transport

Date: 2026-07-27

## Scope and firewall

L0B, L0B1, L0B2, and L0B6 are retained as zero-data HTTP 400 stops and will
not be rerun. The subsequent protocol-only sequence identified the actual
issue: unauthenticated Zenodo `size=100` is rejected, while the corresponding
`size=25` request is accepted. L0B11 is a fresh, independent, no-retry
metadata-discovery run based on that verified transport shape.

It may read only the seven registered public bibliographic/repository metadata
responses. It may not read article full text, supplementary files,
measurements, structures, FORT labels, `.env`, or call a model API.

* `sealed_test_consumed=false`
* `confirmation_labels_read=false`
* `llm_api_called=false`
* `raw_measurement_files_read=false`
* `model_trained=false`

## Frozen discovery protocol

Plan: `manifests/lexor_l0b11_discovery_plan.v1.json`

<!-- LEXOR_L0B_DISCOVERY_PLAN_SHA256: f5715a1052c042ef2331fe224341acf0a7d4cb55d1c89c3ad5080b25d7813ff9 -->

The runner may make exactly seven unauthenticated HTTPS GET requests under this
plan. It records only approved metadata fields and may not infer a
scaffold-diverse query count. Any transport error stops L0B11 without a partial
inventory; the runner now reports the failing URL, HTTP status, and bounded
protocol-error detail.

## Transition

If discovery completes, a separate L0B11 audit preregistration must freeze the
candidate-inventory hash before a local L0 audit runs. No discovery artifact
may be added after that audit starts.
