# LEXOR L0B Transport Attribution Correction

Date: 2026-07-27

## Correction

The original L0B, L0B1, L0B2, and L0B6 discovery runners wrote output only
after all planned requests completed. Their HTTP 400 exceptions contain neither
the URL nor a persisted request index. Accordingly, those zero-data stops
establish only that **at least one** registered metadata request failed. They
do not establish that the failure occurred on the first request, on OpenAlex,
or in an OpenAlex `select` projection.

The original frozen preregistrations remain unchanged. This correction
supersedes any causal reading of their retrospective prose that attributes a
400 to specific OpenAlex fields or endpoint parameters.

## Consequence

L0B3-L0B7 are transport-only diagnostics and do not amend a discovery
inventory. A subsequent discovery run is justified only after a bounded,
per-endpoint probe identifies a compatible request shape and retains no
candidate response body.

## Firewall State

* candidate metadata accepted by this correction: `False`
* raw measurement files read: `False`
* LLM API called: `False`
* model trained: `False`
* FORT labels read: `False`
* sealed test consumed: `False`
