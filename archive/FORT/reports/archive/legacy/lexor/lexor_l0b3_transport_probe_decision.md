# LEXOR L0B3 Transport Probe Decision

## Scope

This was an HTTP-parameter diagnostic only. Successful response bodies were
not read, parsed, retained, or serialized; no candidate metadata was accepted.

## Results

| probe | HTTP status | bounded protocol detail |
| --- | ---: | --- |
| openalex_base | 200 |  |
| openalex_per_page | 200 |  |
| openalex_filter | 200 |  |
| openalex_sort | 200 |  |
| openalex_select | 200 |  |

## Firewall State

* candidate metadata accepted: `False`
* raw measurement files read: `False`
* LLM API called: `False`
* model trained: `False`
* FORT labels read: `False`
* sealed test consumed: `False`
