# LEXOR L0B3 Transport Probe Decision

## Scope

This was an HTTP-parameter diagnostic only. Successful response bodies were
not read, parsed, retained, or serialized; no candidate metadata was accepted.

## Results

| probe | HTTP status | bounded protocol detail |
| --- | ---: | --- |
| zenodo_q | 200 |  |
| zenodo_q_size | 400 | A validation error occurred. | size | Page size cannot be greater than 25. Please use authenticated requests to increase the limit to 100. |
| zenodo_q_size_sort_mostrecent | 400 | A validation error occurred. | size | Page size cannot be greater than 25. Please use authenticated requests to increase the limit to 100. |

## Firewall State

* candidate metadata accepted: `False`
* raw measurement files read: `False`
* LLM API called: `False`
* model trained: `False`
* FORT labels read: `False`
* sealed test consumed: `False`
