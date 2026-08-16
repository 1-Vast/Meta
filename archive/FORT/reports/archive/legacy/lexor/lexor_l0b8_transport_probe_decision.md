# LEXOR L0B3 Transport Probe Decision

## Scope

This was an HTTP-parameter diagnostic only. Successful response bodies were
not read, parsed, retained, or serialized; no candidate metadata was accepted.

## Results

| probe | HTTP status | bounded protocol detail |
| --- | ---: | --- |
| zenodo_kinase_inhibitor_profiling | 400 | A validation error occurred. |
| zenodo_kinome_profiling | 400 | A validation error occurred. |
| zenodo_drug_target_affinity | 400 | A validation error occurred. |

## Firewall State

* candidate metadata accepted: `False`
* raw measurement files read: `False`
* LLM API called: `False`
* model trained: `False`
* FORT labels read: `False`
* sealed test consumed: `False`
