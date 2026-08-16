# LEXOR L0B Discovery Decision

## Scope

This was a bounded public-metadata discovery phase, not LLM extraction. It made
`7` unauthenticated metadata request(s), read no full text or
measurement data, and did not read `.env` or call a model API.

## Result

* deduplicated candidate documents: `766`
* candidates with explicit accepted repository license metadata: `64`
* candidates with explicit post-firewall scaffold-diverse query count: `0`

The discovery artifact cannot make any candidate countable because the frozen
allowed metadata fields do not state the required post-firewall scaffold-diverse
query-ligand count. The candidate inventory is frozen for a separate L0B audit;
it must not be amended after observing that audit.

## Firewall State

* LLM API called: `False`
* raw measurement files read: `False`
* model trained: `False`
* FORT labels read: `False`
