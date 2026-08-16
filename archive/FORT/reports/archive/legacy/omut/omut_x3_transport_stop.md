# OpenMut `OMUT-X3` transport stop

**Recorded:** 2026-07-28.
**Bound preregistration:** `reports/active/omut_x3_preregistration.md`,
SHA-256
`967e5878d82a1e0f47c57f0d7e18e0b384433fb4c8addc9816410cacf917dcda`.

## What happened

The first formal OpenAlex Works request returned HTTP `429`. The response
reported an anonymous daily budget of zero and a retry interval extending to
the next UTC budget window. The run was stopped rather than sleeping for
hours, inventing an API key or email address, or changing the registered
source during execution.

No `omut_x3.json` result was written. No candidate topology, document body,
abstract, snippet, activity outcome, Davis value, or sealed-test value was
read or inferred by this incomplete run.

## Interpretation

This is a transport/budget failure, not evidence that admissible OpenAlex
locations do not exist and not a topology failure. Therefore none of the
three registered X3 verdicts is asserted.

OpenAlex remains a separately registered, incomplete stage. A different
public metadata source may be evaluated only under a new preregistration with
its own source-specific admissibility rules.
