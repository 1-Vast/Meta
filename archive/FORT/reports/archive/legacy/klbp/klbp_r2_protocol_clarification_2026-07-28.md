# K-LBP v2 R2 pre-live protocol clarification

**Recorded:** 2026-07-28, before any live R2 API call.

R2 was hardened after implementation review, but the formal live stage was not run because R3 returned
`R3_ESTIMATOR_INSENSITIVE_NO_DECISION` and blocked R4. API calls consumed by R2: **0**.

## Schema boundary

The model-emitted card contains exactly the 16 Section A-D fields frozen in
`klbp_r2_preregistration.md`. Every field must be present; an abstention is an explicit null. A
non-null field uses the singular `source_id` required by that stage-specific preregistration. Section E
metadata from `task.md` is runner-owned provenance, never model output and never a predictive feature.

## Artifact and budget boundary

Every future live run must bind its report, cards, and per-call ledger with one `run_id` and SHA-256
hashes. Mock paths are resolved and prohibited from aliasing live paths. The client reserves prompt and
output tokens before each request, fails closed when usage is missing, and stops all subsequent calls
after any budget error. Only `MOONSHOT_API_KEY`, `MOONSHOT_BASE_URL`, and `KIMI_MODEL` may be loaded from
the root `.env`; values are never recorded.

## Fake-cell statistic

Fake target-ligand pairs are verified absent from the complete panel using only target/ligand identity
columns. Their `y_ref` is the sampled real ligand cell's value and exists solely to apply the same
numeric-recovery rule to an unmeasured-pair negative control; it is not asserted to be the fake pair's
unknown affinity. This definition must be retained if R2 is ever reauthorized.

No threshold, prompt, call cap, token cap, sample size, or scientific gate changed in this
clarification.
