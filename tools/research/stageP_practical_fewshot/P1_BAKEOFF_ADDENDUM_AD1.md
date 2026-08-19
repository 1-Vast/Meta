# Stage P1 bake-off addendum AD1 — episode bank record granularity (2026-08-19)

Clarification of P1_BAKEOFF_PREREGISTRATION.md §2, frozen BEFORE the bank
artifact is finalized (the first bank build exposed the ambiguity; no arm
has consumed labels under it). Does NOT move any gate, budget, metric, or
arm definition.

- Eligibility rule (prereg §1) is authoritative: a target is eligible at
  k iff it has >= k + Q unique ligands.
- Therefore bank records are emitted PER (split, target, draw, k) for each
  eligible k in {0,1,2,3,5,10,20,40}: support(k) = the first k ligands of
  the draw's ligand-unique rng ordering; query = the next Q cells after k
  (k=0 -> query = first Q cells). The same rng ordering (keyed
  stage|split|target|draw, never by arm or k) underlies every k, so
  support(k1) subset support(k2) for k1<k2 within a target-draw.
- The prereg's "query after max_k, same query across k" clause is replaced
  by the above per-k emission (it conflicted with the frozen eligibility
  rule for targets with fewer than max_k + Q ligands). Cross-k equality of
  query rows is not required for arm pairing: every arm consumes the SAME
  records at the same k, so comparisons remain paired per k.
- All other clauses unchanged. This addendum is SHA-frozen.
