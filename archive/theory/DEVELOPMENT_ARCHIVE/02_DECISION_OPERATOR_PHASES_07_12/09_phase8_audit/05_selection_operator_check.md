# Selection Operator Check

## Verdict

`PASS_WITH_SCOPING_CORRECTION`

Phase 8 supplies both honest forms:

1. return \(\mathcal A^*=\operatorname*{argmin}_{a\in\mathcal A}\rho(a)\), or its declared \(\eta\)-argmin enlargement;
2. return \(\tau(\mathcal A^*)\) only when \(\tau\) is explicit and uses no unidentified structure or undeclared measure.

This correctly rejects a hidden prior, hidden reference measure, or unstated preference. The symmetry obstruction for two freely swapped ranking actions is valid: there is no symmetry-equivariant single-valued output, so the system must return the set or declare asymmetry through \(\tau\).

The approximation lemma is also sound: if \(\sup_a|\widehat\rho(a)-\rho(a)|\le\eta/2\), a minimizer of \(\widehat\rho\) is \(\eta\)-optimal for \(\rho\).

DR-S4 needs narrower hypotheses. A mandatory jump follows for a branch-switching tie between separated action components, including opposing discrete rankings. It is not true for every continuous problem path that happens to contain a nonunique argmin. This is a theorem-scoping defect, not a reason to reintroduce hidden selection.
