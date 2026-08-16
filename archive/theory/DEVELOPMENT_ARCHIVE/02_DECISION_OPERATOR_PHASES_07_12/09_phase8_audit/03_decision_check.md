# Decision Layer Check

## Verdict

`PASS_WITH_SCOPING_CORRECTIONS`

Phase 8 types the decision by the feasible joint object, population information, declared loss/criterion, and an optional declared tie-break. The selection documents correctly require either an argmin set or an explicit selector \(\tau\). The former implicit full-support measure is explicitly reclassified as a declaration and is forbidden when hidden.

## Criterion dependence counterexample

Let \(J=\{0,1\}\) and \(\mathcal A=\mathbb R\).

- Under deterministic minimax absolute loss, the unique action is \(a=1/2\).
- Under Bayes squared loss with declared \(P(V=1)=0.9\), the unique action is \(a=0.9\).
- Under a declared 0-1 classification loss, the Bayes action is \(1\).

Thus the same identified set legitimately produces different decisions under different declared criteria. Identification is unchanged.

## Qualifications

- `honest_selection_operator.md`, DR-S4, overstates discontinuity by saying every single-valued selector jumps along any path crossing a tie. The conclusion requires a branch-switching crossing between separated action components, as in discrete rankings; a path that merely touches a tie can retain a continuous selected branch.
- The decision tuple must include every distributional preference used by a Bayes or ambiguity criterion. Writing only an abstract \(\Delta\) is acceptable only when its semantics and conditioning are declared.

The criterion and selection typing is basically sound, but the robust-floor component used by the final decision contract is not; see `04_information_floor_check.md`.
