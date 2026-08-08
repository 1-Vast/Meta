# Abstention Semantics Audit

## Verdict

`PASS`

The closure repair correctly places abstention inside the action space:

\[
a_{\rm abs}\in\mathcal A,
\qquad
\rho(a_{\rm abs})=\sup_{v\in J}L(a_{\rm abs},v).
\]

For constant minimax abstention cost \(c\), DC-A1 proves exactly

\[
a_{\rm abs}\in\arg\min_{a\in\mathcal A}\rho(a)
\iff
c\le R_{\rm set}(J,\mathcal A_{\rm strict},L).
\]

Thus the only permitted threshold is the derived comparison between the strict-action risk and the declared abstention loss. There is no generic rule of the form “uncertainty exceeds an unrelated threshold, therefore abstain.” The repair explicitly adopts counterexamples to radius thresholds and mismatched tolerance thresholds.

Failure is now a separate, out-of-game report. In particular, when

\[
\min_{a\in\mathcal A}\widehat\rho(a)>T,
\]

the result is tolerance infeasibility, not forced abstention. This correctly handles the previous witness with strict value 1, abstention cost 2, and tolerance \(1/2\).

The outer-envelope version may select abstention more often than the exact problem because it optimizes the certified surrogate, but it remains criterion-optimal for that declared surrogate and its emitted loss is a valid upper guarantee. It is not presented as exact-current-task optimality.
