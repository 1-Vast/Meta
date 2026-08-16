# Task Space and History Audit

## Task object

`PASS`

Phase 9 defines

\[
T=(O,S,Q,\gamma),
\]

with a complete observation record, released support, query points and pushforward, and a declared decision specification. Historical tasks may omit \(Q\) or \(\gamma\); the current task requires all components.

The task space is declared standard Borel, which is an acceptable explicit measurability assumption for regular conditional laws on variables actually included in the probability space.

## Historical sample

`PASS`

The history is correctly typed as an ordered sequence

\[
H_N=(T_1,\ldots,T_N)\in\mathbb T^N,
\]

with multiset and set quotients. The feasibility channel consumes the set of distinct traces; the population-frequency channel consumes the multiset, or the full sequence under declared drift.

Two histories with the same distinct elements can carry different meta-information. For records \(A,B\),

\[
H_2=(A,B),\qquad H_3=(A,A,B)
\]

have the same set quotient but empirical masses \(1/2\) and \(2/3\) on \(A\), different sample sizes, and different confidence radii. Phase 9 correctly preserves this difference.

## Downstream indexing warning

The task definition places the query points and demanded pushforward \(g\) in the separate object \(Q\), while the transferable operator is later indexed only by \((c,\gamma)\). Unless the complete \(Q\) is explicitly duplicated inside \(\gamma\) with a coherence rule, the task information is lost at the meta-operator boundary. This is a downstream interface failure, not a defect in \(T\) or \(H_N\) themselves.
