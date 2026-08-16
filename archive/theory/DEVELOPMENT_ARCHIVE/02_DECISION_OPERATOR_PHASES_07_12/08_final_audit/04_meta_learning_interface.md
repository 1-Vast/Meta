# Meta-Learning Interface Audit

## Correct object choice

The required mathematical interface is candidate D: the combination of A+B+C.

### A. Identification operator

$$
I:O\to I(O),\qquad I_Q:O\to J_Q(O).
$$

This is learned/approximated only within the frozen class, archive, coverage,
capacity, and one-sided certificate constraints.

### B. Population decision object

For a declared decision-relevant functional `psi`, history supplies an estimator
or ambiguity class

$$
H_n\to\Pi_{n,\psi}.
$$

This requires an explicit cross-task relation: IID, a specified exchangeable
model/class, or a declared transport ambiguity class. A stochastic observation
kernel is required for a point posterior; otherwise likelihood ambiguity remains
set-valued.

### C. Decision functional

$$
D:(J_Q(O),\mathcal Q(O),\mathcal A,L,\tau)
\to\arg\min_{a\in\mathcal A}
\sup_{\mu\in\mathcal Q(O)}\int L(a,v)d\mu(v),
$$

where `tau` is an explicit general tie-break, or the output remains the argmin
set.

## Contract boundary

`B` may tilt action choice inside `A`'s admissible set. It cannot narrow `I(O)` or
claim extra current-member identification. `C` must emit criterion-relative risk
separately from the identification set and its original scalar certificate.

No architecture is determined by this interface.
