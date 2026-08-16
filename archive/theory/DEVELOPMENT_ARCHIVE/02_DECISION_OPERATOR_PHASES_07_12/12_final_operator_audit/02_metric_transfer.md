# Route A Metric Transfer Audit

## Result

`FAIL AS STATED; PROBABILITY-COORDINATE TRANSFER PASSES`

## Valid part

Route A explicitly declares:

- `|Omega_Q| <= n_bar`;
- at most `e_bar` event constraints at each index;
- finite VC dimension for the full indicator class; and
- pullback closure handled separately by OM-5.

With finite outcome and event bounds, the 0/1 constraint matrices range over a
finite family, up to harmless outcome relabeling. Each fixed matrix has a finite
Hoffman constant, so their maximum `H_bar` is finite. For two nonempty
polytopes with the same matrix and right-hand sides within `epsilon`,

`d_H^TV(K_1,K_2) <= (H_bar/2) epsilon`

is correctly proved. An infinite query atlas is permissible because its
constraint geometry has finitely many patterns and its statistical class is
separately controlled by finite VC dimension.

The earlier dyadic, unbounded-outcome counterexample violates RA1 and therefore
does not refute this probability-set result.

## Full-metric counterexample

OM-2 states the stronger inequality

`d_M(M_1,M_2) <= (H_bar/2) d_desc(M_1,M_2)`

for arbitrary `M_1,M_2 in M`. That statement is false after OM-3 redefines
`d_M` on probability, confidence, and rung coordinates.

Take a one-index, one-outcome space. Let both operators have the same singleton
probability set and the same rung, but confidence coordinates zero and one:

`M_1=(K,0,r)`, `M_2=(K,1,r)`.

All constraint endpoints coincide, so `d_desc(M_1,M_2)=0`. The complete
product metric gives `d_M(M_1,M_2)=1`. This satisfies the strongest possible
RA1 bound and falsifies OM-2 exactly as written.

The same construction gives a sequence with `d_desc -> 0` while `d_M` does
not converge. Bounded outcome complexity transfers only the probability-set
coordinate. Confidence convergence and rung equality require separate
hypotheses. OM-5-pre and S5 later supply those hypotheses for the canonical
sequence, but they do not make OM-2's universal full-metric statement true.

## Scope

Route A also restricts the learnability theorem to finite outcome spaces. The
frozen interface includes scalar real-valued pushforwards; those require Route B
or another declared stability result. Route A is a valid restricted contract,
not a proof for the entire earlier codomain.

