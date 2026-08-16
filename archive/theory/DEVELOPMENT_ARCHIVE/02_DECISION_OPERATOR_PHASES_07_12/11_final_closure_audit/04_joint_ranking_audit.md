# Joint Ranking Object Audit

## Core object

The population ranking outcome is correctly typed on the complete space

\[
\Omega_m=S_m
\]

or on a declared weak-order analogue when ties are retained. A learned object is a constraint class of laws \(P\in\Delta(\Omega_m)\); queried orders, pairwise comparisons, and top-\(k\) statements are events constraining that law. They are not treated as outcomes themselves.

Pairwise probabilities are correctly identified as marginals

\[
p_{ab}=P\{\pi:a\succ_\pi b\}.
\]

The repair also correctly distinguishes pairwise-decomposable losses, such as Kendall discordance, from non-decomposable listwise losses that require the full law.

## Marginal-inconsistent example

For three items, set

\[
p_{12}=p_{23}=p_{31}=1.
\]

Each pair is individually a valid deterministic marginal, but together they assert the cycle \(1\succ2\succ3\succ1\), which no law on total orders can realize. The DC-R1/DC-R2 polytope rejects this example through transitivity/dicycle constraints. This part of the repair passes.

## Counterexample to DC-R5's “iff”

DC-R5 says Tier-2 decision robustness holds **iff** one action's marginal risk upper bound lies below every other action's marginal risk lower bound. Interval separation is sufficient, but not necessary because risks across a common \(P\) are correlated.

Let the ambiguity class contain all mixtures \(P_p\), \(p\in[0,1]\), on two ordering outcomes \(\omega_0,\omega_1\). Give two ranking actions losses

\[
L(a_0,\omega_0)=0,quad L(a_0,\omega_1)=1,
\]

\[
L(a_1,\omega_0)=0.1,quad L(a_1,\omega_1)=1.1,
\]

and assign all remaining actions larger losses. Then

\[
r_0(P_p)=1-p,qquad r_1(P_p)=1.1-p.
\]

Action \(a_0\) is uniquely Bayes-optimal for every \(P_p\), so it is decision-robust in the DE-R6 sense. Yet its risk interval \([0,1]\) is not below \(a_1\)'s interval \([0.1,1.1]\). The stated necessary-and-sufficient test fails.

The correct Tier-2 test is common-argmin stability over \(P\in\widehat{\mathcal Q}\), equivalently checking the coupled inequalities \(r_{a_0}(P)\le r_a(P)\) for every \(P\) and action \(a\). Separate marginal risk intervals discard the dependence needed for necessity.

## Verdict

`FAIL`

The coherent joint probability object is repaired, but a theorem used to classify robust listwise decisions remains false.
