# Selector Theorem Audit

## DC-S1 check

The corrected discontinuity theorem now explicitly includes all load-bearing hypotheses:

1. closed branches \(A_0,A_1\) at positive distance;
2. confinement \(\mathcal A^*(t)\subseteq A_0\cup A_1\) for every intermediate \(t\);
3. endpoint assignment \(\mathcal A^*(0)\subseteq A_0\), \(\mathcal A^*(1)\subseteq A_1\).

Under those hypotheses, the connectedness proof is valid: a continuous selector would split \([0,1]\) into two disjoint nonempty closed preimages. The oscillation lower bound follows at a common boundary point of the two preimages.

The Berge and Michael conditions in DC-S4(i)-(ii) are legitimate sufficient conditions for continuous selectors: respectively, a unique argmin under the stated compact/strictly-quasiconvex regime, and a lower-semicontinuous correspondence with nonempty closed convex values.

## Counterexample to DC-S4(iii)

The same proposition then claims that continuous selectors exist “whenever any argmin bridge connects the branches.” This is false.

Let \(\mathcal A=[0,1]\), \(t\in[0,1]\), and

\[
\rho_t(a)=(2t-1)a.
\]

Then

\[
\mathcal A^*(t)=
\begin{cases}
\{1\},&t<1/2,\\
[0,1],&t=1/2,\\
\{0\},&t>1/2.
\end{cases}
\]

The full argmin interval at \(t=1/2\) is a bridge connecting the branches, but every selector has left limit 1 and right limit 0. No continuous selector exists. A bridge persisting over a parameter interval may permit a continuous ramp, as in DC-S3; the existence of a bridge at an isolated parameter does not.

## Discrete-action scoping error

DC-S5's final conclusion is correct: any nonconstant map from a connected interval into a finite discrete action space is discontinuous. Its claim that DC-S1's confinement hypothesis holds “automatically for any two actions” is false when the finite action space has a third action that is optimal at an intermediate parameter. The discrete conclusion follows directly from connectedness, not from automatic two-branch confinement.

## Verdict

`FAIL`

DC-S1 is repaired, but the selector file still labels a false continuous-selection converse as proved. The mandate requires failure if any selector counterexample exists.
