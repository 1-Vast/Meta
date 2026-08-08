# Phase 8.2 Decision Operator Closure: Consolidated Final Audit

## Executive verdict

Phase 8.2 correctly repairs several central defects: DC-S1 now states branch confinement explicitly; abstention is selected only through the declared criterion and is separated from failure; ranking probabilities live on the full order space; and the interface distinguishes deterministic/randomized floors and attained/approximate guarantees.

The claim `DECISION_OPERATOR_CLOSED` is nevertheless false. Counterexamples remain inside claimed proved results, and the final population interface destroys the sample multiplicities required by its own learning theorem.

## Audit summary

| Audit area | Result |
|---|---|
| Phase 0-7 freeze | `THEORY_FREEZE_CONFIRMED` |
| Corrected branch-confinement theorem DC-S1 | Pass |
| Complete selector file | Fail: false bridge converse and false automatic-confinement claim |
| Abstention/failure semantics | Pass |
| Full-space coherent ranking law | Pass |
| DC-R5 robustness characterization | Fail: interval separation is sufficient, not necessary |
| Conditional-fiber idea under actual-law sufficiency | Conditionally valid |
| DC-C2/DC-C3 as written | Fail: noise-kernel quantifier and singleton-posterior overclaims |
| Existence vs identification vs learning separation | Incomplete at the terminal theorem |
| Trainable interface | Fail: archive set type discards frequency multiplicity |

## Decisive counterexamples

### Selector

For \(\rho_t(a)=(2t-1)a\) on \(\mathcal A=[0,1]\), the argmin is \(\{1\}\) left of \(1/2\), all of \([0,1]\) at \(1/2\), and \(\{0\}\) right of \(1/2\). An argmin bridge exists, but no continuous selector exists. This refutes DC-S4(iii)'s claim that any bridge suffices. DC-S1 itself remains valid because its confinement premise fails here.

### Ranking robustness

Let \(r_0(P_p)=1-p\) and \(r_1(P_p)=1.1-p\) over an ambiguity class \(p\in[0,1]\). Action 0 is uniquely optimal for every population law, but its interval \([0,1]\) is not separated below action 1's \([0.1,1.1]\). Common-argmin robustness holds while DC-R5's claimed “iff” interval test fails.

### Conditional population law

`SUFF-kappa` is a statement about a particular joint law of latent member and record. DC-C2 assumes it once and concludes the conditional identity for every bounded-support noise kernel. With constant \(\kappa\), an admissible uninformative kernel can satisfy sufficiency while another admissible revealing kernel does not. Uniform validity requires a uniform sufficiency assumption. Moreover, a declared likelihood maps a prior ambiguity class to a posterior ambiguity class; it does not create the single posterior asserted by DC-C3.

### Archive type

The final interface makes \(\mathcal H\) a finite set for both operators. Duplicate invariance is correct for \(I_\theta\)'s feasibility view but invalid for \(M_\phi\)'s frequency view. Deduplication changes empirical frequencies, sample size, and confidence radii, so DR-L3-R cannot apply to the stated input object.

## Conditional population determination

The closure repair does identify a valid route in principle:

\[
\text{query-indexed actual-law }(\mathrm{C\!\!-IID}_\kappa)
+(\mathrm{SUFF}_\kappa)
\Longrightarrow
P(g(f_\beta)\mid O,Q)=P(g(f_\beta)\mid\kappa(O),Q).
\]

That route is a declared-assumption compilation target, not a derivation from identification alone. It becomes usable only after the quantifier is corrected, the historical frequency input is a multiset/sequence, the zero-fiber case has a vacuous fallback, and posterior ambiguity is retained unless the prior is singleton.

## Existence, identification, learning

- **Existence:** available under compactness/lower-semicontinuity or as an explicit \(\eta\)-argmin object.
- **Identification:** set-valued decisions are determined once all valid decision and conditional-population declarations are fixed.
- **Learning:** simultaneous event constraints are estimable under the declared IID/C-IID rungs, but the complete operator is not learnable from the interface as typed.

## Engineering handoff

An engineering agent would still have to change mathematical domains and repair theorem statements concerning selector continuity, common-argmin ranking robustness, sufficiency across noise kernels, posterior uncertainty, and empty conditional fibers. Those are theorem-level choices, so neither readiness verdict nor “closed but not learnable” is appropriate.

PROCESS_VERDICT: `THEORY_FREEZE_CONFIRMED`

MODEL_COMPILATION_VERDICT: `DECISION_OPERATOR_INVALID`
