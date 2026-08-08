# Phase 8.1 Repair Verification

## Process result

`THEORY_FREEZE_CONFIRMED`

The 26 Phase 0-6 files still match the prior SHA-256 snapshot; mismatches: `0`. The latest Phase-7 write is `2026-08-03T08:18:51Z`, before the repaired directory was created. The repair is isolated in eight Markdown files, and all explicit local links in those files resolve.

## Required repair checks

| Check | Result | Finding |
|---|---|---|
| A. Outer-envelope semantics | `PASS` | DR-F4-R correctly makes \(R_{\rm set}(J)\) the exact deterministic floor and \(R_{\rm set}(\widehat J)\) an outer robust surrogate/upper bound, never a lower floor. Inner witnesses are correctly typed as lower-bound devices. |
| B. Selection | `PASS` for the required selector | The output is an argmin set or an explicit \(\tau\)-selection. Hidden measures and implicit priors are forbidden. The separate DR-S4-R repair remains false as stated; see below. |
| C. Abstention | `FAIL` | The basic values are loss-typed, but the terminal failure rule can mandate abstention when it is costlier and not criterion-optimal. |
| D. Population operator typing | `PASS` for the signature | \(M_\phi\) now takes \(\gamma=(g,Q,\text{context},\mathcal A,L,\ldots)\), or returns an indexed family. Conditional-law coverage is not established by DR-L3-R; see the composition audit. |
| E. Learnability assumptions | `PASS` | Rates require IID, conditional IID, or a separately declared concentration inequality. Bare exchangeability is explicitly insufficient. |
| F. Joint-object learning | `FAIL` as a law-valued claim | The union bound proves simultaneous event intervals, but the asserted “polytope of laws on \(S\)” is not defined correctly when \(S\subsetneq S_m\) or when \(S\) is a family of overlapping pairwise events. |
| G. Failure contract | `PASS` on loss-typed infinity | Off-coverage gives \(+\infty\) only when the declared loss is unbounded on the feasible set. The action response still fails Check C. |

## Remaining mathematical failures

### 1. DR-S4-R omits its load-bearing hypothesis

`honest_selection_operator.md` states that endpoint argmins in separated closed sets force every selector to jump. Its proof additionally assumes

\[
\bigcup_t\mathcal A^*(P_t)\subseteq\mathcal A_0\cup\mathcal A_1,
\]

but that assumption is absent from the theorem statement.

Counterexample: let \(\mathcal A=[0,1]\), \(\rho_t(a)=(a-t)^2\), \(\mathcal A_0=\{0\}\), and \(\mathcal A_1=\{1\}\). All stated hypotheses hold, while \(s(t)=t\) is a continuous unique selector. The intended theorem becomes valid only after formally requiring that all minimizers remain in the two separated branches, or equivalently restricting the action space to their union.

### 2. The abstention failure rule contradicts the action semantics

For pairwise 0-1 ranking with both signs feasible, give strict actions worst-case value 1, abstention cost \(c=2\), and declared tolerance \(T=1/2\). Then \(G_{\rm cert}=1>T\). Contract Failure 2 says to abstain, but abstention has loss 2 and is not criterion-optimal. The honest output is a tolerance-failure flag (and, if an action is still required, the declared criterion's minimizer), not forced abstention.

Abstention may be selected only when its declared loss makes it optimal, or through a separately declared refusal policy with its own semantics.

### 3. Simultaneous event coverage is not yet a joint law object

DR-L3-R's Hoeffding/union calculation is valid for the displayed event intervals. However, if \(S=\{123\}\subsetneq S_3\), a “law on \(S\)” assigns mass 1 to 123 and cannot represent a true probability \(p(123)<1\). A residual outcome is required. For multiple pairwise comparisons, the events overlap, so their marginal intervals likewise do not define a categorical law on the set of pair labels.

The ambiguity object must instead be a set of laws on the full outcome space (for listwise ranking, \(S_m\), or \(S\cup\{\text{other}\}\)) constrained by the simultaneous event bounds and joint realizability.

## Additional scope corrections

- DR-F4-R(c)'s statement that “no rule” beats \(R_{\rm set}(\widetilde J)\) is valid only for deterministic rules. Randomized rules are bounded by \(R_{\rm rand}(\widetilde J)\); the fair coin in the two-sign ranking problem attains \(1/2<1=R_{\rm set}\).
- Calling \(G_{\rm cert}=\inf_a\sup_{\widehat J}L(a,v)\) an achieved guarantee requires argmin existence. For \(\mathcal A=(0,1)\), \(J=\{0\}\), and \(L(a,0)=a\), the infimum 0 is not attained. Under approximate selection, the emitted guarantee must be the selected action's outer risk, or at most \(G_{\rm cert}+\eta\), not an unattained \(G_{\rm cert}\).

## Repair-completeness verdict

`FAIL`
