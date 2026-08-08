# Phase 8.1 Repair Verification: Summary

## Verdict

PROCESS_VERDICT: `THEORY_FREEZE_CONFIRMED`

MODEL_COMPILATION_VERDICT: `DECISION_OPERATOR_REPAIR_FAILED`

## What was repaired successfully

- The exact-set minimax value remains the information floor; an outer-envelope value is now correctly typed as an upper robust surrogate, never a lower floor.
- Selection returns an argmin set or uses an explicit tie-break \(\tau\); hidden measures and implicit priors are forbidden.
- \(M_\phi\) is indexed by the query, context, and decision specification.
- Concentration rates now require IID, conditional IID, or an explicit concentration assumption. Bare exchangeability is correctly rejected as a source of Hoeffding/DKW rates.
- The finite-family union bound gives simultaneous event confidence.
- Off-coverage behavior is loss-typed; unbounded query values no longer imply universally infinite risk.

## Why the repair still fails

1. **DR-S4-R is false as stated.** Its proof assumes all intermediate minimizers remain in two separated branches, but the theorem statement omits that condition. The path \(\rho_t(a)=(a-t)^2\) on \([0,1]\) satisfies the written endpoint hypotheses and has the continuous selector \(a=t\).
2. **The abstention contract can select a worse action.** If strict ranking actions have robust loss 1, abstention costs 2, and tolerance is \(1/2\), the contract mandates abstention even though its loss is 2 and it is not criterion-optimal. Failure to meet tolerance should be flagged; abstention should be selected only when its declared criterion makes it optimal.
3. **Simultaneous intervals are mis-typed as a joint law.** A law on a strict subset \(S\subsetneq S_m\) cannot represent probability mass outside \(S\), and overlapping pairwise events are not categorical outcomes. The ambiguity class must live on the full order space, include a residual outcome, and enforce joint constraints.
4. **Current-observation conditioning is not proved.** DR-L3-R estimates transported population event frequencies, but DR-M1-R uses them as a law conditioned on current observations. This requires a joint pushforward, declared likelihood, or a proved conditional-fiber construction.

## Further corrections required

- Use \(R_{\rm set}\) only for deterministic policies and \(R_{\rm rand}\) for randomized policies.
- Require argmin existence before calling \(G_{\rm cert}\) achieved; otherwise report the selected action's outer risk or \(G_{\rm cert}+\eta\).
- Formally enforce that the archive feasibility channel ignores task ordering and duplicate multiplicity, preventing frequency information from entering identification.

## Readiness conclusion

A separate engineering agent would still need to invent theorem-level rules for conditioning, joint ranking probabilities, abstention failure behavior, policy randomization, and approximate attainment. The repaired interface is therefore not ready for model compilation.

The complete evidence and counterexamples are in `FINAL_VERDICT.md` and the five component audit files in this directory.
