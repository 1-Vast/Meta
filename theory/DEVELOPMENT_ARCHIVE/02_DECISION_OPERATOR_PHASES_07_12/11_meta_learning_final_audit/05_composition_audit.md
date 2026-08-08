# Identification and Meta-Composition Audit

## Intended flow

The two-channel design is correct in principle:

\[
H_N\xrightarrow{A_\phi}M
\xrightarrow{\mathrm{eval}}\Delta_{\rm pop},
\qquad
(\mathrm{set}(H_N),O_*,Q_*)\xrightarrow{I_\theta}
(\widehat J,\widetilde J,\mathrm{flags}).
\]

The frequency path preserves multiplicities; the feasibility path is duplicate-invariant. Population information has no edge into identification and can only weight decisions inside the current feasible object. No hidden task identifier is an allowed input.

## Why MC-9/MC-20 do not follow

### Missing query in population evaluation

The composition evaluates \(M\) at \((\kappa(O_*),\gamma_*)\), not at \((\kappa(O_*),Q_*,\gamma_*)\). The identification arm is query-specific, while the population arm is not. Therefore the claimed pair need not concern the same decision target, and it is not generally decision-sufficient.

### Undefined ideal population target

The population target invokes a latent member not contained in the law on observable tasks. Until the joint marked-task law is typed, \(V_A\)'s phrase “true rung-r object” has no unique referent.

### Overstated support falsification

Support restriction is a safe intersection operation. A conservative confidence class can include unsupported candidate laws even when its true member satisfies sufficiency. Such breadth is not a falsification signal unless every candidate compatible with the confidence statement violates support.

## Decision layer

The decision-facing corrections referenced by Phase 9 are otherwise appropriate:

- the input includes joint \(J_Q\), a population ambiguity object, declared loss, and explicit tie-breaking;
- ranking laws live on the full order space and robustness uses the corrected per-law LP test;
- abstention is selected only when loss/criterion-optimal, while failure is a separate report;
- deterministic/randomized floors and approximate attainment are separately typed.

Those valid downstream rules cannot repair an upstream population object with the wrong index and target probability space.

## Verdict

`FAIL`

The information-flow graph is honest, but the composed mathematical objects do not share a complete common domain.
