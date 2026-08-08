# Robust Information Floor Check

## Verdict

`FAIL`

## Sound part

For the exact identified joint object, the definition

\[
R_{\mathrm{set}}(J,\mathcal A,L)=\inf_{a\in\mathcal A}\sup_{v\in J}L(a,v)
\]

is the correct deterministic minimax value. Its randomized analogue is also correctly typed. It supports arbitrary declared losses and structured action spaces by definition. With scalar \(J\subset\mathbb R\), \(\mathcal A=\mathbb R\), and \(L(a,v)=|a-v|\), it recovers

\[
R_{\mathrm{set}}=\tfrac12(\sup J-\inf J).
\]

For a pairwise 0-1 ranking action with both strict signs feasible and no abstention, the deterministic and randomized values are respectively \(1\) and \(1/2\).

## Fatal direction error in DR-F4

`loss_typed_information_floor.md`, DR-F4(i), correctly proves monotonicity but draws the wrong inferential conclusion. If \(J\subseteq\widehat J\), then

\[
R_{\mathrm{set}}(J)\le R_{\mathrm{set}}(\widehat J).
\]

Therefore the outer-set value is an **upper bound on the true minimax value**, not a valid inflated lower information floor.

Counterexample: take scalar absolute loss, \(J=\{0\}\), and \(\widehat J=\{0,100\}\). Then

\[
R_{\mathrm{set}}(J)=0,\qquad R_{\mathrm{set}}(\widehat J)=50.
\]

The action \(a=0\) has true worst-case loss zero, so the assertion that no guarantee below 50 is possible is false. Conversely, an inner set gives a smaller minimax value; as a lower bound on the true minimax value it is weaker but not false. An inner set is unsafe for feasible-set certification, but that is a different claim.

The correct approximation semantics are:

- exact \(R_{\mathrm{set}}(J)\): information lower bound and exact minimax value;
- outer \(R_{\mathrm{set}}(\widehat J)\): conservative surrogate minimax value and an upper risk certificate for the action optimized over \(\widehat J\);
- inner \(R_{\mathrm{set}}(\widetilde J)\): potentially weak lower bound, not a valid representation of feasibility.

## Abstention error

DR-F1(iii) claims the pairwise values \(1\) and \(1/2\) “with or without abstention.” If abstention is an action with declared constant loss \(c<1/2\), both minimax values are at most \(c\), contradicting that claim. Abstention must be included in \(\mathcal A\) with its declared loss and recomputed.

Because the false outer-floor conclusion is reused in DR-M1, \(V_D\), the compilation contract, and the stopping criterion, the Phase-8 bridge fails mathematically even though its exact-set minimax definition is sound.
