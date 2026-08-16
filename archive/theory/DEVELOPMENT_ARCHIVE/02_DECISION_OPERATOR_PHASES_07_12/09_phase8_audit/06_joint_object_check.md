# Joint Query Object Check

## Verdict

`PASS`

Phase 8 correctly uses the joint pushforward \(J_Q(O)\), its difference projections for pairwise decisions, and its jointly realizable order set \(\Sigma(J_Q(O))\) for listwise ranking. It explicitly rejects reconstructing listwise behavior from independent marginals or separately realizable pairwise signs.

## Same marginals, different ordering behavior

Let

\[
J_{\mathrm{diag}}=\{(0,0),(1,1)\},\qquad
J_{\mathrm{anti}}=\{(0,1),(1,0)\}.
\]

Both have coordinatewise marginal intervals \([0,1]\times[0,1]\). In \(J_{\mathrm{diag}}\), the two coordinates are always tied. In \(J_{\mathrm{anti}}\), either strict order is possible and equality is impossible. Marginal intervals therefore cannot determine ranking behavior.

DR-J3's three-item witness is also correct: the two reverse vectors admit only two joint orders, while independently compatible pairwise signs over-admit all six orders. The Kendall minimax values are 2 for the exact order set and 3 for the pairwise outer proxy.

That value 3 is a conservative surrogate robust value, not an inflated lower information floor. The joint-object theorem is correct; its cross-reference to DR-F4 inherits the floor-direction error identified separately.
