# F3K preregistration: source-selected identifiability sample complexity

Date frozen: 2026-08-08, after the empty-domain result for F2C at `k=5` and
before computing any certificate at a larger support budget.  KCGS numeric
outcomes remain unread.

## Hypothesis

The failure may be informational rather than architectural.  The theory only
requires `d_adapt <= k`; it does not require `k=5`.  F3K estimates the minimum
support budget at which the four-coordinate partial atlas section becomes
identifiable.

The candidate ladder is fixed as

\[
    k\in\{5,10,20,40\},\qquad
    z=(\tau,u_1,u_2,c),\quad d_{adapt}=4.
\]

For each `k`, the F2C 14-feature support-only certificate is independently
cross-fitted and one-sided 80% conformally calibrated on the same three
simultaneous scaffold-cold and kinase-group-cold PKIS1 folds.  Five deterministic
random, distinct-scaffold episodes per held target are used.  Query scaffolds
overlapping any support are excluded.

## Source-only budget choice

Choose the smallest `k` satisfying every source criterion:

1. cross-fitted certificate/observed-minimum-margin Pearson correlation at
   least `0.20`;
2. conformally admitted episode rate at least `0.20`;
3. at least 30 held target clusters have an admitted episode; and
4. the target-cluster bootstrap lower 95% bound of the observed minimum margin
   among admitted episodes is strictly positive.

If no candidate satisfies these conditions, stop without reading another
panel: the present observation law is not identifiable up to 40 supports.
No development outcome may choose the budget or the conformal threshold.

## External gate after a source pass

Only the single source-selected budget is evaluated.  PKIS2 requires at least
20% episode coverage over at least 30 targets plus positive lower 95% target-
bootstrap bounds for all six raw controls and the two correctly centred
interaction controls.  Anastassiadis2011 requires at least 10% coverage over at
least 10 targets and positive point estimates for the same eight contrasts.

Passing identifies the minimum tested observation budget for this biological
section.  It does not reinstate a five-shot claim, does not justify a scalar
bypass, and does not alter any frozen theory or production file.
