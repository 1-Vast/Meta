# Blocking Issues

## Minimal obstruction

The closure does not define the learnable implementation family.

It introduces an implementation map:

`F_omega : Z -> C`

but does not define:

- the parameter space containing `omega`;
- the hypothesis class indexed by `omega`; or
- an empirical optimization objective for `omega).

Condition C3:

`sup_z ||F_omega(z)-g_star(z)|| <= epsilon`

is the desired approximation property itself. It is imposed as an
implementation obligation, not derived for a specified trainable class.

The existing ERM/generalization theorem applies to the separate convex
mathematical band family. No theorem connects empirical optimization of
`F_omega` to C3 or to operator-metric approximation.

Therefore the bridge from the meta-learning interface to a trainable deep
operator is undefined.

## Scope note

The scalar and ranking outputs are separately well-typed. Ranking is not
derivable from continuous affinity marginals and is only supported as a
separately supervised Route-A objective. This is an explicit scope boundary,
not the minimal invalidity causing the verdict.

