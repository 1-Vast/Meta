# F6I to frozen-theory bridge

## Admitted experimental coordinates

F6I separates the task statistic into:

- a zero-shot, protein-dependent interaction coordinate `b(P,L)` supplied by
  the KLIFS/KiSSim bioactivity atlas;
- a bounded support-dependent location `tau(S)` that is independent of the
  protein representation; and
- a bounded reliability coordinate `c(S)` used only to control law dispersion.

The affine point `y=clip(base+b+tau,0,1)` is not emitted directly.  It is used
to parameterize a probability law on the frozen seven-point response mesh.

## Concrete operator path

For mesh `t=(0,1/6,...,1)`, `F(z)` is the unique discrete exponential-family
law

\[
 p_j \propto \exp\{\eta t_j-\lambda(c)(t_j-y)^2\}
\]

whose barycentre is exactly `y`.  The scalar `eta` is found by deterministic
bisection.  `B(c)` is a nonnegative, column-stochastic, tridiagonal diffusion
kernel: each interior mesh atom retains mass `1-d` and sends `d/2` to each
neighbour.  Endpoint atoms remain fixed.  Every column preserves its mesh
barycentre, so `B(c)` preserves both total probability and the predicted mean.
Finally `K(beta)` is the categorical law on the fixed mesh.

This realizes the required path

\[
 z\longmapsto F(z)\longmapsto B(z)F(z)\longmapsto K(\beta)
\]

without a scalar output bypass.  Low confidence increases the diffusion
strength and therefore the law's spread while preserving the affinity mean.

## Status boundary

This bridge is a research implementation demonstrating mathematical
compatibility.  It does not edit the frozen theory or production model, and it
does not turn the consumed PKIS2/Anastassiadis development evidence into an
untouched external validation.  A production admission still requires a fresh
endpoint-consistent panel.
