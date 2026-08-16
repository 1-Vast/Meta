# Phase 23 Final Theory Freeze Audit

## Verdict

`THEORY_STILL_INVALID`

## Consolidated audit

### 1. Scope verification

Pass. DT-A explicitly fixes `z_H=z_H^0` per deployment and expressly excludes generalization over varying `z_H`. The statistic remains `z(S_T,Q_T,gamma)`, so the current result is support/query conditioned but deployment-state specific.

### 2. Target validity

Fail on target identity. RP-3 defines an everywhere measurable, unique minimizer in the new barycentric coordinate under A-CONT and genuine strong convexity, and RP-10 distinguishes a ridge-regularized target from an unregularized target. Those local properties are valid.

However, the repair also says the Phase-21 target is unchanged and that only its coordinates changed. That is false when the ridge is active. Phase 21 regularizes the old coordinate `c=(lambda,w)` by `mu/2 ||c||^2`; Phase 22.1 regularizes `p=(1-lambda,lambda*w)` by `mu/2 ||p||^2`. Euclidean norm is not invariant under this nonlinear, non-bijective-at-`lambda=0` reparameterization, so the two objectives and their minimizers generally differ.

A one-anchor counterexample is enough. Let the base loss be identically zero, which is bounded and convex. With `w=1`, the old ridge is proportional to `lambda^2+1` and is minimized at `lambda=0`. The new ridge is proportional to `(1-lambda)^2+lambda^2` and is minimized at `lambda=1/2`. Thus the old `g*_mu` and repaired `g*_mu` are different targets even though both parameterizations have the same assembled-band image. The statement that Bayes optimality transfers with the “same risk” is therefore invalid. The repair neither retracts the old regularized target nor proves equivalence to it.

### 3. Strong convexity

Pass for the new objective. For each fixed `z`, `B(z)p` is linear in `p`; dependence of `B` on the observable context does not create a bilinear parameter term. A convex band loss composed with this linear map remains convex, and adding `mu/2 ||p||^2` makes the new objective at least `mu`-strongly convex. The phrase “modulus exactly mu” is stronger than proved when the base risk is itself strongly convex, but the lower modulus `mu` used later is valid.

### 4. Approximation theorem

Pass after reading RP-4 with the necessary two-sided argument. Strong convexity at both `z` and `z'`, plus the two uniform value-modulus transfers, gives

```text
mu ||g*(z)-g*(z')||^2 <= 2 varpi_ell(d_Z(z,z')),
```

and hence the boxed square-root bound. The displayed proof first derives a looser constant and describes the refinement imprecisely, but the boxed conclusion follows from the stated assumptions. RP-5 uses a square-root/Hölder modulus and makes no default Lipschitz claim; the optional linear rate is correctly conditioned on A-GRAD.

### 5. Statistical consistency

Conditional and not fully frozen. RP-8 correctly separates

```text
d <= Phi(excess risk) + epsilon_design,
epsilon_design = 2h,
```

and states that total error vanishes only when `h -> 0`. RP-9 also states the needed dimension-growth condition and requires approximation mesh, output mesh, generalization, and optimization errors to vanish.

Two typings remain unstated. First, A-SKEL fixes the Route-B value mesh, while RP-9 changes `h_N`; a common limiting operator space, embeddings of the mesh-dependent spaces, and the corresponding single target across this changing sequence are not defined. Second, the example schedule assigns `dim Omega_N` and `h_N` but does not tie the interpolation resolution `r_N` to the actual node-count dimension of `Omega_N`; nor is a sequence `delta_N -> 0` declared for the claim “with probability -> 1.” The result is therefore only a schedule template, not the fully typed consistency theorem claimed.

### 6. DTA scope

Pass. RP-10 expressly limits the theory to continuous point-valued affinity regression for a fixed deployment and expressly disclaims every ranking guarantee. Pairwise, listwise, joint-ordering, and ranking-calibration claims remain out of scope.

## Minimal obstruction

The barycentric change makes strong convexity valid, but applying the ridge in the new coordinates changes the regularized risk target. Because the package simultaneously preserves the old target and asserts a new-coordinate target with the “same risk,” it does not contain one frozen `g*`. This must be resolved before the theory is ready for modelization; the fixed-deployment and regression-only scope limits do not cure that target mismatch.
