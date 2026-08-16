# Phase 22 Target Alignment Audit

## Verdict

`TARGET_ALIGNMENT_INVALID`

The repair removes the Phase-20 canonical/risk-optimal target switch, but its central regularity theorem and its claimed end-to-end scope are not mathematically established.

## Consolidated findings

### 1. Single target

Pass. PT-1, PT-3, PT-9, PT-10, and PT-11 consistently refer to the risk-optimal coefficient map

```text
g*(z) = argmin_{c in C} ell(z,c).
```

The earlier canonical operator is not used as an imitation target. When a ridge is used, the package explicitly changes the sole target to the regularized risk minimizer `g*_mu`; this is a declared bias, not a hidden target switch.

### 2. Meta-learning input

Conditional pass with narrower scope than claimed. The statistic is explicitly

```text
zeta = z(S_T, Q_T, gamma),
```

so the map is support- and query-conditioned and is not ordinary regression that ignores the support set. However, `z_H` is not an argument of `z` or `g*`. The population band and anchors are fixed by the deployment skeleton. Therefore the theorems define a separate target for each frozen historical state/skeleton; they do not establish one operator conditioned on a varying observable input `(z_H,S_T,Q,gamma)`. No equivalence between fixing `z_H` and including it in the statistic is proved.

### 3. Risk target

Existence and uniqueness hold if the strong-convexity statement in A-SC is taken as a direct assumption. The accompanying derivation is not valid as written. The assembly

```text
(1-lambda)b_pop + lambda * sum_j w_j b_j
```

is bilinear, not affine, in `c=(lambda,w)`. Consequently convexity of the band loss does not imply convexity in `c`, and adjoining `mu/2 ||c||^2` does not by itself guarantee `mu`-strong convexity. For example, a scalar squared loss containing `(lambda*w-a)^2` has an indefinite Hessian near the origin for suitable `a`; a ridge must dominate its negative curvature, and its resulting strong-convexity modulus is not automatically the ridge coefficient.

Measurability is also overstated in PT-2. A conditional expectation defines the local risk only almost everywhere until a version is selected. A-CONT can provide a continuous, hence measurable, version, but PT-2 omits A-CONT from its hypotheses while claiming a function on every `z in Z`.

### 4. Calibration and approximation

Fail. From the displayed PT-6 inequalities one obtains only

```text
||g*(z)-g*(z')|| <= 2 * sqrt(varpi_ell(d_Z(z,z')) / mu).
```

The asserted refinement to `varpi_ell/mu`, and hence the claim that `g*` is Lipschitz whenever `varpi_ell(delta)=L delta`, does not follow from uniform function-value continuity plus strong convexity. A linear argmin modulus needs stronger gradient/subgradient regularity. PT-6 itself first derives the square-root bound and then replaces it with the unsupported linear bound. PT-9 subsequently uses that false modulus. Continuity and interpolation approximability can survive with the square-root modulus, but the stated theorem and rate are false.

The coefficient-level calibration inequality

```text
||F-g*||_L2 <= sqrt(2 * excess_risk / mu)
```

is valid under genuine strong convexity. Its transfer to a vanishing operator error is only conditional on a stability inequality without an additive Route-B mesh floor. The text mentions such a floor but omits it from `Phi`; with a fixed positive additive floor, `Phi(t)` does not tend to zero.

PT-11.1 also lets approximation resolution, parameter dimension, and sample size vary simultaneously without imposing the required sieve-growth condition. `epsilon_approx -> 0` and `Gamma_N -> 0` do not follow together unless, for example, the chosen dimensions satisfy an explicit complexity-versus-sample-size rate.

### 5. Continuous regression versus ranking

The conditional strong-convex calibration can apply to continuous point-valued affinity regression only when the observable target, sufficient support/query statistic, genuine strong convexity, risk-field continuity, and operator stability assumptions all hold. A ridge defines a regularized affinity target, not the unregularized regression target.

No ranking target, coherent joint ordering object, ranking loss, or ranking calibration theorem appears in the repair. Pairwise, listwise, and metric-ranking objectives generally do not satisfy the stated uniqueness and strong-convexity conditions. None of PT-1 through PT-11 is therefore established for ranking metrics.

## Minimal obstruction

PT-6 is false as stated: its assumptions imply at most a square-root argmin modulus, not the claimed Lipschitz modulus. This invalidates a theorem used by PT-9 and the asserted complete chain. The omitted `z_H` typing, non-affine assembly/strong-convexity gap, Route-B metric floor, unconstrained sieve growth, and absence of any ranking theorem independently prevent the broader completeness claim.
