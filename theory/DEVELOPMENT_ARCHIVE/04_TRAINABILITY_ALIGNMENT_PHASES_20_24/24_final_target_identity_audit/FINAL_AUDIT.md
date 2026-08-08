# Phase 24 Final Target Identity Audit

## Verdict

`FINAL_THEORY_INVALID`

## Consolidated audit

### 1. Single target consistency

Pass for the fixed-mesh theory. TI-1, TI-2, TI-5, TI-6, TI-7, TI-8, and TI-9 consistently use the barycentric regularized target `g_mu^*`. TI-3 expressly retracts equivalence with the old-coordinate target, and later theorems do not use `g_old_mu`. The unregularized `g_0^*` appears only to disclaim equality and optimality, not as a substituted target.

### 2. Regularized target validity

Pass. The operative objective is explicitly

```text
J_mu(z,p) = L_0(z,Bp) + (mu/2)||p||^2,
g_mu^*(z) = argmin_{p in Delta_m} J_mu(z,p).
```

TI-5 proves pointwise and population optimality only for the matched regularized risk `R_mu`. The package no longer claims that `g_mu^*` equals or minimizes the unregularized risk.

### 3. Strong convexity

Pass under the declared assumptions. For each fixed `z`, `p -> B(z)p` is linear, the conditional base loss is convex in `p`, and the Euclidean ridge makes `J_mu(z,.)` at least `mu`-strongly convex on the simplex. Compactness, continuity, and A-CONT give existence, uniqueness, and an everywhere measurable target.

### 4. Mesh refinement

Fail. TI-10 declares mesh-indexed targets `g*_{mu,h_N}` and a continuum target `g*_{mu,0}`, but it does not prove the asserted convergence between them.

The bound

```text
d_M(K_h(beta), K_0(beta)) <= 2h
```

compares discretized and continuum operators for the same band `beta`. It does not imply

```text
d_M(K_h(beta*_{mu,h}), K_0(beta*_{mu,0})) <= 2h,
```

because the two minimizers can be different. That conclusion needs a uniformly convergent family of objectives `J_{mu,h} -> J_{mu,0}` and an argmin-stability argument. Under strong convexity, a uniform objective error of order `h` would normally yield coefficient displacement of order `sqrt(h/mu)`, not the asserted `2h` target rate.

The target spaces are also under-typed. If the grid refinement changes the band basis, then `B_h`, the number of anchors, the simplex `Delta_{m(h)}`, and maps into a common coefficient or band space must be defined. If `B` and `Delta_m` do not change, then the coefficient target is not mesh-indexed and the separate symbol `g*_{mu,h}` has not been defined. TI-10 specifies neither alternative completely.

TI-9 contains a smaller overstatement: its upper bound converging to `2h` proves only `limsup error <= 2h`, not that the error itself converges to exactly `2h`.

Thus fixed-mesh deployments and mesh-refining deployments are named separately, but the mathematical convergence typing connecting their targets is incomplete.

### 5. Calibration

Pass at fixed mesh. TI-8 measures coefficient excess risk and operator error only against `g_mu^*`, with the Route-B design floor shown separately. No old or unregularized target enters the calibration inequality. The continuum conclusion in TI-10 does not follow, however, because the missing mesh-target convergence is required for its triangle step.

### 6. DTA scope

Pass with declared limits. The package retains the fixed-deployment, continuous point-valued affinity-regression scope and disclaims ranking guarantees. It supplies no pairwise, listwise, joint-ordering, or ranking-calibration theorem.

## Minimal obstruction

The fixed-mesh target identity is repaired, but TI-10 mistakes a same-band discretization estimate for convergence of different mesh-dependent risk minimizers. Without defined mesh-indexed objectives/bases and a valid argmin-stability theorem, `g*_{mu,h} -> g*_{mu,0}` and the final zero-error consistency claim are unproved. Because that theorem is claimed as part of the repair, the final theory is not valid as written.
