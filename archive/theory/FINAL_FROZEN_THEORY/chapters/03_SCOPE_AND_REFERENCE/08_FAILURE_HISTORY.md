# Frozen Theory: Rejected Formulations

This chapter records failed formulations for provenance only. None is part of the active definitions or theorems.

## 1. Negative ridge

**Rejected form**

$$L_0(z,B(z)p)-\frac\mu2\|p\|^2.$$

**Reason rejected:** subtracting the quadratic does not supply $\mu$-strong convexity. When the base risk is flat in a direction, the negative quadratic is strictly concave there. Existence may survive by compactness, but uniqueness, continuous argmin control, coefficient calibration, and the retained consistency chain do not follow.

**Retained form:** the positive ridge in $J_\mu$.

## 2. Old coordinate target

**Rejected form:** applying an Euclidean ridge to the earlier coordinates $(\lambda,w)$ and identifying its minimizer with the barycentric-$p$ minimizer.

**Reason rejected:** Euclidean ridge is not invariant under $p_0=1-\lambda$, $p_j=\lambda w_j$. Equal assembled-band images do not imply equal regularized objectives or minimizers.

**Retained form:** one target defined directly in $p\in\Delta_m$.

## 3. Support-intersection output

**Rejected form:** emitting the class $K(B(z)F(z))$ after an additional support/identified-set intersection while calibrating the unrestricted class.

**Reason rejected:** the declared output and calibrated object were different. Intersection need not preserve nonemptiness or be nonexpansive in Hausdorff distance.

**Retained form:** $\mathsf A(F,z)=K(B(z)F(z))$ only.

## 4. Continuum mesh convergence

**Rejected form:** inferring convergence of mesh-dependent risk minimizers to a continuum target from a same-band discretization estimate.

**Reason rejected:** a bound comparing discretized and continuum classes for the same band does not compare different minimizers. That would require a typed objective family and an argmin-stability theorem not present in the retained scope.

**Retained form:** fixed output mesh $h$ and additive floor $2h$.

## 5. Approximation from dimension growth alone

**Rejected form:** concluding $\varepsilon_{\rm approx}\to0$ merely because parameter dimension tends to infinity.

**Reason rejected:** a growing family may still contain only constant or otherwise nondense maps.

**Retained form:** the hypothesis class is the multilinear witness family itself, and density follows from node refinement plus target continuity.

## 6. Canonical/risk-optimal target switching

**Rejected form:** using a computable canonical operator for approximation or imitation while calling the task-risk minimizer the same target.

**Reason rejected:** the two objects need not coincide, and their supervision and guarantees are different.

**Retained form:** one positive-ridge regularized risk target $g_\mu^\star$ throughout.

## 7. Weak probability schedule

**Rejected forms:** using only $\delta_N\to0$, or only $\log(1/\delta_N)/N\to0$, and claiming both confidence tending to one and $\Gamma_N\to0$; claiming an almost-sure result without summability.

**Reason rejected:** the two high-probability conditions have independent roles, and Borel-Cantelli requires summable failure probabilities.

**Retained form:** both high-probability conditions, plus summability for the optional almost-sure conclusion.
