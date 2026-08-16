# Frozen Theory: Scope and Limitations

## Supported scope

The frozen theory supports exactly:

1. **One fixed deployment**
   $$\mathcal D=(z_H^0,B(\cdot),\Delta_m,\mu,h).$$
   All population bands, anchors, ridge modulus, and output mesh are fixed.

2. **Support-conditioned meta-learning**
   Historical IID tasks train $F_{\hat\omega_N}$. At inference, the current support and query enter through $z(S,Q,\gamma)$.

3. **Continuous point-valued affinity regression**
   The supervised target is $Y\in V\subset\mathbb R$. The emitted object is a fixed-grid CDF-band class of probability laws over $V$.

4. **Fixed-resolution guarantees**
   Approximation resolution for the coefficient map may refine, but the output-grid mesh $h$ remains fixed. Calibration and consistency retain the design floor $2h$.

## Not claimed

The frozen theory provides no theorem for:

- pairwise, listwise, or metric ranking;
- coherent joint-order learning;
- derivation of ranking from affinity regression;
- continuum output-mesh refinement or a zero-mesh target;
- convergence between mesh-indexed targets;
- a model conditioned on varying $z_H$;
- transport to an undeclared task distribution;
- conditional-IID fibers or missing-fiber generalization;
- optimization efficiency or an architecture-specific training algorithm;
- equality between the regularized target and an unregularized Bayes target;
- a support-intersected output class.

## Interpretation of the result

The theorem guarantees that empirical learning approaches the single regularized target operator in the declared operator metric up to the fixed $2h$ floor, with high probability under the full schedule in `06_CALIBRATION_AND_GENERALIZATION.md`. It does not claim exact zero-error recovery at fixed output resolution.

## Final status

Within these explicit limitations, the retained mathematical theory is frozen. The limitations are part of the contract and must not be removed by interpretation.
