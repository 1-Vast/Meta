# Audited Model Contract

## Status

`NOT ESTABLISHED`

This file records the contract surface actually defined by the closure. It does
not add the missing hypothesis class or propose an architecture.

## Defined inference contract

### Inputs

`(z_H,S,Q,gamma)`

- `z_H`: frozen meta-training state containing trained mathematical-family
  parameters, context-fiber counts, empirical population bands with margins,
  and tags.
- `S`: current finite observable support, noise level, and optional declared
  auxiliary label.
- `Q`: declared query index.
- `gamma`: declared decision specification.

### Fixed computations

- exact current identification object `I(S)`;
- current context `kappa(S)`;
- population-band selection from `z_H`;
- convex assembly in the valid-description set;
- support restriction by `I(S)`;
- confidence/rung tags from `z_H`;
- certificate rows and failure flags from canonical rules.

### Outputs

- scalar CDF-band law class under the bounded Route-B scope, or a separate
  Route-A order polytope;
- confidence and rung;
- identification certificates; and
- typed flags/fallbacks.

### Metrics

- scalar head: Hausdorff-`W_1` on law classes;
- ranking head: Hausdorff-TV/Hoffman geometry on order polytopes;
- confidence: absolute distance;
- rung: discrete distance.

### Training data semantics

- point-identified scalar targets may enter the calibration-bearing interval
  loss;
- censored scalar tasks contribute forced/compatible population-band
  estimation only;
- ranking consumes separately identified order supervision;
- task sampling claims require the declared IID/C-IID and transport tags.

## Undefined trainable contract

The closure does not provide a complete tuple:

`(Omega, {F_omega}_{omega in Omega}, R_hat_N, R, d_M)`.

Only `R` for the separate mathematical band family and `d_M` are defined.
For the implementation map, `Omega`, the indexed class, and
`R_hat_N(omega)` are absent.

Uniform approximation C3 and optimization tolerance are stated as obligations.
Accordingly, no exact model-engineering contract with a proved learnability
claim follows from the frozen closure.

## Scope boundary

The defined interface supports bounded scalar affinity uncertainty and a
separately supervised finite-order ranking object. It does not support a
ranking distribution derived from continuous scalar affinity marginals.

CONTRACT_VERDICT: `UNDEFINED_TRAINABLE_OPERATOR`
