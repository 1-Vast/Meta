# MetaSieve Theory

This directory has one current theory entry point and three preserved legacy areas.

## Current theory

Read these two documents for current work:

1. [`CURRENT_THEORY/PURE_MATHEMATICAL_THEORY.md`](CURRENT_THEORY/PURE_MATHEMATICAL_THEORY.md)
   contains the architecture-independent mathematics: information limits, additive quotients,
   integrability, finite-support sections, certificates, component generalization, and transport.
2. [`CURRENT_THEORY/QPSMP_COLD_TARGET_MODEL_THEORY.md`](CURRENT_THEORY/QPSMP_COLD_TARGET_MODEL_THEORY.md)
   instantiates that mathematics as a trainable deep meta-learning model for cold-target few-shot
   drug-target affinity prediction.

The pure theory controls mathematical claims. The model theory may specialize its objects but may
not strengthen its theorems. The learnable QPSMP meta-potential is the candidate innovation; the
analytic centered ridge is only a comparator and certificate helper.

## Preserved legacy material

- `FINAL_FROZEN_THEORY/` is the previous frozen operator theory. It remains immutable provenance;
  it is not the current QPSMP theory and does not establish QPSMP cold-target claims.
- `PROJECT_HANDOFFS/` contains historical implementation handoffs. Handoffs are not theorem sources.
- `DEVELOPMENT_ARCHIVE/` contains derivations, repairs, and audits. It is retained for traceability
  and must not be used to override either current document.

## Reading rule

For a current mathematical question, read the pure theory first. For a model, training, input, or
deployment question, read the model theory second. Consult legacy material only for provenance.
