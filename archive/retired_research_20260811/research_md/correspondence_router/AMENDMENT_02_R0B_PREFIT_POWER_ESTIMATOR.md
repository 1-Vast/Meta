# R0-B Amendment 02: Executable pre-fit power estimator

Status: frozen on 2026-08-11 after immutable exact-geometry serialization but
before any frozen-P1B prior score, additive-oracle score, power result, or R0-B
model fit was read. No affinity label is used.

This amendment resolves one implementation detail left unspecified by
`PREREG_R0_EXACT_DISTANCE_RESIDUAL.md`. It changes no population, arm, effect
threshold, Gate, seed, or authorization boundary.

The pre-fit audit is evaluated on heldout-A. Each system is first reduced to
one mean RPS, then systems are averaged inside their registered protein
homology component, and components receive equal weight.

The fixed null dispersion is the vector of frozen P1B prior RPS values across
heldout-A components. Center this vector at zero and draw 10,000 equal-size
component bootstrap means with NumPy `default_rng(1700)`. For a constant true
improvement `e`, the one-sided alpha 0.05 rejection threshold is the 95th
percentile of the centered null means. The minimum `e` giving 80% rejection is
therefore registered as:

`MDE80 = quantile(null_mean, 0.95) - quantile(null_mean, 0.20)`.

This estimator uses no fitted-model result and cannot be changed after the
ceiling report is read. The existing stopping condition remains
`MDE80 <= delta_star`, where `delta_star = 0.05*S_prior`.

