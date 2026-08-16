# Mamba E0 environment/component decision

Date: 2026-07-27

## Verdict

`MAMBA_E0_PASS_COMPONENTS_AUTHORIZED_ONLY`.

The `drug` environment now imports `mamba_ssm 2.2.4` with `transformers 4.46.3` and executes the
selective-scan `Mamba` kernel on the RTX 4060 GPU. The prior `transformers 5.12.1` dependency was
incompatible with `mamba_ssm` because it removed generation output symbols imported by the package.

The E0 runner verifies a Mamba-as-FFN forward/backward pass, two independently trainable directional
scans in `TrueBiMamba`, and the tied-direction component swap under sequence reversal. This rules out
the one-way scan incorrectly presented as bidirectional in some public DTA implementations.

This is an implementation result only. It provides no evidence of target-specific affinity transfer
and does not authorize any Mamba-DTA training before F0-P and F1 pass.

