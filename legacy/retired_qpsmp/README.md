# Retired QPSMP Paths

This isolated package preserves the pre-BPSF development comparators. None of
these classes are imported by the active `model` or `scripts` packages.

- `PooledInteraction`: ligand-conditioned pooled protein interaction.
- `AtomResiduePool`: one-pass atom-residue attention and pooling.
- `CenteredRidgeSection`: analytic positive-ridge support solve.
- `qpsmp.py`, `metasieve_v1.py`: former analytic Cold Target models.
- `scripts/` and `tests/`: their standalone historical runner and direct tests.

They are retained for historical reproduction only. The active model is the
trainable QPSMP-BPSF implementation.

This directory is intentionally outside both active packages. Its files may
require import-path adjustment to reproduce an old run and must not be copied
back into `model/` or `scripts/` without a new admission review.
