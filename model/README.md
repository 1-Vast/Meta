# Model package

This package retains only verified production primitives:

- `bands.py`, `mathematical.py`, `meta_operator.py`, `config.py` and
  `runtime.py` implement the frozen probability-law operator;
- `encoders.py` and `mechanism.py` implement the passed sequence/graph geometry
  bridge.

The final model task is unseen-target few-shot affinity prediction, but no
assembled DTA pipeline is admitted here yet. The active research design will
meta-learn one low-dimensional biological task basis from source targets and
use a closed-form support-identifiable ridge section for adaptation. It will be
migrated only after partner, affinity, transfer and k-shot Gates pass.

Removed support encoders, arbitrary biological latents and terminal-negative
frontends remain recoverable from Git. Adding another PLM, attention branch or
adapter without a falsified information Gate is outside the model contract.
