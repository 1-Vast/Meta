# Model package

This package retains verified operator primitives and the current research
model implementation:

- `bands.py`, `mathematical.py`, `meta_operator.py`, `config.py` and
  `runtime.py` implement the frozen probability-law operator;
- `encoders.py` and `mechanism.py` implement the passed sequence/graph geometry
  bridge;
- `metasieve_v1.py` implements the development-only Cold Target V1 model,
  including the retained uncentered positive dual-ridge Meta-Section.

The final task is unseen-target few-shot affinity prediction. V1 is assembled,
trainable and CUDA-verified, but is not biologically admitted or production
ready. Its best development result and failed Gates are summarized in
`PROJECT_SUMMARY.md`. Migration remains conditional on partner, affinity,
transfer and k-shot Gates.

Removed support encoders, arbitrary biological latents and terminal-negative
frontends remain recoverable from Git. Adding another PLM, attention branch or
adapter without a falsified information Gate is outside the model contract.
