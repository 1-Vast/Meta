# Model package

This package retains verified operator primitives and the current research
model implementation:

- `bands.py`, `mathematical.py`, `meta_operator.py`, `config.py` and
  `runtime.py` implement the frozen probability-law operator;
- `bpsf.py` implements the active bipartite pair field, pair-to-section latent
  encoder, and learned quotient-preserving support operator;
- `mechanism.py` attaches source-only contact/distance supervision to that same
  pair field. Geometry labels are never deployment inputs;
- `qpsmp_meta.py` integrates the pair-section trunk into the trainable
  cold-target episodic model. Analytic ridge remains a comparator only;

The current primary candidate is named **QPSMP-BPSF**. Its BPSF stage is a
bounded, GPU-friendly rectangular pair field with one-way pair-to-section
latent attention. The support-conditioned quotient operator is learned by
query loss; it is not a closed-form deployment solver. Current validation is
development-only and does not imply that the target-specificity or Cold Target
performance gates have passed.
- `metasieve_v1.py` implements the development-only Cold Target V1 model,
  including the retained uncentered positive dual-ridge Meta-Section.

The final task is unseen-target few-shot affinity prediction. V1 is assembled,
trainable and CUDA-verified, but is not biologically admitted or production
ready. Its best development result and failed Gates are summarized in
`PROJECT_CONTEXT.md`, `README.md` and `history.md`. Migration remains conditional on partner, affinity,
transfer and k-shot Gates.

Removed support encoders, arbitrary biological latents and terminal-negative
frontends remain recoverable from Git. Adding another PLM, attention branch or
adapter without a falsified information Gate is outside the model contract.
