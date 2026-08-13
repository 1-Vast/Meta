# Model package

This package retains verified operator primitives and the current research
model implementation:

- `bands.py`, `mathematical.py`, `meta_operator.py`, `config.py` and
  `runtime.py` implement the frozen probability-law operator;
- `bpsf.py` implements the localized bipartite pair field and support-gated
  low-rank interaction adapters;
- `geometry_supervision.py` attaches source-only contact/distance supervision to that same
  pair field. Geometry labels are never deployment inputs;
- `qpsmp_meta.py` integrates the pair trunk and amortized Set2Code conditioner into the trainable
  cold-target episodic model. The active package exposes no analytic solver;

The current primary candidate is named **QPSMP-HyperSAR**. Its support set is
amortized into a transient reference code. A lightweight Siamese relative
conditioner combines it with each query ligand and protein--ligand endpoint,
then gates shared low-rank bases in the last interaction blocks with a
query-specific code. The code is learned by query loss and discarded after
the episode; active few-shot adaptation has no ridge/linear solve and no
test-time gradient update. Current validation is
development-only and does not imply that the target-specificity or Cold Target
performance gates have passed. Retired pooled, one-pass atom-residue and
analytic-section comparators are historical research only and are
never imported by this package or active training scripts.

The PBCNet2.0 connection is limited to reference--query relative recognition:
the support reference is compared with every query before pair modulation. It
is not a reproduction of PBCNet2.0's coordinate-based Cartesian rank-2 tensor
network, because the governed cold-target corpus has no complex coordinates for
every deployment pair.

Removed support encoders, arbitrary biological latents and terminal-negative
frontends remain recoverable from Git. Adding another PLM, attention branch or
adapter without a falsified information Gate is outside the model contract.
