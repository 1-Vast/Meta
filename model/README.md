# Model package

## Active 2026-08-13 candidate

The active cold-target model is the solver-free
**Difference-constrained Mechanism-Evidence Meta-Transformer (D-MEMT)** in
`qpsmp_meta.py`. It retains aligned BPSF interaction slots, binds support
affinity residuals to slot sensitivities, and uses support-order-invariant
mechanism prompts plus a difference-only query/reference path. It jointly
trains k={0,1,2,3,5}; k=0 is exactly zero-shot and k=1 is not structurally
disabled.

`cartesian.py` supplies an optional sparse O(3)-equivariant scalar/vector/
symmetric-traceless-rank-2 encoder. It is enabled only for declared coordinate
inputs. The current BindingDB main bank has no coordinates, so its active path
is sequence+2D and makes no atomic 3D claim.

The older HyperSAR description below is historical context for checkpoints and
ablations; it no longer defines the active architecture.

This package retains verified operator primitives and the current research
model implementation:

- `bands.py`, `mathematical.py`, `meta_operator.py`, `config.py` and
  `runtime.py` retain mathematical/operator diagnostics;
- `bpsf.py` implements the localized bipartite pair field and retained aligned
  mechanism slots;
- `geometry_supervision.py` attaches source-only contact/distance supervision to that same
  pair field. Geometry labels are never deployment inputs;
- `qpsmp_meta.py` integrates the pair trunk, mechanism-evidence set encoder,
  difference-only query transport, and scalar slot gates. The active package
  exposes no analytic solver;
- `cartesian.py` implements the optional sparse O(3) field encoder. It rejects
  cross-sample edges and must not be called with independently framed protein
  and ligand coordinates as if they were one complex.

Current validation is development-only and does not imply that the
target-specificity or Cold Target performance gates have passed. HyperSAR,
Set2Code, pooled/one-pass atom-residue, and analytic-section variants are
historical ablations only.

The Cartesian branch is PBCNet2.0/TensorNet-inspired, not a reproduction or
performance claim. The governed BindingDB corpus has no complex coordinates
for its deployment pairs.

Removed support encoders, arbitrary biological latents and terminal-negative
frontends remain recoverable from Git. Adding another PLM, attention branch or
adapter without a falsified information Gate is outside the model contract.
