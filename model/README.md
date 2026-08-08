# Model Package

This package contains only two verified boundaries:

- `bands.py`, `mathematical.py`, `meta_operator.py`, `config.py`, and
  `runtime.py` implement and test the frozen mathematical operator contract;
- `encoders.py` and `mechanism.py` implement the local sequence/graph geometry
  bridge that passed P1B.

There is no assembled production DTA pipeline. The previous support encoder,
QPMA, arbitrary 28-dimensional biological state, and end-to-end assembly failed
to establish protein-specific incremental affinity information and were removed.
Their conclusions are preserved in `history.md` and the root evidence-triage
document.

The current integration level is therefore **interface-compatible, not deeply
bio-mechanistically integrated**. The theory accepts an observable bounded
statistic `z(S,Q,gamma)`, but no biological statistic has yet passed the gate
required to enter that coordinate system. Terminal synthetic, structural and
source-affinity evidence remains isolated in `research/e0_identifiability/`.
