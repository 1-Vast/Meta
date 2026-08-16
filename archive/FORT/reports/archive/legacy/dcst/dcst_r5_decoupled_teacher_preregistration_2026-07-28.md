# DCST-R5 decoupled privileged-teacher preregistration

Date: 2026-07-28  
Status: frozen before R5 implementation and training

## Candidate innovation

R5 implements the user's two-stage idea as an explicitly decoupled transfer
operator:

### Stage 1A: high-quality structural teacher

Train the exact-target R4 encoder and structural head for 4,000 steps using
only SIFTS-aligned privileged objectives. Multi-ligand targets use absolute,
centered, counterfactual, and bidirectional-retrieval losses; singleton targets
use only absolute structural cross-entropy. No affinity loss updates this
teacher.

### Stage 1B: source affinity readout

Freeze the teacher's target and ligand encoders, reset the bilinear affinity
matrix to the exact null `theta=0`, and train only `theta` for 4,000 steps on
firewalled PLINDER affinity residuals. The structural head cannot be
overwritten by the affinity objective.

### Stage 2: certified downstream transfer

Only held-source target- and ligand-destruction-certified spectral bands of
the Stage-1B matrix may enter the existing frozen transfer bridge. Stage 2
remains unauthorized until the source-only gate passes.

## Matched source controls

All use the same exact target rows, seed, affinity base, optimizer budget and
architecture:

- `PrivTeacher-Frozen`: structurally pretrained frozen encoders; train only
  `theta`;
- `Random-Frozen`: the teacher's exact pre-training initialization, frozen
  without privileged updates; train only `theta`;
- `NoPriv-Trainable`: no privileged labels; train the complete interaction
  branch on affinity as before.

The privileged certificate count must be strictly greater than both controls.

## Frozen source-only gate

1. the Stage-1A teacher passes the existing complete structural mechanism
   gate;
2. `PrivTeacher-Frozen` certifies at least one source affinity band;
3. it certifies strictly more bands than both `Random-Frozen` and
   `NoPriv-Trainable`.

Seed 1729, 4,000 steps per registered phase, exact target/SIFTS data,
certificate thresholds and confirmation policy are unchanged. Failure stops
before any new ChEMBL affinity-label load.

