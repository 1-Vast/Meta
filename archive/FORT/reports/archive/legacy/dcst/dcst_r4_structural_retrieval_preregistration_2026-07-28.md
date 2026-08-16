# DCST-R4 bidirectional structural-retrieval preregistration

Date: 2026-07-28  
Status: frozen before R4 implementation and training

## Candidate mechanism

R4 retains the exact-target, SIFTS-aligned R3 architecture and adds
**bidirectional within-target structural retrieval**. For a target episode
with at least two privileged ligands:

1. center predicted and observed `32 × 8` distributions across ligands;
2. L2-normalize each ligand's centered map;
3. form the complete predicted-to-observed cosine-similarity matrix;
4. apply cross-entropy with the true ligand pairing on the diagonal;
5. apply the symmetric observed-to-predicted cross-entropy.

The temperature is fixed at `0.10` and the loss weight at `1.00`. Identical or
zero-variation rows contribute no invented identity signal. This objective
asks whether each ligand's structural signature is distinguishable from the
other ligands for the same unseen protein; it directly targets the
wrong-ligand failure that an aggregate episode cosine can hide.

## Frozen controls

All R3 data, exact target keys, SIFTS labels, singleton absolute loss,
architecture, optimizer, seed 1729, 4,000-step budget, `NoPriv`, spectral
bands, mechanism thresholds, downstream stop rule, and confirmation policy
remain unchanged. The existing aggregate centered, absolute, and destruction
losses remain present.

R4 must pass the same complete source-only gate before any new ChEMBL
affinity-label load.

