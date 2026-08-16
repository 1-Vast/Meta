# DCST-R2 singleton structure-sampling amendment

Date: 2026-07-28  
Status: frozen before the 4,000-step R2 source run

## Trigger

Exact-accession grouping leaves 1,663 of 2,106 admissible Stage-1 train rows
inside targets with at least two affinity ligands. Of the remaining singleton
targets, 430 rows have a valid SIFTS-aligned structural interaction map. The
existing affinity-episode sampler would never expose those high-quality
structural labels.

## Frozen sampling correction

Each of the existing 4,000 source interaction steps retains one exact-target
affinity episode with at least two ligands and all registered losses. In
addition, sample one target-singleton row with valid privileged supervision
and add only:

```text
0.25 * absolute 32 × 8 joint interaction cross-entropy
```

A singleton cannot identify within-target ligand deltas, ranking,
orthogonality, or ligand derangement; those losses are therefore not invented
for it. The matched `DCST-NoPriv` arm retains the same affinity rows and
steps, but receives no singleton structural update because the additional
information is precisely the privileged intervention being tested.

No seed, optimizer, affinity step count, certificate, mechanism threshold,
downstream arm, or confirmation policy changes.

