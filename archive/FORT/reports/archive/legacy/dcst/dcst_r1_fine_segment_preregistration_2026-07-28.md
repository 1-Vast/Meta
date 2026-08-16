# DCST-R1 fine-segment source preregistration

Date: 2026-07-28  
Status: frozen before generating or evaluating the 32-segment features

## Motivation

`DCST-P0` compressed each protein into eight equal-length ESM-2 means. On the
source development split, 4,000-step privileged training improved absolute
joint-distribution cross-entropy from the uniform value `4.159` to `3.224`,
and wrong-target replacement worsened it to `3.842`. The registered centered
joint gate nevertheless failed because its eight positional bins did not
provide a stable cross-protein coordinate:

```text
true centered alignment          +0.0033
wrong-target centered alignment  +0.0434
wrong-ligand centered alignment  +0.0063
```

For a median source protein of about 403 residues, one token averages roughly
50 residues. `DCST-R1` tests the pre-specified representation diagnosis that
this pooling discards the local sequence information needed to bind PLINDER's
residue-level interaction annotations to the correct target.

## Frozen representation and data

- frozen language model: `facebook/esm2_t33_650M_UR50D`, unchanged;
- target tokens: 32 equal-length residue intervals, mean-pooled from frozen
  residue representations;
- source keys: PLINDER train and development targets only;
- downstream keys, if and only if the source gate passes: ChEMBL train and
  development targets only;
- maximum sequence length and truncation: unchanged at 1,022 residues;
- ligand input: unchanged Morgan-1024 plus ten train-standardized
  descriptors;
- privileged target: a normalized `32 × 8` segment-by-interaction-type map;
- affinity and counterfactual losses, spectral rank, band grouping, source
  firewall, seed and step counts: unchanged from the P0 amendments.

The ESM feature builder may read only target keys, split labels, cached
sequences, and frozen model weights. It must not use affinity values.

## Source-only gate

Run seed 1729 for 4,000 source-base and 4,000 interaction steps with
`--stage1-only --target-segments 32`. All must hold:

1. true centered joint alignment is positive;
2. true alignment exceeds wrong-target alignment by more than `0.05`;
3. true alignment exceeds within-target wrong-ligand alignment by more than
   `0.05`;
4. true absolute joint cross-entropy is below uniform `log(256)`;
5. at least one two-direction source affinity band is certified;
6. the privileged model certifies strictly more bands than `DCST-NoPriv`.

Failure stops R1 before any new ChEMBL affinity-label load. Passing R1
authorizes the unchanged registered Stage-2 arm comparison; it does not
authorize confirmation or sealed-test scoring.

## Claim boundary

R1 is evidence for or against fine-grained structure-privileged source
distillation. It is not evidence that generic ESM pretraining, generic
fine-tuning, or more target tokens alone solve dual-cold affinity prediction.
