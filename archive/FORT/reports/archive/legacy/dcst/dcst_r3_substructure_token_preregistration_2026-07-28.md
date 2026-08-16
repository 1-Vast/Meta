# DCST-R3 substructure-token interaction preregistration

Date: 2026-07-28  
Status: frozen before R3 implementation and training

## Candidate innovation

R3 adds a **substructure-token interaction cross-attention** module to the
entity-aligned two-stage pipeline:

1. every active Morgan-1024 bit is treated as a learned local-environment
   token;
2. the ten ligand descriptors provide a shared ligand context token;
3. each of the 32 exact-UniProt ESM segments attends only to active
   substructure tokens;
4. the resulting segment-specific ligand contexts predict the SIFTS-aligned
   `32 × 8` structural interaction distribution;
5. the affinity ligand state pools the same trained substructure tokens, so
   the source spectral transfer path shares the privileged representation;
6. only held-source target- and ligand-destruction-certified spectral bands
   are exposed to Stage 2.

This is not an atom-contact claim: Morgan bits are hashed local chemical
environments and PLINDER supplies no atom-index target in the local table.
The claim is that retaining local ligand tokens, rather than collapsing them
before target conditioning, makes the available pair-specific structural
signal identifiable.

## Frozen controls and gate

- exact R2 target registry, SIFTS projection, source firewall and 32 protein
  segments;
- seed 1729, 4,000 source-base and 4,000 source-interaction steps;
- identical R3 architecture and affinity rows for `DCST-NoPriv`;
- exact null at interaction matrix `theta=0`;
- unchanged four two-direction spectral bands and certificate scale;
- unchanged structural gate: true centered alignment positive and more than
  0.05 above both wrong target and wrong ligand, with true CE below
  `log(256)`;
- at least one privileged certified band and strictly more bands than
  `NoPriv`.

Failure stops before any new downstream-label load. Passing authorizes the
unchanged Stage-2 comparisons; confirmation remains sealed.

