# DCST-R3 substructure-token decision

Date: 2026-07-28  
Decision: `STOP_R3__CAPACITY_WITHOUT_PRIVILEGED_ATTRIBUTION`

## Result

R3 reduced correct-pair structural cross-entropy from R2's 4.9405 to 4.6520
and retained target/ligand destruction in the structural head. It did not
make the affinity transfer path privileged-specific:

```text
R2 certified bands: privileged 2/4, NoPriv 0/4
R3 certified bands: privileged 1/4, NoPriv 1/4
```

The R3 structural centered alignments were approximately `+0.021` correct,
`-0.006` wrong target, and `-0.014` wrong ligand. The registered 0.05 margins
and the no-privileged certificate comparison failed.

## Interpretation

Hashed Morgan tokens increased capacity, but the episode-level centered cosine
objective reduces all ligand distinctions to one aggregate scalar. It does not
require each predicted ligand-specific map to identify its own observed map
among other ligands for the same protein. R3 therefore stops before Stage 2.

