# DCST-R2 entity-aligned decision

Date: 2026-07-28  
Decision: `STOP_R2__LIGAND_GLOBAL_BOTTLENECK`

## Result

The exact-target, SIFTS-aligned seed-1729 source run used 2,106 train and 220
development rows. The privileged model certified two of four affinity bands;
the matched exact-target `DCST-NoPriv` model certified zero.

The structural mechanism probe was:

| Probe | Cross-entropy | Centered alignment |
| --- | ---: | ---: |
| correct target and ligand | 4.9405 | +0.0273 |
| wrong target | 6.9982 | +0.0044 |
| wrong ligand | 4.9659 | -0.0144 |
| uniform | 5.5452 | — |

Entity alignment fixed the target-side shortcut: wrong-target absolute loss
increased by 2.058 and its centered alignment fell. The registered centered
margin was nevertheless only 0.0229 for target destruction and 0.0417 for
ligand destruction, below 0.05. R2 therefore stops before Stage 2.

## Interpretation

The residual failure is consistent with the ligand architecture. One global
Morgan-plus-descriptor vector is copied to all protein segments. It can
identify a ligand's broad interaction distribution but cannot express which
local chemical environment should interact with which protein segment.

PLINDER's annotation table does not expose a reliable ligand-atom index for
each interaction, so an atom-contact ground truth must not be fabricated.
The admissible successor is R3 substructure-token cross-attention, which
retains each active Morgan environment as a token and supervises its
protein-segment aggregation with the already valid SIFTS-aligned map.

