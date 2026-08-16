# DCST-R8 frozen protein-language atlas decision

Date: 2026-07-28  
Decision: `STOP_R8_PRESERVE_ATLAS_DECOUPLE_SEGMENT_MAP`

## Result

The frozen atlas was built as registered:

- 767 unique firewalled source-train exact targets;
- 24,544 normalized ESM segment rows;
- all eight clusters nonempty;
- atlas SHA-256
  `018ed822858f79fe7d64d75acda7b0491a590aca4b78ad6067656963f5ab0253`;
- zero affinity, structural, source-development, or downstream input.

`dcst_r8_stage1_seed1729.json` produced `1/4` privileged FPLA bands versus
`0/4` FPLA-NoPriv. The privileged-specific affinity criterion therefore
passed. The complete source gate failed because the segment-level structural
mechanism did not pass:

- true centered alignment: `0.01678`;
- target-destroyed: `-0.00317`;
- ligand-destroyed: `-0.02555`;
- target margin: `0.01995`;
- ligand margin: `0.04233`.

The one active band had true utility `0.11212` and certificate confidence
`0.07234`; all no-privileged confidences were zero.

Wall time was `225.384 s`; peak allocated CUDA memory was `941.1 MiB`. No new
downstream affinity label was loaded.

## Diagnosis

R8 fixed role identity and recovered a privileged-specific atlas direction,
but its atlas affinity loss was still allowed to update the upstream
segment-interaction map. That map lost both registered destruction margins.

R6 provides the complementary fact under identical data and budget: joint
segment-map training passed the structural mechanism and certified 2/4
privileged bands, but absolute segment energies failed downstream transport.
Together the results support a more exact two-step interface:

1. learn and freeze the R6 segment interaction measure that is already
   structurally identifiable;
2. fit only the frozen-atlas energy matrix on top of that measure.

The atlas is retained; joint upstream updating is rejected.

