# MMP-X F0 correction: exact protein-accession firewall

Date: 2026-07-27

## Reason for correction

The original label-blind F0 counted independent
`(transformation, KLIFS family)` units. Before MMP-X1 could read numerical labels, the downstream
dual-cold contract required one stricter check: a cross-source repetition must not be created by
measuring the same exact protein accession in two sources.

This is a firewall correction, not an after-result threshold change. No numerical activity,
inhibition or affinity value was read.

## Corrected support

- original cross-source repeated transformation--family units: **39**;
- units for which every source pair is exact-accession-disjoint: **8**;
- strict units span only **8 transformations and 5 KLIFS families**;
- 31/39 nominal repetitions reuse an accession in at least one source pair;
- disjoint source-pair counts: Papyrus--Reinecke 8, KIRHub--Reinecke 1,
  KIRHub--Papyrus 0.

## Verdict

**`MMPX_F0_ACCESSION_FIREWALL_INSUFFICIENT_STOP`.**

The corrected graph cannot support a protein-conditioned directional-label audit. MMP-X1 is not
authorized, and no numerical label was read. The original F0 result remains a valid statement about
chemical transformation coverage, but not about protein-cold independent replication.

This failure is a data-overlap failure, not a rejection of local chemical-edit primitives in
principle. It means the current three-source kinase corpus cannot prove them without exact-target
reuse.

Firewall status: current-run confirmation labels unread; historical confirmation remains
quarantined; `sealed_test_consumed=false`.

