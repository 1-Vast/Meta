# Stage CIIP-0b Davis panel census preregistration (2026-08-19)

Frozen BEFORE the census computation. Substitute positive-control source
after the KiRHub census returned DATA BLOCKER (stageCIIP_kirhub_census_20260819,
CENSUS.json). Read-only data audit of the LOCAL Davis 2011 supplementary
tables (davis_MOESM3.xls kinase metadata with Mutant YES/NO flags;
davis_MOESM5.xls 442-kinase x 72-inhibitor Kd matrix, nM; provenance:
Nature Biotechnology 2011, Davis et al., supplementary material, kept
under tools/research/stageX_csc_signal/downloads/). No training, no GPU,
no label use beyond auditing.

## The same 10-item checklist as the KiRHub census applies

1. usable WT-variant pair count: (WT row, mutant row) pairs within one
   Entrez Gene Symbol with at least one commonly measured ligand;
   fusion/multi-mutant rows counted separately and EXCLUDED from the
   single-mutant estimand (Mutant == YES with a single mutation tag is
   the admission rule; unparseable mutant strings -> quarantined).
2. identical-ligand count per pair (finite Kd in both rows).
3. parent (gene), mutation, fusion counts.
4. condition completeness: the panel is a single-assay-condition
   competition-binding Kd (fixed ATP per published protocol) -> recorded
   as ONE condition, completeness = measured-cell fraction; no
   per-cell ATP/construct variability table exists in these files ->
   PARTIAL.
5. duplicate and saturation fraction: duplicate (gene, mutant) rows;
   NA fraction and the 10-uM cap semantics recorded (9900 max observed;
   >cap vs not-tested distinction recorded as UNKNOWN where the file
   does not distinguish).
6. endpoint direction/units: Kd [nM], competition binding, larger =
   weaker binding; monotone, but NOT pK/Ki/DTA and never relabeled.
7. parent/pocket-group connectivity: genes with >=1 pair; pairs per gene.
8. held-out parent folds: number of parents with >=2 mutant rows (or
   >=2 pairs) usable for leave-one-parent-out CIIP-1B.
9. per-parent ligand coverage and effective centered-effect variance:
   per pair, variance of log10(Kd_mut/Kd_wt) over common ligands,
   pairs with >=3 common ligands reported.
10. data availability: local files, sizes, and the original SI URLs.

## Stop rules (frozen)

- If usable pairs < 20 or median common ligands < 5 or held-out-parent
  count < 10: record INSUFFICIENT and do not authorize CIIP-1A/1B
  training on Davis alone; alternatives (Anastassiadis, Duong-Ly) get
  the same census before any training.
- The census never authorizes training by itself; CIIP-1A/1B need their
  own preregistrations and the Q2d terminal archival completed first.
