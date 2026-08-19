# Stage CIIP-0 KiRHub census preregistration (2026-08-19)

Frozen BEFORE any census conclusions are recorded. This stage performs
ONLY a data-access and usability audit of the KiRHub profiling resource
(Nature Biotechnology 2026, s41587-026-03090-8); no model training, no
label use beyond what the audit requires, no GPU.

## Objective

Determine whether KiRHub can serve as the CIIP-1A/1B same-parent
WT/variant positive control. The census must produce, per item, either a
number or an explicit MISSING/UNKNOWN mark — the paper's abstract-level
totals are never substituted for usable counts:

1. usable WT-variant pair count (single-point-mutant pairs only;
   fusions and multi-mutants counted separately and excluded from the
   single-mutant estimand);
2. identical-ligand count per pair;
3. parent, mutation, fusion counts;
4. construct/substrate/cofactor/ATP-condition completeness rate;
5. duplicate and saturation fraction (censoring at 0/100);
6. endpoint direction/units/comparability (% inhibition, bounded);
7. parent/pocket-group connectivity;
8. whether complete held-out parent folds are constructible;
9. per-parent ligand coverage and effective centered-effect variance;
10. data availability: repository, license, direct download URL.

## Stop rules (frozen)

- If the raw data are not legally/reachably obtainable from this
  environment, STOP model training for CIIP-1A/1B on KiRHub; record the
  blocker with evidence (URLs, dates) and the alternative local sources:
  Davis Kd panel (local MOESM files), Anastassiadis 2011 (local),
  Duong-Ly (local), Saifudeen Q3 census.
- The census does not authorize any training on any dataset.
- Phase ordering: the Q2d terminal archival (stageX Q2d-1e + span-param
  diagnostic) completes independently; this census is read-only and does
  not depend on it, but no CIIP training starts before Q2d archival.

## Evidence state marks

[V] verified against a primary source (URL recorded); [I] interpretation;
[U] unverified claim from the literature review or secondary sources.
