# R0-B Amendment 01: Full frozen-P1B exposure closure

Status: frozen on 2026-08-11 before exact-geometry generation, ceiling scoring,
or any R0-B fit. No affinity value, exact distance label, model prediction, or
R0-B Gate result was read in making this amendment.

This amendment changes only the label/geometry-free development panel in
`PREREG_R0B_GOVERNED_EXACT_DISTANCE.md`. All mathematical arms, loss functions,
thresholds, seeds, stopping rules and authorization boundaries remain fixed.

## Reason

Panel v1 filtered validation and heldout chemistry against the representatives
selected for residual fitting. That is insufficient because the frozen P1B
prior was trained on every P1B train record and selected using every P1B
validation record. A later audit found 56/241 v1 validation records sharing a
CCD, exact connectivity or Murcko scaffold with the full P1B train exposure,
and 26/170 v1 heldout records sharing one with the full P1B train/validation
exposure. Panel v1 is retired from execution and remains evidence only.

## Replacement panel

The immutable upstream input remains:

`dataset/processed/open_structures/pilot20k_homology_split_v2/complexes.jsonl`

SHA256:
`45907b45b590c6ec27242fc07028444133a5f562f79eff9ba5951cb0b09fae1a`.

The amended label/geometry-free output is:

`dataset/processed/correspondence_router/r0b_governed_panel_v3/panel.jsonl`

SHA256:
`a1f3d29a3b5d876f81a23819f82ac1fa07d681ea500ddb1dbc0f81cc95d89a65`.

Companion hashes:

- `exclusions.jsonl`:
  `b995abff78e879226313f5498860144b70c16aee3bc8efc7a735211b664320c5`;
- `panel_audit.json`:
  `1961d74bb35e6e517c487bac1d5478ca2e5523976a6d7ce8616894a75523d23d`.

Selection still chooses one train representative per registered homology group
with namespace `R0B-TRAIN-REP-v1`. The only amendment is the correct exposure
reference:

1. validation is filtered against every frozen-P1B train record;
2. heldout-A is filtered against every frozen-P1B train and validation record;
3. a dependency match is any shared CCD hash, exact connectivity hash, or
   nonempty Murcko scaffold;
4. the original governed protein-homology split remains unchanged.

Frozen counts are 2,516 train records/components, 185 validation records in 54
components, and 144 heldout-A records in 53 components. The largest heldout-A
component is 17/144 (`11.8056%`), below the fixed 20% limit and above the fixed
30-component minimum. All 2,845 selected records have exact protein mapping and
an admitted exact ligand graph.

## Remaining fail-closed conditions

This amendment does not authorize training. The runner must first complete an
immutable exact-protein bank, exact atom-order mapping, exact residue-by-atom
geometry sidecar, distance-bin census, additive ceiling, MDE/power audit and N3
movable-residue audit. Any shortfall records `R0B_NOT_RUN_FAIL_CLOSED` without
fitting. A development PASS still authorizes only a fresh structural
confirmation cohort, never affinity R1 or V1 integration.
