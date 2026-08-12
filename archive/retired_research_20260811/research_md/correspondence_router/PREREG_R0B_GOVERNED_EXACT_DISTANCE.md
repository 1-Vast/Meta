# R0-B Development Preregistration: Governed Exact Distance Residual

Status: frozen after the label/geometry-free panel selection and before exact
geometry scoring, exact ESM caching, or any R0-B model fit.

R0 on C0 stopped before training because the strict, exact-mapping cohort had a
38.05% largest protein closure component, above the registered 20% limit. R0-B
does not change the scientific axis or the Gates in
`PREREG_R0_EXACT_DISTANCE_RESIDUAL.md`; it replaces only the statistically
invalid development population.

## 1. Frozen panel

Input: governed P1B homology split v2, SHA256
`45907b45b590c6ec27242fc07028444133a5f562f79eff9ba5951cb0b09fae1a`.

Label/geometry-free selection output:
`dataset/processed/correspondence_router/r0b_governed_panel_v1/panel.jsonl`,
SHA256 `de1b60f014a82b8f9d8c475aa765ad556b8a9dbbcb45ec603d1140917aae03ed`.

Selection is fixed as follows:

1. require 100% protein coordinate-to-sequence mapping and a nonempty exact
   canonical ligand graph;
2. choose exactly one train record per registered 40%-identity homology group
   using SHA256 namespace `R0B-TRAIN-REP-v1`, without reading later splits;
3. remove validation records whose exact graph or Murcko scaffold occurs in
   selected train;
4. remove heldout-A records whose exact graph or scaffold occurs in selected
   train or retained validation;
5. preserve the original homology-group train/val/test assignment.

Frozen counts are 2,516 train records/components, 241 validation records in 66
components, and 170 heldout-A records in 62 components. The largest heldout-A
component has 17/170 records (10.0%). No affinity or geometry value was read to
make this selection.

The frozen P1B checkpoint used train to fit, validation to select, and its test
split has already been consumed by the P1B Gate. Consequently heldout-A is an
R0-B development set, not an independent confirmation set. A PASS cannot
authorize affinity R1 without a future fresh structural confirmation.

## 2. Exact geometry and feature contract

- Reconstruct each full ligand-heavy-atom by exact-residue minimum-distance
  matrix from the immutable mmCIF and CCD paths in the governed record.
- Require the mmCIF atom-name set to equal the CCD canonical heavy-atom set and
  the stored ligand graph's `atom_mapping_hash`; quarantine any mismatch.
- Require sequence indices to be exactly `0..L-1`, with no duplicates.
- Use every exact cell; no distance-based candidate selection or negative
  sampling is allowed.
- Generate exact residue states with the same frozen ESM2-t30 revision
  `a695f6045e2e32885fa60af20c13cb35398ce30c` and the frozen P1B
  `bank_proj`. Cache only projected fp16 states plus the recomputed 128-slot
  states. ESM and P1B have zero trainable parameters.
- Generate exact ligand-atom states with the frozen P1B GINE and the existing
  exact CCD-hash ligand bank. No ligand encoder parameter is unfrozen.
- The frozen P1B distance posterior is lifted from slot to exact residue as the
  prior. Contact logits are diagnostic only.

## 3. Arms, fitting, and stopping

The arms and Gates remain N0 PRIOR, N1 SLOT_SHARED, N2 ADDITIVE_EXACT, N3
RES_DERANGE, N4 ATOM_DERANGE, and the full exact bilinear residual as defined in
the parent preregistration. Full, N1 and N2 use the same three seeds
`(1701, 1702, 1703)`, train-record order, number of epochs, early-stopping
metric, and selected checkpoint epoch. Hyperparameters are selected on
component-macro validation RPS once and then frozen across arms.

Before fitting, the exact-geometry runner must finish the parent preregistration
ceiling/power audit. In addition to its conditions, R0-B requires at least 30
scorable heldout components after every atom/mapping quarantine and at least
50% of exact residues to be movable in N3; otherwise verdict
`R0B_NOT_RUN_FAIL_CLOSED`.

The primary score, `delta_star = 0.05*S_prior`, G1/G2/G3, seed-direction rule,
NLL guard, component bootstrap, and terminal verdict names are unchanged. No
threshold may be altered after an arm result is read.

## 4. Authorization boundary

An all-Gate PASS records only
`R0B_DEV_EXACT_DISTANCE_RESIDUAL_IDENTIFIED`. It licenses construction of a
fresh structural confirmation cohort. It does not by itself license affinity
R1, Meta-Section integration, `z`, Q-PMA, CSMO, or production migration.
