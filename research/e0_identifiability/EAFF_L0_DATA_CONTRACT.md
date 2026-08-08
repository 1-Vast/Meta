# E-AFF-L0 Data And Estimand Contract

Status: registered before the identifiability statistics were computed and
before any arm was fitted. Frozen here; violations fail closed.

## Source

- Release-pinned ChEMBL37 SQLite only, verified against `release_manifest.json`.
- Governed E0-Core label-blind row index for activity selection, proteins,
  ligand connectivity, closure components and outer folds.
- `Ki` and `Kd` are separate endpoints and are never pooled.
- Research-only. No code enters `model/` or normal `scripts/`.

## The Assay Stratum Must Be Target-Independent

The corpus field `assay_context_sha256` **cannot** be used as the L0 stratum.
Its payload includes `component_accession`, `variant_id` and
`variant_mutation` ([canonicalize.py:72](../../scripts/source_affinity/canonicalize.py:72)),
so conditioning on it would be a target-identifier shortcut, which this contract
prohibits.

L0 therefore defines the stratum from the target-independent context keys
already used and audited by X0:

```text
kappa_raw = (assay_organism, bao_format, cell_id, tissue_id,
             subcellular_fraction, relationship_type, assay_parameters)
```

Document identifiers, assay identifiers, target accessions, variant fields and
task identifiers are excluded from the stratum. Metadata is obtained through the
existing audited value-free projection `x0_metadata.sql`, which selects no
affinity field.

## Prohibited

- DAVIS data and recipient labels;
- target identifiers of any kind;
- assay-identifier embeddings;
- full task identifiers;
- deranged-protein examples in training;
- a support nuisance available to one arm and not to its control;
- post-hoc selection of thresholds, models, anchors, coordinates or margins.

## Required Governance

- homology-document closure components define the outer OOF folds;
- every learned baseline, prior and calibration is cross-fitted by closure
  component, never fitted on its own evaluation fold;
- ligand connectivity keys and Murcko scaffolds are audited, and fit/evaluation
  ligand partitions inside a task are scaffold-disjoint;
- primary inference is closure-component macro with a closure-component
  bootstrap for 95% intervals;
- evaluation proteins lie in closure components disjoint from every component
  used for fitting, anchor calibration or population-band estimation.

## Location Estimand Identifiability, Frozen Criterion

The location estimand is only meaningful if a protein's affinity level can be
distinguished from the offset of the assay stratum it was measured in. Build the
bipartite design graph with protein nodes and stratum nodes, one edge per
observed (protein, stratum) cell, separately per endpoint. The frozen
requirements are:

| Code | Requirement | Threshold |
|---|---|---|
| C1 | strata containing at least two distinct proteins | `>= 30` |
| C2 | proteins appearing in at least two distinct strata | `>= 30` |
| C3 | protein fraction in the largest connected component | `>= 0.50` |

C1 supplies within-stratum protein contrasts, so a protein effect is not forced
to absorb its stratum's offset. C2 supplies links between strata. C3 ensures the
protein levels that survive are mutually comparable rather than fragmented into
incomparable islands.

An endpoint is admitted only if C1, C2 and C3 all hold for it. If no endpoint is
admitted, L0 stops with

```text
L0_NOT_RUN_LOCATION_ESTIMAND_NOT_IDENTIFIED
```

If exactly one endpoint is admitted, L0 runs on that endpoint alone and reports
the other as not identified. Thresholds are frozen here and may not be revised
after the statistics are seen.

## Label Access

L0 reads governed ChEMBL37 source `Ki`/`Kd` values through the existing
canonical row store, for the selected activity IDs only. The exact label fields
accessed are recorded in the run manifest. DAVIS and recipient reads remain
zero. The identifiability check in this contract is computed from design
metadata only and reads no affinity value.
