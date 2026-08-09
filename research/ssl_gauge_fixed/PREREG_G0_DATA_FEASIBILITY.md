# G0 preregistration: ensemble and quantitative-anchor data feasibility

Status: **CONDITIONAL, METADATA-ONLY, NO MODEL TRAINING**.

This candidate stage may run only after the external S5-S9 evidence-recovery
audit and without changing the active S5 verdict. It does not read numeric
affinity values and does not authorize A0/A1/A2 training.

## Question

Do public, version-pinned sources contain enough independently governed systems
to compare bound-only, matched-apo/relaxed-unbound, and quantitative
same-target difference information under the required closure rules?

## Allowed sources and roles

- MISATO: bound-complex trajectories (`U3`) and metadata only;
- PLINDER v2: pocket/protein/ligand similarities, matched-series identifiers,
  and linked apo structures (`U1` candidates);
- RCSB: raw coordinates and provenance cross-checks;
- archived BindingDB/ChEMBL metadata: endpoint, assay, target, ligand, document,
  and relation fields, without numeric values;
- OpenFF protein-ligand benchmark: graph and licence census only.

PDBbind-derived labels are not allowed without a separate licence and
provenance decision.

## Required immutable outputs

```text
SOURCE_RELEASE_MANIFEST.json
ENSEMBLE_ROLE_CENSUS.json
APO_HOLO_MAPPING_AUDIT.json
QUANTITATIVE_GRAPH_CENSUS.json
CLOSURE_AND_EXPOSURE_AUDIT.json
COMPUTE_FEASIBILITY.json
G0_DATA_FEASIBILITY_REPORT.md
```

Each output binds release IDs, URLs, licences, archive hashes, extractor commit,
environment, and upstream manifest hashes.

## Hard constraints

1. U1, U2, and U3 are never pooled in the census.
2. A linked apo candidate requires compatible construct, chain mapping, pocket
   mapping, and no bound ligand in the evaluated pocket.
3. Matched-series connected components are indivisible split units.
4. Protein homology, pocket similarity, exact ligand, scaffold, structure,
   trajectory, document, source, and endpoint edges are all included in the
   closure graph.
5. Numeric affinity values, DAVIS, recipient labels, and S9 evaluation labels
   are forbidden.
6. No count is inferred from pair/edge rows as if they were independent; report
   connected components and effective units.
7. PLINDER affinity annotations are not trusted as a label source.

## Terminal verdicts

```text
G0_SOURCE_OR_LICENCE_FAIL_CLOSED
G0_MAPPING_CONTRACT_FAIL_CLOSED
G0_INDEPENDENT_COMPONENTS_UNDERDETERMINED
G0_COMPUTE_BUDGET_INFEASIBLE
G0_DATA_FEASIBLE_FOR_SEPARATE_MODEL_PREREGISTRATION
```

The last verdict authorizes only a new A0/A1 model preregistration. It does not
authorize affinity label access, GPU training, production integration, or a
thermodynamic claim.
