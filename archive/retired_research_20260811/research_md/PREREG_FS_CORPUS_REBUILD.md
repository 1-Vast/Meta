# Few-Shot Corpus Rebuild Preregistration

## Scope and authority

This is the O1 follow-up authorized by `task.md` and `project_state.json`.
It is a label-blind structural and power audit only.  It does not read affinity
values, train a model, construct a split, or open confirmation data.

## Estimand and canonical unit

The estimand is quantitative Ki affinity for query ligands of an unseen target
after `k in {1, 2, 3, 5}` distinct support ligands.  The upstream input is the
BindingDB Articles 202608 metadata projection.  A candidate row is admitted
only when it is single-chain and its endpoint list is exactly `["Ki"]`.
Rows whose ligand cannot produce a Murcko scaffold are excluded and counted
before canonicalization; they cannot take part in a scaffold-closed estimand.

The canonical measurement cell is `(target_sequence_sha256, ligand_inchikey)`.
Repeated metadata observations of that pair are merged without reading a label;
their document and assay identifiers are retained as sets.  A `k=5` support set
therefore always means five distinct ligands.  The future numeric aggregation
policy remains deliberately unopened until a label-authorized reliability
protocol is registered.

## Dependency closure

The graph has target sequences as vertices.  It unions targets sharing any
document, any RDKit Murcko scaffold, or a protein pair with local BLOSUM62
identity at least 0.40 (`parasail` Smith-Waterman, gap open 10, extend 1,
matches/alignment length).  MMseqs2 retrieves candidates at 30% identity before
the exact local-alignment confirmation.  This is the same strict family/document/scaffold
closure requested for the O1 reconstruction.  No affinity statistic, model
output, or observed task effect participates in admission or closure.

## Gates

FS-C0 reports total targets/cells, `k=1/2/3/5` eligibility, document depth,
within-target scaffold-disjoint feasibility, component target/cell depths and
largest-component shares.  It fails closed when fewer than five components can
contribute a `k=5` eligible target.

If FS-C0 passes, the split is frozen before FS-C1: the largest dependency
component is source and every remaining complete component is evaluation.  This
rule depends only on closure component target depth (with a lexical component
ID tie-break), not labels or model output.  FS-C1 reports the resulting
target-level evaluation supply.  Its unchanged success requirements are at
least 30 evaluation targets at `k=5`, at least 100 source targets at `k=5`, and
`MDE_d = (1.645 + 0.842) / sqrt(N_eval) <= 0.600`.  Passing FS-C1 does not
authorize training.

## Stop conditions

`FEWSHOT_CORPUS_DEPENDENCY_NOT_IDENTIFIABLE` stops the route at FS-C0.
`FEWSHOT_EVALUATION_POWER_NOT_IDENTIFIABLE` stops it at FS-C1.  A PASS result
is `FEWSHOT_CORPUS_STRUCTURALLY_AND_POWER_FEASIBLE_UNSPLIT`.
