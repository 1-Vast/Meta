# E-AFF-X0-FEAS Unit Feasibility Audit

Status: registered before computing any ceiling, packing or verdict. X0-FEAS is
an audit of a completed stage's estimand, not a new affinity stage and not a new
Gate. It reads only label-blind design metadata.

## Question

The frozen X0 requirement is `245` effective components per endpoint, where an
effective component is a connected component of the panel / D1-closure graph.
Before spending acquisition and governance effort on a new crossed corpus,
X0-FEAS asks a narrower prior question:

> Is `245` attainable under that unit definition at all?

If the requirement exceeds the structural ceiling of the unit, then the X0
`STOP_SOURCE_INTERACTION_UNDERDETERMINED` verdict was fixed before the query ran
and carries no information about how crossed the source actually is. That
distinction changes which repair is correct.

## Inputs And Label Firewall

- `e0_input_v1/rows.label_blind.jsonl` for closure components, proteins,
  ligand connectivity and endpoint family.
- `artifacts/eaff_x0_v1/{cells,dependency_components,report,manifest}.jsonl|json`
  for the already-published label-blind panel geometry.
- `energy_pilot_v1/task_manifest.json` for D0 task depth, endpoint family,
  documents and E0-Core eligibility, so the ceiling can be evaluated on the full
  governed corpus rather than only on the filtered E0 population.
- `energy_pilot_v1_governance/homology_assignments.jsonl` for the 40% identity
  components and the DAVIS-protected exclusions.
- No SQLite query, no ligand identity outside the published X0 geometry, no
  affinity field, no DAVIS or recipient read, no training.
- The runner fails closed if any input row carries a value-like field, and
  admits D0 task-manifest fields only from an explicit whitelist.

## Quantities

1. `closure_component_universe`: distinct closure components in the governed
   corpus. This is the hard upper bound on effective components for every
   endpoint simultaneously.
2. `closure_components_with_endpoint_data`: the per-endpoint ceiling, since a
   component with no rows for that endpoint can never host one of its
   rectangles.
3. `closure_components_rectangle_capable_panel_free`: components where two
   proteins share at least two ligands *anywhere* in the component, ignoring
   document, protocol and context entirely. This is the ceiling that survives
   even if the panel key were abolished.
4. `design_disjoint_rectangles_greedy`: a greedy target- and ligand-disjoint
   rectangle packing per panel. It reports what remains after removing target
   and ligand repetition, and is recorded as a candidate alternative unit only.
5. `independently_derived_required_n`: an independent re-derivation of the
   frozen `245` from its stated one-sided chi-square design.
6. `governed_corpus_population_sweep`: the closure recomputed over the full D0
   governed corpus and over task-depth populations that the E0-Core `>=20`
   compound contract excluded. For each population it reports closure components
   and, per endpoint, components holding at least two distinct proteins. A
   component holding one protein for an endpoint cannot host a rectangle for
   that endpoint, so this is a hard ceiling on effective components at that
   population. The sweep also records whether the ceiling *falls* as the corpus
   grows, which is the signature of a unit whose count is destroyed by the same
   documents that create crossing.

## Structural Claim Under Test

A rectangle requires two distinct proteins measured inside one panel, and a
panel key contains exactly one `document_chembl_id`. D1 closure unions every
pair of targets that share a document. Therefore both proteins of any rectangle
are necessarily already in one closure component, panel-to-closure union can
never merge two closure components through crossing, and effective components
reduce to closure components that happen to contain a crossed panel.

The audit tests this by counting panels that touch more than one closure
component. The claim predicts zero.

## Verdicts

- `X0_UNIT_REQUIREMENT_ATTAINABLE_<ENDPOINT>`: some governed population reaches
  the requirement for that endpoint, so the X0 stop was informative for it and
  the correct repair is a better census population.
- `X0_UNIT_REQUIREMENT_UNATTAINABLE_BY_CONSTRUCTION`: no governed population
  reaches the requirement for any endpoint.

## Scope Limits

X0-FEAS does not weaken, replace or re-register the frozen `245` requirement,
does not authorize X1, X2, angular work, RFSA, DAVIS or production, and does not
claim that protein-by-ligand interaction exists. An unattainable requirement is
a specification finding only. Any replacement unit or design calculation
requires its own separate registration and an explicit human decision.
