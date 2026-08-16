# OpenMut `OMUT-D1` preregistration

**Frozen:** 2026-07-28, before `research/omut_d1.py` existed and before any `OMUT-D1`
probe or formal run.
**Environment:** `D:\anaconda\envs\drug\python.exe` (Python 3.11.15, torch 2.6.0+cu124,
CUDA 12.4, RTX 4060 Laptop GPU).
**Predecessors:** `OMUT-D0` (`OMUT_D0_SOURCE_FREEZE_COMPLETE`, anchor
`1953cff4c1d7301c51d1ef934c0c5c913f7c154022ad515c53e416bbce8f82f9`) and `OMUT-F0`
(`OMUT_F0_DAVIS_ROLE_FROZEN__PRESERVE_CONFIRMATION`, policy `preserve_confirmation`).
**Authority:** `task.md`, Route A, `OMUT-D1` unlock condition "F0 passes and adequate
base-protein/family/source topology exists", restricted under F0 to **non-Davis
sources**.

## 1. Question and claim boundary

D1 asks one question, on exactly the two sources F0 leaves admissible — the verified
local BindingDB 202607 Articles archive and the bounded ChEMBL variant layer:

> Do these sources contain enough document-independent, base-protein-independent
> WT-to-single-substitution rectangles with shared query ligands to make the mechanism
> question in `task.md` section 1.2 answerable at all?

Per `task.md`, "D1 may report only observation topology, canonical mutation/edit
coverage, graph connectivity/effective rank, independent-unit projections, and
optimistic/projected MDE." It explicitly excludes censoring, numeric relations,
replicate covariance, reliability, observed reordering variance, and empirical MDE,
which belong to `OMUT-I0`.

A pass or fail here is **topology evidence only**. It is not evidence about
substitution semantics, ligand reordering, or dual-cold transfer, and it does not by
itself authorise `OMUT-X0`. `OMUT-X0` additionally requires the evidence-bound registry
(exact construct, endpoint, assay, and document evidence) that D1 does not build.

D1 does not authorise reading any numeric affinity, relation, censor, unit,
temperature, or pH value. It authorises reading identity and provenance fields only:
ligand identity, protein accession, free-text construct annotation, and document
identifiers.

## 2. Sources and the label-free projector

### 2.1 BindingDB (`bindingdb_local_archive`, D0-frozen, SHA-256 verified)

A **full stream** of the single TSV member, not a sample. The projector is a *local
field-selection* blind read: because BindingDB is a variable-width tab-separated text
file, locating column boundaries requires splitting each line, so the full line is
transiently split into a list. The projector then assigns **only** the frozen allowed
indices to named variables; forbidden indices are never read, never named, never
aggregated, and the split line is discarded after each row. This is weaker than the
server-side `only=` restriction used for ChEMBL in D0, and amendment wording says so
explicitly rather than calling it equivalent.

Allowed columns (frozen from the D0 schema, `dataset/public/open_s` archive):

```text
Ligand InChI Key, Ligand SMILES, BindingDB MonomerID, Target Name,
Curation/DataSource, Article DOI, BindingDB Entry DOI, PMID, Institution,
UniProt (SwissProt) Primary ID of Target Chain 1,
Number of Protein Chains in Target (>1 implies a multichain complex)
```

Forbidden columns, asserted absent from every projected row and from every output
structure: `Ki (nM)`, `IC50 (nM)`, `Kd (nM)`, `EC50 (nM)`, `kon (M-1-s-1)`,
`koff (s-1)`, `pH`, `Temp (C)`.

**Mutation-mention detection.** BindingDB's `Target Name` field encodes explicit
single-residue substitutions in bracketed notation, e.g.
`Dimer of Gag-Pol polyprotein [489-587,L512I]`. The frozen detector is the regex

```text
(?<![A-Za-z0-9])([ACDEFGHIKLMNPQRSTVWY])(\d{1,4})([ACDEFGHIKLMNPQRSTVWY])(?![A-Za-z0-9])
```

applied to `Target Name` only. This is a **label-free candidate detector**, not a
verified construct call: it can false-positive on gene-name tokens that happen to match
the pattern, and it does not confirm the mutation against a reference sequence. Every
count below is reported as "candidate", and construct verification against UniProt
sequence is `OMUT-X0` work.

A row whose `Target Name` yields exactly one candidate token is a **single-substitution
row**. A row yielding two or more is a **multi-substitution row** and is excluded from
the primary graph (multi-mutations are sensitivity data, per `task.md` section 3). A
row yielding zero tokens is a **WT/base row** for its accession.

**Grouping key.** Rows are grouped by UniProt accession only, not by the numeric
residue range also present in `Target Name`. This is a stated simplification: two
isoform windows of one accession are treated as the same base protein. Resolving
construct-range identity precisely is `OMUT-X0` work.

**Ligand identity.** `Ligand InChI Key`, already a canonical structural identifier;
rows with an empty key are dropped from the graph.

**Document identity.** `Article DOI`, `BindingDB Entry DOI`, and `PMID` are passed to
`research.omut_f0.provenance_unit` (frozen at F0). A row with none of these keeps its
ligand/mutation data but is not counted toward the independent-document total.

Documents are tracked per `(accession, ligand)` on the WT side and per
`(accession, token, ligand)` on the mutant side, **not per accession as a whole**. A
mutation component's reported document count is the union of documents attached only
to its *shared* ligand measurements. A document covering some other ligand of the same
accession, with no bearing on the shared-ligand rectangle, must not inflate that
rectangle's apparent independent replication — this is exactly the over-counting this
programme's promiscuous-ligand-block audit found and rejected elsewhere
(`document-overlap-binding-constraint`). The reported count is still an upper bound: it
credits a document if it names either the WT or the mutant measurement of a shared
ligand, without confirming both appear in the same document (that confirmation is
`OMUT-X0` work).

### 2.2 ChEMBL variant layer (`chembl_activity_variant_projection`, D0-bounded)

A **bounded sample**, explicitly not a census. Server-side `only=` projection, capped
at 5 pages of 1000 rows (5,000 rows maximum), of activities carrying
`assay_variant_mutation`:

```text
only=activity_id,assay_variant_accession,assay_variant_mutation,target_chembl_id,
molecule_chembl_id,document_chembl_id,standard_type
&assay_variant_mutation__isnull=false
```

`standard_type` (e.g. `Ki`, `IC50`) is an endpoint-type label, not a value, relation,
or unit, and is retained because `task.md` requires `Kd`/`Ki` to be analysed
separately. The runner asserts the returned rows carry no field outside this list,
reusing the field/schema assertions frozen at D0.

**Amendment A1 (2026-07-28, after two execution attempts aborted on the firewall and
before any result was written or accepted).** ChEMBL's `activity` resource always
echoes an extra field named `type` outside the requested `only=` list. The first
attempt assumed it mirrors `standard_type` and aborted correctly when a row showed
`type="Dissociation rate constant"` against `standard_type="k_off"` — a different
vocabulary entirely, not a copy. `type` is not an affinity value, relation, unit, or
censor field (checked against the D0 token vocabulary), but its meaning is not
resolved here, so A1 does not attempt to interpret or validate it: every stored row is
filtered to exactly `CHEMBL_ONLY_FIELDS` before it is kept, and `type` is dropped
unread. Both aborted attempts wrote no artefact. Nothing else in section 2.2 changes.

This sample reports **variant-side topology only**: distinct targets, distinct
mutation strings, distinct documents (after collapse), and distinct molecules linked
to a variant assay. It does **not** build a WT-vs-mutant shared-ligand graph for
ChEMBL, because that needs a second, larger query for non-variant activities on the
same targets. That comparison is deferred to `OMUT-X0`, where the full evidence-bound
registry is built; D1 states this rather than approximating it.

## 3. Reported topology

For BindingDB, per accession/mutation-token pair with exactly one candidate token:

```text
shared(accession, token) = wt_ligands(accession)  intersect  mutant_ligands(accession, token)
```

Reported, for `k in {4, 8, 16}`: the count of `(accession, token)` pairs with
`|shared| >= k`, the count of distinct accessions among them, and the count of
distinct provenance units covering their contributing rows.

Also reported: repeated substitution classes (`wt_residue -> mut_residue`, position
ignored) and the number of distinct accessions each class appears in among the
`k=4` set; a per-accession document count; and the multi-substitution row count as a
separate, non-graph statistic.

**Comparison-graph connectivity.** An accession-accession graph with an edge when two
accessions share at least one ligand identity anywhere in the streamed data,
restricted to accessions with at least one `k=4` component; reported as connected-
component count and largest-component size.

**Effective rank.** The binary incidence matrix of `k=4` mutation components (rows)
against the ligands in their shared sets (columns); reported as the count of singular
values exceeding `1e-6` times the largest singular value.

**Projected MDE.** Using the frozen convention already in this program (Davis
registration, `reports/active/panel_davis_registration.md`; OpenMut crossreview
section "DAVIS-Complete"): assumed paired SD `0.10`, `mde_from_spread` from
`research/dualcold_power.py`, seed `1729`, `n` equal to the `k=4` independent
mutation-component count. This is optimistic and projected, not empirical; no real
paired difference is computed because no affinity value is read.

Explicitly **not** reported, and asserted absent from the output by key name: any
censoring rate, replicate covariance, reliability estimate, observed reordering
variance, or empirical MDE. Those are `OMUT-I0`.

## 4. Gates

All nine are evaluated; all nine must pass.

| gate | condition |
| --- | --- |
| `D1_F0_COMPLETE` | `reports/active/omut_f0.json` exists with verdict `OMUT_F0_DAVIS_ROLE_FROZEN__PRESERVE_CONFIRMATION` and a matching preregistration hash |
| `D1_PROJECTOR_BLIND` | no forbidden field name or its value ever appears as a key or a value fingerprint in any aggregation structure or in the output; asserted at runtime and by static test |
| `D1_BINDINGDB_FULL_STREAM` | every line of the local TSV member is scanned (`rows_scanned` equals the archive's true row count, not a sample) |
| `D1_CHEMBL_BOUNDED_SAMPLE` | ChEMBL rows fetched `<= 5000`, fetched only through the frozen `only=` projection |
| `D1_DOI_COLLAPSE_APPLIED` | every reported document count is produced through `research.omut_f0.provenance_unit`, not through a raw source field |
| `D1_TOPOLOGY_REPORTED` | the `k in {4,8,16}` counts, substitution-class table, and per-accession document counts are all present, including the zero case |
| `D1_GRAPH_METRICS_COMPUTED` | connected-component count, largest-component size, and effective rank are computed and internally consistent (largest component `<=` total nodes, effective rank `<=` min(rows, cols)) |
| `D1_MDE_PROJECTED` | the MDE is computed with the frozen SD, seed, and `n`, and is deterministic across two calls |
| `D1_NO_I0_SCOPE_LEAK` | the output contains no key naming censoring, replicate covariance, reliability, observed variance, or empirical MDE |

## 5. Verdicts

```text
OMUT_D1_TOPOLOGY_ADEQUATE       all nine gates pass and the k=4 independent
                                 mutation-component count meets or exceeds 25
OMUT_D1_TOPOLOGY_INADEQUATE     all nine gates pass but the k=4 count is below 25
OMUT_D1_INCOMPLETE_STOP         any gate fails
```

The threshold of 25 independent components is not invented here: it is the low end
of the "~25-30 document-independent components" bound this program has used
repeatedly (`document-overlap-binding-constraint`, CROSSDOC, REWIRE) as the point
below which a component-blocked confidence bound cannot plausibly clear a `0.03`
material effect. Falling below it is not a failure of this stage; it is the expected,
informative outcome the stage is designed to detect.

## 6. Stop rules

- A count below threshold is not rescued by lowering the shared-ligand requirement,
  by relaxing the mutation-token regex, by merging BindingDB and ChEMBL components
  before independence is established, or by treating multi-substitution rows as
  single.
- The mutation-token detector is a candidate detector. `OMUT_D1_TOPOLOGY_ADEQUATE`
  does not certify that any candidate is a real single-amino-acid substitution;
  `OMUT-X0` must verify against reference sequence before any coordinate or estimator
  work.
- No GPU computation is performed at D1.
- Grouping by accession only (not construct range) can only ever inflate the reported
  topology, never deflate it; a pass is not weakened by this simplification, but a
  fail is not strengthened by it either — the true, range-resolved topology is `<=`
  what D1 reports.

## 7. Deliverables

- `research/omut_d1.py` — the projector, the graph/MDE computation, and the gates.
- `tests/test_omut_d1.py` — projector-blindness, regex, graph, MDE, and gate tests.
- `reports/active/omut_d1.json` — machine-readable topology result.
- `reports/active/omut_d1_decision.md` — verdict and what it does or does not unlock.
