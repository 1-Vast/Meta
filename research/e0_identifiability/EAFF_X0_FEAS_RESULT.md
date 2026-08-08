# E-AFF-X0-FEAS Unit Feasibility Audit Result

## Verdict

```text
X0_UNIT_REQUIREMENT_UNATTAINABLE_BY_CONSTRUCTION
```

The frozen X0 requirement of `245` effective components per endpoint exceeds the
structural ceiling of the X0 independence unit at every governed population of
ChEMBL37, including the full corpus. X0's
`STOP_SOURCE_INTERACTION_UNDERDETERMINED` verdict was therefore determined
before the census query ran, and does not measure how crossed ChEMBL37 is.

## Label Firewall

X0-FEAS read `152,737` label-blind rows, the published X0 label-blind panel
geometry, the D0 task manifest under a field whitelist, and the D1 homology
assignments. It opened no SQLite connection, selected zero affinity fields,
trained nothing, and performed zero DAVIS or recipient reads.

## 1. The Unit Cannot Be Created By Crossing

A rectangle requires two distinct proteins inside one panel; a panel key carries
exactly one `document_chembl_id`; D1 closure unions every pair of targets that
share a document. Both proteins of every rectangle are therefore already in one
closure component before the panel graph is built.

Predicted panels touching more than one closure component: `0`.
Observed: `0`.

Panel-to-closure union can never merge components through crossing. Effective X0
components collapse exactly onto closure components that happen to contain a
crossed panel, and are bounded above by the closure-component universe. Every X0
dependency component in the published artifact contains exactly one closure
component, which is the same fact seen from the output side.

## 2. The Requirement Exceeds The Ceiling On The X0 Population

Closure-component universe of the E0 input corpus: `245`.
Frozen requirement: `245` **per endpoint**.

| Endpoint | Components with endpoint data | Panel-free rectangle-capable | X0 observed | Required |
|---|---:|---:|---:|---:|
| Ki | 202 | 57 | 36 | 245 |
| Kd | 72 | 12 | 12 | 245 |

Satisfying the requirement for one endpoint would demand that every component in
the corpus carry that endpoint *and* host a crossed panel. Only `202` components
carry any Ki row and only `72` carry any Kd row, so neither endpoint can reach
`245` at any query, threshold or protocol setting.

The panel-free column removes the document, protocol and context key entirely
and asks only whether two proteins in the component share two ligands anywhere.
Even that maximally permissive relaxation gives `57` and `12`. For Kd the
observed `12` already equals the panel-free ceiling, so ChEMBL Kd crossing is
fully exploited and no better query exists.

## 3. The Ceiling Does Not Improve On The Full Governed Corpus

X0 ran on the E0-Core population, which the D0 task contract restricted to tasks
with at least `20` exact compounds. That contract exists so task-local ranking is
estimable; a rectangle needs only two shared ligands. Removing it multiplies the
corpus roughly tenfold, so the ceiling was recomputed across populations.

`components_with_two_or_more_proteins` is a hard upper bound on effective
components: a component holding a single protein for an endpoint cannot host a
rectangle for it.

| Population | Tasks | Proteins | Closure components | Ki ceiling | Kd ceiling |
|---|---:|---:|---:|---:|---:|
| E0-Core eligible (X0 population) | 3,817 | 697 | 253 | 75 | 14 |
| full D0 governed corpus | 37,783 | 4,787 | 459 | 72 | 56 |
| at least 5 compounds per task | 14,589 | 1,492 | 434 | 92 | 51 |
| at least 2 compounds per task | 23,961 | 3,006 | 494 | 97 | 54 |

Best ceiling over every population: Ki `97`, Kd `56`, against a requirement of
`245`. Opening the full corpus grows tasks by `9.9x` and proteins by `6.9x` while
the Ki ceiling *falls* from `97` to `72`, because the added tasks add documents,
and every added document unions more targets into one component. Kd rises from
`14` to `56` over the same range, so the effect is not universal; the reliable
statement is that the ceiling is governed by document disjointness rather than by
data volume, and is not monotone in corpus size.

## 4. Frozen Power Derivation Reproduced

The `245` was independently re-derived from its stated design: one-sided
chi-square variance test, variance ratio `1.25` (interaction RMS `0.5` assay
noise SD), `alpha=0.05`, power `0.80`, `df=n-1`, giving `n=245`. The frozen
number is arithmetically correct.

The defect is not the number. It is that a sample-size formula assuming
independent identically distributed units was bound to a unit whose universe on
this corpus is `245` in total. The two `245` values are unrelated quantities that
coincide numerically.

## 5. Consequence For The Named Next Step

The X0 result named acquisition of a genuinely crossed source/selectivity corpus
as the next admissible action. Under this unit that action does not repair the
census.

Reaching `245` units requires `245` protein groups that are simultaneously
pairwise below 40% identity **and** document-disjoint, each internally crossed by
at least two ligands. Internal crossing comes from screening several proteins
against shared compounds in one study, which is exactly what a selectivity or
profiling corpus is, and any such study unions all of its targets into one
component. A panel spanning several families merges those families. The unit
reaches its maximum on many small, mutually unrelated two-protein studies and is
reduced by a single kinome-wide or receptor-panel screen. Section 3 shows this
happening on real data for Ki.

## 6. Candidate Alternative Unit, Recorded Only

A target- and ligand-disjoint greedy rectangle packing was computed per panel. It
removes the target-repeat and ligand-repeat pseudoreplication that inflated the
nominal counts, while keeping crossing.

| Endpoint | Design-disjoint rectangles | Clusters | Largest cluster share |
|---|---:|---:|---:|
| Ki | 705 | 36 | 0.444 |
| Kd | 62 | 12 | 0.355 |

These are not certified independent units: packed rectangles share an assay run,
and the double difference cancels only effects additive in target and in ligand
within a panel. The numbers are recorded so that a corrected design calculation
has a concrete starting design, not as a replacement threshold.

## Scope

X0-FEAS is an estimand audit. It does not re-register the frozen requirement,
does not change any Gate, does not authorize X1, X2, angular work, RFSA, DAVIS or
production, and is not evidence that protein-by-ligand affinity interaction
exists. It establishes one fact: the current X0 stop is a specification
consequence rather than a measurement of the source, so no corpus acquisition
should be funded against it until the unit and its design calculation are
re-registered together.
