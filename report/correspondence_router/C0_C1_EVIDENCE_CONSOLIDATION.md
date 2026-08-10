# C0/C1 correspondence information audit — evidence consolidation

## Terminal result

```text
C0 corpus and closure ................ ALL GATES PASS
C0 mapping equivalence ............... PASS after amendment 01 (60/60)
C1a exact coupling vs fixed-degree ... FAIL
C1b unit sufficiency ................. PASS
C1c replicate ceiling defined ........ PASS
terminal verdict ..................... EXACT_EDGE_COUPLING_NOT_SUPPORTED_BY_TEACHER
C2 router ............................ NOT AUTHORIZED, NOT TRAINED
```

Audit only. Zero trainable parameters, zero affinity value reads, heldout-A
never referenced, PLINDER not used, heldout-B not created.

## 1. Exposure and the untouched corpus

Every PDB id consumed by any prior stage in any role was unioned into an
exposure registry: the pilot20k governed, ungoverned and homology-split corpora
(P1B), the preflight and preflight_qc100 QC corpora, the MONN development and
additional-PDB edge corpora (S7, B5, S3R, S4R, S5D) and the ssl_b2 independent
structural set. That union is **24,874** PDB ids.

```text
local raw mmCIF entries                    14,169
minus the exposure registry                 2,836   untouched
with a BioLiP2 biologically relevant ligand 2,509
admissible systems after S1-S7              2,039   over 669 PDB entries
after selection rule S8 (CCD -> Murcko)     1,862   scored panel
distinct receptor sequences                   982
distinct ligand CCD ids                       340
exclusions, with reasons, never replaced    5,289
```

Exclusion reasons: covalent link `3,164`, sequence length `1,350`, resolution
`395`, ligand atom count `372`, not X-ray `4`, no contact within `6.0 A` `4`,
plus `177` dropped by S8 for lacking a parseable CCD scaffold.

BioLiP2 was used only to decide which non-polymer entity is a genuine ligand
and, after amendment 01, to supply the receptor sequence string. No numeric
BioLiP2 field and no affinity column was read. PLINDER was not used: the
standing licence audit already forbids it, and that decision was honoured
rather than reopened.

## 2. The mapping correction, and why it matters

The preregistration asserted that P1B's `record["sequence"]` is the mmCIF
`_entity_poly` canonical sequence, so that `sequence_index = label_seq_id - 1`
would reproduce P1B's slot assignment exactly. It required that claim to be
verified against the parasail mapping and to **fail closed** otherwise.

It failed: `23 / 40`. The premise was factually wrong. P1B takes its sequence
from **BioLiP column 20**, the receptor sequence, which is systematically
shorter than the entity sequence — `162` versus `164` for `182l`, `332` versus
`348` for `1a0i` — and is not even a prefix for `1a0q`.

Amendment 01 replaced `M3`/`M4` with the true P1B path — BioLiP receptor
sequence plus the parasail alignment coordinate — and re-keyed system selection
to BioLiP rows with auth-to-label asym resolution, exactly as
`scripts/build_holo_complex_index` does. Under the corrected path the check
passes **60 / 60** on sequence index and on slot assignment.

The C1 execution that had been started under the rejected mapping was stopped,
and its partial output was discarded without being inspected. No Gate,
threshold, margin, seed or verdict was altered. This was a provenance
correction caught by the audit's own fail-closed check, not a rescued Gate.

## 3. C0 admissibility — all gates pass

The union closure over protein identity `>= 0.40`, shared Murcko scaffold or
CCD id, and exact sequence produced **89** components but its largest component
exceeded the `0.25` cap. That is precisely the giant-component failure the
parent instruction warned against assuming away, and the registered
DataSAIL-style two-dimensional fallback was used instead.

| Gate | observed | required | result |
|---|---:|---:|:---:|
| G0a independent inference components | 496 | >= 60 | PASS |
| G0b largest component fraction | 0.0811 | <= 0.25 | PASS |
| G0c minimum detectable effect at 80% power | 0.00453 | <= 0.05 | PASS |

`G0c` used `sigma = 0.04053` estimated from the degree-preserving null arm
only, never from the contrast under test, over `N = 496` components. The panel
is therefore powered to detect an effect an order of magnitude smaller than the
registered margin, so `C1a` is a fair test rather than an underpowered one.

The registered 3-mer containment prefilter was measured rather than relied on:
exact brute-force all-vs-all alignment found `3,037` true identity edges, of
which the prefilter would have missed `0` — recall `1.0`.

## 4. C1 — the decisive measurement

Deconvolution census over the scored panel: `4,702,002` atom-slot units in
candidate slots, `162,276` positive units, `100,563` valid 2x2 checkerboards
read entirely from raw coordinates, multi-contact rate `0.257`, median 3
residues per candidate slot.

Within a positive unit, candidate residues are ranked by their additive contact
degree and scored by exact tied AP, aggregated complex -> closure component.
Rows were never bootstrapped as IID.

| arm | component-macro within-slot AP |
|---|---:|
| empirical | 0.985611 |
| fixed-degree rewire null | 0.953959 |
| atom shuffle | 0.985611 |
| geometry shuffle | 0.993948 |
| empirical, informative units only | 0.983487 |
| rewire, informative units only | 0.950276 |

| Gate | delta | LCB95 | required | result |
|---|---:|---:|---:|:---:|
| C1a empirical - fixed-degree rewire | +0.031652 | +0.029690 | +0.05 | FAIL |
| C1b positive units / checkerboards | 162,276 / 100,563 | — | 10,000 / 1,000 | PASS |
| C1c replicate ceiling defined | jaccard 0.830 | — | positive, stable | PASS |

Complete-edge AP of the rank-1 additive marginal product: `0.753449`.

## 5. What the numbers actually say

**The within-slot deconvolution task is nearly saturated by additive
marginals.** The empirical AP is `0.9856`. Ranking a slot's candidate residues
by nothing more than how many ligand atoms each residue contacts overall
answers the question almost perfectly. The total headroom available to *any*
predictor above that baseline is `1 - 0.9856 = 0.0144` AP.

That single fact is more decisive than the Gate. The registered `+0.05` margin
cannot be reached on this statistic even by an oracle, because there is only
`0.0144` of room above the marginal predictor. A geometry-gated router competing
here would be fighting for at most one and a half AP points that additive
degree information has already taken.

The `+0.031652 [LCB +0.029690]` empirical-minus-rewire effect is real and
tightly bounded away from zero: the true arrangement is measurably more
marginal-consistent than a random degree-preserving arrangement. But it is
below the practical margin, and the panel was powered to `0.00453`, so this is
a genuine effect-size failure and not a detection failure.

The mechanism is chemistry, not statistics. At the frozen `6.0 A` P1B contact
threshold, a slot holds about three **sequence-adjacent** residues, which are
also spatially adjacent along the backbone. If one of them is within `6.0 A` of
a ligand atom, its neighbours usually are too. The slot partition and the
contact threshold together leave very little to deconvolve.

Two controls are reported honestly rather than as evidence:

- **Atom shuffle is an exact no-op** — `0.985611`, identical to the empirical
  arm to every digit. The statistic ranks candidates by a residue-side column
  sum, which is invariant to permuting atom rows. It is a degenerate control
  for this estimand and carries no information.
- **Geometry shuffle scores higher** than the empirical arm, `0.993948`.
  Permuting slot assignment breaks up sequence-adjacent groupings, and
  arbitrary residue groups are *easier* to rank by degree than genuinely
  adjacent ones. This confirms that the real slot partition is the hard case,
  and that "beating the shuffle" would have been the wrong direction of test.

`C1c` passed on the letter of the rule but rests on only **17** replicate pairs
from distinct PDB entries across 5 components, with mean contacted-residue
Jaccard `0.830`. That denominator is positive but thin, so replicate-normalized
headroom is not reported as a headline quantity, exactly as the preregistration
required.

## 6. Verdict and what it forecloses

```text
EXACT_EDGE_COUPLING_NOT_SUPPORTED_BY_TEACHER
```

This is the earliest failed boundary. C0 identifiability passed, the slot
routing estimand is valid after amendment 01, and the units are abundant — so
the failure is squarely about the teacher's exact coupling content, which is
what the verdict names.

The C2 Geometry-Gated Coarse-to-Exact Correspondence Router was **not
preregistered and not trained**, as the registered stopping rule requires.
Nothing here authorizes attention stacks, a new PLM, a parallel GNN, typed
energy heads, orientation channels, affinity supervision, KG features or
adapters, and nothing authorizes widening the corpus, relaxing the `6.0 A`
threshold or re-running with a different closure to chase the margin.

This result is consistent with, and sharpens, the established fact that MONN
exact atom-residue coupling is weak relative to degree-preserving nulls. It
extends it to a fully independent, never-exposed corpus of 1,862 systems built
from raw RCSB coordinates, and it localises the cause: at a `6.0 A` threshold
with 128 linear sequence slots, exact correspondence is very nearly a function
of residue contact degree.

## 7. Governance and remaining frozen boundaries

- Trains nothing; zero trainable parameters introduced.
- Affinity value reads `0` from BioLiP2, PLINDER, PDBbind, ChEMBL, BindingDB,
  DAVIS, KIBA and every recipient dataset.
- Heldout-A is permanently consumed and was never referenced. Heldout-B was not
  created. R6 was not opened.
- `model/`, `scripts/`, `theory/`, CSMO, Band, the mesh, `z` and
  `A(F,z)=K(B(z)F(z))` are unmodified. Prior-stage code was imported for
  equivalence checking only.
- No selection rule, threshold, margin, seed, closure or exclusion rule was
  changed after any C0 census or C1 statistic was read. The one amendment
  preceded every statistic and corrected a provenance error, not a Gate.
- Affinity direction, few-shot sectioning and biological `z` remain
  unidentified and unopened.
